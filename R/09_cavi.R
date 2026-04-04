# ============================================================
# 09_cavi.R
#
# Mean-field variational pipeline for the recommended single-ordering
# MPCurver fit.
# Current scope:
#   * single ordering
#   * homoskedastic across components, but feature-specific sigma_j^2
#   * feature-specific lambda_j
#   * no relative_lambda scaling
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

.cavi_normalize_responsibilities <- function(R, n, K) {
  R <- as.matrix(R)
  if (!all(dim(R) == c(n, K))) {
    stop("responsibilities_init must be an n x K matrix.")
  }
  R[R < 0] <- 0
  rs <- rowSums(R)
  if (any(rs <= 0)) stop("Each row of responsibilities_init must have positive sum.")
  R / rs
}

.cavi_raw_moments_from_responsibilities <- function(X, R,
                                                    sigma_min = 1e-10,
                                                    sigma_max = 1e10) {
  X <- as.matrix(X)
  R <- as.matrix(R)
  n <- nrow(X)
  d <- ncol(X)
  Nk <- colSums(R)
  Nk[Nk < 1e-8] <- 1e-8

  GX <- crossprod(X, R)
  mu_raw <- sweep(GX, 2L, Nk, "/")
  GX2 <- crossprod(X^2, R)
  mu2_Nk <- mu_raw^2 * rep(Nk, each = d)
  resid_term <- rowSums(GX2 - 2 * mu_raw * GX + mu2_Nk)
  sigma2_raw <- pmin(pmax(resid_term / n, sigma_min), sigma_max)

  list(
    mu_raw = mu_raw,
    sigma2_raw = as.numeric(sigma2_raw),
    Nk = Nk
  )
}

.cavi_lambda_init_from_mu_raw <- function(mu_raw, Q_K, r_rank,
                                          lambda_sd_prior_rate = NULL,
                                          lambda_min = 1e-10,
                                          lambda_max = 1e10) {
  mu_raw <- as.matrix(mu_raw)
  eq_quad_raw <- vapply(seq_len(nrow(mu_raw)), function(j) {
    mj <- mu_raw[j, ]
    as.numeric(crossprod(mj, Q_K %*% mj))
  }, numeric(1))
  denom <- pmax(eq_quad_raw, 1e-12)
  lambda_sd_prior_rate_val <- .lambda_sd_prior_rate_value(lambda_sd_prior_rate)

  if (lambda_sd_prior_rate_val <= 0) {
    return(pmin(pmax(r_rank / denom, lambda_min), lambda_max))
  }

  .optimize_lambda_induced_exp(
    eq_quad = denom,
    r_rank = r_rank,
    rate = lambda_sd_prior_rate_val,
    lambda_min = lambda_min,
    lambda_max = lambda_max
  )
}


#' Simulate a toy dataset from the current single-ordering CAVI model
#'
#' @param n Number of observations.
#' @param d Number of features.
#' @param K Number of latent positions/components.
#' @param rw_q Random-walk order for the GMRF prior.
#' @param lambda_range Length-2 positive range for feature-specific \code{lambda_j}.
#' @param sigma_range Length-2 positive range for feature-specific \code{sigma_j}.
#' @param pi Optional component probabilities; defaults to uniform.
#' @param ridge Small ridge used only for simulation to make the prior proper.
#' @param seed Optional integer seed.
#'
#' @return A list with \code{X}, \code{z}, \code{pi}, \code{mu}, \code{sigma2},
#'   \code{lambda_vec}, and \code{Q_K}.
#' @export
simulate_cavi_toy <- function(n = 150,
                              d = 20,
                              K = 8,
                              rw_q = 2L,
                              lambda_range = c(0.5, 3),
                              sigma_range = c(0.08, 0.18),
                              pi = NULL,
                              ridge = 1e-3,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- as.integer(n)
  d <- as.integer(d)
  K <- as.integer(K)
  rw_q <- as.integer(rw_q)
  if (n < 2L) stop("n must be >= 2.")
  if (d < 1L) stop("d must be >= 1.")
  if (K < 2L) stop("K must be >= 2.")
  if (rw_q < 1L || rw_q >= K) stop("rw_q must be in {1, ..., K-1}.")

  lambda_range <- as.numeric(lambda_range)
  sigma_range <- as.numeric(sigma_range)
  if (length(lambda_range) != 2L || any(!is.finite(lambda_range)) || any(lambda_range <= 0)) {
    stop("lambda_range must contain two positive finite values.")
  }
  if (length(sigma_range) != 2L || any(!is.finite(sigma_range)) || any(sigma_range <= 0)) {
    stop("sigma_range must contain two positive finite values.")
  }

  if (is.null(pi)) {
    pi <- rep(1 / K, K)
  } else {
    pi <- as.numeric(pi)
    if (length(pi) != K || any(pi < 0) || sum(pi) <= 0) {
      stop("pi must be a nonnegative vector of length K with positive sum.")
    }
    pi <- pi / sum(pi)
  }

  Q_K <- make_random_walk_precision(K = K, d = 1, q = rw_q, lambda = 1, ridge = ridge)
  Q_K <- 0.5 * (as.matrix(Q_K) + t(as.matrix(Q_K)))

  lambda_vec <- exp(stats::runif(d, log(min(lambda_range)), log(max(lambda_range))))
  sigma2 <- exp(stats::runif(d, log(min(sigma_range)^2), log(max(sigma_range)^2)))

  mu_mat <- matrix(0, nrow = d, ncol = K)
  for (j in seq_len(d)) {
    A_j <- lambda_vec[j] * Q_K
    L_j <- chol(A_j)
    z_j <- stats::rnorm(K)
    mu_j <- backsolve(L_j, z_j)
    mu_j <- as.numeric(mu_j - mean(mu_j))
    mu_mat[j, ] <- mu_j
  }

  z <- sample.int(K, size = n, replace = TRUE, prob = pi)
  X <- matrix(0, nrow = n, ncol = d)
  for (i in seq_len(n)) {
    X[i, ] <- mu_mat[, z[i]] + stats::rnorm(d, sd = sqrt(sigma2))
  }

  list(
    X = X,
    z = z,
    pi = pi,
    mu = mu_mat,
    sigma2 = sigma2,
    lambda_vec = lambda_vec,
    Q_K = Q_K,
    rw_q = rw_q
  )
}


#' Fit the recommended CAVI model for a single MPCurver ordering
#'
#' @description
#' Recommended variational pipeline for a single-ordering MPCurver model with
#' feature-specific noise levels \eqn{\sigma_j^2} and smoothness parameters
#' \eqn{\lambda_j}. The latent trajectory for feature \eqn{j} is a length-\eqn{K}
#' vector with a GMRF prior along the component axis.
#'
#' This function is now the main single-ordering MPCurver backend. In contrast
#' to the legacy \code{csmooth_em} path, it keeps an explicit Gaussian
#' posterior over feature trajectories, with posterior means \eqn{m_j} and
#' covariances \eqn{S_j}.
#'
#' @param X Numeric matrix (\eqn{n \times d}).
#' @param K Number of components. If NULL, defaults to \code{max(2, min(50, floor(n/5)))}.
#' @param method Initialization ordering method when \code{responsibilities_init} is NULL.
#' @param responsibilities_init Optional initial responsibility matrix (\eqn{n \times K}).
#' @param pi_init Optional initial mixture proportions.
#' @param sigma2_init Optional initial feature-specific variances.
#' @param lambda_init Optional initial feature-specific smoothness vector.
#' @param rw_q Random-walk order for the prior precision along components.
#' @param ridge Ridge added to \code{Q_K}.
#' @param discretization Discretization method for ordering-based initialization.
#' @param fix_lambda Logical. If \code{TRUE}, skip variational lambda updates
#'   and keep \code{lambda_j} fixed at the initial values throughout fitting.
#' @param lambda_sd_prior_rate Optional positive rate for an exponential prior
#'   on \code{1 / sqrt(lambda_j)}. The default \code{NULL} means no lambda prior
#'   penalty. For backward compatibility, an explicit \code{0} is treated the
#'   same way; it is only an alias for "no penalty" and does not correspond to
#'   a literal exponential prior with rate zero.
#' @param lambda_min,lambda_max Bounds for \code{lambda_j}.
#' @param sigma_min,sigma_max Bounds for \code{sigma_j^2}.
#' @param max_iter Maximum number of CAVI sweeps. \code{0} returns the
#'   initialization-only state.
#' @param tol Relative ELBO tolerance.
#' @param verbose Logical.
#'
#' @return An object of class \code{"cavi"} with components including:
#' \itemize{
#'   \item \code{$params}: current \code{pi}, \code{mu}, and \code{sigma2}
#'   \item \code{$posterior}: Gaussian posterior summaries for feature trajectories
#'   \item \code{$gamma}: cell-by-component responsibilities
#'   \item \code{$lambda_vec}: feature-specific smoothness parameters
#'   \item \code{$elbo_trace}: the monotone variational objective used for diagnostics
#'   \item \code{$loglik_trace}: a plug-in observed log-likelihood diagnostic
#' }
#' @export
cavi <- function(X,
                 K = NULL,
                 method = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
                 responsibilities_init = NULL,
                 pi_init = NULL,
                 sigma2_init = NULL,
                 lambda_init = NULL,
                 rw_q = 2L,
                 ridge = 0,
                 discretization = c("quantile", "equal", "kmeans"),
                 fix_lambda = FALSE,
                 lambda_sd_prior_rate = NULL,
                 lambda_min = 1e-10,
                 lambda_max = 1e10,
                 sigma_min = 1e-10,
                 sigma_max = 1e10,
                 max_iter = 100L,
                 tol = 1e-6,
                 verbose = FALSE) {
  method <- match.arg(method)
  discretization <- match.arg(discretization)

  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  if (n < 2L) stop("X must have at least two rows.")
  if (d < 1L) stop("X must have at least one column.")

  if (is.null(K)) {
    K <- max(2L, min(50L, as.integer(floor(n / 5))))
  } else {
    K <- as.integer(K)
  }
  rw_q <- as.integer(rw_q)
  max_iter <- as.integer(max_iter)

  if (length(K) != 1L || is.na(K) || K < 2L) stop("K must be a single integer >= 2.")
  if (length(rw_q) != 1L || is.na(rw_q) || rw_q < 1L || rw_q >= K) {
    stop("rw_q must be an integer in {1, ..., K-1}.")
  }
  if (length(max_iter) != 1L || is.na(max_iter) || max_iter < 0L) {
    stop("max_iter must be an integer >= 0.")
  }

  lambda_min <- as.numeric(lambda_min)
  lambda_max <- as.numeric(lambda_max)
  lambda_sd_prior_rate <- .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
  lambda_sd_prior_rate_val <- .lambda_sd_prior_rate_value(lambda_sd_prior_rate)
  sigma_min <- as.numeric(sigma_min)
  sigma_max <- as.numeric(sigma_max)
  tol <- as.numeric(tol)
  if (!is.finite(lambda_min) || !is.finite(lambda_max) || lambda_min <= 0 || lambda_max <= 0 || lambda_min > lambda_max) {
    stop("lambda_min/lambda_max must be positive finite numbers with lambda_min <= lambda_max.")
  }
  if (!is.finite(sigma_min) || !is.finite(sigma_max) || sigma_min <= 0 || sigma_max <= 0 || sigma_min > sigma_max) {
    stop("sigma_min/sigma_max must be positive finite numbers with sigma_min <= sigma_max.")
  }
  if (!is.finite(tol) || tol < 0) stop("tol must be a nonnegative finite number.")

  Q_K <- make_random_walk_precision(K = K, d = 1, q = rw_q, lambda = 1, ridge = ridge)
  Q_K <- 0.5 * (as.matrix(Q_K) + t(as.matrix(Q_K)))
  prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
  r_rank <- prior_meta$rank
  logdet_Q <- prior_meta$logdet
  skip_raw_gamma_preinit <- FALSE
  raw_gamma_preinit <- FALSE
  raw_init <- NULL

  if (is.null(responsibilities_init)) {
    init <- initialize_ordering_csmooth(
      X = X,
      K = K,
      method = method,
      discretization = discretization,
      modelName = "homoskedastic"
    )
    params0 <- init$params
    K_init <- length(params0$pi)
    if (K_init != K) {
      K <- K_init
      if (rw_q >= K) {
        stop("Initialisation reduced K below rw_q + 1; try a smaller rw_q or different discretization.")
      }
      Q_K <- make_random_walk_precision(K = K, d = 1, q = rw_q, lambda = 1, ridge = ridge)
      Q_K <- 0.5 * (as.matrix(Q_K) + t(as.matrix(Q_K)))
      prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
      r_rank <- prior_meta$rank
      logdet_Q <- prior_meta$logdet
    }

    keep_idx <- init$keep_idx %||% seq_len(n)
    cluster_rank <- init$cluster_rank
    R <- matrix(1 / K, nrow = n, ncol = K)
    R[keep_idx, ] <- 0
    R[cbind(keep_idx, cluster_rank)] <- 1
  } else {
    skip_raw_gamma_preinit <- isTRUE(attr(responsibilities_init, "cavi_skip_raw_preinit"))
    R <- .cavi_normalize_responsibilities(responsibilities_init, n = n, K = K)
    params0 <- NULL
    raw_gamma_preinit <- !skip_raw_gamma_preinit
    if (raw_gamma_preinit) {
      raw_init <- .cavi_raw_moments_from_responsibilities(
        X = X,
        R = R,
        sigma_min = sigma_min,
        sigma_max = sigma_max
      )
    }
  }

  pi_vec <- if (!is.null(pi_init)) {
    tmp <- as.numeric(pi_init)
    if (length(tmp) != K || any(tmp < 0) || sum(tmp) <= 0) {
      stop("pi_init must be a nonnegative vector of length K with positive sum.")
    }
    tmp / sum(tmp)
  } else {
    pmax(colMeans(R), .Machine$double.eps)
  }
  pi_vec <- pi_vec / sum(pi_vec)

  sigma2 <- if (!is.null(sigma2_init)) {
    tmp <- as.numeric(sigma2_init)
    if (length(tmp) == 1L) tmp <- rep(tmp, d)
    if (length(tmp) != d) stop("sigma2_init must have length 1 or ncol(X).")
    pmin(pmax(tmp, sigma_min), sigma_max)
  } else if (raw_gamma_preinit) {
    raw_init$sigma2_raw
  } else if (!is.null(params0)) {
    pmin(pmax(as.numeric(params0$sigma2), sigma_min), sigma_max)
  } else {
    pmin(pmax(apply(X, 2, stats::var), sigma_min), sigma_max)
  }

  lambda_vec <- if (raw_gamma_preinit) {
    if (!is.null(lambda_init)) {
      tmp <- as.numeric(lambda_init)
      if (length(tmp) == 1L) tmp <- rep(tmp, d)
      if (length(tmp) != d) stop("lambda_init must have length 1 or ncol(X).")
      pmin(pmax(tmp, lambda_min), lambda_max)
    } else if (isTRUE(fix_lambda)) {
      rep(1, d)
    } else {
      .cavi_lambda_init_from_mu_raw(
        mu_raw = raw_init$mu_raw,
        Q_K = Q_K,
        r_rank = r_rank,
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        lambda_min = lambda_min,
        lambda_max = lambda_max
      )
    }
  } else if (!is.null(lambda_init)) {
    tmp <- as.numeric(lambda_init)
    if (length(tmp) == 1L) tmp <- rep(tmp, d)
    if (length(tmp) != d) stop("lambda_init must have length 1 or ncol(X).")
    pmin(pmax(tmp, lambda_min), lambda_max)
  } else if (!is.null(params0)) {
    mu_init_mat <- do.call(cbind, params0$mu)
    tmp <- numeric(d)
    for (j in seq_len(d)) {
      quad_j <- as.numeric(crossprod(mu_init_mat[j, ], Q_K %*% mu_init_mat[j, ]))
      tmp[j] <- r_rank / pmax(quad_j, 1e-12)
    }
    pmin(pmax(tmp, lambda_min), lambda_max)
  } else {
    rep(1, d)
  }

  .update_q_u <- function(R, sigma2, lambda_vec) {
    Nk <- colSums(R)
    Nk[Nk < 1e-8] <- 1e-8
    GX <- t(X) %*% R

    m_mat <- matrix(0, nrow = d, ncol = K)
    sdiag_mat <- matrix(0, nrow = d, ncol = K)
    S_list <- vector("list", d)
    logdetS <- numeric(d)
    eq_quad <- numeric(d)

    for (j in seq_len(d)) {
      A_j <- diag(as.numeric(Nk / sigma2[j]), K, K) + lambda_vec[j] * Q_K
      A_j <- 0.5 * (A_j + t(A_j))
      L_j <- chol(A_j)
      S_j <- chol2inv(L_j)
      rhs_j <- as.numeric(GX[j, ]) / sigma2[j]
      y_j <- forwardsolve(t(L_j), rhs_j)
      m_j <- as.numeric(backsolve(L_j, y_j))

      S_list[[j]] <- S_j
      m_mat[j, ] <- m_j
      sdiag_mat[j, ] <- diag(S_j)
      logdetS[j] <- -2 * sum(log(diag(L_j)))
      eq_quad[j] <- as.numeric(crossprod(m_j, Q_K %*% m_j) + sum(Q_K * S_j))
    }

    list(
      m_mat = m_mat,
      sdiag_mat = sdiag_mat,
      S_list = S_list,
      logdetS = logdetS,
      eq_quad = eq_quad
    )
  }

  .update_r <- function(q_u, pi_vec, sigma2) {
    m_mat     <- q_u$m_mat      # d x K
    sdiag_mat <- q_u$sdiag_mat  # d x K
    inv_s2    <- 1 / sigma2     # d-vector

    # Decompose quad_mat[i,k] = Σ_j [(x_ij - m_jk)² + S_jkk] / σ²_j
    Xs2   <- drop(X^2 %*% inv_s2)              # n-vec
    cross <- X %*% (inv_s2 * m_mat)             # n x K
    mu2_s <- colSums(inv_s2 * m_mat^2)          # K-vec
    unc_s <- colSums(inv_s2 * sdiag_mat)        # K-vec

    quad_mat <- Xs2 - 2 * cross + rep(mu2_s + unc_s, each = n)

    const_sigma <- sum(log(2 * base::pi * sigma2))
    log_r <- rep(log(pmax(pi_vec, 1e-300)), each = n) - 0.5 * (const_sigma + quad_mat)

    lse <- matrixStats::rowLogSumExps(log_r)
    exp(log_r - lse)
  }

  .update_sigma2 <- function(R, q_u) {
    m_mat     <- q_u$m_mat      # d x K
    sdiag_mat <- q_u$sdiag_mat  # d x K
    Nk <- colSums(R)            # K-vec

    # Expand Σ_i Σ_k R_ik (x_ij - m_jk)² = GX2 - 2·m·GX + m²·Nk
    GX2    <- crossprod(X^2, R)                      # d x K
    GX     <- crossprod(X, R)                         # d x K
    mu2_Nk <- m_mat^2 * rep(Nk, each = d)            # d x K

    resid_term <- rowSums(GX2 - 2 * m_mat * GX + mu2_Nk)  # d-vec
    unc_term   <- rowSums(sdiag_mat * rep(Nk, each = d))   # d-vec

    sigma2_new <- (resid_term + unc_term) / n
    pmin(pmax(sigma2_new, sigma_min), sigma_max)
  }

  .update_lambda <- function(q_u) {
    denom <- pmax(q_u$eq_quad, 1e-12)
    if (lambda_sd_prior_rate_val <= 0) {
      return(pmin(pmax(r_rank / denom, lambda_min), lambda_max))
    }
    .optimize_lambda_induced_exp(
      eq_quad = denom,
      r_rank = r_rank,
      rate = lambda_sd_prior_rate_val,
      lambda_min = lambda_min,
      lambda_max = lambda_max
    )
  }

  .compute_elbo <- function(R, pi_vec, sigma2, lambda_vec, q_u) {
    m_mat     <- q_u$m_mat      # d x K
    sdiag_mat <- q_u$sdiag_mat  # d x K
    inv_s2    <- 1 / sigma2     # d-vector

    log_pi_term <- sum(R * rep(log(pmax(pi_vec, 1e-300)), each = n))
    entropy_z   <- -sum(R * log(R + 1e-300))

    # Likelihood: same quad_mat decomposition as .update_r()
    Xs2      <- drop(X^2 %*% inv_s2)
    cross    <- X %*% (inv_s2 * m_mat)
    mu2_s    <- colSums(inv_s2 * m_mat^2)
    unc_s    <- colSums(inv_s2 * sdiag_mat)
    quad_mat <- Xs2 - 2 * cross + rep(mu2_s + unc_s, each = n)

    const_sigma <- sum(log(2 * base::pi * sigma2))
    lik_term <- -0.5 * sum(R * (const_sigma + quad_mat))

    # Prior and entropy(q_U) — vectorized over d
    prior_term <- 0.5 * sum(
      r_rank * (log(lambda_vec) - log(2 * base::pi)) +
        logdet_Q - lambda_vec * q_u$eq_quad
    )
    lambda_penalty <- sum(.lambda_sd_prior_terms(
      lambda_vec = lambda_vec,
      rate = lambda_sd_prior_rate,
      include_constant = TRUE
    ))
    entropy_u  <- 0.5 * sum(K * log(2 * base::pi * exp(1)) + q_u$logdetS)

    as.numeric(log_pi_term + entropy_z + lik_term + prior_term + lambda_penalty + entropy_u)
  }

  .compute_observed_loglik <- function(q_u, pi_vec, sigma2) {
    m_mat  <- q_u$m_mat   # d x K
    inv_s2 <- 1 / sigma2  # d-vector

    Xs2   <- drop(X^2 %*% inv_s2)
    cross <- X %*% (inv_s2 * m_mat)
    mu2_s <- colSums(inv_s2 * m_mat^2)

    quad_mat <- Xs2 - 2 * cross + rep(mu2_s, each = n)

    const_sigma <- sum(log(2 * base::pi * sigma2))
    logdens <- rep(log(pmax(pi_vec, 1e-300)), each = n) - 0.5 * (const_sigma + quad_mat)

    sum(matrixStats::rowLogSumExps(logdens))
  }

  q_u <- .update_q_u(R, sigma2, lambda_vec)
  elbo_trace <- .compute_elbo(R, pi_vec, sigma2, lambda_vec, q_u)
  loglik_trace <- .compute_observed_loglik(q_u, pi_vec, sigma2)
  lambda_trace <- list(lambda_vec)
  sigma2_trace <- list(sigma2)
  pi_trace <- list(pi_vec)
  converged <- FALSE

  if (max_iter > 0L) {
    for (iter in seq_len(max_iter)) {
      R <- .update_r(q_u, pi_vec, sigma2)
      pi_vec <- pmax(colMeans(R), .Machine$double.eps)
      pi_vec <- pi_vec / sum(pi_vec)
      sigma2 <- .update_sigma2(R, q_u)
      if (!fix_lambda) lambda_vec <- .update_lambda(q_u)
      q_u <- .update_q_u(R, sigma2, lambda_vec)

      elbo_new <- .compute_elbo(R, pi_vec, sigma2, lambda_vec, q_u)
      loglik_new <- .compute_observed_loglik(q_u, pi_vec, sigma2)
      elbo_trace <- c(elbo_trace, elbo_new)
      loglik_trace <- c(loglik_trace, loglik_new)
      lambda_trace[[length(lambda_trace) + 1L]] <- lambda_vec
      sigma2_trace[[length(sigma2_trace) + 1L]] <- sigma2
      pi_trace[[length(pi_trace) + 1L]] <- pi_vec

      delta <- elbo_trace[length(elbo_trace)] - elbo_trace[length(elbo_trace) - 1L]
      rel_delta <- delta / (abs(elbo_trace[length(elbo_trace) - 1L]) + 1)

      if (verbose) {
        cat(sprintf(
          "[cavi %3d] ELBO=%.6f  delta=%.3e  sigma2=[%.3g, %.3g]  lambda=[%.3g, %.3g]\n",
          iter, elbo_new, delta, min(sigma2), max(sigma2), min(lambda_vec), max(lambda_vec)
        ))
      }

      if (delta < -1e-8) {
        warning(
          sprintf("ELBO decreased by %.3e at iteration %d.", delta, iter),
          call. = FALSE
        )
      }

      if (delta >= 0 && abs(rel_delta) < tol) {
        converged <- TRUE
        break
      }
    }
  }

  params <- list(
    pi = pi_vec,
    mu = lapply(seq_len(K), function(k) q_u$m_mat[, k]),
    sigma2 = sigma2
  )

  structure(
    list(
      params = params,
      gamma = R,
      posterior = list(
        mean = q_u$m_mat,
        cov = q_u$S_list,
        var = q_u$sdiag_mat,
        diag = q_u$sdiag_mat
      ),
      lambda_vec = lambda_vec,
      Q_K = Q_K,
      rw_q = rw_q,
      elbo_trace = elbo_trace,
      loglik_trace = loglik_trace,
      ml_trace = numeric(0),
      lambda_trace = lambda_trace,
      sigma2_trace = sigma2_trace,
      pi_trace = pi_trace,
      iter = length(elbo_trace) - 1L,
      converged = converged,
      control = list(
        modelName = "homoskedastic",
        adaptive = if (fix_lambda) "fixed_lambda" else "variational",
        method = method,
        discretization = discretization,
        fix_lambda = fix_lambda,
        prior_proper = prior_meta$proper,
        prior_rank = prior_meta$rank,
        prior_logdet = prior_meta$logdet,
        gamma_preinit = if (is.null(responsibilities_init)) {
          "ordering_init"
        } else if (raw_gamma_preinit) {
          "raw_unpenalized"
        } else {
          "direct_penalized"
        },
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        sigma_min = sigma_min,
        sigma_max = sigma_max,
        max_iter = max_iter,
        tol = tol,
        ridge = ridge
      ),
      data = X
    ),
    class = "cavi"
  )
}


#' Continue CAVI sweeps on an existing \code{cavi} fit
#'
#' @param object A \code{cavi} object.
#' @param iter Integer maximum number of additional CAVI sweeps.
#' @param tol Optional convergence tolerance. If \code{NULL}, reuse the
#'   original fit's tolerance.
#' @param lambda Optional scalar or d-vector. If provided, overrides
#'   \code{lambda_vec} in the fit and fixes lambda (skips variational
#'   lambda updates).
#' @param lambda_min,lambda_max Optional bounds for \code{lambda_j}. If
#'   \code{NULL}, reuse the values stored in the original fit's
#'   \code{$control}. Ignored when \code{lambda} is provided.
#' @param lambda_sd_prior_rate Optional positive rate for the induced
#'   exponential prior on \code{1 / sqrt(lambda_j)}. If \code{NULL}, reuse the
#'   value stored in the original fit's \code{$control}. An explicit \code{0}
#'   is treated as "no penalty" for backward compatibility and does not
#'   represent a literal exponential prior with rate zero.
#' @param sigma_min,sigma_max Optional bounds for \code{sigma_j^2}. If
#'   \code{NULL}, reuse the values stored in the original fit's
#'   \code{$control}.
#' @param verbose Logical.
#'
#' @return An updated \code{cavi} object with extended traces. The returned
#'   object keeps the full accumulated \code{$elbo_trace}, \code{$loglik_trace},
#'   \code{$lambda_trace}, \code{$sigma2_trace}, and \code{$pi_trace}.
#' @export
do_cavi <- function(object,
                    iter = 1L,
                    tol = NULL,
                    lambda = NULL,
                    lambda_sd_prior_rate = NULL,
                    lambda_min = NULL,
                    lambda_max = NULL,
                    sigma_min = NULL,
                    sigma_max = NULL,
                    verbose = FALSE) {
  if (!inherits(object, "cavi")) stop("object must inherit from class 'cavi'.")

  iter <- as.integer(iter)
  if (length(iter) != 1L || is.na(iter) || iter < 0L) {
    stop("iter must be a single integer >= 0.")
  }

  tol_use <- if (is.null(tol)) {
    object$control$tol %||% 1e-6
  } else {
    as.numeric(tol)
  }

  if (!is.null(lambda)) {
    d <- ncol(object$data)
    lambda_init_use <- rep(as.numeric(lambda), length.out = d)
    fix_lambda_use <- TRUE
    lambda_sd_prior_rate_use <- if (is.null(lambda_sd_prior_rate)) {
      object$control$lambda_sd_prior_rate %||% NULL
    } else {
      .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
    }
    lmin_use <- 1e-10
    lmax_use <- 1e10
  } else {
    lambda_init_use <- object$lambda_vec
    fix_lambda_use <- isTRUE(object$control$fix_lambda)
    lambda_sd_prior_rate_use <- if (is.null(lambda_sd_prior_rate)) {
      object$control$lambda_sd_prior_rate %||% NULL
    } else {
      .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
    }
    lmin_use <- if (is.null(lambda_min)) object$control$lambda_min %||% 1e-10 else as.numeric(lambda_min)
    lmax_use <- if (is.null(lambda_max)) object$control$lambda_max %||% 1e10  else as.numeric(lambda_max)
  }
  smin_use <- if (is.null(sigma_min))  object$control$sigma_min  %||% 1e-10 else as.numeric(sigma_min)
  smax_use <- if (is.null(sigma_max))  object$control$sigma_max  %||% 1e10  else as.numeric(sigma_max)

  updated <- cavi(
    X = object$data,
    K = length(object$params$pi),
    method = object$control$method %||% "PCA",
    responsibilities_init = structure(
      object$gamma,
      cavi_skip_raw_preinit = TRUE
    ),
    pi_init = object$params$pi,
    sigma2_init = object$params$sigma2,
    lambda_init = lambda_init_use,
    rw_q = object$rw_q %||% 2L,
    ridge = object$control$ridge %||% 0,
    discretization = object$control$discretization %||% "quantile",
    fix_lambda = fix_lambda_use,
    lambda_sd_prior_rate = lambda_sd_prior_rate_use,
    lambda_min = lmin_use,
    lambda_max = lmax_use,
    sigma_min = smin_use,
    sigma_max = smax_use,
    max_iter = iter,
    tol = tol_use,
    verbose = verbose
  )

  if (length(updated$elbo_trace) > 0L) {
    updated$elbo_trace <- c(object$elbo_trace %||% numeric(0), updated$elbo_trace[-1L])
  } else {
    updated$elbo_trace <- object$elbo_trace %||% numeric(0)
  }

  if (length(updated$loglik_trace) > 0L) {
    updated$loglik_trace <- c(object$loglik_trace %||% numeric(0), updated$loglik_trace[-1L])
  } else {
    updated$loglik_trace <- object$loglik_trace %||% numeric(0)
  }

  if (length(updated$lambda_trace) > 0L) {
    updated$lambda_trace <- c(object$lambda_trace %||% list(), updated$lambda_trace[-1L])
  } else {
    updated$lambda_trace <- object$lambda_trace %||% list()
  }

  if (length(updated$sigma2_trace) > 0L) {
    updated$sigma2_trace <- c(object$sigma2_trace %||% list(), updated$sigma2_trace[-1L])
  } else {
    updated$sigma2_trace <- object$sigma2_trace %||% list()
  }

  if (length(updated$pi_trace) > 0L) {
    updated$pi_trace <- c(object$pi_trace %||% list(), updated$pi_trace[-1L])
  } else {
    updated$pi_trace <- object$pi_trace %||% list()
  }

  updated$iter <- length(updated$elbo_trace) - 1L
  updated
}


#' @export
print.cavi <- function(x, ...) {
  last_elbo <- if (length(x$elbo_trace)) tail(x$elbo_trace, 1L) else NA_real_
  last_ll <- if (length(x$loglik_trace)) tail(x$loglik_trace, 1L) else NA_real_
  cat("cavi fit\n")
  cat(sprintf("  n=%d, d=%d, K=%d\n", nrow(x$data), ncol(x$data), length(x$params$pi)))
  cat(sprintf("  iter=%d, converged=%s\n", x$iter, if (isTRUE(x$converged)) "yes" else "no"))
  cat(sprintf("  ELBO(last)=%.6f\n", last_elbo))
  if (is.finite(last_ll)) cat(sprintf("  logLik(last)=%.6f\n", last_ll))
  cat(sprintf("  sigma2 range=[%.3g, %.3g]\n", min(x$params$sigma2), max(x$params$sigma2)))
  cat(sprintf("  lambda range=[%.3g, %.3g]\n", min(x$lambda_vec), max(x$lambda_vec)))
  invisible(x)
}


#' Summary method for \code{cavi} objects
#'
#' @param object A \code{cavi} object.
#' @param ... Unused.
#'
#' @return An object of class \code{"summary.cavi"}.
#' @export
summary.cavi <- function(object, ...) {
  if (!inherits(object, "cavi")) stop("object must inherit from class 'cavi'.")

  elbo_trace <- object$elbo_trace %||% numeric(0)
  loglik_trace <- object$loglik_trace %||% numeric(0)
  elbo_last <- if (length(elbo_trace) > 0L) tail(elbo_trace, 1L) else NA_real_
  elbo_diff_last <- if (length(elbo_trace) >= 2L) tail(diff(elbo_trace), 1L) else NA_real_
  loglik_last <- if (length(loglik_trace) > 0L) tail(loglik_trace, 1L) else NA_real_
  loglik_diff_last <- if (length(loglik_trace) >= 2L) tail(diff(loglik_trace), 1L) else NA_real_

  rel_change <- function(v) {
    v <- as.numeric(v %||% numeric(0))
    if (length(v) < 2L) return(NA_real_)
    a <- v[length(v) - 1L]
    b <- v[length(v)]
    if (!is.finite(a) || !is.finite(b)) return(NA_real_)
    abs(b - a) / max(1e-12, abs(a))
  }

  rel_change_vec <- function(v) {
    if (is.null(v) || length(v) < 2L) return(NA_real_)
    a <- as.numeric(v[[length(v) - 1L]])
    b <- as.numeric(v[[length(v)]])
    if (length(a) != length(b) || length(a) == 0L) return(NA_real_)
    sqrt(sum((b - a)^2)) / max(1e-12, sqrt(sum(a^2)))
  }

  pi_hat <- object$params$pi %||% numeric(0)

  out <- list(
    n = if (!is.null(object$data)) nrow(object$data) else NA_integer_,
    d = if (!is.null(object$data)) ncol(object$data) else NA_integer_,
    K = length(object$params$pi %||% numeric(0)),
    iter = object$iter %||% (length(elbo_trace) - 1L),
    converged = isTRUE(object$converged),
    modelName = object$control$modelName %||% "homoskedastic",
    rw_q = object$rw_q %||% NA_integer_,
    init_method = object$control$method %||% NA_character_,
    discretization = object$control$discretization %||% NA_character_,
    adaptive = object$control$adaptive %||% "variational",
    elbo_trace = elbo_trace,
    loglik_trace = loglik_trace,
    elbo_last = elbo_last,
    elbo_diff_last = elbo_diff_last,
    loglik_last = loglik_last,
    loglik_diff_last = loglik_diff_last,
    elbo_rel_change = rel_change(elbo_trace),
    obj_rel_change = rel_change(loglik_trace),
    lambda_rel_change = rel_change_vec(object$lambda_trace),
    pi = pi_hat,
    pi_min = if (length(pi_hat)) min(pi_hat) else NA_real_,
    pi_max = if (length(pi_hat)) max(pi_hat) else NA_real_,
    Nk = if (!is.null(object$gamma)) colSums(object$gamma) else numeric(0),
    sigma2_type = "vector(d)",
    sigma2_range = range(object$params$sigma2, na.rm = TRUE),
    sigma2_min = min(object$params$sigma2, na.rm = TRUE),
    sigma2_max = max(object$params$sigma2, na.rm = TRUE),
    lambda_range = range(object$lambda_vec, na.rm = TRUE),
    lambda_min = min(object$lambda_vec, na.rm = TRUE),
    lambda_max = max(object$lambda_vec, na.rm = TRUE),
    lambda_mean = mean(object$lambda_vec, na.rm = TRUE)
  )

  class(out) <- "summary.cavi"
  out
}


#' @export
print.summary.cavi <- function(x, ...) {
  cat("<summary.cavi>\n")
  cat(sprintf("  n = %s, d = %s, K = %s\n",
              ifelse(is.na(x$n), "?", x$n),
              ifelse(is.na(x$d), "?", x$d),
              ifelse(is.na(x$K), "?", x$K)))
  cat(sprintf("  model = %s, RW(q) = %s\n",
              ifelse(is.na(x$modelName), "?", x$modelName),
              ifelse(is.na(x$rw_q), "?", x$rw_q)))
  cat(sprintf("  init = %s, discretization = %s, adaptive = %s\n",
              ifelse(is.na(x$init_method), "?", x$init_method),
              ifelse(is.na(x$discretization), "?", x$discretization),
              ifelse(is.na(x$adaptive), "?", x$adaptive)))
  cat(sprintf("  iter = %s, converged = %s\n",
              ifelse(is.na(x$iter), "?", x$iter),
              if (isTRUE(x$converged)) "yes" else "no"))
  cat(sprintf("  last ELBO = %s (delta last = %s)\n",
              ifelse(is.na(x$elbo_last), "?", format(x$elbo_last, digits = 8)),
              ifelse(is.na(x$elbo_diff_last), "?", format(x$elbo_diff_last, digits = 4))))
  cat(sprintf("  last logLik = %s (delta last = %s)\n",
              ifelse(is.na(x$loglik_last), "?", format(x$loglik_last, digits = 8)),
              ifelse(is.na(x$loglik_diff_last), "?", format(x$loglik_diff_last, digits = 4))))
  cat(sprintf("  pi range = [%s, %s]\n",
              format(x$pi_min, digits = 4),
              format(x$pi_max, digits = 4)))
  cat(sprintf("  sigma2 range = [%s, %s]\n",
              format(x$sigma2_range[1], digits = 4),
              format(x$sigma2_range[2], digits = 4)))
  cat(sprintf("  lambda range = [%s, %s]\n",
              format(x$lambda_range[1], digits = 4),
              format(x$lambda_range[2], digits = 4)))
  if (is.finite(x$elbo_rel_change) || is.finite(x$obj_rel_change) || is.finite(x$lambda_rel_change)) {
    cat("  relative change (last step):\n")
    if (is.finite(x$elbo_rel_change)) cat(sprintf("    ELBO: %.3g\n", x$elbo_rel_change))
    if (is.finite(x$obj_rel_change)) cat(sprintf("    logLik: %.3g\n", x$obj_rel_change))
    if (is.finite(x$lambda_rel_change)) cat(sprintf("    lambda_vec (L2 rel): %.3g\n", x$lambda_rel_change))
  }
  invisible(x)
}


#' Plot a \code{cavi} fit
#'
#' @param x A \code{cavi} object.
#' @param plot_type One of \code{"scatterplot"}, \code{"elbo"}, or \code{"mu"}.
#' @param dims Integer vector of length 1 or 2 used for plotting.
#' @param data Optional data matrix. Defaults to \code{x$data}.
#' @param two_panel Logical; for \code{plot_type = "elbo"}, show ELBO and
#'   plug-in log-likelihood in separate panels?
#' @param pal Colour palette for \code{"scatterplot"}.
#' @param add_legend Logical; draw the pseudotime legend on scatterplots?
#' @param ... Passed through to low-level plotting functions.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.cavi <- function(
    x,
    plot_type = c("scatterplot", "elbo", "mu"),
    dims = c(1L, 2L),
    data = NULL,
    two_panel = FALSE,
    pal = grDevices::colorRampPalette(
      c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"))(256L),
    add_legend = TRUE,
    ...
) {
  if (!inherits(x, "cavi")) stop("x must inherit from class 'cavi'.")
  plot_type <- match.arg(plot_type)

  if (plot_type == "elbo") {
    elbo <- x$elbo_trace %||% numeric(0)
    ll <- x$loglik_trace %||% numeric(0)

    if (two_panel) {
      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(oldpar), add = TRUE)
      graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
      graphics::plot(seq_along(elbo), elbo, type = "l", lwd = 2,
                     xlab = "Iteration", ylab = "ELBO",
                     main = sprintf("CAVI ELBO trace (K=%d)", length(x$params$pi)),
                     ...)
      if (length(ll) > 0L) {
        graphics::plot(seq_along(ll), ll, type = "l", lwd = 2,
                       xlab = "Iteration", ylab = "Observed log-likelihood",
                       main = "Plug-in observed log-likelihood", ...)
      } else {
        graphics::plot.new()
        graphics::title(main = "Plug-in observed log-likelihood (empty)")
      }
    } else {
      xlim <- c(1, max(length(elbo), length(ll)))
      yall <- c(elbo, ll)
      yall <- yall[is.finite(yall)]
      ylim <- if (length(yall)) range(yall) else c(-1, 1)

      graphics::plot.new()
      graphics::plot.window(xlim = xlim, ylim = ylim)
      graphics::axis(1)
      graphics::axis(2)
      graphics::box()
      graphics::title(main = sprintf("CAVI traces (K=%d)", length(x$params$pi)),
                      xlab = "Iteration", ylab = "Value")
      if (length(elbo) > 0L) graphics::lines(seq_along(elbo), elbo, lty = 1, lwd = 2, ...)
      if (length(ll) > 0L) graphics::lines(seq_along(ll), ll, lty = 2, lwd = 2)
      graphics::legend("bottomright",
                       legend = c(if (length(elbo) > 0L) "ELBO" else NULL,
                                  if (length(ll) > 0L) "logLik" else NULL),
                       lty = c(if (length(elbo) > 0L) 1 else NULL,
                               if (length(ll) > 0L) 2 else NULL),
                       lwd = 2, bty = "n")
    }
    return(invisible(x))
  }

  if (is.null(data)) {
    data <- x$data
    if (is.null(data)) stop("No data found in x$data. Please supply `data` explicitly.")
  }
  data <- as.matrix(data)

  dims <- as.integer(dims)
  if (length(dims) < 1L || length(dims) > 2L) stop("dims must have length 1 or 2.")
  if (any(dims < 1L) || any(dims > ncol(data))) stop("dims out of range for data.")

  K <- length(x$params$pi)
  mu_mat <- x$posterior$mean

  if (plot_type == "mu") {
    if (length(dims) == 2L) {
      graphics::plot(
        mu_mat[dims[1], ], mu_mat[dims[2], ],
        type = "o", pch = 16, col = "orange", lwd = 2,
        xlab = sprintf("dim %d", dims[1]),
        ylab = sprintf("dim %d", dims[2]),
        main = sprintf("CAVI posterior mean path (K=%d)", K),
        ...
      )
    } else {
      positions <- (seq_len(K) - 1L) / (K - 1L)
      graphics::plot(
        positions, mu_mat[dims[1], ],
        type = "o", pch = 16, col = "orange", lwd = 2,
        xlab = "pseudotime",
        ylab = sprintf("dim %d", dims[1]),
        main = sprintf("CAVI posterior mean path (K=%d)", K),
        ...
      )
    }
    return(invisible(x))
  }

  positions <- (seq_len(K) - 1L) / (K - 1L)
  t_pseudo <- as.numeric(x$gamma %*% positions)

  if (length(dims) == 2L) {
    mu_list_dims <- lapply(seq_len(K), function(k) mu_mat[dims, k])
    plot_EM_embedding2D(
      mu_list = mu_list_dims,
      X2 = data[, dims, drop = FALSE],
      t_vec = t_pseudo,
      pal = pal,
      add_legend = FALSE,
      main = sprintf("CAVI pseudotime (K=%d)", K),
      xlab = sprintf("dim %d", dims[1]),
      ylab = sprintf("dim %d", dims[2]),
      ...
    )
  } else {
    idx <- pmax(1L, pmin(length(pal), 1L + floor(t_pseudo * (length(pal) - 1L))))
    graphics::plot(
      t_pseudo, data[, dims[1]],
      pch = 19, col = pal[idx], cex = 0.7,
      xlab = "pseudotime",
      ylab = sprintf("dim %d", dims[1]),
      main = sprintf("CAVI pseudotime (K=%d)", K),
      ...
    )
    graphics::lines(positions, mu_mat[dims[1], ], col = "orange", lwd = 2)
    graphics::points(positions, mu_mat[dims[1], ], pch = 8, col = "orange", cex = 1)
  }

  if (isTRUE(add_legend) && exists(".draw_gradient_legend", mode = "function")) {
    .draw_gradient_legend(pal, title = "pseudotime", lo_label = "0", hi_label = "1")
  }

  invisible(x)
}
