# ============================================================
# 08_SoftSmoothEM.R
#
# Legacy compatibility layer for soft-partition dual-trajectory inference via
# annealed EM. New development should happen in the CAVI partition path in
# 10_partition_cavi.R.
# Depends on local MPCurver internals: score_feature_given_Gamma(),
# do_csmoothEM(), initialize_csmoothEM(), MSTEP_csmooth()
#
# Theory summary (v2 — mathematically correct)
# ---------------------------------------------
# The generative model has a feature-level latent indicator z_j in {A,B}.
# The EM Q-function for trajectory m is:
#
#   Q(mu_j^(m)) = w_{jm} * [ data_lik_j + log P(mu_j^(m) | lambda_j) ]
#
# Taking d/d mu_j: the scalar w_{jm} cancels in the MAP update for mu_j.
# => The mu-refresh step is identical to the unweighted case.
#
# For adaptive = "ml", however, the inner hyper-step optimizes the collapsed
# objective after integrating out mu_j. The Laplace normalization term
# contributes -1/2 log|A_j| without the leading w_{jm}, so lambda_j / sigma_j^2
# generally DO depend on the feature weights.
#
# The soft inner update therefore has two pieces:
#   1. weighted E-step:
#        log gamma_{ik} += sum_j  w_{jm} * log N(x_ij | mu_kj, sigma2_j)
#   2. weighted collapsed hyper-step for lambda / sigma when adaptive = "ml"
#      (mu-refresh remains unchanged).
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b


# ============================================================
# 1.  ESTEP_csmooth_weighted
# ============================================================

#' Weighted E-step for csmooth_em
#'
#' Identical to ESTEP_csmooth() except per-feature log-likelihood contributions
#' are multiplied by feature_weights[j] before summing over features.
#' In the adaptive = "ml" path, this E-step is paired with a weighted
#' collapsed hyper-step for lambda / sigma2.
#'
#' @param X               n-by-d data matrix.
#' @param params          csmooth params list (pi, mu, sigma2).
#' @param modelName       "homoskedastic" or "heteroskedastic".
#' @param feature_weights length-d non-negative weight vector.
#' @return n-by-K responsibility matrix.
#' @keywords internal
ESTEP_csmooth_weighted <- function(X, params, modelName, feature_weights) {

  X  <- as.matrix(X)
  n  <- nrow(X); d <- ncol(X)
  K  <- length(params$pi)

  if (missing(feature_weights) || is.null(feature_weights))
    feature_weights <- rep(1, d)
  w <- pmax(as.numeric(feature_weights), 0)
  if (length(w) != d) stop("feature_weights must have length d = ncol(X).")

  Mu     <- do.call(cbind, params$mu)  # d x K
  log_pi <- log(pmax(params$pi, .Machine$double.eps))

  if (modelName == "homoskedastic") {
    sigma2   <- pmax(as.numeric(params$sigma2), .Machine$double.eps)
    inv_s2_w <- w / sigma2             # length-d: weight / variance

    # quad_w[i,k] = sum_j w_j*(x_ij - mu_jk)^2/sigma2_j
    Mu_s  <- Mu * inv_s2_w             # d x K
    xs2   <- drop(X^2 %*% inv_s2_w)   # n-vector
    cross <- X %*% Mu_s                # n x K
    mu2   <- colSums(Mu * Mu_s)        # K-vector

    quad_w  <- xs2 - 2 * cross + matrix(mu2, n, K, byrow = TRUE)
    const_w <- -0.5 * sum(w * (log(2 * pi) + log(sigma2)))
    log_gamma <- const_w - 0.5 * quad_w + matrix(log_pi, n, K, byrow = TRUE)

  } else {
    sigma2_mat  <- pmax(as.matrix(params$sigma2), .Machine$double.eps)  # d x K
    inv_s2_w_mat <- w / sigma2_mat     # d x K  (broadcast w over columns)

    # quad_w[i,k] = sum_j w_j*(x_ij - mu_jk)^2/sigma2_jk
    Mu_s  <- Mu * inv_s2_w_mat         # d x K
    xs2   <- X^2 %*% inv_s2_w_mat      # n x K
    cross <- X %*% Mu_s                # n x K
    mu2   <- colSums(Mu * Mu_s)        # K-vector

    quad_w     <- xs2 - 2 * cross + matrix(mu2, n, K, byrow = TRUE)
    const_w_vec <- -0.5 * colSums(w * (log(2 * pi) + log(sigma2_mat)))  # K-vector
    log_gamma  <- -0.5 * quad_w + matrix(const_w_vec + log_pi, n, K, byrow = TRUE)
  }

  lse <- matrixStats::rowLogSumExps(log_gamma)
  exp(log_gamma - lse)
}


# ============================================================
# 2.  do_csmoothEM_weighted
# ============================================================

#' Run csmoothEM with weighted inner updates
#'
#' Identical interface to do_csmoothEM(), with one extra argument:
#' \code{feature_weights}. The inner loop reuses the same csmoothEM updates,
#' but replaces the responsibility update with
#' \code{ESTEP_csmooth_weighted()} so all inner iterations use the correct
#' weighted log-likelihood. When \code{adaptive = "ml"}, the collapsed
#' \code{lambda}/\code{sigma2} hyper-step is also reweighted; the final MAP
#' refresh of \code{mu} is unchanged.
#'
#' @param object          A csmooth_em object.
#' @param data            n-by-d data matrix.
#' @param feature_weights length-d weight vector in [0,1].
#' @param iter            Number of inner EM iterations.
#' @param adaptive        Adaptive lambda mode (\code{"ml"} or NULL).
#'   \code{"prior"} is obsolete and retained only for backward compatibility.
#' @param sigma_update    Sigma update mode for \code{adaptive = "ml"}.
#'   Defaults to \code{"ml"}; legacy \code{"mstep"} is retained only for backward compatibility.
#' @param lambda_min,lambda_max  Lambda bounds.
#' @param record          Logical; record traces.
#' @param verbose         Logical.
#' @return Updated csmooth_em object.
#' @noRd
do_csmoothEM_weighted <- function(object,
                                  data,
                                  feature_weights,
                                  iter         = 3L,
                                  adaptive     = "ml",
                                  sigma_update = c("ml", "mstep"),
                                  lambda_min   = 1e-10,
                                  lambda_max   = 1e10,
                                  record       = TRUE,
                                  verbose      = FALSE) {

  if (!inherits(object, "csmooth_em"))
    stop("object must be a 'csmooth_em' object.")

  sigma_update <- .normalize_sigma_update(
    sigma_update = sigma_update,
    default = "ml",
    caller = "do_csmoothEM_weighted()"
  )
  adaptive <- .normalize_csmooth_adaptive(
    adaptive = adaptive,
    default = object$control$adaptive %||% "ml",
    none_as_null = TRUE,
    caller = "do_csmoothEM_weighted()"
  )

  data <- as.matrix(data)
  d    <- ncol(data)
  w    <- pmax(as.numeric(feature_weights), 0)
  if (length(w) != d) stop("feature_weights must have length ncol(data).")

  if (!is.null(adaptive) && identical(adaptive, "ml")) {
    return(do_csmoothEM_ml_collapsed_weighted(
      object = object,
      data = data,
      feature_weights = w,
      iter = iter,
      record = record,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_update = sigma_update,
      verbose = verbose
    ))
  }

  X <- data
  modelName       <- object$control$modelName %||% "homoskedastic"
  relative_lambda <- isTRUE(object$control$relative_lambda)
  nugget          <- object$control$nugget %||% 0
  eigen_tol       <- object$control$eigen_tol
  rw_q            <- as.integer(object$prior$rw_q %||% 0L)

  Q_K        <- object$prior$Q_K
  lambda_vec <- as.numeric(object$prior$lambda_vec)

  K <- length(object$params$pi)
  r_rank <- max(K - rw_q, 1L)
  eps_quad <- 1e-12

  params <- cache_csmooth_params(object$params, modelName)

  .compute_logdetH <- function(Gamma, params, lambda_vec) {
    Nk <- colSums(Gamma)
    Nk[Nk < 1e-8] <- 1e-8

    if (modelName == "homoskedastic") {
      sigma2 <- pmax(as.numeric(params$sigma2), .Machine$double.eps)
    } else {
      sigma2 <- pmax(as.matrix(params$sigma2), .Machine$double.eps)
    }

    out <- 0
    for (j in seq_len(d)) {
      lam <- lambda_vec[j]
      if (!is.finite(lam) || lam <= 0) next

      Qb <- .compute_Qbase_j(params$pi, sigma2, Q_K, relative_lambda, modelName, j)

      if (modelName == "homoskedastic") {
        Dj <- Nk / sigma2[j]
      } else {
        Dj <- Nk / as.numeric(sigma2[j, ])
      }

      A <- diag(as.numeric(Dj), K, K) + lam * Qb
      A <- 0.5 * (A + t(A))
      L <- chol(A)
      out <- out + 2 * sum(log(diag(L)))
    }
    out
  }

  .compute_C <- function(X, Gamma, params, lambda_vec, include_constant = TRUE) {
    pelbo <- compute_penalized_ELBO_csmooth(
      X = X, Gamma = Gamma, params = params, Q_K = Q_K, lambda_vec = lambda_vec,
      modelName = modelName, relative_lambda = relative_lambda,
      eigen_tol = eigen_tol, rw_q = rw_q
    )
    ldH <- .compute_logdetH(Gamma, params, lambda_vec)
    const <- if (isTRUE(include_constant)) 0.5 * d * K * log(2 * base::pi) else 0
    list(C = pelbo + const - 0.5 * ldH, logdetH = ldH, pelbo = pelbo)
  }

  .mu_refresh <- function(X, Gamma, params, lambda_vec) {
    Nk <- colSums(Gamma)
    Nk[Nk < 1e-8] <- 1e-8
    pi_new <- Nk / sum(Nk)

    GX <- t(X) %*% Gamma

    if (modelName == "homoskedastic") {
      sigma2 <- pmax(as.numeric(params$sigma2), .Machine$double.eps)
    } else {
      sigma2 <- pmax(as.matrix(params$sigma2), .Machine$double.eps)
    }

    mu_mat <- matrix(0, nrow = d, ncol = K)

    for (j in seq_len(d)) {
      lam <- lambda_vec[j]
      if (!is.finite(lam) || lam < 0) lam <- 0

      if (modelName == "homoskedastic") {
        s2j <- sigma2[j]
        rhs <- as.numeric(GX[j, ]) / s2j
        A_diag <- Nk / s2j

        if (lam <= 0 || is.null(Q_K)) {
          mu_mat[j, ] <- as.numeric(GX[j, ]) / Nk
        } else {
          Qb <- .compute_Qbase_j(pi_new, sigma2, Q_K, relative_lambda, modelName, j)
          A <- diag(A_diag, K, K) + lam * Qb
          A <- 0.5 * (A + t(A))
          L <- chol(A)
          y <- forwardsolve(t(L), rhs)
          mu_mat[j, ] <- as.numeric(backsolve(L, y))
        }
      } else {
        s2jk <- as.numeric(sigma2[j, ])
        rhs <- as.numeric(GX[j, ]) / s2jk
        A_diag <- Nk / s2jk

        if (lam <= 0 || is.null(Q_K)) {
          mu_mat[j, ] <- as.numeric(GX[j, ]) / Nk
        } else {
          Qb <- .compute_Qbase_j(pi_new, sigma2, Q_K, relative_lambda, modelName, j)
          A <- diag(A_diag, K, K) + lam * Qb
          A <- 0.5 * (A + t(A))
          L <- chol(A)
          y <- forwardsolve(t(L), rhs)
          mu_mat[j, ] <- as.numeric(backsolve(L, y))
        }
      }
    }

    params$pi <- pi_new
    params$mu <- lapply(seq_len(K), function(k) mu_mat[, k])
    params
  }

  for (tt in seq_len(iter)) {
    Gamma <- ESTEP_csmooth_weighted(
      X = X,
      params = params,
      modelName = modelName,
      feature_weights = w
    )

    params <- MSTEP_csmooth(
      X = X, Gamma = Gamma, params = params,
      Q_K = Q_K, lambda_vec = lambda_vec,
      modelName = modelName, relative_lambda = relative_lambda,
      nugget = nugget, rw_q = rw_q, iterate_once = TRUE
    )
    params <- cache_csmooth_params(params, modelName)

    if (!is.null(adaptive)) {
      Nk <- colSums(Gamma)
      Nk[Nk < 1e-8] <- 1e-8
      GX <- t(X) %*% Gamma

      if (adaptive == "prior") {
        mu_mat <- do.call(cbind, params$mu)
        for (j in seq_len(d)) {
          Qb <- .compute_Qbase_j(params$pi, params$sigma2, Q_K, relative_lambda, modelName, j)
          mj <- as.numeric(mu_mat[j, ])
          quad_prior <- as.numeric(crossprod(mj, Qb %*% mj))
          lam_star <- r_rank / pmax(quad_prior, eps_quad)
          lambda_vec[j] <- min(max(lam_star, lambda_min), lambda_max)
        }
      } else if (adaptive == "ml") {
        stop("Internal error: adaptive = 'ml' should have dispatched to do_csmoothEM_ml_collapsed_weighted().")
      }

      params <- .mu_refresh(X, Gamma, params, lambda_vec)
    }

    if (record) {
      ll <- compute_log_joint_observed_csmooth(
        X, params, Q_K, lambda_vec,
        modelName = modelName,
        relative_lambda = relative_lambda,
        eigen_tol = eigen_tol,
        rw_q = rw_q
      )
      pelbo <- compute_penalized_ELBO_csmooth(
        X, Gamma, params, Q_K, lambda_vec,
        modelName = modelName,
        relative_lambda = relative_lambda,
        eigen_tol = eigen_tol,
        rw_q = rw_q
      )
      Cres <- .compute_C(X, Gamma, params, lambda_vec, include_constant = TRUE)

      object$loglik_trace <- c(object$loglik_trace, ll)
      object$elbo_trace <- c(object$elbo_trace, pelbo)

      if (is.null(object$ml_trace)) object$ml_trace <- numeric(0)
      object$ml_trace <- c(object$ml_trace, Cres$C)

      if (is.null(object$logdetH_trace)) object$logdetH_trace <- numeric(0)
      object$logdetH_trace <- c(object$logdetH_trace, Cres$logdetH)

      object$iter <- length(object$elbo_trace)
      if (is.null(object$lambda_trace)) object$lambda_trace <- list()
      object$lambda_trace[[object$iter]] <- lambda_vec

      if (verbose) {
        msg_ad <- if (is.null(adaptive)) "none" else adaptive
        cat(sprintf(
          "csmoothEM_weighted %d/%d: penLogLik=%.6f  penELBO=%.6f  C=%.6f  log|H|=%.6f  lambda(range)=[%.3g, %.3g] (adaptive=%s)\n",
          tt, iter, ll, pelbo, Cres$C, Cres$logdetH, min(lambda_vec), max(lambda_vec), msg_ad
        ))
      }
    }

    object$gamma <- Gamma
  }

  object$params <- params
  object$prior$lambda_vec <- lambda_vec
  object$data <- X
  object
}


#' Collapsed-ML csmoothEM with weighted E-step and weighted hyper-step
#' @keywords internal
do_csmoothEM_ml_collapsed_weighted <- function(object,
                                               data,
                                               feature_weights,
                                               iter = 1,
                                               record = TRUE,
                                               lambda_min = NULL,
                                               lambda_max = NULL,
                                               sigma_update = c("ml", "mstep"),
                                               verbose = FALSE) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!inherits(object, "csmooth_em")) stop("object must be a 'csmooth_em' object.")

  X <- as.matrix(data)
  iter <- as.integer(iter)
  if (length(iter) != 1L || is.na(iter) || iter < 1L) stop("iter must be integer >= 1.")

  modelName <- object$control$modelName %||% "homoskedastic"
  if (modelName != "homoskedastic") {
    stop("do_csmoothEM_ml_collapsed_weighted() currently supports modelName='homoskedastic' only.")
  }

  relative_lambda <- isTRUE(object$control$relative_lambda)
  nugget          <- object$control$nugget %||% 0
  eigen_tol       <- object$control$eigen_tol
  rw_q            <- as.integer(object$prior$rw_q %||% 0L)

  Q_K        <- object$prior$Q_K
  lambda_vec <- as.numeric(object$prior$lambda_vec)

  d <- ncol(X)
  K <- length(object$params$pi)
  w <- pmax(as.numeric(feature_weights), sqrt(.Machine$double.eps))
  if (length(w) != d) stop("feature_weights must have length ncol(data).")

  if (is.null(lambda_min)) lambda_min <- object$control$lambda_min %||% 1e-10
  if (is.null(lambda_max)) lambda_max <- object$control$lambda_max %||% 1e10
  lambda_min <- as.numeric(lambda_min)
  lambda_max <- as.numeric(lambda_max)
  if (!is.finite(lambda_min) || !is.finite(lambda_max) || lambda_min <= 0 || lambda_max <= 0 || lambda_min > lambda_max) {
    stop("lambda_min/lambda_max must be positive finite numbers with lambda_min <= lambda_max.")
  }

  sigma_update <- .normalize_sigma_update(
    sigma_update = sigma_update,
    default = "ml",
    caller = "do_csmoothEM_ml_collapsed_weighted()"
  )
  sigma_min <- as.numeric(object$control$sigma_min %||% 1e-10)
  sigma_max <- as.numeric(object$control$sigma_max %||% 1e10)
  if (!is.finite(sigma_min) || !is.finite(sigma_max) || sigma_min <= 0 || sigma_max <= 0 || sigma_min > sigma_max) {
    stop("sigma_min/sigma_max must be positive finite numbers with sigma_min <= sigma_max.")
  }

  r_rank <- max(K - rw_q, 1L)
  params <- object$params
  params <- cache_csmooth_params(params, "homoskedastic")

  .weighted_collapse_eval <- function(Gamma, params, lambda_vec, include_constant = TRUE) {
    res <- compute_C_by_coord_csmooth(
      X = X,
      Gamma = Gamma,
      params = params,
      Q_K = Q_K,
      lambda_vec = lambda_vec,
      modelName = "homoskedastic",
      relative_lambda = relative_lambda,
      eigen_tol = eigen_tol,
      rw_q = rw_q,
      include_constant = include_constant
    )
    const_j <- if (isTRUE(include_constant)) 0.5 * K * log(2 * base::pi) else 0
    pelbo_coord <- res$C_coord - const_j + 0.5 * res$logdetH_coord
    weighted_pelbo <- res$global_logpi + res$global_entropy + sum(w * pelbo_coord)
    weighted_C <- res$global_logpi + res$global_entropy +
      sum(w * res$C_coord + 0.5 * (w - 1) * res$logdetH_coord)

    list(
      C = as.numeric(weighted_C),
      pelbo = as.numeric(weighted_pelbo),
      logdetH = as.numeric(sum(res$logdetH_coord))
    )
  }

  .ell_j_loglam_weighted <- function(loglam, j, Nk, GX, GX2, pi_vec, sigma2, Q_K, r_rank, wj) {
    lam <- exp(loglam)
    s2j <- pmax(as.numeric(sigma2[j]), .Machine$double.eps)

    Qb <- .compute_Qbase_j(pi_vec, sigma2, Q_K, relative_lambda, "homoskedastic", j)
    C_lik <- -0.5 * sum(Nk) * log(2 * base::pi * s2j) -
      0.5 * sum(as.numeric(GX2[j, ]) / s2j)

    logdet_Qb <- generalized_logdet(Q = Qb, eigen_tol = eigen_tol, rank_deficiency = rw_q)

    Dj <- Nk / s2j
    bj <- as.numeric(GX[j, ]) / s2j

    A <- diag(as.numeric(Dj), K, K) + lam * Qb
    A <- 0.5 * (A + t(A))
    L <- chol(A)
    logdetA <- 2 * sum(log(diag(L)))
    y <- forwardsolve(t(L), bj)
    alpha <- backsolve(L, y)
    quad <- sum(bj * alpha)

    as.numeric(
      wj * (C_lik + 0.5 * (r_rank * loglam + logdet_Qb) + 0.5 * quad) -
        0.5 * logdetA
    )
  }

  .ell_j_logsig2_weighted <- function(logsig2, j, Nk, GX, GX2, pi_vec, lambda_vec, Q_K, r_rank, wj) {
    s2j <- exp(logsig2)

    sigma2_tmp <- pmax(as.numeric(params$sigma2), .Machine$double.eps)
    sigma2_tmp[j] <- s2j
    Qb <- .compute_Qbase_j(pi_vec, sigma2_tmp, Q_K, relative_lambda, "homoskedastic", j)

    C_lik <- -0.5 * sum(Nk) * log(2 * base::pi * s2j) -
      0.5 * sum(as.numeric(GX2[j, ]) / s2j)

    logdet_Qb <- generalized_logdet(Q = Qb, eigen_tol = eigen_tol, rank_deficiency = rw_q)

    Dj <- Nk / s2j
    bj <- as.numeric(GX[j, ]) / s2j
    lam <- lambda_vec[j]

    A <- diag(as.numeric(Dj), K, K) + lam * Qb
    A <- 0.5 * (A + t(A))
    L <- chol(A)
    logdetA <- 2 * sum(log(diag(L)))
    y <- forwardsolve(t(L), bj)
    alpha <- backsolve(L, y)
    quad <- sum(bj * alpha)

    as.numeric(
      wj * (C_lik + 0.5 * (r_rank * log(lam) + logdet_Qb) + 0.5 * quad) -
        0.5 * logdetA
    )
  }

  for (tt in seq_len(iter)) {
    Gamma <- ESTEP_csmooth_weighted(X, params, "homoskedastic", feature_weights = w)

    Nk <- colSums(Gamma)
    Nk[Nk < 1e-8] <- 1e-8
    GX  <- t(X) %*% Gamma
    GX2 <- t(X^2) %*% Gamma

    params$pi <- Nk / sum(Nk)

    if (sigma_update == "mstep") {
      tmp <- MSTEP_csmooth(
        X = X, Gamma = Gamma, params = params,
        Q_K = Q_K, lambda_vec = lambda_vec,
        modelName = "homoskedastic", relative_lambda = relative_lambda,
        nugget = nugget, rw_q = rw_q, iterate_once = TRUE
      )
      params$sigma2 <- tmp$sigma2
    } else {
      sigma2_new <- pmax(as.numeric(params$sigma2), .Machine$double.eps)
      for (j in seq_len(d)) {
        opt_s2 <- optimize(
          f = .ell_j_logsig2_weighted,
          interval = log(c(sigma_min, sigma_max)),
          maximum = TRUE,
          j = j, Nk = Nk, GX = GX, GX2 = GX2,
          pi_vec = params$pi, lambda_vec = lambda_vec,
          Q_K = Q_K, r_rank = r_rank, wj = w[j]
        )
        sigma2_new[j] <- exp(opt_s2$maximum)
      }
      params$sigma2 <- pmax(sigma2_new + nugget, .Machine$double.eps)
    }

    for (j in seq_len(d)) {
      opt_lam <- optimize(
        f = .ell_j_loglam_weighted,
        interval = log(c(lambda_min, lambda_max)),
        maximum = TRUE,
        j = j, Nk = Nk, GX = GX, GX2 = GX2,
        pi_vec = params$pi, sigma2 = params$sigma2,
        Q_K = Q_K, r_rank = r_rank, wj = w[j]
      )
      lambda_vec[j] <- exp(opt_lam$maximum)
    }

    params <- .refresh_map_mu_csmooth(
      X = X, Gamma = Gamma, params = params,
      Q_K = Q_K, lambda_vec = lambda_vec,
      modelName = "homoskedastic", relative_lambda = relative_lambda
    )
    params <- cache_csmooth_params(params, "homoskedastic")

    if (record) {
      Cres <- .weighted_collapse_eval(Gamma, params, lambda_vec, include_constant = TRUE)

      object$loglik_trace <- c(object$loglik_trace, NA_real_)
      object$elbo_trace   <- c(object$elbo_trace, Cres$pelbo)
      if (is.null(object$ml_trace)) object$ml_trace <- numeric(0)
      object$ml_trace <- c(object$ml_trace, Cres$C)
      if (is.null(object$logdetH_trace)) object$logdetH_trace <- numeric(0)
      object$logdetH_trace <- c(object$logdetH_trace, Cres$logdetH)

      object$iter <- length(object$elbo_trace)
      if (is.null(object$lambda_trace)) object$lambda_trace <- list()
      object$lambda_trace[[object$iter]] <- lambda_vec

      if (verbose) {
        cat(sprintf(
          "csmoothEM_ml_collapsed_weighted %d/%d: wC=%.6f  wELBO=%.6f  log|H|=%.6f  lambda(range)=[%.3g, %.3g]  sigma2(range)=[%.3g, %.3g]  sigma_update=%s\n",
          tt, iter, Cres$C, Cres$pelbo, Cres$logdetH,
          min(lambda_vec), max(lambda_vec),
          min(params$sigma2), max(params$sigma2),
          sigma_update
        ))
      }
    }

    object$gamma <- Gamma
  }

  object$params <- params
  object$prior$lambda_vec <- lambda_vec
  object$data <- X
  object
}


# ============================================================
# 3.  Internal helpers
# ============================================================

# Pure evaluation of S_jm using already-optimised inner-EM parameters.
# From the CAVI derivation (soft_partition_nestedEM.pdf, Step 1), S_jm is the
# collapsed marginal of feature j under trajectory m, evaluated at the current
# Gamma_m and the parameters (mu, sigma2, lambda) already found by the inner EM.
# Re-fitting these parameters here (as score_feature_given_Gamma does) would
# break the CAVI coordinate-ascent consistency.
.score_all_features <- function(X, fit,
                                rw_q = 2L, score_mode = "ml") {
  # All parameters are read from the fit object itself — this is the only
  # consistent source of truth after the inner EM has run.
  res <- compute_C_by_coord_csmooth(
    X              = X,
    Gamma          = fit$gamma,
    params         = fit$params,
    Q_K            = fit$prior$Q_K,
    lambda_vec     = as.numeric(fit$prior$lambda_vec),
    modelName      = fit$control$modelName      %||% "homoskedastic",
    relative_lambda = isTRUE(fit$control$relative_lambda),
    eigen_tol      = fit$control$eigen_tol,
    rw_q           = as.integer(fit$prior$rw_q  %||% rw_q),
    include_constant = TRUE
  )

  if (score_mode == "ml") {
    # C_coord already includes: Q_data_j + prior_j + K/2*log(2pi) - 1/2*log|A_j|
    as.numeric(res$C_coord)
  } else {
    # score_mode == "none": remove the Laplace correction terms
    # C_coord = Q_data_j + prior_j + const_j - 0.5*logdetA
    # => Q_data_j + prior_j = C_coord - const_j + 0.5*logdetH_coord
    const_j <- 0.5 * ncol(fit$gamma) * log(2 * base::pi)
    as.numeric(res$C_coord - const_j + 0.5 * res$logdetH_coord)
  }
}

# Compute the Gamma-structure terms of the full ELBO for one trajectory fit:
#   H(Gamma_m) + sum_k N_k * log(alpha_k^(m))
# For the collapsed outer objective used in soft partitioning, the full score is:
#   F = L_total + gamma_elbo_terms(fit1) + gamma_elbo_terms(fit2)
# In practice the inner trajectory blocks are only updated for a finite
# number of weighted EM iterations (`inner_iter`), so F is a useful tracking
# target but not a strict monotonicity certificate.
.gamma_elbo_terms <- function(fit) {
  G    <- fit$gamma
  Nk   <- colSums(G)
  pi_k <- pmax(as.numeric(fit$params$pi), .Machine$double.eps)
  entropy <- -sum(G * log(G + 1e-300))
  logpi   <- sum(Nk * log(pi_k))
  entropy + logpi
}

.softmax_rows <- function(M) {
  M <- M - apply(M, 1, max)
  E <- exp(M)
  E / rowSums(E)
}

.finite_floor <- function(s) {
  f <- s[is.finite(s)]
  if (!length(f)) return(rep(-1e6, length(s)))
  s[!is.finite(s)] <- min(f) - 1
  s
}


# ============================================================
# 4.  init_two_trajectories  (initialisation wrapper)
# ============================================================

#' Initialise two trajectory fits guaranteed to face different directions
#'
#' @description
#' Wrapper that produces a pair of \code{csmooth_em} objects suitable as
#' \code{fit1_init} / \code{fit2_init} for \code{soft_two_trajectory_EM()}.
#' Three strategies are available:
#'
#' \describe{
#'   \item{\strong{"score"}}{fit1 is initialised with PCA on all features.
#'     fit2 is seeded from the single feature that fit1 explains least
#'     (lowest per-feature marginal score).}
#'   \item{\strong{"mincor"}}{fit1 is initialised via \code{method1} (default PCA).
#'     fit2 is seeded from the feature with the smallest absolute correlation
#'     with the feature that drives fit1 most (top PC1 loading). Ensures
#'     fit2 starts from a direction maximally orthogonal to fit1.}
#'   \item{\strong{"random_split"}}{Features are randomly split 50/50.
#'     fit1 is trained on the first half, fit2 on the second half; both are
#'     then re-expanded to the full feature set.  Stochastic: set \code{seed}
#'     for reproducibility.}
#' }
#'
#' @param X          Numeric matrix (n x d).
#' @param method     Initialisation strategy: "score" (default), "mincor", or "random_split".
#' @param method1    Method for fit1 ordering: "PCA" (default), "fiedler", "tSNE", "pcurve", "random".
#' @param modelName  Variance structure: "homoskedastic" (default) or "heteroskedastic".
#' @param K          Number of pseudotime bins. If NULL (default), auto-selected as
#'   \code{max(2, min(50, floor(n/5)))}.
#' @param rw_q       Random-walk order for the smoothness prior (default 2).
#' @param adaptive   Adaptive hyperparameter mode passed to inner EM: \code{"ml"} (default)
#'   or \code{"none"}. \code{"prior"} is obsolete and retained only for backward compatibility.
#' @param relative_lambda Scale smoothness prior by feature variance (default TRUE).
#' @param num_iter   Warm-start EM iterations for each fit (default 5).
#' @param lambda_min,lambda_max  Lambda search bounds.
#' @param seed       Integer seed (used by "random_split"; ignored otherwise).
#' @param verbose    Logical.
#'
#' @return List with \code{fit1}, \code{fit2} (both \code{csmooth_em}),
#'   and \code{seed2} (the column index used to seed fit2, or NA for "random_split").
#' @noRd
init_two_trajectories <- function(X,
                                  method          = c("score", "mincor", "random_split"),
                                  method1         = c("PCA", "fiedler", "tSNE", "pcurve", "random"),
                                  modelName       = c("homoskedastic", "heteroskedastic"),
                                  K               = NULL,
                                  rw_q            = 2L,
                                  adaptive        = "ml",
                                  relative_lambda = TRUE,
                                  num_iter        = 5L,
                                  lambda_min      = 1e-10,
                                  lambda_max      = 1e10,
                                  seed            = 42L,
                                  verbose         = FALSE) {
  method    <- match.arg(method)
  method1   <- match.arg(method1)
  modelName <- match.arg(modelName)
  adaptive  <- .normalize_csmooth_adaptive(
    adaptive = adaptive,
    default = "ml",
    caller = "init_two_trajectories()"
  )
  X  <- as.matrix(X)
  n  <- nrow(X); d <- ncol(X)

  # Resolve K: user-supplied, or auto (min(50, floor(n/5)), at least 2)
  K_use <- if (!is.null(K)) as.integer(K) else max(2L, min(50L, as.integer(floor(n / 5))))

  # ---------- shared helper: build a csmooth_em from a single seed column ----
  # seed_col determines the initial ordering; full X is used for all params.
  .fit_from_seed <- function(seed_col, num_iter) {
    ord  <- rank(X[, seed_col], ties.method = "first")

    # use full X so mu/sigma2 have length d throughout
    init <- make_init_csmooth(
      X              = X,
      ordering_vec   = ord,
      K              = K_use,
      modelName      = modelName,
      discretization = "quantile",
      na_action      = "drop",
      eps            = 1e-12
    )
    Q_K <- make_random_walk_precision(K = length(init$pi), d = 1,
                                      lambda = 1, q = rw_q, ridge = 0)
    fit <- as_csmooth_em(
      params     = list(pi = init$pi, mu = init$mu, sigma2 = init$sigma2),
      data       = X,
      Q_K        = Q_K,
      lambda_vec = rep(1, d),
      rw_q       = rw_q,
      modelName  = modelName,
      relative_lambda = relative_lambda
    )
    # explicitly set adaptive in control so do_csmoothEM inherits it correctly
    fit$control$adaptive   <- adaptive
    fit$control$lambda_min <- lambda_min
    fit$control$lambda_max <- lambda_max

    do_csmoothEM(fit, data = X, iter = num_iter,
                 adaptive   = adaptive,
                 lambda_min = lambda_min,
                 lambda_max = lambda_max,
                 record = TRUE, verbose = FALSE)
  }

  # ---------- fit1: user-specified method (default PCA) ---------------------
  if (verbose) cat(sprintf("[init] fit1: %s init...\n", method1))
  fit1 <- initialize_csmoothEM(
    X = X, method = method1, modelName = modelName,
    K = K_use, adaptive = adaptive, num_iter = num_iter,
    relative_lambda = relative_lambda,
    lambda_min = lambda_min, lambda_max = lambda_max,
    include.data = TRUE
  )
  fit1 <- do_csmoothEM(fit1, data = X, iter = 1L,
                       record = TRUE, verbose = FALSE)

  # ---------- fit2: strategy-dependent seed ----------------------------------
  seed2 <- NA_integer_

  if (method == "score") {
    # Seed = feature fit1 explains least
    if (verbose) cat("[init] fit2: score-based seed...
")
    scores1 <- .score_all_features(X, fit1,
                                   rw_q = rw_q, score_mode = "ml")
    seed2 <- which.min(scores1)
    fit2  <- .fit_from_seed(seed2, num_iter)

  } else if (method == "mincor") {
    # Seed = feature with smallest absolute correlation with the fit1 seed feature.
    # Rationale: we want fit2 to start from a direction maximally orthogonal to
    # fit1, not just negatively correlated (which could still share structure).
    if (verbose) cat("[init] fit2: mincor seed...\n")
    # identify which feature drives fit1 most (highest abs loading on PC1 of X)
    pc1_loadings <- abs(prcomp(X, center = TRUE, scale. = FALSE)$rotation[, 1])
    seed1_col    <- which.max(pc1_loadings)
    # find the feature with the SMALLEST absolute correlation to seed1_col
    # (excluding seed1_col itself)
    abs_cors        <- abs(cor(X[, seed1_col], X))   # 1 x d
    abs_cors[seed1_col] <- Inf                        # exclude self
    seed2 <- which.min(abs_cors)
    if (verbose)
      cat(sprintf("[init]   seed1_col=%d  seed2_col=%d  |cor|=%.3f\n",
                  seed1_col, seed2, abs_cors[seed2]))
    fit2 <- .fit_from_seed(seed2, num_iter)

  } else {
    # random_split: train each fit on a random half, expand params to full d
    if (verbose) cat("[init] fit2: random-split init...
")
    set.seed(seed)
    perm   <- sample(d)
    half1  <- sort(perm[seq_len(floor(d / 2))])
    half2  <- sort(perm[(floor(d / 2) + 1):d])

    # fit each half independently
    fit1_half <- initialize_csmoothEM(
      X = X[, half1, drop = FALSE], method = "PCA",
      modelName = modelName, K = K_use, adaptive = adaptive,
      num_iter = num_iter, relative_lambda = relative_lambda,
      lambda_min = lambda_min, lambda_max = lambda_max, include.data = TRUE
    )
    fit2_half <- initialize_csmoothEM(
      X = X[, half2, drop = FALSE], method = "PCA",
      modelName = modelName, K = K_use, adaptive = adaptive,
      num_iter = num_iter, relative_lambda = relative_lambda,
      lambda_min = lambda_min, lambda_max = lambda_max, include.data = TRUE
    )

    # expand mu and sigma2 to full d by inserting neutral values for missing cols
    .expand_params <- function(fit_half, cols_used, d_full) {
      K       <- length(fit_half$params$pi)
      mu_half <- do.call(cbind, fit_half$params$mu)   # d_half x K

      mu_full <- matrix(0, nrow = d_full, ncol = K)
      mu_full[cols_used, ] <- mu_half
      missing <- setdiff(seq_len(d_full), cols_used)
      if (length(missing) > 0)
        mu_full[missing, ] <- matrix(
          rep(colMeans(X[, missing, drop = FALSE]), K),
          nrow = length(missing), ncol = K)

      if (modelName == "homoskedastic") {
        sig_half <- as.numeric(fit_half$params$sigma2)       # d_half
        sig_full <- rep(1, d_full)
        sig_full[cols_used] <- sig_half
        if (length(missing) > 0)
          sig_full[missing] <- apply(X[, missing, drop = FALSE], 2, var)
      } else {
        # heteroskedastic: sigma2 is d_half x K
        sig_half <- as.matrix(fit_half$params$sigma2)        # d_half x K
        sig_full <- matrix(1, nrow = d_full, ncol = K)
        sig_full[cols_used, ] <- sig_half
        if (length(missing) > 0)
          sig_full[missing, ] <- matrix(
            rep(apply(X[, missing, drop = FALSE], 2, var), K),
            nrow = length(missing), ncol = K)
      }

      list(
        pi     = fit_half$params$pi,
        mu     = lapply(seq_len(K), function(k) mu_full[, k]),
        sigma2 = sig_full
      )
    }

    params1 <- .expand_params(fit1_half, half1, d)
    params2 <- .expand_params(fit2_half, half2, d)

    K1 <- length(params1$pi); K2 <- length(params2$pi)
    Q_K1 <- make_random_walk_precision(K=K1, d=1, lambda=1, q=rw_q, ridge=0)
    Q_K2 <- make_random_walk_precision(K=K2, d=1, lambda=1, q=rw_q, ridge=0)

    fit1 <- as_csmooth_em(params=params1, data=X, Q_K=Q_K1,
                          lambda_vec=rep(1,d), rw_q=rw_q,
                          modelName=modelName, relative_lambda=relative_lambda)
    fit2 <- as_csmooth_em(params=params2, data=X, Q_K=Q_K2,
                          lambda_vec=rep(1,d), rw_q=rw_q,
                          modelName=modelName, relative_lambda=relative_lambda)

    # set adaptive in control so do_csmoothEM inherits it correctly
    fit1$control$adaptive <- adaptive; fit1$control$lambda_min <- lambda_min; fit1$control$lambda_max <- lambda_max
    fit2$control$adaptive <- adaptive; fit2$control$lambda_min <- lambda_min; fit2$control$lambda_max <- lambda_max

    fit1 <- do_csmoothEM(fit1, data=X, iter=num_iter, adaptive=adaptive,
                         lambda_min=lambda_min, lambda_max=lambda_max,
                         record=TRUE, verbose=FALSE)
    fit2 <- do_csmoothEM(fit2, data=X, iter=num_iter, adaptive=adaptive,
                         lambda_min=lambda_min, lambda_max=lambda_max,
                         record=TRUE, verbose=FALSE)
  }

  list(fit1 = fit1, fit2 = fit2, seed2 = seed2)
}


# ============================================================
# 4.  soft_two_trajectory_EM  (main exported function)
# ============================================================

#' Soft dual-trajectory partition via annealed CAVI
#'
#' @description
#' Jointly infers two latent 1-D orderings and a soft feature-level partition
#' via coordinate-ascent variational inference (CAVI). The algorithm alternates
#' between:
#'
#' Outer E-step: score each feature under current Gamma_1, Gamma_2 via
#'   \code{compute_C_by_coord_csmooth()}; compute \code{Softmax(S/T)} to get
#'   soft feature weights \code{w_jm}.
#'
#' Inner EM: refit each trajectory via \code{do_csmoothEM_weighted()}.
#'   The E-step uses weighted log-likelihood contributions, and when
#'   \code{adaptive = "ml"} the collapsed \code{lambda}/\code{sigma2}
#'   hyper-step is reweighted as well. The final MAP refresh of \code{mu}
#'   is unchanged.
#'
#' Phase 1 anneals the temperature T from \code{T_start} to \code{T_end} to
#' escape local optima. Each outer step then performs \code{inner_iter}
#' weighted trajectory updates before recomputing feature weights. Phase 2
#' runs the same untempered updates at \code{T = T_end} until the recorded
#' objective stabilises.
#'
#' This is the legacy \code{csmooth_em}-based soft partition routine. It is
#' retained for benchmark comparison and internal regression testing. For new
#' analyses, prefer \code{\link{soft_two_trajectory_cavi}}, which uses the
#' explicit variational \code{cavi} model throughout and is the recommended
#' dual-ordering soft-partition backend.
#'
#' @param X               Numeric matrix (n x d).
#' @param fit1_init       Pre-built \code{csmooth_em} for trajectory A.
#'   If NULL, auto-initialised via \code{init_method}.
#' @param fit2_init       Pre-built \code{csmooth_em} for trajectory B.
#'   If NULL, auto-initialised via \code{init_method}.
#' @param init_method     Initialisation strategy when fits are NULL:
#'   "score" (default), "mincor", or "random_split".
#' @param init_method1    Ordering method for fit1: "PCA" (default), "fiedler",
#'   "tSNE", "pcurve", or "random".
#' @param init_seed       Integer seed for "random_split" (default 42).
#' @param modelName       Variance structure: "homoskedastic" (default) or
#'   "heteroskedastic".
#' @param K               Number of pseudotime bins. If NULL (default),
#'   auto-selected as \code{max(2, min(50, floor(n/5)))}. Ignored when
#'   \code{fit1_init}/\code{fit2_init} are supplied.
#' @param T_start         Starting annealing temperature (default 5).
#' @param T_end           Final temperature; also the Phase 2 temperature
#'   (default 1.0).
#' @param n_outer         Number of Phase 1 annealing steps (default 25).
#' @param inner_iter      Weighted inner EM sweeps per outer step (default 1).
#' @param max_converge_iter  Maximum Phase 2 iterations (default 100).
#' @param tol_outer       Relative ELBO convergence tolerance for Phase 2
#'   (default 1e-4). Checked as |Delta F| / (|F| + 1).
#' @param score_mode      Outer scoring mode: "ml" (collapsed marginal, default)
#'   or "none" (penalised ELBO only).
#' @param rw_q            Random-walk order for the smoothness prior (default 2).
#' @param relative_lambda Scale smoothness prior by feature variance (default TRUE).
#' @param lambda_min,lambda_max  Lambda search bounds.
#' @param adaptive        Adaptive hyperparameter mode for the inner EM:
#'   \code{"ml"} (collapsed marginal, default) or \code{"none"}.
#'   \code{"prior"} is obsolete and retained only for backward compatibility.
#' @param sigma_update    Sigma estimation method inside inner EM when
#'   \code{adaptive = "ml"}: \code{"ml"} (collapsed marginal, default) or
#'   legacy \code{"mstep"} (weighted SSE). Use \code{"ml"} to jointly optimize
#'   \eqn{\sigma^2} and \eqn{\lambda}; \code{"mstep"} is retained only for backward compatibility.
#' @param hard_assign_final  Snap final soft weights to hard 0/1 (default FALSE).
#' @param verbose         Print iteration log (default TRUE).
#'
#' @return A list with:
#' \describe{
#'   \item{pi_weights}{d x 2 matrix of soft feature assignments for trajectories A and B.}
#'   \item{assign}{Character vector of hard assignments ("A" or "B").}
#'   \item{fit1, fit2}{Final \code{csmooth_em} objects for each trajectory.}
#'   \item{ll_history}{Full ELBO F at each post-M-step iteration.}
#'   \item{score_history, weight_history}{Per-iteration scores and weights.}
#'   \item{T_schedule}{Annealing temperature schedule.}
#'   \item{converged}{Logical: did Phase 2 converge?}
#'   \item{n_anneal}{Number of Phase 1 iterations.}
#' }
#' @seealso \code{\link{soft_two_trajectory_cavi}}
#' @noRd
soft_two_trajectory_EM <- function(
    X,
    fit1_init           = NULL,
    fit2_init           = NULL,
    init_method         = c("score", "mincor", "random_split"),
    init_method1        = c("PCA", "fiedler", "tSNE", "pcurve", "random"),
    init_seed           = 42L,
    modelName           = c("homoskedastic", "heteroskedastic"),
    K                   = NULL,
    T_start             = 5,
    T_end               = 1.0,
    n_outer             = 25,
    inner_iter          = 1L,
    max_converge_iter   = 100L,
    tol_outer           = 1e-4,
    score_mode          = c("ml", "none"),
    rw_q                = 2L,
    relative_lambda     = TRUE,
    lambda_min          = 1e-10,
    lambda_max          = 1e10,
    adaptive            = "ml",
    sigma_update        = c("ml", "mstep"),
    hard_assign_final   = FALSE,
    verbose             = TRUE
) {
  score_mode   <- match.arg(score_mode)
  init_method  <- match.arg(init_method)
  init_method1 <- match.arg(init_method1)
  modelName    <- match.arg(modelName)
  sigma_update <- .normalize_sigma_update(
    sigma_update = sigma_update,
    default = "ml",
    caller = "soft_two_trajectory_EM()"
  )
  adaptive     <- .normalize_csmooth_adaptive(
    adaptive = adaptive,
    default = "ml",
    caller = "soft_two_trajectory_EM()"
  )
  inner_iter   <- as.integer(inner_iter)
  if (length(inner_iter) != 1L || is.na(inner_iter) || inner_iter < 1L)
    stop("inner_iter must be an integer >= 1.")
  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)

  # ---- (0) Auto-initialisation ----------------------------------------------
  if (is.null(fit1_init) || is.null(fit2_init)) {
    if (verbose)
      cat(sprintf("[soft_EM] Auto-init: fit1=%s, fit2-seed=%s\n",
                  init_method1, init_method))
    inits <- init_two_trajectories(
      X               = X,
      method          = init_method,
      method1         = init_method1,
      modelName       = modelName,
      K               = K,
      rw_q            = rw_q,
      adaptive        = adaptive,
      relative_lambda = relative_lambda,
      num_iter        = 5L,
      lambda_min      = lambda_min,
      lambda_max      = lambda_max,
      seed            = init_seed,
      verbose         = verbose
    )
    if (is.null(fit1_init)) fit1_init <- inits$fit1
    if (is.null(fit2_init)) fit2_init <- inits$fit2
  }

  stopifnot(inherits(fit1_init, "csmooth_em"), inherits(fit2_init, "csmooth_em"))

  fit1 <- fit1_init; fit2 <- fit2_init

  if (is.null(fit1$gamma))
    fit1 <- do_csmoothEM(fit1, data = X, iter = 1L, record = TRUE, verbose = FALSE)
  if (is.null(fit2$gamma))
    fit2 <- do_csmoothEM(fit2, data = X, iter = 1L, record = TRUE, verbose = FALSE)

  # ---- (1) Annealing schedule -----------------------------------------------
  T_schedule <- exp(seq(log(T_start), log(T_end), length.out = n_outer))

  # ---- (2) Uniform initial weights ------------------------------------------
  pi_weights <- matrix(0.5, nrow = d, ncol = 2L)
  col_names  <- c("A", "B")
  colnames(pi_weights) <- col_names

  total_iters    <- n_outer + max_converge_iter
  score_history  <- vector("list", total_iters)
  weight_history <- vector("list", total_iters)
  ll_history     <- numeric(total_iters)

  # ---- helpers: full ELBO (CAVI §4) ----------------------------------------
  # F = L_total + H(Gamma_1) + logpi_1 + H(Gamma_2) + logpi_2
  # where L_total = sum_j log( sum_m (1/M) exp(S_jm) )
  #
  # IMPORTANT: S_jm must be evaluated with the CURRENT fit parameters
  # (post-M-step), and Gamma terms also come from the CURRENT fit.
  # Both must be from the same state. With finite inner_iter this trace is
  # informative but not guaranteed to be strictly monotone iteration-by-iteration.
  .L_total <- function(score_mat) {
    log_prior <- log(1 / ncol(score_mat))
    sum(matrixStats::rowLogSumExps(score_mat + log_prior))
  }
  .global_elbo <- function(score_mat, fit1, fit2) {
    .L_total(score_mat) +
      .gamma_elbo_terms(fit1) +
      .gamma_elbo_terms(fit2)
  }

  # ---- helper: one soft-partition refinement step ----------------------------
  # Given current fit1, fit2, temperature T:
  #   1. Outer E-step: compute S_jm, w_jm = softmax(S/T)
  #   2. Inner EM:     update Gamma and params for each trajectory
  #   3. Post-M-step:  recompute S_jm with NEW params, refresh w_jm, compute F
  # Returns: list(fit1, fit2, pi_weights, score_mat, elbo)
  # The recorded score is always based on the post-update state.
  .cavi_step <- function(fit1, fit2, T_now) {

    # Step 1: Outer E-step with PRE-M-step parameters
    S1 <- .finite_floor(.score_all_features(X, fit1,
            rw_q = rw_q, score_mode = score_mode))
    S2 <- .finite_floor(.score_all_features(X, fit2,
            rw_q = rw_q, score_mode = score_mode))
    score_mat_pre <- cbind(S1, S2)
    pi_weights    <- .softmax_rows(score_mat_pre / T_now)
    colnames(pi_weights) <- col_names

    # Step 2: Inner EM (weighted E-step + M-step) for each trajectory
    fit1 <- do_csmoothEM_weighted(
      fit1,
      data = X,
      feature_weights = pi_weights[, 1],
      iter = inner_iter,
      adaptive = adaptive,
      sigma_update = sigma_update,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      record = TRUE,
      verbose = FALSE
    )
    fit2 <- do_csmoothEM_weighted(
      fit2,
      data = X,
      feature_weights = pi_weights[, 2],
      iter = inner_iter,
      adaptive = adaptive,
      sigma_update = sigma_update,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      record = TRUE,
      verbose = FALSE
    )

    # Step 3: Recompute S_jm with POST-M-step parameters and refresh the
    # outer weights so the returned state is internally consistent.
    S1_new <- .finite_floor(.score_all_features(X, fit1,
                rw_q = rw_q, score_mode = score_mode))
    S2_new <- .finite_floor(.score_all_features(X, fit2,
                rw_q = rw_q, score_mode = score_mode))
    score_mat_post <- cbind(S1_new, S2_new)
    pi_weights_post <- .softmax_rows(score_mat_post / T_now)
    colnames(pi_weights_post) <- col_names
    elbo           <- .global_elbo(score_mat_post, fit1, fit2)

    list(fit1       = fit1,
         fit2       = fit2,
         pi_weights = pi_weights_post,
         score_mat  = score_mat_post,   # post-M-step scores, consistent with ELBO
         elbo       = elbo)
  }

  # ---- (3) Phase 1: annealing -----------------------------------------------
  if (verbose)
    cat(sprintf("[soft_EM] Phase 1: annealing  T: %.2f -> %.4f  (%d steps)\n",
                T_start, T_end, n_outer))

  actual_iters <- 0L
  converged    <- FALSE

  for (iter in seq_len(n_outer)) {
    T_now  <- T_schedule[iter]
    step   <- .cavi_step(fit1, fit2, T_now)
    fit1   <- step$fit1;  fit2 <- step$fit2
    pi_weights <- step$pi_weights

    score_history[[iter]]  <- step$score_mat
    weight_history[[iter]] <- pi_weights
    ll_history[iter]       <- step$elbo
    actual_iters           <- iter

    if (verbose) {
      ent   <- mean(-rowSums(pi_weights * log(pi_weights + 1e-300)))
      cat(sprintf("[phase1 %2d/%2d]  T=%.4f  |wA>0.5|=%d  |wB>0.5|=%d  mean_H=%.3f  ELBO=%.2f\n",
          iter, n_outer, T_now,
          sum(pi_weights[,1] > 0.5), sum(pi_weights[,2] > 0.5), ent, step$elbo))
    }
  }

  # ---- (4) Phase 2: refinement at T = T_end until convergence ---------------
  if (verbose)
    cat(sprintf("[soft_EM] Phase 2: refinement at T=%.4f  (max %d iters, tol=%.1e)\n",
                T_end, max_converge_iter, tol_outer))

  elbo_prev <- if (actual_iters > 0) ll_history[actual_iters] else -Inf

  for (iter2 in seq_len(max_converge_iter)) {
    iter   <- n_outer + iter2
    step   <- .cavi_step(fit1, fit2, T_end)
    fit1   <- step$fit1;  fit2 <- step$fit2
    pi_weights <- step$pi_weights

    score_history[[iter]]  <- step$score_mat
    weight_history[[iter]] <- pi_weights
    ll_history[iter]       <- step$elbo
    actual_iters           <- iter

    if (verbose) {
      ent <- mean(-rowSums(pi_weights * log(pi_weights + 1e-300)))
      cat(sprintf("[phase2 %3d    ]  T=%.4f  |wA>0.5|=%d  |wB>0.5|=%d  mean_H=%.3f  ELBO=%.2f\n",
          iter2, T_end,
          sum(pi_weights[,1] > 0.5), sum(pi_weights[,2] > 0.5), ent, step$elbo))
    }

    delta <- step$elbo - elbo_prev
    if (is.finite(elbo_prev) &&
        delta >= 0 &&
        abs(delta) / (abs(elbo_prev) + 1) < tol_outer) {
      if (verbose)
        cat(sprintf("[soft_EM] Phase 2 converged at iter %d  (DeltaELBO_rel=%.2e)\n",
                    iter2, abs(delta) / (abs(elbo_prev) + 1)))
      converged <- TRUE
      break
    }
    elbo_prev <- step$elbo
  }

  # truncate pre-allocated history to actual length
  score_history  <- score_history[seq_len(actual_iters)]
  weight_history <- weight_history[seq_len(actual_iters)]
  ll_history     <- ll_history[seq_len(actual_iters)]

  # ---- (5) Optional hard final assignment -----------------------------------
  if (hard_assign_final) {
    hard <- matrix(0, d, ncol(pi_weights))
    for (i in seq_len(d)) hard[i, which.max(pi_weights[i, ])] <- 1
    pi_weights <- hard; colnames(pi_weights) <- col_names
  }

  list(pi_weights     = pi_weights,
       assign         = col_names[apply(pi_weights, 1, which.max)],
       fit1           = fit1,
       fit2           = fit2,
       score_history  = score_history,
       weight_history = weight_history,
       ll_history     = ll_history,
       T_schedule     = T_schedule,
       n_anneal       = n_outer,
       converged      = converged)
}


# ============================================================
# 5.  Diagnostics
# ============================================================

#' Summarise a soft partition result
#' @param result  Output of soft_two_trajectory_EM().
#' @param feature_names  Optional character vector length d.
#' @return data.frame with feature, w_A, w_B, assign, entropy.
#' @noRd
summarise_soft_partition <- function(result, feature_names = NULL) {
  pw  <- result$pi_weights
  d   <- nrow(pw)
  fn  <- feature_names %||% paste0("V", seq_len(d))
  ent <- -rowSums(pw * log(pw + 1e-300))
  df  <- as.data.frame(pw)
  df$feature <- fn; df$assign <- result$assign; df$entropy <- round(ent, 4)
  df[, c("feature", colnames(pw), "assign", "entropy")]
}

#' Plot soft feature weights (4-panel diagnostic)
#' @param result          Output of soft_two_trajectory_EM().
#' @param feature_names   Optional character vector.
#' @param show_convergence  Show convergence panels 2, 3 & 4.
#' @param sort_by_weight  Sort features by w_A descending in bar chart.
#' @param title           Plot title.
#' @noRd
plot_soft_weights <- function(result, feature_names = NULL,
                              show_convergence = TRUE,
                              sort_by_weight   = TRUE,
                              title = "Soft feature partition") {
  pw        <- result$pi_weights
  d         <- nrow(pw)
  fn        <- feature_names %||% paste0("V", seq_len(d))
  n_anneal  <- result$n_anneal %||% 0L   # Phase 1 / Phase 2 boundary

  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, if (show_convergence) 4L else 1L), mar = c(5, 4, 3, 1))

  # Panel 1: stacked bar chart of soft weights
  ord      <- if (sort_by_weight) order(pw[,1], decreasing = TRUE) else seq_len(d)
  cols_bar <- c("#2166AC", "#D6604D")
  barplot(t(pw[ord,, drop=FALSE]), col=cols_bar, border=NA,
          names.arg=fn[ord], las=2, cex.names=max(0.4, 1-d/200),
          ylim=c(0,1), ylab="Soft weight", main=title)
  legend("topright", legend=colnames(pw), fill=cols_bar, border=NA, bty="n", cex=0.85)

  if (!show_convergence) return(invisible(result))

  .phase_line <- function(n_anneal, iters) {
    if (n_anneal > 0 && n_anneal < iters)
      abline(v = n_anneal + 0.5, lty = 3, lwd = 2, col = "grey40")
  }

  # Panel 2: weight certainty trace
  wh        <- result$weight_history
  certainty <- vapply(wh, function(w) mean(abs(w[,1] - 0.5)), numeric(1))
  plot(seq_along(wh), certainty, type="b", pch=19, col="#2166AC",
       xlab="Iteration", ylab="Mean |w_A - 0.5|",
       main="Weight certainty", ylim=c(0, 0.5))
  abline(h=0.25, lty=2, col="grey60")
  .phase_line(n_anneal, length(wh))

  # Panel 3: mean feature score difference trace
  sh       <- result$score_history
  sd_trace <- vapply(sh, function(s) mean(s[,1] - s[,2]), numeric(1))
  plot(seq_along(sh), sd_trace, type="b", pch=19, col="#D6604D",
       xlab="Iteration", ylab="Mean score(A) - score(B)",
       main="Mean feature score difference")
  abline(h=0, lty=2, col="grey60")
  .phase_line(n_anneal, length(sh))

  # Panel 4: full ELBO trace F = L_total + H(Gamma_1) + H(Gamma_2) + cell priors
  # This is the true monotone objective (CAVI); L_total alone is not monotone.
  ll <- result$ll_history
  plot(seq_along(ll), ll, type="b", pch=19, col="#1B7837",
       xlab="Iteration", ylab="Full ELBO (F)",
       main=expression(F ~ "= " ~ L[total] ~ "+ H(" * Gamma[1] * ") + H(" * Gamma[2] * ")"))
  .phase_line(n_anneal, length(ll))
  if (n_anneal > 0 && n_anneal < length(ll))
    legend("bottomright", legend=c("Annealing", "Exact EM (T=1)"),
           lty=c(1,1), col=c("grey40","#1B7837"), bty="n", cex=0.8)

  invisible(result)
}


# ============================================================
# 6.  Simulation + evaluation
# ============================================================

.simulate_dual_signal_block <- function(t,
                                        d,
                                        family = c("sinusoidal", "linear", "monotone", "quadratic"),
                                        sigma = 0.3,
                                        signal_range = c(0.5, 2),
                                        linear_slope_range = c(0.8, 2),
                                        monotone_power_range = c(0.7, 2.5),
                                        quadratic_curvature_range = c(0.8, 2),
                                        quadratic_center_range = c(0.25, 0.75),
                                        intercept_sd = 0,
                                        sinusoid_freq = 1:4) {
  family <- match.arg(family)
  n <- length(t)
  d <- as.integer(d)
  if (d <= 0L) return(matrix(0, nrow = n, ncol = 0L))

  amps <- stats::runif(d, min(signal_range), max(signal_range))
  intercepts <- stats::rnorm(d, mean = 0, sd = intercept_sd)
  out <- matrix(0, nrow = n, ncol = d)

  for (j in seq_len(d)) {
    base <- switch(
      family,
      sinusoidal = {
        freq_j <- sample(sinusoid_freq, 1L)
        phase_j <- stats::runif(1L, 0, 2 * pi)
        sin(freq_j * pi * t + phase_j)
      },
      linear = {
        slope_j <- stats::runif(1L, min(linear_slope_range), max(linear_slope_range))
        slope_j <- slope_j * sample(c(-1, 1), 1L)
        slope_j * (t - 0.5)
      },
      monotone = {
        power_left <- stats::runif(1L, min(monotone_power_range), max(monotone_power_range))
        power_right <- stats::runif(1L, min(monotone_power_range), max(monotone_power_range))
        logistic_mid <- stats::runif(1L, 0.2, 0.8)
        logistic_steep <- stats::runif(1L, 3, 12)
        exp_rate <- stats::runif(1L, 1, 6)
        hill_left <- stats::runif(1L, min(monotone_power_range), max(monotone_power_range))
        hill_right <- stats::runif(1L, min(monotone_power_range), max(monotone_power_range))

        logistic_basis <- stats::plogis((t - logistic_mid) * logistic_steep)
        logistic_basis <- (logistic_basis - min(logistic_basis)) /
          (max(logistic_basis) - min(logistic_basis) + 1e-12)

        exp_basis <- (exp(exp_rate * t) - 1) / (exp(exp_rate) - 1)

        hill_num <- t^hill_left
        hill_den <- hill_num + (1 - t)^hill_right + 1e-12
        hill_basis <- hill_num / hill_den

        monotone_basis <- cbind(
          t^power_left,
          1 - (1 - t)^power_right,
          logistic_basis,
          exp_basis,
          hill_basis
        )
        mix_w <- stats::rexp(ncol(monotone_basis))
        mix_w <- mix_w / sum(mix_w)
        raw <- as.numeric(monotone_basis %*% mix_w)
        if (sample(c(FALSE, TRUE), 1L)) raw <- 1 - raw
        raw - mean(raw)
      },
      quadratic = {
        center_j <- stats::runif(1L, min(quadratic_center_range), max(quadratic_center_range))
        curv_j <- stats::runif(1L, min(quadratic_curvature_range), max(quadratic_curvature_range))
        curv_j <- curv_j * sample(c(-1, 1), 1L)
        raw <- curv_j * (t - center_j)^2
        raw - mean(raw)
      }
    )

    out[, j] <- intercepts[j] + amps[j] * base + stats::rnorm(n, 0, sigma)
  }

  out
}

#' Simulate a multi-ordering intrinsic-trajectory dataset
#'
#' @description
#' Simulates a feature-partition dataset with an arbitrary number of intrinsic
#' orderings. Each ordering contributes its own signal block, and optional
#' noise features are pure Gaussian noise.
#'
#' @param n Integer number of samples.
#' @param d_signal Integer vector giving the number of signal features assigned
#'   to each intrinsic ordering.
#' @param d_noise Integer number of pure-noise features.
#' @param sigma Observation noise standard deviation added independently to each
#'   simulated feature trajectory.
#' @param seed Optional integer seed.
#' @param trajectory_family Character scalar or length-\code{M} vector
#'   specifying the feature family for each ordering block. Each entry must be
#'   one of \code{"sinusoidal"}, \code{"linear"}, \code{"monotone"}, or
#'   \code{"quadratic"}.
#' @param signal_range Length-2 positive range controlling the per-feature
#'   signal amplitude.
#' @param linear_slope_range Length-2 positive range for linear slopes.
#' @param monotone_power_range Length-2 positive range controlling the power-like
#'   basis shapes used in the monotone family. The monotone generator mixes
#'   several increasing bases internally to create more varied monotone curves.
#' @param quadratic_curvature_range Length-2 positive range for quadratic
#'   curvature magnitudes.
#' @param quadratic_center_range Length-2 range inside \code{[0, 1]} for the
#'   quadratic vertex.
#' @param intercept_sd Standard deviation of feature-specific intercept shifts.
#'   Defaults to \code{0} to preserve the historical simulation baseline.
#' @param sinusoid_freq Integer vector of admissible sinusoid frequencies.
#' @param latent_positions Optional \code{n x M} matrix of true latent sample
#'   positions. If \code{NULL}, each intrinsic ordering is generated
#'   independently from \code{Uniform(0, 1)}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{X}: shuffled simulated data matrix
#'   \item \code{true_assign}: feature labels such as \code{"A"}, \code{"B"},
#'     ..., or \code{"noise"}
#'   \item \code{latent_positions}: \code{n x M} matrix of true latent sample
#'     locations
#'   \item \code{ordering_labels}: ordering labels used in \code{true_assign}
#'   \item \code{original_order}: column permutation applied to the feature blocks
#'   \item \code{trajectory_family}: the realized family labels for each ordering
#' }
#' @export
simulate_intrinsic_trajectories <- function(n = 200,
                                            d_signal = c(30, 30),
                                            d_noise = 20,
                                            sigma = 0.3,
                                            seed = 42L,
                                            trajectory_family = "sinusoidal",
                                            signal_range = c(0.5, 2),
                                            linear_slope_range = c(0.8, 2),
                                            monotone_power_range = c(0.7, 2.5),
                                            quadratic_curvature_range = c(0.8, 2),
                                            quadratic_center_range = c(0.25, 0.75),
                                            intercept_sd = 0,
                                            sinusoid_freq = 1:4,
                                            latent_positions = NULL) {
  n <- as.integer(n)
  d_signal <- as.integer(d_signal)
  d_noise <- as.integer(d_noise)
  if (length(n) != 1L || is.na(n) || n < 2L) stop("n must be a single integer >= 2.")
  if (!length(d_signal) || any(is.na(d_signal)) || any(d_signal < 0L)) {
    stop("d_signal must be a nonnegative integer vector of positive length.")
  }
  if (length(d_noise) != 1L || is.na(d_noise) || d_noise < 0L) {
    stop("d_noise must be a single nonnegative integer.")
  }
  if ((sum(d_signal) + d_noise) < 1L) {
    stop("At least one of the signal or noise feature counts must be positive.")
  }
  if (!is.finite(sigma) || sigma < 0) stop("sigma must be a nonnegative number.")
  if (!is.null(seed)) set.seed(seed)

  M <- length(d_signal)
  trajectory_family <- rep(trajectory_family, length.out = M)
  valid_family <- c("sinusoidal", "linear", "monotone", "quadratic")
  if (!all(trajectory_family %in% valid_family)) {
    stop("trajectory_family must contain only: ",
         paste(valid_family, collapse = ", "), ".")
  }
  if (!is.numeric(signal_range) || length(signal_range) != 2L ||
      any(!is.finite(signal_range)) || any(signal_range <= 0)) {
    stop("signal_range must be a positive numeric vector of length 2.")
  }
  if (!is.numeric(linear_slope_range) || length(linear_slope_range) != 2L ||
      any(!is.finite(linear_slope_range)) || any(linear_slope_range <= 0)) {
    stop("linear_slope_range must be a positive numeric vector of length 2.")
  }
  if (!is.numeric(monotone_power_range) || length(monotone_power_range) != 2L ||
      any(!is.finite(monotone_power_range)) || any(monotone_power_range <= 0)) {
    stop("monotone_power_range must be a positive numeric vector of length 2.")
  }
  if (!is.numeric(quadratic_curvature_range) || length(quadratic_curvature_range) != 2L ||
      any(!is.finite(quadratic_curvature_range)) || any(quadratic_curvature_range <= 0)) {
    stop("quadratic_curvature_range must be a positive numeric vector of length 2.")
  }
  if (!is.numeric(quadratic_center_range) || length(quadratic_center_range) != 2L ||
      any(!is.finite(quadratic_center_range)) ||
      quadratic_center_range[1] < 0 || quadratic_center_range[2] > 1 ||
      quadratic_center_range[1] > quadratic_center_range[2]) {
    stop("quadratic_center_range must be a finite numeric vector of length 2 inside [0, 1].")
  }
  if (!is.finite(intercept_sd) || intercept_sd < 0) {
    stop("intercept_sd must be nonnegative.")
  }
  if (!length(sinusoid_freq) || any(!is.finite(sinusoid_freq))) {
    stop("sinusoid_freq must contain at least one finite value.")
  }

  if (is.null(latent_positions)) {
    latent_positions <- matrix(stats::runif(n * M), nrow = n, ncol = M)
  } else {
    latent_positions <- as.matrix(latent_positions)
    if (!all(dim(latent_positions) == c(n, M))) {
      stop("latent_positions must be an n x length(d_signal) matrix.")
    }
    if (any(!is.finite(latent_positions)) ||
        any(latent_positions < 0) || any(latent_positions > 1)) {
      stop("latent_positions must contain finite values inside [0, 1].")
    }
  }

  ordering_labels <- if (M <= 26L) LETTERS[seq_len(M)] else paste0("ord", seq_len(M))
  signal_blocks <- lapply(seq_len(M), function(m) {
    .simulate_dual_signal_block(
      t = latent_positions[, m],
      d = d_signal[m],
      family = trajectory_family[m],
      sigma = sigma,
      signal_range = signal_range,
      linear_slope_range = linear_slope_range,
      monotone_power_range = monotone_power_range,
      quadratic_curvature_range = quadratic_curvature_range,
      quadratic_center_range = quadratic_center_range,
      intercept_sd = intercept_sd,
      sinusoid_freq = sinusoid_freq
    )
  })
  X_noise <- if (d_noise > 0L) {
    matrix(stats::rnorm(n * d_noise), nrow = n, ncol = d_noise)
  } else {
    matrix(0, nrow = n, ncol = 0L)
  }

  X_full <- do.call(cbind, c(signal_blocks, list(X_noise)))
  assign_full <- unlist(Map(function(lbl, d_m) rep(lbl, d_m), ordering_labels, d_signal))
  assign_full <- c(assign_full, rep("noise", d_noise))
  perm <- sample.int(sum(d_signal) + d_noise)
  colnames(X_full) <- paste0("V", seq_len(ncol(X_full)))

  list(
    X = X_full[, perm, drop = FALSE],
    true_assign = assign_full[perm],
    latent_positions = latent_positions,
    ordering_labels = ordering_labels,
    original_order = perm,
    trajectory_family = trajectory_family
  )
}

#' Simulate a dual-trajectory dataset with configurable trajectory families
#'
#' @description
#' Backward-compatible convenience wrapper around
#' \code{\link{simulate_intrinsic_trajectories}} for the two-ordering case.
#' Preserves the old \code{t1}/\code{t2} return values and the
#' \code{crossing} shortcut.
#'
#' @inheritParams simulate_intrinsic_trajectories
#' @param d1,d2 Integer numbers of signal features assigned to trajectories A
#'   and B.
#' @param crossing Logical; if \code{TRUE}, the second latent ordering is
#'   generated as a noisy reversed version of the first, making the two
#'   orderings harder to disentangle.
#'
#' @return A list with the same fields as
#'   \code{\link{simulate_intrinsic_trajectories}}, plus backward-compatible
#'   aliases \code{t1} and \code{t2}.
#' @export
simulate_dual_trajectory <- function(n = 200,
                                     d1 = 30,
                                     d2 = 30,
                                     d_noise = 20,
                                     sigma = 0.3,
                                     crossing = FALSE,
                                     seed = 42L,
                                     trajectory_family = c("sinusoidal", "sinusoidal"),
                                     signal_range = c(0.5, 2),
                                     linear_slope_range = c(0.8, 2),
                                     monotone_power_range = c(0.7, 2.5),
                                     quadratic_curvature_range = c(0.8, 2),
                                     quadratic_center_range = c(0.25, 0.75),
                                     intercept_sd = 0,
                                     sinusoid_freq = 1:4) {
  n <- as.integer(n)
  d1 <- as.integer(d1)
  d2 <- as.integer(d2)
  d_noise <- as.integer(d_noise)
  if (length(n) != 1L || is.na(n) || n < 2L) stop("n must be a single integer >= 2.")
  if (any(c(d1, d2, d_noise) < 0L) || any(is.na(c(d1, d2, d_noise)))) {
    stop("d1, d2, and d_noise must be nonnegative integers.")
  }
  if ((d1 + d2 + d_noise) < 1L) {
    stop("At least one of d1, d2, or d_noise must be positive.")
  }
  if (!is.null(seed)) set.seed(seed)

  t1 <- stats::runif(n)
  if (crossing) {
    t2_raw <- 1 - t1 + stats::rnorm(n, 0, 0.3)
    t2 <- (t2_raw - min(t2_raw)) / (max(t2_raw) - min(t2_raw))
  } else {
    t2 <- stats::runif(n)
  }

  sim <- simulate_intrinsic_trajectories(
    n = n,
    d_signal = c(d1, d2),
    d_noise = d_noise,
    sigma = sigma,
    seed = NULL,
    trajectory_family = trajectory_family,
    signal_range = signal_range,
    linear_slope_range = linear_slope_range,
    monotone_power_range = monotone_power_range,
    quadratic_curvature_range = quadratic_curvature_range,
    quadratic_center_range = quadratic_center_range,
    intercept_sd = intercept_sd,
    sinusoid_freq = sinusoid_freq,
    latent_positions = cbind(t1, t2)
  )

  sim$t1 <- t1
  sim$t2 <- t2
  sim
}

#' Evaluate soft partition accuracy (handles label-switching)
#' @param result      Output of soft_two_trajectory_EM().
#' @param true_assign Ground-truth vector ("A","B","noise").
#' @return List: accuracy, best_alignment, confusion_table.
#' @noRd
evaluate_partition <- function(result, true_assign) {
  pred <- result$assign
  true_assign <- as.character(true_assign)
  pred <- as.character(pred)
  if (length(pred) != length(true_assign)) {
    stop("length(result$assign) must match length(true_assign).")
  }

  mask <- true_assign %in% c("A","B")
  pred_sub <- pred[mask]
  true_sub <- true_assign[mask]
  if (!length(true_sub)) {
    return(list(
      accuracy = NA_real_,
      best_alignment = "no A/B labels in truth",
      confusion_table = table(predicted = pred_sub, truth = true_sub)
    ))
  }

  acc1 <- mean(pred_sub == true_sub)
  flip <- pred_sub
  flip[pred_sub == "A"] <- "B"
  flip[pred_sub == "B"] <- "A"
  acc2 <- mean(flip == true_sub)
  if (acc1 >= acc2)
    list(accuracy=round(acc1,4), best_alignment="direct (A=A, B=B)",
         confusion_table=table(predicted=pred_sub, truth=true_sub))
  else
    list(accuracy=round(acc2,4), best_alignment="flipped (A=B, B=A)",
         confusion_table=table(predicted=flip, truth=true_sub))
}
