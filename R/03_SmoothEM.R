# =========================
# Convergence / objective utilities
# =========================

#' Internal objective utilities for SmoothEM
#'
#' @description
#' Internal helper functions used by \code{\link{EM_algorithm}} to evaluate the
#' penalized observed-data objective and the penalized ELBO (Q-function + entropy),
#' optionally including a Gaussian prior penalty on the stacked mean vector.
#'
#' @param X Numeric matrix \code{(n x d)} of observations.
#' @param Gamma Numeric matrix \code{(n x K)} of responsibilities (rows sum to 1).
#' @param params List with fields \code{pi}, \code{mu}, \code{sigma} (and possibly cached fields).
#' @param Q_prior Optional precision matrix for a Gaussian prior on stacked \code{mu}.
#' @param eigen_tol Optional tolerance used by \code{generalized_logdet()}.
#' @param rank_deficiency Integer rank deficiency for the pseudo-determinant correction.
#'
#' @return
#' \itemize{
#'   \item \code{compute_log_joint_observed}: scalar penalized observed-data log-likelihood.
#'   \item \code{compute_penalized_ELBO}: scalar penalized ELBO.
#' }
#'
#' @keywords internal
#' @rdname EM_algorithm
compute_log_joint_observed <- function(X, params, Q_prior = NULL, eigen_tol = NULL, rank_deficiency = 0) {
  K <- length(params$pi)

  log_dens_mat <- sapply(seq_len(K), function(k) {
    mvtnorm::dmvnorm(
      X,
      mean  = params$mu[[k]],
      sigma = params$sigma[[k]],
      log   = TRUE
    )
  })

  log_gamma <- sweep(log_dens_mat, 2, log(params$pi), "+")  # n × K
  loglik <- sum(matrixStats::rowLogSumExps(log_gamma))

  if (!is.null(Q_prior)) {
    U_vec   <- unlist(params$mu, use.names = FALSE)
    penalty <- 0.5 * as.numeric(crossprod(U_vec, as.matrix(Q_prior %*% U_vec)))

    logdetQ <- generalized_logdet(Q = Q_prior, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)

    loglik <- loglik - penalty + 0.5 * logdetQ
  }

  loglik
}


#' @keywords internal
#' @rdname EM_algorithm
compute_penalized_ELBO <- function(X, Gamma, params, Q_prior = NULL, eigen_tol = NULL, rank_deficiency = 0) {
  K <- length(params$pi)

  log_dens_mat <- sapply(seq_len(K), function(k) {
    mvtnorm::dmvnorm(
      X,
      mean  = params$mu[[k]],
      sigma = params$sigma[[k]],
      log   = TRUE
    )
  })
  log_gamma <- sweep(log_dens_mat, 2, log(params$pi), "+")  # n × K

  Qval <- sum(Gamma * log_gamma)
  entropy <- -sum(Gamma * log(Gamma + 1e-300))
  elbo <- Qval + entropy

  if (!is.null(Q_prior)) {
    U_vec   <- unlist(params$mu, use.names = FALSE)
    penalty <- 0.5 * as.numeric(crossprod(U_vec, as.matrix(Q_prior %*% U_vec)))

    logdetQ <- generalized_logdet(Q = Q_prior, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)

    elbo <- elbo - penalty + 0.5 * logdetQ
  }

  elbo
}


# =========================
# E-step
# =========================

#' Internal E-step for SmoothEM
#'
#' @description
#' Computes responsibilities \eqn{\gamma_{ik}} given current parameters.
#' Expects \code{params$invSigma} and \code{params$logdet} to be precomputed
#' (see \code{init_cov_cache_fast()}).
#'
#' @param data Numeric matrix \code{(n x d)}.
#' @param params Parameter list including cached \code{invSigma} and \code{logdet}.
#'
#' @return Numeric matrix \code{(n x K)} of responsibilities.
#'
#' @keywords internal
#' @rdname EM_algorithm
ESTEP <- function(data, params) {
  n <- nrow(data)
  K <- length(params$pi)
  d <- ncol(data)

  log_gamma <- matrix(0, n, K)
  Mu_mat <- do.call(cbind, params$mu)  # d x K

  for (k in seq_len(K)) {
    Xc <- sweep(data, 2, Mu_mat[, k], "-")
    Qf <- rowSums((Xc %*% params$invSigma[[k]]) * Xc)
    log_gamma[, k] <- log(params$pi[k]) - 0.5 * (Qf + params$logdet[k] + d * log(2 * pi))
  }

  lse <- matrixStats::rowLogSumExps(log_gamma)
  exp(log_gamma - lse)
}


# =========================
# M-step
# =========================

#' Internal M-step for SmoothEM
#'
#' @description
#' Updates mixture parameters given responsibilities, optionally with a quadratic
#' prior penalty on stacked means via \code{Q_prior}. Supports covariance models:
#' \code{"VVV"}, \code{"VII"}, \code{"EII"}, \code{"EEI"}.
#'
#' @param data Numeric matrix \code{(n x d)}.
#' @param gamma Numeric matrix \code{(n x K)} responsibilities.
#' @param params Current parameter list (used mainly for starting \code{sigma}).
#' @param Q_prior Optional precision matrix on stacked \code{mu}.
#' @param relative_lambda Logical; if TRUE, rescales \code{Q_prior} by current marginal variances
#'   (requires identical covariances across clusters; typically \code{"EEI"} usage).
#' @param modelName Covariance model: one of \code{"VVV"}, \code{"VII"}, \code{"EII"}, \code{"EEI"}.
#' @param iterate_once Logical; if TRUE, do one update pass (no inner alternation).
#' @param nugget Nonnegative diagonal jitter added to covariance estimates.
#' @param rank_deficiency Rank deficiency used in your EEI + relative-lambda variance update.
#' @param tol_inner,max_inner Inner-loop controls when iterating \code{mu <-> sigma}.
#' @param verbose Logical.
#'
#' @return List with updated \code{pi}, \code{mu}, \code{sigma}.
#'
#' @keywords internal
#' @rdname EM_algorithm
MSTEP <- function(
    data, gamma, params, Q_prior = NULL,
    relative_lambda = FALSE,
    modelName    = "VVV",
    iterate_once = TRUE,
    nugget      = 0,
    rank_deficiency = 0,
    tol_inner    = 1e-6,
    max_inner    = 20,
    verbose      = FALSE
) {
  n <- nrow(data)
  d <- ncol(data)
  K <- ncol(gamma)
  DK <- d * K

  modelName <- match.arg(modelName, c("VVV", "VII", "EII", "EEI"))

  # 1) Update mixing proportions
  Nk <- colSums(gamma)
  Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / n

  # 2) Weighted sums
  Wx <- t(data) %*% gamma  # d x K

  # 3) Init mu and sigma
  mu_list <- lapply(seq_len(K), function(k) as.numeric(Wx[, k] / Nk[k]))
  sigma_list <- params$sigma

  # Keep original Q prior for relative scaling
  Q_prior_orig <- Q_prior

  # 4) Relative lambda scaling requires shared diagonal covariance (EEI in your usage)
  if (relative_lambda && !is.null(Q_prior)) {
    same <- all(vapply(sigma_list, function(S) isTRUE(all.equal(S, sigma_list[[1]])), logical(1)))
    if (!same) stop("relative_lambda requires identical covariances across clusters (typically EEI).")

    sigma_vec <- diag(sigma_list[[1]])
    scale_vec <- rep(1 / sqrt(pmax(sigma_vec, 1e-12)), times = K)
    Sscale <- Matrix::Diagonal(x = scale_vec)
    Q_prior <- Sscale %*% Q_prior_orig %*% Sscale
  }

  # ---- helper: penalized mu update given invSigma_list ----
  update_mu_penalized <- function(invSigma_list, Q_prior_use) {
    blocks <- vector("list", K)
    rhs <- numeric(DK)

    for (k in seq_len(K)) {
      Sinv <- invSigma_list[[k]]
      blocks[[k]] <- Nk[k] * Sinv
      rhs[((k - 1) * d + 1):(k * d)] <- Sinv %*% Wx[, k]
    }

    Hmat <- Matrix::bdiag(blocks) + Q_prior_use

    if (DK > 200) {
      fac <- Matrix::Cholesky(Hmat, LDL = FALSE)
      mu_vec <- Matrix::solve(fac, rhs)
    } else {
      mu_vec <- solve(as.matrix(Hmat), rhs)
    }

    mu_vec <- as.numeric(mu_vec)
    split(mu_vec, rep(seq_len(K), each = d))
  }

  # ---- helper: update sigma given mu_list ----
  update_sigma <- function(mu_list_now) {
    if (modelName == "EEI") {
      GX  <- t(data) %*% gamma
      GX2 <- t(data^2) %*% gamma
      mu_mat <- do.call(cbind, mu_list_now)  # d x K

      V <- rowSums(GX2 - 2 * (mu_mat * GX) + sweep(mu_mat^2, 2, Nk, "*"))

      if (relative_lambda && !is.null(Q_prior_orig)) {
        U_vec <- unlist(mu_list_now, use.names = FALSE)
        U_vec_d <- numeric(d)
        for (dd in seq_len(d)) {
          u_index <- seq(from = dd, to = length(U_vec), by = d)
          U_vec_d[dd] <- as.numeric(crossprod(U_vec[u_index], Q_prior_orig[u_index, u_index] %*% U_vec[u_index]))
        }
        denom <- (n + K - rank_deficiency)
        Vs <- (V + U_vec_d) / denom
        sharedSig <- diag(Vs, d) + nugget * diag(d)
      } else {
        Vs <- V / n
        sharedSig <- diag(Vs, d) + nugget * diag(d)
      }

      replicate(K, sharedSig, simplify = FALSE)

    } else if (modelName == "VII") {
      lapply(seq_len(K), function(k) {
        diff <- sweep(data, 2, mu_list_now[[k]], "-")
        lam <- sum((diff^2) * gamma[, k]) / (Nk[k] * d)
        diag(lam, d) + nugget * diag(d)
      })

    } else if (modelName == "VVV") {
      lapply(seq_len(K), function(k) {
        diff <- sweep(data, 2, mu_list_now[[k]], "-")
        W <- sqrt(gamma[, k])
        crossprod(diff * W, diff * W) / Nk[k] + nugget * diag(d)
      })

    } else { # EII (spherical shared across clusters)
      sse_k <- vapply(seq_len(K), function(k) {
        diff <- sweep(data, 2, mu_list_now[[k]], "-")
        sum((diff^2) * gamma[, k])
      }, numeric(1))

      total <- sum(sse_k)
      lam_sh <- total / (sum(Nk) * d)
      replicate(K, diag(lam_sh, d) + nugget * diag(d), simplify = FALSE)
    }
  }

  # ---- 5) single-shot branch ----
  if (iterate_once || is.null(Q_prior)) {
    # 5a) mu update (penalized only if Q_prior provided)
    if (!is.null(Q_prior)) {
      invSigma_list <- params$invSigma
      mu_list <- update_mu_penalized(invSigma_list, Q_prior)
    }

    # 5b) sigma update once
    sigma_list <- update_sigma(mu_list)

    return(list(pi = pi_new, mu = mu_list, sigma = sigma_list))
  }

  # ---- 6) inner coordinate updates (mu <-> sigma) ----
  for (inner in seq_len(max_inner)) {
    mu_old_vec <- unlist(mu_list, use.names = FALSE)
    sigma_old_vec <- unlist(lapply(sigma_list, function(S) as.numeric(S)), use.names = FALSE)

    # 6a) recompute invSigma from current sigma_list
    invSigma_list <- vector("list", K)
    for (k in seq_len(K)) {
      L <- chol(sigma_list[[k]])
      invSigma_list[[k]] <- chol2inv(L)
    }

    # 6b) if relative_lambda, rescale Q_prior using current sigma (EEI)
    Q_prior_use <- Q_prior
    if (relative_lambda && !is.null(Q_prior_orig)) {
      sigma_vec <- diag(sigma_list[[1]])
      scale_vec <- rep(1 / sqrt(pmax(sigma_vec, 1e-12)), times = K)
      Sscale <- Matrix::Diagonal(x = scale_vec)
      Q_prior_use <- Sscale %*% Q_prior_orig %*% Sscale
    }

    # 6c) mu update
    mu_list <- update_mu_penalized(invSigma_list, Q_prior_use)

    # 6d) sigma update
    sigma_list <- update_sigma(mu_list)

    # 6e) convergence check
    mu_new_vec <- unlist(mu_list, use.names = FALSE)
    sigma_new_vec <- unlist(lapply(sigma_list, function(S) as.numeric(S)), use.names = FALSE)

    mu_diff <- max(abs(mu_new_vec - mu_old_vec))
    sigma_diff <- max(abs(sigma_new_vec - sigma_old_vec))

    if (verbose) {
      cat(sprintf("  MSTEP inner iter %d: Δμ=%.3e, ΔΣ=%.3e\n", inner, mu_diff, sigma_diff))
    }
    if (max(mu_diff, sigma_diff) < tol_inner) break
  }

  list(pi = pi_new, mu = mu_list, sigma = sigma_list)
}



# =========================
# EM wrapper (exported)
# =========================

#' Penalized Expectation–Maximization (EM) algorithm for SmoothEM
#'
#' @description
#' Fits a Gaussian mixture model with optional quadratic prior penalty on the stacked
#' mean vector \code{mu}. The algorithm alternates E-steps and penalized M-steps and
#' records both a penalized observed-data objective and a penalized ELBO trace.
#'
#' @param data Numeric matrix \code{(n x d)} of observations.
#' @param init_params Initial parameters as a list with \code{pi}, \code{mu}, \code{sigma}.
#' @param Q_prior Optional precision matrix on stacked means (\code{d*K x d*K}).
#' @param iterate_once Logical; if TRUE, stop after one EM iteration (often used inside wrappers).
#' @param max_inner Maximum number of inner iterations inside the M-step.
#' @param modelName Covariance model: one of \code{"VVV"}, \code{"VII"}, \code{"EII"}, \code{"EEI"}.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance on ELBO changes.
#' @param inner_tol Tolerance for inner iterations inside the M-step.
#' @param eigen_tol Optional tolerance passed to \code{generalized_logdet()}.
#' @param rank_deficiency Rank deficiency for pseudo-determinant correction (commonly \code{q*d}).
#' @param nugget Nonnegative diagonal jitter added to covariance estimates.
#' @param relative_lambda Logical; if TRUE, rescales \code{Q_prior} by current marginal variances.
#' @param verbose Logical; print progress.
#' @param include.data Logical; if TRUE, include \code{data} in the returned list.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{params}: final parameters (including cached \code{invSigma} and \code{logdet}).
#'   \item \code{gamma}: responsibility matrix \code{(n x K)}.
#'   \item \code{elbo_trace}: numeric vector.
#'   \item \code{loglik_trace}: numeric vector.
#' }
#'
#' @export
EM_algorithm <- function(
    data, init_params,
    Q_prior   = NULL,
    iterate_once = TRUE,
    max_inner = 10,
    modelName = "VVV",
    max_iter  = 100,
    tol       = 1e-4,
    inner_tol = 1e-6,
    eigen_tol = NULL,
    rank_deficiency = 0,
    nugget    = 0,
    relative_lambda = FALSE,
    verbose   = TRUE,
    include.data = TRUE
) {
  data <- as.matrix(data)

  params <- init_cov_cache_fast(init_params)
  K <- length(params$pi)

  loglik_trace <- numeric(max_iter)
  elbo_trace   <- numeric(max_iter)

  Q_prior_orig <- Q_prior

  for (iter in seq_len(max_iter)) {

    # E-step
    gamma <- ESTEP(data, params)

    # M-step
    new_params <- MSTEP(
      data = data,
      gamma = gamma,
      params = params,
      Q_prior = Q_prior_orig,
      relative_lambda = relative_lambda,
      modelName = modelName,
      iterate_once = iterate_once,
      nugget = nugget,
      rank_deficiency = rank_deficiency,
      tol_inner = inner_tol,
      max_inner = max_inner,
      verbose = verbose
    )

    last_params <- params
    params <- init_cov_cache_fast(new_params)

    # For objective evaluation, use the *current* scaled Q_prior (if relative_lambda)
    Q_eval <- Q_prior_orig
    if (relative_lambda && !is.null(Q_prior_orig)) {
      # requires EEI-like shared diagonal
      same <- all(vapply(params$sigma, function(S) isTRUE(all.equal(S, params$sigma[[1]])), logical(1)))
      if (!same) stop("relative_lambda requires identical covariances across clusters (typically EEI).")

      sigma_vec <- diag(params$sigma[[1]])
      scale_vec <- rep(1 / sqrt(pmax(sigma_vec, 1e-12)), times = K)
      Sscale <- Matrix::Diagonal(x = scale_vec)
      Q_eval <- Sscale %*% Q_prior_orig %*% Sscale
    }

    # Record log-likelihood and ELBO
    ll   <- compute_log_joint_observed(data, params, Q_eval, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)
    elbo <- compute_penalized_ELBO(data, gamma, params, Q_eval, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)

    loglik_trace[iter] <- ll
    elbo_trace[iter]   <- elbo

    if (verbose) {
      cat(sprintf("Iteration %3d: penLogLik = %.6f, ELBO = %.6f\n", iter, ll, elbo))
    }

    if (iter > 1 && abs(elbo - elbo_trace[iter - 1]) < tol) {
      # If ELBO decreased materially, rollback to last params and stop
      if (elbo_trace[iter - 1] - elbo > 1e-10) {
        if (verbose) {
          cat(sprintf("ELBO decreased at iteration %d: %.6f -> %.6f. Rolling back.\n",
                      iter, elbo_trace[iter - 1], elbo))
        }
        params <- last_params
        gamma  <- ESTEP(data, params)

        loglik_trace <- loglik_trace[1:(iter - 1)]
        elbo_trace   <- elbo_trace[1:(iter - 1)]
        break
      } else {
        if (verbose) {
          cat(sprintf("Converged at iteration %d with ELBO %.6f\n", iter, elbo))
          cat(sprintf("Final penLogLik: %.6f\n", ll))
        }
        loglik_trace <- loglik_trace[1:iter]
        elbo_trace   <- elbo_trace[1:iter]
        break
      }
    }
  }

  list(
    params       = params,
    gamma        = gamma,
    elbo_trace   = elbo_trace,
    loglik_trace = loglik_trace,
    data         = if (include.data) data else NULL
  )
}


