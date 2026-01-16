# =========================
# Small internal helpers
# =========================

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Fast log-density matrix using cached invSigma/logdet (no mvtnorm::dmvnorm)
.log_dens_mat_fast <- function(X, params) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  K <- length(params$pi)

  if (is.null(params$invSigma) || is.null(params$logdet)) {
    stop("params must contain cached fields $invSigma and $logdet.")
  }

  Mu_mat <- do.call(cbind, params$mu)  # d x K
  out <- matrix(0, n, K)

  const <- d * log(2 * pi)

  for (k in seq_len(K)) {
    Xc <- sweep(X, 2, Mu_mat[, k], "-")
    Qf <- rowSums((Xc %*% params$invSigma[[k]]) * Xc)
    out[, k] <- -0.5 * (Qf + params$logdet[k] + const)
  }
  out
}

# Apply relative-lambda scaling to a BASE precision (no lambda factor), if needed.
# Returns Q_base_eval (still base; lambda not multiplied).
.compute_Q_base_eval <- function(params, Q_base, relative_lambda) {
  if (is.null(Q_base)) return(NULL)
  Q_eval <- Q_base

  if (isTRUE(relative_lambda)) {
    # requires EEI-like shared diagonal covariance
    same <- all(vapply(params$sigma, function(S) isTRUE(all.equal(S, params$sigma[[1]])), logical(1)))
    if (!same) stop("relative_lambda requires identical covariances across clusters (typically EEI).")

    K <- length(params$pi)
    sigma_vec <- diag(params$sigma[[1]])
    scale_vec <- rep(1 / sqrt(pmax(sigma_vec, 1e-12)), times = K)
    Sscale <- Matrix::Diagonal(x = scale_vec)

    Q_eval <- Sscale %*% Q_eval %*% Sscale
  }

  Q_eval
}

.stack_mu_vec <- function(params) {
  if (is.list(params$mu)) {
    unlist(lapply(params$mu, function(m) as.numeric(m)), use.names = FALSE)
  } else {
    as.numeric(params$mu)
  }
}

# Profile-style update for lambda given params and Q_base (base precision; not multiplied by lambda)
.estimate_lambda_star <- function(params, Q_base, relative_lambda, rank_deficiency,
                                  eps_quad = 1e-12, lambda_min = 1e-8, lambda_max = 1e8) {
  Qb <- .compute_Q_base_eval(params, Q_base, relative_lambda)
  if (is.null(Qb)) return(NA_real_)

  U <- .stack_mu_vec(params)
  quad <- as.numeric(crossprod(U, Qb %*% U))

  DK <- ncol(Qb)
  r  <- as.integer(DK - (rank_deficiency %||% 0L))
  r  <- max(r, 1L)

  lam <- r / pmax(quad, eps_quad)
  lam <- min(max(lam, lambda_min), lambda_max)
  lam
}


# =========================
# Convergence / objective utilities (internal)
# =========================

#' @keywords internal
#' @rdname EM_algorithm
compute_log_joint_observed <- function(X, params,
                                       Q_prior = NULL,
                                       eigen_tol = NULL,
                                       rank_deficiency = 0,
                                       # optional factorized mode:
                                       Q_base = NULL,
                                       lambda = NULL) {
  X <- as.matrix(X)
  K <- length(params$pi)

  # fast log densities
  log_dens_mat <- .log_dens_mat_fast(X, params)
  log_gamma <- sweep(log_dens_mat, 2, log(params$pi), "+")
  loglik <- sum(matrixStats::rowLogSumExps(log_gamma))

  # ---- prior contribution ----
  if (!is.null(Q_prior)) {
    # legacy/full-Q mode (Q_prior already includes lambda)
    U <- .stack_mu_vec(params)
    penalty <- 0.5 * as.numeric(crossprod(U, Q_prior %*% U))
    logdetQ <- generalized_logdet(Q = Q_prior, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)
    return(loglik - penalty + 0.5 * logdetQ)
  }

  if (!is.null(Q_base)) {
    if (is.null(lambda) || !is.finite(lambda) || lambda <= 0) {
      stop("In factorized mode, lambda must be a positive finite number.")
    }
    U <- .stack_mu_vec(params)

    quad  <- as.numeric(crossprod(U, Q_base %*% U))
    logdetQb <- generalized_logdet(Q = Q_base, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)

    DK <- ncol(Q_base)
    r  <- as.integer(DK - (rank_deficiency %||% 0L))
    r  <- max(r, 1L)

    # logdet(λQb)=logdet(Qb)+r log λ
    return(loglik - 0.5 * lambda * quad + 0.5 * (logdetQb + r * log(lambda)))
  }

  loglik
}

#' @keywords internal
#' @rdname EM_algorithm
compute_penalized_ELBO <- function(X, Gamma, params,
                                   Q_prior = NULL,
                                   eigen_tol = NULL,
                                   rank_deficiency = 0,
                                   # optional factorized mode:
                                   Q_base = NULL,
                                   lambda = NULL) {
  X <- as.matrix(X)
  Gamma <- as.matrix(Gamma)
  K <- length(params$pi)

  log_dens_mat <- .log_dens_mat_fast(X, params)
  log_gamma <- sweep(log_dens_mat, 2, log(params$pi), "+")  # n x K

  Qval <- sum(Gamma * log_gamma)
  entropy <- -sum(Gamma * log(Gamma + 1e-300))
  elbo <- Qval + entropy

  # ---- prior contribution ----
  if (!is.null(Q_prior)) {
    U <- .stack_mu_vec(params)
    penalty <- 0.5 * as.numeric(crossprod(U, Q_prior %*% U))
    logdetQ <- generalized_logdet(Q = Q_prior, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)
    return(elbo - penalty + 0.5 * logdetQ)
  }

  if (!is.null(Q_base)) {
    if (is.null(lambda) || !is.finite(lambda) || lambda <= 0) {
      stop("In factorized mode, lambda must be a positive finite number.")
    }
    U <- .stack_mu_vec(params)

    quad  <- as.numeric(crossprod(U, Q_base %*% U))
    logdetQb <- generalized_logdet(Q = Q_base, eigen_tol = eigen_tol, rank_deficiency = rank_deficiency)

    DK <- ncol(Q_base)
    r  <- as.integer(DK - (rank_deficiency %||% 0L))
    r  <- max(r, 1L)

    return(elbo - 0.5 * lambda * quad + 0.5 * (logdetQb + r * log(lambda)))
  }

  elbo
}

# =========================
# E-M Step
# =========================

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
  pi_new <- Nk / sum(Nk)

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
        q_eff <- rank_deficiency / d
        q_eff <- as.integer(round(q_eff))
        denom <- n + K - q_eff

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


# =========================
# Construct smooth_em object (exported)
# =========================

#' Construct a smooth_em object from an EM_algorithm fit
#'
#' @param fit Output list from \code{EM_algorithm()}.
#' @param Q_prior Optional precision matrix used in fitting (legacy; discouraged for continuing).
#' @param Q_base Optional *base* precision matrix (without lambda). Recommended.
#' @param lambda Optional penalty strength; used with \code{Q_base}.
#' @param q Optional RW order, if relevant.
#' @param ridge Optional ridge used in building \code{Q_base}.
#' @param modelName Covariance model used in M-step (e.g. "EEI").
#' @param relative_lambda Logical; whether relative-lambda scaling is used.
#' @param rank_deficiency Rank deficiency used in generalized logdet / EEI update.
#' @param eigen_tol Optional tolerance for generalized logdet.
#' @param nugget Diagonal jitter used in covariance updates.
#' @param max_inner,inner_tol M-step inner loop controls.
#' @param meta Optional list of extra metadata to store.
#'
#' @return An object of class \code{smooth_em}.
#' @export
as_smooth_em <- function(
    fit,
    Q_prior = NULL,
    Q_base  = NULL,
    lambda  = NULL,
    q = NULL,
    ridge = NULL,
    modelName = NULL,
    relative_lambda = FALSE,
    rank_deficiency = 0,
    eigen_tol = NULL,
    nugget = 0,
    max_inner = 10,
    inner_tol = 1e-6,
    meta = NULL
) {
  if (!is.list(fit) || is.null(fit$params) || is.null(fit$gamma)) {
    stop("fit must look like EM_algorithm output with fields $params and $gamma.")
  }

  params <- fit$params
  if (is.null(params$invSigma) || is.null(params$logdet)) {
    params <- init_cov_cache_fast(params)
  }

  # normalize prior storage: prefer Q_base + lambda
  if (is.null(Q_base) && !is.null(Q_prior)) {
    Q_base <- Q_prior
    if (is.null(lambda)) lambda <- 1
  }
  if (is.null(lambda)) lambda <- 1
  lambda <- as.numeric(lambda)

  iter0 <- length(fit$elbo_trace %||% numeric(0))
  lambda_trace0 <- if (iter0 > 0L) rep(lambda, iter0) else numeric(0)

  obj <- list(
    params = params,
    gamma  = fit$gamma,
    data   = fit$data %||% NULL,

    elbo_trace   = fit$elbo_trace   %||% numeric(0),
    loglik_trace = fit$loglik_trace %||% numeric(0),
    lambda_trace = fit$lambda_trace %||% lambda_trace0,

    iter = iter0,

    prior = list(
      Q_base = Q_base,
      lambda = lambda,
      Q_prior = Q_prior,  # legacy

      rw_q = q,
      ridge = ridge,
      rank_deficiency = rank_deficiency
    ),

    control = list(
      modelName = modelName,
      relative_lambda = isTRUE(relative_lambda),
      eigen_tol = eigen_tol,
      nugget = nugget,
      max_inner = as.integer(max_inner),
      inner_tol = inner_tol,

      # IMPORTANT: let do_smoothEM decide whether to iterate MSTEP
      mstep_iterate_once = TRUE,

      # bounds (used by do_smoothEM adaptive update)
      lambda_min = 1e-8,
      lambda_max = 1e8
    ),

    meta = meta
  )

  class(obj) <- "smooth_em"
  obj
}


# =========================
# do_smoothEM (exported)
# =========================

#' Run SmoothEM for a given number of iterations on a smooth_em object
#'
#' @param object A \code{smooth_em} object created by \code{as_smooth_em()}.
#' @param data Numeric matrix (n x d).
#' @param iter Integer >= 1; number of (E-step + M-step) iterations to run.
#' @param record Logical; whether to append objective values to traces.
#' @param check_decrease Logical; if TRUE, rollback if ELBO decreases materially.
#' @param tol_decrease Numeric; tolerance for considering ELBO decrease (default 1e-10).
#' @param adaptive Logical; if TRUE, update \code{lambda} each iteration (profile-style update).
#' @param lambda_min,lambda_max Bounds for adaptive lambda (ignored if adaptive=FALSE).
#' @param verbose Logical.
#'
#' @return Updated \code{smooth_em} object.
#' @export
do_smoothEM <- function(object,
                        data = NULL,
                        iter = 1,
                        record = TRUE,
                        check_decrease = TRUE,
                        tol_decrease = 1e-10,
                        adaptive = TRUE,
                        lambda_min = NULL,
                        lambda_max = NULL,
                        verbose = FALSE) {

  if (!inherits(object, "smooth_em")) stop("object must be a 'smooth_em' object.")

  if (is.null(data)) {
    data <- object$data
    if (is.null(data)) stop("data must be provided either in the object or as an argument.")
  }
  data <- as.matrix(data)

  iter <- as.integer(iter)
  if (length(iter) != 1L || is.na(iter) || iter < 1L) stop("iter must be an integer >= 1.")

  # pull settings
  Q_base          <- object$prior$Q_base %||% NULL
  lambda          <- as.numeric(object$prior$lambda %||% 1)
  rank_deficiency <- object$prior$rank_deficiency %||% 0

  modelName       <- object$control$modelName %||% "VVV"
  relative_lambda <- isTRUE(object$control$relative_lambda)
  eigen_tol       <- object$control$eigen_tol
  nugget          <- object$control$nugget %||% 0
  max_inner       <- object$control$max_inner %||% 10
  inner_tol       <- object$control$inner_tol %||% 1e-6

  # let the object control MSTEP iterate_once
  iterate_once <- isTRUE(object$control$mstep_iterate_once %||% TRUE)

  # bounds: argument overrides object$control if provided
  if (is.null(lambda_min)) lambda_min <- object$control$lambda_min %||% 1e-8
  if (is.null(lambda_max)) lambda_max <- object$control$lambda_max %||% 1e8
  lambda_min <- as.numeric(lambda_min)
  lambda_max <- as.numeric(lambda_max)

  eps_quad <- 1e-12

  # ensure cached params
  params <- object$params
  if (is.null(params$invSigma) || is.null(params$logdet)) {
    params <- init_cov_cache_fast(params)
  }

  # main loop
  for (tt in seq_len(iter)) {

    gamma <- ESTEP(data, params)

    Q_prior_for_M <- if (is.null(Q_base)) NULL else (lambda * Q_base)

    new_params <- MSTEP(
      data = data,
      gamma = gamma,
      params = params,
      Q_prior = Q_prior_for_M,
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

    # adaptive lambda (opt-in)
    if (isTRUE(adaptive) && !is.null(Q_base)) {
      lam_star <- .estimate_lambda_star(
        params = params,
        Q_base = Q_base,
        relative_lambda = relative_lambda,
        rank_deficiency = rank_deficiency,
        eps_quad = eps_quad,
        lambda_min = lambda_min,
        lambda_max = lambda_max
      )
      if (is.finite(lam_star)) {
        lambda <- lam_star
        object$prior$lambda <- lambda
      }
    }

    if (isTRUE(record)) {
      # factorized evaluation: use Q_base_eval + lambda
      Q_base_eval <- .compute_Q_base_eval(params, Q_base, relative_lambda)

      ll <- compute_log_joint_observed(
        X = data, params = params,
        Q_prior = NULL,
        Q_base = Q_base_eval,
        lambda = lambda,
        eigen_tol = eigen_tol,
        rank_deficiency = rank_deficiency
      )

      elbo <- compute_penalized_ELBO(
        X = data, Gamma = gamma, params = params,
        Q_prior = NULL,
        Q_base = Q_base_eval,
        lambda = lambda,
        eigen_tol = eigen_tol,
        rank_deficiency = rank_deficiency
      )

      prev_elbo <- tail(object$elbo_trace %||% numeric(0), 1)

      if (isTRUE(check_decrease) && length(prev_elbo) == 1L && is.finite(prev_elbo)) {
        if (prev_elbo - elbo > tol_decrease) {
          if (isTRUE(verbose)) {
            cat(sprintf("ELBO decreased at step %d: %.6f -> %.6f. Rolling back.\n",
                        tt, prev_elbo, elbo))
          }
          params <- last_params
          gamma  <- ESTEP(data, params)
          object$params <- params
          object$gamma  <- gamma
          return(object)
        }
      }

      object$loglik_trace <- c(object$loglik_trace %||% numeric(0), ll)
      object$elbo_trace   <- c(object$elbo_trace   %||% numeric(0), elbo)
      object$lambda_trace <- c(object$lambda_trace %||% numeric(0), lambda)
      object$iter <- length(object$elbo_trace)

      if (isTRUE(verbose)) {
        cat(sprintf("do_smoothEM step %d/%d: penLogLik=%.6f, ELBO=%.6f, lambda=%.4g%s\n",
                    tt, iter, ll, elbo, lambda,
                    if (isTRUE(adaptive)) " (adaptive)" else ""))
      }
    }

    object$gamma <- gamma
  }

  object$params <- params
  object$prior$Q_base <- Q_base
  object$prior$lambda <- lambda
  object
}




