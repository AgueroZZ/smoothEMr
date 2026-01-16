# ============================================================
# csmooth_em class: constructor
# ============================================================

#' Construct a csmooth_em object
#'
#' @description
#' Create a coordinate-specific SmoothEM object (\code{csmooth_em}) with
#' diagonal covariance and a separable RW prior along the K dimension.
#' Each coordinate \eqn{j} has its own penalty strength \eqn{\lambda_j}.
#'
#' Supported covariance models (diagonal only):
#' \itemize{
#'   \item \code{"homoskedastic"}: \eqn{\sigma^2_j} shared across clusters (but varies across coordinates).
#'   \item \code{"heteroskedastic"}: \eqn{\sigma^2_{j,k}} varies across both coordinates and clusters.
#' }
#'
#' @param params List with fields:
#'   \itemize{
#'     \item \code{pi}: length-K mixing proportions.
#'     \item \code{mu}: list of length K; each element is a length-d mean vector.
#'     \item \code{sigma2}: either
#'       \itemize{
#'         \item length-d numeric vector (homoskedastic), or
#'         \item d-by-K numeric matrix (heteroskedastic).
#'       }
#'   }
#' @param gamma (optional) n-by-K responsibility matrix.
#' @param data  (optional) n-by-d data matrix.
#' @param Q_K   K-by-K base precision matrix along components (RW prior along K).
#'              Should be built with \code{lambda = 1}. Coordinate penalties are in \code{lambda_vec}.
#' @param lambda_vec length-d nonnegative vector of per-coordinate lambdas.
#' @param rw_q RW order along K (e.g. 1 or 2). Used as rank deficiency in generalized logdet.
#' @param ridge Ridge used in building \code{Q_K} (stored for provenance).
#' @param modelName One of \code{"homoskedastic"} or \code{"heteroskedastic"}.
#' @param relative_lambda Logical; if TRUE, scale the prior for coordinate j by
#'   \eqn{1/\sigma_j^2} (homoskedastic) or \eqn{1/\bar\sigma_j^2} (heteroskedastic; see details).
#' @param nugget Nonnegative jitter added to variances after updates.
#' @param eigen_tol Optional tolerance for generalized logdet.
#' @param meta Optional list of metadata.
#'
#' @details
#' For \code{modelName="heteroskedastic"} and \code{relative_lambda=TRUE}, we use
#' \eqn{\bar\sigma_j^2 = \sum_k \pi_k \sigma^2_{j,k}} as a reference scale for the
#' prior scaling. This is a pragmatic analogue of the EEI-style scaling.
#'
#' @return An object of class \code{csmooth_em}.
#' @export
as_csmooth_em <- function(
    params,
    gamma = NULL,
    data  = NULL,
    Q_K,
    lambda_vec,
    rw_q = 2,
    ridge = 0,
    modelName = c("homoskedastic", "heteroskedastic"),
    relative_lambda = TRUE,
    nugget = 0,
    eigen_tol = NULL,
    meta = NULL
) {
  modelName <- match.arg(modelName)

  if (!is.list(params) || is.null(params$pi) || is.null(params$mu) || is.null(params$sigma2)) {
    stop("params must be a list with fields pi, mu, sigma2.")
  }
  K <- length(params$pi)
  if (length(params$mu) != K) stop("params$mu must be a list of length K.")
  d <- length(params$mu[[1]])

  if (!is.matrix(Q_K) || nrow(Q_K) != K || ncol(Q_K) != K) stop("Q_K must be a K-by-K matrix.")
  lambda_vec <- as.numeric(lambda_vec)
  if (length(lambda_vec) != d) stop("lambda_vec must have length d (ncol(data)).")
  if (any(!is.finite(lambda_vec)) || any(lambda_vec < 0)) stop("lambda_vec must be finite and >= 0.")

  # sigma2 shape checks
  if (modelName == "homoskedastic") {
    sigma2 <- as.numeric(params$sigma2)
    if (length(sigma2) != d) stop("homoskedastic requires sigma2 as length-d vector.")
  } else {
    sigma2 <- as.matrix(params$sigma2)
    if (!all(dim(sigma2) == c(d, K))) stop("heteroskedastic requires sigma2 as d-by-K matrix.")
  }

  obj <- list(
    params = list(
      pi = as.numeric(params$pi),
      mu = lapply(params$mu, function(x) as.numeric(x)),
      sigma2 = sigma2
    ),
    gamma = gamma,
    data  = data,

    elbo_trace   = numeric(0),
    loglik_trace = numeric(0),
    lambda_trace = list(),   # store per-iter lambda_vec if you want (optional)
    iter = 0L,

    prior = list(
      Q_K = Q_K,
      lambda_vec = lambda_vec,
      rw_q = as.integer(rw_q),
      ridge = ridge
    ),

    control = list(
      modelName = modelName,
      relative_lambda = isTRUE(relative_lambda),
      nugget = as.numeric(nugget),
      eigen_tol = eigen_tol,

      # adaptive-lambda knobs
      adapt_lambda = TRUE,
      lambda_min = 1e-10,
      lambda_max = 1e10
    ),

    meta = meta
  )

  class(obj) <- "csmooth_em"
  obj
}


# ============================================================
# csmooth_em: E-step (diagonal covariances)
# ============================================================

#' Internal E-step for csmooth_em
#'
#' @param X n-by-d data matrix
#' @param params csmooth params list (pi, mu, sigma2)
#' @param modelName "homoskedastic" or "heteroskedastic"
#'
#' @return n-by-K responsibility matrix
#' @keywords internal
ESTEP_csmooth <- function(X, params, modelName) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  K <- length(params$pi)

  Mu <- do.call(cbind, params$mu)  # d x K
  pi_vec <- params$pi

  log_gamma <- matrix(0, n, K)

  if (modelName == "homoskedastic") {
    sigma2 <- as.numeric(params$sigma2)
    sigma2 <- pmax(sigma2, .Machine$double.eps)

    # log det is shared across k: sum_j log sigma2_j
    logdet <- sum(log(sigma2))
    const  <- -0.5 * (d * log(2*pi) + logdet)

    for (k in seq_len(K)) {
      diff <- sweep(X, 2, Mu[, k], "-")
      quad <- rowSums((diff^2) / matrix(sigma2, n, d, byrow = TRUE))
      log_gamma[, k] <- log(pi_vec[k]) + const - 0.5 * quad
    }

  } else {
    sigma2_mat <- as.matrix(params$sigma2)  # d x K

    for (k in seq_len(K)) {
      sig2 <- pmax(sigma2_mat[, k], .Machine$double.eps)
      logdet <- sum(log(sig2))
      const  <- -0.5 * (d * log(2*pi) + logdet)

      diff <- sweep(X, 2, Mu[, k], "-")
      quad <- rowSums((diff^2) / matrix(sig2, n, d, byrow = TRUE))
      log_gamma[, k] <- log(pi_vec[k]) + const - 0.5 * quad
    }
  }

  lse <- matrixStats::rowLogSumExps(log_gamma)
  exp(log_gamma - lse)
}


# ============================================================
# csmooth_em: objective utilities (optional but recommended)
# ============================================================

#' Internal: compute Q_base_j for coordinate j under relative_lambda
#' @keywords internal
.compute_Qbase_j <- function(pi, sigma2, Q_K, relative_lambda, modelName, j) {
  if (!relative_lambda) return(Q_K)

  if (modelName == "homoskedastic") {
    s2 <- sigma2[j]
  } else {
    # reference scale: sum_k pi_k * sigma2_{j,k}
    s2 <- sum(pi * sigma2[j, ])
  }
  s2 <- pmax(s2, .Machine$double.eps)
  Q_K / s2
}

#' Penalized observed-data objective for csmooth_em (optional)
#' @keywords internal
compute_log_joint_observed_csmooth <- function(X, params, Q_K, lambda_vec,
                                               modelName, relative_lambda,
                                               eigen_tol = NULL, rw_q = 0L) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  K <- length(params$pi)

  # observed log-likelihood: log sum_k pi_k N(x_i | mu_k, Sigma_k)
  Mu <- do.call(cbind, params$mu)  # d x K
  logdens <- matrix(0, n, K)

  if (modelName == "homoskedastic") {
    sigma2 <- pmax(as.numeric(params$sigma2), .Machine$double.eps)  # length d
    logdet <- sum(log(sigma2))
    const  <- -0.5 * (d * log(2 * base::pi) + logdet)

    invsig2 <- 1 / sigma2
    for (k in seq_len(K)) {
      diff <- sweep(X, 2, Mu[, k], "-")
      quad <- rowSums((diff^2) * matrix(invsig2, n, d, byrow = TRUE))
      logdens[, k] <- log(pmax(params$pi[k], .Machine$double.eps)) + const - 0.5 * quad
    }

  } else {
    sigma2_mat <- pmax(as.matrix(params$sigma2), .Machine$double.eps)  # d x K
    for (k in seq_len(K)) {
      sig2 <- sigma2_mat[, k]
      logdet <- sum(log(sig2))
      const  <- -0.5 * (d * log(2 * base::pi) + logdet)

      invsig2 <- 1 / sig2
      diff <- sweep(X, 2, Mu[, k], "-")
      quad <- rowSums((diff^2) * matrix(invsig2, n, d, byrow = TRUE))
      logdens[, k] <- log(pmax(params$pi[k], .Machine$double.eps)) + const - 0.5 * quad
    }
  }

  ll <- sum(matrixStats::rowLogSumExps(logdens))

  # prior penalty + normalization
  penalty_sum <- 0
  logdet_sum  <- 0
  for (j in seq_len(d)) {
    lam <- lambda_vec[j]
    if (!is.finite(lam) || lam <= 0) next

    Qb  <- .compute_Qbase_j(params$pi, params$sigma2, Q_K, relative_lambda, modelName, j)
    Qj  <- lam * Qb

    mu_j <- vapply(seq_len(K), function(k) params$mu[[k]][j], numeric(1))
    penalty_sum <- penalty_sum + 0.5 * as.numeric(crossprod(mu_j, Qj %*% mu_j))
    logdet_sum  <- logdet_sum + generalized_logdet(Q = Qj, eigen_tol = eigen_tol, rank_deficiency = rw_q)
  }

  ll - penalty_sum + 0.5 * logdet_sum
}

#' Penalized ELBO for csmooth_em (optional)
#' @keywords internal
compute_penalized_ELBO_csmooth <- function(X, Gamma, params, Q_K, lambda_vec,
                                           modelName, relative_lambda,
                                           eigen_tol = NULL, rw_q = 0L) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  K <- ncol(Gamma)

  Mu <- do.call(cbind, params$mu)
  log_gamma <- matrix(0, n, K)

  if (modelName == "homoskedastic") {
    sigma2 <- pmax(as.numeric(params$sigma2), .Machine$double.eps)
    logdet <- sum(log(sigma2))
    const  <- -0.5 * (d * log(2*pi) + logdet)

    for (k in seq_len(K)) {
      diff <- sweep(X, 2, Mu[, k], "-")
      quad <- rowSums((diff^2) / matrix(sigma2, n, d, byrow = TRUE))
      log_gamma[, k] <- log(params$pi[k]) + const - 0.5 * quad
    }
  } else {
    sigma2_mat <- as.matrix(params$sigma2)
    for (k in seq_len(K)) {
      sig2 <- pmax(sigma2_mat[, k], .Machine$double.eps)
      logdet <- sum(log(sig2))
      const  <- -0.5 * (d * log(2*pi) + logdet)

      diff <- sweep(X, 2, Mu[, k], "-")
      quad <- rowSums((diff^2) / matrix(sig2, n, d, byrow = TRUE))
      log_gamma[, k] <- log(params$pi[k]) + const - 0.5 * quad
    }
  }

  Qval <- sum(Gamma * log_gamma)
  entropy <- -sum(Gamma * log(Gamma + 1e-300))
  elbo <- Qval + entropy

  penalty_sum <- 0
  logdet_sum  <- 0
  for (j in seq_len(d)) {
    lam <- lambda_vec[j]
    if (!is.finite(lam) || lam <= 0) next

    Qb  <- .compute_Qbase_j(params$pi, params$sigma2, Q_K, relative_lambda, modelName, j)
    Qj  <- lam * Qb

    mu_j <- vapply(seq_len(K), function(k) params$mu[[k]][j], numeric(1))
    penalty_sum <- penalty_sum + 0.5 * as.numeric(crossprod(mu_j, Qj %*% mu_j))
    logdet_sum  <- logdet_sum + generalized_logdet(Q = Qj, eigen_tol = eigen_tol, rank_deficiency = rw_q)
  }

  elbo - penalty_sum + 0.5 * logdet_sum
}


# ============================================================
# csmooth_em: coordinate-wise M-step (diagonal covariance)
# ============================================================

#' Internal M-step for csmooth_em (coordinate-wise)
#'
#' @param X n-by-d data
#' @param Gamma n-by-K responsibilities
#' @param params current params (pi, mu, sigma2)
#' @param Q_K K-by-K base precision
#' @param lambda_vec length-d vector
#' @param modelName "homoskedastic" or "heteroskedastic"
#' @param relative_lambda logical
#' @param nugget nonnegative
#' @param rw_q rank deficiency along K (RW order)
#' @param iterate_once logical; if FALSE, you could add inner loops later (kept simple here)
#'
#' @return updated params list (pi, mu, sigma2)
#' @keywords internal
MSTEP_csmooth <- function(X, Gamma, params, Q_K, lambda_vec,
                          modelName, relative_lambda,
                          nugget = 0, rw_q = 0L, iterate_once = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  K <- ncol(Gamma)

  Nk <- colSums(Gamma)
  Nk[Nk < 1e-8] <- 1e-8
  pi_new <- Nk / sum(Nk)

  GX  <- t(X) %*% Gamma      # d x K
  GX2 <- t(X^2) %*% Gamma    # d x K

  # current sigma2
  if (modelName == "homoskedastic") {
    sigma2 <- pmax(as.numeric(params$sigma2), .Machine$double.eps)  # length d
  } else {
    sigma2 <- pmax(as.matrix(params$sigma2), .Machine$double.eps)   # d x K
  }

  mu_mat <- matrix(0, nrow = d, ncol = K)

  # ---- mu update: per-coordinate KxK solve
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
        Qb <- .compute_Qbase_j(pi_new, sigma2, Q_K, relative_lambda, modelName, j)  # already /s2j if relative
        A <- diag(A_diag, K, K) + (lam * Qb)
        # SPD solve
        L <- chol(A)
        y <- forwardsolve(t(L), rhs)
        mu_mat[j, ] <- as.numeric(backsolve(L, y))
      }

    } else { # heteroskedastic
      s2jk <- as.numeric(sigma2[j, ]) # length K
      rhs  <- as.numeric(GX[j, ]) / s2jk
      A_diag <- Nk / s2jk

      if (lam <= 0 || is.null(Q_K)) {
        mu_mat[j, ] <- as.numeric(GX[j, ]) / Nk
      } else {
        Qb <- .compute_Qbase_j(pi_new, sigma2, Q_K, relative_lambda, modelName, j)
        A <- diag(A_diag, K, K) + (lam * Qb)
        L <- chol(A)
        y <- forwardsolve(t(L), rhs)
        mu_mat[j, ] <- as.numeric(backsolve(L, y))
      }
    }
  }

  mu_list <- lapply(seq_len(K), function(k) mu_mat[, k])

  # ---- sigma2 update
  if (modelName == "homoskedastic") {
    V <- numeric(d)
    for (j in seq_len(d)) {
      # V_j = sum_k [ GX2[j,k] - 2 mu_{j,k} GX[j,k] + mu_{j,k}^2 Nk[k] ]
      term1 <- sum(GX2[j, ])
      term2 <- 2 * sum(mu_mat[j, ] * GX[j, ])
      term3 <- sum((mu_mat[j, ]^2) * Nk)
      V[j] <- term1 - term2 + term3
    }

    if (relative_lambda && !is.null(Q_K)) {
      # add prior contribution U_j = lambda_j * mu_j' Q_K mu_j
      U <- vapply(seq_len(d), function(j) {
        lam <- lambda_vec[j]
        if (!is.finite(lam) || lam <= 0) return(0)
        mj <- mu_mat[j, ]
        lam * as.numeric(crossprod(mj, Q_K %*% mj))
      }, numeric(1))

      denom <- (n + K - as.integer(rw_q))
      denom <- max(denom, 1L)
      sigma2_new <- (V + U) / denom + nugget
    } else {
      sigma2_new <- V / n + nugget
    }

    sigma2_new <- pmax(sigma2_new, .Machine$double.eps)
    sigma2_out <- sigma2_new

  } else {
    # heteroskedastic: per (j,k) variance
    sigma2_new <- matrix(0, nrow = d, ncol = K)
    for (k in seq_len(K)) {
      # For each coordinate j:
      # sigma2_{j,k} = [sum_i gamma_ik (x_ij - mu_jk)^2] / Nk[k]
      # Use sufficient stats: sum gamma x^2 - 2 mu sum gamma x + mu^2 Nk
      sigma2_new[, k] <- (GX2[, k] - 2 * mu_mat[, k] * GX[, k] + (mu_mat[, k]^2) * Nk[k]) / Nk[k] + nugget
    }
    sigma2_out <- pmax(sigma2_new, .Machine$double.eps)
  }

  list(pi = pi_new, mu = mu_list, sigma2 = sigma2_out)
}




#' Run csmoothEM iterations
#'
#' @description
#' Runs EM iterations for coordinate-specific SmoothEM (\code{csmooth_em}).
#' Supports adaptive updates of the coordinate-wise smoothing parameters \eqn{\lambda_j}.
#'
#' When \code{adaptive = "ml"}, this function dispatches to
#' \code{\link{do_csmoothEM_ml_collapsed}} (collapsed-ML variant). In that case,
#' \code{sigma_update} controls whether \eqn{\sigma^2} is updated via the standard
#' M-step or via collapsed-ML optimization.
#'
#' @param object A \code{csmooth_em} object.
#' @param data Optional n-by-d data matrix. If NULL, uses \code{object$data}.
#' @param iter Integer \eqn{\ge 1}; number of EM iterations to run.
#' @param record Logical; if TRUE, record traces.
#' @param adaptive Adaptive option controlling \eqn{\lambda} updates:
#'   \describe{
#'     \item{\code{NULL}}{Inherit \code{object$control$adaptive}.}
#'     \item{\code{"none"}}{No adaptive update (fixed \code{lambda_vec}).}
#'     \item{\code{"prior"}}{Prior-based update for \eqn{\lambda_j} using the current M-step \eqn{\mu}.}
#'     \item{\code{"ml"}}{Collapsed-ML update for \eqn{\lambda_j} (dispatches to
#'       \code{\link{do_csmoothEM_ml_collapsed}}).}
#'   }
#'   A logical value is accepted for backward compatibility: \code{TRUE} is treated as
#'   \code{"prior"} and \code{FALSE} as \code{"none"}.
#' @param lambda_min,lambda_max Positive bounds used for adaptive \eqn{\lambda} updates.
#' @param sigma_update Character. Only used when \code{adaptive="ml"} (i.e., in the
#'   collapsed-ML path). Either \code{"mstep"} or \code{"ml"}; see
#'   \code{\link{do_csmoothEM_ml_collapsed}}.
#' @param sigma_min,sigma_max Positive bounds for \code{sigma2} when \code{adaptive="ml"}
#'   and \code{sigma_update="ml"}.
#' @param verbose Logical; if TRUE, print a short progress line each iteration.
#'
#' @return Updated \code{csmooth_em} object.
#' @export
do_csmoothEM <- function(object,
                         data = NULL,
                         iter = 1,
                         record = TRUE,
                         adaptive = NULL,          # NULL = inherit; "none" = off; "prior"/"ml" explicit
                         lambda_min = NULL,
                         lambda_max = NULL,
                         sigma_update = c("mstep", "ml"),  # only used when adaptive == "ml"
                         sigma_min = 1e-10,
                         sigma_max = 1e10,
                         verbose = FALSE) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!inherits(object, "csmooth_em")) stop("object must be a 'csmooth_em' object.")

  if (is.null(data)) {
    data <- object$data
    if (is.null(data)) stop("data must be provided either in object$data or as an argument.")
  }
  X <- as.matrix(data)

  iter <- as.integer(iter)
  if (length(iter) != 1L || is.na(iter) || iter < 1L) stop("iter must be integer >= 1.")

  modelName       <- object$control$modelName %||% "homoskedastic"
  relative_lambda <- isTRUE(object$control$relative_lambda)
  nugget          <- object$control$nugget %||% 0
  eigen_tol       <- object$control$eigen_tol
  rw_q            <- as.integer(object$prior$rw_q %||% 0L)

  Q_K        <- object$prior$Q_K
  lambda_vec <- as.numeric(object$prior$lambda_vec)

  d <- ncol(X)
  K <- length(object$params$pi)

  if (is.null(lambda_min)) lambda_min <- object$control$lambda_min %||% 1e-10
  if (is.null(lambda_max)) lambda_max <- object$control$lambda_max %||% 1e10
  lambda_min <- as.numeric(lambda_min)
  lambda_max <- as.numeric(lambda_max)
  if (!is.finite(lambda_min) || !is.finite(lambda_max) || lambda_min <= 0 || lambda_max <= 0 || lambda_min > lambda_max) {
    stop("lambda_min/lambda_max must be positive finite numbers with lambda_min <= lambda_max.")
  }

  # normalize adaptive option (NULL = inherit; "none" = off)
  if (is.null(adaptive)) adaptive <- object$control$adaptive %||% NULL
  if (is.logical(adaptive)) adaptive <- if (isTRUE(adaptive)) "prior" else "none"
  if (is.character(adaptive) && identical(adaptive, "none")) adaptive <- NULL
  if (!is.null(adaptive)) adaptive <- match.arg(adaptive, choices = c("prior", "ml"))


  # Dispatch: collapsed-ML path (this is the only place sigma_update matters)
  if (!is.null(adaptive) && identical(adaptive, "ml")) {
    if (missing(sigma_update)) sigma_update <- object$control$sigma_update %||% "mstep"
    sigma_update <- match.arg(sigma_update, choices = c("mstep", "ml"))

    if (missing(sigma_min)) sigma_min <- object$control$sigma_min %||% 1e-10
    if (missing(sigma_max)) sigma_max <- object$control$sigma_max %||% 1e10

    return(do_csmoothEM_ml_collapsed(
      object = object,
      data = data,
      iter = iter,
      record = record,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_update = sigma_update,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      verbose = verbose
    ))
  }

  eps_quad <- 1e-12
  r_rank <- max(K - rw_q, 1L)

  params <- object$params

  # ------------------------------------------------------------------
  # Helpers for collapsed objective C = penELBO + const - 0.5*log|H|
  # ------------------------------------------------------------------

  .compute_logdetH <- function(Gamma, params, lambda_vec) {
    Nk <- colSums(Gamma)
    Nk[Nk < 1e-8] <- 1e-8

    if (modelName == "homoskedastic") {
      sigma2 <- pmax(as.numeric(params$sigma2), .Machine$double.eps)  # length d
    } else {
      sigma2 <- pmax(as.matrix(params$sigma2), .Machine$double.eps)   # d x K
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

  # ------------------------------------------------------------------
  # Helper: refresh mu to be MAP under current Gamma, sigma2, lambda_vec
  # (keeps sigma2 fixed; refreshes pi from Gamma)
  # ------------------------------------------------------------------
  .mu_refresh <- function(X, Gamma, params, lambda_vec) {
    Nk <- colSums(Gamma)
    Nk[Nk < 1e-8] <- 1e-8
    pi_new <- Nk / sum(Nk)

    GX <- t(X) %*% Gamma  # d x K

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
          A  <- diag(A_diag, K, K) + lam * Qb
          A  <- 0.5 * (A + t(A))
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
          A  <- diag(A_diag, K, K) + lam * Qb
          A  <- 0.5 * (A + t(A))
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

  # ------------------------------------------------------------------
  # Helper used only for adaptive="ml": per-coordinate reduced objective in log-lambda
  # (no likelihood constants; we only need it for optimization)
  # ------------------------------------------------------------------
  .ml_reduced_obj_loglam_j <- function(loglam, Dj, bj, Qb, logdet_Qb, r_rank, K) {
    lam <- exp(loglam)
    A <- diag(Dj, K, K) + lam * Qb
    A <- 0.5 * (A + t(A))
    L <- chol(A)
    logdetA <- 2 * sum(log(diag(L)))
    y <- forwardsolve(t(L), bj)
    alpha <- backsolve(L, y)
    quad_term <- sum(bj * alpha)
    0.5 * (r_rank * loglam + logdet_Qb) - 0.5 * logdetA + 0.5 * quad_term
  }

  # helper to build Dj/bj given Gamma and sigma2
  .build_D_b <- function(j, params, Nk, GX) {
    if (modelName == "homoskedastic") {
      s2j <- pmax(as.numeric(params$sigma2[j]), .Machine$double.eps)
      Dj <- Nk / s2j
      bj <- as.numeric(GX[j, ]) / s2j
    } else {
      s2jk <- pmax(as.numeric(params$sigma2[j, ]), .Machine$double.eps)
      Dj <- Nk / s2jk
      bj <- as.numeric(GX[j, ]) / s2jk
    }
    list(Dj = as.numeric(Dj), bj = as.numeric(bj))
  }

  # ------------------------------------------------------------------
  # Main loop
  # ------------------------------------------------------------------
  for (tt in seq_len(iter)) {

    # E-step
    Gamma <- ESTEP_csmooth(X, params, modelName)

    # M-step (updates pi, mu, sigma2 given current lambda_vec)
    params <- MSTEP_csmooth(
      X = X, Gamma = Gamma, params = params,
      Q_K = Q_K, lambda_vec = lambda_vec,
      modelName = modelName, relative_lambda = relative_lambda,
      nugget = nugget, rw_q = rw_q, iterate_once = TRUE
    )

    # Adaptive lambda update
    if (!is.null(adaptive)) {

      Nk <- colSums(Gamma)
      Nk[Nk < 1e-8] <- 1e-8
      GX <- t(X) %*% Gamma

      if (adaptive == "prior") {
        mu_mat <- do.call(cbind, params$mu)  # d x K
        for (j in seq_len(d)) {
          Qb <- .compute_Qbase_j(params$pi, params$sigma2, Q_K, relative_lambda, modelName, j)
          mj <- as.numeric(mu_mat[j, ])
          quad_prior <- as.numeric(crossprod(mj, Qb %*% mj))
          lam_star <- r_rank / pmax(quad_prior, eps_quad)
          lambda_vec[j] <- min(max(lam_star, lambda_min), lambda_max)
        }

      } else if (adaptive == "ml") {

        for (j in seq_len(d)) {
          Qb <- .compute_Qbase_j(params$pi, params$sigma2, Q_K, relative_lambda, modelName, j)
          logdet_Qb <- generalized_logdet(Q = Qb, eigen_tol = eigen_tol, rank_deficiency = rw_q)

          tmp <- .build_D_b(j, params, Nk, GX)
          opt <- optimize(
            f = .ml_reduced_obj_loglam_j,
            interval = log(c(lambda_min, lambda_max)),
            maximum = TRUE,
            Dj = tmp$Dj, bj = tmp$bj, Qb = Qb,
            logdet_Qb = logdet_Qb, r_rank = r_rank, K = K
          )
          lambda_vec[j] <- exp(opt$maximum)
        }
      }

      # After changing lambda, refresh mu so that mu = MAP under current lambda (needed for C/ml_trace)
      params <- .mu_refresh(X, Gamma, params, lambda_vec)
    }

    # Record traces
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
      object$elbo_trace   <- c(object$elbo_trace,   pelbo)

      # unified "ml_trace" := collapsed objective C
      if (is.null(object$ml_trace)) object$ml_trace <- numeric(0)
      object$ml_trace <- c(object$ml_trace, Cres$C)

      # also store log|H| for diagnostics
      if (is.null(object$logdetH_trace)) object$logdetH_trace <- numeric(0)
      object$logdetH_trace <- c(object$logdetH_trace, Cres$logdetH)

      object$iter <- length(object$elbo_trace)

      if (is.null(object$lambda_trace)) object$lambda_trace <- list()
      object$lambda_trace[[object$iter]] <- lambda_vec

      if (verbose) {
        msg_ad <- if (is.null(adaptive)) "none" else adaptive
        cat(sprintf(
          "csmoothEM %d/%d: penLogLik=%.6f  penELBO=%.6f  C=%.6f  log|H|=%.6f  lambda(range)=[%.3g, %.3g] (adaptive=%s)\n",
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


# ============================================================
# do_csmoothEM_ml_collapsed: unified tracing + consistent names
# ============================================================
#' Run collapsed-ML csmoothEM iterations (homoskedastic only)
#'
#' @description
#' Iteration order (collapsed-ML):
#' 1) E-step: update responsibilities using current (MAP) \eqn{\mu}.
#' 2) Hyper-step: update \eqn{\pi}, \eqn{\lambda}, and optionally \eqn{\sigma^2}.
#' 3) M-step: update MAP \eqn{\mu} given current responsibilities and hyperparameters.
#' 4) Record: penalized observed objective (\code{loglik_trace}), penalized ELBO
#'    (\code{elbo_trace}), collapsed objective \eqn{\mathcal{C}} (\code{ml_trace}),
#'    and \code{logdetH_trace}.
#'
#' @details
#' The recorded \code{ml_trace} is the collapsed objective
#' \deqn{\mathcal{C} \;=\; \text{penELBO} \;+\; \frac{dK}{2}\log(2\pi) \;-\; \frac12 \log|H|,}
#' where \eqn{H} is the Hessian/precision w.r.t. \eqn{\mu} (block-diagonal over coordinates).
#'
#' @param object A \code{csmooth_em} object.
#' @param data Optional n-by-d data matrix. If NULL, uses \code{object$data}.
#' @param iter Integer >= 1 number of iterations.
#' @param record Logical; if TRUE, record traces.
#' @param lambda_min,lambda_max Bounds for lambda optimization.
#' @param sigma_update Either \code{"mstep"} (standard M-step update) or \code{"ml"} (optimize collapsed objective).
#' @param sigma_min,sigma_max Bounds for sigma2 when \code{sigma_update="ml"}.
#' @param verbose Logical.
#'
#' @return Updated \code{csmooth_em} object.
#' @export
do_csmoothEM_ml_collapsed <- function(object,
                                      data = NULL,
                                      iter = 1,
                                      record = TRUE,
                                      lambda_min = NULL,
                                      lambda_max = NULL,
                                      sigma_update = c("mstep", "ml"),
                                      sigma_min = 1e-10,
                                      sigma_max = 1e10,
                                      verbose = FALSE) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!inherits(object, "csmooth_em")) stop("object must be a 'csmooth_em' object.")

  if (is.null(data)) {
    data <- object$data
    if (is.null(data)) stop("data must be provided either in object$data or as an argument.")
  }
  X <- as.matrix(data)

  iter <- as.integer(iter)
  if (length(iter) != 1L || is.na(iter) || iter < 1L) stop("iter must be integer >= 1.")

  modelName       <- object$control$modelName %||% "homoskedastic"
  if (modelName != "homoskedastic") {
    stop("do_csmoothEM_ml_collapsed currently supports modelName='homoskedastic' only.")
  }

  relative_lambda <- isTRUE(object$control$relative_lambda)
  nugget          <- object$control$nugget %||% 0
  eigen_tol       <- object$control$eigen_tol
  rw_q            <- as.integer(object$prior$rw_q %||% 0L)

  Q_K        <- object$prior$Q_K
  lambda_vec <- as.numeric(object$prior$lambda_vec)

  d <- ncol(X)
  K <- length(object$params$pi)

  if (is.null(lambda_min)) lambda_min <- object$control$lambda_min %||% 1e-10
  if (is.null(lambda_max)) lambda_max <- object$control$lambda_max %||% 1e10
  lambda_min <- as.numeric(lambda_min)
  lambda_max <- as.numeric(lambda_max)
  if (!is.finite(lambda_min) || !is.finite(lambda_max) || lambda_min <= 0 || lambda_max <= 0 || lambda_min > lambda_max) {
    stop("lambda_min/lambda_max must be positive finite numbers with lambda_min <= lambda_max.")
  }

  sigma_update <- match.arg(sigma_update)
  sigma_min <- as.numeric(sigma_min)
  sigma_max <- as.numeric(sigma_max)
  if (!is.finite(sigma_min) || !is.finite(sigma_max) || sigma_min <= 0 || sigma_max <= 0 || sigma_min > sigma_max) {
    stop("sigma_min/sigma_max must be positive finite numbers with sigma_min <= sigma_max.")
  }

  r_rank <- max(K - rw_q, 1L)
  params <- object$params

  # ---- helpers ------------------------------------------------

  # log|H| where H is block diagonal over coordinates, with blocks A_j = diag(Nk/s2j) + lam_j Qb_j
  .logdetH <- function(Gamma, params, lambda_vec) {
    Nk <- colSums(Gamma)
    Nk[Nk < 1e-8] <- 1e-8
    sigma2 <- pmax(as.numeric(params$sigma2), .Machine$double.eps)

    out <- 0
    for (j in seq_len(d)) {
      lam <- lambda_vec[j]
      if (!is.finite(lam) || lam <= 0) next
      Qb <- .compute_Qbase_j(params$pi, sigma2, Q_K, relative_lambda, "homoskedastic", j)
      A  <- diag(Nk / sigma2[j], K, K) + lam * Qb
      A  <- 0.5 * (A + t(A))
      L  <- chol(A)
      out <- out + 2 * sum(log(diag(L)))
    }
    out
  }

  # Collapsed objective C = penELBO + (dK/2)log(2pi) - 0.5 log|H|
  .compute_C <- function(X, Gamma, params, lambda_vec, include_constant = TRUE) {
    pelbo <- compute_penalized_ELBO_csmooth(
      X = X, Gamma = Gamma, params = params, Q_K = Q_K, lambda_vec = lambda_vec,
      modelName = "homoskedastic", relative_lambda = relative_lambda,
      eigen_tol = eigen_tol, rw_q = rw_q
    )
    ldH <- .logdetH(Gamma, params, lambda_vec)
    const <- if (isTRUE(include_constant)) 0.5 * d * K * log(2 * base::pi) else 0
    list(C = pelbo + const - 0.5 * ldH, logdetH = ldH, pelbo = pelbo)
  }

  # per-j collapsed objective contribution for optimizing lambda_j, given Gamma and sigma2
  .ell_j_loglam <- function(loglam, j, Nk, GX, GX2, pi_vec, sigma2, Q_K, r_rank) {
    lam <- exp(loglam)
    s2j <- pmax(as.numeric(sigma2[j]), .Machine$double.eps)

    Qb <- .compute_Qbase_j(pi_vec, sigma2, Q_K, relative_lambda, "homoskedastic", j)

    # likelihood constants (scalar!)
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

    as.numeric(C_lik + 0.5 * (r_rank * loglam + logdet_Qb) - 0.5 * logdetA + 0.5 * quad)
  }

  # per-j collapsed objective for optimizing sigma2_j (homoskedastic), given Gamma and lambda_j
  .ell_j_logsig2 <- function(logsig2, j, Nk, GX, GX2, pi_vec, lambda_vec, Q_K, r_rank) {
    s2j <- exp(logsig2)

    sigma2_tmp <- rep(1, d)
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

    as.numeric(C_lik + 0.5 * (r_rank * log(lam) + logdet_Qb) - 0.5 * logdetA + 0.5 * quad)
  }

  # ------------------------------------------------------------

  for (tt in seq_len(iter)) {

    # 1) E-step (based on current MAP mu in params)
    Gamma <- ESTEP_csmooth(X, params, "homoskedastic")

    # sufficient stats for hyper-step
    Nk <- colSums(Gamma); Nk[Nk < 1e-8] <- 1e-8
    GX  <- t(X)   %*% Gamma
    GX2 <- t(X^2) %*% Gamma

    # 2) Hyper-step: update pi (closed form)
    params$pi <- Nk / sum(Nk)

    # 2a) optional sigma2 update
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
          f = .ell_j_logsig2,
          interval = log(c(sigma_min, sigma_max)),
          maximum = TRUE,
          j = j, Nk = Nk, GX = GX, GX2 = GX2,
          pi_vec = params$pi, lambda_vec = lambda_vec,
          Q_K = Q_K, r_rank = r_rank
        )
        sigma2_new[j] <- exp(opt_s2$maximum)
      }
      params$sigma2 <- pmax(sigma2_new + nugget, .Machine$double.eps)
    }

    # 2b) lambda update (per coordinate)
    for (j in seq_len(d)) {
      opt_lam <- optimize(
        f = .ell_j_loglam,
        interval = log(c(lambda_min, lambda_max)),
        maximum = TRUE,
        j = j, Nk = Nk, GX = GX, GX2 = GX2,
        pi_vec = params$pi, sigma2 = params$sigma2,
        Q_K = Q_K, r_rank = r_rank
      )
      lambda_vec[j] <- exp(opt_lam$maximum)
    }

    # 3) M-step: update MAP mu given current Gamma and updated hyperparameters
    params <- MSTEP_csmooth(
      X = X, Gamma = Gamma, params = params,
      Q_K = Q_K, lambda_vec = lambda_vec,
      modelName = "homoskedastic", relative_lambda = relative_lambda,
      nugget = nugget, rw_q = rw_q, iterate_once = TRUE
    )

    # 4) unified record
    if (record) {
      ll <- compute_log_joint_observed_csmooth(
        X, params, Q_K, lambda_vec,
        modelName = "homoskedastic",
        relative_lambda = relative_lambda,
        eigen_tol = eigen_tol,
        rw_q = rw_q
      )

      pelbo <- compute_penalized_ELBO_csmooth(
        X, Gamma, params, Q_K, lambda_vec,
        modelName = "homoskedastic",
        relative_lambda = relative_lambda,
        eigen_tol = eigen_tol,
        rw_q = rw_q
      )

      Cres <- .compute_C(X, Gamma, params, lambda_vec, include_constant = TRUE)

      object$loglik_trace <- c(object$loglik_trace, ll)
      object$elbo_trace   <- c(object$elbo_trace,  pelbo)

      # unified "ml_trace" := collapsed objective C
      if (is.null(object$ml_trace)) object$ml_trace <- numeric(0)
      object$ml_trace <- c(object$ml_trace, Cres$C)

      if (is.null(object$logdetH_trace)) object$logdetH_trace <- numeric(0)
      object$logdetH_trace <- c(object$logdetH_trace, Cres$logdetH)

      object$iter <- length(object$elbo_trace)

      if (is.null(object$lambda_trace)) object$lambda_trace <- list()
      object$lambda_trace[[object$iter]] <- lambda_vec

      if (verbose) {
        cat(sprintf(
          "csmoothEM_ml_collapsed %d/%d: penLogLik=%.6f  penELBO=%.6f  C=%.6f  log|H|=%.6f  lambda(range)=[%.3g, %.3g]  sigma2(range)=[%.3g, %.3g]  sigma_update=%s\n",
          tt, iter, ll, pelbo, Cres$C, Cres$logdetH,
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
# initialize_csmoothEM: user-facing initializer (same feel as initialize_smoothEM)
# ============================================================


#' Initialization from an ordering vector for csmoothEM
#'
#' @description
#' Initializes csmoothEM parameters from a 1D ordering score by:
#' \enumerate{
#'   \item Discretizing \code{ordering_vec} into \code{K} ordered groups (equal/quantile/kmeans).
#'   \item Computing component means \code{mu}.
#'   \item Estimating diagonal variances \code{sigma2} under either:
#'     \itemize{
#'       \item \code{modelName="homoskedastic"}: one variance per coordinate shared across clusters (length \code{d}).
#'       \item \code{modelName="heteroskedastic"}: one variance per coordinate per cluster (\code{d x K}).
#'     }
#' }
#'
#' This function mirrors the discretization and mean construction logic of \code{make_init()},
#' but returns \code{sigma2} (diagonal variances) instead of a list of full covariance matrices.
#'
#' @param X Numeric matrix \code{(n x d)}.
#' @param ordering_vec Numeric vector of length \code{n} (can contain NA).
#' @param K Integer >= 2; number of mixture components.
#' @param modelName Either \code{"homoskedastic"} or \code{"heteroskedastic"}.
#' @param nugget Nonnegative scalar added to variance estimates.
#' @param discretization One of \code{"equal"}, \code{"quantile"}, \code{"kmeans"}.
#' @param na_action How to handle NA in \code{ordering_vec}: \code{"drop"} or \code{"error"}.
#' @param eps Small positive floor for \code{pi} and \code{sigma2}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{pi}: length-\code{K} mixing proportions.
#'   \item \code{mu}: list of length \code{K}, each a length-\code{d} mean vector.
#'   \item \code{sigma2}: diagonal variances; length-\code{d} vector (homoskedastic) or \code{d x K} matrix (heteroskedastic).
#'   \item \code{keep_idx}: row indices kept after NA handling.
#'   \item \code{cluster_rank}: integer vector in \code{1..K} for each kept row.
#' }
#'
#' @export
make_init_csmooth <- function(
    X, ordering_vec, K,
    modelName = c("homoskedastic", "heteroskedastic"),
    nugget = 0,
    discretization = c("equal", "quantile", "kmeans"),
    na_action = c("drop", "error"),
    eps = 1e-12
) {
  X <- as.matrix(X)
  modelName <- match.arg(modelName)
  discretization <- match.arg(discretization)
  na_action <- match.arg(na_action)

  K <- as.integer(K)
  if (length(K) != 1L || is.na(K) || K < 2L) stop("K must be an integer >= 2.")
  if (length(ordering_vec) != nrow(X)) stop("ordering_vec must have length nrow(X).")

  keep_idx <- seq_len(nrow(X))
  if (anyNA(ordering_vec)) {
    if (na_action == "error") stop("ordering_vec contains NA.")
    keep_idx <- which(!is.na(ordering_vec))
    ordering_vec <- ordering_vec[keep_idx]
    X <- X[keep_idx, , drop = FALSE]
  }

  n <- nrow(X)
  d <- ncol(X)
  if (n < K) stop("After NA handling, n must be >= K.")

  # ---- discretize ordering_vec into K bins / clusters ----
  if (discretization == "kmeans") {
    cl <- stats::kmeans(ordering_vec, centers = K)
    centers <- as.numeric(cl$centers)
    cluster_ids <- cl$cluster

  } else if (discretization == "quantile") {
    probs <- seq(0, 1, length.out = K + 1)
    brks <- as.numeric(stats::quantile(ordering_vec, probs = probs, na.rm = TRUE))
    brks <- unique(brks)
    if (length(brks) < 2) stop("ordering_vec has too few unique values for quantile discretization.")

    cluster_ids <- cut(ordering_vec, breaks = brks, include.lowest = TRUE, labels = FALSE)

    if (all(is.na(cluster_ids))) stop("Quantile discretization failed (all clusters NA).")
    if (anyNA(cluster_ids)) {
      warning("Some values were not assigned to any quantile bin; dropping them.")
      ok <- which(!is.na(cluster_ids))
      cluster_ids <- cluster_ids[ok]
      ordering_vec <- ordering_vec[ok]
      X <- X[ok, , drop = FALSE]
      keep_idx <- keep_idx[ok]
      n <- nrow(X)
    }

    K_eff <- max(cluster_ids, na.rm = TRUE)
    if (K_eff != K) {
      warning("Quantile discretization produced ", K_eff, " groups (requested K=", K, ").")
      K <- K_eff
    }

    centers <- tapply(ordering_vec, cluster_ids, mean)
    centers <- as.numeric(centers)

  } else { # "equal"
    brks <- seq(min(ordering_vec), max(ordering_vec), length.out = K + 1)
    cluster_ids <- cut(ordering_vec, breaks = brks, include.lowest = TRUE, labels = FALSE)
    if (all(is.na(cluster_ids))) stop("Equal discretization failed (all clusters NA).")

    centers <- tapply(ordering_vec, cluster_ids, mean)
    centers <- as.numeric(centers)
  }

  # reorder clusters from left to right by centers
  cluster_order <- order(centers)
  cluster_rank <- match(cluster_ids, cluster_order)

  # ---- mu ----
  mu_list <- lapply(seq_len(K), function(k) {
    idx <- which(cluster_rank == k)
    if (length(idx) == 0) rep(0, d) else colMeans(X[idx, , drop = FALSE])
  })

  # # ---- pi ----
  Nk <- tabulate(cluster_rank, nbins = K)
  pi_vec <- Nk / n

  # ---- sigma2 ----
  if (modelName == "homoskedastic") {
    # pooled within-cluster residuals, shared across clusters
    diffs <- do.call(rbind, lapply(seq_len(K), function(k) {
      idx <- which(cluster_rank == k)
      if (length(idx) == 0) return(NULL)
      sweep(X[idx, , drop = FALSE], 2, mu_list[[k]])
    }))

    if (is.null(diffs) || nrow(diffs) == 0) {
      sigma2 <- rep(1, d)
    } else {
      sigma2 <- colMeans(diffs^2) + nugget
    }
    sigma2 <- pmax(as.numeric(sigma2), eps)

  } else {
    # heteroskedastic: per (j,k) within-cluster variance
    sigma2 <- matrix(NA_real_, nrow = d, ncol = K)
    global_var <- apply(X, 2, stats::var)
    global_var[!is.finite(global_var)] <- 1
    global_var <- pmax(global_var, eps)

    for (k in seq_len(K)) {
      idx <- which(cluster_rank == k)
      if (length(idx) <= 1L) {
        sigma2[, k] <- global_var + nugget
      } else {
        dif <- sweep(X[idx, , drop = FALSE], 2, mu_list[[k]])
        sigma2[, k] <- colMeans(dif^2) + nugget
      }
    }
    sigma2 <- pmax(sigma2, eps)
  }

  list(
    pi = pi_vec,
    mu = mu_list,
    sigma2 = sigma2,
    keep_idx = keep_idx,
    cluster_rank = cluster_rank
  )
}




#' @keywords internal
.convert_initparams_to_csmooth <- function(init_params,
                                           modelName = c("homoskedastic", "heteroskedastic"),
                                           eps = 1e-12) {
  modelName <- match.arg(modelName)

  pi <- as.numeric(init_params$pi)
  mu <- init_params$mu
  sigma_list <- init_params$sigma
  K <- length(pi)
  d <- length(mu[[1]])

  # extract diag variances (d x K)
  sig2_mat <- vapply(seq_len(K), function(k) diag(as.matrix(sigma_list[[k]])), numeric(d))
  sig2_mat <- matrix(sig2_mat, nrow = d, ncol = K)
  # vapply gives d x K already

  if (modelName == "homoskedastic") {
    sigma2 <- rowMeans(sig2_mat)
  } else {
    sigma2 <- sig2_mat
  }

  pi <- pmax(pi, eps)
  pi <- pi / sum(pi)

  sigma2 <- pmax(sigma2, eps)

  list(pi = pi, mu = mu, sigma2 = sigma2)
}





#' Initialize ordering for csmoothEM (diagonal-variance version)
#'
#' @description
#' Same ordering logic as \code{initialize_ordering()}, but returns csmooth-style
#' parameters with diagonal variances \code{sigma2}.
#'
#' @param X Numeric matrix (n x d).
#' @param K Integer >= 2; number of mixture components.
#' @param method One of \code{"PCA","fiedler","pcurve","tSNE","random"}.
#' @param discretization One of \code{"equal","quantile","kmeans"}.
#' @param modelName Either \code{"homoskedastic"} or \code{"heteroskedastic"}.
#' @param nugget Nonnegative scalar added to variance estimates.
#' @param eps Small positive floor for \code{pi} and \code{sigma2}.
#' @param ... Extra arguments passed to the ordering method (PCA/tSNE/pcurve/fiedler).
#'
#' @return A list with fields:
#' \itemize{
#'   \item \code{params}: list(pi, mu, sigma2)
#'   \item \code{keep_idx}: indices kept after NA handling
#'   \item \code{ordering}: ordering metadata including score \code{t} and method name
#' }
#'
#' @export
initialize_ordering_csmooth <- function(
    X, K,
    method = c("PCA", "fiedler", "pcurve", "tSNE", "random"),
    discretization = c("equal", "quantile", "kmeans"),
    modelName = c("homoskedastic", "heteroskedastic"),
    nugget = 0,
    eps = 1e-12,
    ...
) {
  method <- match.arg(method)
  discretization <- match.arg(discretization)
  modelName <- match.arg(modelName)

  X <- as.matrix(X)
  if (nrow(X) < 1L) stop("X must have at least one row.")
  K <- as.integer(K)
  if (length(K) != 1L || is.na(K) || K < 2L) stop("K must be a single integer >= 2.")

  nugget <- as.numeric(nugget)
  if (!is.finite(nugget) || nugget < 0) stop("nugget must be nonnegative.")
  eps <- as.numeric(eps)
  if (!is.finite(eps) || eps <= 0) stop("eps must be positive.")

  if (method == "random") {
    # IMPORTANT: keep your existing random logic
    init <- make_default_init(X = X, K = K, ordering = TRUE)
    init$keep_idx <- seq_len(nrow(X))
    init$ordering <- list(method = method, t = NA_real_, keep_idx = init$keep_idx)

    params_cs <- .convert_initparams_to_csmooth(
      init_params = init,
      modelName = modelName,
      eps = eps
    )

    # nugget applies after conversion
    if (modelName == "homoskedastic") {
      params_cs$sigma2 <- pmax(as.numeric(params_cs$sigma2) + nugget, eps)
    } else {
      params_cs$sigma2 <- pmax(as.matrix(params_cs$sigma2) + nugget, eps)
    }

    return(list(
      params = params_cs,
      keep_idx = init$keep_idx,
      ordering = init$ordering
    ))
  }

  # non-random: ordering score + make_init_csmooth()
  ordering_result <- switch(
    method,
    PCA     = PCA_ordering(X, ...),
    tSNE    = tSNE_ordering(X, ...),
    pcurve  = pcurve_ordering(X, ...),
    fiedler = fiedler_ordering(X, ...)
  )

  init <- make_init_csmooth(
    X = X,
    ordering_vec = ordering_result$t,
    K = K,
    modelName = modelName,
    nugget = nugget,
    discretization = discretization,
    na_action = "drop",
    eps = eps
  )

  ordering <- c(ordering_result, list(method = method))
  list(
    params = list(pi = init$pi, mu = init$mu, sigma2 = init$sigma2),
    keep_idx = init$keep_idx,
    ordering = ordering
  )
}



#' Initialize csmoothEM (coordinate-specific SmoothEM)
#'
#' @description
#' Create a \code{csmooth_em} object using the same ordering-initialization ideas as
#' \code{initialize_smoothEM()}, but with coordinate-specific penalties \code{lambda_vec}
#' and diagonal covariance structures only.
#'
#' Workflow:
#' \enumerate{
#'   \item Initialize an ordering and discretize into \code{K} components using
#'   \code{initialize_ordering_csmooth()}.
#'   \item Build the base random-walk precision \code{Q_K} along components (with \code{lambda=1}).
#'   \item Initialize \code{lambda_vec}. If \code{adaptive} is not \code{"none"}, estimate an initial
#'   \code{lambda_vec} from the initialized parameters (and clamp to \code{[lambda_min, lambda_max]}).
#'   \item Construct a \code{csmooth_em} object via \code{as_csmooth_em()}, then run a warm start via
#'   \code{do_csmoothEM()} for \code{num_iter} iterations.
#' }
#'
#' @param X n-by-d numeric matrix.
#' @param method One of \code{"tSNE","PCA","random","fiedler","multi_scale"}.
#'   Currently \code{"multi_scale"} is not implemented for csmoothEM and will error.
#' @param rw_q Integer RW order along K for \code{Q_K}.
#' @param lambda Scalar or length-d vector. If scalar, recycled to length d.
#'   Used as an initial value when \code{adaptive="none"}.
#' @param relative_lambda Logical; if TRUE, scale the base prior for coordinate \code{j} by a
#'   coordinate-specific variance scale (see Details).
#' @param K Number of mixture components. If NULL, defaults to \code{min(50, floor(n/5))} (at least 2).
#' @param num_iter Integer >= 1; number of warm-start iterations to run immediately.
#' @param modelName Either \code{"homoskedastic"} or \code{"heteroskedastic"}.
#' @param ridge Nonnegative ridge added when building RW precision \code{Q_K}.
#' @param nugget Nonnegative nugget used in M-step variance updates.
#' @param eigen_tol Optional eigen tolerance used for generalized log-determinants (if needed downstream).
#' @param include.data Logical; store the (kept) data matrix in the returned object.
#' @param adaptive Character specifying how to initialize/update \code{lambda_vec}:
#'   \describe{
#'     \item{\code{"none"}}{Do not estimate \code{lambda_vec} from the initialization; use \code{lambda}.}
#'     \item{\code{"prior"}}{Prior-based estimate using the initialized component means \eqn{\mu}.}
#'     \item{\code{"ml"}}{Collapsed-ML warm start (dispatches to \code{do_csmoothEM_ml_collapsed()} inside
#'       \code{do_csmoothEM()}), allowing optional ML updates of \eqn{\sigma^2}.}
#'   }
#'   Logical values are accepted for backward compatibility: \code{TRUE} is treated as \code{"prior"}
#'   and \code{FALSE} as \code{"none"}.
#' @param lambda_min,lambda_max Positive bounds for \code{lambda_vec} (used when \code{adaptive!="none"}).
#' @param sigma_update Character. Only used when \code{adaptive="ml"} during the warm start.
#'   Either \code{"mstep"} or \code{"ml"}; see \code{\link{do_csmoothEM_ml_collapsed}}.
#' @param sigma_min,sigma_max Positive bounds for \code{sigma2} when \code{adaptive="ml"} and
#'   \code{sigma_update="ml"}.
#' @param discretization Discretization method passed to \code{initialize_ordering_csmooth()}.
#' @param ... Passed to the ordering method (e.g. PCA/tSNE/pcurve/fiedler).
#'
#' @details
#' Let \eqn{Q_K} denote the RW(q) precision along components. For a separable prior
#' \eqn{Q_{\mathrm{full}} = I_d \otimes Q_K}, the rank deficiency is \eqn{rw_q} per coordinate,
#' so a convenient degrees-of-freedom term is \eqn{r = K - rw_q}.
#'
#' When \code{adaptive = "prior"}, the initializer estimates each coordinate penalty by:
#' \deqn{\lambda_j = r \big/ \left( \mu_{j\cdot}^\top Q_{j,\mathrm{base}} \mu_{j\cdot} \right),}
#' where \eqn{\mu_{j\cdot}} is the length-K vector of component means for coordinate \eqn{j}, and
#' \eqn{Q_{j,\mathrm{base}}} is the coordinate-specific base precision returned by
#' \code{.compute_Qbase_j()} (equivalently, \eqn{Q_K} possibly rescaled when
#' \code{relative_lambda = TRUE}). The estimate is clamped to \eqn{[\text{lambda_min},\text{lambda_max}]}.
#'
#' When \code{adaptive = "ml"}, the warm start uses the collapsed-ML routine. In this mode,
#' \code{do_csmoothEM()} dispatches to \code{do_csmoothEM_ml_collapsed()}, which optimizes a
#' collapsed (Laplace-exact, Gaussian) objective and records \code{ml_trace} as the collapsed
#' objective \eqn{\mathcal{C}}.
#'
#' When \code{relative_lambda = TRUE}, the base precision \eqn{Q_{j,\mathrm{base}}} is scaled by a
#' coordinate variance proxy:
#' \itemize{
#'   \item \code{modelName = "homoskedastic"}: scale by \code{1 / sigma2[j]}.
#'   \item \code{modelName = "heteroskedastic"}: scale by \code{1 / sum_k pi_k sigma2[j,k]}.
#' }
#'
#' @return A \code{csmooth_em} object.
#' @export
initialize_csmoothEM <- function(
    X,
    method = c("tSNE", "PCA", "pcurve", "random", "fiedler", "multi_scale"),
    rw_q = 2,
    lambda = 1,
    relative_lambda = TRUE,
    K = NULL,
    num_iter = 1,
    modelName = c("homoskedastic", "heteroskedastic"),
    ridge = 0,
    nugget = 0,
    eigen_tol = NULL,
    include.data = TRUE,
    adaptive = "ml",          # "none" / "prior" / "ml" (logical also accepted)
    lambda_min = 1e-8,
    lambda_max = 1e8,
    sigma_update = c("ml", "mstep"),  # only used when adaptive="ml"
    sigma_min = 1e-10,
    sigma_max = 1e10,
    discretization = c("equal", "quantile", "kmeans"),
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)

  method <- match.arg(method)
  modelName <- match.arg(modelName)
  discretization <- match.arg(discretization)

  if (method == "multi_scale") {
    stop("initialize_csmoothEM: method='multi_scale' is not implemented yet for csmoothEM.")
  }

  rw_q <- as.integer(rw_q)
  if (length(rw_q) != 1L || is.na(rw_q) || rw_q < 0L) stop("rw_q must be a single integer >= 0.")

  num_iter <- as.integer(num_iter)
  if (length(num_iter) != 1L || is.na(num_iter) || num_iter < 1L) stop("num_iter must be integer >= 1.")

  lambda_min <- as.numeric(lambda_min)
  lambda_max <- as.numeric(lambda_max)
  if (!is.finite(lambda_min) || !is.finite(lambda_max) || lambda_min <= 0 || lambda_max <= 0 || lambda_min > lambda_max) {
    stop("lambda_min/lambda_max must be positive finite numbers with lambda_min <= lambda_max.")
  }

  sigma_update <- match.arg(sigma_update)
  sigma_min <- as.numeric(sigma_min)
  sigma_max <- as.numeric(sigma_max)
  if (!is.finite(sigma_min) || !is.finite(sigma_max) || sigma_min <= 0 || sigma_max <= 0 || sigma_min > sigma_max) {
    stop("sigma_min/sigma_max must be positive finite numbers with sigma_min <= sigma_max.")
  }

  # --- normalize adaptive option (match do_csmoothEM semantics)
  if (is.logical(adaptive)) adaptive <- if (isTRUE(adaptive)) "ml" else "none"
  adaptive <- match.arg(adaptive, choices = c("none", "prior", "ml"))

  # choose K if not provided
  if (is.null(K)) {
    K <- min(50L, floor(n / 5))
    K <- max(2L, as.integer(K))
  } else {
    K <- as.integer(K)
    if (length(K) != 1L || is.na(K) || K < 2L) stop("K must be a single integer >= 2.")
  }

  # ------------------------------------------------------------
  # 1) ordering init (csmooth-style: returns sigma2 directly)
  # ------------------------------------------------------------
  init <- initialize_ordering_csmooth(
    X = X, K = K,
    method = method,
    discretization = discretization,
    modelName = modelName,
    ...
  )

  keep_idx <- init$keep_idx %||% seq_len(nrow(X))
  X_use <- X  # currently not subsetting (keep_idx retained for metadata)

  params0 <- init$params
  if (is.null(params0$pi) || is.null(params0$mu) || is.null(params0$sigma2)) {
    stop("initialize_ordering_csmooth must return $params with fields pi, mu, sigma2.")
  }

  K_eff <- length(params0$pi)
  if (K_eff < 2L) stop("Initialization produced K<2 components; cannot proceed.")

  # ------------------------------------------------------------
  # 2) build Q_K base precision (lambda=1)
  # ------------------------------------------------------------
  Q_K <- make_random_walk_precision(K = K_eff, d = 1, lambda = 1, q = rw_q, ridge = ridge)
  if (nrow(Q_K) != K_eff || ncol(Q_K) != K_eff) {
    stop("make_random_walk_precision(K, d=1, ...) must return a K x K matrix/Matrix.")
  }

  # ------------------------------------------------------------
  # 3) initialize lambda_vec: user-specified vs estimated
  # ------------------------------------------------------------
  if (length(lambda) == 1L) {
    lambda_vec <- rep(as.numeric(lambda), d)
  } else {
    lambda_vec <- as.numeric(lambda)
    if (length(lambda_vec) != d) stop("lambda must be scalar or length d.")
  }

  clamp_lambda <- function(lam) {
    lam <- as.numeric(lam)
    if (!is.finite(lam)) return(NA_real_)
    min(max(lam, lambda_min), lambda_max)
  }

  compute_Qbase_j <- function(params, Q_K, relative_lambda, modelName, j, eps = 1e-12) {
    if (!isTRUE(relative_lambda)) return(Q_K)
    if (modelName == "homoskedastic") {
      s2j <- as.numeric(params$sigma2[j])
    } else {
      s2j <- sum(as.numeric(params$pi) * as.numeric(params$sigma2[j, ]))
    }
    s2j <- pmax(s2j, eps)
    Q_K / s2j
  }

  estimate_lambda_vec_prior <- function(params, Q_K, relative_lambda, modelName, rw_q,
                                        eps_quad = 1e-12) {
    K <- length(params$pi)
    r <- max(as.integer(K - rw_q), 1L)
    mu_mat <- do.call(cbind, params$mu)  # d x K
    lam_hat <- rep(NA_real_, d)
    for (j in seq_len(d)) {
      Qb <- compute_Qbase_j(params, Q_K, relative_lambda, modelName, j)
      mj <- as.numeric(mu_mat[j, ])
      quad <- as.numeric(crossprod(mj, Qb %*% mj))
      lam_hat[j] <- clamp_lambda(r / pmax(quad, eps_quad))
    }
    lam_hat
  }

  # Note: for adaptive="ml" we do not need a special pre-estimate of lambda;
  # the warm-start routine will update lambda anyway. We keep a lightweight initialization:
  lambda_init_est <- rep(NA_real_, d)
  lambda_init_method <- adaptive

  if (adaptive == "prior") {
    lambda_init_est <- estimate_lambda_vec_prior(
      params = params0, Q_K = Q_K,
      relative_lambda = relative_lambda,
      modelName = modelName,
      rw_q = rw_q
    )
    ok <- is.finite(lambda_init_est)
    if (any(ok)) lambda_vec <- lambda_init_est
  }

  # ------------------------------------------------------------
  # 4) build csmooth_em object
  # ------------------------------------------------------------
  obj <- as_csmooth_em(
    params = params0,
    gamma = NULL,
    data  = if (include.data) X_use else NULL,
    Q_K   = Q_K,
    lambda_vec = lambda_vec,
    rw_q  = rw_q,
    ridge = ridge,
    modelName = modelName,
    relative_lambda = relative_lambda,
    nugget = nugget,
    eigen_tol = eigen_tol,
    meta = list(
      init = list(
        method = method,
        details = list(
          K = K_eff,
          keep_idx = keep_idx,
          lambda_init_method = lambda_init_method,
          lambda_init_est = lambda_init_est
        )
      )
    )
  )

  obj$control <- obj$control %||% list()
  obj$control$adaptive <- adaptive
  obj$control$lambda_min <- lambda_min
  obj$control$lambda_max <- lambda_max

  obj$control$sigma_update <- sigma_update
  obj$control$sigma_min <- sigma_min
  obj$control$sigma_max <- sigma_max

  # ------------------------------------------------------------
  # 5) warm start: run num_iter iterations
  # ------------------------------------------------------------
  obj <- do_csmoothEM(
    object = obj,
    data = X_use,
    iter = num_iter,
    record = TRUE,
    adaptive = adaptive,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_update = sigma_update,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    verbose = FALSE
  )

  obj
}





#' Plot a csmooth_em object
#'
#' @description
#' Visualization for a \code{csmooth_em} object.
#' \itemize{
#'   \item \code{plot_type="scatterplot"}: calls \code{plot_EM_embedding()} (same behavior as \code{smooth_em}).
#'   \item \code{plot_type="elbo"}: plots ELBO and penalized observed-data objective traces.
#'   \item \code{plot_type="mu"}: plots only the component means (with arrows for 2D).
#' }
#'
#' @param x A \code{csmooth_em} object.
#' @param data Numeric matrix (n x d). Required when \code{plot_type="scatterplot"} unless stored in \code{x$data}.
#' @param plot_type One of \code{"scatterplot"}, \code{"elbo"}, \code{"mu"}.
#' @param dims Integer vector of length 1 or 2. Used for \code{"scatterplot"} and \code{"mu"}.
#' @param two_panel Logical; if TRUE and \code{plot_type="elbo"}, draw ELBO and objective in two panels.
#' @param verbose Logical; not passed to graphics (prevents "not a graphical parameter" warnings).
#' @param ... Passed to the underlying plotting functions.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.csmooth_em <- function(
    x,
    data = NULL,
    plot_type = c("scatterplot", "elbo", "mu"),
    dims = c(1, 2),
    two_panel = FALSE,
    verbose = FALSE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!inherits(x, "csmooth_em")) stop("x must be a 'csmooth_em' object.")
  plot_type <- match.arg(plot_type)

  # prevent accidental verbose passing to base graphics
  dots <- list(...)
  if ("verbose" %in% names(dots)) {
    warning("Argument 'verbose' was provided in ...; ignoring it for plotting.")
    dots$verbose <- NULL
  }

  dims <- as.integer(dims)
  if (!(length(dims) %in% c(1L, 2L))) stop("dims must have length 1 or 2.")

  # -------------------------
  # scatterplot: reuse your existing helper
  # -------------------------
  if (plot_type == "scatterplot") {
    if (is.null(data)) {
      data <- x$data
      if (is.null(data)) stop("data must be provided either as an argument or stored in the object.")
    }
    data <- as.matrix(data)

    if (any(is.na(dims)) || any(dims < 1L) || any(dims > ncol(data))) stop("dims out of range.")

    do.call(plot_EM_embedding, c(list(fit = x, X = data, dims = dims), dots))
    return(invisible(x))
  }

  # -------------------------
  # elbo: match smooth_em style (optionally two-panel)
  # -------------------------
  if (plot_type == "elbo") {
    elbo <- x$elbo_trace %||% numeric(0)
    ll   <- x$loglik_trace %||% numeric(0)
    ml   <- x$ml_trace %||% numeric(0)

    if (length(elbo) == 0L && length(ll) == 0L && length(ml) == 0L) {
      warning("No traces found: elbo_trace/loglik_trace/ml_trace are empty.")
      return(invisible(x))
    }

    if (two_panel) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar), add = TRUE)
      par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

      # Panel 1: ELBO (and ML if available)
      if (length(elbo) > 0L || length(ml) > 0L) {
        xlim <- c(1, max(length(elbo), length(ml)))
        yall <- c(elbo, ml)
        yall <- yall[is.finite(yall)]
        ylim <- if (length(yall)) range(yall) else c(-1, 1)

        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
        axis(1); axis(2)
        box()
        title(main = "Penalized ELBO / collapsed objective traces",
              xlab = "Iteration", ylab = "Value")

        if (length(elbo) > 0L) lines(seq_along(elbo), elbo, lty = 1, lwd = 2, ...)
        if (length(ml)   > 0L) lines(seq_along(ml),   ml,   lty = 3, lwd = 2)

        legend("bottomright",
               legend = c(if (length(elbo) > 0L) "ELBO" else NULL,
                          if (length(ml)   > 0L) "ML (collapsed C)" else NULL),
               lty = c(if (length(elbo) > 0L) 1 else NULL,
                       if (length(ml)   > 0L) 3 else NULL),
               lwd = 2, bty = "n")
      } else {
        plot.new(); title("Penalized ELBO / collapsed objective traces (empty)")
      }

      # Panel 2: penalized observed objective
      if (length(ll) > 0L) {
        plot(seq_along(ll), ll, type = "l",
             xlab = "Iteration", ylab = "Penalized log-likelihood",
             main = "Penalized observed-data objective trace", ...)
      } else {
        plot.new(); title("Penalized observed-data objective trace (empty)")
      }

    } else {
      # Overlaid in one panel: ELBO (solid), objective (dashed), ML (dotted)
      xlim <- c(1, max(length(elbo), length(ll), length(ml)))
      yall <- c(elbo, ll, ml)
      yall <- yall[is.finite(yall)]
      ylim <- if (length(yall)) range(yall) else c(-1, 1)

      plot.new()
      plot.window(xlim = xlim, ylim = ylim)
      axis(1); axis(2)
      box()
      title(main = "Traces (overlaid)", xlab = "Iteration", ylab = "Value")

      if (length(elbo) > 0L) lines(seq_along(elbo), elbo, lty = 1, lwd = 2, ...)
      if (length(ll)   > 0L) lines(seq_along(ll),   ll,   lty = 2, lwd = 2)
      if (length(ml)   > 0L) lines(seq_along(ml),   ml,   lty = 3, lwd = 2)

      legend("bottomright",
             legend = c(if (length(elbo) > 0L) "ELBO" else NULL,
                        if (length(ll)   > 0L) "penalized objective" else NULL,
                        if (length(ml)   > 0L) "ML (collapsed C)" else NULL),
             lty = c(if (length(elbo) > 0L) 1 else NULL,
                     if (length(ll)   > 0L) 2 else NULL,
                     if (length(ml)   > 0L) 3 else NULL),
             lwd = 2, bty = "n")
    }

    return(invisible(x))
  }

  # -------------------------
  # mu only (no data)
  # -------------------------
  if (is.null(x$params) || is.null(x$params$mu)) stop("x$params$mu is missing.")
  mu_mat_full <- do.call(rbind, x$params$mu)  # K x d
  K <- nrow(mu_mat_full)
  d <- ncol(mu_mat_full)

  if (any(is.na(dims)) || any(dims < 1L) || any(dims > d)) stop("dims out of range for mu.")

  mup <- mu_mat_full[, dims, drop = FALSE]

  if (length(dims) == 2L) {
    plot(mup[, 1], mup[, 2], type = "n",
         xlab = sprintf("mu[, %d]", dims[1]),
         ylab = sprintf("mu[, %d]", dims[2]),
         main = sprintf("csmoothEM: component means on dims (%d, %d)", dims[1], dims[2]),
         ...)
    for (k in seq_len(K - 1L)) {
      dx <- mup[k + 1, 1] - mup[k, 1]
      dy <- mup[k + 1, 2] - mup[k, 2]
      if (sqrt(dx^2 + dy^2) > 1e-12) {
        arrows(mup[k, 1], mup[k, 2], mup[k + 1, 1], mup[k + 1, 2],
               lwd = 3, length = 0.08)
      }
    }
    points(mup[, 1], mup[, 2], pch = 8, cex = 1.0, lwd = 1)
  } else {
    j <- dims[1]
    plot(seq_len(K), mup[, 1], type = "l",
         xlab = "k (component)", ylab = sprintf("mu[, %d]", j),
         main = sprintf("csmoothEM: mean curve on dim %d", j),
         lwd = 2, ...)
    points(seq_len(K), mup[, 1], pch = 8, cex = 1.0)
  }

  invisible(x)
}



#' Print csmooth_em object
#'
#' @param x A \code{csmooth_em} object.
#' @param ... Unused.
#'
#' @export
print.csmooth_em <- function(x, ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!inherits(x, "csmooth_em")) stop("x must be a 'csmooth_em' object.")

  K <- length(x$params$pi %||% numeric(0))
  d <- if (!is.null(x$params$mu) && length(x$params$mu) > 0) length(x$params$mu[[1]]) else NA_integer_
  n <- if (!is.null(x$data)) nrow(x$data) else NA_integer_

  modelName <- x$control$modelName %||% NA_character_
  it <- x$iter %||% length(x$elbo_trace %||% numeric(0))

  init_method <- x$meta$init$method %||% NA_character_
  adaptive_mode <- x$control$adaptive %||% NA_character_

  rw_q <- x$prior$rw_q %||% NA_integer_
  prior_desc <- sprintf("RW(%s) prior along K = %s",
                        ifelse(is.na(rw_q), "?", as.integer(rw_q)),
                        ifelse(is.na(K), "NA", as.integer(K)))

  lam <- x$prior$lambda_vec
  lam_line <- "lambda_vec: NA"
  if (!is.null(lam) && length(lam) > 0) {
    lam_line <- sprintf("lambda_vec: range=[%.3g, %.3g], mean=%.3g, relative = %s",
                        min(lam), max(lam), mean(lam),
                        ifelse(isTRUE(x$control$relative_lambda), "TRUE", "FALSE"))
  }

  get_last <- function(v) {
    v <- as.numeric(v %||% numeric(0))
    if (length(v) == 0L) return(NA_real_)
    tail(v, 1)
  }

  elbo_last <- get_last(x$elbo_trace)
  obj_last  <- get_last(x$loglik_trace)
  ml_last   <- get_last(x$ml_trace)

  cat(sprintf("Fitted csmoothEM object with %s\n", prior_desc))
  cat("-----\n")
  cat(sprintf("  n = %s, d = %s, modelName = %s\n",
              ifelse(is.na(n), "NA", n),
              ifelse(is.na(d), "NA", d),
              modelName))
  cat(sprintf("  iter = %s; init_method = %s; adaptive = %s;\n",
              as.integer(it), init_method, adaptive_mode))
  cat(sprintf("  %s\n", lam_line))

  cat(sprintf("  ELBO last = %s; penLogLik last = %s; ML/C last = %s\n",
              ifelse(is.finite(elbo_last), sprintf("%.6f", elbo_last), "NA"),
              ifelse(is.finite(obj_last),  sprintf("%.6f", obj_last),  "NA"),
              ifelse(is.finite(ml_last),   sprintf("%.6f", ml_last),   "NA")))

  invisible(x)
}




#' Summary for csmooth_em object
#'
#' @description
#' Returns a detailed summary of a fitted \code{csmooth_em} object, including
#' dimensions, initialization/adaptive settings, parameter ranges, last trace
#' values, and last-step relative changes.
#'
#' @param object A \code{csmooth_em} object.
#' @param ... Unused.
#'
#' @return A list with class \code{"summary.csmooth_em"}.
#' @export
summary.csmooth_em <- function(object, ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!inherits(object, "csmooth_em")) stop("object must be a 'csmooth_em' object.")

  K <- length(object$params$pi %||% numeric(0))
  d <- if (!is.null(object$params$mu) && length(object$params$mu) > 0) length(object$params$mu[[1]]) else NA_integer_
  n <- if (!is.null(object$data)) nrow(object$data) else NA_integer_
  it <- object$iter %||% length(object$elbo_trace %||% numeric(0))

  init_method <- object$meta$init$method %||% NA_character_
  adaptive_mode <- object$control$adaptive %||% NA_character_
  sigma_update <- object$control$sigma_update %||% NA_character_

  # pi
  pi <- object$params$pi %||% rep(NA_real_, K)
  pi_min <- if (K > 0) min(pi) else NA_real_
  pi_max <- if (K > 0) max(pi) else NA_real_

  # sigma2
  sig2 <- object$params$sigma2
  sig2_rng <- c(NA_real_, NA_real_)
  sig2_type <- NA_character_
  if (!is.null(sig2)) {
    sig2_type <- if (is.matrix(sig2)) "matrix(d×K)" else "vector(d)"
    sig2_rng <- c(min(sig2), max(sig2))
  }

  # lambda
  lam <- object$prior$lambda_vec
  lam_rng <- c(NA_real_, NA_real_)
  lam_mean <- NA_real_
  if (!is.null(lam) && length(lam) > 0L) {
    lam_rng <- c(min(lam), max(lam))
    lam_mean <- mean(lam)
  }

  get_last <- function(v) {
    v <- as.numeric(v %||% numeric(0))
    if (length(v) == 0L) return(NA_real_)
    tail(v, 1)
  }

  rel_change <- function(v) {
    v <- as.numeric(v %||% numeric(0))
    if (length(v) < 2) return(NA_real_)
    a <- v[length(v) - 1]
    b <- v[length(v)]
    if (!is.finite(a) || !is.finite(b)) return(NA_real_)
    abs(b - a) / max(1e-12, abs(a))
  }

  rel_change_vec <- function(v) {
    if (is.null(v) || length(v) < 2) return(NA_real_)
    a <- as.numeric(v[[length(v) - 1]])
    b <- as.numeric(v[[length(v)]])
    if (length(a) != length(b) || length(a) == 0) return(NA_real_)
    sqrt(sum((b - a)^2)) / max(1e-12, sqrt(sum(a^2)))
  }

  out <- list(
    # core dimensions
    n = n, d = d, K = K,
    iter = as.integer(it),

    # configuration
    modelName = object$control$modelName %||% NA_character_,
    relative_lambda = isTRUE(object$control$relative_lambda),
    init_method = init_method,
    adaptive = adaptive_mode,
    sigma_update = sigma_update,

    # parameter ranges
    pi_min = pi_min,
    pi_max = pi_max,
    sigma2_type = sig2_type,
    sigma2_min = sig2_rng[1],
    sigma2_max = sig2_rng[2],
    lambda_min = lam_rng[1],
    lambda_max = lam_rng[2],
    lambda_mean = lam_mean,

    # last values
    elbo_last = get_last(object$elbo_trace),
    obj_last  = get_last(object$loglik_trace),
    ml_last   = get_last(object$ml_trace),
    logdetH_last = get_last(object$logdetH_trace),

    # last-step relative changes
    elbo_rel_change = rel_change(object$elbo_trace),
    obj_rel_change  = rel_change(object$loglik_trace),
    ml_rel_change   = rel_change(object$ml_trace),
    lambda_rel_change = rel_change_vec(object$lambda_trace)
  )

  class(out) <- "summary.csmooth_em"
  out
}

#' Print summary of csmooth_em
#'
#' @param x A \code{summary.csmooth_em} object.
#' @param ... Unused.
#'
#' @export
print.summary.csmooth_em <- function(x, ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  cat("<summary.csmooth_em>\n")
  cat(sprintf("  n=%s, d=%s, K=%s, iter=%d\n",
              ifelse(is.na(x$n), "NA", x$n),
              ifelse(is.na(x$d), "NA", x$d),
              ifelse(is.na(x$K), "NA", x$K),
              x$iter))
  cat(sprintf("  modelName=%s, relative_lambda=%s\n",
              x$modelName, ifelse(isTRUE(x$relative_lambda), "TRUE", "FALSE")))

  cat(sprintf("  init_method=%s, adaptive=%s\n",
              x$init_method %||% "NA", x$adaptive %||% "NA"))
  if (!is.na(x$sigma_update) && nzchar(x$sigma_update)) {
    cat(sprintf("  sigma_update (ml path)=%s\n", x$sigma_update))
  }

  cat(sprintf("  pi range: [%.3g, %.3g]\n", x$pi_min, x$pi_max))
  cat(sprintf("  sigma2 (%s) range: [%.3g, %.3g]\n",
              x$sigma2_type %||% "NA", x$sigma2_min, x$sigma2_max))
  cat(sprintf("  lambda range: [%.3g, %.3g], mean=%.3g\n",
              x$lambda_min, x$lambda_max, x$lambda_mean))

  cat(sprintf("  ELBO last: %.6f\n", x$elbo_last))
  cat(sprintf("  penLogLik last: %.6f\n", x$obj_last))
  if (!is.na(x$ml_last)) cat(sprintf("  ML/C last: %.6f\n", x$ml_last))
  if (!is.na(x$logdetH_last)) cat(sprintf("  log|H| last: %.6f\n", x$logdetH_last))

  any_rc <- any(is.finite(c(x$elbo_rel_change, x$obj_rel_change, x$ml_rel_change, x$lambda_rel_change)))
  if (any_rc) {
    cat("  relative change (last step):\n")
    if (is.finite(x$elbo_rel_change))   cat(sprintf("    ELBO: %.3g\n", x$elbo_rel_change))
    if (is.finite(x$obj_rel_change))    cat(sprintf("    penLogLik: %.3g\n", x$obj_rel_change))
    if (is.finite(x$ml_rel_change))     cat(sprintf("    ML/C: %.3g\n", x$ml_rel_change))
    if (is.finite(x$lambda_rel_change)) cat(sprintf("    lambda_vec (L2 rel): %.3g\n", x$lambda_rel_change))
  }

  invisible(x)
}



#' Run multiple csmoothEM initializations in parallel and summarize results
#'
#' @description
#' Runs \code{\link{initialize_csmoothEM}} with multiple initialization \code{methods},
#' each for \code{num_iter} warm-start iterations, and returns a named list of fits.
#' A per-method summary table is attached as an attribute \code{"summary"}.
#'
#' The summary table reports the last values of \code{elbo_trace}, \code{loglik_trace},
#' and \code{ml_trace}. In the current codebase, \code{ml_trace} is the collapsed
#' objective \eqn{\mathcal{C}} (Laplace-exact in the Gaussian csmooth setting).
#'
#' @param X Numeric matrix \code{(n x d)}.
#' @param methods Character vector of initialization methods. Default includes
#'   \code{"PCA"}, \code{"tSNE"}, \code{"random"}, \code{"fiedler"}, \code{"pcurve"}.
#' @param num_iter Integer \eqn{\ge 1}. Number of warm-start iterations to run inside
#'   \code{initialize_csmoothEM} for each method.
#' @param num_cores Integer \eqn{\ge 1}. Number of cores to use. On non-Windows systems,
#'   uses \code{parallel::mclapply}; on Windows, falls back to a PSOCK cluster.
#' @param K Optional integer \eqn{\ge 2}. Number of mixture components. If \code{NULL},
#'   the default logic in \code{initialize_csmoothEM} is used.
#' @param adaptive Character. One of \code{"none"}, \code{"prior"}, \code{"ml"}.
#'   Logical values are accepted for backward compatibility: \code{TRUE} is treated as
#'   \code{"prior"} and \code{FALSE} as \code{"none"}.
#' @param lambda_min,lambda_max Positive bounds for \code{lambda_vec} (passed to \code{initialize_csmoothEM}).
#' @param sigma_update Character. Only used when \code{adaptive="ml"} (passed to \code{initialize_csmoothEM}).
#'   Either \code{"mstep"} or \code{"ml"}; see \code{\link{do_csmoothEM_ml_collapsed}}.
#' @param sigma_min,sigma_max Positive bounds for \code{sigma2} when \code{adaptive="ml"} and
#'   \code{sigma_update="ml"}.
#' @param seed Optional integer seed. If provided, a different derived seed is used per method.
#' @param quiet Logical; reserved for future use.
#' @param ... Additional arguments passed to \code{initialize_csmoothEM} (and downstream ordering routines).
#'
#' @return A named list of fits, with names equal to \code{methods}. Each entry is either
#'   a \code{csmooth_em} object (on success) or \code{NULL} (on failure).
#'   The list has an attribute \code{"summary"} which is a data.frame with one row per method:
#'   \itemize{
#'     \item \code{method}: method name
#'     \item \code{success}: logical
#'     \item \code{n,d,K}: scalar integers from \code{summary(csmooth_em)} when available
#'     \item \code{elbo_last}: last penalized ELBO (or \code{NA})
#'     \item \code{obj_last}: last penalized observed objective (or \code{NA})
#'     \item \code{ml_last}: last collapsed objective \eqn{\mathcal{C}} (or \code{NA})
#'     \item \code{error}: error message if failed
#'   }
#'
#' @seealso \code{\link{optimize_initial_csmoothEM}}, \code{\link{initialize_csmoothEM}}
#' @export
parallel_initial_csmoothEM <- function(
    X,
    methods = c("PCA", "tSNE", "random", "fiedler", "pcurve"),
    num_iter = 1,
    num_cores = 2,
    K = NULL,
    adaptive = "prior",
    lambda_min = 1e-8,
    lambda_max = 1e8,
    sigma_update = c("mstep", "ml"),
    sigma_min = 1e-10,
    sigma_max = 1e10,
    seed = NULL,
    quiet = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  X <- as.matrix(X)
  methods <- unique(methods)

  # normalize adaptive option (match initialize semantics)
  if (is.null(adaptive)) adaptive <- "none"
  if (is.logical(adaptive)) adaptive <- if (isTRUE(adaptive)) "prior" else "none"
  adaptive <- match.arg(adaptive, choices = c("none", "prior", "ml"))

  sigma_update <- match.arg(sigma_update)

  if (!is.null(seed)) {
    set.seed(seed)
    seeds <- sample.int(1e9, length(methods))
  } else {
    seeds <- rep(NA_integer_, length(methods))
  }

  one_run <- function(i) {
    mm <- methods[i]
    if (!is.na(seeds[i])) set.seed(seeds[i])

    obj <- tryCatch(
      initialize_csmoothEM(
        X = X,
        method = mm,
        K = K,
        num_iter = num_iter,
        adaptive = adaptive,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        sigma_update = sigma_update,
        sigma_min = sigma_min,
        sigma_max = sigma_max,
        ...
      ),
      error = function(e) e
    )

    if (inherits(obj, "error")) {
      return(list(method = mm, fit = NULL, summary = NULL, error = conditionMessage(obj)))
    }
    if (!inherits(obj, "csmooth_em")) {
      return(list(method = mm, fit = NULL, summary = NULL,
                  error = sprintf("Returned object is not 'csmooth_em' (class=%s).",
                                  paste(class(obj), collapse = ","))))
    }

    # Prefer summary(obj) for n/d/K if available, but use traces directly for last values.
    s <- tryCatch(summary(obj), error = function(e) e)
    if (inherits(s, "error")) s <- NULL

    get_last <- function(v) {
      v <- as.numeric(v %||% numeric(0))
      if (length(v) == 0L) return(NA_real_)
      tail(v, 1)
    }

    elbo_last <- get_last(obj$elbo_trace)
    obj_last  <- get_last(obj$loglik_trace)
    ml_last   <- get_last(obj$ml_trace)

    # validate scalars
    if (!(is.na(elbo_last) || (length(elbo_last) == 1L && is.finite(elbo_last)))) elbo_last <- NA_real_
    if (!(is.na(obj_last)  || (length(obj_last)  == 1L && is.finite(obj_last))))  obj_last  <- NA_real_
    if (!(is.na(ml_last)   || (length(ml_last)   == 1L && is.finite(ml_last))))   ml_last   <- NA_real_

    list(method = mm, fit = obj, summary = s,
         elbo_last = elbo_last, obj_last = obj_last, ml_last = ml_last,
         error = NA_character_)
  }

  # run (prefer mclapply on mac/linux)
  if (num_cores > 1 && .Platform$OS.type != "windows") {
    res <- parallel::mclapply(seq_along(methods), one_run, mc.cores = num_cores)
  } else if (num_cores > 1 && .Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    res <- parallel::parLapply(cl, seq_along(methods), one_run)
  } else {
    res <- lapply(seq_along(methods), one_run)
  }

  fits <- setNames(vector("list", length(methods)), methods)
  sum_df <- data.frame(
    method = methods,
    success = FALSE,
    n = NA_integer_,
    d = NA_integer_,
    K = NA_integer_,
    elbo_last = NA_real_,
    obj_last  = NA_real_,
    ml_last   = NA_real_,
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(res)) {
    rr <- res[[i]]
    mm <- rr$method
    fits[[mm]] <- rr$fit

    sum_df$success[i] <- !is.null(rr$fit) && is.na(rr$error)
    sum_df$error[i] <- rr$error %||% NA_character_

    if (!is.null(rr$fit) && is.na(rr$error)) {
      # n/d/K from summary if available (best effort)
      if (!is.null(rr$summary)) {
        sum_df$n[i] <- as.integer(rr$summary$n %||% NA_integer_)
        sum_df$d[i] <- as.integer(rr$summary$d %||% NA_integer_)
        sum_df$K[i] <- as.integer(rr$summary$K %||% NA_integer_)
      }
      sum_df$elbo_last[i] <- rr$elbo_last
      sum_df$obj_last[i]  <- rr$obj_last
      sum_df$ml_last[i]   <- rr$ml_last
    }
  }

  attr(fits, "summary") <- sum_df
  fits
}



#' Optimize csmoothEM initialization by comparing multiple methods after a warm start
#'
#' @description
#' Runs several initialization methods (via \code{\link{parallel_initial_csmoothEM}}),
#' each for \code{num_iter} warm-start iterations, then selects the best fit.
#'
#' Selection criterion:
#' \itemize{
#'   \item If \code{adaptive == "ml"}, selects the fit with the largest \code{ml_last}
#'         (collapsed objective \eqn{\mathcal{C}}).
#'   \item Otherwise, selects the fit with the largest \code{elbo_last}.
#' }
#'
#' If \code{plot=TRUE}, plots traces across methods and highlights the selected best method.
#' When \code{adaptive=="ml"}, the primary plot shows \code{ml_trace}; otherwise it shows
#' \code{elbo_trace}. If \code{two_panel=TRUE}, a second panel shows \code{loglik_trace}.
#'
#' @param X Numeric matrix \code{(n x d)}.
#' @param methods Character vector of initialization methods. Default includes
#'   \code{"PCA"}, \code{"tSNE"}, \code{"random"}, \code{"fiedler"}, \code{"pcurve"}.
#' @param num_iter Integer \eqn{\ge 1}. Number of warm-start iterations to run per method.
#' @param num_cores Integer \eqn{\ge 1}. Number of cores to use (see \code{\link{parallel_initial_csmoothEM}}).
#' @param K Optional integer \eqn{\ge 2}. Number of mixture components. If \code{NULL},
#'   the default logic in \code{initialize_csmoothEM} is used.
#' @param adaptive Character. One of \code{"none"}, \code{"prior"}, \code{"ml"}.
#'   Logical values are accepted for backward compatibility: \code{TRUE} is treated as
#'   \code{"prior"} and \code{FALSE} as \code{"none"}.
#' @param lambda_min,lambda_max Positive bounds for \code{lambda_vec} when \code{adaptive!="none"}.
#' @param sigma_update Character. Only used when \code{adaptive="ml"} (passed through).
#'   Either \code{"mstep"} or \code{"ml"}.
#' @param sigma_min,sigma_max Positive bounds for \code{sigma2} when \code{adaptive="ml"} and
#'   \code{sigma_update="ml"}.
#' @param plot Logical; if \code{TRUE}, plot traces across methods.
#' @param two_panel Logical; if \code{TRUE}, also plot penalized objective traces in a second panel.
#' @param seed Optional integer seed. If provided, a different derived seed is used per method.
#' @param quiet Logical; reserved for future use.
#' @param ... Additional arguments passed to \code{initialize_csmoothEM} (and downstream ordering routines).
#'
#' @return A \code{csmooth_em} object corresponding to the selected best method.
#'   The returned object has an additional \code{meta$initial_search} list containing:
#'   \itemize{
#'     \item \code{best_method}: selected method name
#'     \item \code{criterion}: selection criterion (\code{"ml_last"} or \code{"elbo_last"})
#'     \item \code{summary}: the per-method summary data.frame
#'     \item \code{fits}: the full named list of fits
#'     \item \code{options}: options used for the search
#'   }
#'
#' @seealso \code{\link{parallel_initial_csmoothEM}}, \code{\link{initialize_csmoothEM}}
#' @export
optimize_initial_csmoothEM <- function(
    X,
    methods = c("PCA", "tSNE", "random", "fiedler", "pcurve"),
    num_iter = 5,
    num_cores = 2,
    K = NULL,
    adaptive = "prior",
    lambda_min = 1e-8,
    lambda_max = 1e8,
    sigma_update = c("mstep", "ml"),
    sigma_min = 1e-10,
    sigma_max = 1e10,
    plot = FALSE,
    two_panel = FALSE,
    seed = NULL,
    quiet = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (is.null(adaptive)) adaptive <- "none"
  if (is.logical(adaptive)) adaptive <- if (isTRUE(adaptive)) "prior" else "none"
  adaptive <- match.arg(adaptive, choices = c("none", "prior", "ml"))
  sigma_update <- match.arg(sigma_update)

  fits <- parallel_initial_csmoothEM(
    X = X,
    methods = methods,
    num_iter = num_iter,
    num_cores = num_cores,
    K = K,
    adaptive = adaptive,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_update = sigma_update,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    seed = seed,
    quiet = quiet,
    ...
  )
  sum_df <- attr(fits, "summary")

  if (adaptive == "ml") {
    score <- sum_df$ml_last
    score[!is.finite(score)] <- -Inf
    criterion <- "ml_last"
  } else {
    score <- sum_df$elbo_last
    score[!is.finite(score)] <- -Inf
    criterion <- "elbo_last"
  }

  if (all(score == -Inf)) {
    stop("All initializations failed or produced empty traces under the selection criterion.\n",
         paste(sprintf("  %s: %s", sum_df$method, sum_df$error), collapse = "\n"))
  }

  best_idx <- which.max(score)
  best_method <- sum_df$method[best_idx]
  best <- fits[[best_method]]
  if (is.null(best)) stop("Best method returned NULL (unexpected).")

  best$meta <- best$meta %||% list()
  best$meta$initial_search <- list(
    best_method = best_method,
    criterion = criterion,
    summary = sum_df,
    fits = fits,
    options = list(
      methods = methods,
      num_iter = num_iter,
      num_cores = num_cores,
      K = K,
      adaptive = adaptive,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_update = sigma_update,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      seed = seed
    )
  )

  if (isTRUE(plot)) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    methods_use <- methods
    nM <- length(methods_use)
    cols <- grDevices::hcl.colors(nM, palette = "Dark 3")
    names(cols) <- methods_use

    get_trace <- function(obj, what = c("elbo_trace", "loglik_trace", "ml_trace")) {
      what <- match.arg(what)
      if (is.null(obj)) return(numeric(0))
      tr <- obj[[what]]
      if (is.null(tr)) numeric(0) else as.numeric(tr)
    }

    plot_family <- function(traces, ylab, main, best_method = NULL) {
      first_ok <- NULL
      for (mm in methods_use) {
        if (length(traces[[mm]]) > 0) { first_ok <- mm; break }
      }
      if (is.null(first_ok)) {
        plot.new(); title(paste0(main, " (empty)")); return(invisible(NULL))
      }

      y <- unlist(traces, use.names = FALSE)
      y <- y[is.finite(y)]
      ylim <- if (length(y)) range(y) else c(-1, 1)

      tr0 <- traces[[first_ok]]
      plot(seq_along(tr0), tr0, type = "l",
           col = cols[first_ok], lwd = 2,
           xlab = "Iteration", ylab = ylab, main = main, ylim = ylim)

      for (mm in methods_use) {
        tr <- traces[[mm]]
        if (length(tr) > 0) {
          is_best <- (!is.null(best_method) && identical(mm, best_method))
          lines(seq_along(tr), tr,
                col = cols[mm],
                lwd = if (is_best) 4 else 2,
                lty = if (is_best) 1 else 2)
        }
      }
      legend("bottomright",
             legend = methods_use,
             col = cols[methods_use],
             lwd = ifelse(methods_use == best_method, 4, 2),
             lty = ifelse(methods_use == best_method, 1, 2),
             bty = "n")
      invisible(NULL)
    }

    if (two_panel) par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

    if (adaptive == "ml") {
      ml_traces <- setNames(lapply(methods_use, function(mm) get_trace(fits[[mm]], "ml_trace")), methods_use)
      plot_family(ml_traces, ylab = "Collapsed objective (C)", main = "Collapsed objective traces across initializations",
                  best_method = best_method)
      mtext(sprintf("Best: %s (last C = %.6f)", best_method, sum_df$ml_last[best_idx]),
            side = 3, line = 0.2)
    } else {
      elbo_traces <- setNames(lapply(methods_use, function(mm) get_trace(fits[[mm]], "elbo_trace")), methods_use)
      plot_family(elbo_traces, ylab = "Penalized ELBO", main = "ELBO traces across initializations",
                  best_method = best_method)
      mtext(sprintf("Best: %s (last ELBO = %.6f)", best_method, sum_df$elbo_last[best_idx]),
            side = 3, line = 0.2)
    }

    if (two_panel) {
      obj_traces <- setNames(lapply(methods_use, function(mm) get_trace(fits[[mm]], "loglik_trace")), methods_use)
      plot_family(obj_traces, ylab = "Penalized observed objective",
                  main = "Penalized objective traces across initializations",
                  best_method = best_method)
    }
  }

  best
}
