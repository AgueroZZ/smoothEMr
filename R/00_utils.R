# internal helper
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Generalized log-determinant (sum of log positive eigenvalues)
#'
#' Computes \eqn{\sum_i \log(\lambda_i)} over eigenvalues \eqn{\lambda_i}
#' of a symmetric matrix \code{Q} that exceed \code{eigen_tol}.
#'
#' @param Q Symmetric matrix.
#' @param eigen_tol Nonnegative threshold; eigenvalues <= eigen_tol are ignored.
#'   If \code{NULL}, uses a heuristic based on machine precision.
#' @param rank_deficiency Optional integer. If > 0, drops the smallest
#'   \code{rank_deficiency} eigenvalues (useful when you know the null-space
#'   dimension a priori).
#' @return Scalar generalized log-determinant.
generalized_logdet <- function(Q, eigen_tol = NULL, rank_deficiency = 0) {
  # Avoid unnecessary densification of sparse matrices
  if (!is.matrix(Q)) Q <- as.matrix(Q)

  # eigen() returns values in decreasing order for symmetric matrices
  vals <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values

  if (rank_deficiency < 0) stop("rank_deficiency must be >= 0.")
  if (rank_deficiency > 0) {
    if (rank_deficiency >= length(vals)) {
      stop("rank_deficiency must be < ncol(Q).")
    }
    # drop the smallest 'rank_deficiency' eigenvalues (tail)
    vals <- vals[seq_len(length(vals) - rank_deficiency)]
  }

  if (is.null(eigen_tol)) {
    # Heuristic tolerance: scale by matrix size and largest eigenvalue magnitude
    eigen_tol <- length(vals) * max(abs(vals)) * .Machine$double.eps
  }
  if (eigen_tol < 0) stop("eigen_tol must be >= 0 (or NULL).")

  positive <- vals[vals > eigen_tol]
  if (!length(positive)) {
    warning("No eigenvalues above tolerance; returning -Inf.")
    return(-Inf)
  }

  sum(log(positive))
}


#' Internal metadata for RW precision normalizing constants
#'
#' Distinguishes the intrinsic RW(q) case from a properized/full-rank precision
#' (for example after adding a ridge nugget). In the proper case the Gaussian
#' prior normalizer uses the full matrix rank; in the intrinsic case we retain
#' the known RW(q) null-space dimension.
#'
#' @param Q Symmetric precision matrix.
#' @param rw_q Expected RW(q) nullity in the intrinsic case.
#' @param eigen_tol Optional eigenvalue threshold used to decide whether
#'   \code{Q} is numerically full-rank.
#' @return A list with \code{rank}, \code{logdet}, \code{proper},
#'   \code{nullity}, and \code{eigen_tol}.
#' @keywords internal
.rw_precision_metadata <- function(Q, rw_q = 0L, eigen_tol = NULL) {
  if (!is.matrix(Q)) Q <- as.matrix(Q)

  vals <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  if (is.null(eigen_tol)) {
    eigen_tol <- length(vals) * max(abs(vals)) * .Machine$double.eps
  }
  if (eigen_tol < 0) stop("eigen_tol must be >= 0 (or NULL).")

  positive <- vals[vals > eigen_tol]
  nullity <- length(vals) - length(positive)
  proper <- nullity == 0L

  if (proper) {
    return(list(
      rank = length(vals),
      logdet = sum(log(vals)),
      proper = TRUE,
      nullity = 0L,
      eigen_tol = eigen_tol
    ))
  }

  rw_q <- as.integer(rw_q)
  if (length(rw_q) != 1L || is.na(rw_q) || rw_q < 0L || rw_q >= length(vals)) {
    stop("rw_q must be a single integer in {0, ..., ncol(Q)-1}.")
  }

  list(
    rank = max(length(vals) - rw_q, 1L),
    logdet = generalized_logdet(Q, eigen_tol = eigen_tol, rank_deficiency = rw_q),
    proper = FALSE,
    nullity = nullity,
    eigen_tol = eigen_tol
  )
}


#' Normalize the optional exponential rate for the lambda SD prior
#'
#' Public APIs use \code{NULL} to mean "no lambda prior penalty". For backward
#' compatibility, an explicit numeric \code{0} is treated the same way.
#'
#' @param rate Optional non-negative scalar.
#' @param arg_name Argument name used in error messages.
#' @return Either \code{NULL} (penalty off) or a strictly positive numeric
#'   scalar.
#' @keywords internal
.normalize_lambda_sd_prior_rate <- function(rate,
                                            arg_name = "lambda_sd_prior_rate") {
  if (is.null(rate)) return(NULL)
  rate <- as.numeric(rate)[1]
  if (!is.finite(rate) || rate < 0) {
    stop(arg_name, " must be NULL or a single finite non-negative number.")
  }
  if (rate == 0) return(NULL)
  rate
}


#' Numeric value of the optional lambda SD prior rate
#'
#' Converts the public \code{NULL}-means-off convention into a numeric value
#' used internally by ELBO bookkeeping helpers.
#'
#' @param rate Optional rate.
#' @return Zero when the penalty is off, otherwise the positive numeric rate.
#' @keywords internal
.lambda_sd_prior_rate_value <- function(rate) {
  rate <- .normalize_lambda_sd_prior_rate(rate)
  if (is.null(rate)) 0 else rate
}


#' Log-prior contribution induced by an exponential prior on 1/sqrt(lambda)
#'
#' For \eqn{\tau = 1 / \sqrt{\lambda}} with
#' \eqn{\tau \sim \mathrm{Exp}(\rho)}, the induced prior on \eqn{\lambda > 0}
#' contributes
#' \deqn{
#' \log p(\lambda)
#' =
#' \log(\rho / 2)
#' - \frac{3}{2}\log \lambda
#' - \rho \lambda^{-1/2}.
#' }
#'
#' @param lambda_vec Positive numeric vector.
#' @param rate Non-negative scalar rate parameter \eqn{\rho}.
#' @param include_constant Logical; include \eqn{\log(\rho / 2)} when
#'   \code{rate > 0}?
#' @return Numeric vector of log-prior contributions.
#' @keywords internal
.lambda_sd_prior_terms <- function(lambda_vec, rate, include_constant = TRUE) {
  lambda_vec <- as.numeric(lambda_vec)
  rate <- .lambda_sd_prior_rate_value(rate)
  if (!length(lambda_vec)) return(numeric(0))
  if (any(!is.finite(lambda_vec)) || any(lambda_vec <= 0)) {
    stop("lambda_vec must contain positive finite values.")
  }
  if (rate == 0) return(rep(0, length(lambda_vec)))

  out <- -1.5 * log(lambda_vec) - rate / sqrt(lambda_vec)
  if (isTRUE(include_constant)) {
    out <- out + log(rate / 2)
  }
  out
}


#' Optimize lambda under the induced exponential prior on 1/sqrt(lambda)
#'
#' Solves the penalized one-dimensional problem
#' \deqn{
#' \max_{\lambda \in [\lambda_{\min}, \lambda_{\max}]}
#' \frac{r_Q - 3}{2}\log \lambda
#' - \frac{q}{2}\lambda
#' - \rho \lambda^{-1/2},
#' }
#' where \eqn{q} is the current posterior second moment
#' \eqn{m^\top Q m + \mathrm{tr}(Q S)}.
#'
#' The optimization is carried out in \eqn{\eta = \log \lambda}, for which the
#' objective is strictly concave and the derivative is monotone decreasing.
#'
#' @param eq_quad Positive numeric vector of posterior second moments.
#' @param r_rank Positive scalar RW precision rank.
#' @param rate Positive scalar rate parameter \eqn{\rho}.
#' @param lambda_min,lambda_max Positive bounds with \code{lambda_min <= lambda_max}.
#' @return Numeric vector of optimized \code{lambda}.
#' @keywords internal
.optimize_lambda_induced_exp <- function(eq_quad, r_rank, rate,
                                         lambda_min, lambda_max) {
  eq_quad <- as.numeric(eq_quad)
  r_rank <- as.numeric(r_rank)[1]
  rate <- as.numeric(rate)[1]
  lambda_min <- as.numeric(lambda_min)[1]
  lambda_max <- as.numeric(lambda_max)[1]

  if (!is.finite(r_rank) || r_rank <= 0) {
    stop("r_rank must be a single positive finite number.")
  }
  if (!is.finite(rate) || rate <= 0) {
    stop("rate must be a single positive finite number.")
  }
  if (!is.finite(lambda_min) || !is.finite(lambda_max) ||
      lambda_min <= 0 || lambda_max <= 0 || lambda_min > lambda_max) {
    stop("lambda_min/lambda_max must be positive finite numbers with lambda_min <= lambda_max.")
  }
  if (!length(eq_quad)) return(numeric(0))
  if (any(!is.finite(eq_quad)) || any(eq_quad <= 0)) {
    stop("eq_quad must contain positive finite values.")
  }
  if (lambda_min == lambda_max) {
    return(rep(lambda_min, length(eq_quad)))
  }

  eta_lo <- log(lambda_min)
  eta_hi <- log(lambda_max)

  solve_one <- function(qj) {
    deriv_eta <- function(eta) {
      0.5 * (r_rank - 3) -
        0.5 * qj * exp(eta) +
        0.5 * rate * exp(-eta / 2)
    }

    d_lo <- deriv_eta(eta_lo)
    d_hi <- deriv_eta(eta_hi)

    if (!is.finite(d_lo) || !is.finite(d_hi)) {
      stop("Non-finite derivative encountered while updating lambda.")
    }
    if (d_lo <= 0) return(lambda_min)
    if (d_hi >= 0) return(lambda_max)

    exp(stats::uniroot(deriv_eta, lower = eta_lo, upper = eta_hi,
                       tol = 1e-10)$root)
  }

  vapply(eq_quad, solve_one, numeric(1))
}


#' Numerically stable log-sum-exp
#'
#' @param x Numeric vector.
#' @return log(sum(exp(x))) computed stably.
logsumexp <- function(x) {
  if (!length(x)) return(-Inf)
  m <- max(x)
  if (!is.finite(m)) return(m)  # handles all -Inf, or Inf present
  m + log(sum(exp(x - m)))
}


#' Cache covariance inverses and log-determinants
#'
#' Adds \code{invSigma} and \code{logdet} to a parameter list.
#'
#' @param params A list with at least \code{pi} and \code{sigma}, where
#'   \code{sigma} is a list of covariance matrices.
#' @param jitter Nonnegative diagonal jitter added if Cholesky fails.
#' @return Updated params list with \code{invSigma} and \code{logdet}.
init_cov_cache_fast <- function(params, jitter = 0) {
  K <- length(params$pi)
  if (is.null(params$sigma) || length(params$sigma) != K) {
    stop("params$sigma must be a list of length length(params$pi).")
  }
  if (jitter < 0) stop("jitter must be >= 0.")

  invSigma <- vector("list", K)
  logdet   <- numeric(K)

  for (k in seq_len(K)) {
    S <- as.matrix(params$sigma[[k]])

    # Add jitter if requested
    if (jitter > 0) {
      S <- S + diag(jitter, nrow(S))
    }

    # Cholesky (with a clear error if not PD)
    L <- tryCatch(chol(S), error = function(e) {
      stop("Cholesky failed for sigma[[", k, "]] (matrix not PD). Consider adding jitter.")
    })

    invSigma[[k]] <- chol2inv(L)
    logdet[k]     <- 2 * sum(log(diag(L)))
  }

  params$invSigma <- invSigma
  params$logdet   <- logdet
  params
}


#' Cache inverse-sigma2 and log-determinant for csmooth_em (diagonal covariance)
#'
#' For \code{csmooth_em} the covariance is diagonal, so no Cholesky is needed.
#' Stores \code{invsig2} and \code{logdet} into \code{params} so that
#' \code{ESTEP_csmooth()} can skip recomputing them on every call.
#'
#' @param params A \code{csmooth_em} parameter list with a \code{sigma2} field.
#' @param modelName One of \code{"homoskedastic"} or \code{"heteroskedastic"}.
#' @return Updated \code{params} with \code{$invsig2} and \code{$logdet}.
cache_csmooth_params <- function(params,
    modelName = c("homoskedastic", "heteroskedastic")) {
  modelName <- match.arg(modelName)
  if (modelName == "homoskedastic") {
    s2 <- pmax(as.numeric(params$sigma2), .Machine$double.eps)
    params$invsig2 <- 1 / s2          # d-vector
    params$logdet  <- sum(log(s2))    # scalar
  } else {
    s2 <- pmax(as.matrix(params$sigma2), .Machine$double.eps)
    params$invsig2 <- 1 / s2          # d x K matrix
    params$logdet  <- colSums(log(s2)) # K-vector
  }
  params
}
