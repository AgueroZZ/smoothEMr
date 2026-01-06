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
  Q <- as.matrix(Q)

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


