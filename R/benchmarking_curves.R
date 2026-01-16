#' Simulate a 2D Archimedean spiral (helix-like) with noise
#'
#' @description
#' Simulates a 2D Archimedean spiral curve parameterized by a latent ordering
#' variable \eqn{t \in [0,1]}. The radius increases linearly with \eqn{t} and the
#' angle makes \code{turns} revolutions. Independent Gaussian noise is added to
#' each coordinate.
#'
#' @param n Integer \eqn{\ge 1}. Number of samples.
#' @param turns Positive numeric. Number of revolutions over \eqn{t \in [0,1]}.
#' @param noise Nonnegative numeric. Standard deviation of Gaussian noise added
#'   to both coordinates.
#' @param r0,r1 Nonnegative numerics with \code{r1 >= r0}. Start/end radius.
#' @param seed Optional integer. If not NULL, sets the RNG seed via
#'   \code{set.seed()}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{obs}: numeric matrix \code{(n x 2)} with column names \code{c("x1","x2")}.
#'   \item \code{t}: numeric vector of length \code{n}, the latent ordering in \eqn{[0,1]}.
#' }
#'
#' @examples
#' sim <- simulate_spiral2d(n = 200, turns = 4, noise = 0.05, seed = 1)
#' head(sim$obs)
#' head(sim$t)
#'
#' @export
simulate_spiral2d <- function(n = 500, turns = 3, noise = 0.05,
                              r0 = 0.2, r1 = 1.0, seed = 1) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) stop("n must be a single integer >= 1.")
  if (!is.finite(turns) || turns <= 0) stop("turns must be a positive number.")
  if (!is.finite(noise) || noise < 0) stop("noise must be a nonnegative number.")
  if (!is.finite(r0) || r0 < 0) stop("r0 must be a nonnegative number.")
  if (!is.finite(r1) || r1 < 0 || r1 < r0) stop("r1 must be >= r0 and nonnegative.")
  if (!is.null(seed)) {
    if (length(seed) != 1L || is.na(seed)) stop("seed must be a single integer or NULL.")
    set.seed(seed)
  }

  t <- sort(stats::runif(n))                 # latent ordering in [0,1]
  theta <- 2 * base::pi * turns * t          # angle
  r <- r0 + (r1 - r0) * t                    # radius increases with t

  x1 <- r * cos(theta) + stats::rnorm(n, 0, noise)
  x2 <- r * sin(theta) + stats::rnorm(n, 0, noise)

  X <- cbind(x1, x2)
  colnames(X) <- c("x1", "x2")
  list(obs = X, t = t)
}


#' Simulate a 2D "swiss roll" spiral with a 1D latent parameter
#'
#' @description
#' Simulates a 2D spiral (often used as a 1D analogue of a swiss roll) parameterized
#' by \eqn{t \in [t_{\min}, t_{\max}]}. The noiseless curve is
#' \eqn{(x(t),y(t)) = (t\cos t, t\sin t)}. Independent Gaussian noise is added to
#' each coordinate.
#'
#' @param n Integer \eqn{\ge 1}. Number of samples.
#' @param t_range Numeric vector of length 2 specifying \eqn{[t_{\min}, t_{\max}]}.
#' @param sigma Nonnegative numeric. Standard deviation of Gaussian noise. May be
#'   a scalar (applied to both coordinates) or a length-2 vector for coordinate-specific
#'   noise.
#' @param seed Optional integer. If not NULL, sets the RNG seed via \code{set.seed()}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{t}: numeric vector of length \code{n}, the latent parameter.
#'   \item \code{truth}: data.frame with columns \code{x}, \code{y} for the noiseless curve.
#'   \item \code{obs}: data.frame with columns \code{x}, \code{y} for the noisy observations.
#' }
#'
#' @examples
#' sim <- simulate_swiss_roll_1d_2d(n = 300, sigma = 0.1, seed = 1)
#' head(sim$obs)
#'
#' @export
simulate_swiss_roll_1d_2d <- function(n = 800,
                                      t_range = c(1.5 * base::pi, 6 * base::pi),
                                      sigma = 0.15,
                                      seed = 1) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L) stop("n must be a single integer >= 1.")
  if (!is.numeric(t_range) || length(t_range) != 2L || any(!is.finite(t_range))) {
    stop("t_range must be a finite numeric vector of length 2.")
  }
  if (t_range[2] <= t_range[1]) stop("t_range[2] must be > t_range[1].")

  if (!is.numeric(sigma) || !(length(sigma) == 1L || length(sigma) == 2L) || any(!is.finite(sigma)) || any(sigma < 0)) {
    stop("sigma must be a nonnegative numeric scalar or a length-2 vector.")
  }
  if (length(sigma) == 1L) sigma <- rep(sigma, 2)

  if (!is.null(seed)) {
    if (length(seed) != 1L || is.na(seed)) stop("seed must be a single integer or NULL.")
    set.seed(seed)
  }

  t <- stats::runif(n, min = t_range[1], max = t_range[2])

  x_true <- t * cos(t)
  y_true <- t * sin(t)

  x_obs <- x_true + stats::rnorm(n, 0, sigma[1])
  y_obs <- y_true + stats::rnorm(n, 0, sigma[2])

  list(
    t = t,
    truth = data.frame(x = x_true, y = y_true),
    obs   = data.frame(x = x_obs,  y = y_obs)
  )
}


#' Simulate a two-ordering GP dataset (Matern) for feature partitioning
#'
#' @description
#' Simulates an \eqn{N \times D} dataset with two latent sample orderings \eqn{t_1} and \eqn{t_2}.
#' The first \eqn{D/2} features are generated from a Gaussian process over \eqn{t_1}, and the
#' remaining \eqn{D/2} features are generated from an independent Gaussian process over \eqn{t_2}.
#'
#' Each feature is a GP draw evaluated at the \eqn{N} sample locations. Optionally:
#' \itemize{
#'   \item permute feature columns (\code{permute_cols}),
#'   \item permute the sample order within block 2 (\code{permute_rows_block2}),
#'   \item shift the data to be positive (\code{shift_positive}),
#'   \item add i.i.d. Gaussian noise (\code{noise_sd}).
#' }
#'
#' Column names are assigned *before* permutation and then permuted consistently with the columns,
#' so that names remain aligned with the returned \code{true_group}.
#'
#' @param N Integer \eqn{\ge 2}. Number of samples (rows).
#' @param D Integer \eqn{\ge 2} and even. Number of features (columns).
#' @param t_range Numeric length-2 vector. Range for sampling \code{t1} and \code{t2}.
#' @param range,smoothness,variance Matern GP hyperparameters.
#' @param noise_sd Nonnegative numeric. Standard deviation of i.i.d. Gaussian noise added to \code{X}.
#' @param shift_positive Logical; if TRUE, shift each GP block so its minimum is 1.
#' @param permute_cols Logical; if TRUE, permute feature columns.
#' @param permute_rows_block2 Logical; if TRUE, permute rows of the second GP block before combining.
#' @param seed Optional integer seed. If not NULL, sets \code{set.seed(seed)}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{X}: numeric matrix \eqn{N \times D}.
#'   \item \code{t1}, \code{t2}: numeric vectors length N (latent orderings).
#'   \item \code{permut_cols}: integer vector length D (the applied column permutation; identity if \code{permute_cols=FALSE}).
#'   \item \code{true_group}: integer vector length D in \code{1,2} indicating which ordering generated each feature (after permutation).
#'   \item \code{row_perm_block2}: integer vector length N giving the row permutation applied to block 2 (or \code{NULL}).
#' }
#'
#' @details
#' Requires \code{fields::rdist}, \code{fields::Matern}, and \code{MASS::mvrnorm}.
#'
#' @export
simulate_two_order_gp_dataset <- function(
    N = 1000,
    D = 16,
    t_range = c(0, 10),
    range = 5,
    smoothness = 2.5,
    variance = 3.0,
    noise_sd = 0.05,
    shift_positive = TRUE,
    permute_cols = FALSE,
    permute_rows_block2 = TRUE,
    seed = NULL
) {
  N <- as.integer(N)
  D <- as.integer(D)
  if (length(N) != 1L || is.na(N) || N < 2L) stop("N must be a single integer >= 2.")
  if (length(D) != 1L || is.na(D) || D < 2L || (D %% 2L) != 0L) stop("D must be an even integer >= 2.")
  if (!is.numeric(t_range) || length(t_range) != 2L || any(!is.finite(t_range)) || t_range[2] <= t_range[1]) {
    stop("t_range must be a finite numeric vector length 2 with t_range[2] > t_range[1].")
  }
  if (!is.finite(range) || range <= 0) stop("range must be positive.")
  if (!is.finite(smoothness) || smoothness <= 0) stop("smoothness must be positive.")
  if (!is.finite(variance) || variance <= 0) stop("variance must be positive.")
  if (!is.finite(noise_sd) || noise_sd < 0) stop("noise_sd must be nonnegative.")

  if (!is.null(seed)) {
    if (length(seed) != 1L || is.na(seed)) stop("seed must be a single integer or NULL.")
    set.seed(seed)
  }

  half <- D %/% 2L

  # latent orderings
  t1 <- stats::runif(N, min = t_range[1], max = t_range[2])
  t2 <- stats::runif(N, min = t_range[1], max = t_range[2])

  matern_cov <- function(t) {
    dists <- fields::rdist(t, t)
    cov_mat <- fields::Matern(dists, range = range, smoothness = smoothness)
    cov_mat * variance
  }

  Sigma1 <- matern_cov(matrix(t1, ncol = 1))
  Sigma2 <- matern_cov(matrix(t2, ncol = 1))

  # each column is one GP draw across N samples
  X1 <- t(MASS::mvrnorm(n = half, mu = rep(0, N), Sigma = Sigma1))  # N x half
  X2 <- t(MASS::mvrnorm(n = half, mu = rep(0, N), Sigma = Sigma2))  # N x half

  row_perm_block2 <- NULL
  if (isTRUE(permute_rows_block2)) {
    row_perm_block2 <- sample.int(N)
    X2 <- X2[row_perm_block2, , drop = FALSE]
  }

  if (isTRUE(shift_positive)) {
    X1 <- X1 + abs(min(X1)) + 1
    X2 <- X2 + abs(min(X2)) + 1
  }

  X <- cbind(X1, X2)

  if (noise_sd > 0) {
    X <- X + matrix(stats::rnorm(N * D, mean = 0, sd = noise_sd), nrow = N, ncol = D)
  }

  # names and truth BEFORE permutation
  colnames0 <- c(paste0("ord1_", seq_len(half)), paste0("ord2_", seq_len(half)))
  true_group0 <- c(rep(1L, half), rep(2L, half))

  colnames(X) <- colnames0

  permut_cols <- seq_len(D)
  true_group <- true_group0

  if (isTRUE(permute_cols)) {
    permut_cols <- sample.int(D)
    X <- X[, permut_cols, drop = FALSE]
    colnames(X) <- colnames0[permut_cols]   # IMPORTANT: permute ORIGINAL names once
    true_group <- true_group0[permut_cols]
  }

  list(
    X = X,
    t1 = t1,
    t2 = t2,
    permut_cols = permut_cols,
    true_group = true_group,
    row_perm_block2 = row_perm_block2
  )
}
