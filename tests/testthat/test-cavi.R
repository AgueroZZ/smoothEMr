test_that("simulate_cavi_toy returns a coherent model-matched dataset", {
  sim <- simulate_cavi_toy(
    n = 60,
    d = 10,
    K = 6,
    rw_q = 2,
    seed = 1
  )

  expect_equal(dim(sim$X), c(60, 10))
  expect_length(sim$z, 60)
  expect_equal(dim(sim$mu), c(10, 6))
  expect_length(sim$sigma2, 10)
  expect_length(sim$lambda_vec, 10)
  expect_equal(dim(sim$Q_K), c(6, 6))
  expect_true(all(sim$sigma2 > 0))
  expect_true(all(sim$lambda_vec > 0))
})

test_that("cavi returns valid responsibilities and nondecreasing ELBO on a toy example", {
  sim <- simulate_cavi_toy(
    n = 120,
    d = 12,
    K = 6,
    rw_q = 2,
    seed = 2
  )

  fit <- cavi(
    sim$X,
    K = 6,
    method = "PCA",
    rw_q = 2,
    ridge = 0,
    max_iter = 25,
    tol = 0,
    verbose = FALSE
  )

  expect_s3_class(fit, "cavi")
  expect_equal(dim(fit$gamma), c(120, 6))
  expect_equal(rowSums(fit$gamma), rep(1, 120), tolerance = 1e-8)
  expect_equal(length(fit$elbo_trace), fit$iter + 1L)
  expect_equal(length(fit$loglik_trace), fit$iter + 1L)
  expect_gte(min(diff(fit$elbo_trace)), -1e-8)
  expect_true(all(fit$params$sigma2 > 0))
  expect_true(all(fit$lambda_vec > 0))
  expect_equal(dim(fit$posterior$mean), c(12, 6))
  expect_length(fit$posterior$cov, 12)
})

test_that("lambda induced-prior penalty is opt-in and keeps a monotone penalized ELBO", {
  sim <- simulate_cavi_toy(
    n = 100,
    d = 10,
    K = 6,
    rw_q = 2,
    seed = 21
  )

  fit_base <- cavi(
    sim$X,
    K = 6,
    method = "PCA",
    rw_q = 2,
    max_iter = 5,
    tol = 0,
    verbose = FALSE
  )
  fit_null <- cavi(
    sim$X,
    K = 6,
    method = "PCA",
    rw_q = 2,
    lambda_sd_prior_rate = NULL,
    max_iter = 5,
    tol = 0,
    verbose = FALSE
  )
  fit_zero <- cavi(
    sim$X,
    K = 6,
    method = "PCA",
    rw_q = 2,
    lambda_sd_prior_rate = 0,
    max_iter = 5,
    tol = 0,
    verbose = FALSE
  )
  fit_pen <- cavi(
    sim$X,
    K = 6,
    method = "PCA",
    rw_q = 2,
    lambda_sd_prior_rate = 0.8,
    max_iter = 5,
    tol = 0,
    verbose = FALSE
  )

  expect_equal(fit_null$gamma, fit_base$gamma, tolerance = 1e-10)
  expect_equal(fit_null$lambda_vec, fit_base$lambda_vec, tolerance = 1e-10)
  expect_equal(fit_null$elbo_trace, fit_base$elbo_trace, tolerance = 1e-10)
  expect_equal(fit_zero$gamma, fit_base$gamma, tolerance = 1e-10)
  expect_equal(fit_zero$lambda_vec, fit_base$lambda_vec, tolerance = 1e-10)
  expect_equal(fit_zero$elbo_trace, fit_base$elbo_trace, tolerance = 1e-10)
  expect_null(fit_null$control$lambda_sd_prior_rate)
  expect_null(fit_zero$control$lambda_sd_prior_rate)
  expect_equal(fit_pen$control$lambda_sd_prior_rate, 0.8)
  expect_gte(min(diff(fit_pen$elbo_trace)), -1e-8)
  expect_false(isTRUE(all.equal(fit_pen$lambda_vec, fit_base$lambda_vec, tolerance = 1e-8)))
})

test_that("cavi supports zero-iteration initialization and continuation", {
  sim <- simulate_cavi_toy(
    n = 80,
    d = 8,
    K = 5,
    rw_q = 2,
    seed = 3
  )

  fit0 <- cavi(
    sim$X,
    K = 5,
    method = "PCA",
    rw_q = 2,
    max_iter = 0,
    verbose = FALSE
  )

  expect_s3_class(fit0, "cavi")
  expect_equal(fit0$iter, 0L)
  expect_equal(length(fit0$elbo_trace), 1L)
  expect_equal(length(fit0$loglik_trace), 1L)
  expect_equal(dim(fit0$gamma), c(80, 5))
  expect_equal(rowSums(fit0$gamma), rep(1, 80), tolerance = 1e-8)

  fit1 <- do_cavi(fit0, iter = 1L, verbose = FALSE)
  expect_equal(fit1$iter, 1L)
  expect_equal(length(fit1$elbo_trace), 2L)
  expect_gte(tail(fit1$elbo_trace, 1L), fit0$elbo_trace[1L] - 1e-8)
  expect_equal(fit1$control$gamma_preinit, "direct_penalized")
})

test_that("raw gamma pre-init moments match responsibility-weighted formulas", {
  X <- matrix(c(
    1, 4,
    2, 5,
    3, 6,
    4, 7
  ), nrow = 4, byrow = TRUE)
  R <- matrix(c(
    1.0, 0.0,
    0.7, 0.3,
    0.2, 0.8,
    0.0, 1.0
  ), nrow = 4, byrow = TRUE)

  raw <- getFromNamespace(".cavi_raw_moments_from_responsibilities", "MPCurver")(
    X = X,
    R = R,
    sigma_min = 1e-12,
    sigma_max = 1e12
  )

  Nk_expected <- colSums(R)
  mu_expected <- vapply(seq_len(ncol(R)), function(k) {
    colSums(R[, k] * X) / Nk_expected[k]
  }, numeric(ncol(X)))

  sigma_expected <- vapply(seq_len(ncol(X)), function(j) {
    centered <- sweep(
      matrix(X[, j], nrow = nrow(X), ncol = ncol(R)),
      2L,
      mu_expected[j, ],
      "-"
    )
    sum(R * centered^2) / nrow(X)
  }, numeric(1))

  expect_equal(raw$Nk, Nk_expected, tolerance = 1e-10)
  expect_equal(raw$mu_raw, mu_expected, tolerance = 1e-10)
  expect_equal(raw$sigma2_raw, sigma_expected, tolerance = 1e-10)
})

test_that("gamma-based cavi starts from raw moments and first records a penalized state", {
  sim <- simulate_cavi_toy(
    n = 70,
    d = 7,
    K = 5,
    rw_q = 2,
    seed = 31
  )
  gamma_init <- matrix(0, nrow = nrow(sim$X), ncol = 5)
  gamma_init[cbind(seq_len(nrow(sim$X)), sim$z)] <- 1

  raw <- getFromNamespace(".cavi_raw_moments_from_responsibilities", "MPCurver")(
    X = sim$X,
    R = gamma_init,
    sigma_min = 1e-10,
    sigma_max = 1e10
  )

  fit <- cavi(
    sim$X,
    K = 5,
    responsibilities_init = gamma_init,
    lambda_init = 6,
    fix_lambda = TRUE,
    rw_q = 2,
    max_iter = 0,
    verbose = FALSE
  )

  expect_s3_class(fit, "cavi")
  expect_equal(fit$iter, 0L)
  expect_equal(length(fit$elbo_trace), 1L)
  expect_equal(length(fit$loglik_trace), 1L)
  expect_equal(fit$control$gamma_preinit, "raw_unpenalized")
  expect_equal(fit$gamma, gamma_init, tolerance = 1e-10)
  expect_equal(fit$params$sigma2, raw$sigma2_raw, tolerance = 1e-10)
  expect_equal(fit$lambda_vec, rep(6, ncol(sim$X)), tolerance = 1e-10)
  expect_gt(max(abs(fit$posterior$mean - raw$mu_raw)), 1e-6)
})

test_that("explicit sigma2_init and lambda_init are honored with responsibilities_init", {
  sim <- simulate_cavi_toy(
    n = 70,
    d = 5,
    K = 5,
    rw_q = 2,
    seed = 34
  )
  gamma_init <- matrix(0, nrow = nrow(sim$X), ncol = 5)
  gamma_init[cbind(seq_len(nrow(sim$X)), sim$z)] <- 1

  sigma2_init <- rep(0.25, ncol(sim$X))
  lambda_init <- rep(7, ncol(sim$X))

  fit <- cavi(
    sim$X,
    K = 5,
    responsibilities_init = gamma_init,
    sigma2_init = sigma2_init,
    lambda_init = lambda_init,
    rw_q = 2,
    fix_lambda = FALSE,
    max_iter = 0,
    verbose = FALSE
  )

  expect_equal(fit$params$sigma2, sigma2_init, tolerance = 1e-10)
  expect_equal(fit$lambda_vec, lambda_init, tolerance = 1e-10)
})

test_that("adaptive lambda initialization is derived from raw unpenalized trajectories", {
  sim <- simulate_cavi_toy(
    n = 60,
    d = 6,
    K = 5,
    rw_q = 2,
    seed = 32
  )
  gamma_init <- matrix(0, nrow = nrow(sim$X), ncol = 5)
  gamma_init[cbind(seq_len(nrow(sim$X)), sim$z)] <- 1

  raw <- getFromNamespace(".cavi_raw_moments_from_responsibilities", "MPCurver")(
    X = sim$X,
    R = gamma_init,
    sigma_min = 1e-10,
    sigma_max = 1e10
  )
  qk <- make_random_walk_precision(K = 5, d = 1, q = 2, lambda = 1, ridge = 0)
  qk <- 0.5 * (as.matrix(qk) + t(as.matrix(qk)))
  r_rank <- getFromNamespace(".rw_precision_metadata", "MPCurver")(qk, rw_q = 2)$rank

  lambda_expected <- getFromNamespace(".cavi_lambda_init_from_mu_raw", "MPCurver")(
    mu_raw = raw$mu_raw,
    Q_K = qk,
    r_rank = r_rank,
    lambda_sd_prior_rate = NULL,
    lambda_min = 1e-10,
    lambda_max = 1e10
  )

  fit <- cavi(
    sim$X,
    K = 5,
    responsibilities_init = gamma_init,
    rw_q = 2,
    fix_lambda = FALSE,
    max_iter = 0,
    verbose = FALSE
  )

  expect_equal(fit$control$gamma_preinit, "raw_unpenalized")
  expect_true(all(is.finite(fit$lambda_vec)))
  expect_true(all(fit$lambda_vec >= fit$control$lambda_min))
  expect_true(all(fit$lambda_vec <= fit$control$lambda_max))
  expect_equal(fit$lambda_vec, lambda_expected, tolerance = 1e-10)
})

test_that("adaptive lambda initialization respects the induced prior when requested", {
  sim <- simulate_cavi_toy(
    n = 60,
    d = 6,
    K = 5,
    rw_q = 2,
    seed = 33
  )
  gamma_init <- matrix(0, nrow = nrow(sim$X), ncol = 5)
  gamma_init[cbind(seq_len(nrow(sim$X)), sim$z)] <- 1

  raw <- getFromNamespace(".cavi_raw_moments_from_responsibilities", "MPCurver")(
    X = sim$X,
    R = gamma_init,
    sigma_min = 1e-10,
    sigma_max = 1e10
  )
  qk <- make_random_walk_precision(K = 5, d = 1, q = 2, lambda = 1, ridge = 0)
  qk <- 0.5 * (as.matrix(qk) + t(as.matrix(qk)))
  r_rank <- getFromNamespace(".rw_precision_metadata", "MPCurver")(qk, rw_q = 2)$rank

  lambda_expected <- getFromNamespace(".cavi_lambda_init_from_mu_raw", "MPCurver")(
    mu_raw = raw$mu_raw,
    Q_K = qk,
    r_rank = r_rank,
    lambda_sd_prior_rate = 0.8,
    lambda_min = 1e-10,
    lambda_max = 1e10
  )

  fit <- cavi(
    sim$X,
    K = 5,
    responsibilities_init = gamma_init,
    rw_q = 2,
    fix_lambda = FALSE,
    lambda_sd_prior_rate = 0.8,
    max_iter = 0,
    verbose = FALSE
  )

  expect_equal(fit$control$gamma_preinit, "raw_unpenalized")
  expect_equal(fit$control$lambda_sd_prior_rate, 0.8)
  expect_equal(fit$lambda_vec, lambda_expected, tolerance = 1e-10)
})

test_that("cavi supports known feature-wise measurement sd", {
  sim <- simulate_cavi_toy(
    n = 90,
    d = 7,
    K = 5,
    rw_q = 2,
    seed = 41
  )
  S <- seq(0.09, 0.15, length.out = ncol(sim$X))

  fit <- cavi(
    sim$X,
    K = 5,
    method = "PCA",
    S = S,
    rw_q = 2,
    max_iter = 8,
    tol = 0,
    verbose = FALSE
  )

  expect_s3_class(fit, "cavi")
  expect_equal(fit$control$noise_model, "known_feature_sd")
  expect_equal(fit$measurement_sd, S, tolerance = 1e-12)
  expect_null(fit$params$sigma2)
  expect_length(fit$sigma2_trace, 0L)
  expect_equal(rowSums(fit$gamma), rep(1, nrow(sim$X)), tolerance = 1e-8)
  expect_gte(min(diff(fit$elbo_trace)), -1e-8)
})

test_that("known feature-sd vector and repeated matrix inputs are equivalent", {
  sim <- simulate_cavi_toy(
    n = 60,
    d = 6,
    K = 5,
    rw_q = 2,
    seed = 42
  )
  gamma_init <- matrix(0, nrow = nrow(sim$X), ncol = 5)
  gamma_init[cbind(seq_len(nrow(sim$X)), sim$z)] <- 1
  S_vec <- seq(0.08, 0.14, length.out = ncol(sim$X))
  S_mat <- matrix(rep(S_vec, each = nrow(sim$X)), nrow = nrow(sim$X), ncol = ncol(sim$X))

  fit_vec <- cavi(
    sim$X,
    K = 5,
    responsibilities_init = gamma_init,
    pi_init = colMeans(gamma_init),
    S = S_vec,
    lambda_init = rep(3, ncol(sim$X)),
    fix_lambda = TRUE,
    rw_q = 2,
    max_iter = 3,
    tol = 0,
    verbose = FALSE
  )
  fit_mat <- cavi(
    sim$X,
    K = 5,
    responsibilities_init = gamma_init,
    pi_init = colMeans(gamma_init),
    S = S_mat,
    lambda_init = rep(3, ncol(sim$X)),
    fix_lambda = TRUE,
    rw_q = 2,
    max_iter = 3,
    tol = 0,
    verbose = FALSE
  )

  expect_equal(fit_vec$gamma, fit_mat$gamma, tolerance = 1e-10)
  expect_equal(fit_vec$posterior$mean, fit_mat$posterior$mean, tolerance = 1e-10)
  expect_equal(fit_vec$elbo_trace, fit_mat$elbo_trace, tolerance = 1e-10)
})

test_that("do_cavi reuses known observation-level measurement sd", {
  sim <- simulate_cavi_toy(
    n = 50,
    d = 5,
    K = 4,
    rw_q = 2,
    seed = 43
  )
  base_sd <- seq(0.07, 0.11, length.out = ncol(sim$X))
  row_scale <- seq(1, 1.15, length.out = nrow(sim$X))
  S <- outer(row_scale, base_sd)

  fit0 <- cavi(
    sim$X,
    K = 4,
    method = "PCA",
    S = S,
    rw_q = 2,
    max_iter = 0,
    verbose = FALSE
  )
  fit1 <- do_cavi(fit0, iter = 2, tol = 0, verbose = FALSE)

  expect_equal(dim(fit1$measurement_sd), dim(S))
  expect_equal(fit1$measurement_sd, S, tolerance = 1e-12)
  expect_equal(fit1$control$noise_model, "known_observation_sd")
  expect_null(fit1$params$sigma2)
  expect_gte(min(diff(fit1$elbo_trace)), -1e-8)
  expect_error(
    do_cavi(fit0, iter = 1, S = S * 1.01, verbose = FALSE),
    "does not match"
  )
})
