test_that("soft_two_trajectory_cavi returns a stable A/B interface", {
  sim <- simulate_dual_trajectory(
    n = 100,
    d1 = 5,
    d2 = 5,
    d_noise = 0,
    sigma = 0.1,
    seed = 2
  )

  res <- soft_two_trajectory_cavi(
    sim$X,
    K = 10,
    n_outer = 4,
    max_converge_iter = 6,
    inner_iter = 1L,
    verbose = FALSE
  )

  expect_equal(dim(res$pi_weights), c(ncol(sim$X), 2L))
  expect_identical(colnames(res$pi_weights), c("A", "B"))
  expect_true(all(res$assign %in% c("A", "B")))
  expect_true(all(is.finite(res$objective_history)))
  expect_true(all(diff(res$objective_history) >= -1e-8))
  expect_gte(evaluate_partition(list(assign = res$assign), sim$true_assign)$accuracy, 0.8)
})

test_that("soft partition accepts full-data fits rebuilt from subset gamma starts", {
  sim <- simulate_dual_trajectory(
    n = 70,
    d1 = 4,
    d2 = 4,
    d_noise = 0,
    sigma = 0.1,
    seed = 7
  )

  idx_a <- which(sim$true_assign == "A")
  idx_b <- which(sim$true_assign == "B")

  fit_a_sub <- cavi(
    sim$X[, idx_a, drop = FALSE],
    K = 7,
    method = "PCA",
    rw_q = 2,
    max_iter = 2,
    verbose = FALSE
  )
  fit_b_sub <- cavi(
    sim$X[, idx_b, drop = FALSE],
    K = 7,
    method = "PCA",
    rw_q = 2,
    max_iter = 2,
    verbose = FALSE
  )

  fit_a_full <- cavi(
    sim$X,
    K = 7,
    responsibilities_init = fit_a_sub$gamma,
    rw_q = 2,
    max_iter = 0,
    verbose = FALSE
  )
  fit_b_full <- cavi(
    sim$X,
    K = 7,
    responsibilities_init = fit_b_sub$gamma,
    rw_q = 2,
    max_iter = 0,
    verbose = FALSE
  )

  res <- soft_partition_cavi(
    sim$X,
    M = 2,
    fits_init = list(fit_a_full, fit_b_full),
    n_outer = 2,
    inner_iter = 1L,
    max_converge_iter = 2L,
    verbose = FALSE
  )

  expect_equal(fit_a_full$control$gamma_preinit, "raw_unpenalized")
  expect_equal(fit_b_full$control$gamma_preinit, "raw_unpenalized")
  expect_s3_class(res, "soft_partition_cavi")
  expect_equal(length(res$fits), 2L)
  expect_true(all(is.finite(res$objective_history)))
  expect_true(all(diff(res$objective_history) >= -1e-8))
})

test_that("soft_partition_cavi supports similarity-based initialization metadata", {
  sim <- simulate_intrinsic_trajectories(
    n = 90,
    d_signal = c(4, 4, 4),
    d_noise = 0,
    sigma = 0.08,
    seed = 11,
    trajectory_family = rep("linear", 3)
  )

  res <- suppressWarnings(soft_partition_cavi(
    sim$X,
    M = 3,
    init_methods = c("PCA", "random", "PCA"),
    partition_init = "similarity",
    similarity_metric = "spearman",
    cluster_linkage = "single",
    similarity_min_feature_sd = 1e-8,
    K = 7,
    n_outer = 2,
    inner_iter = 1L,
    max_converge_iter = 2L,
    verbose = FALSE
  ))

  expect_s3_class(res, "soft_partition_cavi")
  expect_equal(res$control$partition_init, "similarity")
  expect_equal(res$control$similarity_metric, "spearman")
  expect_equal(res$control$cluster_linkage, "single")
  expect_true(!is.null(res$similarity_init))
  expect_true(is.matrix(res$similarity_init$S))
  expect_true(is.matrix(res$similarity_init$distance))
  expect_equal(dim(res$similarity_init$S), c(ncol(sim$X), ncol(sim$X)))
  expect_equal(dim(res$similarity_init$distance), c(ncol(sim$X), ncol(sim$X)))
  expect_equal(unname(diag(res$similarity_init$S)), rep(1, ncol(sim$X)))
  expect_length(res$similarity_init$feature_cluster, ncol(sim$X))
  expect_equal(sum(res$similarity_init$cluster_sizes), ncol(sim$X))
  expect_setequal(unique(unname(res$similarity_init$feature_cluster)), 1:3)
  expect_null(res$similarity_init$directional_score)
  expect_null(res$similarity_init$directional_delta)
  expect_null(res$similarity_init$directional_success)
  expect_equal(
    unname(vapply(res$init_info, `[[`, character(1), "method_requested")),
    c("PCA", "random", "PCA")
  )
})

test_that("soft_partition_cavi stores smooth-fit similarity diagnostics", {
  sim <- simulate_intrinsic_trajectories(
    n = 50,
    d_signal = c(3, 3),
    d_noise = 0,
    sigma = 0.08,
    seed = 21,
    trajectory_family = c("monotone", "quadratic")
  )

  res <- suppressWarnings(soft_partition_cavi(
    sim$X,
    M = 2,
    partition_init = "similarity",
    similarity_metric = "smooth_fit",
    cluster_linkage = "single",
    K = 6,
    ridge = 1e-6,
    n_outer = 2,
    inner_iter = 1L,
    max_converge_iter = 2L,
    verbose = FALSE
  ))

  S <- res$similarity_init$S
  expect_true(is.matrix(S))
  expect_equal(dim(S), c(ncol(sim$X), ncol(sim$X)))
  expect_equal(unname(diag(S)), rep(1, ncol(sim$X)))
  expect_true(all(is.finite(S)))
  expect_true(all(S >= 0 & S <= 1))
  expect_equal(S, t(S), tolerance = 1e-10)
  expect_true(is.matrix(res$similarity_init$distance))
  expect_true(is.matrix(res$similarity_init$directional_score))
  expect_true(is.matrix(res$similarity_init$directional_delta))
  expect_true(is.matrix(res$similarity_init$directional_success))
  expect_true(is.matrix(res$similarity_init$directional_lambda))
  expect_true(is.matrix(res$similarity_init$directional_sigma2))
  expect_true(all(is.finite(res$similarity_init$directional_score)))
  expect_true(all(is.finite(res$similarity_init$directional_delta)))
  expect_true(all(is.finite(res$similarity_init$directional_lambda)))
  expect_true(all(is.finite(res$similarity_init$directional_sigma2)))
  expect_equal(dim(res$similarity_init$directional_score), dim(S))
  expect_equal(dim(res$similarity_init$directional_delta), dim(S))
  expect_equal(dim(res$similarity_init$directional_success), dim(S))
  expect_equal(dim(res$similarity_init$directional_lambda), dim(S))
  expect_equal(dim(res$similarity_init$directional_sigma2), dim(S))
})

test_that("smooth-fit similarity supports fixed lambda mode and warns for ridge = 0", {
  sim <- simulate_intrinsic_trajectories(
    n = 50,
    d_signal = c(3, 3),
    d_noise = 0,
    sigma = 0.08,
    seed = 22,
    trajectory_family = c("monotone", "quadratic")
  )

  expect_warning(
    res_intrinsic <- soft_partition_cavi(
      sim$X,
      M = 2,
      partition_init = "similarity",
      similarity_metric = "smooth_fit",
      K = 6,
      ridge = 0,
      n_outer = 2,
      inner_iter = 1L,
      max_converge_iter = 2L,
      verbose = FALSE
    ),
    "pseudo-evidence"
  )
  expect_true(is.matrix(res_intrinsic$similarity_init$S))
  expect_true(all(is.finite(res_intrinsic$similarity_init$directional_score)))
  expect_true(all(is.finite(res_intrinsic$similarity_init$directional_delta)))
  expect_true(all(is.finite(res_intrinsic$similarity_init$directional_lambda)))
  expect_true(all(is.finite(res_intrinsic$similarity_init$directional_sigma2)))

  res_fixed <- suppressWarnings(soft_partition_cavi(
    sim$X,
    M = 2,
    partition_init = "similarity",
    similarity_metric = "smooth_fit",
    smooth_fit_lambda_mode = "fixed",
    smooth_fit_lambda_value = 2,
    K = 6,
    ridge = 1e-6,
    n_outer = 2,
    inner_iter = 1L,
    max_converge_iter = 2L,
    verbose = FALSE
  ))

  expect_equal(res_fixed$control$smooth_fit_lambda_mode, "fixed")
  expect_equal(res_fixed$control$smooth_fit_lambda_value, 2)
  expect_true(all(res_fixed$similarity_init$directional_lambda > 0))
})

test_that("soft_partition_cavi supports known measurement sd with smooth-fit initialization", {
  sim <- simulate_intrinsic_trajectories(
    n = 45,
    d_signal = c(3, 3),
    d_noise = 0,
    sigma = 0.08,
    seed = 45,
    trajectory_family = c("monotone", "quadratic")
  )
  base_sd <- seq(0.07, 0.12, length.out = ncol(sim$X))
  row_scale <- seq(1, 1.1, length.out = nrow(sim$X))
  S <- outer(row_scale, base_sd)

  res <- suppressWarnings(soft_partition_cavi(
    sim$X,
    S = S,
    M = 2,
    partition_init = "similarity",
    similarity_metric = "smooth_fit",
    K = 6,
    ridge = 1e-6,
    n_outer = 2,
    inner_iter = 1L,
    max_converge_iter = 2L,
    verbose = FALSE
  ))

  expect_s3_class(res, "soft_partition_cavi")
  expect_equal(res$control$noise_model, "known_observation_sd")
  expect_equal(dim(res$measurement_sd), dim(S))
  expect_equal(res$measurement_sd, S, tolerance = 1e-12)
  expect_true(all(vapply(res$fits, function(f) !is.null(f$measurement_sd), logical(1))))
  expect_true(all(vapply(res$fits, function(f) {
    identical(f$control$noise_model, "known_observation_sd")
  }, logical(1))))
  expect_true(all(is.finite(res$objective_history)))
  # Only check monotonicity in convergence phase (post-annealing).
  # During annealing the temperature changes, so the objective function
  # itself changes and monotonicity is not guaranteed.
  converge_idx <- seq(from = res$n_anneal, to = length(res$objective_history))
  if (length(converge_idx) >= 2L) {
    expect_true(all(diff(res$objective_history[converge_idx]) >= -1e-8))
  }
})

test_that("fit_mpcurve can use similarity initialization for intrinsic_dim > 1", {
  sim <- simulate_intrinsic_trajectories(
    n = 80,
    d_signal = c(4, 4, 4),
    d_noise = 0,
    sigma = 0.08,
    seed = 12,
    trajectory_family = rep("monotone", 3)
  )

  fit <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 3,
    partition_init = "similarity",
    similarity_metric = "spearman",
    cluster_linkage = "single",
    K = 7,
    n_outer = 2,
    inner_iter = 1L,
    max_converge_iter = 2L,
    verbose = FALSE
  ))

  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$algorithm, "cavi")
  expect_equal(fit$requested_intrinsic_dim, 3L)
  expect_equal(fit$fit$control$partition_init, "similarity")
  expect_true(!is.null(fit$similarity_init))
  expect_true(is.matrix(fit$similarity_init$S))
  expect_equal(length(fit$fits), 3L)
})
