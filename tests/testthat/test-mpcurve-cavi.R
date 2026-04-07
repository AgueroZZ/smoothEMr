test_that("fit_mpcurve defaults to cavi and returns a valid mpcurve object", {
  sim <- simulate_cavi_toy(
    n = 100,
    d = 10,
    K = 6,
    rw_q = 2,
    seed = 11
  )

  fit <- fit_mpcurve(
    sim$X,
    K = 6,
    iter = 40,
    num_cores = 1
  )

  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$algorithm, "cavi")
  expect_equal(dim(fit$params$mu), c(10, 6))
  expect_equal(dim(fit$gamma), c(100, 6))
  expect_true(!is.null(fit$locations))
  expect_length(fit$locations$mean$index, 100)
  expect_length(fit$locations$map$index, 100)
  expect_true(all(fit$locations$mean$pseudotime >= 0 & fit$locations$mean$pseudotime <= 1))
  expect_true(all(fit$locations$map$pseudotime >= 0 & fit$locations$map$pseudotime <= 1))
  expect_equal(length(fit$elbo_trace), fit$iter + 1L)
  expect_gte(min(diff(fit$elbo_trace)), -1e-8)

  s <- summary(fit)
  expect_s3_class(s, "summary.mpcurve")
  expect_s3_class(s$underlying, "summary.cavi")
})

test_that("do_cavi and do_mpcurve extend cavi traces", {
  sim <- simulate_cavi_toy(
    n = 90,
    d = 8,
    K = 5,
    rw_q = 2,
    seed = 12
  )

  raw_fit <- cavi(
    sim$X,
    K = 5,
    method = "PCA",
    rw_q = 2,
    max_iter = 4,
    tol = 0,
    verbose = FALSE
  )
  raw_more <- do_cavi(raw_fit, iter = 3, tol = 0, verbose = FALSE)

  expect_s3_class(raw_more, "cavi")
  expect_gt(length(raw_more$elbo_trace), length(raw_fit$elbo_trace))
  expect_equal(length(raw_more$elbo_trace), raw_more$iter + 1L)
  expect_gte(min(diff(raw_more$elbo_trace)), -1e-8)

  mp_fit <- as_mpcurve(raw_fit)
  mp_more <- do_mpcurve(mp_fit, iter = 3, verbose = FALSE)

  expect_s3_class(mp_more, "mpcurve")
  expect_equal(mp_more$algorithm, "cavi")
  expect_gt(length(mp_more$elbo_trace), length(mp_fit$elbo_trace))
  expect_equal(length(mp_more$elbo_trace), mp_more$iter + 1L)
  expect_equal(length(mp_more$locations$mean$index), 90)
  expect_gte(min(diff(mp_more$elbo_trace)), -1e-8)
})

test_that("fit_mpcurve and do_mpcurve expose a cavi-only public wrapper", {
  fit_formals <- names(formals(fit_mpcurve))
  do_formals <- names(formals(do_mpcurve))

  expect_true("greedy" %in% fit_formals)
  expect_true("S" %in% fit_formals)
  expect_false("algorithm" %in% fit_formals)
  expect_false("relative_lambda" %in% fit_formals)
  expect_false("adaptive" %in% fit_formals)
  expect_false("sigma_update" %in% fit_formals)
  expect_false("check_decrease" %in% fit_formals)
  expect_false("tol_decrease" %in% fit_formals)

  expect_false("adaptive" %in% do_formals)
  expect_true("S" %in% do_formals)
  expect_false("sigma_update" %in% do_formals)
  expect_false("check_decrease" %in% do_formals)
  expect_false("tol_decrease" %in% do_formals)
})

test_that("fit_mpcurve and do_mpcurve preserve known measurement sd", {
  sim <- simulate_cavi_toy(
    n = 70,
    d = 6,
    K = 5,
    rw_q = 2,
    seed = 44
  )
  S <- seq(0.09, 0.14, length.out = ncol(sim$X))

  fit <- fit_mpcurve(
    sim$X,
    K = 5,
    S = S,
    iter = 4,
    tol = 0,
    verbose = FALSE
  )
  fit_more <- do_mpcurve(fit, iter = 2, tol = 0, verbose = FALSE)

  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$measurement_sd, S, tolerance = 1e-12)
  expect_equal(fit$fit$control$noise_model, "known_feature_sd")
  expect_null(fit$params$sigma2)
  expect_equal(fit_more$measurement_sd, S, tolerance = 1e-12)
  expect_equal(fit_more$fit$control$noise_model, "known_feature_sd")
  expect_null(fit_more$params$sigma2)
  expect_gte(min(diff(fit_more$elbo_trace)), -1e-8)
})

test_that("greedy forward selection uses fixed-M comparison semantics", {
  sim <- simulate_intrinsic_trajectories(
    n = 100,
    d_signal = c(4, 4),
    d_noise = 0,
    sigma = 0.08,
    seed = 21,
    trajectory_family = c("linear", "linear")
  )

  fit <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 4,
    greedy = "forward",
    partition_init = "similarity",
    method = "isomap",
    K = 8,
    n_outer = 3L,
    inner_iter = 1L,
    max_converge_iter = 4L,
    drop_unused_ordering = FALSE,
    verbose = FALSE
  ))

  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$greedy_selection$mode, "forward")
  expect_equal(fit$greedy_selection$requested_upper_bound, 4L)
  expect_equal(fit$greedy_selection$selected_intrinsic_dim, 2L)
  expect_true(nrow(fit$greedy_selection$history) >= 1L)
  expect_false(isTRUE(fit$control$drop_unused_ordering))
})

test_that("greedy backward can compact after internal comparison search", {
  sim <- simulate_intrinsic_trajectories(
    n = 100,
    d_signal = c(4, 4),
    d_noise = 0,
    sigma = 0.08,
    seed = 22,
    trajectory_family = c("monotone", "monotone")
  )

  fit <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 4,
    greedy = "backward",
    partition_init = "similarity",
    method = "isomap",
    K = 8,
    n_outer = 3L,
    inner_iter = 1L,
    max_converge_iter = 4L,
    drop_unused_ordering = TRUE,
    verbose = FALSE
  ))

  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$greedy_selection$mode, "backward")
  expect_equal(fit$greedy_selection$requested_upper_bound, 4L)
  expect_true(fit$greedy_selection$starting_active_intrinsic_dim <= 4L)
  expect_true(fit$greedy_selection$selected_intrinsic_dim >= 1L)
  expect_true(isTRUE(fit$fit$control$drop_unused_ordering))
  expect_equal(fit$requested_intrinsic_dim, fit$greedy_selection$selected_intrinsic_dim)
})

test_that("legacy hard CAVI partition selectors are no longer exported", {
  ns_exports <- getNamespaceExports("MPCurver")
  expect_false("forward_two_ordering_partition_cavi" %in% ns_exports)
  expect_false("backward_two_ordering_partition_cavi" %in% ns_exports)
})

test_that("legacy wrapper controls are rejected cleanly", {
  sim <- simulate_cavi_toy(
    n = 60,
    d = 6,
    K = 4,
    rw_q = 2,
    seed = 13
  )

  expect_warning(
    fit_mpcurve(
      sim$X,
      algorithm = "cavi",
      method = "PCA",
      K = 4,
      iter = 2,
      verbose = FALSE
    ),
    "no longer part of the public fit_mpcurve"
  )

  expect_error(
    fit_mpcurve(
      sim$X,
      algorithm = "csmooth_em",
      method = "PCA",
      K = 4,
      iter = 2
    ),
    "CAVI-only public wrapper"
  )

  expect_error(
    fit_mpcurve(
      sim$X,
      K = 4,
      iter = 2,
      adaptive = "ml"
    ),
    "Legacy smoothEM/csmoothEM controls"
  )

  legacy_raw <- initialize_csmoothEM(sim$X, method = "PCA", K = 4, num_iter = 1)
  legacy_mp <- as_mpcurve(legacy_raw)
  expect_error(
    do_mpcurve(legacy_mp, iter = 1),
    "CAVI-only public wrapper"
  )
})
