test_that("public namespace exports the curated mpcurve-first API", {
  expect_setequal(
    getNamespaceExports("MPCurver"),
    c(
      "PCA_ordering",
      "do_mpcurve",
      "fiedler_ordering",
      "fit_mpcurve",
      "fitted_prior",
      "isomap_ordering",
      "pcurve_ordering",
      "simulate_cavi_toy",
      "simulate_dual_trajectory",
      "simulate_intrinsic_trajectories",
      "simulate_spiral2d",
      "simulate_swiss_roll_1d_2d",
      "simulate_two_order_gp_dataset",
      "tSNE_ordering"
    )
  )

  namespace_text <- paste(
    readLines(testthat::test_path("..", "..", "NAMESPACE")),
    collapse = "\n"
  )

  expect_false(grepl("S3method\\(print,cavi\\)", namespace_text, fixed = FALSE))
  expect_false(grepl("S3method\\(summary,cavi\\)", namespace_text, fixed = FALSE))
  expect_false(grepl("S3method\\(plot,cavi\\)", namespace_text, fixed = FALSE))
  expect_true(is.function(utils::getS3method("print", "mpcurve")))
  expect_true(is.function(utils::getS3method("summary", "mpcurve")))
  expect_true(is.function(utils::getS3method("plot", "mpcurve")))
})

test_that("fit_mpcurve keeps algorithm formal but only accepts cavi", {
  fit_formals <- names(formals(fit_mpcurve))
  do_formals <- names(formals(do_mpcurve))

  expect_true("algorithm" %in% fit_formals)
  expect_true("greedy" %in% fit_formals)
  expect_true("S" %in% fit_formals)
  expect_true("position_prior" %in% fit_formals)
  expect_true("position_prior_init" %in% fit_formals)
  expect_true("partition_prior" %in% fit_formals)
  expect_true("partition_prior_init" %in% fit_formals)
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

test_that("fit_mpcurve rejects legacy backend requests explicitly", {
  sim <- simulate_cavi_toy(
    n = 60,
    d = 6,
    K = 4,
    rw_q = 2,
    seed = 13
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
      algorithm = "smooth_em",
      method = "PCA",
      K = 4,
      iter = 2
    ),
    "CAVI-only public wrapper"
  )
})

test_that("print.mpcurve is shorter than summary for single and partition fits", {
  sim_single <- simulate_cavi_toy(
    n = 60,
    d = 8,
    K = 5,
    rw_q = 2,
    seed = 101
  )
  fit_single <- fit_mpcurve(
    sim_single$X,
    method = "PCA",
    K = 5,
    iter = 4,
    tol = 0,
    verbose = FALSE
  )

  single_print <- paste(capture.output(print(fit_single)), collapse = "\n")
  single_summary <- paste(capture.output(summary(fit_single)), collapse = "\n")
  single_summary_obj <- summary(fit_single)
  expect_lt(length(strsplit(single_print, "\n", fixed = TRUE)[[1]]),
            length(strsplit(single_summary, "\n", fixed = TRUE)[[1]]))
  expect_true(all(c("priors", "converged", "underlying") %in% names(single_summary_obj)))
  expect_snapshot_output(print(fit_single))

  sim_partition <- simulate_dual_trajectory(
    n = 60,
    d1 = 4,
    d2 = 4,
    d_noise = 0,
    sigma = 0.15,
    seed = 102
  )
  fit_partition <- suppressWarnings(fit_mpcurve(
    sim_partition$X,
    intrinsic_dim = 2,
    method = "PCA",
    K = 6,
    n_outer = 2L,
    inner_iter = 1L,
    max_converge_iter = 2L,
    verbose = FALSE
  ))

  partition_print <- paste(capture.output(print(fit_partition)), collapse = "\n")
  partition_summary <- paste(capture.output(summary(fit_partition)), collapse = "\n")
  partition_summary_obj <- summary(fit_partition)
  expect_lt(length(strsplit(partition_print, "\n", fixed = TRUE)[[1]]),
            length(strsplit(partition_summary, "\n", fixed = TRUE)[[1]]))
  expect_true(all(c(
    "priors",
    "partition",
    "objective_history",
    "convergence_info",
    "greedy_selection"
  ) %in% names(partition_summary_obj)))
  expect_snapshot_output(print(fit_partition))
})
