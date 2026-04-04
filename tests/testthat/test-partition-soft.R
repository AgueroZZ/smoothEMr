test_that("forward greedy partition no longer crashes with refinement enabled", {
  set.seed(1)
  X <- matrix(rnorm(60 * 8), nrow = 60, ncol = 8)

  res <- forward_two_ordering_partition_csmooth(
    X,
    K = 5,
    greedy_em_refine = TRUE,
    greedy_em_max_iter = 1,
    verbose = 0
  )

  expect_length(res$coord_assign, ncol(X))
  expect_setequal(unique(res$coord_assign), c(1L, 2L))
})

test_that("soft two-trajectory output stays on the A/B interface", {
  set.seed(1)
  sim <- simulate_dual_trajectory(
    n = 30,
    d1 = 3,
    d2 = 3,
    d_noise = 1,
    sigma = 0.2,
    seed = 1
  )

  res <- soft_two_trajectory_EM(
    sim$X,
    init_method1 = "PCA",
    K = 4,
    n_outer = 1,
    inner_iter = 1,
    max_converge_iter = 1,
    verbose = FALSE
  )

  expect_equal(dim(res$pi_weights), c(ncol(sim$X), 2L))
  expect_identical(colnames(res$pi_weights), c("A", "B"))
  expect_true(all(res$assign %in% c("A", "B")))
})

test_that("soft score and weight histories are synchronised to the post-update state", {
  softmax_rows <- function(M) {
    M <- M - apply(M, 1, max)
    E <- exp(M)
    E / rowSums(E)
  }

  set.seed(2)
  sim <- simulate_dual_trajectory(
    n = 40,
    d1 = 4,
    d2 = 4,
    d_noise = 0,
    sigma = 0.15,
    seed = 2
  )

  res <- soft_two_trajectory_EM(
    sim$X,
    init_method1 = "PCA",
    K = 5,
    n_outer = 2,
    inner_iter = 1,
    max_converge_iter = 2,
    tol_outer = 0,
    verbose = FALSE
  )

  last_idx <- length(res$score_history)
  expected_last_weights <- softmax_rows(res$score_history[[last_idx]])
  colnames(expected_last_weights) <- c("A", "B")

  expect_equal(res$weight_history[[last_idx]], expected_last_weights, tolerance = 1e-10)
  expect_equal(res$pi_weights, expected_last_weights, tolerance = 1e-10)
})

test_that("weighted collapsed ml path matches the unweighted fit when all weights are one", {
  set.seed(11)
  sim <- simulate_dual_trajectory(
    n = 40,
    d1 = 4,
    d2 = 4,
    d_noise = 0,
    sigma = 0.15,
    seed = 11
  )

  fit0 <- initialize_csmoothEM(
    sim$X,
    method = "PCA",
    K = 5,
    adaptive = "ml",
    num_iter = 0,
    include.data = TRUE
  )

  ref <- do_csmoothEM_ml_collapsed(
    fit0,
    data = sim$X,
    iter = 2,
    sigma_update = "ml",
    verbose = FALSE
  )
  wres <- do_csmoothEM_weighted(
    fit0,
    data = sim$X,
    feature_weights = rep(1, ncol(sim$X)),
    iter = 2,
    adaptive = "ml",
    sigma_update = "ml",
    verbose = FALSE
  )

  expect_equal(do.call(cbind, ref$params$mu), do.call(cbind, wres$params$mu), tolerance = 1e-5)
  expect_equal(ref$params$sigma2, wres$params$sigma2, tolerance = 1e-5)
  expect_equal(ref$gamma, wres$gamma, tolerance = 1e-5)
  expect_equal(tail(ref$ml_trace, 1), tail(wres$ml_trace, 1), tolerance = 1e-5)
})

test_that("weighted collapsed ml trace is nondecreasing for fixed soft weights", {
  set.seed(12)
  sim <- simulate_dual_trajectory(
    n = 40,
    d1 = 4,
    d2 = 4,
    d_noise = 0,
    sigma = 0.15,
    seed = 12
  )

  fit0 <- initialize_csmoothEM(
    sim$X,
    method = "PCA",
    K = 5,
    adaptive = "ml",
    num_iter = 0,
    include.data = TRUE
  )

  w <- c(0.95, 0.9, 0.8, 0.65, 0.35, 0.2, 0.1, 0.05)
  res <- do_csmoothEM_weighted(
    fit0,
    data = sim$X,
    feature_weights = w,
    iter = 5,
    adaptive = "ml",
    sigma_update = "ml",
    verbose = FALSE
  )

  expect_gte(min(diff(res$ml_trace)), -1e-8)
})

test_that("simulate_dual_trajectory supports multiple trajectory families", {
  families <- c("linear", "monotone", "quadratic", "sinusoidal")

  for (fam in families) {
    sim <- simulate_dual_trajectory(
      n = 40,
      d1 = 3,
      d2 = 2,
      d_noise = 1,
      sigma = 0.05,
      seed = 10,
      trajectory_family = fam
    )

    expect_equal(dim(sim$X), c(40, 6))
    expect_equal(length(sim$true_assign), 6)
    expect_equal(sim$trajectory_family, c(fam, fam))
  }
})

test_that("simulate_intrinsic_trajectories supports multiple intrinsic dimensions", {
  latent <- cbind(
    seq(0, 1, length.out = 30),
    seq(1, 0, length.out = 30),
    rep(0.5, 30)
  )
  sim <- simulate_intrinsic_trajectories(
    n = 30,
    d_signal = c(2, 3, 1),
    d_noise = 2,
    sigma = 0.05,
    seed = 11,
    trajectory_family = c("linear", "quadratic", "monotone"),
    latent_positions = latent
  )

  expect_equal(dim(sim$X), c(30, 8))
  expect_equal(dim(sim$latent_positions), c(30, 3))
  expect_equal(sim$ordering_labels, c("A", "B", "C"))
  expect_equal(as.integer(sort(table(sim$true_assign)[c("A", "B", "C", "noise")])),
               c(1, 2, 2, 3))
  expect_equal(sim$trajectory_family, c("linear", "quadratic", "monotone"))
  expect_equal(sim$latent_positions, latent)
})

test_that("evaluate_partition does not treat noise predictions as correct", {
  res <- evaluate_partition(list(assign = "noise"), true_assign = "A")

  expect_identical(res$accuracy, 0)
  expect_identical(res$best_alignment, "direct (A=A, B=B)")
  expect_equal(unname(res$confusion_table["noise", "A"]), 1)
})

test_that("obsolete prior adaptive mode warns and public defaults moved to ml", {
  set.seed(7)
  X <- matrix(rnorm(20 * 4), nrow = 20, ncol = 4)

  expect_warning(
    initialize_csmoothEM(X, method = "PCA", K = 4, num_iter = 0, adaptive = TRUE),
    'adaptive = "ml"'
  )

  fit0 <- initialize_csmoothEM(X, method = "PCA", K = 4, num_iter = 0, adaptive = "none")
  expect_warning(
    do_csmoothEM(fit0, data = X, iter = 1, adaptive = "prior", verbose = FALSE),
    'adaptive = "prior" is obsolete'
  )

  expect_warning(
    initialize_csmoothEM(
      X,
      method = "PCA",
      K = 4,
      num_iter = 0,
      adaptive = "ml",
      sigma_update = "mstep"
    ),
    'sigma_update = "mstep" is deprecated'
  )

  expect_identical(deparse(formals(init_two_trajectories)$adaptive), '"ml"')
  expect_identical(deparse(formals(greedy_backward_filter_csmooth)$adaptive), '"ml"')
  expect_identical(deparse(formals(parallel_initial_csmoothEM)$adaptive), '"ml"')
  expect_identical(deparse(formals(optimize_initial_csmoothEM)$adaptive), '"ml"')
  expect_match(deparse(formals(do_csmoothEM)$sigma_update), 'c\\("ml", "mstep"\\)')
  expect_false("sigma_update" %in% names(formals(do_mpcurve)))
})
