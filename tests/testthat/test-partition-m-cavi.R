test_that("PCA_ordering with component parameter works", {
  set.seed(1)
  X <- matrix(rnorm(500), 100, 5)

  o1 <- PCA_ordering(X, component = 1)
  o2 <- PCA_ordering(X, component = 2)
  o3 <- PCA_ordering(X, component = 3)


  # Different components should give different orderings
  expect_true(abs(cor(o1$t, o2$t)) < 0.3)
  expect_true(abs(cor(o1$t, o3$t)) < 0.3)

  # Component 1 default
  o_default <- PCA_ordering(X)
  expect_equal(o1$t, o_default$t)

  # Out of range

  expect_error(PCA_ordering(X, component = 0))
  expect_error(PCA_ordering(X, component = 100))
})

test_that("soft_partition_cavi with M=2 produces correct structure", {
  set.seed(1)
  sim <- simulate_dual_trajectory(n = 80, d1 = 4, d2 = 4,
                                   d_noise = 0, sigma = 0.2, seed = 1)

  res <- suppressWarnings(soft_partition_cavi(
    sim$X, M = 2, K = 8,
    n_outer = 3, max_converge_iter = 3,
    verbose = FALSE
  ))

  expect_s3_class(res, "soft_partition_cavi")
  expect_equal(res$M, 2L)
  expect_equal(dim(res$pi_weights), c(ncol(sim$X), 2L))
  expect_equal(length(res$fits), 2L)
  expect_s3_class(res$fits[[1]], "cavi")
  expect_s3_class(res$fits[[2]], "cavi")
  expect_equal(length(res$assign), ncol(sim$X))
  expect_true(all(res$assign %in% c("A", "B")))
  expect_true(length(res$objective_history) > 0L)
  expect_equal(res$control$freeze_unused_ordering_threshold, 0.5)
})

test_that("fit_mpcurve uses the updated partition freeze threshold default", {
  set.seed(11)
  sim <- simulate_dual_trajectory(
    n = 60, d1 = 4, d2 = 4,
    d_noise = 0, sigma = 0.2, seed = 11
  )

  fit <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 2,
    method = "PCA",
    K = 8,
    n_outer = 2,
    max_converge_iter = 2,
    verbose = FALSE
  ))

  expect_equal(fit$fit$control$freeze_unused_ordering_threshold, 0.5)
})

test_that("single-fit feature-score decomposition matches the cavi ELBO", {
  set.seed(7)
  X <- matrix(rnorm(240), 40, 6)

  fit <- cavi(X, K = 6, max_iter = 4, verbose = FALSE)
  scored <- score_features_onefit_cavi(fit, X = X, include_prior = TRUE)

  expect_equal(
    sum(scored$feature_score) + scored$global_terms,
    tail(fit$elbo_trace, 1L),
    tolerance = 1e-8
  )
})

test_that("weighted q(U) update uses the coherent partition precision", {
  set.seed(5)
  X <- matrix(rnorm(300), 50, 6)
  fit <- cavi(X, K = 5, max_iter = 3, verbose = FALSE)

  j <- 2L
  w_j <- 0.2
  q_u <- getFromNamespace(".cavi_update_q_u_weighted", "MPCurver")(
    X = X,
    R = fit$gamma,
    sigma2 = fit$params$sigma2,
    lambda_vec = fit$lambda_vec,
    Q_K = fit$Q_K,
    feature_weights = replace(rep(1, ncol(X)), j, w_j),
    rw_q = if (is.null(fit$rw_q)) 2L else fit$rw_q
  )

  Nk <- colSums(fit$gamma)
  GX <- t(X) %*% fit$gamma
  A_expected <- fit$lambda_vec[j] * fit$Q_K +
    w_j * diag(as.numeric(Nk / fit$params$sigma2[j]), ncol(fit$gamma), ncol(fit$gamma))
  A_expected <- 0.5 * (A_expected + t(A_expected))
  S_expected <- solve(A_expected)
  m_expected <- as.numeric(S_expected %*% ((w_j * GX[j, ]) / fit$params$sigma2[j]))

  expect_equal(q_u$S_list[[j]], S_expected, tolerance = 1e-8)
  expect_equal(q_u$m_mat[j, ], m_expected, tolerance = 1e-8)
})

test_that("soft partition objective follows the coherent ELBO decomposition", {
  set.seed(11)
  X <- matrix(rnorm(320), 40, 8)

  inits <- suppressWarnings(init_m_trajectories_cavi(
    X,
    M = 3,
    methods = c("PCA", "PCA", "random"),
    pca_components = c(1L, 2L, 1L),
    K = 6,
    ridge = 1e-3,
    num_iter = 2L,
    verbose = FALSE
  ))

  step <- getFromNamespace(".soft_partition_step", "MPCurver")(
    fits = inits$fits,
    X = X,
    T_now = 1,
    inner_iter = 1L,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    freeze_unused_ordering = FALSE
  )

  obj_terms <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = step$fits,
    X = X,
    weights = step$pi_weights,
    T_now = 1,
    active_orderings = step$active_orderings,
    active_feature_pairs = step$active_feature_pairs
  )

  expect_equal(
    step$objective,
    obj_terms$objective,
    tolerance = 1e-8
  )
  weights_assignment <- obj_terms$weights_assignment
  omega_hat <- colSums(weights_assignment)
  omega_hat <- omega_hat[obj_terms$assignment_info$active_idx]
  omega_hat <- omega_hat / sum(omega_hat)
  expected_assign <- sum(
    weights_assignment[, obj_terms$assignment_info$active_idx, drop = FALSE] *
      rep(log(omega_hat), each = ncol(X))
  )
  expect_equal(obj_terms$assignment_log_prior, expected_assign, tolerance = 1e-8)
})

test_that("fixed partition_prior retains the uniform assignment term", {
  set.seed(15)
  X <- matrix(rnorm(240), 40, 6)

  inits <- suppressWarnings(init_m_trajectories_cavi(
    X,
    M = 3,
    methods = c("PCA", "PCA", "random"),
    pca_components = c(1L, 2L, 1L),
    K = 6,
    ridge = 1e-3,
    num_iter = 2L,
    verbose = FALSE
  ))

  weights <- matrix(1 / 3, nrow = ncol(X), ncol = 3)
  obj_terms <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = inits$fits,
    X = X,
    weights = weights,
    T_now = 1,
    partition_prior = "fixed"
  )

  expect_equal(obj_terms$assignment_log_prior, -ncol(X) * log(3), tolerance = 1e-8)
  expect_equal(obj_terms$assignment_info$partition_prior, "fixed")
})

test_that("dirichlet assignment prior is included in the partition objective", {
  set.seed(12)
  X <- matrix(rnorm(320), 40, 8)

  inits <- suppressWarnings(init_m_trajectories_cavi(
    X,
    M = 3,
    methods = c("PCA", "PCA", "random"),
    pca_components = c(1L, 2L, 1L),
    K = 6,
    num_iter = 2L,
    verbose = FALSE
  ))

  weights <- matrix(c(
    0.8, 0.15, 0.05,
    0.7, 0.2, 0.1,
    0.6, 0.25, 0.15,
    0.75, 0.15, 0.1,
    0.1, 0.8, 0.1,
    0.15, 0.7, 0.15,
    0.2, 0.6, 0.2,
    0.1, 0.2, 0.7
  ), ncol = 3, byrow = TRUE)

  assignment_info <- getFromNamespace(".cavi_partition_assignment_info", "MPCurver")(
    weights = weights,
    T_now = 1,
    assignment_prior = "dirichlet",
    ordering_alpha = 0.5
  )

  obj_terms <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = inits$fits,
    X = X,
    weights = weights,
    T_now = 1,
    assignment_prior = "dirichlet",
    ordering_alpha = 0.5
  )

  expect_equal(obj_terms$assignment_info$assignment_prior, "dirichlet")
  expect_equal(obj_terms$assignment_info$e_log_omega, assignment_info$e_log_omega, tolerance = 1e-8)
  expect_equal(obj_terms$assignment_info$posterior_alpha, 0.5 + colSums(weights), tolerance = 1e-8)
  expect_true(is.finite(obj_terms$assignment_info$objective))
})

test_that("soft partition validates freeze/drop controls", {
  X <- matrix(rnorm(120), 30, 4)

  expect_error(
    soft_partition_cavi(
      X,
      M = 2,
      K = 5,
      freeze_unused_ordering_threshold = -1,
      n_outer = 1,
      max_converge_iter = 1,
      verbose = FALSE
    ),
    "freeze_unused_ordering_threshold"
  )
  expect_error(
    soft_partition_cavi(
      X,
      M = 2,
      K = 5,
      freeze_unused_ordering = NA,
      n_outer = 1,
      max_converge_iter = 1,
      verbose = FALSE
    ),
    "freeze_unused_ordering"
  )
  expect_error(
    soft_partition_cavi(
      X,
      M = 2,
      K = 5,
      drop_unused_ordering = NA,
      n_outer = 1,
      max_converge_iter = 1,
      verbose = FALSE
    ),
    "drop_unused_ordering"
  )
})

test_that("inactive ordering is pre-frozen before weighted update and stays frozen", {
  set.seed(23)
  X <- matrix(rnorm(240), 40, 6)
  fit_good <- cavi(X, K = 5, ridge = 0, max_iter = 3, verbose = FALSE)
  fit_bad <- fit_good
  fit_bad$posterior$mean[] <- 1e6
  fit_bad_gamma0 <- fit_bad$gamma

  step1 <- suppressWarnings(getFromNamespace(".soft_partition_step", "MPCurver")(
    fits = list(fit_good, fit_good, fit_bad),
    X = X,
    T_now = 1,
    inner_iter = 1L,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    active_orderings = c(TRUE, TRUE, TRUE),
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    iter_index = 1L,
    ordering_labels = c("A", "B", "C")
  ))

  expect_false(step1$active_orderings[3])
  expect_true(all(step1$effective_pi_weights[, 3] == 0))
  expect_true(all(step1$pi_weights[, 3] >= 0))
  expect_equal(step1$fits[[3]]$gamma, fit_bad_gamma0, tolerance = 0)

  frozen_gamma <- step1$fits[[3]]$gamma
  frozen_weights <- step1$pi_weights[, 3]
  step2 <- suppressWarnings(getFromNamespace(".soft_partition_step", "MPCurver")(
    fits = step1$fits,
    X = X,
    T_now = 1,
    inner_iter = 1L,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    active_orderings = step1$active_orderings,
    weights_prev = step1$pi_weights,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    iter_index = 2L,
    ordering_labels = c("A", "B", "C")
  ))

  expect_equal(step2$fits[[3]]$gamma, frozen_gamma, tolerance = 0)
  expect_true(all(step2$effective_pi_weights[, 3] == 0))
  expect_equal(step2$pi_weights[, 3], frozen_weights, tolerance = 0)
  event_types <- vapply(step1$ordering_events, `[[`, character(1), "event")
  event_stages <- vapply(step1$ordering_events, `[[`, character(1), "stage")
  expect_true(any(event_types %in% c("freeze", "drop_zero_mass")))
  expect_true(any(event_stages == "pre_update"))
})

test_that("comparison objective keeps frozen C blocks but zeroes frozen U blocks", {
  set.seed(29)
  X <- matrix(rnorm(240), 40, 6)
  fit_good <- cavi(X, K = 5, ridge = 0, max_iter = 3, verbose = FALSE)
  fit_bad <- fit_good
  fit_bad$posterior$mean[] <- 1e6

  step <- suppressWarnings(getFromNamespace(".soft_partition_step", "MPCurver")(
    fits = list(fit_good, fit_good, fit_bad),
    X = X,
    T_now = 1,
    inner_iter = 1L,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    active_orderings = c(TRUE, TRUE, TRUE),
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = FALSE,
    iter_index = 1L,
    ordering_labels = c("A", "B", "C")
  ))

  obj_compare <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = step$fits,
    X = X,
    weights = step$pi_weights,
    T_now = 1,
    active_orderings = step$active_orderings,
    drop_unused_ordering = FALSE
  )
  obj_drop <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = step$fits,
    X = X,
    weights = step$effective_pi_weights,
    T_now = 1,
    active_orderings = step$active_orderings,
    drop_unused_ordering = TRUE
  )

  expect_equal(step$objective, obj_compare$objective, tolerance = 1e-8)
  expect_equal(sum(obj_compare$prior_entropy_mat[, !step$active_orderings, drop = FALSE]), 0, tolerance = 1e-12)
  expect_equal(obj_compare$cell_terms, sum(obj_compare$cell_terms_by_fit), tolerance = 1e-8)
  expect_equal(obj_drop$cell_terms, sum(obj_drop$cell_terms_by_fit[step$active_orderings]), tolerance = 1e-8)
})

test_that("feature-level freezing suppresses pairwise score blow-ups for nearly identical orderings", {
  sim <- simulate_dual_trajectory(
    n = 500,
    d1 = 10,
    d2 = 10,
    d_noise = 0,
    sigma = 0.05,
    trajectory_family = c("quadratic", "monotone"),
    seed = 1
  )
  X <- sim$X

  inits <- suppressWarnings(init_m_trajectories_cavi(
    X,
    M = 2,
    methods = c("fiedler", "PCA"),
    pca_components = c(1L, 1L),
    K = 16,
    num_iter = 5L,
    verbose = FALSE
  ))

  step_no <- suppressWarnings(getFromNamespace(".soft_partition_step", "MPCurver")(
    fits = inits$fits,
    X = X,
    T_now = 1,
    inner_iter = 1L,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    active_orderings = c(TRUE, TRUE),
    active_feature_pairs = matrix(TRUE, nrow = ncol(X), ncol = 2),
    weights_prev = matrix(0.5, nrow = ncol(X), ncol = 2),
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    freeze_feature = FALSE,
    drop_unused_ordering = FALSE,
    iter_index = 1L,
    ordering_labels = c("A", "B")
  ))
  step_yes <- suppressWarnings(getFromNamespace(".soft_partition_step", "MPCurver")(
    fits = inits$fits,
    X = X,
    T_now = 1,
    inner_iter = 1L,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    active_orderings = c(TRUE, TRUE),
    active_feature_pairs = matrix(TRUE, nrow = ncol(X), ncol = 2),
    weights_prev = matrix(0.5, nrow = ncol(X), ncol = 2),
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    freeze_feature = TRUE,
    freeze_feature_weight_threshold = 0.1,
    drop_unused_ordering = FALSE,
    iter_index = 1L,
    ordering_labels = c("A", "B")
  ))

  expect_true(length(step_yes$feature_events) > 0L)
  expect_gt(max(abs(step_no$score_mat)), 1e6)
  expect_lt(max(abs(step_yes$score_mat)), 1e4)
})

test_that("drop=false preserves frozen reported weights while drop=true uses effective weights", {
  sim <- simulate_intrinsic_trajectories(
    n = 140,
    d_signal = c(6, 6),
    d_noise = 0,
    sigma = 0.08,
    seed = 2,
    trajectory_family = c("monotone", "monotone")
  )

  fit_keep <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 6,
    partition_init = "similarity",
    method = "isomap",
    K = 25,
    rw_q = 2,
    ridge = 0,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = FALSE,
    assignment_prior = "uniform",
    verbose = FALSE
  ))
  fit_drop <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 6,
    partition_init = "similarity",
    method = "isomap",
    K = 25,
    rw_q = 2,
    ridge = 0,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = TRUE,
    assignment_prior = "uniform",
    verbose = FALSE
  ))

  frozen_keep <- !fit_keep$fit$active_orderings
  expect_true(any(frozen_keep))
  expect_true(any(fit_keep$fit$pi_weights[, frozen_keep, drop = FALSE] > 0))
  expect_true(all(fit_keep$fit$effective_pi_weights[, frozen_keep, drop = FALSE] == 0))
  expect_equal(rowSums(fit_keep$fit$pi_weights), rep(1, ncol(sim$X)), tolerance = 1e-8)
  expect_equal(fit_drop$fit$pi_weights, fit_drop$fit$effective_pi_weights, tolerance = 1e-12)
})

test_that("fixed-M frozen objective decreases with larger M on a simple linear example", {
  sim <- simulate_dual_trajectory(
    n = 120,
    d1 = 6,
    d2 = 6,
    d_noise = 0,
    sigma = 0.15,
    seed = 1,
    trajectory_family = c("linear", "linear")
  )

  fit2 <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 2,
    method = "isomap",
    K = 30,
    rw_q = 2,
    n_outer = 10,
    max_converge_iter = 15,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = FALSE,
    verbose = FALSE
  ))
  fit4 <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 4,
    method = "isomap",
    K = 30,
    rw_q = 2,
    n_outer = 10,
    max_converge_iter = 15,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = FALSE,
    verbose = FALSE
  ))

  expect_equal(sum(fit4$partition$active_orderings), 2L)
  expect_lt(tail(fit4$objective_history, 1L), tail(fit2$objective_history, 1L))
})

test_that("simple monotone similarity example concentrates objective at the true dimension", {
  sim <- simulate_intrinsic_trajectories(
    n = 140,
    d_signal = c(6, 6),
    d_noise = 0,
    sigma = 0.08,
    seed = 1,
    trajectory_family = c("monotone", "monotone")
  )

  fit2 <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 2,
    partition_init = "similarity",
    method = "isomap",
    K = 25,
    rw_q = 2,
    ridge = 0,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = FALSE,
    assignment_prior = "uniform",
    verbose = FALSE
  ))
  fit6_keep <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 6,
    partition_init = "similarity",
    method = "isomap",
    K = 25,
    rw_q = 2,
    ridge = 0,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = FALSE,
    assignment_prior = "uniform",
    verbose = FALSE
  ))
  fit6_drop <- suppressWarnings(fit_mpcurve(
    sim$X,
    intrinsic_dim = 6,
    partition_init = "similarity",
    method = "isomap",
    K = 25,
    rw_q = 2,
    ridge = 0,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = TRUE,
    assignment_prior = "uniform",
    verbose = FALSE
  ))

  expect_lt(tail(fit6_keep$objective_history, 1L), tail(fit2$objective_history, 1L))
  expect_equal(sum(fit6_keep$fit$active_orderings), 2L)
  expect_equal(fit6_drop$intrinsic_dim, 6L)
  expect_equal(fit6_drop$active_intrinsic_dim, 2L)
  expect_equal(fit6_drop$displayed_intrinsic_dim, 2L)
})

test_that("weighted cavi with unit feature weights matches standard cavi exactly", {
  set.seed(13)
  sim <- simulate_dual_trajectory(
    n = 80,
    d1 = 6,
    d2 = 6,
    d_noise = 0,
    sigma = 0.2,
    seed = 13
  )
  X <- sim$X

  fit0 <- cavi(X, K = 10, max_iter = 2, tol = 0, verbose = FALSE)
  fit_std <- do_cavi(fit0, iter = 2, tol = 0, verbose = FALSE)
  fit_w <- getFromNamespace(".do_cavi_weighted", "MPCurver")(
    object = fit0,
    data = X,
    feature_weights = rep(1, ncol(X)),
    iter = 2,
    tol = 0,
    verbose = FALSE
  )
  obj1 <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = list(fit_std),
    X = X,
    weights = matrix(1, nrow = ncol(X), ncol = 1),
    T_now = 1
  )

  expect_equal(fit_std$gamma, fit_w$gamma, tolerance = 1e-10)
  expect_equal(fit_std$posterior$mean, fit_w$posterior$mean, tolerance = 1e-10)
  expect_equal(fit_std$posterior$var, fit_w$posterior$var, tolerance = 1e-10)
  expect_equal(unname(fit_std$params$sigma2), unname(fit_w$params$sigma2), tolerance = 1e-10)
  expect_equal(fit_std$lambda_vec, fit_w$lambda_vec, tolerance = 1e-10)
  expect_equal(obj1$objective, tail(fit_std$elbo_trace, 1L), tolerance = 1e-10)
})

test_that("penalized weighted cavi stays aligned with the penalized ELBO", {
  set.seed(18)
  X <- matrix(rnorm(320), 40, 8)

  fit0 <- cavi(
    X,
    K = 6,
    lambda_sd_prior_rate = 0.7,
    max_iter = 2,
    tol = 0,
    verbose = FALSE
  )
  fit_more <- getFromNamespace(".do_cavi_weighted", "MPCurver")(
    object = fit0,
    data = X,
    feature_weights = rep(1, ncol(X)),
    iter = 2,
    tol = 0,
    verbose = FALSE
  )
  obj1 <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = list(fit_more),
    X = X,
    weights = matrix(1, nrow = ncol(X), ncol = 1),
    T_now = 1
  )

  expect_gte(min(diff(fit_more$elbo_trace)), -1e-8)
  expect_equal(obj1$objective, tail(fit_more$elbo_trace, 1L), tolerance = 1e-8)
})

test_that("NULL and zero lambda prior rates both mean no partition penalty", {
  set.seed(19)
  X <- matrix(rnorm(240), 40, 6)

  fit_base <- suppressWarnings(soft_partition_cavi(
    X,
    M = 2,
    K = 6,
    n_outer = 2,
    max_converge_iter = 2,
    lambda_sd_prior_rate = NULL,
    verbose = FALSE
  ))
  fit_zero <- suppressWarnings(soft_partition_cavi(
    X,
    M = 2,
    K = 6,
    n_outer = 2,
    max_converge_iter = 2,
    lambda_sd_prior_rate = 0,
    verbose = FALSE
  ))

  expect_equal(fit_base$objective_history, fit_zero$objective_history, tolerance = 1e-8)
  expect_equal(fit_base$pi_weights, fit_zero$pi_weights, tolerance = 1e-8)
  expect_null(fit_base$control$lambda_sd_prior_rate)
  expect_null(fit_zero$control$lambda_sd_prior_rate)
})

test_that("unused properized ordering contributes zero prior-entropy mass", {
  set.seed(17)
  X <- matrix(rnorm(80 * 12), 80, 12)

  fit_used <- cavi(
    X,
    K = 10,
    rw_q = 2,
    ridge = 1e-3,
    max_iter = 4,
    tol = 0,
    verbose = FALSE
  )
  fit_unused <- getFromNamespace(".do_cavi_weighted", "MPCurver")(
    object = fit_used,
    data = X,
    feature_weights = rep(0, ncol(X)),
    iter = 8,
    tol = 0,
    verbose = FALSE
  )

  obj_single <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = list(fit_used),
    X = X,
    weights = matrix(1, nrow = ncol(X), ncol = 1),
    T_now = 1
  )
  obj_pair <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = list(fit_used, fit_unused),
    X = X,
    weights = cbind(rep(1, ncol(X)), rep(0, ncol(X))),
    T_now = 1
  )
  unused_terms <- getFromNamespace(".cavi_partition_terms_from_fit", "MPCurver")(
    fit_unused,
    data = X
  )

  expect_equal(unused_terms$cell_terms, 0, tolerance = 1e-10)
  expect_equal(sum(unused_terms$prior_entropy), 0, tolerance = 1e-6)
  expect_equal(
    obj_pair$objective - obj_single$objective,
    0,
    tolerance = 1e-6
  )

  obj_pair_fixed <- getFromNamespace(".cavi_partition_objective_from_fits", "MPCurver")(
    fits = list(fit_used, fit_unused),
    X = X,
    weights = cbind(rep(1, ncol(X)), rep(0, ncol(X))),
    T_now = 1,
    partition_prior = "fixed"
  )
  expect_equal(
    obj_pair_fixed$objective - obj_single$objective,
    -ncol(X) * log(2),
    tolerance = 1e-6
  )
})

test_that("soft_partition_cavi with M=3 produces correct structure", {
  set.seed(1)
  X <- matrix(rnorm(900), 100, 9)

  res <- suppressWarnings(soft_partition_cavi(
    X,
    M = 3,
    K = 8,
    partition_init = "ordering_methods",
    n_outer = 3,
    max_converge_iter = 3,
    verbose = FALSE
  ))

  expect_s3_class(res, "soft_partition_cavi")
  expect_equal(res$M, 3L)
  expect_equal(dim(res$pi_weights), c(9L, 3L))
  expect_equal(length(res$fits), 3L)
  expect_true(all(res$assign %in% c("A", "B", "C")))
  expect_equal(colnames(res$pi_weights), c("A", "B", "C"))
})

test_that("fit_mpcurve with intrinsic_dim=2 works", {
  set.seed(1)
  sim <- simulate_dual_trajectory(n = 80, d1 = 4, d2 = 4,
                                   d_noise = 0, sigma = 0.2, seed = 1)

  fit <- suppressWarnings(fit_mpcurve(
    sim$X, K = 8, intrinsic_dim = 2, method = "PCA",
    n_outer = 3, max_converge_iter = 3,
    freeze_unused_ordering_threshold = 1.0,
    verbose = FALSE
  ))

  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$intrinsic_dim, 2L)
  expect_equal(length(fit$fits), 2L)
  expect_true(!is.null(fit$partition))
  expect_equal(length(fit$partition$assign), ncol(sim$X))
  expect_true(!is.null(fit$locations))
  expect_named(fit$locations, c("A", "B"))
  expect_equal(length(fit$locations$A$mean$index), nrow(sim$X))
})

test_that("fit_mpcurve uses PC2 as the default second ordering after a non-PCA init", {
  set.seed(1)
  sim <- simulate_dual_trajectory(n = 80, d1 = 4, d2 = 4,
                                   d_noise = 0, sigma = 0.2, seed = 1)

  fit <- suppressWarnings(fit_mpcurve(
    sim$X,
    K = 8,
    intrinsic_dim = 2,
    partition_init = "ordering_methods",
    method = "isomap",
    n_outer = 2,
    max_converge_iter = 2,
    freeze_unused_ordering_threshold = 1.0,
    verbose = FALSE
  ))

  expect_true(!is.null(fit$init_info))
  expect_equal(fit$init_info$A$method, "isomap")
  expect_equal(fit$init_info$B$method, "PCA")
  expect_equal(fit$init_info$B$pca_component, 2L)
})

test_that("explicit method vectors default PCA components by PCA occurrence", {
  set.seed(1)
  sim <- simulate_dual_trajectory(n = 80, d1 = 4, d2 = 4,
                                   d_noise = 0, sigma = 0.2, seed = 1)

  fit <- suppressWarnings(fit_mpcurve(
    sim$X,
    K = 8,
    intrinsic_dim = 2,
    partition_init = "ordering_methods",
    method = c("fiedler", "PCA"),
    n_outer = 2,
    max_converge_iter = 2,
    freeze_unused_ordering_threshold = 1.0,
    verbose = FALSE
  ))

  expect_true(!is.null(fit$init_info))
  expect_equal(fit$init_info$A$method, "fiedler")
  expect_equal(fit$init_info$B$method, "PCA")
  expect_equal(fit$init_info$B$pca_component, 1L)
})

test_that("fit_mpcurve with intrinsic_dim=3 and method vector works", {
  set.seed(1)
  X <- matrix(rnorm(900), 100, 9)

  fit <- suppressWarnings(fit_mpcurve(
    X,
    K = 8,
    intrinsic_dim = 3,
    partition_init = "ordering_methods",
    method = c("PCA", "PCA", "PCA"),
    n_outer = 3,
    max_converge_iter = 3,
    verbose = FALSE
  ))

  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$intrinsic_dim, 3L)
  expect_equal(length(fit$fits), 3L)
  expect_named(fit$locations, c("A", "B", "C"))
})

test_that("partition path forwards ridge to all ordering fits", {
  set.seed(19)
  X <- matrix(rnorm(720), 80, 9)

  fit <- suppressWarnings(fit_mpcurve(
    X,
    K = 8,
    intrinsic_dim = 3,
    partition_init = "ordering_methods",
    method = "PCA",
    ridge = 1e-3,
    n_outer = 2,
    max_converge_iter = 2,
    verbose = FALSE
  ))

  expect_true(all(vapply(fit$fits, function(fm) {
    isTRUE(all.equal(fm$fit$control$ridge, 1e-3))
  }, logical(1))))
})

test_that("do_mpcurve continues partition fit", {
  set.seed(1)
  sim <- simulate_dual_trajectory(n = 80, d1 = 4, d2 = 4,
                                   d_noise = 0, sigma = 0.2, seed = 1)

  fit <- suppressWarnings(fit_mpcurve(
    sim$X, K = 8, intrinsic_dim = 2, method = "PCA",
    n_outer = 3, max_converge_iter = 3,
    freeze_unused_ordering_threshold = 1.0,
    verbose = FALSE
  ))

  n_before <- length(fit$objective_history)
  fit2 <- do_mpcurve(fit, iter = 3, verbose = FALSE)

  expect_equal(length(fit2$objective_history), n_before + 3L)
  expect_s3_class(fit2, "mpcurve")
  expect_equal(fit2$intrinsic_dim, 2L)
})

test_that("partition wrappers preserve freeze controls and optional dropped view", {
  set.seed(37)
  X <- matrix(rnorm(300), 50, 6)
  fit_good <- cavi(X, K = 5, ridge = 1e-3, max_iter = 3, verbose = FALSE)
  fit_bad <- fit_good
  fit_bad$posterior$mean[] <- 1e6

  res <- suppressWarnings(soft_partition_cavi(
    X,
    M = 3,
    fits_init = list(fit_good, fit_good, fit_bad),
    K = 5,
    lambda_sd_prior_rate = 0.6,
    assignment_prior = "dirichlet",
    ordering_alpha = 0.5,
    n_outer = 1,
    max_converge_iter = 1,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = TRUE,
    verbose = FALSE
  ))

  expect_true(res$frozen_orderings[3])
  expect_true(all(res$weight_history[[length(res$weight_history)]][, 3] == 0))
  expect_equal(res$control$lambda_sd_prior_rate, 0.6)
  expect_equal(res$control$assignment_prior, "dirichlet")
  expect_equal(res$control$ordering_alpha, 0.5)
  expect_true(isTRUE(res$control$freeze_unused_ordering))
  expect_equal(res$control$freeze_unused_ordering_threshold, 1.0)
  expect_true(isTRUE(res$control$drop_unused_ordering))
  active_idx <- which(res$active_orderings)
  expect_equal(
    unname(res$assignment_posterior$alpha[active_idx]),
    0.5 + unname(colSums(res$pi_weights)[active_idx]),
    tolerance = 1e-8
  )
  expect_true(all(is.na(res$assignment_posterior$alpha[-active_idx])))
  expect_null(res$assignment_posterior$omega)
  expect_null(res$assignment_posterior$partition_prior)
  expect_null(res$assignment_posterior$assignment_prior)

  fit <- as_mpcurve(res)
  expect_equal(fit$intrinsic_dim, 3L)
  expect_equal(fit$requested_intrinsic_dim, 3L)
  expect_equal(fit$active_intrinsic_dim, 2L)
  expect_equal(fit$displayed_intrinsic_dim, 2L)
  expect_equal(fit$partition$dropped_labels, "C")
  expect_equal(fit$partition$frozen_labels, "C")
  expect_true(isTRUE(fit$partition$compacted))
  expect_true(length(fit$ordering_events) >= 1L)
  expect_equal(names(fitted_prior(fit, type = "partition")$omega), c("A", "B", "C"))

  fit_wrap <- suppressWarnings(fit_mpcurve(
    X,
    K = 5,
    intrinsic_dim = 3,
    method = "PCA",
    fits_init = list(fit_good, fit_good, fit_bad),
    lambda_sd_prior_rate = 0.6,
    assignment_prior = "dirichlet",
    ordering_alpha = 0.5,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = TRUE,
    n_outer = 1,
    max_converge_iter = 1,
    verbose = FALSE
  ))
  expect_equal(fit_wrap$intrinsic_dim, 3L)
  expect_equal(fit_wrap$requested_intrinsic_dim, 3L)
  expect_equal(fit_wrap$active_intrinsic_dim, 2L)
  expect_equal(fit_wrap$displayed_intrinsic_dim, 2L)
  expect_equal(fit_wrap$partition$dropped_labels, "C")
  expect_equal(fit_wrap$partition$frozen_labels, "C")
  expect_true(isTRUE(fit_wrap$partition$compacted))
  expect_equal(fit_wrap$fit$control$lambda_sd_prior_rate, 0.6)
  expect_equal(fit_wrap$fit$control$assignment_prior, "dirichlet")
  expect_equal(fit_wrap$fit$control$ordering_alpha, 0.5)
  expect_true(isTRUE(fit_wrap$fit$control$freeze_unused_ordering))
  expect_equal(fit_wrap$fit$control$freeze_unused_ordering_threshold, 1.0)
  expect_true(isTRUE(fit_wrap$fit$control$drop_unused_ordering))
  expect_true(is.matrix(fit_wrap$data))
  expect_null(fit_wrap$fit$assignment_posterior$omega)

  tmp_plot <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_plot)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(tmp_plot)
  }, add = TRUE)
  expect_no_error(plot(fit_wrap, dims = 1))
  grDevices::dev.off()

  fit2 <- suppressWarnings(do_mpcurve(fit, iter = 1, verbose = FALSE))
  expect_equal(fit2$intrinsic_dim, 3L)
  expect_equal(fit2$requested_intrinsic_dim, 3L)
  expect_equal(fit2$active_intrinsic_dim, 2L)
  expect_equal(fit2$displayed_intrinsic_dim, 2L)
  expect_equal(fit2$partition$dropped_labels, "C")
  expect_equal(fit2$partition$frozen_labels, "C")
  expect_true(isTRUE(fit2$partition$compacted))
  expect_equal(fit2$fit$control$lambda_sd_prior_rate, 0.6)
  expect_equal(fit2$fit$control$assignment_prior, "dirichlet")
  expect_equal(fit2$fit$control$ordering_alpha, 0.5)
  expect_true(isTRUE(fit2$fit$control$freeze_unused_ordering))
  expect_equal(fit2$fit$control$freeze_unused_ordering_threshold, 1.0)
  expect_true(isTRUE(fit2$fit$control$drop_unused_ordering))
  expect_true(all(fit2$fit$weight_history[[length(fit2$fit$weight_history)]][, 3] == 0))
})

test_that("partition initialization enforces a common K when quantile bins collapse", {
  X <- cbind(
    rep(seq_len(5), each = 10),
    rep(seq_len(5), each = 10)
  )

  init <- suppressWarnings(init_m_trajectories_cavi(
    X,
    M = 2,
    methods = c("PCA", "PCA"),
    K = 12,
    discretization = "quantile",
    num_iter = 0L,
    verbose = FALSE
  ))

  expect_equal(vapply(init$fits, function(f) length(f$params$pi), integer(1)), c(12L, 12L))
  expect_true(!is.null(init$init_info))
})

test_that("plot.mpcurve works for partition fits", {
  set.seed(1)
  sim <- simulate_dual_trajectory(n = 80, d1 = 4, d2 = 4,
                                   d_noise = 0, sigma = 0.2, seed = 1)
  fit <- suppressWarnings(fit_mpcurve(
    sim$X, K = 8, intrinsic_dim = 2, method = "PCA",
    n_outer = 3, max_converge_iter = 3,
    freeze_unused_ordering_threshold = 1.0,
    verbose = FALSE
  ))

  pdf(tempfile(fileext = ".pdf"))
  expect_no_error(plot(fit, dims = c(1, 2)))
  expect_no_error(plot(fit, dims = 1))
  expect_no_error(plot(fit, plot_type = "elbo"))
  dev.off()
})

test_that("soft_two_trajectory_cavi wrapper delegates correctly", {
  set.seed(1)
  sim <- simulate_dual_trajectory(n = 80, d1 = 4, d2 = 4,
                                   d_noise = 0, sigma = 0.2, seed = 1)

  res <- suppressWarnings(soft_two_trajectory_cavi(
    sim$X, K = 8,
    n_outer = 3, max_converge_iter = 3,
    freeze_unused_ordering_threshold = 1.0,
    verbose = FALSE
  ))

  expect_s3_class(res, "soft_partition_cavi")
  expect_equal(res$M, 2L)
  expect_true(!is.null(res$fits))
  expect_equal(length(res$fits), 2L)
  expect_equal(dim(res$pi_weights), c(ncol(sim$X), 2L))
})

test_that("fit_mpcurve intrinsic_dim=1 is backward compatible", {
  set.seed(1)
  X <- matrix(rnorm(500), 100, 5)

  fit <- fit_mpcurve(X, K = 10, iter = 10, verbose = FALSE)
  expect_s3_class(fit, "mpcurve")
  expect_equal(fit$intrinsic_dim, 1L)
  expect_true(!is.null(fit$gamma))
  expect_true(!is.null(fit$params$mu))
})

test_that("partition convergence helper is stricter than the old one-step 1e-4 rule", {
  init_info <- getFromNamespace(".partition_init_convergence_info", "MPCurver")
  update_info <- getFromNamespace(".partition_update_convergence_info", "MPCurver")

  info <- init_info(tol_outer = 1e-5)
  obj_prev <- 12775.909030
  deltas <- c(3.470178, 3.673248, 3.402639, 2.690425, 2.005244, 1.541720, 1.276872)

  for (delta in deltas) {
    obj_now <- obj_prev + delta
    info <- update_info(info, obj_prev = obj_prev, obj_now = obj_now)
    obj_prev <- obj_now
  }

  expect_false(info$converged)
  expect_gt(info$last_rel_delta, 1e-5)
  expect_lt(info$last_rel_delta, 1e-4)
  expect_equal(info$consecutive_small_steps, 0L)
  expect_equal(info$phase2_iters, length(deltas))
})

test_that("partition convergence info and frozen wording are exposed to users", {
  set.seed(71)
  X <- matrix(rnorm(300), 50, 6)
  fit_good <- cavi(X, K = 5, ridge = 1e-3, max_iter = 3, verbose = FALSE)
  fit_bad <- fit_good
  fit_bad$posterior$mean[] <- 1e6

  fit <- suppressWarnings(fit_mpcurve(
    X,
    K = 5,
    intrinsic_dim = 3,
    method = "PCA",
    fits_init = list(fit_good, fit_good, fit_bad),
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 1.0,
    drop_unused_ordering = FALSE,
    n_outer = 1,
    max_converge_iter = 3,
    verbose = FALSE
  ))

  expect_true(!is.null(fit$convergence_info))
  expect_true(!is.null(fit$fit$convergence_info))
  expect_equal(fit$fit$convergence_info$consecutive_required, 3L)
  expect_equal(fit$fit$convergence_info$min_phase2_iter, 3L)
  expect_true("C" %in% fit$partition$frozen_labels)
  expect_identical(fit$partition$dropped_labels, character(0))

  printed <- capture.output(print(fit))
  expect_true(any(grepl("Frozen", printed, fixed = TRUE)))
  expect_false(any(grepl("Dropped        :", printed, fixed = TRUE)))

  summarized <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Reason", summarized, fixed = TRUE)))
})

test_that("fit_mpcurve errors for non-cavi partition", {
  X <- matrix(rnorm(200), 50, 4)
  expect_error(
    fit_mpcurve(X, algorithm = "csmooth_em", intrinsic_dim = 2),
    "CAVI-only public wrapper"
  )
})
