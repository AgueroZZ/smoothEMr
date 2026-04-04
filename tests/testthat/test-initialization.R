test_that("fiedler_ordering auto-increases the default k until the graph is connected", {
  centers <- rbind(
    c(-8, 0),
    c(0, 0),
    c(8, 0)
  )
  X <- do.call(
    rbind,
    lapply(seq_len(nrow(centers)), function(g) {
      matrix(
        stats::rnorm(40, mean = rep(centers[g, ], each = 20), sd = 0.15),
        ncol = 2,
        byrow = FALSE
      )
    })
  )

  fit_default <- expect_no_warning(fiedler_ordering(X))
  expect_equal(fit_default$n_components, 1L)
  expect_equal(length(fit_default$keep_idx), nrow(X))
  expect_equal(sum(!is.na(fit_default$t)), nrow(X))
  expect_gt(fit_default$k_used, 15)

  expect_warning(
    fit_explicit <- fiedler_ordering(X, k = 15),
    "connected components"
  )
  expect_gt(fit_explicit$n_components, 1L)
  expect_lt(length(fit_explicit$keep_idx), nrow(X))
  expect_equal(fit_explicit$k_used, 15)
})
