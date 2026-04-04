source("R/00_utils.R")
source("R/01_initialization.R")
source("R/02_prior.R")
source("R/06_cSmoothEM.R")
source("R/09_cavi.R")
source("R/benchmarking_curves.R")

args <- commandArgs(trailingOnly = TRUE)
plot_out <- if (length(args) >= 1L) args[[1L]] else NULL

make_truth_bins <- function(sim, K) {
  ord <- order(sim$t)
  truth <- as.matrix(sim$truth[ord, c("x", "y")])
  bin_id <- cut(seq_along(ord), breaks = K, labels = FALSE)
  t(vapply(seq_len(K), function(k) {
    colMeans(truth[bin_id == k, , drop = FALSE])
  }, numeric(2)))
}

score_fit <- function(fit, sim, K) {
  pt_hat <- as.numeric(fit$gamma %*% seq_len(K))
  rho_forward <- suppressWarnings(stats::cor(pt_hat, sim$t, method = "spearman"))
  rho_reverse <- suppressWarnings(stats::cor((K + 1 - pt_hat), sim$t, method = "spearman"))
  use_reverse <- is.finite(rho_reverse) && abs(rho_reverse) > abs(rho_forward)
  rho_best <- if (use_reverse) rho_reverse else rho_forward

  est_mu <- t(fit$posterior$mean)
  if (use_reverse) {
    est_mu <- est_mu[K:1, , drop = FALSE]
  }

  truth_bin <- make_truth_bins(sim, K)
  cor_x <- suppressWarnings(stats::cor(est_mu[, 1], truth_bin[, 1], method = "spearman"))
  cor_y <- suppressWarnings(stats::cor(est_mu[, 2], truth_bin[, 2], method = "spearman"))

  list(
    reverse = use_reverse,
    rho_best = rho_best,
    est_mu = est_mu,
    truth_bin = truth_bin,
    cor_x = cor_x,
    cor_y = cor_y
  )
}

run_method <- function(method, discretization = "quantile", sim, K = 100L) {
  fit <- cavi(
    as.matrix(sim$obs),
    K = K,
    method = method,
    discretization = discretization,
    rw_q = 2,
    max_iter = 250,
    tol = 1e-6,
    verbose = FALSE
  )

  score <- score_fit(fit, sim, K)
  list(
    method = method,
    discretization = discretization,
    fit = fit,
    score = score,
    row = data.frame(
      method = method,
      discretization = discretization,
      K = K,
      iter = fit$iter,
      converged = fit$converged,
      min_delta = min(diff(fit$elbo_trace)),
      rho_t = score$rho_best,
      cor_x = score$cor_x,
      cor_y = score$cor_y,
      reversed = score$reverse,
      stringsAsFactors = FALSE
    )
  )
}

sim <- simulate_swiss_roll_1d_2d(
  n = 1500,
  t_range = c(1.5 * base::pi, 4.5 * base::pi),
  sigma = 0.12,
  seed = 123
)

runs <- list(
  run_method("PCA", "quantile", sim = sim, K = 100L),
  run_method("fiedler", "quantile", sim = sim, K = 100L),
  run_method("pcurve", "kmeans", sim = sim, K = 100L)
)
results <- do.call(rbind, lapply(runs, `[[`, "row"))

cat("=== cavi easier spiral, larger K ===\n")
print(results, digits = 6, row.names = FALSE)
cat(sprintf(
  "\nall tested initialisations nondecreasing ELBO: %s\n",
  if (all(results$min_delta >= -1e-8)) "TRUE" else "FALSE"
))
cat("\nnotes\n")
cat("- This uses n = 1500, K = 100, sigma = 0.12, and t_range = [1.5*pi, 4.5*pi].\n")
cat("- pcurve is run with discretization = 'kmeans' so the initial clustering still uses all 100 components.\n")

if (!is.null(plot_out)) {
  grDevices::png(filename = plot_out, width = 1200, height = 900, res = 140)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  pal <- grDevices::colorRampPalette(c("#2b6cb0", "#5cb85c", "#f6c85f", "#d9534f"))(100)
  cols <- pal[cut(sim$t, breaks = 100, labels = FALSE)]

  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 1.1))
  graphics::plot(
    sim$obs$x, sim$obs$y,
    col = cols, pch = 16, cex = 0.35,
    xlab = "x", ylab = "y",
    main = "Truth"
  )
  graphics::lines(sim$truth$x[order(sim$t)], sim$truth$y[order(sim$t)], lwd = 2.2, col = "#1f2430")

  for (run in runs) {
    score <- run$score
    graphics::plot(
      sim$obs$x, sim$obs$y,
      col = "grey85", pch = 16, cex = 0.25,
      xlab = "x", ylab = "y",
      main = sprintf("%s (%s)  rho=%.3f", run$method, run$discretization, score$rho_best)
    )
    graphics::lines(sim$truth$x[order(sim$t)], sim$truth$y[order(sim$t)], lwd = 1.8, col = "#9aa5b1")
    graphics::lines(score$est_mu[, 1], score$est_mu[, 2], lwd = 2.2, col = "#c7522a")
    graphics::points(score$est_mu[, 1], score$est_mu[, 2], pch = 16, cex = 0.4, col = "#c7522a")
    graphics::legend(
      "topleft",
      legend = c("Truth curve", "Posterior mean curve"),
      col = c("#9aa5b1", "#c7522a"),
      lty = 1, pch = c(NA, 16), bty = "n", cex = 0.82
    )
  }
}
