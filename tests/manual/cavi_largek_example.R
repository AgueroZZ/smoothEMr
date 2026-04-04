source("R/00_utils.R")
source("R/01_initialization.R")
source("R/02_prior.R")
source("R/06_cSmoothEM.R")
source("R/09_cavi.R")

args <- commandArgs(trailingOnly = TRUE)
plot_out <- if (length(args) >= 1L) args[[1L]] else NULL

score_orientation <- function(est_mu, true_mu) {
  cors_forward <- vapply(seq_len(nrow(true_mu)), function(j) {
    stats::cor(est_mu[j, ], true_mu[j, ], method = "spearman")
  }, numeric(1))

  est_flip <- est_mu[, ncol(est_mu):1, drop = FALSE]
  cors_flip <- vapply(seq_len(nrow(true_mu)), function(j) {
    stats::cor(est_flip[j, ], true_mu[j, ], method = "spearman")
  }, numeric(1))

  mean_forward <- mean(cors_forward)
  mean_flip <- mean(cors_flip)
  use_flip <- is.finite(mean_flip) && mean_flip > mean_forward

  if (use_flip) {
    list(
      flip = 1L,
      est_mu = est_flip,
      curve_cor = cors_flip,
      mean_curve_cor = mean_flip
    )
  } else {
    list(
      flip = 0L,
      est_mu = est_mu,
      curve_cor = cors_forward,
      mean_curve_cor = mean_forward
    )
  }
}

evaluate_run <- function(seed) {
  sim <- simulate_cavi_toy(
    n = 220,
    d = 24,
    K = 12,
    rw_q = 2,
    lambda_range = c(0.8, 3.5),
    sigma_range = c(0.07, 0.18),
    seed = seed
  )

  fit <- cavi(
    sim$X,
    K = 12,
    method = "PCA",
    rw_q = 2,
    max_iter = 60,
    tol = 1e-6,
    verbose = FALSE
  )

  orient <- score_orientation(fit$posterior$mean, sim$mu)
  z_hat <- max.col(fit$gamma, ties.method = "first")
  z_acc <- max(
    mean(z_hat == sim$z),
    mean((ncol(fit$gamma) + 1L - z_hat) == sim$z)
  )

  list(
    seed = seed,
    sim = sim,
    fit = fit,
    orient = orient,
    row = data.frame(
      seed = seed,
      iter = fit$iter,
      converged = fit$converged,
      elbo_start = fit$elbo_trace[1],
      elbo_end = tail(fit$elbo_trace, 1L),
      min_delta = min(diff(fit$elbo_trace)),
      z_acc = z_acc,
      mean_curve_cor = mean(orient$curve_cor),
      median_curve_cor = stats::median(orient$curve_cor),
      min_curve_cor = min(orient$curve_cor),
      flip = orient$flip,
      stringsAsFactors = FALSE
    )
  )
}

seeds <- 101:105
runs <- lapply(seeds, evaluate_run)
results <- do.call(rbind, lapply(runs, `[[`, "row"))

cat("=== cavi larger-K stress check ===\n")
print(results, digits = 6, row.names = FALSE)
cat("\nsummary\n")
print(summary(results[, c("iter", "converged", "min_delta", "z_acc",
                          "mean_curve_cor", "median_curve_cor", "min_curve_cor")]))
cat(sprintf(
  "\nall runs nondecreasing ELBO: %s\n",
  if (all(results$min_delta >= -1e-8)) "TRUE" else "FALSE"
))

if (!is.null(plot_out)) {
  rep_idx <- which.max(results$mean_curve_cor)
  rep_run <- runs[[rep_idx]]
  est_mu <- rep_run$orient$est_mu
  true_mu <- rep_run$sim$mu
  feat_idx <- c(1L, ceiling(nrow(true_mu) / 2), nrow(true_mu))

  grDevices::png(filename = plot_out, width = 1260, height = 420, res = 140)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3.4, 1.2), oma = c(0, 0, 1, 0))
  for (j in feat_idx) {
    ylim_j <- range(c(true_mu[j, ], est_mu[j, ]))
    graphics::plot(
      seq_len(ncol(true_mu)), true_mu[j, ],
      type = "o", pch = 16, lwd = 2, col = "#2f6db3",
      ylim = ylim_j, xlab = "Component index", ylab = "Feature mean",
      main = sprintf("Feature %d", j)
    )
    graphics::lines(
      seq_len(ncol(est_mu)), est_mu[j, ],
      type = "o", pch = 17, lwd = 2, col = "#c7522a"
    )
    graphics::legend(
      "topright",
      legend = c("True", "Posterior mean"),
      col = c("#2f6db3", "#c7522a"),
      lty = 1, pch = c(16, 17), bty = "n", cex = 0.9
    )
  }
  graphics::mtext(
    sprintf("Representative K=12 run (seed=%d, flip=%d)", rep_run$seed, rep_run$orient$flip),
    outer = TRUE, cex = 1
  )
}
