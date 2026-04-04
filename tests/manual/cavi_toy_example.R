source("R/00_utils.R")
source("R/01_initialization.R")
source("R/02_prior.R")
source("R/06_cSmoothEM.R")
source("R/09_cavi.R")

set.seed(1)

sim <- simulate_cavi_toy(
  n = 150,
  d = 20,
  K = 8,
  rw_q = 2,
  lambda_range = c(0.8, 3),
  sigma_range = c(0.08, 0.18),
  seed = 1
)

fit <- cavi(
  sim$X,
  K = 8,
  method = "PCA",
  rw_q = 2,
  max_iter = 50,
  tol = 1e-6,
  verbose = TRUE
)

z_hat <- max.col(fit$gamma, ties.method = "first")

cat("\n=== cavi toy example ===\n")
cat(sprintf("iter: %d\n", fit$iter))
cat(sprintf("converged: %s\n", fit$converged))
cat(sprintf("ELBO start/end: %.6f -> %.6f\n", fit$elbo_trace[1], tail(fit$elbo_trace, 1)))
cat(sprintf("ELBO min delta: %.3e\n", min(diff(fit$elbo_trace))))
cat(sprintf("direct MAP assignment accuracy: %.3f\n", mean(z_hat == sim$z)))
cat(sprintf("sigma2 range: [%.3g, %.3g]\n", min(fit$params$sigma2), max(fit$params$sigma2)))
cat(sprintf("lambda range: [%.3g, %.3g]\n", min(fit$lambda_vec), max(fit$lambda_vec)))
