devtools::load_all(".")

simulate_dual_funbank <- function(fun_bank, n = 120, d1 = 8, d2 = 8, sigma = 0.05, seed = 1L) {
  set.seed(seed)
  t1 <- runif(n)
  t2 <- runif(n)
  fun_names <- names(fun_bank)

  make_block <- function(t, d) {
    idx <- sample(seq_along(fun_bank), d, replace = TRUE)
    amp <- runif(d, 0.8, 1.5)
    shift <- runif(d, -0.4, 0.4)
    out <- sapply(seq_len(d), function(j) {
      shift[j] + amp[j] * fun_bank[[idx[j]]](t) + rnorm(length(t), 0, sigma)
    })
    list(X = matrix(out, nrow = length(t), ncol = d), fun_type = fun_names[idx])
  }

  A <- make_block(t1, d1)
  B <- make_block(t2, d2)
  X_full <- cbind(A$X, B$X)
  truth <- c(rep("A", d1), rep("B", d2))
  perm <- sample(ncol(X_full))

  colnames(X_full) <- paste0("V", seq_len(ncol(X_full)))
  list(
    X = X_full[, perm, drop = FALSE],
    true_assign = truth[perm],
    t1 = t1,
    t2 = t2
  )
}

mixed_bank <- list(
  linear = function(t) t,
  sqrt = function(t) sqrt(pmax(t, 0)),
  log = function(t) log1p(4 * t) / log(5),
  exp = function(t) (exp(2 * t) - 1) / (exp(2) - 1),
  logistic = function(t) plogis(8 * (t - 0.5)),
  quad_u = function(t) 4 * (t - 0.5)^2,
  quad_hump = function(t) 1 - 4 * (t - 0.5)^2,
  quad_shift = function(t) pmax(0, 1 - 6 * (t - 0.3)^2)
)

score_fit <- function(fit, X_sub, t_true) {
  mu <- fit$posterior$mean
  pt <- max.col(fit$gamma, ties.method = "first")
  rho <- suppressWarnings(cor(pt, t_true, method = "spearman"))
  if (!is.finite(rho)) rho <- 0
  if (rho < 0) {
    mu <- mu[, ncol(mu):1, drop = FALSE]
    rho <- -rho
  }

  K <- ncol(mu)
  bins <- cut(rank(t_true, ties.method = "first"), breaks = K, labels = FALSE)
  true_curve <- t(vapply(seq_len(ncol(X_sub)), function(j) {
    out <- numeric(K)
    for (k in seq_len(K)) out[k] <- mean(X_sub[bins == k, j])
    out
  }, numeric(K)))

  curve_cor <- mean(vapply(seq_len(nrow(mu)), function(j) {
    suppressWarnings(cor(mu[j, ], true_curve[j, ]))
  }, numeric(1)), na.rm = TRUE)

  data.frame(
    rho = rho,
    curve_cor = curve_cor,
    elbo_last = tail(fit$elbo_trace, 1),
    loglik_last = tail(fit$loglik_trace, 1),
    iter = fit$iter,
    converged = isTRUE(fit$converged)
  )
}

run_method <- function(X_sub, t_true, method, K = 8L) {
  fit <- tryCatch(
    cavi(
      X_sub,
      K = K,
      method = method,
      rw_q = 2,
      max_iter = 150,
      tol = 1e-8,
      verbose = FALSE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      method = method,
      rho = NA_real_,
      curve_cor = NA_real_,
      elbo_last = NA_real_,
      loglik_last = NA_real_,
      iter = NA_integer_,
      converged = FALSE,
      error = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }

  cbind(
    data.frame(method = method, error = NA_character_, stringsAsFactors = FALSE),
    score_fit(fit, X_sub, t_true)
  )
}

run_seed <- function(seed, methods = c("PCA", "fiedler", "pcurve", "isomap"), K = 8L) {
  sim <- simulate_dual_funbank(mixed_bank, n = 120, d1 = 8, d2 = 8, sigma = 0.05, seed = seed)
  X <- sim$X
  truth <- sim$true_assign

  XA <- X[, truth == "A", drop = FALSE]
  XB <- X[, truth == "B", drop = FALSE]

  outA <- do.call(rbind, lapply(methods, function(m) run_method(XA, sim$t1, m, K = K)))
  outA$block <- "A"
  outA$seed <- seed

  outB <- do.call(rbind, lapply(methods, function(m) run_method(XB, sim$t2, m, K = K)))
  outB$block <- "B"
  outB$seed <- seed

  rbind(outA, outB)
}

res <- do.call(rbind, lapply(1:5, run_seed))
print(res, row.names = FALSE)

cat("\nAverage by block/method:\n")
agg <- aggregate(
  cbind(rho, curve_cor, elbo_last, loglik_last, iter) ~ block + method,
  data = res,
  FUN = function(x) mean(x, na.rm = TRUE)
)
print(agg, row.names = FALSE)

cat("\nBest method by mean curve correlation within each block:\n")
for (blk in unique(agg$block)) {
  sub <- agg[agg$block == blk, ]
  best <- sub[which.max(sub$curve_cor), ]
  cat(sprintf(
    "block %s: %s (curve_cor=%.3f, rho=%.3f, elbo=%.2f)\n",
    blk, best$method, best$curve_cor, best$rho, best$elbo_last
  ))
}
