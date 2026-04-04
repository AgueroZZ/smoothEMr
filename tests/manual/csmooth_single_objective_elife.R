suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

x <- read.csv("elife-61271-fig2-data1-v2.csv", check.names = FALSE)
fitness_cols <- grep("_fitness$", names(x), value = TRUE)
X <- as.matrix(x[, fitness_cols, drop = FALSE])
mode(X) <- "numeric"
colnames(X) <- sub("_fitness$", "", fitness_cols)

methods_single <- c("PCA", "fiedler", "pcurve", "tSNE", "random")
res_single <- lapply(methods_single, function(m) {
  fit <- fit_mpcurve(
    X,
    algorithm = "csmooth_em",
    method = m,
    K = 30,
    rw_q = 2,
    iter = 100,
    adaptive = "ml",
    verbose = FALSE
  )
  ml_last <- if (!is.null(fit$fit$ml_trace) && length(fit$fit$ml_trace)) {
    tail(fit$fit$ml_trace, 1)
  } else {
    NA_real_
  }
  data.frame(
    method = m,
    ml_last = ml_last,
    elbo_last = tail(fit$elbo_trace, 1),
    loglik_last = tail(fit$loglik_trace, 1),
    iter = fit$iter,
    stringsAsFactors = FALSE
  )
})

print(do.call(rbind, res_single), row.names = FALSE)
