#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
if (!length(file_arg)) {
  stop("This script should be run with Rscript so --file= is available.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, export_all = FALSE, helpers = FALSE, quiet = TRUE)
  cat("Loaded package with pkgload::load_all().\n\n")
} else if (requireNamespace("MPCurver", quietly = TRUE)) {
  library(MPCurver)
  cat("Loaded installed MPCurver package.\n\n")
} else {
  stop("Need either pkgload or an installed MPCurver package to run these experiments.")
}

run_experiment <- function(name, expr) {
  cat("==", name, "==\n")
  started <- Sys.time()
  result <- tryCatch(force(expr), error = function(e) e)
  elapsed <- round(as.numeric(Sys.time() - started, units = "secs"), 2)

  if (inherits(result, "error")) {
    cat("Status: FAIL\n")
    cat("Error :", conditionMessage(result), "\n\n")
    return(invisible(FALSE))
  }

  cat("Status: PASS\n")
  if (!is.null(result) && length(result)) {
    print(result)
  }
  cat("Elapsed:", elapsed, "sec\n\n")
  invisible(TRUE)
}

ok <- TRUE

ok <- run_experiment("Experiment 1: forward greedy default path no longer crashes", {
  set.seed(1)
  X <- matrix(rnorm(60 * 8), nrow = 60, ncol = 8)

  res <- forward_two_ordering_partition_csmooth(
    X,
    K = 5,
    greedy_em_refine = TRUE,
    greedy_em_max_iter = 1,
    verbose = 0
  )

  stopifnot(length(res$coord_assign) == ncol(X))
  stopifnot(setequal(unique(res$coord_assign), c(1L, 2L)))

  data.frame(
    n_features = ncol(X),
    group1 = length(res$J1),
    group2 = length(res$J2)
  )
}) && ok

ok <- run_experiment("Experiment 2: backward greedy still runs after cache invalidation changes", {
  set.seed(2)
  sim <- simulate_dual_trajectory(
    n = 50,
    d1 = 4,
    d2 = 4,
    d_noise = 2,
    sigma = 0.2,
    seed = 2
  )

  res <- backward_two_ordering_partition_csmooth(
    sim$X,
    K = 5,
    init_method1 = "PCA",
    warm_iter_init = 2,
    warm_iter_refit = 1,
    max_steps = 4,
    verbose = FALSE
  )

  stopifnot(length(res$assign) == ncol(sim$X))
  stopifnot(all(res$assign %in% c("A", "B")))

  data.frame(
    n_features = ncol(sim$X),
    groupA = length(res$J1),
    groupB = length(res$J2)
  )
}) && ok

ok <- run_experiment("Experiment 3: soft two-trajectory output stays on the A/B interface", {
  set.seed(3)
  sim <- simulate_dual_trajectory(
    n = 50,
    d1 = 4,
    d2 = 4,
    d_noise = 2,
    sigma = 0.2,
    seed = 3
  )

  res <- soft_two_trajectory_EM(
    sim$X,
    init_method1 = "PCA",
    K = 5,
    n_outer = 2,
    inner_iter = 1,
    max_converge_iter = 2,
    verbose = FALSE
  )

  stopifnot(identical(colnames(res$pi_weights), c("A", "B")))
  stopifnot(all(res$assign %in% c("A", "B")))

  eval <- evaluate_partition(res, sim$true_assign)

  data.frame(
    n_features = ncol(sim$X),
    weight_cols = ncol(res$pi_weights),
    mean_entropy = round(mean(-rowSums(res$pi_weights * log(res$pi_weights + 1e-300))), 4),
    accuracy_on_AB = eval$accuracy
  )
}) && ok

ok <- run_experiment("Experiment 4: evaluate_partition no longer rewards noise predictions", {
  res <- evaluate_partition(list(assign = "noise"), true_assign = "A")
  stopifnot(identical(res$accuracy, 0))

  data.frame(
    predicted = "noise",
    truth = "A",
    reported_accuracy = res$accuracy,
    alignment = res$best_alignment
  )
}) && ok

if (!ok) {
  quit(status = 1)
}

cat("All manual validation experiments passed.\n")
