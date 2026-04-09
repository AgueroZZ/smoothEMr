`%||%` <- function(a, b) if (!is.null(a)) a else b


# ============================================================
# Dual-trajectory / feature-partition helpers for the recommended CAVI path.
# The exported user-facing routines are:
#   * soft_two_trajectory_cavi()
# ============================================================

.cavi_partition_weight_floor <- function() sqrt(.Machine$double.eps)
.cavi_partition_feature_freeze_threshold_default <- function() 0.1

.cavi_resolve_K <- function(X, K) {
  if (!is.null(K)) return(as.integer(K))
  max(2L, min(50L, as.integer(floor(nrow(X) / 5))))
}

.cavi_resolve_feature_weights <- function(feature_weights, d,
                                          floor = .cavi_partition_weight_floor()) {
  if (is.null(feature_weights)) return(rep(1, d))
  w <- pmax(as.numeric(feature_weights), floor)
  if (length(w) != d) stop("feature_weights must have length d.")
  w
}

.cavi_subset_measurement_sd <- function(S, cols) {
  if (is.null(S)) return(NULL)
  if (is.matrix(S)) {
    return(as.matrix(S)[, cols, drop = FALSE])
  }
  as.numeric(S)[cols]
}

.cavi_partition_order_labels <- function(M) {
  if (M <= 26L) LETTERS[seq_len(M)] else paste0("ord", seq_len(M))
}

.cavi_default_pca_components <- function(methods) {
  methods <- as.character(methods)
  out <- rep(NA_integer_, length(methods))
  pca_idx <- which(methods == "PCA")
  if (length(pca_idx) > 0L) {
    out[pca_idx] <- seq_along(pca_idx)
  }
  out
}

.cavi_validate_cluster_linkage <- function(cluster_linkage = "single") {
  valid <- c("ward.D", "ward.D2", "single", "complete",
             "average", "mcquitty", "median", "centroid")
  cluster_linkage <- as.character(cluster_linkage)[1]
  if (is.na(cluster_linkage) || !nzchar(cluster_linkage) || !cluster_linkage %in% valid) {
    stop(
      "cluster_linkage must be one of: ",
      paste(valid, collapse = ", "),
      "."
    )
  }
  cluster_linkage
}

.cavi_constant_null_score <- function(x,
                                      sigma_min = 1e-10,
                                      sigma_max = 1e10) {
  x <- as.numeric(x)
  sigma2 <- mean((x - mean(x))^2)
  sigma2 <- pmin(pmax(sigma2, sigma_min), sigma_max)
  -0.5 * length(x) * (log(2 * base::pi * sigma2) + 1)
}

.cavi_constant_null_sigma2 <- function(x,
                                       sigma_min = 1e-10,
                                       sigma_max = 1e10) {
  x <- as.numeric(x)
  sigma2 <- mean((x - mean(x))^2)
  pmin(pmax(sigma2, sigma_min), sigma_max)
}

.cavi_constant_null_score_known <- function(x, measurement_sd) {
  x <- as.numeric(x)
  sd_vec <- pmax(as.numeric(measurement_sd), sqrt(.Machine$double.eps))
  if (length(sd_vec) == 1L) sd_vec <- rep(sd_vec, length(x))
  if (length(sd_vec) != length(x)) {
    stop("measurement_sd must have length equal to length(x).")
  }
  var_vec <- sd_vec^2
  prec_vec <- 1 / var_vec
  mu_hat <- sum(prec_vec * x) / sum(prec_vec)
  as.numeric(
    -0.5 * (
      sum(log(2 * base::pi * var_vec)) +
        sum((x - mu_hat)^2 * prec_vec)
    )
  )
}

.cavi_constant_null_sigma2_known <- function(measurement_sd) {
  mean(pmax(as.numeric(measurement_sd), sqrt(.Machine$double.eps))^2)
}

.compute_same_ordering_similarity_cor <- function(X,
                                                  metric = c("spearman", "pearson"),
                                                  use = "pairwise.complete.obs",
                                                  abs_value = TRUE,
                                                  min_feature_sd = 1e-8) {
  X <- as.matrix(X)
  metric <- match.arg(metric)
  min_feature_sd <- as.numeric(min_feature_sd)[1]
  if (!is.finite(min_feature_sd) || min_feature_sd < 0) {
    stop("min_feature_sd must be a single finite nonnegative number.")
  }

  d <- ncol(X)
  feature_names <- colnames(X) %||% paste0("V", seq_len(d))
  feature_sd <- apply(X, 2, stats::sd)
  low_variance <- !is.finite(feature_sd) | feature_sd < min_feature_sd

  if (d == 1L) {
    S <- matrix(1, nrow = 1L, ncol = 1L)
  } else {
    S <- suppressWarnings(stats::cor(X, method = metric, use = use))
    if (!is.matrix(S)) {
      S <- matrix(S, nrow = d, ncol = d)
    }
    S[!is.finite(S)] <- 0
    if (isTRUE(abs_value)) {
      S <- abs(S)
    }
    S <- pmin(pmax(S, 0), 1)
    S <- 0.5 * (S + t(S))
  }

  dimnames(S) <- list(feature_names, feature_names)
  if (any(low_variance)) {
    S[low_variance, ] <- 0
    S[, low_variance] <- 0
  }
  diag(S) <- 1

  list(
    S = S,
    distance = 1 - S,
    metric = metric,
    feature_info = data.frame(
      feature = feature_names,
      sd = as.numeric(feature_sd),
      low_variance = low_variance,
      stringsAsFactors = FALSE
    )
  )
}

.cavi_equal_width_bin_index <- function(x, K) {
  x <- as.numeric(x)
  K <- as.integer(K)[1]
  xr <- range(x, na.rm = TRUE)
  if (!is.finite(xr[1]) || !is.finite(xr[2])) {
    stop("x must contain only finite values.")
  }
  if (K < 2L) stop("K must be >= 2.")
  if (xr[1] == xr[2]) {
    return(list(
      bin_index = rep(as.integer(ceiling(K / 2)), length(x)),
      breaks = c(xr[1], rep(xr[2], K))
    ))
  }
  breaks <- seq(xr[1], xr[2], length.out = K + 1L)
  idx <- findInterval(x, vec = breaks, rightmost.closed = TRUE, all.inside = TRUE)
  idx <- pmin(pmax(as.integer(idx), 1L), K)
  list(bin_index = idx, breaks = breaks)
}

.cavi_binned_directional_stats <- function(x_cov, y, K) {
  bin_info <- .cavi_equal_width_bin_index(x_cov, K = K)
  idx <- bin_info$bin_index
  counts <- tabulate(idx, nbins = K)
  y_sum <- numeric(K)
  grouped <- rowsum(matrix(as.numeric(y), ncol = 1L), group = idx, reorder = FALSE)
  y_sum[as.integer(rownames(grouped))] <- grouped[, 1]
  list(
    bin_index = idx,
    breaks = bin_info$breaks,
    counts = as.numeric(counts),
    y_sum = as.numeric(y_sum),
    y_sq = sum(as.numeric(y)^2),
    n = length(y)
  )
}

.cavi_binned_directional_stats_known <- function(x_cov, y, measurement_sd, K) {
  bin_info <- .cavi_equal_width_bin_index(x_cov, K = K)
  idx <- bin_info$bin_index
  counts <- tabulate(idx, nbins = K)
  sd_vec <- pmax(as.numeric(measurement_sd), sqrt(.Machine$double.eps))
  if (length(sd_vec) == 1L) sd_vec <- rep(sd_vec, length(y))
  if (length(sd_vec) != length(y)) {
    stop("measurement_sd must have length equal to length(y).")
  }
  var_vec <- sd_vec^2
  prec_vec <- 1 / var_vec

  prec_sum <- numeric(K)
  wy_sum <- numeric(K)
  prec_grouped <- rowsum(matrix(prec_vec, ncol = 1L), group = idx, reorder = FALSE)
  wy_grouped <- rowsum(matrix(prec_vec * as.numeric(y), ncol = 1L), group = idx, reorder = FALSE)
  prec_sum[as.integer(rownames(prec_grouped))] <- prec_grouped[, 1]
  wy_sum[as.integer(rownames(wy_grouped))] <- wy_grouped[, 1]

  list(
    bin_index = idx,
    breaks = bin_info$breaks,
    counts = as.numeric(counts),
    prec_sum = as.numeric(prec_sum),
    wy_sum = as.numeric(wy_sum),
    wy_sq = sum((as.numeric(y)^2) * prec_vec),
    log_const = sum(log(2 * base::pi * var_vec)),
    n = length(y)
  )
}

.cavi_smooth_fit_log_evidence <- function(stats,
                                          lambda,
                                          sigma2,
                                          Q_K,
                                          logdet_Q,
                                          rank_Q) {
  K <- nrow(Q_K)
  b <- stats$y_sum / sigma2
  P <- lambda * Q_K + diag(stats$counts / sigma2, K)
  chol_info <- .cavi_chol_with_jitter(P)
  L <- chol_info$chol
  alpha <- backsolve(L, forwardsolve(t(L), b))
  logdet_P <- 2 * sum(log(diag(L)))
  as.numeric(
    -0.5 * stats$n * log(2 * base::pi * sigma2) +
      0.5 * (rank_Q * log(lambda) + logdet_Q) -
      0.5 * logdet_P -
      0.5 * stats$y_sq / sigma2 +
      0.5 * sum(b * alpha)
  )
}

.cavi_smooth_fit_log_evidence_known <- function(stats,
                                                lambda,
                                                Q_K,
                                                logdet_Q,
                                                rank_Q) {
  K <- nrow(Q_K)
  b <- stats$wy_sum
  P <- lambda * Q_K + diag(stats$prec_sum, K)
  chol_info <- .cavi_chol_with_jitter(P)
  L <- chol_info$chol
  alpha <- backsolve(L, forwardsolve(t(L), b))
  logdet_P <- 2 * sum(log(diag(L)))
  as.numeric(
    -0.5 * stats$log_const +
      0.5 * (rank_Q * log(lambda) + logdet_Q) -
      0.5 * logdet_P -
      0.5 * stats$wy_sq +
      0.5 * sum(b * alpha)
  )
}

.cavi_fit_directional_smooth_evidence <- function(x_cov,
                                                  y,
                                                  measurement_sd = NULL,
                                                  K,
                                                  rw_q = 2L,
                                                  ridge = 0,
                                                  lambda_sd_prior_rate = NULL,
                                                  lambda_mode = c("optimize", "fixed"),
                                                  lambda_value = 1,
                                                  lambda_min = 1e-10,
                                                  lambda_max = 1e10,
                                                  sigma_min = 1e-10,
                                                  sigma_max = 1e10) {
  lambda_mode <- match.arg(lambda_mode)
  lambda_value <- as.numeric(lambda_value)[1]
  if (!is.finite(lambda_value) || lambda_value <= 0) {
    stop("smooth_fit_lambda_value must be a single finite positive number.")
  }
  if (!is.finite(sigma_min) || !is.finite(sigma_max) || sigma_min <= 0 || sigma_max <= 0 || sigma_min > sigma_max) {
    stop("sigma_min/sigma_max must be positive finite numbers with sigma_min <= sigma_max.")
  }
  Q_K <- make_random_walk_precision(K = K, d = 1, q = rw_q, lambda = 1, ridge = ridge)
  Q_K <- 0.5 * (as.matrix(Q_K) + t(as.matrix(Q_K)))
  prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
  has_known_sd <- !is.null(measurement_sd)
  if (has_known_sd) {
    measurement_sd <- pmax(as.numeric(measurement_sd), sqrt(.Machine$double.eps))
    if (length(measurement_sd) == 1L) measurement_sd <- rep(measurement_sd, length(y))
    if (length(measurement_sd) != length(y)) {
      stop("measurement_sd must have length equal to length(y).")
    }
    null_sigma2 <- .cavi_constant_null_sigma2_known(measurement_sd)
    null_score <- .cavi_constant_null_score_known(y, measurement_sd = measurement_sd)
    stats <- .cavi_binned_directional_stats_known(
      x_cov = x_cov,
      y = y,
      measurement_sd = measurement_sd,
      K = K
    )
  } else {
    null_sigma2 <- .cavi_constant_null_sigma2(y, sigma_min = sigma_min, sigma_max = sigma_max)
    null_score <- .cavi_constant_null_score(y, sigma_min = sigma_min, sigma_max = sigma_max)
    stats <- .cavi_binned_directional_stats(x_cov = x_cov, y = y, K = K)
  }
  rate_val <- .lambda_sd_prior_rate_value(lambda_sd_prior_rate)

  score_at <- function(lambda, sigma2) {
    if (has_known_sd) {
      .cavi_smooth_fit_log_evidence_known(
        stats = stats,
        lambda = lambda,
        Q_K = Q_K,
        logdet_Q = prior_meta$logdet,
        rank_Q = prior_meta$rank
      )
    } else {
      .cavi_smooth_fit_log_evidence(
        stats = stats,
        lambda = lambda,
        sigma2 = sigma2,
        Q_K = Q_K,
        logdet_Q = prior_meta$logdet,
        rank_Q = prior_meta$rank
      )
    }
  }

  objective_at <- function(lambda, sigma2) {
    score <- score_at(lambda, sigma2)
    if (rate_val > 0) {
      score <- score + .lambda_sd_prior_terms(
        lambda_vec = lambda,
        rate = rate_val,
        include_constant = TRUE
      )
    }
    score
  }

  lambda_init <- pmin(pmax(lambda_value, lambda_min), lambda_max)
  sigma_init <- null_sigma2

  if (has_known_sd) {
    if (identical(lambda_mode, "fixed")) {
      log_evidence <- score_at(lambda = lambda_init, sigma2 = null_sigma2)
      return(list(
        log_evidence = as.numeric(log_evidence),
        null_score = as.numeric(null_score),
        delta = as.numeric(max(log_evidence - null_score, 0)),
        lambda = as.numeric(lambda_init),
        sigma2 = as.numeric(null_sigma2),
        null_sigma2 = as.numeric(null_sigma2),
        success = TRUE,
        breaks = stats$breaks,
        bin_counts = stats$counts
      ))
    }

    opt_lambda <- tryCatch(
      stats::optimize(
        f = function(log_lambda) {
          -objective_at(lambda = exp(log_lambda), sigma2 = null_sigma2)
        },
        interval = c(log(lambda_min), log(lambda_max))
      ),
      error = function(e) NULL
    )
    if (is.null(opt_lambda)) {
      return(list(
        log_evidence = as.numeric(null_score),
        null_score = as.numeric(null_score),
        delta = 0,
        lambda = as.numeric(lambda_init),
        sigma2 = as.numeric(null_sigma2),
        null_sigma2 = as.numeric(null_sigma2),
        success = FALSE,
        breaks = stats$breaks,
        bin_counts = stats$counts
      ))
    }

    lambda_hat <- exp(opt_lambda$minimum)
    log_evidence <- score_at(lambda = lambda_hat, sigma2 = null_sigma2)
    return(list(
      log_evidence = as.numeric(log_evidence),
      null_score = as.numeric(null_score),
      delta = as.numeric(max(log_evidence - null_score, 0)),
      lambda = as.numeric(lambda_hat),
      sigma2 = as.numeric(null_sigma2),
      null_sigma2 = as.numeric(null_sigma2),
      success = TRUE,
      breaks = stats$breaks,
      bin_counts = stats$counts
    ))
  }

  if (identical(lambda_mode, "fixed")) {
    opt_sigma <- stats::optimize(
      f = function(log_sigma2) {
        -score_at(lambda = lambda_init, sigma2 = exp(log_sigma2))
      },
      interval = c(log(sigma_min), log(sigma_max))
    )
    sigma_hat <- exp(opt_sigma$minimum)
    log_evidence <- score_at(lambda = lambda_init, sigma2 = sigma_hat)
    return(list(
      log_evidence = as.numeric(log_evidence),
      null_score = as.numeric(null_score),
      delta = as.numeric(max(log_evidence - null_score, 0)),
      lambda = as.numeric(lambda_init),
      sigma2 = as.numeric(sigma_hat),
      null_sigma2 = as.numeric(null_sigma2),
      success = TRUE,
      breaks = stats$breaks,
      bin_counts = stats$counts
    ))
  }

  opt_res <- tryCatch(
    stats::optim(
      par = c(log(lambda_init), log(sigma_init)),
      fn = function(theta) {
        lambda <- exp(theta[1])
        sigma2 <- exp(theta[2])
        -objective_at(lambda = lambda, sigma2 = sigma2)
      },
      method = "L-BFGS-B",
      lower = c(log(lambda_min), log(sigma_min)),
      upper = c(log(lambda_max), log(sigma_max))
    ),
    error = function(e) NULL
  )

  if (is.null(opt_res)) {
    return(list(
      log_evidence = as.numeric(null_score),
      null_score = as.numeric(null_score),
      delta = 0,
      lambda = as.numeric(lambda_init),
      sigma2 = as.numeric(null_sigma2),
      null_sigma2 = as.numeric(null_sigma2),
      success = FALSE,
      breaks = stats$breaks,
      bin_counts = stats$counts
    ))
  }

  lambda_hat <- exp(opt_res$par[1])
  sigma_hat <- exp(opt_res$par[2])
  log_evidence <- score_at(lambda = lambda_hat, sigma2 = sigma_hat)

  list(
    log_evidence = as.numeric(log_evidence),
    null_score = as.numeric(null_score),
    delta = as.numeric(max(log_evidence - null_score, 0)),
    lambda = as.numeric(lambda_hat),
    sigma2 = as.numeric(sigma_hat),
    null_sigma2 = as.numeric(null_sigma2),
    success = TRUE,
    breaks = stats$breaks,
    bin_counts = stats$counts
  )
}

.compute_same_ordering_similarity_smooth_fit <- function(X,
                                                         S = NULL,
                                                         K,
                                                         rw_q = 2L,
                                                         ridge = 0,
                                                         lambda_sd_prior_rate = NULL,
                                                         smooth_fit_lambda_mode = c("optimize", "fixed"),
                                                         smooth_fit_lambda_value = 1,
                                                         lambda_min = 1e-10,
                                                         lambda_max = 1e10,
                                                         min_feature_sd = 1e-8,
                                                         sigma_min = 1e-10,
                                                         sigma_max = 1e10) {
  X <- as.matrix(X)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)
  min_feature_sd <- as.numeric(min_feature_sd)[1]
  if (!is.finite(min_feature_sd) || min_feature_sd < 0) {
    stop("min_feature_sd must be a single finite nonnegative number.")
  }
  if (!is.finite(as.numeric(smooth_fit_lambda_value)[1]) ||
      as.numeric(smooth_fit_lambda_value)[1] <= 0) {
    stop("smooth_fit_lambda_value must be a single finite positive number.")
  }

  n <- nrow(X)
  d <- ncol(X)
  noise_info <- .cavi_resolve_measurement_sd(S = S, X = X,
                                             caller = ".compute_same_ordering_similarity_smooth_fit()")
  feature_names <- colnames(X) %||% paste0("V", seq_len(d))
  feature_sd <- apply(X, 2, stats::sd)
  low_variance <- !is.finite(feature_sd) | feature_sd < min_feature_sd
  if (is.null(noise_info$measurement_sd)) {
    null_sigma2 <- vapply(seq_len(d), function(k) {
      .cavi_constant_null_sigma2(X[, k], sigma_min = sigma_min, sigma_max = sigma_max)
    }, numeric(1))
    null_scores <- vapply(seq_len(d), function(k) {
      .cavi_constant_null_score(X[, k], sigma_min = sigma_min, sigma_max = sigma_max)
    }, numeric(1))
  } else {
    null_sigma2 <- vapply(seq_len(d), function(k) {
      .cavi_constant_null_sigma2_known(noise_info$sd_mat[, k])
    }, numeric(1))
    null_scores <- vapply(seq_len(d), function(k) {
      .cavi_constant_null_score_known(X[, k], measurement_sd = noise_info$sd_mat[, k])
    }, numeric(1))
  }
  lambda_default <- pmin(pmax(as.numeric(smooth_fit_lambda_value)[1], lambda_min), lambda_max)

  if (d == 1L) {
    S <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    distance <- matrix(0, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_score <- matrix(null_scores, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_delta <- matrix(0, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_success <- matrix(FALSE, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_lambda <- matrix(lambda_default, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    directional_sigma2 <- matrix(null_sigma2, nrow = 1L, ncol = 1L, dimnames = list(feature_names, feature_names))
    feature_info <- data.frame(
      feature = feature_names,
      sd = as.numeric(feature_sd),
      low_variance = low_variance,
      null_score = as.numeric(null_scores),
      null_sigma2 = as.numeric(null_sigma2),
      stringsAsFactors = FALSE
    )
    return(list(
      S = S,
      distance = distance,
      metric = "smooth_fit",
      feature_info = feature_info,
      directional_score = directional_score,
      directional_delta = directional_delta,
      directional_success = directional_success,
      directional_lambda = directional_lambda,
      directional_sigma2 = directional_sigma2
    ))
  }

  directional_score <- matrix(
    rep(null_scores, each = d),
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )
  directional_delta <- matrix(
    0,
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )
  directional_success <- matrix(
    FALSE,
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )
  directional_lambda <- matrix(
    lambda_default,
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )
  directional_sigma2 <- matrix(
    rep(null_sigma2, each = d),
    nrow = d,
    ncol = d,
    dimnames = list(feature_names, feature_names)
  )

  for (j in seq_len(d)) {
    if (low_variance[j]) next
    for (k in seq_len(d)) {
      if (j == k || low_variance[k]) next
      fit_jk <- tryCatch(
        .cavi_fit_directional_smooth_evidence(
          x_cov = X[, j],
          y = X[, k],
          measurement_sd = if (is.null(noise_info$measurement_sd)) NULL else noise_info$sd_mat[, k],
          K = K,
          rw_q = rw_q,
          ridge = ridge,
          lambda_sd_prior_rate = lambda_sd_prior_rate,
          lambda_mode = smooth_fit_lambda_mode,
          lambda_value = smooth_fit_lambda_value,
          lambda_min = lambda_min,
          lambda_max = lambda_max,
          sigma_min = sigma_min,
          sigma_max = sigma_max
        ),
        error = function(e) NULL
      )
      if (is.null(fit_jk)) {
        directional_success[j, k] <- FALSE
        next
      }
      directional_score[j, k] <- fit_jk$log_evidence
      directional_delta[j, k] <- fit_jk$delta
      directional_lambda[j, k] <- fit_jk$lambda
      directional_sigma2[j, k] <- fit_jk$sigma2
      directional_success[j, k] <- isTRUE(fit_jk$success)
    }
  }

  S_raw <- pmax(directional_delta, t(directional_delta))
  diag(S_raw) <- 0
  off_diag <- S_raw[row(S_raw) != col(S_raw)]
  max_off_diag <- if (length(off_diag)) max(off_diag, na.rm = TRUE) else 0
  if (!is.finite(max_off_diag) || max_off_diag <= 0) {
    S <- diag(d)
  } else {
    S <- S_raw / max_off_diag
    S <- pmin(pmax(S, 0), 1)
    diag(S) <- 1
  }

  dimnames(S) <- list(feature_names, feature_names)
  distance <- 1 - S
  dimnames(distance) <- dimnames(S)

  feature_info <- data.frame(
    feature = feature_names,
    sd = as.numeric(feature_sd),
    low_variance = low_variance,
    null_score = as.numeric(null_scores),
    null_sigma2 = as.numeric(null_sigma2),
    stringsAsFactors = FALSE
  )

  list(
    S = S,
    distance = distance,
    metric = "smooth_fit",
    feature_info = feature_info,
    directional_score = directional_score,
    directional_delta = directional_delta,
    directional_success = directional_success,
    directional_lambda = directional_lambda,
    directional_sigma2 = directional_sigma2
  )
}

.compute_same_ordering_similarity <- function(X,
                                              S = NULL,
                                              metric = c("spearman", "pearson", "smooth_fit"),
                                              use = "pairwise.complete.obs",
                                              abs_value = TRUE,
                                              min_feature_sd = 1e-8,
                                              K = NULL,
                                              rw_q = 2L,
                                              ridge = 0,
                                              lambda_sd_prior_rate = NULL,
                                              smooth_fit_lambda_mode = c("optimize", "fixed"),
                                              smooth_fit_lambda_value = 1,
                                              lambda_min = 1e-10,
                                              lambda_max = 1e10,
                                              sigma_min = 1e-10,
                                              sigma_max = 1e10,
                                              discretization = c("quantile", "equal", "kmeans")) {
  metric <- match.arg(metric)
  if (metric %in% c("spearman", "pearson")) {
    return(.compute_same_ordering_similarity_cor(
      X = X,
      metric = metric,
      use = use,
      abs_value = abs_value,
      min_feature_sd = min_feature_sd
    ))
  }

  if (is.null(K)) {
    stop("K must be supplied when similarity_metric = 'smooth_fit'.")
  }
  .compute_same_ordering_similarity_smooth_fit(
    X = X,
    S = S,
    K = as.integer(K)[1],
    rw_q = rw_q,
    ridge = ridge,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    smooth_fit_lambda_mode = smooth_fit_lambda_mode,
    smooth_fit_lambda_value = smooth_fit_lambda_value,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    min_feature_sd = min_feature_sd,
    sigma_min = sigma_min,
    sigma_max = sigma_max
  )
}

.cavi_canonicalize_feature_clusters <- function(cluster_assign) {
  cluster_assign <- as.integer(cluster_assign)
  raw_ids <- sort(unique(cluster_assign))
  min_feature_idx <- vapply(raw_ids, function(id) {
    min(which(cluster_assign == id))
  }, integer(1))
  cluster_order <- raw_ids[order(min_feature_idx)]
  canonical_assign <- match(cluster_assign, cluster_order)
  cluster_members <- split(seq_along(canonical_assign), canonical_assign)

  list(
    feature_cluster = canonical_assign,
    cluster_order = cluster_order,
    cluster_members = cluster_members,
    cluster_sizes = vapply(cluster_members, length, integer(1))
  )
}

.cavi_resolve_similarity_methods <- function(methods = NULL,
                                             pca_components = NULL,
                                             M) {
  valid_methods <- c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap")

  if (is.null(methods)) {
    methods <- rep("PCA", M)
  } else {
    methods <- as.character(methods)
    if (length(methods) == 1L) {
      methods <- rep(methods, M)
    } else if (length(methods) != M) {
      stop(sprintf("init_methods must have length 1 or M=%d under similarity initialization.", M))
    }
    methods <- vapply(methods, match.arg, character(1), choices = valid_methods)
  }

  if (is.null(pca_components)) {
    pca_components <- .cavi_default_pca_components(methods)
  } else {
    pca_components <- as.integer(pca_components)
    if (length(pca_components) == 1L) {
      pca_components <- rep(pca_components, M)
    } else if (length(pca_components) != M) {
      stop(sprintf("pca_components must have length 1 or M=%d under similarity initialization.", M))
    }
  }

  bad_pca <- methods == "PCA" & (is.na(pca_components) | pca_components < 1L)
  if (any(bad_pca)) {
    stop("pca_components for PCA-based methods must contain positive integers.")
  }

  list(
    methods = methods,
    pca_components = pca_components
  )
}

.cavi_similarity_subset_fit <- function(X_sub,
                                        S = NULL,
                                        K,
                                        method,
                                        pca_component = 1L,
                                        rw_q = 2L,
                                        ridge = 0,
                                        lambda_sd_prior_rate = NULL,
                                        lambda_min = 1e-10,
                                        lambda_max = 1e10,
                                        sigma_min = 1e-10,
                                        sigma_max = 1e10,
                                        max_iter = 5L,
                                        tol = 1e-6,
                                        discretization = c("quantile", "equal", "kmeans"),
                                        cluster_label = NULL,
                                        verbose = FALSE) {
  X_sub <- as.matrix(X_sub)
  method <- match.arg(method, c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"))
  pca_component <- as.integer(pca_component)[1]

  method_requested <- method
  method_used <- method
  pca_requested <- if (method == "PCA") pca_component else NA_integer_
  pca_used <- pca_requested
  fallback <- FALSE
  fallback_reason <- NULL

  if (ncol(X_sub) == 1L) {
    ordering_vec <- rank(X_sub[, 1], ties.method = "first")
    fit <- .cavi_build_from_ordering(
      X = X_sub,
      ordering_vec = ordering_vec,
      S = S,
      K = K,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = max_iter,
      tol = tol,
      discretization = discretization,
      strict_K = TRUE,
      ordering_label = cluster_label %||% "single-feature cluster",
      verbose = verbose
    )
    method_used <- "single_feature_rank"
    pca_used <- NA_integer_
    fallback <- !identical(method_requested, method_used)
    fallback_reason <- if (fallback) {
      "Cluster contained a single feature; initialized by ranking that feature directly."
    } else {
      NULL
    }

    return(list(
      fit = fit,
      method_requested = method_requested,
      method_used = method_used,
      pca_component_requested = pca_requested,
      pca_component_used = pca_used,
      fallback = fallback,
      fallback_reason = fallback_reason
    ))
  }

  fit_try <- tryCatch(
    .cavi_fit_from_method(
      X = X_sub,
      S = S,
      method = method,
      pca_component = pca_component,
      K = K,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = max_iter,
      tol = tol,
      discretization = discretization,
      strict_K = TRUE,
      verbose = verbose
    ),
    error = identity
  )

  if (inherits(fit_try, "error")) {
    fallback <- TRUE
    fallback_reason <- conditionMessage(fit_try)
    fit_try <- tryCatch(
      .cavi_fit_from_method(
        X = X_sub,
        S = S,
        method = "PCA",
        pca_component = 1L,
        K = K,
        rw_q = rw_q,
        ridge = ridge,
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        sigma_min = sigma_min,
        sigma_max = sigma_max,
        max_iter = max_iter,
        tol = tol,
        discretization = discretization,
        strict_K = TRUE,
        verbose = verbose
      ),
      error = identity
    )
    if (inherits(fit_try, "error")) {
      ordering_vec <- rank(X_sub[, 1], ties.method = "first")
      fit_try <- .cavi_build_from_ordering(
        X = X_sub,
        ordering_vec = ordering_vec,
        S = S,
        K = K,
        rw_q = rw_q,
        ridge = ridge,
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        max_iter = max_iter,
        tol = tol,
        discretization = discretization,
        strict_K = TRUE,
        ordering_label = cluster_label %||% "fallback single-feature rank",
        verbose = verbose
      )
      method_used <- "single_feature_rank"
      pca_used <- NA_integer_
      fallback_reason <- paste(
        fallback_reason,
        "Fallback PCA also failed; used the first feature rank as a final fallback."
      )
    } else {
      method_used <- "PCA"
      pca_used <- 1L
    }
  }

  list(
    fit = fit_try,
    method_requested = method_requested,
    method_used = method_used,
    pca_component_requested = pca_requested,
    pca_component_used = pca_used,
    fallback = fallback,
    fallback_reason = fallback_reason
  )
}

.cavi_init_m_trajectories_similarity <- function(X,
                                                 S = NULL,
                                                 M = 2L,
                                                 methods = NULL,
                                                 pca_components = NULL,
                                                 K = NULL,
                                                 rw_q = 2L,
                                                 ridge = 0,
                                                 lambda_sd_prior_rate = NULL,
                                                 smooth_fit_lambda_mode = c("optimize", "fixed"),
                                                 smooth_fit_lambda_value = 1,
                                                 lambda_min = 1e-10,
                                                 lambda_max = 1e10,
                                                 sigma_min = 1e-10,
                                                 sigma_max = 1e10,
                                                 discretization = c("quantile", "equal", "kmeans"),
                                                 num_iter = 0L,
                                                 similarity_metric = c("spearman", "pearson", "smooth_fit"),
                                                 cluster_linkage = "single",
                                                 similarity_min_feature_sd = 1e-8,
                                                 verbose = FALSE) {
  X <- as.matrix(X)
  M <- as.integer(M)
  d <- ncol(X)
  if (M > d) {
    stop(sprintf("similarity initialization requires M=%d to be <= ncol(X)=%d.", M, d))
  }

  discretization <- match.arg(discretization)
  similarity_metric <- match.arg(similarity_metric)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)
  cluster_linkage <- .cavi_validate_cluster_linkage(cluster_linkage)
  K_use <- .cavi_resolve_K(X, K)
  if (identical(similarity_metric, "smooth_fit") && isTRUE(M > 1L) &&
      is.finite(ridge) && ridge <= 0) {
    warning(
      paste0(
        "similarity_metric = 'smooth_fit' with ridge = 0 uses an intrinsic RW prior. ",
        "This is fine for pairwise initialization, but the directional scores are pseudo-evidence ",
        "rather than fully proper marginal likelihoods; use a small positive ridge if you want ",
        "the smoother evidence to be theoretically proper."
      ),
      call. = FALSE
    )
  }
  method_info <- .cavi_resolve_similarity_methods(
    methods = methods,
    pca_components = pca_components,
    M = M
  )

  similarity <- .compute_same_ordering_similarity(
    X = X,
    S = S,
    metric = similarity_metric,
    min_feature_sd = similarity_min_feature_sd,
    K = K_use,
    rw_q = rw_q,
    ridge = ridge,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    smooth_fit_lambda_mode = smooth_fit_lambda_mode,
    smooth_fit_lambda_value = smooth_fit_lambda_value,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    discretization = discretization
  )
  hc <- stats::hclust(stats::as.dist(similarity$distance), method = cluster_linkage)
  raw_cluster <- stats::cutree(hc, k = M)
  cluster_info <- .cavi_canonicalize_feature_clusters(raw_cluster)
  ord_labels <- .cavi_partition_order_labels(M)

  fits <- vector("list", M)
  init_info <- vector("list", M)

  for (m in seq_len(M)) {
    cols_m <- cluster_info$cluster_members[[m]]
    subset_res <- .cavi_similarity_subset_fit(
      X_sub = X[, cols_m, drop = FALSE],
      S = .cavi_subset_measurement_sd(S, cols_m),
      K = K_use,
      method = method_info$methods[m],
      pca_component = method_info$pca_components[m],
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = num_iter,
      tol = 1e-6,
      discretization = discretization,
      cluster_label = sprintf("cluster %s", ord_labels[m]),
      verbose = FALSE
    )

    # Use max_iter = 0 to preserve the subset-specialized cell ordering R.
    # Iterating on all d features without feature weights would corrupt the
    # ordering (especially with known S where all features contribute equally).
    # The partition CAVI loop (with feature weights) handles refinement.
    fits[[m]] <- cavi(
      X = X,
      K = length(subset_res$fit$params$pi),
      responsibilities_init = subset_res$fit$gamma,
      position_prior_init = colMeans(subset_res$fit$gamma),
      S = S,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = 0L,
      tol = 1e-6,
      verbose = FALSE
    )

    init_info[[m]] <- list(
      label = ord_labels[m],
      method_requested = subset_res$method_requested,
      method_used = subset_res$method_used,
      pca_component_requested = subset_res$pca_component_requested,
      pca_component_used = subset_res$pca_component_used,
      K = K_use,
      discretization = discretization,
      feature_idx = cols_m,
      cluster_size = length(cols_m),
      cluster_min_feature = min(cols_m),
      fallback = subset_res$fallback,
      fallback_reason = subset_res$fallback_reason
    )
  }

  names(init_info) <- ord_labels
  names(cluster_info$cluster_sizes) <- ord_labels
  names(cluster_info$cluster_members) <- ord_labels
  feature_cluster_named <- cluster_info$feature_cluster
  names(feature_cluster_named) <- colnames(X) %||% paste0("V", seq_len(ncol(X)))

  similarity_init <- list(
    partition_init = "similarity",
    similarity_metric = similarity_metric,
    cluster_linkage = cluster_linkage,
    similarity_min_feature_sd = as.numeric(similarity_min_feature_sd)[1],
    smooth_fit_lambda_mode = smooth_fit_lambda_mode,
    smooth_fit_lambda_value = as.numeric(smooth_fit_lambda_value)[1],
    S = similarity$S,
    distance = similarity$distance,
    feature_cluster = feature_cluster_named,
    cluster_order = cluster_info$cluster_order,
    cluster_members = cluster_info$cluster_members,
    cluster_sizes = cluster_info$cluster_sizes,
    feature_info = similarity$feature_info,
    init_methods = method_info$methods,
    pca_components = method_info$pca_components,
    hclust_order = hc$order,
    directional_score = similarity$directional_score %||% NULL,
    directional_delta = similarity$directional_delta %||% NULL,
    directional_success = similarity$directional_success %||% NULL,
    directional_lambda = similarity$directional_lambda %||% NULL,
    directional_sigma2 = similarity$directional_sigma2 %||% NULL
  )

  if (isTRUE(verbose)) {
    cat(sprintf(
      "[init_m_cavi similarity] M=%d metric=%s linkage=%s sizes=%s\n",
      M,
      similarity_metric,
      cluster_linkage,
      paste(cluster_info$cluster_sizes, collapse = ",")
    ))
  }

  list(
    fits = fits,
    init_info = init_info,
    ordering_similarity = NULL,
    similarity_init = similarity_init
  )
}

.cavi_partition_uniform_log_prior <- function(d, M) {
  -as.numeric(d) * log(as.numeric(M))
}

.warn_assignment_prior_deprecated <- local({
  warned <- FALSE
  function(caller = "soft_partition_cavi()") {
    if (warned) return(invisible(NULL))
    warned <<- TRUE
    warning(
      "`assignment_prior` is deprecated; use `partition_prior` and ",
      "`partition_prior_init` instead.",
      call. = FALSE
    )
    invisible(NULL)
  }
})

.warn_ordering_alpha_deprecated <- local({
  warned <- FALSE
  function(caller = "soft_partition_cavi()") {
    if (warned) return(invisible(NULL))
    warned <<- TRUE
    warning(
      "`ordering_alpha` is deprecated and only used by the legacy ",
      "`assignment_prior = \"dirichlet\"` compatibility path.",
      call. = FALSE
    )
    invisible(NULL)
  }
})

.validate_partition_assignment_controls <- function(partition_prior = c("adaptive", "fixed"),
                                                    partition_prior_init = NULL,
                                                    assignment_prior = NULL,
                                                    ordering_alpha = NULL,
                                                    M,
                                                    partition_prior_missing = FALSE,
                                                    caller = "soft_partition_cavi()") {
  partition_prior <- match.arg(partition_prior)
  assignment_mode <- partition_prior
  legacy_assignment_prior <- NULL
  ordering_alpha_use <- NA_real_

  if (!is.null(assignment_prior)) {
    assignment_prior <- as.character(assignment_prior)[1]
    if (!assignment_prior %in% c("uniform", "dirichlet")) {
      stop("assignment_prior must be one of \"uniform\" or \"dirichlet\".", call. = FALSE)
    }
    .warn_assignment_prior_deprecated(caller = caller)
    if (identical(assignment_prior, "uniform")) {
      if (!isTRUE(partition_prior_missing) && !identical(partition_prior, "fixed")) {
        stop(
          caller,
          ": deprecated `assignment_prior = \"uniform\"` conflicts with ",
          "`partition_prior = \"", partition_prior, "\"`.",
          call. = FALSE
        )
      }
      partition_prior <- "fixed"
      assignment_mode <- "fixed"
      legacy_assignment_prior <- "uniform"
      if (!is.null(ordering_alpha)) {
        .warn_ordering_alpha_deprecated(caller = caller)
      }
    } else {
      if (!isTRUE(partition_prior_missing) && identical(partition_prior, "fixed")) {
        stop(
          caller,
          ": deprecated `assignment_prior = \"dirichlet\"` conflicts with ",
          "`partition_prior = \"fixed\"`.",
          call. = FALSE
        )
      }
      partition_prior <- "adaptive"
      assignment_mode <- "legacy_dirichlet"
      legacy_assignment_prior <- "dirichlet"
      if (!is.null(ordering_alpha)) {
        .warn_ordering_alpha_deprecated(caller = caller)
      }
      ordering_alpha_use <- if (is.null(ordering_alpha)) 0.5 else as.numeric(ordering_alpha)[1]
      if (!is.finite(ordering_alpha_use) || ordering_alpha_use <= 0) {
        stop("ordering_alpha must be a single finite positive number.", call. = FALSE)
      }
    }
  }

  if (!is.null(partition_prior_init)) {
    init <- as.numeric(partition_prior_init)
    if (length(init) != M || any(!is.finite(init)) || any(init < 0) || sum(init) <= 0) {
      stop(
        caller,
        ": `partition_prior_init` must be a nonnegative vector of length M with positive sum.",
        call. = FALSE
      )
    }
    partition_prior_init <- init / sum(init)
  }

  if (identical(assignment_mode, "legacy_dirichlet") && !is.null(partition_prior_init)) {
    warning(
      "`partition_prior_init` is ignored by the deprecated ",
      "`assignment_prior = \"dirichlet\"` compatibility path.",
      call. = FALSE
    )
  }

  list(
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    assignment_mode = assignment_mode,
    legacy_assignment_prior = legacy_assignment_prior,
    ordering_alpha = ordering_alpha_use
  )
}

.partition_assignment_weights <- function(weights,
                                          active_orderings,
                                          active_feature_pairs = NULL) {
  weights <- as.matrix(weights)
  M <- ncol(weights)
  d <- nrow(weights)
  active_orderings <- as.logical(active_orderings)
  if (length(active_orderings) != M) {
    stop("active_orderings must have length equal to ncol(weights).")
  }
  if (is.null(active_feature_pairs)) {
    active_feature_pairs <- matrix(rep(active_orderings, each = d), nrow = d, ncol = M)
  } else {
    active_feature_pairs <- as.matrix(active_feature_pairs)
    if (!all(dim(active_feature_pairs) == c(d, M))) {
      stop("active_feature_pairs must be a d x M matrix.")
    }
  }
  active_feature_pairs[, !active_orderings] <- FALSE
  weights_use <- weights
  weights_use[!active_feature_pairs] <- 0
  weights_use[, !active_orderings] <- 0
  weights_use
}

.partition_zero_mass_ordering_drop <- function(weights,
                                               active_orderings,
                                               active_feature_pairs = NULL,
                                               iter_index = NA_integer_,
                                               T_now = NA_real_,
                                               ordering_labels = NULL,
                                               stage = c("pre_update", "post_update")) {
  stage <- match.arg(stage)
  weights_use <- .partition_assignment_weights(
    weights = weights,
    active_orderings = active_orderings,
    active_feature_pairs = active_feature_pairs
  )
  M <- ncol(weights_use)
  if (is.null(ordering_labels)) {
    ordering_labels <- .cavi_partition_order_labels(M)
  }
  active_now <- as.logical(active_orderings)
  usage_mass <- colSums(weights_use)
  drop_idx <- integer(0)
  ordering_events <- list()

  if (sum(active_now) > 1L) {
    drop_idx <- which(active_now & usage_mass <= 0)
    if (length(drop_idx) >= sum(active_now)) {
      keep_idx <- which.max(ifelse(active_now, usage_mass, -Inf))
      drop_idx <- setdiff(drop_idx, keep_idx)
    }
    for (m in drop_idx) {
      active_now[m] <- FALSE
      ordering_events[[length(ordering_events) + 1L]] <- .make_partition_ordering_event(
        event = "drop_zero_mass",
        ordering = m,
        ordering_label = ordering_labels[m],
        step = iter_index,
        T_now = T_now,
        usage_mass = usage_mass[m],
        threshold = 0,
        stage = stage
      )
    }
  }

  list(
    active_orderings = active_now,
    drop_idx = drop_idx,
    usage_mass = usage_mass,
    ordering_events = ordering_events
  )
}

.validate_partition_feature_controls <- function(freeze_feature = TRUE,
                                                 freeze_feature_weight_threshold = .cavi_partition_feature_freeze_threshold_default()) {
  if (!is.logical(freeze_feature) || length(freeze_feature) != 1L ||
      is.na(freeze_feature)) {
    stop("freeze_feature must be a single TRUE/FALSE value.")
  }
  if (!is.numeric(freeze_feature_weight_threshold) ||
      length(freeze_feature_weight_threshold) != 1L ||
      !is.finite(freeze_feature_weight_threshold) ||
      freeze_feature_weight_threshold < 0) {
    stop("freeze_feature_weight_threshold must be a single finite non-negative number.")
  }
  list(
    freeze_feature = isTRUE(freeze_feature),
    freeze_feature_weight_threshold = as.numeric(freeze_feature_weight_threshold)
  )
}

.cavi_dirichlet_kl <- function(alpha_q, alpha_p) {
  alpha_q <- as.numeric(alpha_q)
  alpha_p <- as.numeric(alpha_p)
  if (length(alpha_q) != length(alpha_p)) {
    stop("alpha_q and alpha_p must have the same length.")
  }
  aq0 <- sum(alpha_q)
  ap0 <- sum(alpha_p)
  as.numeric(
    lgamma(aq0) - sum(lgamma(alpha_q)) -
      lgamma(ap0) + sum(lgamma(alpha_p)) +
      sum((alpha_q - alpha_p) * (digamma(alpha_q) - digamma(aq0)))
  )
}

.cavi_partition_assignment_info <- function(weights,
                                            T_now = 1,
                                            partition_prior = c("adaptive", "fixed"),
                                            partition_prior_init = NULL,
                                            assignment_prior = NULL,
                                            ordering_alpha = NULL,
                                            assignment_M = NULL,
                                            active_orderings = NULL,
                                            active_feature_pairs = NULL,
                                            drop_unused_ordering = FALSE,
                                            partition_prior_missing = FALSE,
                                            caller = ".cavi_partition_assignment_info()") {
  weights <- as.matrix(weights)
  ctl <- .validate_partition_assignment_controls(
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    M = ncol(weights),
    partition_prior_missing = partition_prior_missing,
    caller = caller
  )
  M <- ncol(weights)
  d <- nrow(weights)
  if (is.null(active_orderings)) {
    active_orderings <- rep(TRUE, M)
  }
  active_orderings <- as.logical(active_orderings)
  if (length(active_orderings) != M) {
    stop("active_orderings must have length equal to ncol(weights).")
  }
  active_idx <- if (identical(ctl$assignment_mode, "legacy_dirichlet")) {
    if (isTRUE(drop_unused_ordering)) which(active_orderings) else seq_len(M)
  } else {
    which(active_orderings)
  }
  if (!length(active_idx)) {
    stop("At least one ordering must remain active.")
  }
  weights_assignment <- .partition_assignment_weights(
    weights = weights,
    active_orderings = active_orderings,
    active_feature_pairs = active_feature_pairs
  )
  weights_use <- weights_assignment[, active_idx, drop = FALSE]
  if (identical(ctl$assignment_mode, "legacy_dirichlet")) {
    if (isTRUE(drop_unused_ordering)) {
      assignment_M <- length(active_idx)
    } else if (is.null(assignment_M)) {
      assignment_M <- M
    }
  } else {
    assignment_M <- length(active_idx)
  }
  if (!is.numeric(assignment_M) || length(assignment_M) != 1L ||
      !is.finite(assignment_M) || assignment_M < 1) {
    stop("assignment_M must be a single finite number >= 1.")
  }

  z_entropy <- -sum(weights_use * log(weights_use + 1e-300))
  e_log_omega_full <- rep(0, M)
  posterior_alpha_full <- rep(NA_real_, M)
  omega_full <- rep(0, M)
  assignment_prior_label <- ctl$legacy_assignment_prior %||%
    if (identical(ctl$assignment_mode, "fixed") && is.null(ctl$partition_prior_init)) {
      "uniform"
    } else {
      ctl$partition_prior
    }

  if (identical(ctl$assignment_mode, "fixed")) {
    if (is.null(ctl$partition_prior_init)) {
      omega_use <- rep(1 / length(active_idx), length(active_idx))
    } else {
      omega_use <- ctl$partition_prior_init[active_idx]
      if (sum(omega_use) <= 0) {
        stop(
          "partition_prior_init must place positive mass on the currently active orderings.",
          call. = FALSE
        )
      }
      omega_use <- omega_use / sum(omega_use)
    }
    e_log_omega_use <- log(pmax(omega_use, 1e-300))
    e_log_omega_full[active_idx] <- e_log_omega_use
    omega_full[active_idx] <- omega_use
    expected_log_assign <- sum(weights_use * rep(e_log_omega_use, each = d))
    return(list(
      assignment_prior = assignment_prior_label,
      partition_prior = ctl$partition_prior,
      assignment_mode = ctl$assignment_mode,
      legacy_assignment_prior = ctl$legacy_assignment_prior,
      ordering_alpha = ctl$ordering_alpha,
      omega = omega_full,
      e_log_omega = e_log_omega_full,
      posterior_alpha = NULL,
      expected_log_assign = as.numeric(expected_log_assign),
      kl_q_p = 0,
      z_entropy = as.numeric(z_entropy),
      objective = as.numeric(expected_log_assign + T_now * z_entropy),
      active_idx = active_idx,
      assignment_M = as.numeric(assignment_M)
    ))
  }

  if (identical(ctl$assignment_mode, "adaptive")) {
    omega_use <- colSums(weights_use)
    omega_use <- omega_use / sum(omega_use)
    e_log_omega_use <- log(pmax(omega_use, 1e-300))
    expected_log_assign <- sum(weights_use * rep(e_log_omega_use, each = d))
    e_log_omega_full[active_idx] <- e_log_omega_use
    omega_full[active_idx] <- omega_use
    return(list(
      assignment_prior = assignment_prior_label,
      partition_prior = ctl$partition_prior,
      assignment_mode = ctl$assignment_mode,
      legacy_assignment_prior = ctl$legacy_assignment_prior,
      ordering_alpha = ctl$ordering_alpha,
      omega = omega_full,
      e_log_omega = e_log_omega_full,
      posterior_alpha = NULL,
      expected_log_assign = as.numeric(expected_log_assign),
      kl_q_p = 0,
      z_entropy = as.numeric(z_entropy),
      objective = as.numeric(expected_log_assign + T_now * z_entropy),
      active_idx = active_idx,
      assignment_M = as.numeric(assignment_M)
    ))
  }

  prior_alpha <- rep(ctl$ordering_alpha, length(active_idx))
  posterior_alpha <- prior_alpha + colSums(weights_use)
  e_log_omega_use <- digamma(posterior_alpha) - digamma(sum(posterior_alpha))
  expected_log_assign <- sum(weights_use * rep(e_log_omega_use, each = d))
  kl_q_p <- .cavi_dirichlet_kl(posterior_alpha, prior_alpha)
  omega_use <- posterior_alpha / sum(posterior_alpha)
  e_log_omega_full[active_idx] <- e_log_omega_use
  posterior_alpha_full[active_idx] <- posterior_alpha
  omega_full[active_idx] <- omega_use

  list(
    assignment_prior = assignment_prior_label,
    partition_prior = ctl$partition_prior,
    assignment_mode = ctl$assignment_mode,
    legacy_assignment_prior = ctl$legacy_assignment_prior,
    ordering_alpha = ctl$ordering_alpha,
    omega = omega_full,
    e_log_omega = e_log_omega_full,
    posterior_alpha = posterior_alpha_full,
    expected_log_assign = as.numeric(expected_log_assign),
    kl_q_p = as.numeric(kl_q_p),
    z_entropy = as.numeric(z_entropy),
    objective = as.numeric(expected_log_assign + T_now * z_entropy - kl_q_p),
    active_idx = active_idx,
    assignment_M = as.numeric(assignment_M)
  )
}

.cavi_partition_assignment_state <- function(assignment_info) {
  list(
    assignment_mode = assignment_info$assignment_mode,
    legacy_assignment_prior = assignment_info$legacy_assignment_prior,
    alpha = assignment_info$posterior_alpha,
    e_log_omega = assignment_info$e_log_omega,
    assignment_M = assignment_info$assignment_M,
    active_idx = assignment_info$active_idx
  )
}

.cavi_validate_partition_methods <- function(methods, pca_components) {
  pca_idx <- which(methods == "PCA")
  if (length(pca_idx) > 1L && anyDuplicated(pca_components[pca_idx])) {
    stop("For intrinsic_dim > 1, PCA-based partition initializations must use distinct PCA components.")
  }

  repeated_det <- unique(methods[methods %in% c("fiedler", "pcurve", "isomap") &
                                   duplicated(methods)])
  if (length(repeated_det) > 0L) {
    warning(
      "Repeated deterministic partition initialization methods can produce nearly identical orderings: ",
      paste(repeated_det, collapse = ", "),
      ". Consider mixing methods, using distinct PCA components, or supplying fits_init explicitly.",
      call. = FALSE
    )
  }

  if (sum(methods == "tSNE") > 1L) {
    warning(
      "Repeated tSNE-based partition initializations rely on stochastic variation rather than guaranteed orthogonal orderings. ",
      "If you want clearly distinct starts, prefer distinct PCA components or mixed methods.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

.cavi_warn_if_similar_orderings <- function(ordering_mat, labels,
                                            threshold = 0.98) {
  if (is.null(ordering_mat) || ncol(ordering_mat) < 2L) return(invisible(NULL))

  corr_mat <- suppressWarnings(stats::cor(ordering_mat, use = "pairwise.complete.obs"))
  if (!is.matrix(corr_mat)) return(invisible(NULL))
  diag(corr_mat) <- 0

  hit <- which(abs(corr_mat) >= threshold, arr.ind = TRUE)
  if (!nrow(hit)) return(invisible(NULL))

  hit <- hit[hit[, 1] < hit[, 2], , drop = FALSE]
  if (!nrow(hit)) return(invisible(NULL))

  msgs <- apply(hit, 1, function(idx) {
    sprintf("%s vs %s (|cor|=%.3f)",
            labels[idx[1]], labels[idx[2]],
            abs(corr_mat[idx[1], idx[2]]))
  })
  warning(
    "Some partition initializations are nearly identical: ",
    paste(msgs, collapse = "; "),
    ". This can make intrinsic_dim > 1 fits hard to interpret.",
    call. = FALSE
  )

  invisible(NULL)
}

.cavi_hard_gamma_from_cluster_rank <- function(cluster_rank, n, K) {
  gamma <- matrix(0, nrow = n, ncol = K)
  gamma[cbind(seq_len(length(cluster_rank)), cluster_rank)] <- 1
  gamma
}

.cavi_exact_k_init_from_ordering <- function(X,
                                             ordering_vec,
                                             K,
                                             discretization = c("quantile", "equal", "kmeans"),
                                             strict_K = TRUE) {
  discretization <- match.arg(discretization)

  capture_init <- function(disc) {
    warnings <- character(0)
    init <- withCallingHandlers(
      make_init_csmooth(
        X = X,
        ordering_vec = ordering_vec,
        K = K,
        modelName = "homoskedastic",
        discretization = disc,
        na_action = "drop",
        eps = 1e-12
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(init = init, warnings = unique(warnings))
  }

  primary <- capture_init(discretization)
  K_primary <- length(primary$init$pi)
  if (!strict_K || K_primary == K) {
    attr(primary$init, "discretization_used") <- discretization
    attr(primary$init, "suppressed_warnings") <- primary$warnings
    return(primary$init)
  }

  if (discretization != "equal") {
    fallback <- capture_init("equal")
    if (length(fallback$init$pi) == K) {
      attr(fallback$init, "discretization_used") <- "equal"
      attr(fallback$init, "suppressed_warnings") <- c(primary$warnings, fallback$warnings)
      attr(fallback$init, "fallback_from") <- discretization
      attr(fallback$init, "fallback_reason") <- sprintf(
        "Requested discretization '%s' produced K=%d instead of K=%d.",
        discretization, K_primary, K
      )
      return(fallback$init)
    }
  }

  stop(
    "Unable to construct a partition initialization with the requested K=",
    K, ". For intrinsic_dim > 1, all orderings must use the same K."
  )
}

.cavi_get_ordering_result <- function(X, method, pca_component = 1L) {
  method <- match.arg(method, c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"))
  if (method == "random") return(NULL)

  switch(
    method,
    PCA = PCA_ordering(X, component = pca_component),
    fiedler = fiedler_ordering(X),
    pcurve = pcurve_ordering(X),
    tSNE = tSNE_ordering(X),
    isomap = isomap_ordering(X)
  )
}

.cavi_build_from_ordering <- function(X,
                                      ordering_vec,
                                      S = NULL,
                                      K,
                                      rw_q = 2L,
                                      ridge = 0,
                                      lambda_init = 1,
                                      lambda_sd_prior_rate = NULL,
                                      lambda_min = 1e-10,
                                      lambda_max = 1e10,
                                      sigma_min = 1e-10,
                                      sigma_max = 1e10,
                                      max_iter = 5L,
                                      tol = 1e-6,
                                      discretization = c("quantile", "equal", "kmeans"),
                                      strict_K = FALSE,
                                      ordering_label = NULL,
                                      verbose = FALSE) {
  discretization <- match.arg(discretization)
  init <- .cavi_exact_k_init_from_ordering(
    X = X,
    ordering_vec = ordering_vec,
    K = K,
    discretization = discretization,
    strict_K = strict_K
  )
  disc_used <- attr(init, "discretization_used") %||% discretization
  fallback_from <- attr(init, "fallback_from")
  if (!is.null(fallback_from)) {
    prefix <- if (is.null(ordering_label) || !nzchar(ordering_label)) "" else paste0(ordering_label, ": ")
    warning(
      prefix,
      attr(init, "fallback_reason"),
      " Falling back to equal-width bins so all orderings use the same K.",
      call. = FALSE
    )
  }
  gamma0 <- .cavi_hard_gamma_from_cluster_rank(init$cluster_rank, nrow(X), length(init$pi))
  cavi_args <- list(
    X = X,
    K = length(init$pi),
    responsibilities_init = gamma0,
    position_prior_init = init$pi,
    S = S,
    lambda_init = rep(lambda_init, ncol(X)),
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    rw_q = rw_q,
    ridge = ridge,
    discretization = disc_used,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )
  if (is.null(S)) {
    cavi_args$sigma2_init <- init$sigma2
  }
  do.call(cavi, cavi_args)
}

.fit_cavi_from_seed_col <- function(X,
                                    seed_col,
                                    S = NULL,
                                    K,
                                    rw_q = 2L,
                                    ridge = 0,
                                    lambda_sd_prior_rate = NULL,
                                    lambda_min = 1e-10,
                                    lambda_max = 1e10,
                                    sigma_min = 1e-10,
                                    sigma_max = 1e10,
                                    max_iter = 5L,
                                    tol = 1e-6,
                                    discretization = c("quantile", "equal", "kmeans"),
                                    strict_K = FALSE,
                                    verbose = FALSE) {
  ordering_vec <- rank(X[, seed_col], ties.method = "first")
  .cavi_build_from_ordering(
    X = X,
    ordering_vec = ordering_vec,
    S = S,
    K = K,
    rw_q = rw_q,
    ridge = ridge,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    max_iter = max_iter,
    tol = tol,
    discretization = discretization,
    strict_K = strict_K,
    ordering_label = sprintf("seed feature %d", seed_col),
    verbose = verbose
  )
}

.cavi_rebuild_subset <- function(X_sub,
                                 S = NULL,
                                 K = NULL,
                                 rw_q = 2L,
                                 ridge = 0,
                                 gamma_init = NULL,
                                 method = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
                                 max_iter = 5L,
                                 tol = 1e-6,
                                 lambda_init = 1,
                                 lambda_sd_prior_rate = NULL,
                                 verbose = FALSE) {
  X_sub <- as.matrix(X_sub)
  method <- match.arg(method)

  if (is.null(gamma_init)) {
    return(cavi(
      X = X_sub,
      K = K,
      method = method,
      S = S,
      rw_q = rw_q,
      ridge = ridge,
      lambda_init = rep(lambda_init, ncol(X_sub)),
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose
    ))
  }

  gamma_init <- as.matrix(gamma_init)
  K_use <- ncol(gamma_init)
  cavi_args <- list(
    X = X_sub,
    K = K_use,
    responsibilities_init = gamma_init,
    position_prior_init = colMeans(gamma_init),
    S = S,
    lambda_init = rep(lambda_init, ncol(X_sub)),
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    rw_q = rw_q,
    ridge = ridge,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )
  if (is.null(S)) {
    cavi_args$sigma2_init <- pmax(apply(X_sub, 2, stats::var), 1e-10)
  }
  do.call(cavi, cavi_args)
}

.cavi_prepare_weighted_metadata <- function(K, rw_q, Q_K) {
  prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
  list(
    r_rank = prior_meta$rank,
    logdet_Q = prior_meta$logdet,
    prior_proper = prior_meta$proper
  )
}

.cavi_chol_with_jitter <- function(A,
                                   base_jitter = 1e-8,
                                   max_tries = 6L) {
  A <- 0.5 * (A + t(A))
  L <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(L)) {
    return(list(chol = L, matrix = A, jitter = 0))
  }

  jitter0 <- max(base_jitter, 1e-10 * max(1, mean(diag(A))))
  for (i in seq_len(max_tries)) {
    jitter_i <- jitter0 * (10 ^ (i - 1L))
    A_try <- A + diag(jitter_i, nrow(A))
    L <- tryCatch(chol(A_try), error = function(e) NULL)
    if (!is.null(L)) {
      return(list(chol = L, matrix = A_try, jitter = jitter_i))
    }
  }

  stop("Unable to stabilize the weighted q(U) precision matrix for Cholesky factorization.")
}

.cavi_update_q_u_weighted <- function(X, R, sigma2, lambda_vec, Q_K,
                                      measurement_sd = NULL,
                                      feature_weights, rw_q = 2L,
                                      feature_active = NULL,
                                      q_u_current = NULL) {
  X <- as.matrix(X)
  R <- as.matrix(R)
  d <- ncol(X)
  K <- ncol(R)
  w <- .cavi_resolve_feature_weights(feature_weights, d)
  if (is.null(feature_active)) {
    feature_active <- rep(TRUE, d)
  }
  feature_active <- as.logical(feature_active)
  if (length(feature_active) != d) {
    stop("feature_active must have length d.")
  }
  if (any(!feature_active) && is.null(q_u_current)) {
    stop("q_u_current must be supplied when feature_active contains FALSE.")
  }

  noise_info <- .cavi_resolve_measurement_sd(
    S = measurement_sd,
    X = X,
    caller = ".cavi_update_q_u_weighted()"
  )
  use_known_sd <- !is.null(noise_info$measurement_sd)
  if (!use_known_sd) {
    Nk <- colSums(R)
    Nk[Nk < 1e-8] <- 1e-8
    GX <- t(X) %*% R
  }

  m_mat <- matrix(0, nrow = d, ncol = K)
  sdiag_mat <- matrix(0, nrow = d, ncol = K)
  S_list <- vector("list", d)
  logdetS <- numeric(d)
  eq_quad <- numeric(d)

  for (j in seq_len(d)) {
    if (!feature_active[j]) {
      S_list[[j]] <- q_u_current$S_list[[j]]
      m_mat[j, ] <- q_u_current$m_mat[j, ]
      sdiag_mat[j, ] <- q_u_current$sdiag_mat[j, ]
      logdetS[j] <- q_u_current$logdetS[j]
      eq_quad[j] <- q_u_current$eq_quad[j]
      next
    }
    if (use_known_sd) {
      inv_vj <- w[j] / noise_info$var_mat[, j]
      diag_j <- as.numeric(crossprod(inv_vj, R))
      A_j <- diag(diag_j, K, K) + lambda_vec[j] * Q_K
      rhs_j <- as.numeric(crossprod(inv_vj * X[, j], R))
    } else {
      A_j <- w[j] * diag(as.numeric(Nk / sigma2[j]), K, K) +
        lambda_vec[j] * Q_K
      rhs_j <- (w[j] * as.numeric(GX[j, ])) / sigma2[j]
    }
    chol_info <- .cavi_chol_with_jitter(A_j)
    A_j <- chol_info$matrix
    L_j <- chol_info$chol
    S_j <- chol2inv(L_j)
    y_j <- forwardsolve(t(L_j), rhs_j)
    m_j <- as.numeric(backsolve(L_j, y_j))

    S_list[[j]] <- S_j
    m_mat[j, ] <- m_j
    sdiag_mat[j, ] <- diag(S_j)
    logdetS[j] <- -2 * sum(log(diag(L_j)))
    eq_quad[j] <- as.numeric(crossprod(m_j, Q_K %*% m_j) + sum(Q_K * S_j))
  }

  list(
    m_mat = m_mat,
    sdiag_mat = sdiag_mat,
    S_list = S_list,
    logdetS = logdetS,
    eq_quad = eq_quad,
    weights = w
  )
}

.cavi_update_r_weighted <- function(X, q_u, pi_vec, sigma2,
                                    measurement_sd = NULL,
                                    feature_weights) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  K <- length(pi_vec)
  w <- .cavi_resolve_feature_weights(feature_weights, d)
  noise_info <- .cavi_resolve_measurement_sd(
    S = measurement_sd,
    X = X,
    caller = ".cavi_update_r_weighted()"
  )

  m_mat <- q_u$m_mat
  sdiag_mat <- q_u$sdiag_mat

  log_r <- matrix(0, nrow = n, ncol = K)
  if (is.null(noise_info$measurement_sd)) {
    const_sigma <- sum(w * log(2 * base::pi * sigma2))
    inv_sigma2_w <- matrix(w / sigma2, n, d, byrow = TRUE)
  } else {
    inv_sigma2_w <- matrix(w, n, d, byrow = TRUE) / noise_info$var_mat
    const_sigma <- rowSums(matrix(w, n, d, byrow = TRUE) *
                             log(2 * base::pi * noise_info$var_mat))
  }

  for (k in seq_len(K)) {
    diff_k <- sweep(X, 2, m_mat[, k], "-")
    quad_det <- rowSums(inv_sigma2_w * diff_k^2)
    if (is.null(noise_info$measurement_sd)) {
      quad_unc <- sum((w * sdiag_mat[, k]) / sigma2)
      log_r[, k] <- log(pmax(pi_vec[k], .Machine$double.eps)) -
        0.5 * (const_sigma + quad_det + quad_unc)
    } else {
      quad_unc <- drop(inv_sigma2_w %*% sdiag_mat[, k])
      log_r[, k] <- log(pmax(pi_vec[k], .Machine$double.eps)) -
        0.5 * (const_sigma + quad_det + quad_unc)
    }
  }

  lse <- matrixStats::rowLogSumExps(log_r)
  exp(log_r - lse)
}

.cavi_update_sigma2_weighted <- function(X, R, q_u,
                                         sigma2_current = NULL,
                                         measurement_sd = NULL,
                                         sigma_min = 1e-10,
                                         sigma_max = 1e10,
                                         feature_weights = NULL,
                                         feature_active = NULL) {
  X <- as.matrix(X)
  R <- as.matrix(R)
  d <- ncol(X)
  K <- ncol(R)
  n <- nrow(X)
  w <- .cavi_resolve_feature_weights(feature_weights, d)
  if (is.null(feature_active)) {
    feature_active <- rep(TRUE, d)
  }
  feature_active <- as.logical(feature_active)
  if (length(feature_active) != d) {
    stop("feature_active must have length d.")
  }
  if (!is.null(measurement_sd)) {
    return(sigma2_current)
  }

  m_mat <- q_u$m_mat
  sdiag_mat <- q_u$sdiag_mat
  Nk <- colSums(R)

  sigma2_new <- numeric(d)
  for (j in seq_len(d)) {
    if (!feature_active[j] ||
        (!is.null(sigma2_current) && w[j] <= .cavi_partition_weight_floor())) {
      sigma2_new[j] <- sigma2_current[j]
      next
    }
    accum_j <- 0
    for (k in seq_len(K)) {
      diff2_jk <- (X[, j] - m_mat[j, k])^2
      accum_j <- accum_j + sum(R[, k] * diff2_jk) + Nk[k] * sdiag_mat[j, k]
    }
    sigma2_new[j] <- accum_j / n
  }
  pmin(pmax(sigma2_new, sigma_min), sigma_max)
}

.cavi_update_lambda_weighted <- function(q_u, r_rank,
                                         lambda_sd_prior_rate = NULL,
                                         lambda_min = 1e-10,
                                         lambda_max = 1e10,
                                         lambda_current = NULL,
                                         feature_weights = NULL,
                                         feature_active = NULL) {
  d <- nrow(q_u$m_mat)
  w <- .cavi_resolve_feature_weights(feature_weights, d)
  if (is.null(feature_active)) {
    feature_active <- rep(TRUE, d)
  }
  feature_active <- as.logical(feature_active)
  if (length(feature_active) != d) {
    stop("feature_active must have length d.")
  }
  denom <- pmax(q_u$eq_quad, 1e-12)
  lambda_sd_prior_rate_val <- .lambda_sd_prior_rate_value(lambda_sd_prior_rate)
  if (lambda_sd_prior_rate_val <= 0) {
    lambda_new <- pmin(pmax(r_rank / denom, lambda_min), lambda_max)
  } else {
    lambda_new <- .optimize_lambda_induced_exp(
      eq_quad = denom,
      r_rank = r_rank,
      rate = lambda_sd_prior_rate_val,
      lambda_min = lambda_min,
      lambda_max = lambda_max
    )
  }
  if (!is.null(lambda_current)) {
    hold <- !feature_active | w <= .cavi_partition_weight_floor()
    lambda_new[hold] <- lambda_current[hold]
  }
  lambda_new
}

.cavi_entropy_u <- function(q_u, K) {
  sum(0.5 * (K * log(2 * base::pi * exp(1)) + q_u$logdetS))
}

.cavi_entropy_u_terms <- function(q_u, K) {
  0.5 * (K * log(2 * base::pi * exp(1)) + q_u$logdetS)
}

.cavi_prior_terms_from_state <- function(q_u, lambda_vec, Q_K, rw_q = 2L,
                                         lambda_sd_prior_rate = NULL) {
  prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
  r_rank <- prior_meta$rank
  logdet_Q <- prior_meta$logdet
  0.5 * (r_rank * (log(lambda_vec) - log(2 * base::pi)) + logdet_Q) -
    0.5 * lambda_vec * q_u$eq_quad +
    .lambda_sd_prior_terms(lambda_vec, rate = lambda_sd_prior_rate, include_constant = TRUE)
}

.cavi_cell_terms_from_state <- function(R, pi_vec) {
  log_pi_term <- sum(R * rep(log(pmax(pi_vec, .Machine$double.eps)), each = nrow(R)))
  entropy_c <- -sum(R * log(R + 1e-300))
  as.numeric(log_pi_term + entropy_c)
}

.cavi_posterior_state_from_fit <- function(fit) {
  q_u <- list(
    m_mat = fit$posterior$mean,
    sdiag_mat = fit$posterior$var,
    S_list = fit$posterior$cov,
    logdetS = vapply(
      fit$posterior$cov,
      function(S) as.numeric(determinant(S, logarithm = TRUE)$modulus),
      numeric(1)
    )
  )
  q_u$eq_quad <- vapply(seq_len(nrow(q_u$m_mat)), function(j) {
    mj <- q_u$m_mat[j, ]
    Sj <- q_u$S_list[[j]]
    as.numeric(crossprod(mj, fit$Q_K %*% mj) + sum(fit$Q_K * Sj))
  }, numeric(1))
  q_u
}

.cavi_partition_terms_from_fit <- function(fit, data = NULL) {
  if (!inherits(fit, "cavi")) stop("fit must inherit from class 'cavi'.")
  X <- as.matrix(data %||% fit$data)
  if (is.null(X)) stop("X must be provided or stored in fit$data.")
  feature_active <- (fit$control %||% list())$feature_active %||% rep(TRUE, ncol(X))
  feature_active <- as.logical(feature_active)
  if (length(feature_active) != ncol(X)) {
    stop("fit$control$feature_active must have length d.")
  }

  q_u <- .cavi_posterior_state_from_fit(fit)
  K <- length(fit$params$pi)
  rw_q <- fit$rw_q %||% 2L
  lambda_sd_prior_rate <- (fit$control %||% list())$lambda_sd_prior_rate %||% NULL

  like_score <- .cavi_feature_scores_from_state(
    X = X,
    R = fit$gamma,
    q_u = q_u,
    sigma2 = fit$params$sigma2,
    measurement_sd = fit$measurement_sd %||% NULL,
    lambda_vec = fit$lambda_vec,
    Q_K = fit$Q_K,
    rw_q = rw_q,
    feature_active = feature_active,
    include_prior = FALSE
  )

  prior_entropy <- .cavi_prior_terms_from_state(
    q_u = q_u,
    lambda_vec = fit$lambda_vec,
    Q_K = fit$Q_K,
    rw_q = rw_q,
    lambda_sd_prior_rate = lambda_sd_prior_rate
  ) + .cavi_entropy_u_terms(q_u, K = K)
  prior_entropy[!feature_active] <- 0

  list(
    like_score = like_score,
    prior_entropy = prior_entropy,
    cell_terms = .cavi_cell_terms_from_state(fit$gamma, fit$params$pi)
  )
}

.cavi_partition_objective_from_fits <- function(fits, X, weights,
                                                T_now = 1,
                                                active_orderings = NULL,
                                                active_feature_pairs = NULL,
                                                assignment_M = NULL,
                                                partition_prior = c("adaptive", "fixed"),
                                                partition_prior_init = NULL,
                                                assignment_prior = NULL,
                                                ordering_alpha = NULL,
                                                drop_unused_ordering = FALSE,
                                                partition_prior_missing = FALSE,
                                                caller = ".cavi_partition_objective_from_fits()") {
  X <- as.matrix(X)
  weights <- as.matrix(weights)
  M <- length(fits)
  d <- ncol(X)
  assignment_ctl <- .validate_partition_assignment_controls(
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    M = M,
    partition_prior_missing = partition_prior_missing,
    caller = caller
  )
  if (is.null(active_orderings)) {
    active_orderings <- rep(TRUE, M)
  }
  active_orderings <- as.logical(active_orderings)
  if (length(active_orderings) != M) {
    stop("active_orderings must have length equal to length(fits).")
  }
  if (is.null(active_feature_pairs)) {
    active_feature_pairs <- matrix(rep(active_orderings, each = d), nrow = d, ncol = M)
  } else {
    active_feature_pairs <- as.matrix(active_feature_pairs)
    if (!all(dim(active_feature_pairs) == c(d, M))) {
      stop("active_feature_pairs must be a d x M matrix.")
    }
  }
  active_feature_pairs[, !active_orderings] <- FALSE
  if (!any(active_orderings)) {
    stop("At least one ordering must remain active.")
  }
  if (is.null(assignment_M)) {
    assignment_M <- M
  }
  if (!is.numeric(assignment_M) || length(assignment_M) != 1L ||
      !is.finite(assignment_M) || assignment_M < 1) {
    stop("assignment_M must be a single finite number >= 1.")
  }

  term_list <- lapply(seq_along(fits), function(m) {
    fit_m <- fits[[m]]
    fit_m$control <- utils::modifyList(
      fit_m$control %||% list(),
      list(feature_active = as.logical(active_feature_pairs[, m]))
    )
    .cavi_partition_terms_from_fit(fit_m, data = X)
  })
  like_mat <- vapply(term_list, `[[`, numeric(d), "like_score")
  prior_entropy_raw <- vapply(term_list, `[[`, numeric(d), "prior_entropy")
  prior_entropy_mat <- prior_entropy_raw
  prior_entropy_mat[, !active_orderings] <- 0
  cell_terms_by_fit <- vapply(term_list, `[[`, numeric(1), "cell_terms")
  cell_terms <- if (isTRUE(drop_unused_ordering)) {
    sum(cell_terms_by_fit[active_orderings])
  } else {
    sum(cell_terms_by_fit)
  }
  weights_assignment <- .partition_assignment_weights(
    weights = weights,
    active_orderings = active_orderings,
    active_feature_pairs = active_feature_pairs
  )
  assignment_info <- .cavi_partition_assignment_info(
    weights = weights,
    T_now = T_now,
    partition_prior = assignment_ctl$partition_prior,
    partition_prior_init = assignment_ctl$partition_prior_init,
    assignment_prior = assignment_ctl$legacy_assignment_prior,
    ordering_alpha = assignment_ctl$ordering_alpha,
    assignment_M = assignment_M,
    active_orderings = active_orderings,
    active_feature_pairs = active_feature_pairs,
    drop_unused_ordering = drop_unused_ordering
  )

  list(
    objective = as.numeric(
      sum(weights * like_mat) +
        sum(prior_entropy_mat[, active_orderings, drop = FALSE]) +
        cell_terms +
        assignment_info$objective
    ),
    like_mat = like_mat,
    weights_assignment = weights_assignment,
    prior_entropy_mat = prior_entropy_mat,
    prior_entropy_raw = prior_entropy_raw,
    cell_terms = as.numeric(cell_terms),
    cell_terms_by_fit = as.numeric(cell_terms_by_fit),
    assignment_log_prior = as.numeric(assignment_info$expected_log_assign),
    assignment_info = assignment_info,
    z_entropy = as.numeric(assignment_info$z_entropy),
    active_orderings = active_orderings,
    assignment_M = as.numeric(assignment_info$assignment_M)
  )
}

.validate_partition_inactive_controls <- function(freeze_unused_ordering = TRUE,
                                                  freeze_unused_ordering_threshold = 0.5,
                                                  drop_unused_ordering = FALSE) {
  if (!is.logical(freeze_unused_ordering) || length(freeze_unused_ordering) != 1L ||
      is.na(freeze_unused_ordering)) {
    stop("freeze_unused_ordering must be a single TRUE/FALSE value.")
  }
  if (!is.logical(drop_unused_ordering) || length(drop_unused_ordering) != 1L ||
      is.na(drop_unused_ordering)) {
    stop("drop_unused_ordering must be a single TRUE/FALSE value.")
  }
  if (!is.numeric(freeze_unused_ordering_threshold) ||
      length(freeze_unused_ordering_threshold) != 1L ||
      !is.finite(freeze_unused_ordering_threshold) ||
      freeze_unused_ordering_threshold < 0) {
    stop("freeze_unused_ordering_threshold must be a single finite non-negative number.")
  }
  list(
    freeze_unused_ordering = isTRUE(freeze_unused_ordering),
    freeze_unused_ordering_threshold = as.numeric(freeze_unused_ordering_threshold),
    drop_unused_ordering = isTRUE(drop_unused_ordering)
  )
}

.soft_partition_softmax <- function(score_mat, T_now = 1,
                                    active_orderings = NULL,
                                    active_feature_pairs = NULL) {
  score_mat <- as.matrix(score_mat)
  d <- nrow(score_mat)
  M <- ncol(score_mat)
  if (is.null(active_orderings)) {
    active_orderings <- rep(TRUE, M)
  }
  active_orderings <- as.logical(active_orderings)
  if (length(active_orderings) != M) {
    stop("active_orderings must have length equal to ncol(score_mat).")
  }
  if (!any(active_orderings)) {
    stop("At least one ordering must remain active.")
  }
  if (is.null(active_feature_pairs)) {
    active_feature_pairs <- matrix(rep(active_orderings, each = d), nrow = d, ncol = M)
  } else {
    active_feature_pairs <- as.matrix(active_feature_pairs)
    if (!all(dim(active_feature_pairs) == c(d, M))) {
      stop("active_feature_pairs must be a d x M matrix.")
    }
  }
  active_feature_pairs[, !active_orderings] <- FALSE

  weights <- matrix(0, nrow = d, ncol = M)
  for (j in seq_len(d)) {
    active_idx <- which(active_feature_pairs[j, ])
    if (!length(active_idx)) {
      active_idx <- which(active_orderings)
      active_idx <- active_idx[which.max(score_mat[j, active_idx])]
    }
    if (length(active_idx) == 1L) {
      weights[j, active_idx] <- 1
      next
    }
    log_w <- score_mat[j, active_idx] / T_now
    log_w <- log_w - max(log_w)
    w <- exp(log_w)
    weights[j, active_idx] <- w / sum(w)
  }
  colnames(weights) <- colnames(score_mat)
  weights
}

.partition_preserved_frozen_weights <- function(reported_weights, active_feature_pairs) {
  reported_weights <- as.matrix(reported_weights)
  active_feature_pairs <- as.matrix(active_feature_pairs)
  frozen <- reported_weights
  frozen[active_feature_pairs] <- 0
  colnames(frozen) <- colnames(reported_weights)
  frozen
}

.partition_compose_reported_weights <- function(effective_weights,
                                                active_feature_pairs,
                                                reported_prev,
                                                drop_unused_ordering = FALSE) {
  effective_weights <- as.matrix(effective_weights)
  if (isTRUE(drop_unused_ordering)) {
    out <- effective_weights
    colnames(out) <- colnames(effective_weights)
    return(out)
  }
  active_feature_pairs <- as.matrix(active_feature_pairs)
  preserved <- .partition_preserved_frozen_weights(reported_prev, active_feature_pairs)
  remaining_mass <- pmax(1 - rowSums(preserved), 0)
  out <- preserved
  if (any(active_feature_pairs)) {
    row_idx <- row(effective_weights)[active_feature_pairs]
    out[active_feature_pairs] <- effective_weights[active_feature_pairs] * remaining_mass[row_idx]
  }
  colnames(out) <- colnames(effective_weights)
  out
}

.partition_register_new_frozen_weights <- function(reported_prev,
                                                   snapshot_weights,
                                                   drop_idx = integer(0),
                                                   frozen_mask = NULL) {
  reported_prev <- as.matrix(reported_prev)
  snapshot_weights <- as.matrix(snapshot_weights)
  if (length(drop_idx) > 0L) {
    reported_prev[, drop_idx] <- snapshot_weights[, drop_idx, drop = FALSE]
  }
  if (!is.null(frozen_mask)) {
    frozen_mask <- as.matrix(frozen_mask)
    reported_prev[frozen_mask] <- snapshot_weights[frozen_mask]
  }
  reported_prev
}

.make_partition_ordering_event <- function(event, ordering, ordering_label,
                                           step = NA_integer_, T_now = NA_real_,
                                           usage_mass = NA_real_,
                                           threshold = NA_real_,
                                           stage = NA_character_) {
  list(
    event = as.character(event),
    ordering = as.integer(ordering),
    label = as.character(ordering_label),
    step = as.integer(step),
    T_now = as.numeric(T_now),
    usage_mass = as.numeric(usage_mass),
    threshold = as.numeric(threshold),
    stage = as.character(stage)
  )
}

.make_partition_feature_event <- function(event, feature, ordering, ordering_label,
                                          step = NA_integer_, T_now = NA_real_,
                                          weight = NA_real_,
                                          threshold = NA_real_,
                                          stage = NA_character_) {
  list(
    event = as.character(event),
    feature = as.integer(feature),
    ordering = as.integer(ordering),
    label = as.character(ordering_label),
    step = as.integer(step),
    T_now = as.numeric(T_now),
    weight = as.numeric(weight),
    threshold = as.numeric(threshold),
    stage = as.character(stage)
  )
}

.partition_freeze_orderings <- function(weights,
                                        active_orderings,
                                        threshold,
                                        iter_index = NA_integer_,
                                        T_now = NA_real_,
                                        ordering_labels = NULL,
                                        stage = c("pre_update", "post_update")) {
  stage <- match.arg(stage)
  weights <- as.matrix(weights)
  M <- ncol(weights)

  if (is.null(ordering_labels)) {
    ordering_labels <- .cavi_partition_order_labels(M)
  }
  if (length(active_orderings) != M) {
    stop("active_orderings must have length equal to ncol(weights).")
  }

  active_now <- as.logical(active_orderings)
  usage_mass <- colSums(weights)
  ordering_events <- list()
  drop_idx <- integer(0)

  if (sum(active_now) > 1L) {
    drop_idx <- which(active_now & usage_mass <= threshold)
    if (length(drop_idx) >= sum(active_now)) {
      keep_idx <- which.max(ifelse(active_now, usage_mass, -Inf))
      drop_idx <- setdiff(drop_idx, keep_idx)
    }
    if (length(drop_idx) > 0L) {
      stage_text <- if (identical(stage, "pre_update")) {
        "before weighted update"
      } else {
        "after weighted update"
      }
      for (m in drop_idx) {
        active_now[m] <- FALSE
        ordering_events[[length(ordering_events) + 1L]] <- .make_partition_ordering_event(
          event = "freeze",
          ordering = m,
          ordering_label = ordering_labels[m],
          step = iter_index,
          T_now = T_now,
          usage_mass = usage_mass[m],
          threshold = threshold,
          stage = stage
        )
        warning(
          sprintf(
            "soft_partition_cavi froze ordering %s %s at step %d: posterior feature mass %.6f was at or below threshold %.6f.",
            ordering_labels[m],
            stage_text,
            as.integer(iter_index),
            usage_mass[m],
            threshold
          ),
          call. = FALSE
        )
      }
    }
  }

  list(
    active_orderings = active_now,
    drop_idx = drop_idx,
    usage_mass = usage_mass,
    ordering_events = ordering_events
  )
}

.partition_freeze_features <- function(weights,
                                       active_feature_pairs,
                                       active_orderings,
                                       threshold,
                                       iter_index = NA_integer_,
                                       T_now = NA_real_,
                                       ordering_labels = NULL,
                                       stage = c("pre_update", "post_update")) {
  stage <- match.arg(stage)
  weights <- as.matrix(weights)
  active_feature_pairs <- as.matrix(active_feature_pairs)
  d <- nrow(weights)
  M <- ncol(weights)
  if (!all(dim(active_feature_pairs) == c(d, M))) {
    stop("active_feature_pairs must be a d x M matrix.")
  }
  active_now <- active_feature_pairs
  active_now[, !active_orderings] <- FALSE
  if (is.null(ordering_labels)) {
    ordering_labels <- .cavi_partition_order_labels(M)
  }

  feature_events <- list()
  new_frozen_mask <- matrix(FALSE, nrow = d, ncol = M)
  for (j in seq_len(d)) {
    active_idx <- which(active_now[j, ] & active_orderings)
    if (length(active_idx) <= 1L) next
    drop_idx <- active_idx[weights[j, active_idx] <= threshold]
    if (length(drop_idx) >= length(active_idx)) {
      keep_idx <- active_idx[which.max(weights[j, active_idx])]
      drop_idx <- setdiff(drop_idx, keep_idx)
    }
    if (!length(drop_idx)) next
    active_now[j, drop_idx] <- FALSE
    new_frozen_mask[j, drop_idx] <- TRUE
    for (m in drop_idx) {
      feature_events[[length(feature_events) + 1L]] <- .make_partition_feature_event(
        event = "freeze_feature",
        feature = j,
        ordering = m,
        ordering_label = ordering_labels[m],
        step = iter_index,
        T_now = T_now,
        weight = weights[j, m],
        threshold = threshold,
        stage = stage
      )
    }
  }

  list(
    active_feature_pairs = active_now,
    new_frozen_mask = new_frozen_mask,
    feature_events = feature_events
  )
}

.cavi_feature_scores_from_state <- function(X, R, q_u, sigma2, lambda_vec,
                                            measurement_sd = NULL,
                                            Q_K, rw_q = 2L,
                                            lambda_sd_prior_rate = NULL,
                                            include_prior = TRUE,
                                            feature_active = NULL) {
  X <- as.matrix(X)
  R <- as.matrix(R)
  n <- nrow(X)
  d <- ncol(X)
  K <- ncol(R)
  if (is.null(feature_active)) {
    feature_active <- rep(TRUE, d)
  }
  feature_active <- as.logical(feature_active)
  if (length(feature_active) != d) {
    stop("feature_active must have length d.")
  }
  prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
  r_rank <- prior_meta$rank
  logdet_Q <- prior_meta$logdet
  noise_info <- .cavi_resolve_measurement_sd(
    S = measurement_sd,
    X = X,
    caller = ".cavi_feature_scores_from_state()"
  )

  scores <- numeric(d)
  m_mat <- q_u$m_mat
  sdiag_mat <- q_u$sdiag_mat

  for (j in seq_len(d)) {
    if (!feature_active[j]) {
      scores[j] <- 0
      next
    }
    diff_j <- matrix(X[, j], nrow = n, ncol = K) -
      matrix(m_mat[j, ], nrow = n, ncol = K, byrow = TRUE)
    quad_jk <- colSums(R * (diff_j^2))
    unc_jk <- colSums(R) * sdiag_mat[j, ]
    if (is.null(noise_info$measurement_sd)) {
      data_j <- -0.5 * sum(
        colSums(R) * log(2 * base::pi * sigma2[j]) +
          (quad_jk + unc_jk) / sigma2[j]
      )
    } else {
      var_j <- noise_info$var_mat[, j]
      const_j <- sum(log(2 * base::pi * var_j))
      unc_row_j <- rowSums(R * rep(sdiag_mat[j, ], each = n))
      data_j <- -0.5 * (
        const_j +
          sum(rowSums(R * (diff_j^2) / var_j)) +
          sum(unc_row_j / var_j)
      )
    }

    prior_j <- 0
    if (isTRUE(include_prior)) {
      prior_j <- 0.5 * (r_rank * (log(lambda_vec[j]) - log(2 * base::pi)) + logdet_Q) -
        0.5 * lambda_vec[j] * q_u$eq_quad[j] +
        .lambda_sd_prior_terms(
          lambda_vec = lambda_vec[j],
          rate = lambda_sd_prior_rate,
          include_constant = TRUE
        )
    }
    scores[j] <- data_j + prior_j
  }

  scores
}

.cavi_global_terms_from_fit <- function(fit) {
  if (!inherits(fit, "cavi")) stop("fit must inherit from class 'cavi'.")
  q_u <- .cavi_posterior_state_from_fit(fit)
  K <- length(fit$params$pi)
  .cavi_cell_terms_from_state(fit$gamma, fit$params$pi) +
    .cavi_entropy_u(q_u = q_u, K = K)
}

.cavi_fixed_R_local_feature_fit <- function(xj,
                                            R,
                                            Q_K,
                                            measurement_sd = NULL,
                                            rw_q = 2L,
                                            lambda_init = 1,
                                            lambda_sd_prior_rate = NULL,
                                            sigma2_init = NULL,
                                            lambda_min = 1e-10,
                                            lambda_max = 1e10,
                                            sigma_min = 1e-10,
                                            sigma_max = 1e10,
                                            max_iter = 25L,
                                            tol = 1e-6) {
  X1 <- matrix(as.numeric(xj), ncol = 1L)
  K <- ncol(R)
  prior_meta <- .rw_precision_metadata(Q_K, rw_q = rw_q)
  r_rank <- prior_meta$rank
  logdet_Q <- prior_meta$logdet

  sigma2 <- if (!is.null(measurement_sd)) {
    NULL
  } else if (is.null(sigma2_init)) {
    pmin(pmax(stats::var(as.numeric(xj)), sigma_min), sigma_max)
  } else {
    pmin(pmax(as.numeric(sigma2_init)[1], sigma_min), sigma_max)
  }
  lambda_j <- pmin(pmax(as.numeric(lambda_init)[1], lambda_min), lambda_max)

  compute_local_elbo <- function(q_u, sigma2, lambda_j) {
    scores <- .cavi_feature_scores_from_state(
      X = X1,
      R = R,
      q_u = q_u,
      sigma2 = sigma2,
      measurement_sd = measurement_sd,
      lambda_vec = lambda_j,
      Q_K = Q_K,
      rw_q = rw_q,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      include_prior = TRUE
    )
    data_prior <- scores[1]
    entropy_u <- .cavi_entropy_u(q_u, K = K)
    as.numeric(data_prior + entropy_u)
  }

  q_u <- .cavi_update_q_u_weighted(
    X = X1,
    R = R,
    sigma2 = sigma2,
    lambda_vec = lambda_j,
    Q_K = Q_K,
    measurement_sd = measurement_sd,
    feature_weights = 1,
    rw_q = rw_q
  )
  elbo_prev <- compute_local_elbo(q_u, sigma2, lambda_j)

  for (iter in seq_len(as.integer(max_iter))) {
    sigma2 <- .cavi_update_sigma2_weighted(
      X = X1,
      R = R,
      q_u = q_u,
      sigma2_current = sigma2,
      measurement_sd = measurement_sd,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      feature_weights = 1
    )
    lambda_j <- .cavi_update_lambda_weighted(
      q_u = q_u,
      r_rank = r_rank,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      lambda_current = lambda_j,
      feature_weights = 1
    )
    q_u <- .cavi_update_q_u_weighted(
      X = X1,
      R = R,
      sigma2 = sigma2,
      lambda_vec = lambda_j,
      Q_K = Q_K,
      measurement_sd = measurement_sd,
      feature_weights = 1,
      rw_q = rw_q
    )
    elbo_new <- compute_local_elbo(q_u, sigma2, lambda_j)
    if (abs(elbo_new - elbo_prev) / (abs(elbo_prev) + 1) < tol) break
    elbo_prev <- elbo_new
  }

  score <- .cavi_feature_scores_from_state(
    X = X1,
    R = R,
    q_u = q_u,
    sigma2 = sigma2,
    measurement_sd = measurement_sd,
    lambda_vec = lambda_j,
    Q_K = Q_K,
    rw_q = rw_q,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    include_prior = TRUE
  )[1]

  list(
    score = as.numeric(score),
    sigma2 = if (is.null(sigma2)) NA_real_ else as.numeric(sigma2),
    lambda = as.numeric(lambda_j),
    q_u = q_u,
    elbo = as.numeric(elbo_prev),
    logdet_Q = as.numeric(logdet_Q),
    r_rank = as.integer(r_rank)
  )
}

.do_cavi_weighted <- function(object,
                              data = NULL,
                              feature_weights = NULL,
                              feature_active = NULL,
                              iter = 1L,
                              tol = 0,
                              lambda_sd_prior_rate = NULL,
                              sigma_min = 1e-10,
                              sigma_max = 1e10,
                              lambda_min = 1e-10,
                              lambda_max = 1e10,
                              verbose = FALSE) {
  if (!inherits(object, "cavi")) stop("object must inherit from class 'cavi'.")

  X <- as.matrix(data %||% object$data)
  if (is.null(X)) stop("data must be supplied for weighted CAVI updates.")
  iter <- as.integer(iter)
  if (iter < 1L) stop("iter must be >= 1.")

  d <- ncol(X)
  K <- length(object$params$pi)
  w <- .cavi_resolve_feature_weights(feature_weights, d)
  if (is.null(feature_active)) {
    feature_active <- rep(TRUE, d)
  }
  feature_active <- as.logical(feature_active)
  if (length(feature_active) != d) {
    stop("feature_active must have length d.")
  }
  w_eff <- w
  w_eff[!feature_active] <- 0
  R <- as.matrix(object$gamma)
  pi_vec <- as.numeric(object$params$pi)
  position_prior <- (object$control %||% list())$position_prior %||% "adaptive"
  measurement_sd <- object$measurement_sd %||% NULL
  sigma2 <- if (is.null(measurement_sd)) as.numeric(object$params$sigma2) else NULL
  lambda_vec <- as.numeric(object$lambda_vec)
  Q_K <- object$Q_K
  rw_q <- object$rw_q %||% 2L
  lambda_sd_prior_rate_use <- if (is.null(lambda_sd_prior_rate)) {
    (object$control %||% list())$lambda_sd_prior_rate %||% NULL
  } else {
    .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
  }
  meta <- .cavi_prepare_weighted_metadata(K = K, rw_q = rw_q, Q_K = Q_K)
  q_u_current <- .cavi_posterior_state_from_fit(object)

  if (!any(feature_active)) {
    out <- object
    out$control$feature_weights <- w_eff
    out$control$feature_active <- feature_active
    out$control$adaptive <- "variational"
    out$control$position_prior <- position_prior
    out$control$lambda_sd_prior_rate <- lambda_sd_prior_rate_use
    out$data <- X
    return(out)
  }

  elbo_trace <- numeric(iter + 1L)
  q_u <- .cavi_update_q_u_weighted(
    X = X,
    R = R,
    sigma2 = sigma2,
    lambda_vec = lambda_vec,
    Q_K = Q_K,
    measurement_sd = measurement_sd,
    feature_weights = w_eff,
    rw_q = rw_q,
    feature_active = feature_active,
    q_u_current = q_u_current
  )

  compute_fit_elbo <- function(R, pi_vec, sigma2, lambda_vec, q_u) {
    like_scores <- .cavi_feature_scores_from_state(
      X = X,
      R = R,
      q_u = q_u,
      sigma2 = sigma2,
      measurement_sd = measurement_sd,
      lambda_vec = lambda_vec,
      Q_K = Q_K,
      rw_q = rw_q,
      lambda_sd_prior_rate = lambda_sd_prior_rate_use,
      include_prior = FALSE,
      feature_active = feature_active
    )
    prior_terms <- .cavi_prior_terms_from_state(
      q_u = q_u,
      lambda_vec = lambda_vec,
      Q_K = Q_K,
      rw_q = rw_q,
      lambda_sd_prior_rate = lambda_sd_prior_rate_use
    )
    entropy_terms <- .cavi_entropy_u_terms(q_u, K = K)
    prior_terms[!feature_active] <- 0
    entropy_terms[!feature_active] <- 0
    cell_terms <- .cavi_cell_terms_from_state(R, pi_vec)

    as.numeric(
      cell_terms +
        sum(w_eff * like_scores + prior_terms + entropy_terms)
    )
  }

  elbo_trace[1] <- compute_fit_elbo(R, pi_vec, sigma2, lambda_vec, q_u)
  converged <- FALSE

  for (step in seq_len(iter)) {
    R <- .cavi_update_r_weighted(
      X = X,
      q_u = q_u,
      pi_vec = pi_vec,
      sigma2 = sigma2,
      measurement_sd = measurement_sd,
      feature_weights = w_eff
    )
    if (identical(position_prior, "adaptive")) {
      pi_vec <- pmax(colMeans(R), .Machine$double.eps)
      pi_vec <- pi_vec / sum(pi_vec)
    }
    sigma2 <- .cavi_update_sigma2_weighted(
      X = X,
      R = R,
      q_u = q_u,
      sigma2_current = sigma2,
      measurement_sd = measurement_sd,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      feature_weights = w_eff,
      feature_active = feature_active
    )
    lambda_vec <- .cavi_update_lambda_weighted(
      q_u = q_u,
      r_rank = meta$r_rank,
      lambda_sd_prior_rate = lambda_sd_prior_rate_use,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      lambda_current = lambda_vec,
      feature_weights = w_eff,
      feature_active = feature_active
    )
    q_u <- .cavi_update_q_u_weighted(
      X = X,
      R = R,
      sigma2 = sigma2,
      lambda_vec = lambda_vec,
      Q_K = Q_K,
      measurement_sd = measurement_sd,
      feature_weights = w_eff,
      rw_q = rw_q,
      feature_active = feature_active,
      q_u_current = q_u
    )
    elbo_trace[step + 1L] <- compute_fit_elbo(R, pi_vec, sigma2, lambda_vec, q_u)

    delta <- elbo_trace[step + 1L] - elbo_trace[step]
    if (verbose) {
      cat(sprintf("[cavi-weighted %3d] ELBO=%.6f delta=%.3e\n",
                  step, elbo_trace[step + 1L], delta))
    }
    if (delta < -1e-8) {
      warning(sprintf("Weighted CAVI ELBO decreased by %.3e at step %d.", delta, step),
              call. = FALSE)
    }
    if (delta >= 0 && abs(delta) / (abs(elbo_trace[step]) + 1) < tol) {
      converged <- TRUE
      elbo_trace <- elbo_trace[seq_len(step + 1L)]
      break
    }
  }

  params <- list(
    pi = pi_vec,
    mu = lapply(seq_len(K), function(k) q_u$m_mat[, k]),
    sigma2 = sigma2
  )

  out <- object
  out$params <- params
  out$gamma <- R
  out$measurement_sd <- measurement_sd
  out$posterior <- list(
    mean = q_u$m_mat,
    cov = q_u$S_list,
    var = q_u$sdiag_mat
  )
  out$lambda_vec <- lambda_vec
  out$elbo_trace <- c(out$elbo_trace %||% numeric(0), elbo_trace[-1L])
  out$iter <- length(out$elbo_trace)
  out$converged <- converged
  out$control$feature_weights <- w_eff
  out$control$feature_active <- feature_active
  out$control$adaptive <- "variational"
  out$control$position_prior <- position_prior
  out$control$lambda_sd_prior_rate <- lambda_sd_prior_rate_use
  out$data <- X
  out
}


# ============================================================
# Feature scoring and soft partition based on CAVI
# ============================================================

#' Score features under a fitted \code{cavi} model
#'
#' @param fit A \code{cavi} object.
#' @param X Optional data matrix. If NULL, uses \code{fit$data}.
#' @param include_prior Logical; include the GMRF prior contribution.
#'
#' @return A list with \code{feature_score}, \code{global_terms}, and
#'   \code{entropy_u}.
#' @noRd
score_features_onefit_cavi <- function(fit, X = NULL, include_prior = TRUE) {
  if (!inherits(fit, "cavi")) stop("fit must be a 'cavi' object.")
  X <- as.matrix(X %||% fit$data)
  if (is.null(X)) stop("X must be provided or stored in fit$data.")
  q_u <- .cavi_posterior_state_from_fit(fit)

  scores <- .cavi_feature_scores_from_state(
    X = X,
    R = fit$gamma,
    q_u = q_u,
    sigma2 = fit$params$sigma2,
    measurement_sd = fit$measurement_sd %||% NULL,
    lambda_vec = fit$lambda_vec,
    Q_K = fit$Q_K,
    rw_q = fit$rw_q %||% 2L,
    lambda_sd_prior_rate = (fit$control %||% list())$lambda_sd_prior_rate %||% NULL,
    include_prior = include_prior
  )

  list(
    feature_score = scores,
    global_terms = .cavi_global_terms_from_fit(fit),
    entropy_u = .cavi_entropy_u(q_u = q_u, K = length(fit$params$pi))
  )
}

#' Partition features by comparing CAVI feature scores from two fits
#'
#' @param fitA,fitB \code{cavi} objects.
#' @param X Optional data matrix.
#' @param delta Nonnegative assignment margin.
#' @param include_prior Logical.
#'
#' @return A list with \code{assign}, \code{score_diff}, \code{CA}, and \code{CB}.
#' @noRd
partition_features_twofits_cavi <- function(fitA, fitB, X = NULL, delta = 0,
                                            include_prior = TRUE) {
  if (!inherits(fitA, "cavi") || !inherits(fitB, "cavi")) {
    stop("fitA and fitB must both be 'cavi' objects.")
  }
  X <- as.matrix(X %||% fitA$data)
  if (is.null(X)) stop("X must be provided or stored in fitA$data.")

  CA <- score_features_onefit_cavi(fitA, X = X, include_prior = include_prior)$feature_score
  CB <- score_features_onefit_cavi(fitB, X = X, include_prior = include_prior)$feature_score
  score_diff <- CA - CB
  assign <- ifelse(score_diff > delta, "A", "B")

  list(assign = assign, score_diff = score_diff, CA = CA, CB = CB)
}

.cavi_fit_from_method <- function(X,
                                  S = NULL,
                                  method,
                                  K,
                                  rw_q = 2L,
                                  ridge = 0,
                                  lambda_sd_prior_rate = NULL,
                                  lambda_min = 1e-10,
                                  lambda_max = 1e10,
                                  sigma_min = 1e-10,
                                  sigma_max = 1e10,
                                  max_iter = 5L,
                                  tol = 1e-6,
                                  discretization = c("quantile", "equal", "kmeans"),
                                  pca_component = 1L,
                                  strict_K = FALSE,
                                  verbose = FALSE) {
  method <- match.arg(method, c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"))
  discretization <- match.arg(discretization)

  if (method == "random") {
    return(cavi(
      X = X,
      K = K,
      method = "random",
      S = S,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose
    ))
  }

  ordering_result <- .cavi_get_ordering_result(X, method = method, pca_component = pca_component)
  .cavi_build_from_ordering(
    X = X,
    ordering_vec = ordering_result$t,
    S = S,
    K = K,
    rw_q = rw_q,
    ridge = ridge,
    lambda_init = 1,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    max_iter = max_iter,
    tol = tol,
    discretization = discretization,
    strict_K = strict_K,
    ordering_label = if (method == "PCA") {
      sprintf("%s(PC%d)", method, pca_component)
    } else {
      method
    },
    verbose = verbose
  )
}

#' Initialize two CAVI trajectories for dual-ordering partitioning
#'
#' @param X Numeric matrix (\code{n x d}).
#' @param method Initialisation strategy for the second trajectory.
#' @param method1 Initialisation method for the first trajectory.
#' @param K Number of pseudotime bins.
#' @param rw_q Random-walk order.
#' @param ridge Optional nugget added to the RW precision. Mainly useful for
#'   internal diagnostics that compare intrinsic and properized priors.
#' @param lambda_sd_prior_rate Optional positive rate for an exponential prior
#'   on \code{1 / sqrt(lambda_j^{(m)})}. The default \code{NULL} means no lambda
#'   prior penalty. For backward compatibility, an explicit \code{0} is treated
#'   the same way; it is only an alias for "no penalty" and does not
#'   correspond to a literal exponential prior with rate zero.
#' @param lambda_min,lambda_max Bounds used when clipping the initial
#'   feature-specific \code{lambda_j} values for each warm-start fit.
#' @param sigma_min,sigma_max Bounds for the feature-specific noise variance
#'   \code{sigma2_j}. Forwarded to internal \code{cavi()} calls.
#' @param discretization Initial discretization scheme used when converting
#'   ordering scores into \code{K} bins. For dual-ordering initialisation the
#'   requested \code{K} is enforced strictly; if quantile cuts collapse, the
#'   code falls back to equal-width bins so both orderings remain comparable.
#' @param num_iter Warm-start CAVI sweeps per fit.
#' @param seed Random seed used by \code{"random_split"}.
#' @param verbose Logical.
#'
#' @return A list with \code{fit1}, \code{fit2}, and \code{seed2}.
#' @noRd
init_two_trajectories_cavi <- function(X,
                                       S = NULL,
                                       method = c("score", "mincor", "random_split"),
                                       method1 = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
                                       K = NULL,
                                       rw_q = 2L,
                                       ridge = 0,
                                       lambda_sd_prior_rate = NULL,
                                       lambda_min = 1e-10,
                                       lambda_max = 1e10,
                                       sigma_min = 1e-10,
                                       sigma_max = 1e10,
                                       discretization = c("quantile", "equal", "kmeans"),
                                       num_iter = 5L,
                                       seed = 42L,
                                       verbose = FALSE) {
  method <- match.arg(method)
  method1 <- match.arg(method1)
  discretization <- match.arg(discretization)
  X <- as.matrix(X)
  d <- ncol(X)
  K_use <- .cavi_resolve_K(X, K)

  fit1 <- .cavi_fit_from_method(
    X = X,
    S = S,
    K = K_use,
    method = method1,
    rw_q = rw_q,
    ridge = ridge,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    max_iter = num_iter,
    discretization = discretization,
    strict_K = TRUE,
    verbose = FALSE
  )

  seed2 <- NA_integer_
  if (method == "score") {
    scores1 <- score_features_onefit_cavi(fit1, X = X, include_prior = TRUE)$feature_score
    seed2 <- which.min(scores1)
    fit2 <- .fit_cavi_from_seed_col(
      X = X,
      S = S,
      seed_col = seed2,
      K = K_use,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = num_iter,
      discretization = discretization,
      strict_K = TRUE,
      verbose = FALSE
    )
  } else if (method == "mincor") {
    pc1_loadings <- abs(prcomp(X, center = TRUE, scale. = FALSE)$rotation[, 1])
    seed1_col <- which.max(pc1_loadings)
    abs_cors <- abs(cor(X[, seed1_col], X))
    abs_cors[seed1_col] <- Inf
    seed2 <- which.min(abs_cors)
    fit2 <- .fit_cavi_from_seed_col(
      X = X,
      S = S,
      seed_col = seed2,
      K = K_use,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = num_iter,
      discretization = discretization,
      strict_K = TRUE,
      verbose = FALSE
    )
  } else {
    set.seed(seed)
    perm <- sample(d)
    half1 <- sort(perm[seq_len(floor(d / 2))])
    half2 <- sort(perm[(floor(d / 2) + 1):d])
    fit1_half <- .cavi_fit_from_method(
      X = X[, half1, drop = FALSE],
      S = .cavi_subset_measurement_sd(S, half1),
      method = "PCA",
      pca_component = 1L,
      K = K_use,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = num_iter,
      discretization = discretization,
      strict_K = TRUE,
      verbose = FALSE
    )
    fit2_half <- .cavi_fit_from_method(
      X = X[, half2, drop = FALSE],
      S = .cavi_subset_measurement_sd(S, half2),
      method = "PCA",
      pca_component = 1L,
      K = K_use,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = num_iter,
      discretization = discretization,
      strict_K = TRUE,
      verbose = FALSE
    )
    fit1_args <- list(
      X = X,
      K = K_use,
      responsibilities_init = fit1_half$gamma,
      position_prior_init = fit1_half$params$pi,
      S = S,
      lambda_init = rep(1, d),
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      rw_q = rw_q,
      ridge = ridge,
      discretization = discretization,
      max_iter = num_iter,
      verbose = FALSE
    )
    fit2_args <- fit1_args
    fit2_args$responsibilities_init <- fit2_half$gamma
    fit2_args$position_prior_init <- fit2_half$params$pi
    if (is.null(S)) {
      fit1_args$sigma2_init <- pmax(apply(X, 2, stats::var), 1e-10)
      fit2_args$sigma2_init <- fit1_args$sigma2_init
    }
    fit1 <- do.call(cavi, fit1_args)
    fit2 <- do.call(cavi, fit2_args)
  }

  if (isTRUE(verbose)) {
    cat(sprintf("[init_cavi] fit1=%s fit2_strategy=%s seed2=%s\n",
                method1, method, ifelse(is.na(seed2), "NA", as.character(seed2))))
  }

  list(fit1 = fit1, fit2 = fit2, seed2 = seed2)
}


# ============================================================
# Generalized M-ordering soft partition
# ============================================================

#' Initialize M CAVI trajectories for multi-ordering partitioning
#'
#' @param X Numeric matrix (\code{n x d}).
#' @param M Integer >= 2. Number of orderings.
#' @param methods Character vector of length M. Each element is an ordering
#'   method (e.g., \code{"PCA"}, \code{"fiedler"}, \code{"random"}).
#' @param pca_components Integer vector of length M. For PCA-based methods,
#'   which PC component to use. Ignored for non-PCA methods.
#' @param partition_init Either \code{"similarity"} for feature-similarity
#'   clustering followed by block-specific ordering initialization or
#'   \code{"ordering_methods"} for the existing per-ordering warm starts. The
#'   default is \code{"similarity"}.
#' @param similarity_metric For \code{partition_init = "similarity"} only:
#'   feature-similarity metric used to construct the clustering. One of
#'   \code{"spearman"}, \code{"pearson"}, or \code{"smooth_fit"}.
#'   \code{"smooth_fit"} can be used with \code{ridge = 0}, but for
#'   \code{M > 1} this uses an intrinsic-RW pseudo-evidence rather than a
#'   fully proper marginal likelihood; use a small positive \code{ridge} if
#'   you want the smoother evidence to be theoretically proper.
#' @param smooth_fit_lambda_mode For \code{similarity_metric = "smooth_fit"}
#'   only: whether the directional smoother optimizes \code{lambda} or keeps it
#'   fixed at \code{smooth_fit_lambda_value}.
#' @param smooth_fit_lambda_value For \code{similarity_metric = "smooth_fit"}
#'   only: fixed \code{lambda} value used when
#'   \code{smooth_fit_lambda_mode = "fixed"}, and the starting value when
#'   \code{smooth_fit_lambda_mode = "optimize"}.
#' @param cluster_linkage For \code{partition_init = "similarity"} only:
#'   hierarchical-clustering linkage applied to \code{1 - S(X)}.
#' @param similarity_min_feature_sd For \code{partition_init = "similarity"}
#'   only: features below this standard-deviation threshold are treated as
#'   low-information and get zero off-diagonal similarity.
#' @param K Number of pseudotime bins.
#' @param rw_q Random-walk order.
#' @param ridge Optional nugget added to the RW precision. Mainly useful for
#'   internal diagnostics that compare intrinsic and properized priors.
#' @param lambda_sd_prior_rate Optional positive rate for an exponential prior
#'   on \code{1 / sqrt(lambda_j^{(m)})}. The default \code{NULL} means no lambda
#'   prior penalty. For backward compatibility, an explicit \code{0} is treated
#'   the same way; it is only an alias for "no penalty" and does not
#'   correspond to a literal exponential prior with rate zero.
#' @param lambda_min,lambda_max Lambda bounds.
#' @param sigma_min,sigma_max Bounds used by the \code{"smooth_fit"}
#'   similarity metric and the warm-start single-ordering CAVI fits.
#' @param discretization Initial discretization scheme used when converting
#'   ordering scores into \code{K} bins. For partition fits, MPCurver enforces a
#'   common \code{K} across orderings and may fall back to equal-width bins if
#'   quantile discretization collapses.
#' @param num_iter Warm-start CAVI sweeps per fit.
#' @param verbose Logical.
#'
#' @return A list with \code{$fits} (list of M cavi objects),
#'   \code{$init_info} (per-ordering method metadata), and
#'   \code{$ordering_similarity} (pairwise ordering correlations for
#'   \code{partition_init = "ordering_methods"}) plus optional
#'   \code{$similarity_init} metadata for \code{partition_init = "similarity"}.
#'   When similarity initialization is used, \code{$similarity_init} stores the
#'   full symmetric similarity matrix \code{S}, the corresponding
#'   \code{distance = 1 - S}, and for \code{similarity_metric = "smooth_fit"}
#'   the raw directional score diagnostics used to build \code{S}, including
#'   directional \code{lambda} and \code{sigma^2} estimates.
#' @noRd
init_m_trajectories_cavi <- function(X,
                                     S = NULL,
                                     M = 2L,
                                     methods = NULL,
                                     pca_components = NULL,
                                     K = NULL,
                                     rw_q = 2L,
                                     ridge = 0,
                                     lambda_sd_prior_rate = NULL,
                                     smooth_fit_lambda_mode = c("optimize", "fixed"),
                                     smooth_fit_lambda_value = 1,
                                     lambda_min = 1e-10,
                                     lambda_max = 1e10,
                                     sigma_min = 1e-10,
                                     sigma_max = 1e10,
                                     discretization = c("quantile", "equal", "kmeans"),
                                     partition_init = c("similarity", "ordering_methods"),
                                     similarity_metric = c("spearman", "pearson", "smooth_fit"),
                                     cluster_linkage = "single",
                                     similarity_min_feature_sd = 1e-8,
                                     num_iter = 5L,
                                     verbose = FALSE) {
  X <- as.matrix(X)
  M <- as.integer(M)
  K_use <- .cavi_resolve_K(X, K)
  discretization <- match.arg(discretization)
  partition_init <- match.arg(partition_init)
  similarity_metric <- match.arg(similarity_metric)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)

  if (identical(partition_init, "similarity")) {
    return(.cavi_init_m_trajectories_similarity(
      X = X,
      S = S,
      M = M,
      methods = methods,
      pca_components = pca_components,
      K = K_use,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      smooth_fit_lambda_mode = smooth_fit_lambda_mode,
      smooth_fit_lambda_value = smooth_fit_lambda_value,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      discretization = discretization,
      num_iter = num_iter,
      similarity_metric = similarity_metric,
      cluster_linkage = cluster_linkage,
      similarity_min_feature_sd = similarity_min_feature_sd,
      verbose = verbose
    ))
  }

  # Default methods: first = PCA, rest = PCA (with higher components)
  if (is.null(methods)) {
    methods <- rep("PCA", M)
  }
  if (length(methods) != M)
    stop(sprintf("methods must have length M=%d.", M))

  # Default PCA components follow ordering slots, so method vectors like
  # c("fiedler", "PCA") use PC1 for the single PCA ordering, while
  # repeated PCA entries use PC1, PC2, ... in order of appearance.
  if (is.null(pca_components)) {
    pca_components <- .cavi_default_pca_components(methods)
  }
  if (length(pca_components) != M)
    stop(sprintf("pca_components must have length M=%d.", M))

  .cavi_validate_partition_methods(methods, pca_components)

  fits <- vector("list", M)
  ordering_mat <- matrix(NA_real_, nrow(X), M)
  ord_labels <- .cavi_partition_order_labels(M)
  init_info <- vector("list", M)

  for (m in seq_len(M)) {
    method_m <- methods[m]
    fits[[m]] <- .cavi_fit_from_method(
      X = X,
      S = S,
      method = method_m,
      pca_component = pca_components[m],
      K = K_use,
      rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      max_iter = num_iter,
      discretization = discretization,
      strict_K = TRUE,
      verbose = FALSE
    )

    ordering_result <- .cavi_get_ordering_result(X, method = method_m,
                                                 pca_component = pca_components[m])
    if (!is.null(ordering_result)) {
      ordering_mat[, m] <- ordering_result$t
    }
    init_info[[m]] <- list(
      label = ord_labels[m],
      method = method_m,
      pca_component = if (method_m == "PCA") pca_components[m] else NA_integer_,
      K = K_use,
      discretization = discretization
    )
  }

  names(init_info) <- ord_labels
  colnames(ordering_mat) <- ord_labels
  .cavi_warn_if_similar_orderings(ordering_mat, ord_labels)
  ordering_similarity <- suppressWarnings(stats::cor(ordering_mat, use = "pairwise.complete.obs"))

  if (isTRUE(verbose)) {
    cat(sprintf("[init_m_cavi] M=%d methods=%s\n",
                M, paste(methods, collapse = ",")))
  }

  list(
    fits = fits,
    init_info = init_info,
    ordering_similarity = ordering_similarity,
    similarity_init = NULL
  )
}


# Internal: one step of soft partition (score → softmax → weighted CAVI → rescore)
.soft_partition_step <- function(fits, X, T_now, inner_iter,
                                 lambda_min, lambda_max,
                                 sigma_min = 1e-10,
                                 sigma_max = 1e10,
                                 active_orderings = NULL,
                                 active_feature_pairs = NULL,
                                 weights_prev = NULL,
                                 position_prior = c("adaptive", "fixed"),
                                 freeze_unused_ordering = TRUE,
                                 freeze_unused_ordering_threshold = 0.5,
                                 freeze_feature = TRUE,
                                 freeze_feature_weight_threshold = .cavi_partition_feature_freeze_threshold_default(),
                                 drop_unused_ordering = FALSE,
                                 partition_prior = c("adaptive", "fixed"),
                                 partition_prior_init = NULL,
                                 assignment_prior = NULL,
                                 ordering_alpha = NULL,
                                 iter_index = NA_integer_,
                                 ordering_labels = NULL) {
  M <- length(fits)
  d <- ncol(X)
  if (is.null(active_orderings)) {
    active_orderings <- rep(TRUE, M)
  }
  active_orderings <- as.logical(active_orderings)
  if (length(active_orderings) != M) {
    stop("active_orderings must have length equal to length(fits).")
  }
  if (!any(active_orderings)) {
    stop("At least one ordering must remain active.")
  }
  if (is.null(active_feature_pairs)) {
    active_feature_pairs <- matrix(rep(active_orderings, each = d), nrow = d, ncol = M)
  } else {
    active_feature_pairs <- as.matrix(active_feature_pairs)
    if (!all(dim(active_feature_pairs) == c(d, M))) {
      stop("active_feature_pairs must be a d x M matrix.")
    }
  }
  active_feature_pairs <- active_feature_pairs & matrix(rep(active_orderings, each = d), nrow = d, ncol = M)
  if (is.null(ordering_labels)) {
    ordering_labels <- .cavi_partition_order_labels(M)
  }
  if (length(ordering_labels) != M) {
    stop("ordering_labels must have length equal to length(fits).")
  }
  position_prior <- match.arg(position_prior)
  assignment_ctl <- .validate_partition_assignment_controls(
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    M = M,
    caller = ".soft_partition_step()"
  )
  inactive_ctl <- .validate_partition_inactive_controls(
    freeze_unused_ordering = freeze_unused_ordering,
    freeze_unused_ordering_threshold = freeze_unused_ordering_threshold,
    drop_unused_ordering = drop_unused_ordering
  )
  feature_ctl <- .validate_partition_feature_controls(
    freeze_feature = freeze_feature,
    freeze_feature_weight_threshold = freeze_feature_weight_threshold
  )
  if (is.null(weights_prev)) {
    weights_prev <- matrix(1 / M, nrow = d, ncol = M)
  } else {
    weights_prev <- as.matrix(weights_prev)
    if (!all(dim(weights_prev) == c(d, M))) {
      stop("weights_prev must be a d x M matrix.")
    }
  }
  colnames(weights_prev) <- ordering_labels
  active_now <- active_orderings
  active_pairs_now <- active_feature_pairs
  ordering_events <- list()
  feature_events <- list()
  weights_prev_work <- weights_prev

  if (identical(assignment_ctl$assignment_mode, "adaptive")) {
    zero_drop_pre_seed <- .partition_zero_mass_ordering_drop(
      weights = weights_prev_work,
      active_orderings = active_now,
      active_feature_pairs = active_pairs_now,
      iter_index = iter_index,
      T_now = T_now,
      ordering_labels = ordering_labels,
      stage = "pre_update"
    )
    active_now <- zero_drop_pre_seed$active_orderings
    active_pairs_now[, !active_now] <- FALSE
    ordering_events <- c(ordering_events, zero_drop_pre_seed$ordering_events)
    if (!inactive_ctl$drop_unused_ordering && length(zero_drop_pre_seed$drop_idx) > 0L) {
      weights_prev_work <- .partition_register_new_frozen_weights(
        reported_prev = weights_prev_work,
        snapshot_weights = weights_prev,
        drop_idx = zero_drop_pre_seed$drop_idx
      )
    }
  }

  # Score each fit using the expected log-likelihood block and the current
  # assignment prior state implied by the previous weights.
  objective_pre <- .cavi_partition_objective_from_fits(
    fits = fits,
    X = X,
    weights = weights_prev_work,
    T_now = T_now,
    active_orderings = active_now,
    active_feature_pairs = active_pairs_now,
    assignment_M = M,
    partition_prior = assignment_ctl$partition_prior,
    partition_prior_init = assignment_ctl$partition_prior_init,
    assignment_prior = assignment_ctl$legacy_assignment_prior,
    ordering_alpha = assignment_ctl$ordering_alpha,
    drop_unused_ordering = inactive_ctl$drop_unused_ordering
  )
  score_mat <- objective_pre$like_mat
  score_mat_aug <- sweep(score_mat, 2L, objective_pre$assignment_info$e_log_omega, `+`)

  # Softmax at temperature using only the active orderings.
  effective_weights_candidate <- .soft_partition_softmax(
    score_mat_aug,
    T_now = T_now,
    active_orderings = active_now,
    active_feature_pairs = active_pairs_now
  )
  reported_weights_candidate <- .partition_compose_reported_weights(
    effective_weights = effective_weights_candidate,
    active_feature_pairs = active_pairs_now,
    reported_prev = weights_prev_work,
    drop_unused_ordering = inactive_ctl$drop_unused_ordering
  )

  if (feature_ctl$freeze_feature) {
    freeze_feature_pre <- .partition_freeze_features(
      weights = effective_weights_candidate,
      active_feature_pairs = active_pairs_now,
      active_orderings = active_now,
      threshold = feature_ctl$freeze_feature_weight_threshold,
      iter_index = iter_index,
      T_now = T_now,
      ordering_labels = ordering_labels,
      stage = "pre_update"
    )
    active_pairs_now <- freeze_feature_pre$active_feature_pairs
    feature_events <- c(feature_events, freeze_feature_pre$feature_events)
    if (!inactive_ctl$drop_unused_ordering && any(freeze_feature_pre$new_frozen_mask)) {
      weights_prev_work <- .partition_register_new_frozen_weights(
        reported_prev = weights_prev_work,
        snapshot_weights = reported_weights_candidate,
        frozen_mask = freeze_feature_pre$new_frozen_mask
      )
    }
    if (any(freeze_feature_pre$new_frozen_mask)) {
      effective_weights_candidate <- .soft_partition_softmax(
        score_mat_aug,
        T_now = T_now,
        active_orderings = active_now,
        active_feature_pairs = active_pairs_now
      )
      reported_weights_candidate <- .partition_compose_reported_weights(
        effective_weights = effective_weights_candidate,
        active_feature_pairs = active_pairs_now,
        reported_prev = weights_prev_work,
        drop_unused_ordering = inactive_ctl$drop_unused_ordering
      )
    }
  }

  if (identical(assignment_ctl$assignment_mode, "adaptive")) {
    zero_drop_pre <- .partition_zero_mass_ordering_drop(
      weights = effective_weights_candidate,
      active_orderings = active_now,
      active_feature_pairs = active_pairs_now,
      iter_index = iter_index,
      T_now = T_now,
      ordering_labels = ordering_labels,
      stage = "pre_update"
    )
    active_now <- zero_drop_pre$active_orderings
    active_pairs_now[, !active_now] <- FALSE
    ordering_events <- c(ordering_events, zero_drop_pre$ordering_events)
    if (!inactive_ctl$drop_unused_ordering && length(zero_drop_pre$drop_idx) > 0L) {
      weights_prev_work <- .partition_register_new_frozen_weights(
        reported_prev = weights_prev_work,
        snapshot_weights = reported_weights_candidate,
        drop_idx = zero_drop_pre$drop_idx
      )
    }
    if (length(zero_drop_pre$drop_idx) > 0L) {
      effective_weights_candidate <- .soft_partition_softmax(
        score_mat_aug,
        T_now = T_now,
        active_orderings = active_now,
        active_feature_pairs = active_pairs_now
      )
      reported_weights_candidate <- .partition_compose_reported_weights(
        effective_weights = effective_weights_candidate,
        active_feature_pairs = active_pairs_now,
        reported_prev = weights_prev_work,
        drop_unused_ordering = inactive_ctl$drop_unused_ordering
      )
    }
  }

  if (inactive_ctl$freeze_unused_ordering) {
    freeze_pre <- .partition_freeze_orderings(
      weights = effective_weights_candidate,
      active_orderings = active_now,
      threshold = inactive_ctl$freeze_unused_ordering_threshold,
      iter_index = iter_index,
      T_now = T_now,
      ordering_labels = ordering_labels,
      stage = "pre_update"
    )
    active_now <- freeze_pre$active_orderings
    ordering_events <- c(ordering_events, freeze_pre$ordering_events)
    active_pairs_now[, !active_now] <- FALSE
    weights_prev_pre <- weights_prev_work
    if (!inactive_ctl$drop_unused_ordering && length(freeze_pre$drop_idx) > 0L) {
      weights_prev_pre <- .partition_register_new_frozen_weights(
        reported_prev = weights_prev_pre,
        snapshot_weights = reported_weights_candidate,
        drop_idx = freeze_pre$drop_idx
      )
    }
    if (length(freeze_pre$drop_idx) > 0L) {
      effective_weights <- .soft_partition_softmax(
        score_mat_aug,
        T_now = T_now,
        active_orderings = active_now,
        active_feature_pairs = active_pairs_now
      )
    } else {
      effective_weights <- effective_weights_candidate
    }
    weights <- .partition_compose_reported_weights(
      effective_weights = effective_weights,
      active_feature_pairs = active_pairs_now,
      reported_prev = weights_prev_pre,
      drop_unused_ordering = inactive_ctl$drop_unused_ordering
    )
  } else {
    effective_weights <- effective_weights_candidate
    weights <- reported_weights_candidate
  }

  # Weighted CAVI update on each fit
  for (m in which(active_now)) {
    fits[[m]]$control$position_prior <- position_prior
    fits[[m]] <- .do_cavi_weighted(
      object = fits[[m]],
      data = X,
      feature_weights = effective_weights[, m],
      feature_active = active_pairs_now[, m],
      iter = inner_iter,
      lambda_sd_prior_rate = (fits[[m]]$control %||% list())$lambda_sd_prior_rate %||% NULL,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      tol = 0,
      verbose = FALSE
    )
  }

  # Re-score and re-weight post-update
  objective_post_seed <- .cavi_partition_objective_from_fits(
    fits = fits,
    X = X,
    weights = weights,
    T_now = T_now,
    active_orderings = active_now,
    active_feature_pairs = active_pairs_now,
    assignment_M = M,
    partition_prior = assignment_ctl$partition_prior,
    partition_prior_init = assignment_ctl$partition_prior_init,
    assignment_prior = assignment_ctl$legacy_assignment_prior,
    ordering_alpha = assignment_ctl$ordering_alpha,
    drop_unused_ordering = inactive_ctl$drop_unused_ordering
  )
  score_mat_post <- objective_post_seed$like_mat
  score_mat_post_aug <- sweep(
    score_mat_post,
    2L,
    objective_post_seed$assignment_info$e_log_omega,
    `+`
  )
  effective_weights_post_candidate <- .soft_partition_softmax(
    score_mat_post_aug,
    T_now = T_now,
    active_orderings = active_now,
    active_feature_pairs = active_pairs_now
  )
  reported_weights_post_candidate <- .partition_compose_reported_weights(
    effective_weights = effective_weights_post_candidate,
    active_feature_pairs = active_pairs_now,
    reported_prev = weights,
    drop_unused_ordering = inactive_ctl$drop_unused_ordering
  )
  weights_prev_post_work <- weights

  if (feature_ctl$freeze_feature) {
    freeze_feature_post <- .partition_freeze_features(
      weights = effective_weights_post_candidate,
      active_feature_pairs = active_pairs_now,
      active_orderings = active_now,
      threshold = feature_ctl$freeze_feature_weight_threshold,
      iter_index = iter_index,
      T_now = T_now,
      ordering_labels = ordering_labels,
      stage = "post_update"
    )
    active_pairs_now <- freeze_feature_post$active_feature_pairs
    feature_events <- c(feature_events, freeze_feature_post$feature_events)
    if (!inactive_ctl$drop_unused_ordering && any(freeze_feature_post$new_frozen_mask)) {
      weights_prev_post_work <- .partition_register_new_frozen_weights(
        reported_prev = weights_prev_post_work,
        snapshot_weights = reported_weights_post_candidate,
        frozen_mask = freeze_feature_post$new_frozen_mask
      )
    }
    if (any(freeze_feature_post$new_frozen_mask)) {
      effective_weights_post_candidate <- .soft_partition_softmax(
        score_mat_post_aug,
        T_now = T_now,
        active_orderings = active_now,
        active_feature_pairs = active_pairs_now
      )
      reported_weights_post_candidate <- .partition_compose_reported_weights(
        effective_weights = effective_weights_post_candidate,
        active_feature_pairs = active_pairs_now,
        reported_prev = weights_prev_post_work,
        drop_unused_ordering = inactive_ctl$drop_unused_ordering
      )
    }
  }

  if (identical(assignment_ctl$assignment_mode, "adaptive")) {
    zero_drop_post <- .partition_zero_mass_ordering_drop(
      weights = effective_weights_post_candidate,
      active_orderings = active_now,
      active_feature_pairs = active_pairs_now,
      iter_index = iter_index,
      T_now = T_now,
      ordering_labels = ordering_labels,
      stage = "post_update"
    )
    active_now <- zero_drop_post$active_orderings
    active_pairs_now[, !active_now] <- FALSE
    ordering_events <- c(ordering_events, zero_drop_post$ordering_events)
    if (!inactive_ctl$drop_unused_ordering && length(zero_drop_post$drop_idx) > 0L) {
      weights_prev_post_work <- .partition_register_new_frozen_weights(
        reported_prev = weights_prev_post_work,
        snapshot_weights = reported_weights_post_candidate,
        drop_idx = zero_drop_post$drop_idx
      )
    }
    if (length(zero_drop_post$drop_idx) > 0L) {
      effective_weights_post_candidate <- .soft_partition_softmax(
        score_mat_post_aug,
        T_now = T_now,
        active_orderings = active_now,
        active_feature_pairs = active_pairs_now
      )
      reported_weights_post_candidate <- .partition_compose_reported_weights(
        effective_weights = effective_weights_post_candidate,
        active_feature_pairs = active_pairs_now,
        reported_prev = weights_prev_post_work,
        drop_unused_ordering = inactive_ctl$drop_unused_ordering
      )
    }
  }

  if (inactive_ctl$freeze_unused_ordering) {
    freeze_post <- .partition_freeze_orderings(
      weights = effective_weights_post_candidate,
      active_orderings = active_now,
      threshold = inactive_ctl$freeze_unused_ordering_threshold,
      iter_index = iter_index,
      T_now = T_now,
      ordering_labels = ordering_labels,
      stage = "post_update"
    )
    active_now <- freeze_post$active_orderings
    ordering_events <- c(ordering_events, freeze_post$ordering_events)
    active_pairs_now[, !active_now] <- FALSE
    weights_prev_post <- weights_prev_post_work
    if (!inactive_ctl$drop_unused_ordering && length(freeze_post$drop_idx) > 0L) {
      weights_prev_post <- .partition_register_new_frozen_weights(
        reported_prev = weights_prev_post,
        snapshot_weights = reported_weights_post_candidate,
        drop_idx = freeze_post$drop_idx
      )
    }
    if (length(freeze_post$drop_idx) > 0L) {
      effective_weights_post <- .soft_partition_softmax(
        score_mat_post_aug,
        T_now = T_now,
        active_orderings = active_now,
        active_feature_pairs = active_pairs_now
      )
    } else {
      effective_weights_post <- effective_weights_post_candidate
    }
    weights_post <- .partition_compose_reported_weights(
      effective_weights = effective_weights_post,
      active_feature_pairs = active_pairs_now,
      reported_prev = weights_prev_post,
      drop_unused_ordering = inactive_ctl$drop_unused_ordering
    )
  } else {
    effective_weights_post <- effective_weights_post_candidate
    weights_post <- reported_weights_post_candidate
  }

  # Refresh the assignment state so the returned weights/objective are aligned
  # with the final prior-augmented coordinate updates.
  objective_final_seed <- .cavi_partition_objective_from_fits(
    fits = fits,
    X = X,
    weights = weights_post,
    T_now = T_now,
    active_orderings = active_now,
    active_feature_pairs = active_pairs_now,
    assignment_M = M,
    partition_prior = assignment_ctl$partition_prior,
    partition_prior_init = assignment_ctl$partition_prior_init,
    assignment_prior = assignment_ctl$legacy_assignment_prior,
    ordering_alpha = assignment_ctl$ordering_alpha,
    drop_unused_ordering = inactive_ctl$drop_unused_ordering
  )
  score_mat_final_aug <- sweep(
    objective_final_seed$like_mat,
    2L,
    objective_final_seed$assignment_info$e_log_omega,
    `+`
  )
  effective_weights_final <- .soft_partition_softmax(
    score_mat_final_aug,
    T_now = T_now,
    active_orderings = active_now,
    active_feature_pairs = active_pairs_now
  )
  weights_final <- .partition_compose_reported_weights(
    effective_weights = effective_weights_final,
    active_feature_pairs = active_pairs_now,
    reported_prev = weights_post,
    drop_unused_ordering = inactive_ctl$drop_unused_ordering
  )
  objective_post <- .cavi_partition_objective_from_fits(
    fits = fits,
    X = X,
    weights = weights_final,
    T_now = T_now,
    active_orderings = active_now,
    active_feature_pairs = active_pairs_now,
    assignment_M = M,
    partition_prior = assignment_ctl$partition_prior,
    partition_prior_init = assignment_ctl$partition_prior_init,
    assignment_prior = assignment_ctl$legacy_assignment_prior,
    ordering_alpha = assignment_ctl$ordering_alpha,
    drop_unused_ordering = inactive_ctl$drop_unused_ordering
  )
  freeze_happened <- (length(ordering_events) + length(feature_events)) > 0L

  list(
    fits = fits,
    score_mat = objective_post$like_mat,
    pi_weights = weights_final,
    effective_pi_weights = effective_weights_final,
    objective = objective_post$objective,
    active_orderings = active_now,
    active_feature_pairs = active_pairs_now,
    frozen_orderings = !active_now,
    ordering_events = ordering_events,
    feature_events = feature_events,
    freeze_happened = freeze_happened,
    assignment_info = objective_post$assignment_info,
    assignment_posterior = .cavi_partition_assignment_state(objective_post$assignment_info)
  )
}

.partition_convergence_defaults <- function() {
  list(
    tol_outer = 1e-5,
    min_phase2_iter = 3L,
    consecutive_small_steps = 3L
  )
}

.partition_init_convergence_info <- function(tol_outer = .partition_convergence_defaults()$tol_outer,
                                             phase2_iters = 0L,
                                             consecutive_small_steps = 0L,
                                             converged = FALSE,
                                             last_delta = NA_real_,
                                             last_rel_delta = NA_real_,
                                             reason = "") {
  defaults <- .partition_convergence_defaults()
  tol_outer <- as.numeric(tol_outer)[1]
  if (!is.finite(tol_outer)) {
    stop("tol_outer must be a finite numeric value.")
  }
  list(
    criterion = sprintf(
      "phase-2 rel_delta < %.1e for %d consecutive nonnegative steps after at least %d phase-2 iterations",
      tol_outer,
      defaults$consecutive_small_steps,
      defaults$min_phase2_iter
    ),
    tol_outer = tol_outer,
    phase2_iters = as.integer(phase2_iters),
    last_delta = as.numeric(last_delta),
    last_rel_delta = as.numeric(last_rel_delta),
    consecutive_small_steps = as.integer(consecutive_small_steps),
    min_phase2_iter = defaults$min_phase2_iter,
    consecutive_required = defaults$consecutive_small_steps,
    converged = isTRUE(converged),
    reason = as.character(reason %||% "")[1]
  )
}

.partition_update_convergence_info <- function(info, obj_prev, obj_now,
                                              freeze_step = FALSE) {
  if (is.null(info)) {
    info <- .partition_init_convergence_info()
  }
  delta <- obj_now - obj_prev
  rel_delta <- if (is.finite(obj_prev)) delta / (abs(obj_prev) + 1) else NA_real_

  info$last_delta <- delta
  info$last_rel_delta <- rel_delta
  info$converged <- FALSE

  if (isTRUE(freeze_step)) {
    info$reason <- "Phase-2 step skipped for convergence because freezing occurred."
    return(info)
  }

  info$phase2_iters <- as.integer(info$phase2_iters %||% 0L) + 1L

  if (!is.finite(obj_prev)) {
    info$consecutive_small_steps <- 0L
    info$reason <- "Phase-2 convergence check waiting for a finite previous objective."
    return(info)
  }

  if (info$tol_outer <= 0) {
    info$consecutive_small_steps <- 0L
    info$reason <- "tol_outer <= 0; early stopping disabled."
    return(info)
  }

  is_small_step <- is.finite(rel_delta) && delta >= 0 && abs(rel_delta) < info$tol_outer
  if (is_small_step) {
    info$consecutive_small_steps <- as.integer(info$consecutive_small_steps %||% 0L) + 1L
  } else {
    info$consecutive_small_steps <- 0L
  }

  if (delta < 0) {
    info$reason <- sprintf("Objective decreased by %.3e at T=1; convergence not declared.", delta)
    return(info)
  }

  if (info$phase2_iters >= info$min_phase2_iter &&
      info$consecutive_small_steps >= info$consecutive_required) {
    info$converged <- TRUE
    info$reason <- sprintf(
      "Phase-2 rel_delta stayed below %.1e for %d consecutive nonnegative steps.",
      info$tol_outer,
      info$consecutive_required
    )
  } else {
    info$reason <- sprintf(
      "Phase-2 rel_delta=%.3e with streak %d/%d.",
      rel_delta,
      info$consecutive_small_steps,
      info$consecutive_required
    )
  }

  info
}


#' Soft M-ordering partition using mean-field CAVI
#'
#' @description
#' Generalised soft partition routine supporting \eqn{M \ge 2} orderings.
#' Features are assigned probabilistic weights across orderings via
#' temperature-annealed softmax on per-feature ELBO scores.
#'
#' @param X Numeric matrix (\code{n x d}).
#' @param M Integer >= 2. Number of orderings.
#' @param fits_init Optional list of M \code{cavi} fits. If NULL, initialized
#'   via \code{init_m_trajectories_cavi()}.
#' @param init_methods Character vector of length M. Ordering initialisation
#'   methods (e.g., \code{c("fiedler", "PCA", "PCA")}). If NULL, defaults to
#'   \code{rep("PCA", M)}.
#' @param pca_components Integer vector of length M. Which PCA component to
#'   use for PCA-based methods. If NULL, auto-assigned sequentially.
#' @param partition_init Partition initialisation strategy. Use
#'   \code{"similarity"} to cluster features into \code{M} blocks via
#'   \code{1 - S(X)} before building ordering-specific warm starts, or
#'   \code{"ordering_methods"} for the existing per-ordering method starts.
#'   The default is \code{"similarity"}.
#' @param similarity_metric For \code{partition_init = "similarity"} only:
#'   pairwise feature similarity metric used to form \code{S(X)}. One of
#'   \code{"spearman"}, \code{"pearson"}, or \code{"smooth_fit"}.
#'   \code{"smooth_fit"} can be used with \code{ridge = 0}, but for
#'   \code{M > 1} this uses an intrinsic-RW pseudo-evidence rather than a
#'   fully proper marginal likelihood; use a small positive \code{ridge} if
#'   you want the smoother evidence to be theoretically proper.
#' @param smooth_fit_lambda_mode For \code{similarity_metric = "smooth_fit"}
#'   only: whether the directional smoother optimizes \code{lambda} or keeps it
#'   fixed at \code{smooth_fit_lambda_value}.
#' @param smooth_fit_lambda_value For \code{similarity_metric = "smooth_fit"}
#'   only: fixed \code{lambda} value used when
#'   \code{smooth_fit_lambda_mode = "fixed"}, and the starting value when
#'   \code{smooth_fit_lambda_mode = "optimize"}.
#' @param cluster_linkage For \code{partition_init = "similarity"} only:
#'   hierarchical-clustering linkage passed to \code{\link[stats]{hclust}}.
#'   The default is \code{"single"}.
#' @param similarity_min_feature_sd For \code{partition_init = "similarity"}
#'   only: features with standard deviation below this threshold are treated as
#'   low-information and assigned zero off-diagonal similarity.
#' @param K Number of pseudotime bins.
#' @param discretization Initial discretization scheme used when converting
#'   ordering scores into \code{K} bins. Partition fits require a common
#'   \code{K} across orderings, so if quantile discretization collapses the
#'   number of groups, MPCurver falls back to equal-width bins.
#' @param T_start,T_end Annealing temperatures.
#' @param n_outer Number of annealing steps.
#' @param inner_iter Number of weighted CAVI sweeps per outer step.
#' @param max_converge_iter Maximum number of post-annealing iterations at
#'   \code{T = T_end}.
#' @param tol_outer Relative objective tolerance for the phase-2 convergence
#'   rule. Convergence is declared only after at least three exact
#'   \code{T = T_end} iterations and three consecutive nonnegative relative
#'   objective improvements below this threshold.
#' @param rw_q Random-walk order.
#' @param ridge Optional nugget added to the RW precision for all orderings.
#' @param lambda_sd_prior_rate Optional positive rate for an exponential prior
#'   on \code{1 / sqrt(lambda_j^{(m)})}. The default \code{NULL} means no lambda
#'   prior penalty. For backward compatibility, an explicit \code{0} is treated
#'   the same way; it is only an alias for "no penalty" and does not
#'   correspond to a literal exponential prior with rate zero.
#' @param lambda_min,lambda_max Lambda bounds.
#' @param sigma_min,sigma_max Bounds for \code{sigma_j^2}. These are used both
#'   in the weighted CAVI updates and by the \code{"smooth_fit"} similarity
#'   metric when \code{partition_init = "similarity"}.
#' @param position_prior Either \code{"adaptive"} or \code{"fixed"} for the
#'   within-ordering component prior \code{pi}. When fixed, the weighted CAVI
#'   updates keep each ordering's \code{pi} constant.
#' @param position_prior_init Optional length-\code{K} initial/fixed vector for
#'   the within-ordering component prior \code{pi}. When
#'   \code{position_prior = "fixed"} and this is \code{NULL}, a uniform prior
#'   over the \code{K} positions is used. A single vector is broadcast to all
#'   orderings.
#' @param partition_prior Either \code{"adaptive"} or \code{"fixed"} for the
#'   ordering-usage prior \code{omega}. The adaptive mode uses an empirical-
#'   Bayes update over the active orderings; exact zero-mass orderings are
#'   dropped from the active set before the next prior/objective evaluation.
#' @param partition_prior_init Optional length-\code{M} initial/fixed vector for
#'   the ordering prior \code{omega}. Used only when
#'   \code{partition_prior = "fixed"}. When omitted, the fixed prior defaults
#'   to uniform over the active orderings.
#' @param assignment_prior Deprecated compatibility argument for the old
#'   partition-prior API. \code{"uniform"} maps to
#'   \code{partition_prior = "fixed"} with a uniform prior. \code{"dirichlet"}
#'   remains available only as a deprecated compatibility path.
#' @param ordering_alpha Deprecated compatibility argument used only for the
#'   legacy \code{assignment_prior = "dirichlet"} path.
#' @param hard_assign_final Logical; make hard assignments at the end?
#' @param freeze_unused_ordering Logical; if \code{TRUE}, orderings whose
#'   feature-posterior mass falls below
#'   \code{freeze_unused_ordering_threshold} are frozen and no longer updated.
#' @param freeze_unused_ordering_threshold Non-negative threshold on the summed
#'   feature-posterior mass \code{sum(pi_weights[, m])} used to decide whether
#'   an ordering is effectively unused.
#' @param freeze_feature Logical; if \code{TRUE}, feature-ordering pairs with
#'   posterior weight below \code{freeze_feature_weight_threshold} are frozen,
#'   their trajectory \code{U}-block is no longer updated, and their
#'   feature-trajectory contribution is neutralized to zero.
#' @param freeze_feature_weight_threshold Non-negative threshold on
#'   \code{w_{jm}} used to decide whether a feature-ordering pair is
#'   effectively unused.
#' @param drop_unused_ordering Logical; if \code{TRUE}, the user-facing
#'   \code{mpcurve} wrapper returned by \code{fit_mpcurve()} drops frozen
#'   orderings from its displayed partition view. Internally, fitting still
#'   proceeds by freezing rather than deleting orderings. When
#'   \code{drop_unused_ordering = FALSE}, the exposed \code{$objective_history}
#'   is a fixed-requested-\code{M} comparison objective: frozen orderings keep
#'   their preserved feature-assignment weights and cell-ordering block, while
#'   their trajectory \code{U}-block contribution is neutralized to zero. When
#'   \code{drop_unused_ordering = TRUE}, the exposed \code{$objective_history}
#'   is the post-drop fitting objective of the active model and should not be
#'   compared across requested values of \code{M}. Fixed-\code{M} fits may
#'   still contain frozen orderings; greedy model selection is responsible for
#'   normalizing cross-\code{M} comparisons to the active dimension.
#' @param verbose Logical.
#'
#' @return An object of class \code{"soft_partition_cavi"} with fields:
#'   \code{$fits} (list of M cavi fits), \code{$pi_weights} (d x M matrix),
#'   \code{$assign} (character vector), \code{$M}, \code{$objective_history},
#'   \code{$score_history}, \code{$weight_history}, \code{$T_schedule},
#'   \code{$n_anneal}, \code{$converged}, \code{$convergence_info},
#'   \code{$active_orderings},
#'   \code{$active_feature_pairs}, \code{$frozen_orderings},
#'   \code{$ordering_events}, \code{$feature_events}, and \code{$control},
#'   plus optional initialisation diagnostics \code{$init_info} and
#'   \code{$ordering_similarity}, and optional similarity metadata
#'   \code{$similarity_init}. When similarity initialisation is used,
#'   \code{$similarity_init} stores the full symmetric similarity matrix
#'   \code{S}, \code{distance = 1 - S}, and for
#'   \code{similarity_metric = "smooth_fit"} the raw directional score
#'   diagnostics used to construct \code{S}, including directional
#'   \code{lambda} and \code{sigma^2} estimates. The exposed
#'   \code{$objective_history} has two user-facing semantics. If
#'   \code{drop_unused_ordering = FALSE}, it is the fixed-requested-\code{M}
#'   comparison objective intended for comparing fits across requested
#'   \code{M}; frozen orderings may still be present in that fixed-\code{M}
#'   state. If \code{drop_unused_ordering = TRUE}, it is the post-drop fitting
#'   objective of the active model and is not intended for cross-\code{M}
#'   comparison.
#' @noRd
soft_partition_cavi <- function(X,
                                S = NULL,
                                M = 2L,
                                fits_init = NULL,
                                init_methods = NULL,
                                pca_components = NULL,
                                partition_init = c("similarity", "ordering_methods"),
                                similarity_metric = c("spearman", "pearson", "smooth_fit"),
                                smooth_fit_lambda_mode = c("optimize", "fixed"),
                                smooth_fit_lambda_value = 1,
                                cluster_linkage = "single",
                                similarity_min_feature_sd = 1e-8,
                                K = NULL,
                                discretization = c("quantile", "equal", "kmeans"),
                                T_start = 5,
                                T_end = 1,
                                n_outer = 25L,
                                inner_iter = 1L,
                                max_converge_iter = 100L,
                                tol_outer = 1e-5,
                                rw_q = 2L,
                                ridge = 0,
                                lambda_sd_prior_rate = NULL,
                                lambda_min = 1e-10,
                                lambda_max = 1e10,
                                sigma_min = 1e-10,
                                sigma_max = 1e10,
                                position_prior = c("adaptive", "fixed"),
                                position_prior_init = NULL,
                                partition_prior = c("adaptive", "fixed"),
                                partition_prior_init = NULL,
                                assignment_prior = NULL,
                                ordering_alpha = NULL,
                                hard_assign_final = FALSE,
                                freeze_unused_ordering = TRUE,
                                freeze_unused_ordering_threshold = 0.5,
                                freeze_feature = TRUE,
                                freeze_feature_weight_threshold = 0.1,
                                drop_unused_ordering = FALSE,
                                verbose = TRUE) {
  X <- as.matrix(X)
  M <- as.integer(M)
  d <- ncol(X)
  partition_init <- match.arg(partition_init)
  discretization <- match.arg(discretization)
  similarity_metric <- match.arg(similarity_metric)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)
  position_prior <- match.arg(position_prior)
  partition_prior_missing <- missing(partition_prior)
  partition_prior <- match.arg(partition_prior)
  cluster_linkage <- .cavi_validate_cluster_linkage(cluster_linkage)
  similarity_min_feature_sd <- as.numeric(similarity_min_feature_sd)[1]
  lambda_sd_prior_rate <- .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
  tol_outer <- as.numeric(tol_outer)[1]
  if (M < 2L) stop("M must be >= 2 for partition models.")
  if (!is.finite(tol_outer)) stop("tol_outer must be a finite numeric value.")
  if (!is.finite(similarity_min_feature_sd) || similarity_min_feature_sd < 0) {
    stop("similarity_min_feature_sd must be a single finite nonnegative number.")
  }
  assignment_ctl <- .validate_partition_assignment_controls(
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    M = M,
    partition_prior_missing = partition_prior_missing,
    caller = "soft_partition_cavi()"
  )
  inactive_ctl <- .validate_partition_inactive_controls(
    freeze_unused_ordering = freeze_unused_ordering,
    freeze_unused_ordering_threshold = freeze_unused_ordering_threshold,
    drop_unused_ordering = drop_unused_ordering
  )
  feature_ctl <- .validate_partition_feature_controls(
    freeze_feature = freeze_feature,
    freeze_feature_weight_threshold = freeze_feature_weight_threshold
  )

  # Ordering labels
  ord_labels <- .cavi_partition_order_labels(M)
  init_info <- NULL
  ordering_similarity <- NULL
  similarity_init <- NULL
  applied_partition_init <- if (is.null(fits_init)) partition_init else "fits_init"

  # Initialize fits
  if (is.null(fits_init)) {
    inits <- init_m_trajectories_cavi(
      X = X, S = S, M = M,
      methods = init_methods,
      pca_components = pca_components,
      partition_init = partition_init,
      similarity_metric = similarity_metric,
      smooth_fit_lambda_mode = smooth_fit_lambda_mode,
      smooth_fit_lambda_value = smooth_fit_lambda_value,
      cluster_linkage = cluster_linkage,
      similarity_min_feature_sd = similarity_min_feature_sd,
      K = K, rw_q = rw_q,
      ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      sigma_min = sigma_min,
      sigma_max = sigma_max,
      discretization = discretization,
      num_iter = 5L,
      verbose = verbose
    )
    fits <- inits$fits
    init_info <- inits$init_info
    ordering_similarity <- inits$ordering_similarity
    similarity_init <- inits$similarity_init %||% NULL
  } else {
    if (length(fits_init) != M)
      stop(sprintf("fits_init must be a list of length M=%d.", M))
    for (m in seq_len(M)) {
      if (!inherits(fits_init[[m]], "cavi"))
        stop(sprintf("fits_init[[%d]] must be a cavi object.", m))
    }
    fits <- fits_init
  }

  stored_partition_S <- fits[[1]]$measurement_sd %||% NULL
  if (is.null(S)) {
    partition_measurement_sd <- stored_partition_S
  } else {
    if (!is.null(stored_partition_S) && !.cavi_same_measurement_sd(stored_partition_S, S)) {
      stop("Supplied S does not match the measurement_sd stored on fits_init.", call. = FALSE)
    }
    partition_measurement_sd <- S
  }
  for (m in seq_len(M)) {
    fit_S <- fits[[m]]$measurement_sd %||% NULL
    if (!.cavi_same_measurement_sd(fit_S, partition_measurement_sd)) {
      stop("All fits in fits_init must share the same measurement_sd specification.", call. = FALSE)
    }
  }
  partition_noise_info <- .cavi_resolve_measurement_sd(
    S = partition_measurement_sd,
    X = X,
    caller = "soft_partition_cavi()"
  )

  fits <- lapply(fits, function(fit) {
    fit$measurement_sd <- partition_measurement_sd
    fit$control <- utils::modifyList(
      fit$control %||% list(),
      list(
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        noise_model = partition_noise_info$noise_model
      )
    )
    fit
  })

  K_vals <- vapply(fits, function(f) length(f$params$pi), integer(1))
  if (length(unique(K_vals)) != 1L) {
    stop("All orderings in a partition fit must use the same K. ",
         "Please re-initialize with a common K or let soft_partition_cavi() build the fits.")
  }
  if (!is.null(K) && unique(K_vals) != as.integer(K)) {
    warning("fits_init already determines K=", unique(K_vals),
            "; ignoring the requested K=", as.integer(K), ".", call. = FALSE)
  }
  K_common <- unique(K_vals)

  position_prior_ctl <- .cavi_validate_position_prior(
    position_prior = position_prior,
    position_prior_init = position_prior_init,
    K = K_common,
    caller = "soft_partition_cavi()"
  )

  fits <- lapply(fits, function(fit) {
    fit$control <- utils::modifyList(
      fit$control %||% list(),
      list(position_prior = position_prior_ctl$position_prior)
    )
    if (!is.null(position_prior_ctl$position_prior_init) ||
        identical(position_prior_ctl$position_prior, "fixed")) {
      fit$params$pi <- .cavi_initialize_position_prior(
        R = fit$gamma,
        K = K_common,
        position_prior = position_prior_ctl$position_prior,
        position_prior_init = position_prior_ctl$position_prior_init
      )
      if (length(fit$pi_trace %||% list())) {
        fit$pi_trace[[length(fit$pi_trace)]] <- fit$params$pi
      }
    }
    fit
  })

  # Initial weights: uniform 1/M
  pi_weights <- matrix(1 / M, nrow = d, ncol = M)
  colnames(pi_weights) <- ord_labels
  effective_pi_weights <- pi_weights

  # Annealing schedule
  T_schedule <- exp(seq(log(T_start), log(T_end), length.out = as.integer(n_outer)))
  total_max <- as.integer(n_outer + max_converge_iter)
  score_history <- vector("list", total_max)
  weight_history <- vector("list", total_max)
  effective_weight_history <- vector("list", total_max)
  objective_history <- numeric(total_max)
  actual_iters <- 0L
  converged <- FALSE
  convergence_info <- .partition_init_convergence_info(tol_outer = tol_outer)
  active_orderings <- rep(TRUE, M)
  active_feature_pairs <- matrix(TRUE, nrow = d, ncol = M)
  ordering_events <- list()
  feature_events <- list()
  assignment_posterior <- NULL

  # Annealing phase
  for (iter in seq_len(as.integer(n_outer))) {
    step <- .soft_partition_step(fits, X, T_schedule[iter], inner_iter,
                                 lambda_min, lambda_max,
                                 sigma_min = sigma_min,
                                 sigma_max = sigma_max,
                                 active_orderings = active_orderings,
                                 active_feature_pairs = active_feature_pairs,
                                 weights_prev = pi_weights,
                                 position_prior = position_prior_ctl$position_prior,
                                 freeze_unused_ordering = inactive_ctl$freeze_unused_ordering,
                                 freeze_unused_ordering_threshold = inactive_ctl$freeze_unused_ordering_threshold,
                                 freeze_feature = feature_ctl$freeze_feature,
                                 freeze_feature_weight_threshold = feature_ctl$freeze_feature_weight_threshold,
                                 drop_unused_ordering = inactive_ctl$drop_unused_ordering,
                                 partition_prior = assignment_ctl$partition_prior,
                                 partition_prior_init = assignment_ctl$partition_prior_init,
                                 assignment_prior = assignment_ctl$legacy_assignment_prior,
                                 ordering_alpha = assignment_ctl$ordering_alpha,
                                 iter_index = iter,
                                 ordering_labels = ord_labels)
    fits <- step$fits
    pi_weights <- step$pi_weights
    effective_pi_weights <- step$effective_pi_weights
    active_orderings <- step$active_orderings
    active_feature_pairs <- step$active_feature_pairs
    colnames(pi_weights) <- ord_labels
    colnames(effective_pi_weights) <- ord_labels
    score_history[[iter]] <- step$score_mat
    weight_history[[iter]] <- pi_weights
    effective_weight_history[[iter]] <- effective_pi_weights
    objective_history[iter] <- step$objective
    actual_iters <- iter
    ordering_events <- c(ordering_events, step$ordering_events)
    feature_events <- c(feature_events, step$feature_events)
    assignment_posterior <- step$assignment_posterior

    if (verbose) {
      cat(sprintf("[soft-partition %3d/%d] T=%.3f obj=%.6f\n",
                  iter, n_outer, T_schedule[iter], step$objective))
    }
  }

  # Convergence phase at T_end
  obj_prev <- if (actual_iters > 0L) objective_history[actual_iters] else -Inf
  for (iter2 in seq_len(as.integer(max_converge_iter))) {
    idx <- as.integer(n_outer) + iter2
    step <- .soft_partition_step(fits, X, T_end, inner_iter,
                                 lambda_min, lambda_max,
                                 sigma_min = sigma_min,
                                 sigma_max = sigma_max,
                                 active_orderings = active_orderings,
                                 active_feature_pairs = active_feature_pairs,
                                 weights_prev = pi_weights,
                                 position_prior = position_prior_ctl$position_prior,
                                 freeze_unused_ordering = inactive_ctl$freeze_unused_ordering,
                                 freeze_unused_ordering_threshold = inactive_ctl$freeze_unused_ordering_threshold,
                                 freeze_feature = feature_ctl$freeze_feature,
                                 freeze_feature_weight_threshold = feature_ctl$freeze_feature_weight_threshold,
                                 drop_unused_ordering = inactive_ctl$drop_unused_ordering,
                                 partition_prior = assignment_ctl$partition_prior,
                                 partition_prior_init = assignment_ctl$partition_prior_init,
                                 assignment_prior = assignment_ctl$legacy_assignment_prior,
                                 ordering_alpha = assignment_ctl$ordering_alpha,
                                 iter_index = idx,
                                 ordering_labels = ord_labels)
    fits <- step$fits
    pi_weights <- step$pi_weights
    effective_pi_weights <- step$effective_pi_weights
    active_orderings <- step$active_orderings
    active_feature_pairs <- step$active_feature_pairs
    colnames(pi_weights) <- ord_labels
    colnames(effective_pi_weights) <- ord_labels
    score_history[[idx]] <- step$score_mat
    weight_history[[idx]] <- pi_weights
    effective_weight_history[[idx]] <- effective_pi_weights
    objective_history[idx] <- step$objective
    actual_iters <- idx
    ordering_events <- c(ordering_events, step$ordering_events)
    feature_events <- c(feature_events, step$feature_events)
    assignment_posterior <- step$assignment_posterior

    delta <- step$objective - obj_prev
    convergence_info <- .partition_update_convergence_info(
      info = convergence_info,
      obj_prev = obj_prev,
      obj_now = step$objective,
      freeze_step = step$freeze_happened
    )
    if (verbose) {
      cat(sprintf(
        "[soft-partition conv %3d] T=%.3f obj=%.6f delta=%.3e rel_delta=%.3e streak=%d/%d\n",
        iter2, T_end, step$objective, delta,
        convergence_info$last_rel_delta,
        convergence_info$consecutive_small_steps,
        convergence_info$consecutive_required
      ))
    }
    if (is.finite(obj_prev) && delta < -1e-8) {
      warning(sprintf("soft_partition_cavi objective decreased by %.3e at phase-2 iter %d.",
                      delta, iter2), call. = FALSE)
    }
    if (isTRUE(convergence_info$converged)) {
      converged <- TRUE
      break
    }
    obj_prev <- step$objective
  }
  if (!isTRUE(converged) && convergence_info$tol_outer > 0) {
    convergence_info$reason <- sprintf(
      "Stopped after %d phase-2 iterations without meeting the %d-step convergence rule (last rel_delta=%.3e).",
      convergence_info$phase2_iters,
      convergence_info$consecutive_required,
      convergence_info$last_rel_delta
    )
  }

  # Trim histories
  score_history <- score_history[seq_len(actual_iters)]
  weight_history <- weight_history[seq_len(actual_iters)]
  effective_weight_history <- effective_weight_history[seq_len(actual_iters)]
  objective_history <- objective_history[seq_len(actual_iters)]

  # Hard assignment if requested
  if (hard_assign_final) {
    hard <- matrix(0, nrow = d, ncol = M)
    hard[cbind(seq_len(d), max.col(pi_weights, ties.method = "first"))] <- 1
    pi_weights <- hard
    colnames(pi_weights) <- ord_labels
  }

  assignment_info_final <- .cavi_partition_assignment_info(
    weights = pi_weights,
    T_now = T_end,
    partition_prior = assignment_ctl$partition_prior,
    partition_prior_init = assignment_ctl$partition_prior_init,
    assignment_prior = assignment_ctl$legacy_assignment_prior,
    ordering_alpha = assignment_ctl$ordering_alpha,
    assignment_M = M,
    active_orderings = active_orderings,
    active_feature_pairs = active_feature_pairs,
    drop_unused_ordering = inactive_ctl$drop_unused_ordering
  )
  assignment_posterior <- .cavi_partition_assignment_state(assignment_info_final)

  out <- structure(
    list(
      data = X,
      measurement_sd = partition_measurement_sd,
      pi_weights = pi_weights,
      effective_pi_weights = effective_pi_weights,
      assign = ord_labels[max.col(pi_weights, ties.method = "first")],
      fits = fits,
      score_history = score_history,
      weight_history = weight_history,
      effective_weight_history = effective_weight_history,
      objective_history = objective_history,
      init_info = init_info,
      ordering_similarity = ordering_similarity,
      similarity_init = similarity_init,
      T_schedule = T_schedule,
      n_anneal = as.integer(n_outer),
      converged = converged,
      convergence_info = convergence_info,
      M = M,
      active_orderings = active_orderings,
      active_feature_pairs = active_feature_pairs,
      frozen_orderings = !active_orderings,
      ordering_events = ordering_events,
      feature_events = feature_events,
      assignment_posterior = assignment_posterior,
      control = list(
        noise_model = partition_noise_info$noise_model,
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        position_prior = position_prior_ctl$position_prior,
        partition_prior = assignment_ctl$partition_prior,
        partition_prior_init = assignment_ctl$partition_prior_init,
        assignment_prior = assignment_info_final$assignment_prior,
        assignment_mode = assignment_ctl$assignment_mode,
        ordering_alpha = assignment_ctl$ordering_alpha,
        tol_outer = tol_outer,
        partition_init = applied_partition_init,
        similarity_metric = similarity_metric,
        smooth_fit_lambda_mode = smooth_fit_lambda_mode,
        smooth_fit_lambda_value = smooth_fit_lambda_value,
        cluster_linkage = cluster_linkage,
        similarity_min_feature_sd = similarity_min_feature_sd,
        sigma_min = sigma_min,
        sigma_max = sigma_max,
        freeze_unused_ordering = inactive_ctl$freeze_unused_ordering,
        freeze_unused_ordering_threshold = inactive_ctl$freeze_unused_ordering_threshold,
        freeze_feature = feature_ctl$freeze_feature,
        freeze_feature_weight_threshold = feature_ctl$freeze_feature_weight_threshold,
        drop_unused_ordering = inactive_ctl$drop_unused_ordering
      )
    ),
    class = "soft_partition_cavi"
  )
  out$priors <- .mpcurve_priors_from_partition_fit(out)
  out
}


#' Soft dual-trajectory partition using explicit mean-field CAVI
#'
#' @description
#' Wrapper around \code{\link{soft_partition_cavi}} with \code{M = 2}. Provides
#' the original dual-trajectory interface for backward compatibility.
#'
#' @param X Numeric matrix (\code{n x d}).
#' @param fit1_init,fit2_init Optional initial \code{cavi} fits. If NULL,
#'   they are initialized via \code{init_two_trajectories_cavi()}.
#' @param init_method,init_method1 Initialization methods.
#' @param init_seed Random seed for \code{"random_split"}.
#' @param K Number of pseudotime bins.
#' @param discretization Initial discretization scheme used when converting
#'   ordering scores into \code{K} bins. Dual-ordering initialisation enforces a
#'   common \code{K} across the two orderings.
#' @param T_start,T_end Annealing temperatures.
#' @param n_outer Number of annealing steps.
#' @param inner_iter Number of weighted CAVI sweeps per outer step.
#' @param max_converge_iter Maximum number of exact \code{T=1} iterations.
#' @param tol_outer Relative objective tolerance for the stricter \code{T=1}
#'   convergence rule used after annealing.
#' @param rw_q Random-walk order.
#' @param ridge Optional nugget added to the RW precision for both orderings.
#' @param lambda_sd_prior_rate Optional positive rate for an exponential prior
#'   on \code{1 / sqrt(lambda_j^{(m)})}. The default \code{NULL} means no lambda
#'   prior penalty. For backward compatibility, an explicit \code{0} is treated
#'   the same way; it is only an alias for "no penalty" and does not
#'   correspond to a literal exponential prior with rate zero.
#' @param lambda_min,lambda_max Bounds passed through to the weighted
#'   single-ordering \code{cavi} updates.
#' @param sigma_min,sigma_max Bounds for the feature-specific noise variance
#'   \code{sigma2_j}. Forwarded to internal \code{cavi()} calls.
#' @param position_prior Either \code{"adaptive"} or \code{"fixed"} for the
#'   within-ordering component prior \code{pi}.
#' @param position_prior_init Optional length-\code{K} initial/fixed vector for
#'   the within-ordering component prior \code{pi}. When
#'   \code{position_prior = "fixed"} and this is \code{NULL}, a uniform prior
#'   over the \code{K} positions is used.
#' @param partition_prior Either \code{"adaptive"} or \code{"fixed"} for the
#'   ordering-usage prior \code{omega}.
#' @param partition_prior_init Optional length-2 initial/fixed vector for the
#'   ordering prior \code{omega}. Used only when
#'   \code{partition_prior = "fixed"}.
#' @param assignment_prior Deprecated compatibility argument for the old
#'   partition-prior API.
#' @param ordering_alpha Deprecated compatibility argument used only for the
#'   legacy \code{assignment_prior = "dirichlet"} path.
#' @param hard_assign_final Logical.
#' @param freeze_unused_ordering Logical; if \code{TRUE}, unused orderings may
#'   be frozen and skipped in subsequent updates.
#' @param freeze_unused_ordering_threshold Non-negative feature-mass threshold
#'   used when \code{freeze_unused_ordering = TRUE}.
#' @param freeze_feature Logical; if \code{TRUE}, feature-ordering pairs with
#'   sufficiently small posterior weight are frozen and their trajectory
#'   contribution is neutralized.
#' @param freeze_feature_weight_threshold Non-negative threshold on
#'   \code{w_{jm}} used when \code{freeze_feature = TRUE}.
#' @param drop_unused_ordering Logical; if \code{TRUE}, user-facing wrappers may
#'   compact the returned partition view to the active orderings after fitting.
#'   With \code{FALSE}, the exposed \code{$objective_history} keeps the
#'   fixed-requested-\code{M} comparison semantics of
#'   \code{\link{soft_partition_cavi}}; with \code{TRUE}, it becomes the
#'   post-drop fitting objective of the active model.
#' @param verbose Logical.
#'
#' @return An object of class \code{"soft_partition_cavi"} with fields:
#'   \code{$fits} (list of 2 cavi fits), \code{$pi_weights} (d x 2),
#'   \code{$assign}, \code{$M}, \code{$objective_history}, etc. The
#'   \code{$objective_history} semantics match
#'   \code{\link{soft_partition_cavi}}.
#' @noRd
soft_two_trajectory_cavi <- function(X,
                                     S = NULL,
                                     fit1_init = NULL,
                                     fit2_init = NULL,
                                     init_method = c("score", "mincor", "random_split"),
                                     init_method1 = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
                                     init_seed = 42L,
                                     K = NULL,
                                     discretization = c("quantile", "equal", "kmeans"),
                                     T_start = 5,
                                     T_end = 1,
                                     n_outer = 25L,
                                     inner_iter = 1L,
                                     max_converge_iter = 100L,
                                     tol_outer = 1e-5,
                                     rw_q = 2L,
                                     ridge = 0,
                                     lambda_sd_prior_rate = NULL,
                                     lambda_min = 1e-10,
                                     lambda_max = 1e10,
                                     sigma_min = 1e-10,
                                     sigma_max = 1e10,
                                     position_prior = c("adaptive", "fixed"),
                                     position_prior_init = NULL,
                                     partition_prior = c("adaptive", "fixed"),
                                     partition_prior_init = NULL,
                                     assignment_prior = NULL,
                                     ordering_alpha = NULL,
                                     hard_assign_final = FALSE,
                                     freeze_unused_ordering = TRUE,
                                     freeze_unused_ordering_threshold = 0.5,
                                     freeze_feature = TRUE,
                                     freeze_feature_weight_threshold = 0.1,
                                     drop_unused_ordering = FALSE,
                                     verbose = TRUE) {
  init_method <- match.arg(init_method)
  init_method1 <- match.arg(init_method1)
  discretization <- match.arg(discretization)
  position_prior <- match.arg(position_prior)
  partition_prior_missing <- missing(partition_prior)
  partition_prior <- match.arg(partition_prior)
  if (isTRUE(partition_prior_missing) && !is.null(assignment_prior)) {
    assignment_prior_chr <- as.character(assignment_prior)[1]
    if (assignment_prior_chr %in% c("uniform", "dirichlet")) {
      partition_prior <- if (identical(assignment_prior_chr, "uniform")) "fixed" else "adaptive"
    }
  }
  X <- as.matrix(X)

  # Build fits_init from old-style parameters
  fits_init <- NULL
  if (!is.null(fit1_init) && !is.null(fit2_init)) {
    fits_init <- list(fit1_init, fit2_init)
  } else if (!is.null(fit1_init) || !is.null(fit2_init)) {
    # One provided, one not — use old init path
    inits <- init_two_trajectories_cavi(
      X = X, S = S, method = init_method, method1 = init_method1,
      K = K, rw_q = rw_q, ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min, lambda_max = lambda_max,
      sigma_min = sigma_min, sigma_max = sigma_max,
      discretization = discretization,
      num_iter = 5L, seed = init_seed, verbose = verbose
    )
    fits_init <- list(
      fit1_init %||% inits$fit1,
      fit2_init %||% inits$fit2
    )
  } else {
    # Both NULL — use old init strategy then pass to soft_partition_cavi
    inits <- init_two_trajectories_cavi(
      X = X, S = S, method = init_method, method1 = init_method1,
      K = K, rw_q = rw_q, ridge = ridge,
      lambda_sd_prior_rate = lambda_sd_prior_rate,
      lambda_min = lambda_min, lambda_max = lambda_max,
      sigma_min = sigma_min, sigma_max = sigma_max,
      discretization = discretization,
      num_iter = 5L, seed = init_seed, verbose = verbose
    )
    fits_init <- list(inits$fit1, inits$fit2)
  }

  soft_partition_cavi(
    X = X, S = S, M = 2L,
    fits_init = fits_init,
    K = K,
    discretization = discretization,
    T_start = T_start, T_end = T_end,
    n_outer = n_outer, inner_iter = inner_iter,
    max_converge_iter = max_converge_iter, tol_outer = tol_outer,
    rw_q = rw_q, ridge = ridge,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min, lambda_max = lambda_max,
    sigma_min = sigma_min, sigma_max = sigma_max,
    position_prior = position_prior,
    position_prior_init = position_prior_init,
    partition_prior = partition_prior,
    partition_prior_init = partition_prior_init,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    hard_assign_final = hard_assign_final,
    freeze_unused_ordering = freeze_unused_ordering,
    freeze_unused_ordering_threshold = freeze_unused_ordering_threshold,
    freeze_feature = freeze_feature,
    freeze_feature_weight_threshold = freeze_feature_weight_threshold,
    drop_unused_ordering = drop_unused_ordering,
    verbose = verbose
  )
}
