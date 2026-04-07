# ============================================================
# mpcurve S3 class - unified model object for MPCurve fits
# ============================================================
# MPCurve is the statistical model (GMM + GMRF prior on mu);
# cavi, csmooth_em, and smooth_em are the currently supported fitting
# backends. The recommended/default path is now cavi.
#
# as_mpcurve()   - convert a raw algorithm fit to an mpcurve object
# fit_mpcurve()  - run the full pipeline and return an mpcurve
#
# mpcurve object fields:
#   $params$pi      K-vector   mixture weights
#   $params$mu      d x K      cluster means (unified matrix format)
#   $params$sigma2  d-vec or d x K  per-feature variances
#   $gamma          n x K      responsibility matrix
#   $locations      inferred cell locations from gamma
#   $elbo_trace     numeric    penalized ELBO over iterations
#   $loglik_trace   numeric    observed log-likelihood over iterations
#   $lambda_trace   list       lambda_vec per iteration
#   $iter           integer    number of EM iterations run
#   $K, $n, $d      integers   dimensions
#   $algorithm      character  "cavi", "smooth_em", or "csmooth_em"
#   $modelName      character  covariance structure string
#   $intrinsic_dim  integer    1 = single ordering, 2 = dual-ordering partition
#   $fit            raw        underlying cavi / smooth_em / csmooth_em object
#
# For intrinsic_dim = 2 (soft partition), additional fields:
#   $fits           list(fit1, fit2)  mpcurve objects for each ordering
#   $partition      list with:
#     $pi_weights   d x 2 matrix of feature-to-ordering weights
#     $assign       d-length character vector ("A" / "B")
#   $objective_history  numeric  variational objective trace
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

.mpcurve_ordering_labels <- function(M) {
  if (M <= 26L) LETTERS[seq_len(M)] else paste0("ord", seq_len(M))
}

.mpcurve_extract_algorithm_arg <- function(dots, caller = "fit_mpcurve()") {
  if (!("algorithm" %in% names(dots))) {
    return(list(dots = dots, warned = FALSE))
  }

  alg <- as.character(dots$algorithm %||% "cavi")[1]
  dots$algorithm <- NULL
  if (!identical(alg, "cavi")) {
    stop(
      caller,
      " is now a CAVI-only public wrapper. For legacy fits, call cavi(), ",
      "do_cavi(), do_csmoothEM(), or do_smoothEM() directly.",
      call. = FALSE
    )
  }

  warning(
    "`algorithm` is no longer part of the public ", caller,
    " interface; MPCurver now always uses the CAVI backend here.",
    call. = FALSE
  )
  list(dots = dots, warned = TRUE)
}

.mpcurve_check_legacy_public_args <- function(dots,
                                              caller = "fit_mpcurve()") {
  legacy_args <- c(
    "relative_lambda",
    "adaptive",
    "sigma_update",
    "check_decrease",
    "tol_decrease"
  )
  bad <- intersect(names(dots), legacy_args)
  if (!length(bad)) return(invisible(NULL))

  stop(
    caller,
    " is now a CAVI-only public wrapper. Legacy smoothEM/csmoothEM controls ",
    "are no longer supported here: ",
    paste(sprintf("`%s`", bad), collapse = ", "),
    ". Use the lower-level legacy functions directly if you need those paths.",
    call. = FALSE
  )
}

.mpcurve_update_named_list <- function(x, updates) {
  if (!length(updates)) return(x)
  for (nm in names(updates)) x[[nm]] <- updates[[nm]]
  x
}

.mpcurve_last_objective <- function(x) {
  obj <- if (inherits(x, "mpcurve")) {
    if (!is.null(x$partition) || inherits(x$fit, "soft_partition_cavi")) {
      x$objective_history %||% numeric(0)
    } else {
      x$elbo_trace %||% numeric(0)
    }
  } else if (inherits(x, "soft_partition_cavi")) {
    x$objective_history %||% numeric(0)
  } else if (inherits(x, "cavi")) {
    x$elbo_trace %||% numeric(0)
  } else {
    numeric(0)
  }
  if (!length(obj)) -Inf else as.numeric(tail(obj, 1L))
}

.mpcurve_greedy_method_for_dim <- function(method,
                                           target_dim,
                                           method_missing = FALSE) {
  target_dim <- as.integer(target_dim)[1]
  if (target_dim <= 1L) {
    if (isTRUE(method_missing)) return("PCA")
    return(method[[1L]])
  }
  if (isTRUE(method_missing)) return("PCA")
  if (length(method) == 1L) return(method)
  if (length(method) >= target_dim) return(method[seq_len(target_dim)])
  method
}

.mpcurve_fit_once <- function(fit_args,
                              dots = list(),
                              fit_overrides = list(),
                              dot_overrides = list()) {
  fit_args_use <- .mpcurve_update_named_list(fit_args, fit_overrides)
  dots_use <- .mpcurve_update_named_list(dots, dot_overrides)
  do.call(
    fit_mpcurve,
    c(fit_args_use, dots_use, list(greedy = "none"))
  )
}

.mpcurve_partition_active_idx <- function(raw_fit) {
  M <- raw_fit$M %||% length(raw_fit$fits)
  which(raw_fit$active_orderings %||% rep(TRUE, M))
}

.mpcurve_partition_ordering_contributions <- function(raw_fit) {
  if (!inherits(raw_fit, "soft_partition_cavi")) {
    stop("raw_fit must inherit from 'soft_partition_cavi'.")
  }
  X <- raw_fit$fits[[1]]$data
  ctl <- raw_fit$control %||% list()
  assignment_prior <- ctl$assignment_prior %||% "uniform"
  ordering_alpha <- ctl$ordering_alpha %||% 0.5
  obj_terms <- .cavi_partition_objective_from_fits(
    fits = raw_fit$fits,
    X = X,
    weights = raw_fit$pi_weights,
    T_now = 1,
    active_orderings = raw_fit$active_orderings,
    active_feature_pairs = raw_fit$active_feature_pairs,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    drop_unused_ordering = FALSE
  )
  contrib <- colSums(raw_fit$pi_weights * obj_terms$like_mat) +
    colSums(obj_terms$prior_entropy_mat) +
    obj_terms$cell_terms_by_fit
  ord_labels <- colnames(raw_fit$pi_weights) %||%
    .mpcurve_ordering_labels(length(contrib))
  stats::setNames(as.numeric(contrib), ord_labels)
}

.mpcurve_partition_gamma_correlation <- function(raw_fit) {
  active_idx <- .mpcurve_partition_active_idx(raw_fit)
  if (length(active_idx) < 2L) {
    return(list(
      active_idx = active_idx,
      corr_mat = matrix(NA_real_, nrow = length(active_idx), ncol = length(active_idx))
    ))
  }
  gamma_vecs <- lapply(raw_fit$fits[active_idx], function(fit) as.numeric(fit$gamma))
  corr_mat <- matrix(NA_real_, nrow = length(active_idx), ncol = length(active_idx))
  diag(corr_mat) <- 1
  for (i in seq_along(active_idx)) {
    for (j in seq_len(i - 1L)) {
      cij <- suppressWarnings(stats::cor(gamma_vecs[[i]], gamma_vecs[[j]]))
      if (!is.finite(cij)) cij <- 0
      corr_mat[i, j] <- abs(cij)
      corr_mat[j, i] <- abs(cij)
    }
  }
  list(active_idx = active_idx, corr_mat = corr_mat)
}

.mpcurve_partition_backward_drop <- function(raw_fit) {
  corr_info <- .mpcurve_partition_gamma_correlation(raw_fit)
  active_idx <- corr_info$active_idx
  if (length(active_idx) < 2L) return(NULL)

  corr_mat <- corr_info$corr_mat
  diag(corr_mat) <- -Inf
  max_corr <- max(corr_mat, na.rm = TRUE)
  if (!is.finite(max_corr)) return(NULL)

  pair_rows <- which(corr_mat == max_corr, arr.ind = TRUE)
  pair_rows <- pair_rows[pair_rows[, 1] < pair_rows[, 2], , drop = FALSE]
  if (!nrow(pair_rows)) return(NULL)

  contrib <- .mpcurve_partition_ordering_contributions(raw_fit)
  ord_labels <- names(contrib)
  pair_tbl <- data.frame(
    i = integer(0),
    j = integer(0),
    drop_idx = integer(0),
    drop_contrib = numeric(0),
    stringsAsFactors = FALSE
  )
  for (r in seq_len(nrow(pair_rows))) {
    i_loc <- pair_rows[r, 1]
    j_loc <- pair_rows[r, 2]
    i_idx <- active_idx[i_loc]
    j_idx <- active_idx[j_loc]
    ci <- contrib[[ord_labels[i_idx]]]
    cj <- contrib[[ord_labels[j_idx]]]
    drop_idx <- if (ci <= cj) i_idx else j_idx
    drop_contrib <- min(ci, cj)
    pair_tbl <- rbind(
      pair_tbl,
      data.frame(
        i = i_idx,
        j = j_idx,
        drop_idx = drop_idx,
        drop_contrib = drop_contrib,
        stringsAsFactors = FALSE
      )
    )
  }
  best <- pair_tbl[order(pair_tbl$drop_contrib, pair_tbl$drop_idx), , drop = FALSE][1L, , drop = FALSE]
  keep_idx <- setdiff(active_idx, best$drop_idx)
  list(
    drop_idx = as.integer(best$drop_idx),
    drop_label = ord_labels[best$drop_idx],
    pair_idx = c(as.integer(best$i), as.integer(best$j)),
    pair_label = paste(ord_labels[c(best$i, best$j)], collapse = " vs "),
    corr = as.numeric(max_corr),
    keep_idx = keep_idx,
    contributions = contrib
  )
}

.mpcurve_compact_partition_fit <- function(current_fit,
                                           fit_args,
                                           dots,
                                           method_missing = FALSE) {
  raw_fit <- current_fit$fit
  active_idx <- .mpcurve_partition_active_idx(raw_fit)
  active_dim <- length(active_idx)
  current_requested <- current_fit$requested_intrinsic_dim %||% current_fit$intrinsic_dim

  if (!inherits(raw_fit, "soft_partition_cavi") ||
      active_dim == current_requested ||
      active_dim < 1L) {
    return(list(fit = current_fit, compacted = FALSE, active_dim = active_dim))
  }

  if (active_dim == 1L) {
    compact_fit <- as_mpcurve(raw_fit$fits[[active_idx]])
    return(list(fit = compact_fit, compacted = TRUE, active_dim = 1L))
  }

  compact_fit <- .mpcurve_fit_once(
    fit_args = fit_args,
    dots = dots,
    fit_overrides = list(
      intrinsic_dim = active_dim,
      method = .mpcurve_greedy_method_for_dim(
        fit_args$method,
        active_dim,
        method_missing = method_missing
      ),
      drop_unused_ordering = FALSE
    ),
    dot_overrides = list(
      fits_init = raw_fit$fits[active_idx]
    )
  )
  list(fit = compact_fit, compacted = TRUE, active_dim = active_dim)
}

.mpcurve_finalize_greedy_return <- function(selected_fit,
                                            fit_args,
                                            dots,
                                            method_missing = FALSE,
                                            greedy_info = NULL) {
  want_drop <- isTRUE(fit_args$drop_unused_ordering)
  final_fit <- selected_fit

  if (want_drop && inherits(selected_fit$fit, "soft_partition_cavi")) {
    selected_requested <- selected_fit$requested_intrinsic_dim %||% selected_fit$intrinsic_dim
    final_fit <- .mpcurve_fit_once(
      fit_args = fit_args,
      dots = dots,
      fit_overrides = list(
        intrinsic_dim = selected_requested,
        method = .mpcurve_greedy_method_for_dim(
          fit_args$method,
          selected_requested,
          method_missing = method_missing
        ),
        drop_unused_ordering = TRUE
      ),
      dot_overrides = list(
        fits_init = selected_fit$fit$fits
      )
    )
  }

  final_fit$greedy_selection <- greedy_info
  final_fit
}

.mpcurve_greedy_forward <- function(fit_args, dots, method_missing = FALSE) {
  upper_bound <- as.integer(fit_args$intrinsic_dim)[1]
  history <- data.frame(
    step = integer(0),
    direction = character(0),
    current_M = integer(0),
    candidate_M = integer(0),
    current_objective = numeric(0),
    candidate_objective = numeric(0),
    accepted = logical(0),
    dropped_label = character(0),
    pair_label = character(0),
    note = character(0),
    stringsAsFactors = FALSE
  )

  current_fit <- .mpcurve_fit_once(
    fit_args = fit_args,
    dots = dots,
    fit_overrides = list(
      intrinsic_dim = 1L,
      method = .mpcurve_greedy_method_for_dim(
        fit_args$method,
        1L,
        method_missing = method_missing
      ),
      num_cores = 1L,
      drop_unused_ordering = FALSE
    )
  )
  current_obj <- .mpcurve_last_objective(current_fit)
  current_M <- 1L
  stop_reason <- if (upper_bound <= 1L) "reached_upper_bound" else "smaller_model_preferred"

  if (upper_bound >= 2L) {
    for (step in seq_len(upper_bound - 1L)) {
      candidate_M <- current_M + 1L
      candidate_fit <- .mpcurve_fit_once(
        fit_args = fit_args,
        dots = dots,
        fit_overrides = list(
          intrinsic_dim = candidate_M,
          method = .mpcurve_greedy_method_for_dim(
            fit_args$method,
            candidate_M,
            method_missing = method_missing
          ),
          drop_unused_ordering = FALSE
        )
      )
      candidate_obj <- .mpcurve_last_objective(candidate_fit)
      accepted <- is.finite(candidate_obj) && candidate_obj > current_obj
      history <- rbind(
        history,
        data.frame(
          step = step,
          direction = "forward",
          current_M = current_M,
          candidate_M = candidate_M,
          current_objective = current_obj,
          candidate_objective = candidate_obj,
          accepted = accepted,
          dropped_label = NA_character_,
          pair_label = NA_character_,
          note = if (accepted) "larger_model_preferred" else "smaller_model_preferred",
          stringsAsFactors = FALSE
        )
      )
      if (!accepted) {
        stop_reason <- "smaller_model_preferred"
        break
      }
      current_fit <- candidate_fit
      current_obj <- candidate_obj
      current_M <- candidate_M
      stop_reason <- if (current_M >= upper_bound) "reached_upper_bound" else "smaller_model_preferred"
    }
  }

  greedy_info <- list(
    mode = "forward",
    requested_upper_bound = upper_bound,
    selected_intrinsic_dim = current_M,
    stop_reason = stop_reason,
    history = history
  )
  .mpcurve_finalize_greedy_return(
    selected_fit = current_fit,
    fit_args = fit_args,
    dots = dots,
    method_missing = method_missing,
    greedy_info = greedy_info
  )
}

.mpcurve_greedy_backward <- function(fit_args, dots, method_missing = FALSE) {
  upper_bound <- as.integer(fit_args$intrinsic_dim)[1]
  if (upper_bound <= 1L) {
    fit1 <- .mpcurve_fit_once(
      fit_args = fit_args,
      dots = dots,
      fit_overrides = list(
        intrinsic_dim = 1L,
        method = .mpcurve_greedy_method_for_dim(
          fit_args$method,
          1L,
          method_missing = method_missing
        ),
        num_cores = 1L,
        drop_unused_ordering = FALSE
      )
    )
    greedy_info <- list(
      mode = "backward",
      requested_upper_bound = upper_bound,
      selected_intrinsic_dim = 1L,
      starting_active_intrinsic_dim = 1L,
      stop_reason = "reached_dimension_1",
      history = data.frame(
        step = integer(0),
        direction = character(0),
        current_M = integer(0),
        candidate_M = integer(0),
        current_objective = numeric(0),
        candidate_objective = numeric(0),
        accepted = logical(0),
        dropped_label = character(0),
        pair_label = character(0),
        note = character(0),
        stringsAsFactors = FALSE
      )
    )
    return(.mpcurve_finalize_greedy_return(
      selected_fit = fit1,
      fit_args = fit_args,
      dots = dots,
      method_missing = method_missing,
      greedy_info = greedy_info
    ))
  }

  history <- data.frame(
    step = integer(0),
    direction = character(0),
    current_M = integer(0),
    candidate_M = integer(0),
    current_objective = numeric(0),
    candidate_objective = numeric(0),
    accepted = logical(0),
    dropped_label = character(0),
    pair_label = character(0),
    note = character(0),
    stringsAsFactors = FALSE
  )

  current_fit <- .mpcurve_fit_once(
    fit_args = fit_args,
    dots = dots,
    fit_overrides = list(
      intrinsic_dim = upper_bound,
      method = .mpcurve_greedy_method_for_dim(
        fit_args$method,
        upper_bound,
        method_missing = method_missing
      ),
      drop_unused_ordering = FALSE
    )
  )
  start_active_dim <- if (inherits(current_fit$fit, "soft_partition_cavi")) {
    length(.mpcurve_partition_active_idx(current_fit$fit))
  } else {
    1L
  }
  compact_start <- .mpcurve_compact_partition_fit(
    current_fit = current_fit,
    fit_args = fit_args,
    dots = dots,
    method_missing = method_missing
  )
  current_fit <- compact_start$fit
  current_obj <- .mpcurve_last_objective(current_fit)
  current_M <- current_fit$requested_intrinsic_dim %||% current_fit$intrinsic_dim
  stop_reason <- if (current_M <= 1L) "reached_dimension_1" else "larger_model_preferred"

  step_id <- 0L
  while (current_M > 1L && inherits(current_fit$fit, "soft_partition_cavi")) {
    compact_now <- .mpcurve_compact_partition_fit(
      current_fit = current_fit,
      fit_args = fit_args,
      dots = dots,
      method_missing = method_missing
    )
    current_fit <- compact_now$fit
    current_obj <- .mpcurve_last_objective(current_fit)
    current_M <- current_fit$requested_intrinsic_dim %||% current_fit$intrinsic_dim
    if (!inherits(current_fit$fit, "soft_partition_cavi") || current_M <= 1L) {
      stop_reason <- "reached_dimension_1"
      break
    }

    drop_info <- .mpcurve_partition_backward_drop(current_fit$fit)
    if (is.null(drop_info)) {
      stop_reason <- "reached_dimension_1"
      break
    }

    candidate_M <- current_M - 1L
    if (candidate_M <= 1L) {
      candidate_fit <- .mpcurve_fit_once(
        fit_args = fit_args,
        dots = dots,
        fit_overrides = list(
          intrinsic_dim = 1L,
          method = .mpcurve_greedy_method_for_dim(
            fit_args$method,
            1L,
            method_missing = method_missing
          ),
          num_cores = 1L,
          drop_unused_ordering = FALSE
        ),
        dot_overrides = list(
          responsibilities_init = current_fit$fit$fits[[drop_info$keep_idx]]$gamma
        )
      )
    } else {
      candidate_fit <- .mpcurve_fit_once(
        fit_args = fit_args,
        dots = dots,
        fit_overrides = list(
          intrinsic_dim = candidate_M,
          method = .mpcurve_greedy_method_for_dim(
            fit_args$method,
            candidate_M,
            method_missing = method_missing
          ),
          drop_unused_ordering = FALSE
        ),
        dot_overrides = list(
          fits_init = current_fit$fit$fits[drop_info$keep_idx]
        )
      )
    }
    candidate_obj <- .mpcurve_last_objective(candidate_fit)
    accepted <- is.finite(candidate_obj) && candidate_obj > current_obj
    step_id <- step_id + 1L
    history <- rbind(
      history,
      data.frame(
        step = step_id,
        direction = "backward",
        current_M = current_M,
        candidate_M = candidate_M,
        current_objective = current_obj,
        candidate_objective = candidate_obj,
        accepted = accepted,
        dropped_label = drop_info$drop_label,
        pair_label = drop_info$pair_label,
        note = if (accepted) "smaller_model_preferred" else "larger_model_preferred",
        stringsAsFactors = FALSE
      )
    )
    if (!accepted) {
      stop_reason <- "larger_model_preferred"
      break
    }
    current_fit <- candidate_fit
    current_obj <- candidate_obj
    current_M <- current_fit$requested_intrinsic_dim %||% current_fit$intrinsic_dim
    stop_reason <- if (current_M <= 1L) "reached_dimension_1" else "larger_model_preferred"
  }

  greedy_info <- list(
    mode = "backward",
    requested_upper_bound = upper_bound,
    selected_intrinsic_dim = current_M,
    starting_active_intrinsic_dim = start_active_dim,
    stop_reason = stop_reason,
    history = history
  )
  .mpcurve_finalize_greedy_return(
    selected_fit = current_fit,
    fit_args = fit_args,
    dots = dots,
    method_missing = method_missing,
    greedy_info = greedy_info
  )
}

.mpcurve_locations_from_gamma <- function(gamma) {
  if (is.null(gamma)) return(NULL)

  gamma <- as.matrix(gamma)
  if (!nrow(gamma) || !ncol(gamma)) return(NULL)

  K <- ncol(gamma)
  comp_grid <- seq_len(K)
  pseudo_grid <- if (K <= 1L) rep(0, K) else (comp_grid - 1L) / (K - 1L)

  map_index <- max.col(gamma, ties.method = "first")
  mean_index <- as.numeric(gamma %*% comp_grid)

  list(
    mean = list(
      index = mean_index,
      pseudotime = as.numeric(gamma %*% pseudo_grid)
    ),
    map = list(
      index = map_index,
      pseudotime = pseudo_grid[map_index]
    ),
    component_grid = comp_grid,
    pseudotime_grid = pseudo_grid
  )
}


# ---- as_mpcurve generic ----------------------------------------

#' Convert an algorithm fit object to an \code{mpcurve} model object
#'
#' @param x A \code{cavi}, \code{smooth_em}, or \code{csmooth_em} object.
#' @param ... Ignored.
#' @return An object of class \code{"mpcurve"}.
#'
#' For single-ordering fits, the returned object stores a unified \code{$params}
#' block, the responsibility matrix \code{$gamma}, and a \code{$locations} field
#' containing inferred cell locations derived from \code{$gamma}:
#' \itemize{
#'   \item \code{$locations$mean$index}: posterior-mean component index in \code{[1, K]}
#'   \item \code{$locations$mean$pseudotime}: responsibility-weighted pseudotime in \code{[0, 1]}
#'   \item \code{$locations$map$index}: MAP component index from \code{max.col(gamma)}
#'   \item \code{$locations$map$pseudotime}: pseudotime corresponding to the MAP component
#' }
#'
#' For partition fits (\code{intrinsic_dim >= 2}), \code{$locations} is a named
#' list with one such location object per ordering.
#' @export
as_mpcurve <- function(x, ...) UseMethod("as_mpcurve")


#' @export
as_mpcurve.csmooth_em <- function(x, ...) {
  params <- x$params

  # Unified mu: d x K matrix
  mu_mat <- do.call(cbind, params$mu)

  K <- length(params$pi)
  n <- if (!is.null(x$data)) nrow(x$data) else
    if (!is.null(x$gamma)) nrow(x$gamma) else NA_integer_
  d <- nrow(mu_mat)

  modelName <- x$control$modelName %||% "homoskedastic"

  structure(
    list(
      params = list(
        pi     = as.numeric(params$pi),
        mu     = mu_mat,
        sigma2 = params$sigma2   # d-vec (homo) or d x K matrix (hetero)
      ),
      gamma        = x$gamma,
      locations    = .mpcurve_locations_from_gamma(x$gamma),
      elbo_trace   = x$elbo_trace   %||% numeric(0),
      loglik_trace = x$loglik_trace %||% numeric(0),
      lambda_trace = x$lambda_trace %||% list(),
      iter         = as.integer(x$iter %||% length(x$elbo_trace %||% numeric(0))),
      K            = as.integer(K),
      n            = as.integer(n),
      d            = as.integer(d),
      algorithm      = "csmooth_em",
      modelName      = modelName,
      intrinsic_dim  = 1L,
      fit            = x
    ),
    class = "mpcurve"
  )
}


#' @export
as_mpcurve.cavi <- function(x, ...) {
  params <- x$params
  mu_mat <- x$posterior$mean

  K <- length(params$pi)
  n <- if (!is.null(x$data)) nrow(x$data) else
    if (!is.null(x$gamma)) nrow(x$gamma) else NA_integer_
  d <- nrow(mu_mat)

  structure(
    list(
      params = list(
        pi     = as.numeric(params$pi),
        mu     = mu_mat,
        sigma2 = params$sigma2
      ),
      gamma        = x$gamma,
      measurement_sd = x$measurement_sd %||% NULL,
      locations    = .mpcurve_locations_from_gamma(x$gamma),
      elbo_trace   = x$elbo_trace   %||% numeric(0),
      loglik_trace = x$loglik_trace %||% numeric(0),
      lambda_trace = x$lambda_trace %||% list(),
      iter         = as.integer(x$iter %||% length(x$elbo_trace %||% numeric(0))),
      K            = as.integer(K),
      n            = as.integer(n),
      d            = as.integer(d),
      algorithm      = "cavi",
      modelName      = x$control$modelName %||% "homoskedastic",
      intrinsic_dim  = 1L,
      fit            = x
    ),
    class = "mpcurve"
  )
}


#' @export
as_mpcurve.smooth_em <- function(x, ...) {
  params <- x$params

  # Unified mu: d x K matrix
  mu_mat <- do.call(cbind, params$mu)

  K <- length(params$pi)
  n <- if (!is.null(x$data)) nrow(x$data) else
    if (!is.null(x$gamma)) nrow(x$gamma) else NA_integer_
  d <- nrow(mu_mat)

  # Extract sigma2: diagonal of each sigma[[k]] -> d x K matrix
  sigma2_mat <- if (!is.null(params$sigma) && length(params$sigma) == K) {
    vapply(params$sigma, function(S) diag(as.matrix(S)), numeric(d))
  } else {
    matrix(NA_real_, d, K)
  }

  modelName <- x$control$modelName %||% (x$modelName %||% "unknown")

  structure(
    list(
      params = list(
        pi     = as.numeric(params$pi),
        mu     = mu_mat,
        sigma2 = sigma2_mat   # d x K  (one column per cluster)
      ),
      gamma        = x$gamma,
      locations    = .mpcurve_locations_from_gamma(x$gamma),
      elbo_trace   = x$elbo_trace   %||% numeric(0),
      loglik_trace = x$loglik_trace %||% numeric(0),
      lambda_trace = x$lambda_trace %||% list(),
      iter         = as.integer(x$iter %||% length(x$elbo_trace %||% numeric(0))),
      K            = as.integer(K),
      n            = as.integer(n),
      d            = as.integer(d),
      algorithm      = "smooth_em",
      modelName      = modelName,
      intrinsic_dim  = 1L,
      fit            = x
    ),
    class = "mpcurve"
  )
}


#' @export
as_mpcurve.soft_partition_cavi <- function(x, ...) {
  M <- x$M %||% 2L
  ord_labels_full <- colnames(x$pi_weights) %||% .mpcurve_ordering_labels(M)
  active_full <- x$active_orderings %||% rep(TRUE, M)
  frozen_full <- x$frozen_orderings %||% !active_full
  drop_view <- isTRUE((x$control %||% list())$drop_unused_ordering)
  keep_idx <- if (drop_view) which(active_full) else seq_len(M)
  if (!length(keep_idx)) {
    keep_idx <- which.max(ifelse(active_full, 1, 0))
  }
  ord_labels <- ord_labels_full[keep_idx]

  # Support unified $fits list
  raw_fits <- x$fits
  if (is.null(raw_fits)) {
    # Legacy fallback (should not happen with new code)
    raw_fits <- list(x$fit1, x$fit2)
  }

  fits_mp <- lapply(raw_fits[keep_idx], as_mpcurve)
  ref <- fits_mp[[1]]
  visible_weights <- x$pi_weights[, keep_idx, drop = FALSE]
  colnames(visible_weights) <- ord_labels
  visible_active <- if (drop_view) rep(TRUE, length(keep_idx)) else active_full[keep_idx]
  visible_frozen <- if (drop_view) rep(FALSE, length(keep_idx)) else frozen_full[keep_idx]

  structure(
    list(
      fits      = fits_mp,
      measurement_sd = ref$measurement_sd %||% NULL,
      locations = stats::setNames(lapply(fits_mp, `[[`, "locations"), ord_labels),
      partition = list(
        pi_weights = visible_weights,
        assign     = x$assign,
        active_orderings = visible_active,
        frozen_orderings = visible_frozen,
        frozen_labels = ord_labels_full[frozen_full],
        ordering_events = x$ordering_events %||% list(),
        assignment_posterior = x$assignment_posterior,
        requested_intrinsic_dim = as.integer(M),
        displayed_intrinsic_dim = as.integer(length(keep_idx)),
        dropped_labels = if (drop_view) ord_labels_full[frozen_full] else character(0),
        kept_labels = ord_labels,
        compacted = drop_view
      ),
      objective_history = x$objective_history %||% numeric(0),
      K             = ref$K,
      n             = ref$n,
      d             = ref$d,
      algorithm     = "cavi",
      intrinsic_dim = as.integer(length(keep_idx)),
      requested_intrinsic_dim = as.integer(M),
      init_info     = x$init_info,
      ordering_similarity = x$ordering_similarity,
      similarity_init = x$similarity_init %||% NULL,
      converged     = x$converged,
      convergence_info = x$convergence_info,
      ordering_events = x$ordering_events %||% list(),
      control       = x$control %||% list(),
      fit           = x
    ),
    class = "mpcurve"
  )
}


# ---- S3 methods ------------------------------------------------

#' @export
print.mpcurve <- function(x, ...) {
  idim <- x$intrinsic_dim %||% 1L
  is_partition <- !is.null(x$partition) || inherits(x$fit, "soft_partition_cavi")

  if (is_partition) {
    cat(sprintf("MPCurve model fit (%d-ordering partition)\n", idim))
    cat(sprintf("  Algorithm      : %s\n", x$algorithm))
    greedy_info <- x$greedy_selection %||% NULL
    if (!is.null(greedy_info)) {
      cat(sprintf(
        "  Greedy search  : %s (upper bound %d -> selected %d)\n",
        greedy_info$mode,
        greedy_info$requested_upper_bound,
        greedy_info$selected_intrinsic_dim
      ))
    }
    requested_idim <- x$requested_intrinsic_dim %||% idim
    if (!identical(as.integer(requested_idim), as.integer(idim))) {
      cat(sprintf("  Intrinsic dim  : %d (requested %d)\n", idim, requested_idim))
    } else {
      cat(sprintf("  Intrinsic dim  : %d\n", idim))
    }
    cat(sprintf("  n / d / K      : %d / %d / %d\n", x$n, x$d, x$K))
    part <- x$partition
    if (!is.null(part)) {
      tbl <- table(part$assign)
      part_str <- paste(sprintf("%s=%d", names(tbl), as.integer(tbl)), collapse = ", ")
      cat(sprintf("  Partition      : %s\n", part_str))
      active <- part$active_orderings
      dropped <- part$dropped_labels
      ord_labels <- colnames(part$pi_weights) %||% .mpcurve_ordering_labels(idim)
      if (!is.null(active)) {
        cat(sprintf("  Active         : %s\n",
                    paste(ord_labels[active], collapse = ", ")))
      }
      frozen <- part$frozen_labels
      if (!is.null(frozen) && length(frozen) > 0L) {
        if (isTRUE(part$compacted)) {
          cat(sprintf("  Dropped from view: %s\n",
                      paste(frozen, collapse = ", ")))
        } else {
          cat(sprintf("  Frozen         : %s\n",
                      paste(frozen, collapse = ", ")))
        }
      }
      events <- part$ordering_events %||% x$ordering_events %||% list()
      if (length(events) > 0L) {
        event_types <- vapply(events, `[[`, character(1), "event")
        cat(sprintf("  Events         : %d freeze\n",
                    sum(event_types == "freeze")))
      }
    }
    if (length(x$objective_history) > 0L)
      cat(sprintf("  Objective (last): %.6f\n",
                  tail(x$objective_history, 1L)))
    drop_view <- isTRUE((x$control %||% list())$drop_unused_ordering)
    cat(sprintf("  Objective type : %s\n",
                if (drop_view) {
                  "post-drop fit (not comparable across requested dimensions)"
                } else {
                  "fixed-requested-M comparison"
                }))
    if (!is.null(x$converged)) {
      cat(sprintf("  Converged      : %s\n", x$converged))
    }
    conv_info <- x$convergence_info %||% x$fit$convergence_info
    if (!is.null(conv_info)) {
      if (!is.null(conv_info$last_rel_delta) && is.finite(conv_info$last_rel_delta)) {
        cat(sprintf("  Phase-2 rel_delta   : %.3e\n", conv_info$last_rel_delta))
      }
      if (!is.null(conv_info$consecutive_small_steps)) {
        cat(sprintf("  Small-step streak: %d\n", conv_info$consecutive_small_steps))
      }
      if (!is.null(conv_info$reason) && nzchar(conv_info$reason)) {
        cat(sprintf("  Reason         : %s\n", conv_info$reason))
      }
    }
  } else {
    cat("MPCurve model fit\n")
    cat(sprintf("  Algorithm : %s\n",  x$algorithm))
    greedy_info <- x$greedy_selection %||% NULL
    if (!is.null(greedy_info)) {
      cat(sprintf(
        "  Greedy    : %s (upper bound %d -> selected %d)\n",
        greedy_info$mode,
        greedy_info$requested_upper_bound,
        greedy_info$selected_intrinsic_dim
      ))
    }
    cat(sprintf("  Model     : %s\n",  x$modelName))
    cat(sprintf("  n / d / K : %d / %d / %d\n", x$n, x$d, x$K))
    cat(sprintf("  Iterations: %d\n",  x$iter))
    if (length(x$elbo_trace) > 0L)
      cat(sprintf("  ELBO (last)  : %.6f\n", tail(x$elbo_trace, 1L)))
    ml_trace <- x$fit$ml_trace
    if (!is.null(ml_trace) && length(ml_trace) > 0L)
      cat(sprintf("  ML   (last)  : %.6f\n", tail(ml_trace, 1L)))
    if (length(x$loglik_trace) > 0L)
      cat(sprintf("  logLik (last): %.6f\n", tail(x$loglik_trace, 1L)))
  }
  invisible(x)
}


#' Summarise an \code{mpcurve} model fit
#'
#' Delegates to the underlying \code{summary.cavi},
#' \code{summary.csmooth_em}, or \code{summary.smooth_em}, wrapping the result in a
#' \code{summary.mpcurve} object so users can access all fields
#' (\code{$ml_last}, \code{$elbo_last}, convergence diagnostics, etc.)
#' that the underlying class exposes.
#'
#' @param object An \code{mpcurve} object.
#' @param ... Passed to the underlying summary method.
#' @return An object of class \code{summary.mpcurve} with fields:
#'   \code{$algorithm}, \code{$modelName}, \code{$K}, \code{$n},
#'   \code{$d}, and \code{$underlying} (the full underlying summary).
#' @export
summary.mpcurve <- function(object, ...) {
  idim <- object$intrinsic_dim %||% 1L
  is_partition <- !is.null(object$partition) || inherits(object$fit, "soft_partition_cavi")

  if (is_partition) {
    result <- list(
      algorithm     = object$algorithm,
      intrinsic_dim = idim,
      requested_intrinsic_dim = object$requested_intrinsic_dim %||% idim,
      K             = object$K,
      n             = object$n,
      d             = object$d,
      partition     = object$partition,
      objective_history = object$objective_history,
      converged     = object$converged,
      convergence_info = object$convergence_info %||% object$fit$convergence_info,
      ordering_events = object$ordering_events %||% list(),
      control       = object$control %||% list(),
      greedy_selection = object$greedy_selection %||% NULL
    )
  } else {
    underlying <- summary(object$fit, ...)
    result <- list(
      algorithm     = object$algorithm,
      modelName     = object$modelName,
      intrinsic_dim = idim,
      K             = object$K,
      n             = object$n,
      d             = object$d,
      underlying    = underlying,
      greedy_selection = object$greedy_selection %||% NULL
    )
  }
  class(result) <- "summary.mpcurve"
  result
}


#' @export
print.summary.mpcurve <- function(x, ...) {
  idim <- x$intrinsic_dim %||% 1L
  is_partition <- !is.null(x$partition)

  if (is_partition) {
    cat("MPCurve Partition Summary\n")
    greedy_info <- x$greedy_selection %||% NULL
    if (!is.null(greedy_info)) {
      cat(sprintf("Greedy search : %s  |  upper bound = %d  |  selected = %d\n",
                  greedy_info$mode,
                  greedy_info$requested_upper_bound,
                  greedy_info$selected_intrinsic_dim))
    }
    requested_idim <- x$requested_intrinsic_dim %||% idim
    if (!identical(as.integer(requested_idim), as.integer(idim))) {
      cat(sprintf("Algorithm : %s  |  Intrinsic dim : %d (requested %d)  |  n=%d  d=%d  K=%d\n",
                  x$algorithm, idim, requested_idim, x$n, x$d, x$K))
    } else {
      cat(sprintf("Algorithm : %s  |  Intrinsic dim : %d  |  n=%d  d=%d  K=%d\n",
                  x$algorithm, idim, x$n, x$d, x$K))
    }
    if (!is.null(x$partition)) {
      tbl <- table(x$partition$assign)
      part_str <- paste(sprintf("%s=%d", names(tbl), as.integer(tbl)), collapse = "  ")
      cat(sprintf("Partition : %s\n", part_str))
      active <- x$partition$active_orderings
      ord_labels <- colnames(x$partition$pi_weights) %||% .mpcurve_ordering_labels(idim)
      if (!is.null(active)) {
        cat(sprintf("Active    : %s\n", paste(ord_labels[active], collapse = ", ")))
      }
      frozen <- x$partition$frozen_labels
      if (!is.null(frozen) && length(frozen) > 0L) {
        if (isTRUE(x$partition$compacted)) {
          cat(sprintf("Dropped from view : %s\n", paste(frozen, collapse = ", ")))
        } else {
          cat(sprintf("Frozen    : %s\n", paste(frozen, collapse = ", ")))
        }
      }
    }
    events <- x$ordering_events %||% list()
    if (length(events) > 0L) {
      event_types <- vapply(events, `[[`, character(1), "event")
      cat(sprintf("Events    : %d freeze\n",
                  sum(event_types == "freeze")))
    }
    if (length(x$objective_history) > 0L)
      cat(sprintf("Objective (last) : %.6f\n", tail(x$objective_history, 1L)))
    drop_view <- isTRUE((x$control %||% list())$drop_unused_ordering)
    cat(sprintf("Objective type : %s\n",
                if (drop_view) {
                  "post-drop fit (not comparable across requested dimensions)"
                } else {
                  "fixed-requested-M comparison"
                }))
    if (!is.null(x$converged))
      cat(sprintf("Converged : %s\n", x$converged))
    conv_info <- x$convergence_info
    if (!is.null(conv_info)) {
      if (!is.null(conv_info$last_rel_delta) && is.finite(conv_info$last_rel_delta)) {
        cat(sprintf("Phase-2 rel_delta : %.3e\n", conv_info$last_rel_delta))
      }
      if (!is.null(conv_info$consecutive_small_steps)) {
        cat(sprintf("Small-step streak : %d\n", conv_info$consecutive_small_steps))
      }
      if (!is.null(conv_info$reason) && nzchar(conv_info$reason)) {
        cat(sprintf("Reason : %s\n", conv_info$reason))
      }
    }
  } else {
    cat("MPCurve Model Summary\n")
    greedy_info <- x$greedy_selection %||% NULL
    if (!is.null(greedy_info)) {
      cat(sprintf("Greedy search : %s  |  upper bound = %d  |  selected = %d\n",
                  greedy_info$mode,
                  greedy_info$requested_upper_bound,
                  greedy_info$selected_intrinsic_dim))
    }
    cat(sprintf("Algorithm : %s  |  Model : %s  |  n=%d  d=%d  K=%d\n",
                x$algorithm, x$modelName, x$n, x$d, x$K))
    cat("Underlying fit summary\n")
    print(x$underlying, ...)
  }
  invisible(x)
}


# ---- internal gradient colour-bar legend ------------------------------------

# Draws a vertical gradient bar in the top-right corner of the current plot.
# Must be called after the main plot so par("usr") reflects real axis limits.
.draw_gradient_legend <- function(
    pal,
    title    = "pseudotime",
    lo_label = "0",
    hi_label = "1",
    n_rect   = 200L
) {
  usr <- par("usr")
  pw  <- usr[2] - usr[1]
  ph  <- usr[4] - usr[3]

  # bar: 3 % wide, 35 % tall, top-right with a 2 % margin
  bar_w  <- 0.03 * pw
  bar_h  <- 0.35 * ph
  margin <- 0.02

  xr <- usr[2] - margin * pw
  xl <- xr - bar_w
  yt <- usr[4] - margin * ph
  yb <- yt - bar_h

  # draw gradient (bottom = pseudotime 0, top = pseudotime 1)
  ys   <- seq(yb, yt, length.out = n_rect + 1L)
  cols <- pal[pmax(1L, round(seq(1L, length(pal), length.out = n_rect)))]

  old_xpd <- par(xpd = FALSE)
  on.exit(par(old_xpd), add = TRUE)

  for (i in seq_len(n_rect))
    rect(xl, ys[i], xr, ys[i + 1L], col = cols[i], border = NA)
  rect(xl, yb, xr, yt, border = "black", lwd = 0.5)

  # endpoint labels and title
  gap <- 0.008 * pw
  text(xr + gap, yb, lo_label, adj = c(0, 0.5), cex = 0.65)
  text(xr + gap, yt, hi_label, adj = c(0, 0.5), cex = 0.65)
  text((xl + xr) / 2, yt + 0.018 * ph, title, adj = c(0.5, 0), cex = 0.7)
}


#' Plot an \code{mpcurve} model fit
#'
#' @description
#' Visualise an \code{mpcurve} fit. The default (\code{plot_type="scatterplot"})
#' colors each observation by its estimated pseudotime position along the
#' inferred trajectory, mapped to a \eqn{[0, 1]} scale where component \eqn{k}
#' occupies position \eqn{(k-1)/(K-1)}.  The pseudotime for observation \eqn{i}
#' is the responsibility-weighted average
#' \deqn{t_i = \sum_{k=1}^{K} \gamma_{ik} \cdot \frac{k-1}{K-1}.}
#' Component means are overlaid as orange stars connected by arrows.
#' A vertical gradient color bar is drawn in the top-right corner of the plot.
#'
#' \describe{
#'   \item{\code{"scatterplot"}}{Scatter of two chosen dimensions, colored by
#'     pseudotime (default).}
#'   \item{\code{"elbo"}}{ELBO / log-likelihood / collapsed-ML traces.
#'     Delegates to the underlying \code{plot.cavi}, \code{plot.csmooth_em},
#'     or \code{plot.smooth_em}.}
#'   \item{\code{"mu"}}{Component-means-only plot.
#'     Delegates to the underlying fit's plot method.}
#' }
#'
#' @param x An \code{mpcurve} object.
#' @param plot_type One of \code{"scatterplot"} (default), \code{"elbo"},
#'   \code{"mu"}.
#' @param dims Integer vector of length 2 (or 1). Columns of the data matrix to
#'   use for the scatter plot.  Defaults to \code{c(1, 2)}.
#' @param data Optional numeric matrix (n x d).  If \code{NULL} (default), uses
#'   \code{x$fit$data}.  Must be supplied if the data were not stored at fit
#'   time.
#' @param pal Colour palette (length-256 character vector) used to map
#'   pseudotime to point colors.  Defaults to the classic pseudotime rainbow
#'   (navy \eqn{\to} cyan \eqn{\to} green \eqn{\to} gold \eqn{\to} red),
#'   which avoids white and is easy to read against a white background.
#' @param add_legend Logical; draw a gradient color-bar legend?  Default
#'   \code{TRUE}.
#' @param ... Further arguments forwarded to \code{plot_EM_embedding2D} (for
#'   \code{"scatterplot"}) or the underlying algorithm's \code{plot} method
#'   (for \code{"elbo"} / \code{"mu"}).
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.mpcurve <- function(
    x,
    plot_type  = c("scatterplot", "elbo", "mu"),
    dims       = c(1L, 2L),
    data       = NULL,
    pal        = grDevices::colorRampPalette(
                   c("#0000FF", "#00FFFF", "#00FF00",
                     "#FFFF00", "#FF0000"))(256L),
    add_legend = TRUE,
    ...
) {
  if (!inherits(x, "mpcurve")) stop("x must be an 'mpcurve' object.")
  plot_type <- match.arg(plot_type)
  idim <- x$intrinsic_dim %||% 1L

  # ---- Multi-ordering partition plot ----
  if (idim >= 2L && plot_type == "scatterplot") {
    M <- idim
    fits_list <- x$fits
    pi_w <- x$partition$pi_weights   # d x M

    # Resolve data
    if (is.null(data)) {
      data <- fits_list[[1]]$fit$data
      if (is.null(data))
        stop("No data found. Please supply `data` explicitly.")
    }
    data <- as.matrix(data)

    dims <- as.integer(dims)
    if (length(dims) < 1L || length(dims) > 2L) stop("dims must have length 1 or 2.")
    if (any(dims < 1L) || any(dims > ncol(data)))
      stop("dims out of range for the number of columns in data.")

    ord_labels <- colnames(pi_w)
    if (is.null(ord_labels)) {
      ord_labels <- if (M <= 26L) LETTERS[seq_len(M)] else paste0("ord", seq_len(M))
    }

    old_par <- par(mfrow = c(1, M), mar = c(4, 4, 3.5, 1))
    on.exit(par(old_par), add = TRUE)

    for (m in seq_len(M)) {
      sub_fit <- fits_list[[m]]
      K_m <- sub_fit$K
      gamma_m <- sub_fit$gamma
      positions_m <- (seq_len(K_m) - 1L) / (K_m - 1L)
      t_pseudo_m <- as.numeric(gamma_m %*% positions_m)

      # Colour by pseudotime
      rng <- range(t_pseudo_m, na.rm = TRUE)
      z <- (t_pseudo_m - rng[1]) / (rng[2] - rng[1] + 1e-12)
      idx <- pmax(1L, pmin(length(pal), 1L + floor(z * (length(pal) - 1L))))
      pt_col <- pal[idx]

      # Weight annotation for plotted dims
      w_dims <- pi_w[dims, m]
      w_str <- paste(sprintf("d%d:w=%.2f", dims, w_dims), collapse = ", ")
      main_m <- sprintf("Ordering %s  (%s)", ord_labels[m], w_str)

      mu_m <- sub_fit$params$mu   # d x K

      if (length(dims) == 2L) {
        plot(data[, dims[1]], data[, dims[2]],
             pch = 19, col = pt_col, cex = 0.6,
             main = main_m,
             xlab = sprintf("dim %d", dims[1]),
             ylab = sprintf("dim %d", dims[2]),
             ...)
        # Overlay component means with arrows
        for (k in 2:K_m) {
          arrows(mu_m[dims[1], k - 1], mu_m[dims[2], k - 1],
                 mu_m[dims[1], k], mu_m[dims[2], k],
                 col = "orange", lwd = 1.5, length = 0.06)
        }
        points(mu_m[dims[1], ], mu_m[dims[2], ],
               pch = 8, col = "orange", cex = 0.9)
      } else {
        j <- dims[1]
        plot(t_pseudo_m, data[, j],
             pch = 19, col = pt_col, cex = 0.6,
             main = main_m,
             xlab = "pseudotime",
             ylab = sprintf("dim %d", j),
             ...)
        mu_j <- mu_m[j, ]
        lines(positions_m, mu_j, col = "orange", lwd = 2)
        points(positions_m, mu_j, pch = 8, col = "orange", cex = 0.9)
      }

      if (add_legend)
        .draw_gradient_legend(pal, title = "pseudotime",
                              lo_label = "0", hi_label = "1")
    }

    return(invisible(x))
  }

  # ---- delegate non-scatterplot types to the underlying fit ----
  if (plot_type != "scatterplot") {
    if (idim >= 2L) {
      # For partition objects, plot objective history
      obj <- x$objective_history
      if (length(obj) > 0L) {
        plot(obj, type = "b", pch = 19, cex = 0.7,
             xlab = "Iteration", ylab = "Variational objective",
             main = "Partition objective trace")
      }
    } else {
      plot(x$fit, plot_type = plot_type, ...)
    }
    return(invisible(x))
  }

  # ---- resolve data ----
  if (is.null(data)) {
    data <- x$fit$data
    if (is.null(data))
      stop("No data found in x$fit$data. Please supply `data` explicitly.")
  }
  data <- as.matrix(data)

  dims <- as.integer(dims)
  if (length(dims) < 1L || length(dims) > 2L) stop("dims must have length 1 or 2.")
  if (any(dims < 1L) || any(dims > ncol(data)))
    stop("dims out of range for the number of columns in data.")

  # ---- pseudotime: responsibility-weighted component position in [0, 1] ----
  K     <- x$K
  gamma <- x$gamma   # n x K
  if (is.null(gamma)) stop("x$gamma is NULL; cannot compute pseudotime.")

  positions <- (seq_len(K) - 1L) / (K - 1L)   # 0, 1/(K-1), ..., 1
  t_pseudo  <- as.numeric(gamma %*% positions)  # length n

  # ---- component means projected onto dims ----
  mu_mat <- x$params$mu   # d x K

  # ---- title ----
  main <- sprintf("MPCurve pseudotime  (algorithm: %s, K=%d)",
                  x$algorithm, K)

  # ---- map pseudotime -> color ----
  rng <- range(t_pseudo, na.rm = TRUE)
  z   <- (t_pseudo - rng[1]) / (rng[2] - rng[1] + 1e-12)
  idx <- pmax(1L, pmin(length(pal), 1L + floor(z * (length(pal) - 1L))))
  pt_col <- pal[idx]

  if (length(dims) == 2L) {
    # ===== 2D scatter =====
    mu_list_dims <- lapply(seq_len(K), function(k) mu_mat[dims, k])

    plot_EM_embedding2D(
      mu_list    = mu_list_dims,
      X2         = data[, dims, drop = FALSE],
      t_vec      = t_pseudo,
      pal        = pal,
      add_legend = FALSE,
      main       = main,
      xlab       = sprintf("dim %d", dims[1]),
      ylab       = sprintf("dim %d", dims[2]),
      ...
    )

  } else {
    # ===== 1D: pseudotime (x) vs selected dimension (y) =====
    j    <- dims[1]
    xlab <- "pseudotime"
    ylab <- sprintf("dim %d", j)

    mu_j <- mu_mat[j, ]               # length K
    mu_t <- positions                  # pseudotime of each component (0..1)

    graphics::plot(t_pseudo, data[, j],
                   pch = 19, col = pt_col, cex = 0.7,
                   xlab = xlab, ylab = ylab, main = main, ...)
    graphics::lines(mu_t, mu_j, col = "orange", lwd = 2)
    graphics::points(mu_t, mu_j, pch = 8, col = "orange", cex = 1)
  }

  if (add_legend)
    .draw_gradient_legend(pal, title = "pseudotime", lo_label = "0", hi_label = "1")

  invisible(x)
}


# ---- do_mpcurve ------------------------------------------------

#' Continue CAVI iterations on an existing \code{mpcurve} fit
#'
#' Public continuation wrapper for MPCurve's recommended CAVI paths. For
#' single-ordering fits, this delegates to \code{\link{do_cavi}}. For partition
#' fits, it continues the exact \eqn{T = 1} phase of
#' \code{\link{soft_partition_cavi}} using the stored frozen-ordering and
#' assignment-prior controls.
#'
#' Legacy smoothEM/csmoothEM continuation remains available through the lower
#' level \code{do_smoothEM()} and \code{do_csmoothEM()} functions, not through
#' this high-level wrapper.
#'
#' @param object An \code{mpcurve} object wrapping a \code{cavi} or
#'   \code{soft_partition_cavi} fit.
#' @param iter Integer >= 1. Maximum number of additional CAVI sweeps (single
#'   ordering) or exact \code{T = 1} partition steps.
#' @param lambda Optional scalar or d-vector. For single-ordering CAVI only:
#'   fix \code{lambda_j} at this value for the continuation run.
#' @param S Optional known measurement standard deviations. If \code{NULL},
#'   reuse the value stored on \code{object}. Supplying a new \code{S} is only
#'   allowed when it matches the stored specification.
#' @param tol Optional relative ELBO tolerance for single-ordering CAVI. If
#'   \code{NULL}, reuse the stored value.
#' @param lambda_sd_prior_rate Optional positive rate for the induced
#'   exponential prior on \code{1 / sqrt(lambda_j)}. If \code{NULL}, the stored
#'   control value is reused. An explicit \code{0} is treated as "no penalty"
#'   for backward compatibility and does not represent a literal exponential
#'   prior with rate zero.
#' @param lambda_min,lambda_max Optional positive bounds for \code{lambda_j}. If
#'   \code{NULL}, reuse the stored values.
#' @param sigma_min,sigma_max Optional positive bounds for \code{sigma_j^2}. If
#'   \code{NULL}, reuse the stored values.
#' @param tol_outer For partition fits only: relative objective tolerance used
#'   by the phase-2 convergence rule. If \code{NULL}, reuse the stored value.
#' @param freeze_unused_ordering For partition fits only: if \code{TRUE},
#'   inactive orderings remain frozen and skipped in subsequent weighted
#'   updates. If \code{NULL}, reuse the stored value.
#' @param freeze_unused_ordering_threshold For partition fits only: non-negative
#'   feature-mass threshold used to decide when an ordering is effectively
#'   unused. If \code{NULL}, reuse the stored value.
#' @param freeze_feature For partition fits only: if \code{TRUE},
#'   feature-ordering pairs with sufficiently small posterior weight are frozen
#'   and their trajectory contribution is neutralized. If \code{NULL}, reuse
#'   the stored value.
#' @param freeze_feature_weight_threshold For partition fits only: non-negative
#'   threshold on \code{w_{jm}} used when \code{freeze_feature = TRUE}. If
#'   \code{NULL}, reuse the stored value.
#' @param drop_unused_ordering For partition fits only: if \code{TRUE}, the
#'   returned \code{mpcurve} view may compact away frozen orderings. If
#'   \code{NULL}, reuse the stored value. When \code{FALSE}, the reported
#'   partition \code{$objective_history} retains fixed-requested-\code{M}
#'   comparison semantics; when \code{TRUE}, it becomes the post-drop fitting
#'   objective of the active model.
#' @param assignment_prior For partition fits only: either \code{"uniform"} or
#'   \code{"dirichlet"}. If \code{NULL}, reuse the stored value.
#' @param ordering_alpha For partition fits only: positive scalar concentration
#'   used when \code{assignment_prior = "dirichlet"}. If \code{NULL}, reuse the
#'   stored value.
#' @param verbose Logical. Print per-iteration progress?
#'
#' @return An updated \code{mpcurve} object.
#' @export
do_mpcurve <- function(object,
                       iter = 1,
                       lambda = NULL,
                       S = NULL,
                       tol = NULL,
                       lambda_sd_prior_rate = NULL,
                       lambda_min = NULL,
                       lambda_max = NULL,
                       sigma_min = NULL,
                       sigma_max = NULL,
                       tol_outer = NULL,
                       freeze_unused_ordering = NULL,
                       freeze_unused_ordering_threshold = NULL,
                       freeze_feature = NULL,
                       freeze_feature_weight_threshold = NULL,
                       drop_unused_ordering = NULL,
                       assignment_prior = NULL,
                       ordering_alpha = NULL,
                       verbose = FALSE) {
  if (!inherits(object, "mpcurve")) {
    stop("object must be an 'mpcurve' object.")
  }
  iter <- as.integer(iter)
  if (length(iter) != 1L || is.na(iter) || iter < 1L) {
    stop("iter must be a single integer >= 1.")
  }
  if (!identical(object$algorithm, "cavi")) {
    stop(
      "do_mpcurve() is now a CAVI-only public wrapper. ",
      "For legacy fits, call do_csmoothEM() or do_smoothEM() directly.",
      call. = FALSE
    )
  }

  if (inherits(object$fit, "soft_partition_cavi")) {
    sp <- object$fit
    M <- sp$M %||% 2L
    fits <- sp$fits
    if (!length(fits)) {
      stop("No ordering fits found in object$fit$fits.")
    }
    stored_S <- object$measurement_sd %||% sp$measurement_sd %||% (fits[[1]]$measurement_sd %||% NULL)
    if (!is.null(S)) {
      if (is.null(stored_S)) {
        stop("This fit was created without S; refit with S instead of adding it in do_mpcurve().",
             call. = FALSE)
      }
      if (!.cavi_same_measurement_sd(stored_S, S)) {
        stop("Supplied S does not match the measurement_sd stored on object.", call. = FALSE)
      }
    }

    X <- fits[[1]]$data
    if (is.null(X)) stop("No data found in fits[[1]]$data.")

    control <- sp$control %||% list()
    fit_control <- fits[[1]]$control %||% list()
    lmin <- lambda_min %||% fit_control$lambda_min %||% 1e-10
    lmax <- lambda_max %||% fit_control$lambda_max %||% 1e10
    tol_outer_use <- tol_outer %||% control$tol_outer %||%
      .partition_convergence_defaults()$tol_outer
    tol_outer_use <- as.numeric(tol_outer_use)[1]
    ord_labels <- .mpcurve_ordering_labels(M)
    active_orderings <- sp$active_orderings %||% rep(TRUE, M)
    active_feature_pairs <- sp$active_feature_pairs %||%
      matrix(rep(active_orderings, each = ncol(X)), nrow = ncol(X), ncol = M)

    if (is.null(freeze_unused_ordering)) {
      freeze_unused_ordering <- control$freeze_unused_ordering %||% TRUE
    }
    if (is.null(freeze_unused_ordering_threshold)) {
      freeze_unused_ordering_threshold <- control$freeze_unused_ordering_threshold %||% 0.5
    }
    if (is.null(drop_unused_ordering)) {
      drop_unused_ordering <- control$drop_unused_ordering %||% FALSE
    }
    if (is.null(freeze_feature)) {
      freeze_feature <- control$freeze_feature %||% TRUE
    }
    if (is.null(freeze_feature_weight_threshold)) {
      freeze_feature_weight_threshold <- control$freeze_feature_weight_threshold %||% 0.1
    }
    if (is.null(lambda_sd_prior_rate)) {
      lambda_sd_prior_rate <- control$lambda_sd_prior_rate %||% NULL
    } else {
      lambda_sd_prior_rate <- .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
    }
    if (is.null(assignment_prior)) {
      assignment_prior <- control$assignment_prior %||% "uniform"
    }
    if (is.null(ordering_alpha)) {
      ordering_alpha <- control$ordering_alpha %||% 0.5
    }

    fits <- lapply(fits, function(fit) {
      fit$control <- utils::modifyList(
        fit$control %||% list(),
        list(
          lambda_sd_prior_rate = lambda_sd_prior_rate,
          lambda_min = lmin,
          lambda_max = lmax
        )
      )
      fit
    })

    convergence_info <- sp$convergence_info
    if (is.null(convergence_info) ||
        !isTRUE(all.equal(convergence_info$tol_outer, tol_outer_use))) {
      convergence_info <- .partition_init_convergence_info(
        tol_outer = tol_outer_use,
        phase2_iters = sp$convergence_info$phase2_iters %||% 0L
      )
    }

    new_obj <- numeric(iter)
    new_scores <- vector("list", iter)
    new_weights <- vector("list", iter)
    new_effective_weights <- vector("list", iter)
    new_events <- list()
    new_feature_events <- list()
    assignment_posterior <- sp$assignment_posterior %||% NULL
    obj_prev <- if (length(sp$objective_history %||% numeric(0)) > 0L) {
      tail(sp$objective_history, 1L)
    } else {
      -Inf
    }
    n_run <- 0L

    for (i in seq_len(iter)) {
      step <- .soft_partition_step(
        fits, X, T_now = 1, inner_iter = 1L,
        lambda_min = lmin, lambda_max = lmax,
        active_orderings = active_orderings,
        active_feature_pairs = active_feature_pairs,
        weights_prev = sp$pi_weights,
        freeze_unused_ordering = freeze_unused_ordering,
        freeze_unused_ordering_threshold = freeze_unused_ordering_threshold,
        freeze_feature = freeze_feature,
        freeze_feature_weight_threshold = freeze_feature_weight_threshold,
        drop_unused_ordering = drop_unused_ordering,
        assignment_prior = assignment_prior,
        ordering_alpha = ordering_alpha,
        iter_index = length(sp$objective_history %||% numeric(0)) + i,
        ordering_labels = ord_labels
      )
      fits <- step$fits
      active_orderings <- step$active_orderings
      active_feature_pairs <- step$active_feature_pairs
      sp$pi_weights <- step$pi_weights
      sp$effective_pi_weights <- step$effective_pi_weights
      new_obj[i] <- step$objective
      new_scores[[i]] <- step$score_mat
      new_weights[[i]] <- step$pi_weights
      new_effective_weights[[i]] <- step$effective_pi_weights
      new_events <- c(new_events, step$ordering_events)
      new_feature_events <- c(new_feature_events, step$feature_events)
      assignment_posterior <- step$assignment_posterior
      n_run <- i

      delta <- step$objective - obj_prev
      convergence_info <- .partition_update_convergence_info(
        info = convergence_info,
        obj_prev = obj_prev,
        obj_now = step$objective,
        freeze_step = step$freeze_happened
      )
      if (verbose) {
        cat(sprintf(
          "[do_mpcurve conv %3d] obj=%.6f delta=%.3e rel_delta=%.3e streak=%d/%d\n",
          i,
          step$objective,
          delta,
          convergence_info$last_rel_delta,
          convergence_info$consecutive_small_steps,
          convergence_info$consecutive_required
        ))
      }
      if (is.finite(obj_prev) && delta < -1e-8) {
        warning(
          sprintf("soft_partition_cavi objective decreased by %.3e during do_mpcurve() step %d.",
                  delta, i),
          call. = FALSE
        )
      }
      if (isTRUE(convergence_info$converged)) {
        break
      }
      obj_prev <- step$objective
    }

    if (n_run < iter) {
      new_obj <- new_obj[seq_len(n_run)]
      new_scores <- new_scores[seq_len(n_run)]
      new_weights <- new_weights[seq_len(n_run)]
      new_effective_weights <- new_effective_weights[seq_len(n_run)]
    }
    if (!isTRUE(convergence_info$converged) && convergence_info$tol_outer > 0) {
      convergence_info$reason <- sprintf(
        "do_mpcurve() ran %d exact T=1 steps without meeting the %d-step convergence rule (last rel_delta=%.3e).",
        convergence_info$phase2_iters,
        convergence_info$consecutive_required,
        convergence_info$last_rel_delta
      )
    }

    sp$fits <- fits
    sp$pi_weights <- if (n_run > 0L) step$pi_weights else sp$pi_weights
    sp$effective_pi_weights <- if (n_run > 0L) step$effective_pi_weights else (sp$effective_pi_weights %||% sp$pi_weights)
    colnames(sp$pi_weights) <- ord_labels
    colnames(sp$effective_pi_weights) <- ord_labels
    sp$assign <- ord_labels[max.col(sp$pi_weights, ties.method = "first")]
    sp$objective_history <- c(sp$objective_history, new_obj)
    sp$score_history <- c(sp$score_history, new_scores)
    sp$weight_history <- c(sp$weight_history, new_weights)
    sp$effective_weight_history <- c(sp$effective_weight_history %||% list(), new_effective_weights)
    sp$active_orderings <- active_orderings
    sp$active_feature_pairs <- active_feature_pairs
    sp$frozen_orderings <- !active_orderings
    sp$ordering_events <- c(sp$ordering_events %||% list(), new_events)
    sp$feature_events <- c(sp$feature_events %||% list(), new_feature_events)
    sp$assignment_posterior <- assignment_posterior
    sp$converged <- isTRUE(convergence_info$converged)
    sp$convergence_info <- convergence_info
    sp$control <- utils::modifyList(
      control,
      list(
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        assignment_prior = as.character(assignment_prior)[1],
        ordering_alpha = as.numeric(ordering_alpha)[1],
        tol_outer = tol_outer_use,
        freeze_unused_ordering = isTRUE(freeze_unused_ordering),
        freeze_unused_ordering_threshold = as.numeric(freeze_unused_ordering_threshold),
        freeze_feature = isTRUE(freeze_feature),
        freeze_feature_weight_threshold = as.numeric(freeze_feature_weight_threshold),
        drop_unused_ordering = isTRUE(drop_unused_ordering)
      )
    )

    out <- as_mpcurve(sp)
    out$greedy_selection <- object$greedy_selection %||% NULL
    return(out)
  }

  if (!inherits(object$fit, "cavi")) {
    stop("object$fit must be a cavi or soft_partition_cavi object.")
  }

  tol_use <- tol %||% (object$fit$control %||% list())$tol %||% 1e-6
    new_fit <- do_cavi(
    object = object$fit,
    iter = iter,
    lambda = lambda,
    S = S,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    tol = tol_use,
    verbose = verbose
  )

  as_mpcurve(new_fit)
}


# ---- fit_mpcurve -----------------------------------------------

#' Fit an MPCurve model with the public CAVI interface
#'
#' High-level user-facing wrapper around MPCurve's recommended CAVI fitting
#' paths. Use \code{intrinsic_dim = 1} for a standard single-ordering fit and
#' \code{intrinsic_dim >= 2} for the partition-CAVI model with one ordering per
#' latent dimension.
#'
#' Legacy \code{smooth_em} and \code{csmooth_em} algorithms remain available as
#' lower-level compatibility functions, but they are no longer part of the
#' public \code{fit_mpcurve()} interface and are not the active development
#' path for the package.
#'
#' @param X Numeric matrix (n x d) of observations.
#' @param method Initialisation method(s) for the trajectory ordering. If
#'   \code{num_cores > 1} and \code{length(method) > 1} in the single-ordering
#'   case, all methods are run in parallel and a named list of fits is
#'   returned. For partition fits, \code{length(method)} must be either 1 or
#'   \code{intrinsic_dim}.
#' @param K Integer number of mixture components (grid knots).
#' @param rw_q Integer random-walk order for the GMRF prior.
#' @param lambda Positive scalar initial value for \code{lambda_j} in the
#'   single-ordering CAVI fit.
#' @param S Optional known measurement standard deviations. If a length-\code{d}
#'   vector, each feature uses a known shared standard deviation across
#'   observations. If an \code{n x d} matrix, each observation-feature pair uses
#'   its own known standard deviation.
#' @param fix_lambda Logical; if \code{TRUE}, keep \code{lambda_j} fixed at the
#'   supplied initial value in the single-ordering CAVI fit.
#' @param iter Maximum number of CAVI sweeps for the single-ordering fit.
#' @param tol Relative ELBO tolerance for the single-ordering CAVI fit.
#' @param num_cores Integer >= 1. Workers for parallel multi-method
#'   single-ordering runs.
#' @param intrinsic_dim Integer intrinsic dimensionality of the latent ordering
#'   system. \code{intrinsic_dim = 1} fits the standard single-ordering model.
#'   Values \code{>= 2} fit the partition-CAVI model.
#' @param greedy Dimension-selection mode. \code{"none"} keeps the requested
#'   \code{intrinsic_dim}. \code{"forward"} treats \code{intrinsic_dim} as an
#'   upper bound and compares \eqn{M} versus \eqn{M+1} sequentially starting
#'   from \eqn{M=1}. \code{"backward"} fits the upper bound first, then
#'   compares \eqn{M} versus \eqn{M-1} while greedily removing the most
#'   correlated active ordering.
#' @param partition_init For partition fits only: either
#'   \code{"similarity"} for feature-similarity-driven block initialization or
#'   \code{"ordering_methods"} for the existing per-ordering warm starts. The
#'   default is \code{"similarity"}.
#' @param discretization Optional discretization method passed to ordering-based
#'   initialization. For partition fits, MPCurver enforces a common \code{K}
#'   across orderings; if quantile cuts collapse, it falls back to equal-width
#'   bins.
#' @param ridge Optional nugget added to the RW precision.
#' @param lambda_sd_prior_rate Optional positive rate for an exponential prior
#'   on \code{1 / sqrt(lambda_j)}. Passed to the CAVI backend and the partition
#'   CAVI path. The default \code{NULL} means no lambda prior penalty. For
#'   backward compatibility, an explicit \code{0} is treated the same way; it
#'   is only an alias for "no penalty" and does not correspond to a literal
#'   exponential prior with rate zero.
#' @param lambda_min,lambda_max Positive bounds for \code{lambda_j}.
#' @param sigma_min,sigma_max Positive bounds for \code{sigma_j^2}. These apply
#'   to the single-ordering CAVI path and to the \code{"smooth_fit"}
#'   similarity metric when \code{partition_init = "similarity"}.
#' @param assignment_prior For partition fits only: either \code{"uniform"} or
#'   \code{"dirichlet"}.
#' @param ordering_alpha For partition fits only: positive scalar concentration
#'   used when \code{assignment_prior = "dirichlet"}.
#' @param similarity_metric For \code{partition_init = "similarity"} only:
#'   feature-similarity metric used to construct feature blocks. One of
#'   \code{"spearman"}, \code{"pearson"}, or \code{"smooth_fit"}. The
#'   \code{"smooth_fit"} metric can be used with \code{ridge = 0}, but for
#'   \code{intrinsic_dim > 1} this uses an intrinsic-RW pseudo-evidence rather
#'   than a fully proper marginal likelihood; use a small positive
#'   \code{ridge} if you want the smoother evidence to be theoretically proper.
#' @param smooth_fit_lambda_mode For \code{similarity_metric = "smooth_fit"}
#'   only: whether the directional smoother optimizes \code{lambda} or keeps it
#'   fixed at \code{smooth_fit_lambda_value}.
#' @param smooth_fit_lambda_value For \code{similarity_metric = "smooth_fit"}
#'   only: fixed \code{lambda} value used when
#'   \code{smooth_fit_lambda_mode = "fixed"}, and the starting value when
#'   \code{smooth_fit_lambda_mode = "optimize"}.
#' @param cluster_linkage For \code{partition_init = "similarity"} only:
#'   hierarchical-clustering linkage applied to \code{1 - S(X)}. Defaults to
#'   \code{"single"}.
#' @param similarity_min_feature_sd For \code{partition_init = "similarity"}
#'   only: low-variance feature threshold used when building \code{S(X)}.
#' @param T_start,T_end Annealing temperatures for partition CAVI.
#' @param n_outer Number of annealing steps for partition CAVI.
#' @param inner_iter Number of weighted CAVI sweeps per annealing step.
#' @param max_converge_iter Maximum number of exact \code{T = 1} partition
#'   iterations after annealing. Defaults to \code{iter} when \code{NULL}.
#' @param tol_outer Relative objective tolerance for partition phase-2
#'   convergence.
#' @param freeze_unused_ordering For partition fits only: if \code{TRUE},
#'   inactive orderings are frozen and skipped in subsequent weighted updates.
#' @param freeze_unused_ordering_threshold For partition fits only: non-negative
#'   feature-mass threshold used to decide when an ordering is effectively
#'   unused.
#' @param freeze_feature For partition fits only: if \code{TRUE},
#'   feature-ordering pairs with sufficiently small posterior weight are frozen
#'   and their trajectory contribution is neutralized.
#' @param freeze_feature_weight_threshold For partition fits only: non-negative
#'   threshold on \code{w_{jm}} used when \code{freeze_feature = TRUE}.
#' @param drop_unused_ordering For partition fits only: if \code{TRUE}, the
#'   returned \code{mpcurve} view may compact away frozen orderings. With
#'   \code{FALSE}, the reported partition \code{$objective_history} is the
#'   fixed-requested-\code{M} comparison objective; with \code{TRUE}, it is
#'   the post-drop fitting objective of the active model.
#' @param verbose Logical; print per-iteration progress?
#' @param ... Advanced CAVI / partition-CAVI options forwarded to the underlying
#'   backend. This is intended for advanced initialization controls such as
#'   \code{responsibilities_init}, \code{pi_init}, \code{sigma2_init},
#'   \code{fits_init}, \code{init_methods}, \code{pca_components}, or
#'   \code{hard_assign_final}. Legacy smoothEM/csmoothEM controls are rejected.
#'
#' @return An \code{\link{mpcurve}} object, or (for parallel multi-method
#'   single-ordering runs) a named list of \code{mpcurve} objects with a
#'   \code{summary} attribute. Every returned fit stores inferred cell
#'   locations in \code{$locations}; for single-ordering fits this includes both
#'   posterior-mean and MAP locations derived from \code{$gamma}. Greedy runs
#'   additionally attach a \code{$greedy_selection} metadata block recording
#'   the stepwise search history and stop reason.
#'
#' @export
fit_mpcurve <- function(
    X,
    method = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
    K = NULL,
    rw_q = 2,
    lambda = 1,
    S = NULL,
    fix_lambda = FALSE,
    iter = 100,
    tol = 1e-6,
    num_cores = 1L,
    intrinsic_dim = 1L,
    greedy = c("none", "forward", "backward"),
    partition_init = c("similarity", "ordering_methods"),
    discretization = NULL,
    ridge = 0,
    lambda_sd_prior_rate = NULL,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    sigma_min = 1e-10,
    sigma_max = 1e10,
    assignment_prior = c("uniform", "dirichlet"),
    ordering_alpha = 0.5,
    similarity_metric = c("spearman", "pearson", "smooth_fit"),
    smooth_fit_lambda_mode = c("optimize", "fixed"),
    smooth_fit_lambda_value = 1,
    cluster_linkage = "single",
    similarity_min_feature_sd = 1e-8,
    T_start = 5,
    T_end = 1,
    n_outer = 25L,
    inner_iter = 1L,
    max_converge_iter = NULL,
    tol_outer = 1e-5,
    freeze_unused_ordering = TRUE,
    freeze_unused_ordering_threshold = 0.5,
    freeze_feature = TRUE,
    freeze_feature_weight_threshold = 0.1,
    drop_unused_ordering = FALSE,
    verbose = FALSE,
    ...
) {
  method_missing <- missing(method)
  num_cores <- as.integer(num_cores)
  intrinsic_dim <- as.integer(intrinsic_dim)
  greedy <- match.arg(greedy)
  partition_init <- match.arg(partition_init)
  assignment_prior <- match.arg(assignment_prior)
  similarity_metric <- match.arg(similarity_metric)
  smooth_fit_lambda_mode <- match.arg(smooth_fit_lambda_mode)
  lambda_sd_prior_rate <- .normalize_lambda_sd_prior_rate(lambda_sd_prior_rate)
  dots <- list(...)
  algo_check <- .mpcurve_extract_algorithm_arg(dots, caller = "fit_mpcurve()")
  dots <- algo_check$dots
  .mpcurve_check_legacy_public_args(dots, caller = "fit_mpcurve()")
  max_converge_iter <- max_converge_iter %||% as.integer(iter)

  fit_args <- list(
    X = X,
    method = method,
    K = K,
    rw_q = rw_q,
    lambda = lambda,
    S = S,
    fix_lambda = fix_lambda,
    iter = iter,
    tol = tol,
    num_cores = num_cores,
    intrinsic_dim = intrinsic_dim,
    partition_init = partition_init,
    discretization = discretization,
    ridge = ridge,
    lambda_sd_prior_rate = lambda_sd_prior_rate,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    sigma_min = sigma_min,
    sigma_max = sigma_max,
    assignment_prior = assignment_prior,
    ordering_alpha = ordering_alpha,
    similarity_metric = similarity_metric,
    smooth_fit_lambda_mode = smooth_fit_lambda_mode,
    smooth_fit_lambda_value = smooth_fit_lambda_value,
    cluster_linkage = cluster_linkage,
    similarity_min_feature_sd = similarity_min_feature_sd,
    T_start = T_start,
    T_end = T_end,
    n_outer = n_outer,
    inner_iter = inner_iter,
    max_converge_iter = max_converge_iter,
    tol_outer = tol_outer,
    freeze_unused_ordering = freeze_unused_ordering,
    freeze_unused_ordering_threshold = freeze_unused_ordering_threshold,
    freeze_feature = freeze_feature,
    freeze_feature_weight_threshold = freeze_feature_weight_threshold,
    drop_unused_ordering = drop_unused_ordering,
    verbose = verbose
  )

  if (!identical(greedy, "none")) {
    return(
      switch(
        greedy,
        forward = .mpcurve_greedy_forward(
          fit_args = fit_args,
          dots = dots,
          method_missing = method_missing
        ),
        backward = .mpcurve_greedy_backward(
          fit_args = fit_args,
          dots = dots,
          method_missing = method_missing
        )
      )
    )
  }

  # ---- Partition model (intrinsic_dim >= 2) ----
  if (intrinsic_dim >= 2L) {
    M <- intrinsic_dim

    if (method_missing) {
      method <- "PCA"
    }

    if (identical(partition_init, "ordering_methods")) {
      if (length(method) == 1L) {
        init_methods <- c(method[1], rep("PCA", M - 1L))
        pca_counter <- 1L
        pca_components <- rep(NA_integer_, M)
        if (method[1] == "PCA") {
          pca_components[1] <- 1L
        }
        for (m in 2:M) {
          pca_counter <- pca_counter + 1L
          pca_components[m] <- pca_counter
        }
      } else if (length(method) == M) {
        init_methods <- method
        pca_components <- NULL
      } else {
        stop(sprintf("method must have length 1 or intrinsic_dim=%d.", M))
      }
    } else {
      if (!(length(method) %in% c(1L, M))) {
        stop(sprintf("method must have length 1 or intrinsic_dim=%d.", M))
      }
      init_methods <- method
      pca_components <- NULL
    }

    if ("init_methods" %in% names(dots)) {
      init_methods <- dots$init_methods
      dots$init_methods <- NULL
    }
    if ("pca_components" %in% names(dots)) {
      pca_components <- dots$pca_components
      dots$pca_components <- NULL
    }

    sp_args <- c(
      list(
        X = X,
        S = S,
        M = M,
        init_methods = init_methods,
        pca_components = pca_components,
        partition_init = partition_init,
        similarity_metric = similarity_metric,
        smooth_fit_lambda_mode = smooth_fit_lambda_mode,
        smooth_fit_lambda_value = smooth_fit_lambda_value,
        cluster_linkage = cluster_linkage,
        similarity_min_feature_sd = similarity_min_feature_sd,
        K = K,
        rw_q = rw_q,
        discretization = discretization %||% c("quantile", "equal", "kmeans")[1L],
        T_start = T_start,
        T_end = T_end,
        n_outer = n_outer,
        inner_iter = inner_iter,
        max_converge_iter = max_converge_iter,
        tol_outer = tol_outer,
        ridge = ridge,
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        sigma_min = sigma_min,
        sigma_max = sigma_max,
        assignment_prior = assignment_prior,
        ordering_alpha = ordering_alpha,
        freeze_unused_ordering = freeze_unused_ordering,
        freeze_unused_ordering_threshold = freeze_unused_ordering_threshold,
        freeze_feature = freeze_feature,
        freeze_feature_weight_threshold = freeze_feature_weight_threshold,
        drop_unused_ordering = drop_unused_ordering,
        verbose = verbose
      ),
      dots
    )
    raw <- do.call(soft_partition_cavi, sp_args)
    return(as_mpcurve(raw))
  }

  parallel <- num_cores > 1L && length(method) > 1L
  run_one_cavi <- function(method_i) {
    cavi_args <- c(
      list(
        X = X,
        K = K,
        method = method_i,
        rw_q = rw_q,
        ridge = ridge,
        discretization = discretization,
        S = S,
        fix_lambda = fix_lambda,
        lambda_sd_prior_rate = lambda_sd_prior_rate,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        sigma_min = sigma_min,
        sigma_max = sigma_max,
        max_iter = iter,
        tol = tol,
        verbose = verbose
      ),
      dots
    )
    if (!("lambda_init" %in% names(cavi_args))) {
      cavi_args$lambda_init <- lambda
    }
    tryCatch(
      do.call(cavi, cavi_args),
      error = function(e) NULL
    )
  }

  if (parallel) {
    raw_list <- if (.Platform$OS.type == "unix") {
      parallel::mclapply(method, run_one_cavi, mc.cores = num_cores)
    } else {
      lapply(method, run_one_cavi)
    }
    names(raw_list) <- method

    smry <- data.frame(
      method = method,
      K = vapply(raw_list, function(r) if (is.null(r)) NA_integer_ else length(r$params$pi), integer(1)),
      success = vapply(raw_list, function(r) !is.null(r), logical(1)),
      elbo_last = vapply(raw_list, function(r) if (is.null(r)) NA_real_ else tail(r$elbo_trace, 1L), numeric(1)),
      converged = vapply(raw_list, function(r) if (is.null(r)) FALSE else isTRUE(r$converged), logical(1)),
      stringsAsFactors = FALSE
    )

    out <- lapply(raw_list, function(r) if (!is.null(r)) as_mpcurve(r) else NULL)
    attr(out, "summary") <- smry
    return(out)
  }

  raw <- run_one_cavi(method[[1L]])
  if (is.null(raw)) {
    stop("fit_mpcurve() failed to produce a valid cavi fit.", call. = FALSE)
  }

  as_mpcurve(raw)
}
