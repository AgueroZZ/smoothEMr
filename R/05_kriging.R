#' Match locations on an old grid to indices on a new grid
#'
#' @param loc_old Numeric vector of locations to be matched (e.g. u_obs).
#' @param loc_new Numeric vector of candidate locations (e.g. u_final).
#' @param tol Nonnegative numeric tolerance for matching.
#'
#' @return Integer vector of indices into loc_new.
#' @export
match_locations_to_grid <- function(loc_old, loc_new, tol = 1e-8) {
  loc_old <- as.numeric(loc_old)
  loc_new <- as.numeric(loc_new)

  if (anyNA(loc_old) || anyNA(loc_new)) stop("loc_old/loc_new must not contain NA.")
  if (!is.numeric(tol) || length(tol) != 1L || tol < 0) stop("tol must be a single nonnegative number.")
  if (length(loc_new) < 1L) stop("loc_new must have length >= 1.")

  is_sorted <- all(diff(loc_new) >= 0)

  if (!is_sorted) {
    idx_obs <- vapply(loc_old, function(x) {
      j <- which.min(abs(loc_new - x))   # integer
      if (abs(loc_new[j] - x) > tol) {
        stop(sprintf("Cannot match location %.17g to new grid within tol=%.3e", x, tol))
      }
      as.integer(j)
    }, integer(1))
  } else {
    n_new <- length(loc_new)  # integer
    idx_obs <- vapply(loc_old, function(x) {
      j <- findInterval(x, loc_new)      # integer in 0..n_new
      cand <- c(j, j + 1L)               # keep integer
      cand <- cand[cand >= 1L & cand <= n_new]
      cand <- unique(as.integer(cand))

      dd <- abs(loc_new[cand] - x)
      jbest <- cand[which.min(dd)]

      if (abs(loc_new[jbest] - x) > tol) {
        stop(sprintf("Cannot match location %.17g to new grid within tol=%.3e", x, tol))
      }
      as.integer(jbest)
    }, integer(1))
  }

  if (any(duplicated(idx_obs))) stop("Duplicated indices in matching (loc_old not unique on loc_new?).")
  idx_obs
}


#' Construct a hierarchy of nested grids inside a final grid
#'
#' @description
#' Build nested grid levels indexed into a final grid.
#' Final grid has \code{K_final = 2^m_max + 1} equally-spaced points on [0,1].
#' Level \code{m} has \code{K_m = 2^m + 1} points, appearing as a subset of final indices.
#'
#' @param m_max Nonnegative integer; the finest level exponent.
#'
#' @return List with \code{m_max}, \code{K_final}, \code{u_final}, \code{idx_levels}, \code{K_levels}.
#' @export
make_hierarchical_levels <- function(m_max = 6) {
  if (!is.numeric(m_max) || length(m_max) != 1 || m_max < 0 || m_max != as.integer(m_max)) {
    stop("m_max must be a nonnegative integer.")
  }

  K_final <- 2^m_max + 1L
  u_final <- seq(0, 1, length.out = K_final)

  idx_levels <- lapply(0:m_max, function(m) {
    step <- 2^(m_max - m)
    seq.int(1L, K_final, by = step)
  })
  K_levels <- vapply(idx_levels, length, integer(1))

  list(
    m_max = m_max,
    K_final = K_final,
    u_final = u_final,
    idx_levels = idx_levels,
    K_levels = K_levels
  )
}


#' Scale random-walk penalty strength to account for grid spacing
#'
#' @description
#' Rescale \code{lambda} across nested equally-spaced grids on [0,1].
#'
#' Using the common heuristic for RW(q):
#' \deqn{\lambda(h) \propto h^{-(2q-1)}}
#' and \code{h = 1/(K-1)}, we use:
#' \deqn{\lambda_level = \lambda_final * ((K_level-1)/(K_final-1))^{(2q-1)}}.
#'
#' @param lambda_final Scalar lambda used on final grid.
#' @param K_final Number of points on the final grid.
#' @param K_level Number of points on the current level grid.
#' @param q RW order (e.g. 2 for RW2).
#'
#' @return Scalar \code{lambda_level}.
#' @export
lambda_scale_for_spacing <- function(lambda_final, K_final, K_level, q = 2) {
  if (!is.numeric(lambda_final) || length(lambda_final) != 1) stop("lambda_final must be a scalar.")
  if (!is.numeric(K_final) || length(K_final) != 1 || K_final < 2) stop("K_final must be >= 2.")
  if (!is.numeric(K_level) || length(K_level) != 1 || K_level < 2) stop("K_level must be >= 2.")
  if (!is.numeric(q) || length(q) != 1 || q < 1 || q != as.integer(q)) stop("q must be a positive integer.")

  ratio <- (K_level - 1) / (K_final - 1)
  as.numeric(lambda_final * (ratio^(2*q - 1)))
}


#' Krige/interpolate from observed nodes using a Gaussian precision matrix
#'
#' @description
#' Conditional mean under a Gaussian model specified by a precision matrix \code{Q_new}.
#'
#' @param f_obs Numeric vector length |O| or matrix |O| x p.
#' @param idx_obs Integer vector of observed indices (subset of 1:K_new).
#' @param Q_new K_new x K_new precision matrix (dense or sparse).
#' @param nugget Nonnegative scalar: add nugget*I to Q_UU for stability.
#' @param keep_obs If TRUE, returns full length K_new (or K_new x p) with obs filled in.
#'
#' @return Full vector/matrix if keep_obs=TRUE; else list(pred_idx, f_pred).
#' @export
kriging_from_precision <- function(f_obs, idx_obs, Q_new,
                                   nugget = 0,
                                   keep_obs = TRUE) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for kriging_from_precision().")
  }

  Q_new <- if (inherits(Q_new, "Matrix")) Q_new else Matrix::Matrix(Q_new, sparse = TRUE)

  K_new <- nrow(Q_new)
  if (ncol(Q_new) != K_new) stop("Q_new must be square.")

  idx_obs <- as.integer(idx_obs)
  if (any(idx_obs < 1L | idx_obs > K_new)) stop("idx_obs out of range.")
  if (anyDuplicated(idx_obs)) stop("idx_obs contains duplicates.")

  idx_pred <- setdiff(seq_len(K_new), idx_obs)

  if (length(idx_pred) == 0L) {
    return(if (keep_obs) f_obs else list(pred_idx = integer(0), f_pred = NULL))
  }

  f_obs_mat <- if (is.null(dim(f_obs))) matrix(as.numeric(f_obs), ncol = 1) else as.matrix(f_obs)
  if (nrow(f_obs_mat) != length(idx_obs)) stop("nrow(f_obs) must equal length(idx_obs).")
  p <- ncol(f_obs_mat)

  Q_UU <- Q_new[idx_pred, idx_pred, drop = FALSE]
  Q_UO <- Q_new[idx_pred, idx_obs,  drop = FALSE]

  if (nugget > 0) {
    Q_UU <- Q_UU + nugget * Matrix::Diagonal(n = nrow(Q_UU))
  }

  rhs <- Q_UO %*% f_obs_mat

  # Solve robustly
  f_pred <- if (inherits(Q_UU, "Matrix")) {
    fac <- Matrix::Cholesky(Q_UU, LDL = FALSE)
    -Matrix::solve(fac, rhs)
  } else {
    -solve(Q_UU, rhs)
  }

  # Force strict base matrix with correct dims
  if (is.null(dim(f_pred))) {
    f_pred_mat <- matrix(as.numeric(f_pred), nrow = length(idx_pred), ncol = p)
  } else {
    f_pred_mat <- as.matrix(f_pred)
  }

  if (nrow(f_pred_mat) != length(idx_pred) || ncol(f_pred_mat) != p) {
    stop(sprintf("Internal dimension mismatch: f_pred is %d x %d but expected %d x %d.",
                 nrow(f_pred_mat), ncol(f_pred_mat), length(idx_pred), p))
  }

  if (!keep_obs) {
    out_pred <- if (p == 1) as.numeric(f_pred_mat[, 1]) else f_pred_mat
    return(list(pred_idx = idx_pred, f_pred = out_pred))
  }

  if (p == 1) {
    f_full <- rep(NA_real_, K_new)
    f_full[idx_obs]  <- f_obs_mat[, 1]
    f_full[idx_pred] <- f_pred_mat[, 1]
    return(f_full)
  } else {
    f_full <- matrix(NA_real_, nrow = K_new, ncol = p)
    f_full[idx_obs, ]  <- f_obs_mat
    f_full[idx_pred, ] <- f_pred_mat
    return(f_full)
  }
}


#' Krige SmoothEM mean functions from an observed grid to the final grid
#'
#' @param mu_list_obs List length K_obs; each element a numeric vector length d.
#' @param u_obs Numeric length K_obs.
#' @param u_final Numeric length K_final.
#' @param Q_final_1d K_final x K_final precision acting along the grid.
#' @param nugget Nonnegative scalar passed to kriging_from_precision().
#'
#' @return List(Mu_full, mu_full_list, idx_obs).
#' @export
krige_mu_list_to_full_grid <- function(mu_list_obs, u_obs, u_final, Q_final_1d,
                                       nugget = 0) {
  if (!is.list(mu_list_obs) || length(mu_list_obs) < 1) stop("mu_list_obs must be a nonempty list.")
  if (length(mu_list_obs) != length(u_obs)) stop("length(mu_list_obs) must equal length(u_obs).")

  idx_obs <- match_locations_to_grid(u_obs, u_final)

  Mu_obs <- do.call(rbind, lapply(mu_list_obs, function(v) as.numeric(v)))
  Mu_full <- kriging_from_precision(
    f_obs    = Mu_obs,
    idx_obs  = idx_obs,
    Q_new    = Q_final_1d,
    nugget   = nugget,
    keep_obs = TRUE
  )

  mu_full_list <- lapply(seq_len(nrow(Mu_full)), function(k) as.numeric(Mu_full[k, ]))
  list(Mu_full = Mu_full, mu_full_list = mu_full_list, idx_obs = idx_obs)
}


#' Progressive-resolution initialization for SmoothEM via kriging on a final grid
#'
#' @description
#' Coarse-to-fine continuation scheme for initializing SmoothEM means over nested grids.
#'
#' @inheritParams EM_algorithm
#' @param data Numeric matrix n x d.
#' @param m_max Finest grid exponent; final grid size is 2^m_max + 1.
#' @param lambda_final Penalty strength on final grid.
#' @param q RW order (e.g. 2 for RW2).
#' @param ridge Ridge added in RW precision construction.
#' @param nugget_kriging Nugget added to Q_UU during kriging solve.
#' @param coords_show Coords to visualize when plot_each_stage=TRUE.
#' @param plot_each_stage If TRUE, plot kriged mean curves each stage.
#' @param verbose If TRUE, print stage summaries and forward verbose to EM_algorithm().
#' @param include.data If TRUE, include data in returned fits.
#'
#' @return List(grid, Q_final_1d, fits, mu_full_history, mu_full_list_final, meta_history, mclust).
#' @export
progressive_smoothEM <- function(
    data,
    m_max = 6,
    lambda_final = 500,
    q = 2,
    ridge = 0,
    nugget_kriging = 0,
    tol = 1e-3,
    max_iter = 1000,
    relative_lambda = TRUE,
    modelName = "EEI",
    coords_show = c(1, 2, 3),
    plot_each_stage = TRUE,
    verbose = TRUE,
    include.data = TRUE
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for progressive_smoothEM().")
  }
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for progressive_smoothEM().")
  }

  data <- as.matrix(data)
  d <- ncol(data)
  if (d < 1) stop("data must have at least 1 column.")
  coords_show <- as.integer(coords_show)
  coords_show <- coords_show[coords_show >= 1 & coords_show <= d]

  meta_history <- list()

  # 1) hierarchical grids
  grid <- make_hierarchical_levels(m_max = m_max)
  u_final <- grid$u_final
  K_final <- grid$K_final

  # 2) final-grid precision for kriging (1D)
  # NOTE: if you have a sparse constructor, swap it in here.
  Q_final_1d <- make_random_walk_precision(
    K = K_final, d = 1, lambda = lambda_final, q = q, ridge = ridge
  )
  if (inherits(Q_final_1d, "Matrix")) {
    # fine
  } else {
    Q_final_1d <- Matrix::Matrix(Q_final_1d, sparse = TRUE)
  }

  # 3) stage 0: K=2 using Mclust
  u_obs0 <- c(0, 1)

  mclustBIC <- mclust::mclustBIC
  mcl0 <- mclust::Mclust(data, G = 2, verbose = FALSE)

  mu_list_obs0 <- list(
    as.numeric(mcl0$parameters$mean[, 1]),
    as.numeric(mcl0$parameters$mean[, 2])
  )

  # krige K=2 means to full grid
  kr0 <- krige_mu_list_to_full_grid(mu_list_obs0, u_obs0, u_final, Q_final_1d, nugget = nugget_kriging)
  mu_full_list <- kr0$mu_full_list

  fits <- list()
  mu_full_history <- list()
  mu_full_history[[1]] <- kr0$Mu_full

  meta_history[[1]] <- list(
    stage = 1,
    u_obs = u_obs0,
    mu_obs_list = mu_list_obs0,
    K_obs = length(u_obs0),
    lambda = NA_real_
  )

  if (plot_each_stage && length(coords_show) > 0) {
    graphics::matplot(u_final, kr0$Mu_full[, coords_show, drop = FALSE], type = "l", lty = 1,
                      xlab = "u", ylab = expression(mu[j](u)),
                      main = "Stage 1 (K=2): kriged means on final grid")
    for (jj in coords_show) {
      graphics::points(u_obs0, c(mu_list_obs0[[1]][jj], mu_list_obs0[[2]][jj]), pch = 16)
    }
  }

  # 4) stages m=1..m_max
  for (m in 1:m_max) {
    idx_fit <- grid$idx_levels[[m + 1]]
    u_fit   <- u_final[idx_fit]
    K_next  <- length(u_fit)

    lambda_next <- lambda_scale_for_spacing(lambda_final, K_final, K_next, q = q)

    Q_next <- make_random_walk_precision(
      K = K_next, d = d, lambda = lambda_next, q = q, ridge = ridge
    )
    if (!inherits(Q_next, "Matrix")) Q_next <- Matrix::Matrix(Q_next, sparse = TRUE)

    init_params <- make_default_init(data, K = K_next, ordering = TRUE)
    init_params$mu <- mu_full_list[idx_fit]

    rank_def <- q * d

    if (verbose) {
      cat(sprintf("\n[Stage %d] K=%d, lambda_next=%.3g (lambda_final=%.3g)\n",
                  m, K_next, lambda_next, lambda_final))
    }

    fit_m <- EM_algorithm(
      data = data,
      init_params = init_params,
      Q_prior = Q_next,
      rank_deficiency = rank_def,
      max_iter = max_iter,
      modelName = modelName,
      tol = tol,
      relative_lambda = relative_lambda,
      verbose = verbose
    )

    fits[[paste0("K_", K_next)]] <- fit_m

    mu_list_obs <- fit_m$params$mu
    kr_m <- krige_mu_list_to_full_grid(mu_list_obs, u_fit, u_final, Q_final_1d, nugget = nugget_kriging)

    mu_full_list <- kr_m$mu_full_list
    mu_full_history[[m + 1]] <- kr_m$Mu_full

    meta_history[[m + 1]] <- list(
      stage = (m + 1),
      u_obs = u_fit,
      mu_obs_list = mu_list_obs,
      K_obs = K_next,
      lambda = lambda_next
    )

    if (plot_each_stage && length(coords_show) > 0) {
      graphics::matplot(u_final, kr_m$Mu_full[, coords_show, drop = FALSE], type = "l", lty = 1,
                        xlab = "u", ylab = expression(mu[j](u)),
                        main = sprintf("Stage %d (K=%d): kriged means on final grid", m, K_next))
      for (jj in coords_show) {
        graphics::points(u_fit, vapply(mu_list_obs, function(mu_k) mu_k[jj], numeric(1)), pch = 16)
      }
    }
  }

  # if include data in fits, add it into the last fit
  fits[[length(fits)]]$data <- if (include.data) data else NULL

  list(
    grid = grid,
    Q_final_1d = Q_final_1d,
    fits = fits,
    mu_full_history = mu_full_history,
    mu_full_list_final = mu_full_list,
    meta_history = meta_history,
    mclust = mcl0
  )
}


#' Plot kriged mean curves from mu_full_history
#'
#' @param mu_full_history List of matrices (K_final x d).
#' @param u_final Numeric length K_final.
#' @param history_i Which history element to plot.
#' @param coords Which coordinates (columns) to plot.
#' @param main Plot title.
#' @param add_points Overlay observed means at design points.
#' @param u_obs,mu_obs_list Optional overrides for points.
#' @param res Optional progressive_smoothEM() result used to auto-fill points/title.
#' @return Invisibly returns a list with plot inputs.
#' @export
plot_mu_history <- function(mu_full_history,
                            u_final,
                            history_i = 1,
                            coords = c(1, 2, 3),
                            main = NULL,
                            xlab = "u (grid position)",
                            ylab = expression(mu[j](u)),
                            lty = 1,
                            add_points = FALSE,
                            u_obs = NULL,
                            mu_obs_list = NULL,
                            res = NULL,
                            pch = 16,
                            cex = 0.8,
                            legend_loc = "topright",
                            legend_prefix = "coord ",
                            bty = "n") {

  if (history_i < 1 || history_i > length(mu_full_history)) {
    stop("history_i out of range.")
  }
  Mu <- mu_full_history[[history_i]]
  if (!is.matrix(Mu)) stop("mu_full_history[[history_i]] must be a matrix (K_final x d).")

  K_final <- nrow(Mu)
  d <- ncol(Mu)
  if (length(u_final) != K_final) stop("length(u_final) must equal nrow(Mu).")

  coords <- as.integer(coords)
  if (any(coords < 1 | coords > d)) stop("coords out of range.")

  if (add_points && (is.null(u_obs) || is.null(mu_obs_list))) {
    if (!is.null(res) && !is.null(res$meta_history) &&
        history_i <= length(res$meta_history) &&
        !is.null(res$meta_history[[history_i]]$u_obs) &&
        !is.null(res$meta_history[[history_i]]$mu_obs_list)) {

      u_obs <- res$meta_history[[history_i]]$u_obs
      mu_obs_list <- res$meta_history[[history_i]]$mu_obs_list
    } else {
      warning("add_points=TRUE but u_obs/mu_obs_list not available; skipping points.")
      add_points <- FALSE
    }
  }

  if (is.null(main)) {
    if (!is.null(res) && !is.null(res$meta_history) && history_i <= length(res$meta_history)) {
      mh <- res$meta_history[[history_i]]
      lam_txt <- if (!is.null(mh$lambda) && is.finite(mh$lambda)) sprintf(", lambda=%.3g", mh$lambda) else ""
      main <- sprintf("history %d (K=%d%s)", history_i, mh$K_obs, lam_txt)
    } else {
      main <- sprintf("mu_full_history[[%d]] (selected coordinates)", history_i)
    }
  }

  graphics::matplot(u_final, Mu[, coords, drop = FALSE],
                    type = "l", lty = lty, xlab = xlab, ylab = ylab, main = main)

  if (add_points) {
    for (j in coords) {
      yj <- vapply(mu_obs_list, function(mu_k) mu_k[j], numeric(1))
      graphics::points(u_obs, yj, pch = pch, cex = cex)
    }
  }

  graphics::legend(legend_loc,
                   legend = paste0(legend_prefix, coords),
                   lty = rep(lty, length(coords)),
                   bty = bty)

  invisible(list(Mu = Mu, coords = coords, history_i = history_i,
                 u_obs = u_obs, mu_obs_list = mu_obs_list))
}


#' Plot change in a single coordinate between two histories
#'
#' @export
plot_coordinate_change <- function(mu_full_history,
                                   u_final,
                                   coord,
                                   history_i,
                                   history_j,
                                   type = c("both", "overlay", "diff"),
                                   highlight = c("both", "i", "j", "none"),
                                   highlight_u = NULL,
                                   res = NULL,
                                   main = NULL,
                                   xlab = "u (grid position)",
                                   ylab_overlay = expression(mu[j](u)),
                                   ylab_diff = expression(Delta*mu[j](u)),
                                   lty_i = 1,
                                   lty_j = 2,
                                   lty_diff = 1,
                                   pch_i = 16,
                                   pch_j = 17,
                                   pch_diff = 16,
                                   cex_pts = 0.9,
                                   add_zero_line = TRUE,
                                   legend_loc = "topright",
                                   bty = "n",
                                   ylim_overlay = NULL,
                                   pad_ylim = 0.04,
                                   include_highlight_in_ylim = TRUE,
                                   label_points = FALSE) {

  type <- match.arg(type)
  highlight <- match.arg(highlight)

  H <- length(mu_full_history)
  if (history_i < 1 || history_i > H) stop("history_i out of range.")
  if (history_j < 1 || history_j > H) stop("history_j out of range.")
  if (history_i == history_j) stop("history_i and history_j must be different.")

  Mu_i <- mu_full_history[[history_i]]
  Mu_j <- mu_full_history[[history_j]]
  if (!is.matrix(Mu_i) || !is.matrix(Mu_j)) stop("mu_full_history entries must be matrices (K_final x d).")
  if (nrow(Mu_i) != nrow(Mu_j)) stop("Different K_final across histories.")
  if (length(u_final) != nrow(Mu_i)) stop("length(u_final) must equal nrow(Mu_i).")

  d <- ncol(Mu_i)
  if (coord < 1 || coord > d) stop("coord out of range.")

  yi <- Mu_i[, coord]
  yj <- Mu_j[, coord]
  dy <- yj - yi

  get_u_from_res <- function(h) {
    if (!is.null(res) && !is.null(res$meta_history) &&
        h <= length(res$meta_history) &&
        !is.null(res$meta_history[[h]]$u_obs)) {
      return(res$meta_history[[h]]$u_obs)
    }
    NULL
  }

  u_hi_i <- u_hi_j <- NULL
  if (!is.null(highlight_u)) {
    u_hi_i <- highlight_u
    u_hi_j <- highlight_u
  } else {
    if (highlight %in% c("i", "both")) u_hi_i <- get_u_from_res(history_i)
    if (highlight %in% c("j", "both")) u_hi_j <- get_u_from_res(history_j)
  }

  idx_hi_i <- if (!is.null(u_hi_i)) match_locations_to_grid(u_hi_i, u_final) else integer(0)
  idx_hi_j <- if (!is.null(u_hi_j)) match_locations_to_grid(u_hi_j, u_final) else integer(0)

  if (is.null(main)) {
    if (!is.null(res) && !is.null(res$meta_history) &&
        history_i <= length(res$meta_history) &&
        history_j <= length(res$meta_history)) {
      Ki <- res$meta_history[[history_i]]$K_obs
      Kj <- res$meta_history[[history_j]]$K_obs
      main <- sprintf("coord %d: history %d (K=%s) → history %d (K=%s)",
                      coord, history_i, ifelse(is.null(Ki), "?", Ki),
                      history_j, ifelse(is.null(Kj), "?", Kj))
    } else {
      main <- sprintf("coord %d: history %d → %d", coord, history_i, history_j)
    }
  }

  if (type %in% c("overlay", "both")) {
    if (is.null(ylim_overlay)) {
      y_all <- c(yi, yj)
      if (include_highlight_in_ylim) {
        if (length(idx_hi_i) > 0) y_all <- c(y_all, yi[idx_hi_i])
        if (length(idx_hi_j) > 0) y_all <- c(y_all, yj[idx_hi_j])
      }
      rng <- range(y_all, finite = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2])) stop("Non-finite values in curves; cannot set ylim.")
      if (rng[1] == rng[2]) {
        eps <- ifelse(rng[1] == 0, 1, abs(rng[1])) * 0.05
        rng <- rng + c(-eps, eps)
      } else {
        pad <- diff(rng) * pad_ylim
        rng <- rng + c(-pad, pad)
      }
      ylim_use <- rng
    } else {
      ylim_use <- ylim_overlay
    }

    graphics::plot(u_final, yi, type = "l", lty = lty_i,
                   xlab = xlab, ylab = ylab_overlay,
                   ylim = ylim_use,
                   main = if (type == "overlay") main else paste0(main, " (overlay)"))
    graphics::lines(u_final, yj, lty = lty_j)

    if (highlight %in% c("i", "both") && length(idx_hi_i) > 0) {
      graphics::points(u_final[idx_hi_i], yi[idx_hi_i], pch = pch_i, cex = cex_pts)
      if (label_points) graphics::text(u_final[idx_hi_i], yi[idx_hi_i], labels = idx_hi_i, pos = 3, cex = 0.7)
    }
    if (highlight %in% c("j", "both") && length(idx_hi_j) > 0) {
      graphics::points(u_final[idx_hi_j], yj[idx_hi_j], pch = pch_j, cex = cex_pts)
      if (label_points) graphics::text(u_final[idx_hi_j], yj[idx_hi_j], labels = idx_hi_j, pos = 3, cex = 0.7)
    }

    graphics::legend(legend_loc,
                     legend = c(paste0("history ", history_i),
                                paste0("history ", history_j)),
                     lty = c(lty_i, lty_j),
                     bty = bty)
  }

  if (type %in% c("diff", "both")) {
    graphics::plot(u_final, dy, type = "l", lty = lty_diff,
                   xlab = xlab, ylab = ylab_diff,
                   main = if (type == "diff") main else paste0(main, " (difference)"))
    if (add_zero_line) graphics::abline(h = 0, lty = 2)

    idx_hi <- sort(unique(c(idx_hi_i, idx_hi_j)))
    if (length(idx_hi) > 0) {
      graphics::points(u_final[idx_hi], dy[idx_hi], pch = pch_diff, cex = cex_pts)
      if (label_points) graphics::text(u_final[idx_hi], dy[idx_hi], labels = idx_hi, pos = 3, cex = 0.7)
    }
  }

  invisible(list(coord = coord, history_i = history_i, history_j = history_j,
                 mu_i = yi, mu_j = yj, delta = dy,
                 idx_hi_i = idx_hi_i, idx_hi_j = idx_hi_j))
}


#' Extract a (u, gamma) block for a given progressive history index
#'
#' @export
get_history_block <- function(res, h, normalize_gamma = TRUE) {
  if (is.null(res$meta_history)) stop("res must contain meta_history.")
  if (h < 1 || h > length(res$meta_history)) stop("h out of range.")

  mh <- res$meta_history[[h]]
  if (is.null(mh$K_obs) || is.null(mh$u_obs)) stop("meta_history[[h]] must contain K_obs and u_obs.")

  K <- mh$K_obs
  u <- mh$u_obs

  gamma <- NULL
  if (h == 1) {
    if (!is.null(res$mclust) && !is.null(res$mclust$z)) {
      gamma <- res$mclust$z
    } else {
      stop("history 1: cannot find gamma. Expect res$mclust$z.")
    }
  } else {
    key <- paste0("K_", K)
    if (!is.null(res$fits) && !is.null(res$fits[[key]]) && !is.null(res$fits[[key]]$gamma)) {
      gamma <- res$fits[[key]]$gamma
    } else if (!is.null(res$fits_history) && !is.null(res$fits_history[[h]]) && !is.null(res$fits_history[[h]]$gamma)) {
      gamma <- res$fits_history[[h]]$gamma
    } else {
      stop(sprintf("history %d: cannot find gamma. Expected res$fits[['%s']]$gamma.", h, key))
    }
  }

  if (ncol(gamma) != length(u)) {
    stop(sprintf("gamma ncol (%d) != length(u_obs) (%d) at history %d.", ncol(gamma), length(u), h))
  }

  if (normalize_gamma) {
    rs <- rowSums(gamma)
    gamma <- gamma / pmax(rs, 1e-12)
  }

  list(h = h, K = K, u = u, gamma = gamma)
}


#' Compare posterior position summaries across two histories
#'
#' @export
compare_positions_by_history <- function(res, h1, h2,
                                         type = c("mean", "max"),
                                         add_identity = TRUE,
                                         main = NULL) {
  type <- match.arg(type)

  b1 <- get_history_block(res, h1, normalize_gamma = TRUE)
  b2 <- get_history_block(res, h2, normalize_gamma = TRUE)

  if (type == "mean") {
    pos1 <- as.numeric(b1$gamma %*% b1$u)
    pos2 <- as.numeric(b2$gamma %*% b2$u)
  } else {
    pos1 <- b1$u[max.col(b1$gamma, ties.method = "first")]
    pos2 <- b2$u[max.col(b2$gamma, ties.method = "first")]
  }

  if (is.null(main)) {
    main <- sprintf("%s position: history %d (K=%d) vs history %d (K=%d)",
                    type, h1, b1$K, h2, b2$K)
  }

  graphics::plot(pos2 ~ pos1,
                 xlab = sprintf("%s position (h=%d, K=%d)", type, h1, b1$K),
                 ylab = sprintf("%s position (h=%d, K=%d)", type, h2, b2$K),
                 main = main)

  if (add_identity) graphics::abline(0, 1, lty = 3)
  graphics::abline(stats::lm(pos2 ~ pos1), lty = 2)

  invisible(list(pos1 = pos1, pos2 = pos2, cor = stats::cor(pos1, pos2, use = "complete.obs")))
}
