#' Plot data embedding with SmoothEM component means
#'
#' Visualize an EM fit (from \code{EM_algorithm()}) on 1D or 2D selected coordinates.
#' For 2D, draws arrows connecting component means in order. For 1D, draws a scatter
#' of posterior position vs the selected coordinate and overlays the mean curve.
#'
#' @param fit A fitted object returned by \code{EM_algorithm()}.
#' @param X Numeric matrix (n x d). The data used to fit the model.
#' @param dims Integer vector of length 1 or 2 indicating which coordinates (columns of X) to plot.
#' @param position Optional numeric vector of length n giving x-axis "position".
#'   If NULL, uses \code{fit$position} if present; otherwise uses posterior mean position computed as
#'   \code{fit$gamma \%*\% seq_len(K)}.
#' @param use_posterior_mean Logical; if TRUE (default), use posterior mean position.
#'   If FALSE, use MAP position via \code{max.col(fit$gamma)}.
#' @param pch,col,cex Point style for data scatter.
#' @param mu_pch,mu_col,mu_cex Mean point style.
#' @param arrow_col,arrow_lwd,arrow_len Arrow style (2D only).
#' @param line_col,line_lwd Mean curve style (1D only).
#' @param add Logical; if TRUE, add to existing plot.
#' @param xlab,ylab,main Labels; if NULL, auto-generated.
#' @param ... Passed to \code{plot()} for the scatter.
#'
#' @return Invisibly returns a list with \code{dims}, \code{pos}, \code{mu_mat}.
#' @export
plot_EM_embedding <- function(
    fit,
    X,
    dims = c(1, 2),
    position = NULL,
    use_posterior_mean = TRUE,
    pch = 19,
    col = "grey60",
    cex = 0.7,
    mu_pch = 8,
    mu_col = "orange",
    mu_cex = 1,
    arrow_col = "orange",
    arrow_lwd = 3,
    arrow_len = 0.08,
    line_col = "orange",
    line_lwd = 2,
    add = FALSE,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    ...
) {
  X <- as.matrix(X)
  d <- ncol(X)

  dims <- as.integer(dims)
  if (length(dims) < 1 || length(dims) > 2) {
    stop("dims must have length 1 or 2.")
  }
  if (any(dims < 1 | dims > d)) stop("dims out of range for columns of X.")

  if (is.null(fit$params) || is.null(fit$params$mu)) stop("fit must contain fit$params$mu.")
  if (is.null(fit$gamma)) stop("fit must contain fit$gamma.")

  K <- length(fit$params$mu)
  mu_mat_full <- do.call(rbind, fit$params$mu)  # K x d
  if (!is.matrix(mu_mat_full) || nrow(mu_mat_full) != K) stop("Cannot form mu matrix from fit$params$mu.")

  # ----- position (length n)
  n <- nrow(X)
  if (is.null(position)) {
    if (!is.null(fit$position)) {
      position <- fit$position
    } else {
      if (use_posterior_mean) {
        position <- as.numeric(fit$gamma %*% seq_len(K))
      } else {
        position <- max.col(fit$gamma, ties.method = "first")
      }
    }
  }
  if (length(position) != n) stop("position must have length nrow(X).")

  # select dims
  Xp <- X[, dims, drop = FALSE]
  mup <- mu_mat_full[, dims, drop = FALSE]

  # ----- labels
  if (is.null(main)) {
    main <- if (length(dims) == 2) {
      sprintf("EM embedding on dims (%d, %d)", dims[1], dims[2])
    } else {
      sprintf("EM embedding on dim %d", dims[1])
    }
  }

  if (length(dims) == 2) {
    # ===== 2D plot =====
    if (is.null(xlab)) xlab <- sprintf("X[, %d]", dims[1])
    if (is.null(ylab)) ylab <- sprintf("X[, %d]", dims[2])

    if (!add) {
      plot(Xp[, 1], Xp[, 2],
           pch = pch, col = col, cex = cex,
           xlab = xlab, ylab = ylab, main = main, ...)
    } else {
      points(Xp[, 1], Xp[, 2], pch = pch, col = col, cex = cex, ...)
    }

    # arrows between consecutive means (in component order)
    for (k in seq_len(nrow(mup) - 1)) {
      dx <- mup[k + 1, 1] - mup[k, 1]
      dy <- mup[k + 1, 2] - mup[k, 2]
      if (sqrt(dx^2 + dy^2) > 1e-8) {
        arrows(mup[k, 1], mup[k, 2],
               mup[k + 1, 1], mup[k + 1, 2],
               col = arrow_col, lwd = arrow_lwd, length = arrow_len)
      }
    }
    points(mup[, 1], mup[, 2], pch = mu_pch, col = mu_col, cex = mu_cex, lwd = 1)

  } else {
    # ===== 1D plot =====
    j <- dims[1]
    if (is.null(xlab)) xlab <- "posterior position"
    if (is.null(ylab)) ylab <- sprintf("X[, %d]", j)

    if (!add) {
      plot(position, Xp[, 1],
           pch = pch, col = col, cex = cex,
           xlab = xlab, ylab = ylab, main = main, ...)
    } else {
      points(position, Xp[, 1], pch = pch, col = col, cex = cex, ...)
    }

    # overlay mean curve: x = component index (or component "position"), y = mean on dim j
    # Here we use component index 1..K to match position scale you use.
    lines(seq_len(K), mup[, 1], col = line_col, lwd = line_lwd)
    points(seq_len(K), mup[, 1], pch = mu_pch, col = mu_col, cex = mu_cex)
  }

  invisible(list(dims = dims, pos = position, mu_mat = mup))
}














