#' Print method for smooth_em objects
#'
#' @param x A \code{smooth_em} object.
#' @param ... Unused.
#'
#' @export
print.smooth_em <- function(x, ...) {
  if (!inherits(x, "smooth_em")) stop("print.smooth_em(): x must inherit from class 'smooth_em'.")

  # ---- infer dimensions safely
  K <- NA_integer_
  n <- NA_integer_
  d <- NA_integer_

  if (!is.null(x$params) && !is.null(x$params$pi)) K <- length(x$params$pi)
  if (!is.null(x$gamma)) {
    n <- nrow(x$gamma)
    if (is.na(K)) K <- ncol(x$gamma)
  }
  if (!is.null(x$params) && !is.null(x$params$mu) && length(x$params$mu) > 0) {
    d <- length(x$params$mu[[1]])
  }

  # ---- key settings
  modelName <- x$control$modelName %||% NA_character_
  rel_lam   <- x$control$relative_lambda %||% NA
  lambda    <- x$prior$lambda %||% NA_real_
  rw_q      <- x$prior$rw_q %||% NA_integer_

  it <- x$iter %||% (if (!is.null(x$elbo_trace)) length(x$elbo_trace) else NA_integer_)

  elbo_last <- if (!is.null(x$elbo_trace) && length(x$elbo_trace) > 0) tail(x$elbo_trace, 1) else NA_real_
  ll_last   <- if (!is.null(x$loglik_trace) && length(x$loglik_trace) > 0) tail(x$loglik_trace, 1) else NA_real_

  init_method <- x$meta$init$method %||% NA_character_

  cat("<smooth_em>\n")
  cat(sprintf("  n = %s, d = %s, K = %s\n",
              ifelse(is.na(n), "?", n),
              ifelse(is.na(d), "?", d),
              ifelse(is.na(K), "?", K)))
  cat(sprintf("  model = %s, RW(q) = %s, lambda = %s, relative_lambda = %s\n",
              ifelse(is.na(modelName), "?", modelName),
              ifelse(is.na(rw_q), "?", rw_q),
              ifelse(is.na(lambda), "?", format(lambda, digits = 6)),
              ifelse(is.na(rel_lam), "?", as.character(rel_lam))))
  cat(sprintf("  iter = %s, last ELBO = %s, last penLogLik = %s\n",
              ifelse(is.na(it), "?", it),
              ifelse(is.na(elbo_last), "?", format(elbo_last, digits = 8)),
              ifelse(is.na(ll_last), "?", format(ll_last, digits = 8))))
  cat(sprintf("  init = %s\n", ifelse(is.na(init_method), "?", init_method)))

  invisible(x)
}

#' Summary method for smooth_em objects
#'
#' @param object A \code{smooth_em} object.
#' @param ... Passed through (unused).
#'
#' @return An object of class \code{summary.smooth_em}.
#' @export
summary.smooth_em <- function(object, ...) {
  if (!inherits(object, "smooth_em")) stop("summary.smooth_em(): object must inherit from class 'smooth_em'.")

  x <- object

  # ---- infer dimensions
  K <- if (!is.null(x$params$pi)) length(x$params$pi) else if (!is.null(x$gamma)) ncol(x$gamma) else NA_integer_
  n <- if (!is.null(x$gamma)) nrow(x$gamma) else NA_integer_
  d <- if (!is.null(x$params$mu) && length(x$params$mu) > 0) length(x$params$mu[[1]]) else NA_integer_

  # ---- traces
  elbo_trace <- x$elbo_trace %||% numeric(0)
  ll_trace   <- x$loglik_trace %||% numeric(0)

  elbo_last <- if (length(elbo_trace) > 0) tail(elbo_trace, 1) else NA_real_
  ll_last   <- if (length(ll_trace) > 0) tail(ll_trace, 1) else NA_real_

  elbo_diff_last <- if (length(elbo_trace) >= 2) tail(diff(elbo_trace), 1) else NA_real_

  # ---- mixture summaries
  pi_hat <- x$params$pi %||% rep(NA_real_, K %||% 0L)

  Nk <- NA
  if (!is.null(x$gamma)) {
    Nk <- colSums(x$gamma)
    names(Nk) <- paste0("k", seq_along(Nk))
  }

  # ---- ordering / grid info (if present)
  lambda    <- x$prior$lambda %||% NA_real_
  rw_q      <- x$prior$rw_q %||% NA_integer_
  ridge     <- x$prior$ridge %||% NA_real_
  rank_def  <- x$prior$rank_deficiency %||% NA_integer_
  rel_lam   <- x$control$relative_lambda %||% NA
  modelName <- x$control$modelName %||% NA_character_
  eigen_tol <- x$control$eigen_tol %||% NA_real_

  init_method <- x$meta$init$method %||% NA_character_
  init_details <- x$meta$init$details %||% list()

  out <- list(
    n = n, d = d, K = K,
    iter = x$iter %||% (if (length(elbo_trace) > 0) length(elbo_trace) else NA_integer_),
    modelName = modelName,
    lambda = lambda,
    rw_q = rw_q,
    ridge = ridge,
    relative_lambda = rel_lam,
    rank_deficiency = rank_def,
    eigen_tol = eigen_tol,
    elbo_trace = elbo_trace,
    loglik_trace = ll_trace,
    elbo_last = elbo_last,
    loglik_last = ll_last,
    elbo_diff_last = elbo_diff_last,
    pi = pi_hat,
    Nk = Nk,
    init_method = init_method,
    init_details = init_details
  )

  class(out) <- "summary.smooth_em"
  out
}

#' Print method for summary.smooth_em objects
#'
#' @param x A \code{summary.smooth_em} object.
#' @param ... Unused.
#'
#' @export
print.summary.smooth_em <- function(x, ...) {
  cat("<summary.smooth_em>\n")
  cat(sprintf("  n = %s, d = %s, K = %s\n",
              ifelse(is.na(x$n), "?", x$n),
              ifelse(is.na(x$d), "?", x$d),
              ifelse(is.na(x$K), "?", x$K)))

  cat(sprintf("  model = %s\n", ifelse(is.na(x$modelName), "?", x$modelName)))
  cat(sprintf("  RW(q) = %s, lambda = %s, ridge = %s, rank_def = %s\n",
              ifelse(is.na(x$rw_q), "?", x$rw_q),
              ifelse(is.na(x$lambda), "?", format(x$lambda, digits = 6)),
              ifelse(is.na(x$ridge), "?", format(x$ridge, digits = 6)),
              ifelse(is.na(x$rank_deficiency), "?", x$rank_deficiency)))
  cat(sprintf("  relative_lambda = %s, eigen_tol = %s\n",
              ifelse(is.na(x$relative_lambda), "?", as.character(x$relative_lambda)),
              ifelse(is.na(x$eigen_tol), "?", format(x$eigen_tol, digits = 6))))

  cat(sprintf("  iter = %s\n", ifelse(is.na(x$iter), "?", x$iter)))
  cat(sprintf("  last ELBO = %s (Δ last = %s)\n",
              ifelse(is.na(x$elbo_last), "?", format(x$elbo_last, digits = 8)),
              ifelse(is.na(x$elbo_diff_last), "?", format(x$elbo_diff_last, digits = 4))))
  cat(sprintf("  last penLogLik = %s\n",
              ifelse(is.na(x$loglik_last), "?", format(x$loglik_last, digits = 8))))

  # mixture weights
  if (!is.null(x$pi) && length(x$pi) > 0) {
    pi_show <- x$pi
    names(pi_show) <- paste0("k", seq_along(pi_show))
    cat("  pi:\n")
    cat("   ", paste(sprintf("%s=%.4f", names(pi_show), pi_show), collapse = ", "), "\n")
  }

  # effective cluster sizes
  if (!is.null(x$Nk) && all(is.finite(x$Nk))) {
    cat("  Nk (from gamma):\n")
    cat("   ", paste(sprintf("%s=%.1f", names(x$Nk), x$Nk), collapse = ", "), "\n")
  }

  cat(sprintf("  init = %s\n", ifelse(is.na(x$init_method), "?", x$init_method)))
  if (length(x$init_details) > 0) {
    # compact one-liner for common fields
    dd <- x$init_details
    keys <- intersect(names(dd), c("K", "m_max", "K_final"))
    if (length(keys) > 0) {
      cat("  init details:", paste(sprintf("%s=%s", keys, vapply(dd[keys], as.character, "")), collapse = ", "), "\n")
    }
  }

  invisible(x)
}


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



#' Plot a smooth_em object
#'
#' @description
#' Visualization for a \code{smooth_em} object.
#' \itemize{
#'   \item \code{plot_type="scatterplot"}: calls \code{plot_EM_embedding()}.
#'   \item \code{plot_type="elbo"}: plots ELBO and penalized observed-data objective traces.
#' }
#'
#' @param x A \code{smooth_em} object.
#' @param data Numeric matrix (n x d). Required when \code{plot_type="scatterplot"}.
#' @param plot_type One of \code{"scatterplot"}, \code{"elbo"}.
#' @param dims Integer vector of length 1 or 2 (only used for \code{"scatterplot"}).
#' @param two_panel Logical; if TRUE and \code{plot_type="elbo"}, draw ELBO and objective in two panels.
#' @param verbose Logical; not passed to graphics (prevents "not a graphical parameter" warnings).
#' @param ... Passed to the underlying plotting functions.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.smooth_em <- function(x,
                           data = NULL,
                           plot_type = c("scatterplot", "elbo"),
                           dims = c(1, 2),
                           two_panel = FALSE,
                           verbose = FALSE,
                           ...) {

  if (!inherits(x, "smooth_em")) stop("x must be a 'smooth_em' object.")
  plot_type <- match.arg(plot_type)

  # prevent accidental verbose passing to base graphics
  dots <- list(...)
  if ("verbose" %in% names(dots)) {
    warning("Argument 'verbose' was provided in ...; ignoring it for plotting.")
    dots$verbose <- NULL
  }

  if (plot_type == "scatterplot") {
    if (is.null(data)){
      data <- x$data
      if (is.null(data)) stop("data must be provided either as an argument or stored in the object.")
    }
    dims <- as.integer(dims)
    if (!(length(dims) %in% c(1L, 2L))) stop("dims must have length 1 or 2.")
    if (any(is.na(dims)) || any(dims < 1L) || any(dims > ncol(data))) stop("dims out of range.")

    # Use your existing visualization helper
    # (assumes plot_EM_embedding can take smooth_em directly; if not, pass list(params=..., gamma=...))
    do.call(plot_EM_embedding, c(list(fit = x, X = data, dims = dims), dots))
    return(invisible(x))
  }

  # plot_type == "elbo"
  elbo <- x$elbo_trace %||% numeric(0)
  ll   <- x$loglik_trace %||% numeric(0)

  if (length(elbo) == 0L && length(ll) == 0L) {
    warning("No traces found: elbo_trace/loglik_trace are empty.")
    return(invisible(x))
  }

  if (two_panel) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

    if (length(elbo) > 0L) {
      plot(seq_along(elbo), elbo, type = "l",
           xlab = "Iteration", ylab = "ELBO",
           main = "Penalized ELBO trace", ...)
    } else {
      plot.new(); title("Penalized ELBO trace (empty)")
    }

    if (length(ll) > 0L) {
      plot(seq_along(ll), ll, type = "l",
           xlab = "Iteration", ylab = "Penalized log-likelihood",
           main = "Penalized observed-data objective trace", ...)
    } else {
      plot.new(); title("Penalized observed-data objective trace (empty)")
    }

  } else {
    # single-panel fallback: overlay (note scales may differ)
    it_elbo <- seq_along(elbo)
    plot(it_elbo, elbo, type = "l",
         xlab = "Iteration", ylab = "Value",
         main = "ELBO (and objective, overlaid)", ...)
    if (length(ll) > 0L) lines(seq_along(ll), ll, lty = 2)
    legend("bottomright",
           legend = c("ELBO", "penalized objective"),
           lty = c(1, 2), bty = "n")
  }

  invisible(x)
}











