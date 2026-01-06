
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

# helper: a tiny null-coalescing operator (internal)
`%||%` <- function(a, b) if (!is.null(a)) a else b


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


#' Construct a smooth_em object from an EM_algorithm fit
#'
#' @param fit Output list from \code{EM_algorithm()}.
#' @param Q_prior Optional precision matrix used in fitting (stored for continuing).
#' @param lambda Optional penalty strength used to build \code{Q_prior}.
#' @param q Optional RW order, if relevant.
#' @param ridge Optional ridge used in building \code{Q_prior}.
#' @param modelName Covariance model used in M-step (e.g. "EEI").
#' @param relative_lambda Logical; whether relative-lambda scaling is used.
#' @param rank_deficiency Rank deficiency used in generalized logdet / EEI update.
#' @param eigen_tol Optional tolerance for generalized logdet.
#' @param nugget Diagonal jitter used in covariance updates.
#' @param max_inner,inner_tol M-step inner loop controls.
#' @param meta Optional list of extra metadata to store.
#'
#' @return An object of class \code{smooth_em}.
#' @export
as_smooth_em <- function(
    fit,
    Q_prior = NULL,
    lambda = NULL,
    q = NULL,
    ridge = NULL,
    modelName = NULL,
    relative_lambda = FALSE,
    rank_deficiency = 0,
    eigen_tol = NULL,
    nugget = 0,
    max_inner = 10,
    inner_tol = 1e-6,
    meta = NULL
) {
  if (!is.list(fit) || is.null(fit$params) || is.null(fit$gamma)) {
    stop("fit must look like EM_algorithm output with fields $params and $gamma.")
  }

  params <- fit$params
  if (is.null(params$invSigma) || is.null(params$logdet)) {
    params <- init_cov_cache_fast(params)
  }

  obj <- list(
    params = params,
    gamma  = fit$gamma,
    data = fit$data %||% NULL,
    elbo_trace   = fit$elbo_trace   %||% numeric(0),
    loglik_trace = fit$loglik_trace %||% numeric(0),
    iter = length(fit$elbo_trace %||% numeric(0)),
    prior = list(
      Q_prior = Q_prior,
      lambda = lambda,
      rw_q = q,
      ridge = ridge,
      rank_deficiency = rank_deficiency
    ),
    control = list(
      modelName = modelName,
      relative_lambda = isTRUE(relative_lambda),
      eigen_tol = eigen_tol,
      nugget = nugget,
      max_inner = as.integer(max_inner),
      inner_tol = inner_tol,
      mstep_iterate_once = FALSE
    ),
    meta = meta
  )

  class(obj) <- "smooth_em"
  obj
}

# internal helper
`%||%` <- function(x, y) if (is.null(x)) y else x




#' Run SmoothEM for a given number of iterations on a smooth_em object
#'
#' @param object A \code{smooth_em} object created by \code{as_smooth_em()}.
#' @param data Numeric matrix (n x d).
#' @param iter Integer >= 1; number of (E-step + M-step) iterations to run.
#' @param record Logical; whether to append objective values to traces.
#' @param check_decrease Logical; if TRUE, rollback if ELBO decreases materially.
#' @param tol_decrease Numeric; tolerance for considering ELBO decrease (default 1e-10).
#' @param verbose Logical.
#'
#' @return Updated \code{smooth_em} object.
#' @export
do_smoothEM <- function(object,
                        data = NULL,
                        iter = 1,
                        record = TRUE,
                        check_decrease = TRUE,
                        tol_decrease = 1e-10,
                        verbose = FALSE) {

  if (!inherits(object, "smooth_em")) stop("object must be a 'smooth_em' object.")

  # if data is provided, use it; else use stored data
  if (is.null(data)) {
    data <- object$data
    if (is.null(data)) stop("data must be provided either in the object or as an argument.")
  }

  iter <- as.integer(iter)
  if (length(iter) != 1L || is.na(iter) || iter < 1L) stop("iter must be an integer >= 1.")

  # ---- pull settings ----
  Q_prior_orig    <- object$prior$Q_prior
  rank_deficiency <- object$prior$rank_deficiency %||% 0

  modelName       <- object$control$modelName %||% "VVV"
  relative_lambda <- isTRUE(object$control$relative_lambda)
  eigen_tol       <- object$control$eigen_tol
  nugget          <- object$control$nugget %||% 0
  max_inner       <- object$control$max_inner %||% 10
  inner_tol       <- object$control$inner_tol %||% 1e-6

  # ---- ensure cached params ----
  params <- object$params
  if (is.null(params$invSigma) || is.null(params$logdet)) {
    params <- init_cov_cache_fast(params)
  }

  # ---- helpers ----
  compute_Q_eval <- function(params, Q_prior_orig, relative_lambda) {
    Q_eval <- Q_prior_orig
    if (relative_lambda && !is.null(Q_prior_orig)) {
      same <- all(vapply(params$sigma, function(S) isTRUE(all.equal(S, params$sigma[[1]])), logical(1)))
      if (!same) stop("relative_lambda requires identical covariances across clusters (typically EEI).")

      K <- length(params$pi)
      sigma_vec <- diag(params$sigma[[1]])
      scale_vec <- rep(1 / sqrt(pmax(sigma_vec, 1e-12)), times = K)
      Sscale <- Matrix::Diagonal(x = scale_vec)
      Q_eval <- Sscale %*% Q_prior_orig %*% Sscale
    }
    Q_eval
  }

  # ---- main loop ----
  for (tt in seq_len(iter)) {

    # E-step
    gamma <- ESTEP(data, params)

    # M-step (single pass; can later add iterate_once = FALSE option if you want)
    new_params <- MSTEP(
      data = data,
      gamma = gamma,
      params = params,
      Q_prior = Q_prior_orig,
      relative_lambda = relative_lambda,
      modelName = modelName,
      iterate_once = TRUE,
      nugget = nugget,
      rank_deficiency = rank_deficiency,
      tol_inner = inner_tol,
      max_inner = max_inner,
      verbose = verbose
    )

    last_params <- params
    params <- init_cov_cache_fast(new_params)

    if (record) {
      Q_eval <- compute_Q_eval(params, Q_prior_orig, relative_lambda)

      ll <- compute_log_joint_observed(
        data, params, Q_eval,
        eigen_tol = eigen_tol,
        rank_deficiency = rank_deficiency
      )
      elbo <- compute_penalized_ELBO(
        data, gamma, params, Q_eval,
        eigen_tol = eigen_tol,
        rank_deficiency = rank_deficiency
      )

      prev_elbo <- tail(object$elbo_trace %||% numeric(0), 1)

      if (check_decrease && length(prev_elbo) == 1L && is.finite(prev_elbo)) {
        if (prev_elbo - elbo > tol_decrease) {
          if (verbose) {
            cat(sprintf("ELBO decreased at inner step %d: %.6f -> %.6f. Rolling back.\n",
                        tt, prev_elbo, elbo))
          }
          # rollback this iteration
          params <- last_params
          gamma  <- ESTEP(data, params)

          object$params <- params
          object$gamma  <- gamma
          return(object)
        }
      }

      object$loglik_trace <- c(object$loglik_trace %||% numeric(0), ll)
      object$elbo_trace   <- c(object$elbo_trace   %||% numeric(0), elbo)
      object$iter <- length(object$elbo_trace)

      if (verbose) {
        cat(sprintf("do_smoothEM step %d/%d: penLogLik=%.6f, ELBO=%.6f\n", tt, iter, ll, elbo))
      }
    }

    # keep latest gamma each step
    object$gamma <- gamma
  }

  object$params <- params
  object
}


#' Initialize SmoothEM (single- or multi-scale)
#'
#' @param X Numeric matrix (n x d).
#' @param method Initialization method. If not "multi_scale", it is passed to
#'   \code{initialize_ordering()}.
#' @param rw_q Random-walk order q (default 2).
#' @param lambda Penalty strength (default 10). For multi-scale, this is lambda_final.
#' @param relative_lambda Logical; default TRUE.
#' @param K Number of grid points (only used when method != "multi_scale").
#'   Default: min(30, floor(nrow(X)/5)), with a minimum of 2.
#' @param m_max Finest exponent for multi-scale grid (K_final = 2^m_max + 1).
#' @param num_iter Number of EM iterations to run (default 1).
#' @param modelName Covariance model passed to \code{EM_algorithm()}.
#' @param ridge Ridge passed to \code{make_random_walk_precision()}.
#' @param nugget Nugget passed to \code{EM_algorithm()}.
#' @param eigen_tol eigen_tol passed to objective evaluation inside \code{EM_algorithm()}.
#' @param keep_history Logical; if TRUE and method=="multi_scale", attach the full
#'   progressive result in the returned object.
#' @param include.data Logical; if TRUE, store the data matrix in the returned object.
#' @param ... Extra args passed to \code{initialize_ordering()} (when method != "multi_scale"),
#'   or to \code{progressive_smoothEM()} (when method == "multi_scale").
#'
#' @return A \code{smooth_em} object.
#' @export
initialize_smoothEM <- function(
    X,
    method = c("tSNE", "PCA", "random", "multi_scale", "fiedler"),
    rw_q = 2,
    lambda = 10,
    relative_lambda = TRUE,
    K = NULL,
    m_max = 6,
    num_iter = 1,
    modelName = "EEI",
    ridge = 0,
    nugget = 0,
    eigen_tol = NULL,
    keep_history = FALSE,
    include.data = TRUE,
    ...
) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)

  method <- match.arg(method)

  # rank deficiency for RW(q) separable prior: q per coordinate
  rank_def <- rw_q * d

  if (method != "multi_scale") {

    if (is.null(K)) {
      K <- min(30L, floor(n / 5))
      K <- max(2L, as.integer(K))
    } else {
      K <- as.integer(K)
      if (K < 2L) stop("K must be >= 2.")
    }

    init_params <- initialize_ordering(X = X, K = K, method = method, ...)

    Q_prior <- make_random_walk_precision(
      K = K, d = d, lambda = lambda, q = rw_q, ridge = ridge
    )

    fit <- EM_algorithm(
      data = X,
      init_params = init_params,
      Q_prior = Q_prior,
      max_iter = as.integer(num_iter),
      tol = 0,                     # <- avoid early stopping; run exactly num_iter (unless numerical rollback)
      modelName = modelName,
      eigen_tol = eigen_tol,
      rank_deficiency = rank_def,
      nugget = nugget,
      relative_lambda = relative_lambda,
      verbose = FALSE,
      include.data = include.data
    )

    return(as_smooth_em(
      fit = fit,
      Q_prior = Q_prior,
      lambda = lambda,
      q = rw_q,
      ridge = ridge,
      relative_lambda = relative_lambda,
      modelName = modelName,
      eigen_tol = eigen_tol,
      rank_deficiency = rank_def,
      meta = list(
        init = list(
          method  = method,
          details = list(K = K)
        )
      )
    ))

  }

  # ---- multi-scale branch ----
  prog <- progressive_smoothEM(
    data = X,
    m_max = as.integer(m_max),
    lambda_final = lambda,
    q = rw_q,
    ridge = ridge,
    tol = 0,                 # <- same idea: don’t stop early inside each stage
    max_iter = as.integer(num_iter),
    relative_lambda = relative_lambda,
    modelName = modelName,
    plot_each_stage = FALSE,
    verbose = FALSE,
    include.data = include.data,
    ...
  )

  K_final <- prog$grid$K_final
  key_final <- paste0("K_", K_final)

  fit_final <- prog$fits[[key_final]]
  if (is.null(fit_final)) {
    # fallback: take the last fit in fits list
    fit_final <- prog$fits[[length(prog$fits)]]
  }
  if (is.null(fit_final)) stop("progressive_smoothEM did not return a final-stage fit in $fits.")

  Q_prior_final <- make_random_walk_precision(
    K = K_final, d = d, lambda = lambda, q = rw_q, ridge = ridge
  )

  obj <- as_smooth_em(
    fit = fit_final,
    Q_prior = Q_prior_final,
    lambda = lambda,
    q = rw_q,
    ridge = ridge,
    relative_lambda = relative_lambda,
    modelName = modelName,
    eigen_tol = eigen_tol,
    rank_deficiency = rank_def,
    meta = list(
      init = list(
        method  = "multi_scale",
        details = list(m_max = m_max, K_final = K_final)
      )
    )
  )

  if (keep_history) obj$init$details$progressive <- prog
  obj
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

#' Run multiple SmoothEM initializations in parallel
#'
#' @description
#' Runs \code{initialize_smoothEM()} for a set of initialization methods in parallel.
#' Each method is fit for \code{num_iter} EM iterations (or per-stage iterations for
#' \code{method="multi_scale"}).
#'
#' For compatibility with \code{method="multi_scale"}, if \code{K} is not provided,
#' then non-multi-scale methods default to \code{K = 2^m_max + 1}.
#'
#' @param X Numeric matrix (n x d).
#' @param methods Character vector of methods to try. Defaults to
#'   \code{c("PCA","tSNE","random","fiedler","multi_scale")}.
#' @param num_iter Integer >= 1. Number of EM iterations to run for each method.
#' @param num_cores Integer >= 1. Number of cores for parallel execution.
#' @param m_max Integer >= 1. Used by \code{multi_scale} and to set default \code{K}.
#' @param K Optional integer grid size for non-multi-scale methods. If NULL, uses \code{2^m_max+1}.
#' @param seed Optional base seed for reproducibility. If provided, each method gets a deterministic
#'   derived seed.
#' @param quiet Logical; suppress messages from workers.
#' @param ... Extra args passed to \code{initialize_smoothEM()}.
#'
#' @return A named list of \code{smooth_em} objects (or \code{NULL} for failed fits),
#' with attributes:
#' \itemize{
#'   \item \code{summary}: a data.frame summarizing last ELBO / last objective for each method.
#' }
#'
#' @export
parallel_initial <- function(
    X,
    methods = c("PCA", "tSNE", "random", "fiedler", "multi_scale"),
    num_iter = 1,
    num_cores = 2,
    m_max = 6,
    K = NULL,
    seed = NULL,
    quiet = TRUE,
    ...
) {
  X <- as.matrix(X)
  methods <- unique(as.character(methods))

  num_iter  <- as.integer(num_iter)
  num_cores <- as.integer(num_cores)
  m_max     <- as.integer(m_max)

  if (length(num_iter) != 1L || is.na(num_iter) || num_iter < 1L) stop("num_iter must be integer >= 1.")
  if (length(num_cores) != 1L || is.na(num_cores) || num_cores < 1L) stop("num_cores must be integer >= 1.")
  if (length(m_max) != 1L || is.na(m_max) || m_max < 1L) stop("m_max must be integer >= 1.")

  # capture ...
  dots <- list(...)

  # default K consistent with multi_scale final grid
  if (is.null(K)) {
    K_default <- 2L^m_max + 1L
  } else {
    K_default <- as.integer(K)
    if (length(K_default) != 1L || is.na(K_default) || K_default < 2L) stop("K must be integer >= 2.")
  }

  # deterministic per-method seed
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (length(seed) != 1L || is.na(seed)) stop("seed must be a single integer.")
    method_seeds <- setNames(seed + seq_along(methods) * 1000L, methods)
  } else {
    method_seeds <- setNames(rep(NA_integer_, length(methods)), methods)
  }

  # worker
  worker_one <- function(method) {
    s <- method_seeds[[method]]
    if (!is.na(s)) set.seed(s)

    if (method != "multi_scale") {
      do.call(
        initialize_smoothEM,
        c(list(
          X = X,
          method = method,
          K = K_default,
          m_max = m_max,      # ignored if not multi_scale
          num_iter = num_iter
        ), dots)
      )
    } else {
      do.call(
        initialize_smoothEM,
        c(list(
          X = X,
          method = "multi_scale",
          m_max = m_max,
          num_iter = num_iter
        ), dots)
      )
    }
  }

  results <- vector("list", length(methods))
  names(results) <- methods

  # ---- backend selection ----
  if (num_cores == 1L) {

    for (mm in methods) {
      results[[mm]] <- tryCatch(worker_one(mm), error = function(e) {
        if (!quiet) message(sprintf("[parallel_initial] %s failed: %s", mm, e$message))
        NULL
      })
    }

  } else if (.Platform$OS.type != "windows") {

    # mac/linux: fork
    out <- parallel::mclapply(methods, function(mm) {
      tryCatch(worker_one(mm), error = function(e) {
        if (!quiet) message(sprintf("[parallel_initial] %s failed: %s", mm, e$message))
        NULL
      })
    }, mc.cores = num_cores)

    names(out) <- methods
    results <- out

  } else {

    # windows: PSOCK cluster
    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # IMPORTANT: do NOT use `pkg` inside worker unless you export it.
    # Just load your package explicitly.
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(smoothEMr))
      NULL
    })

    parallel::clusterExport(
      cl,
      varlist = c("X", "K_default", "m_max", "num_iter", "method_seeds", "quiet", "dots", "worker_one"),
      envir = environment()
    )

    out <- parallel::parLapply(cl, methods, function(mm) {
      tryCatch(worker_one(mm), error = function(e) {
        if (!quiet) message(sprintf("[parallel_initial] %s failed: %s", mm, e$message))
        NULL
      })
    })

    names(out) <- methods
    results <- out
  }

  # ---- build summary ----
  get_last <- function(x, field) {
    if (is.null(x) || is.null(x[[field]]) || length(x[[field]]) == 0) return(NA_real_)
    tail(x[[field]], 1)
  }

  sum_df <- data.frame(
    method = methods,
    K = vapply(results, function(obj) {
      if (is.null(obj) || is.null(obj$params$pi)) NA_integer_ else length(obj$params$pi)
    }, integer(1)),
    iter = vapply(results, function(obj) {
      if (is.null(obj)) NA_integer_ else (obj$iter %||% length(obj$elbo_trace %||% numeric(0)))
    }, integer(1)),
    elbo_last = vapply(results, get_last, numeric(1), field = "elbo_trace"),
    obj_last  = vapply(results, get_last, numeric(1), field = "loglik_trace"),
    ok = vapply(results, function(obj) !is.null(obj), logical(1)),
    stringsAsFactors = FALSE
  )

  attr(results, "summary") <- sum_df
  results
}

#' Choose the best SmoothEM initialization by ELBO
#'
#' @description
#' Runs \code{parallel_initial()} and returns the \code{smooth_em} object with the
#' largest last-iteration ELBO.
#'
#' If \code{plot=TRUE}, plots all ELBO traces (different line types) and optionally
#' overlays the penalized objective traces in a second panel.
#'
#' @param X Numeric matrix (n x d).
#' @param methods Methods to try (passed to \code{parallel_initial()}).
#' @param num_iter Number of iterations for each fit.
#' @param num_cores Number of cores.
#' @param m_max Used for multi_scale and default K for others.
#' @param K Optional K for non-multi-scale methods.
#' @param plot Logical; if TRUE, plot traces.
#' @param two_panel Logical; if TRUE, show ELBO and objective in 2 panels.
#' @param seed Optional base seed.
#' @param quiet Logical.
#' @param ... Passed to \code{initialize_smoothEM()} via \code{parallel_initial()}.
#'
#' @return A \code{smooth_em} object (best by last ELBO). The returned object gains:
#' \itemize{
#'   \item \code{$meta$initial_search}: list with the full fits and summary table.
#' }
#'
#' @export
optimize_initial <- function(
    X,
    methods = c("PCA", "tSNE", "random", "fiedler", "multi_scale"),
    num_iter = 1,
    num_cores = 2,
    m_max = 6,
    K = NULL,
    plot = FALSE,
    two_panel = FALSE,
    seed = NULL,
    quiet = TRUE,
    ...
) {
  fits <- parallel_initial(
    X = X,
    methods = methods,
    num_iter = num_iter,
    num_cores = num_cores,
    m_max = m_max,
    K = K,
    seed = seed,
    quiet = quiet,
    ...
  )
  sum_df <- attr(fits, "summary")

  # pick best by last ELBO
  elbo_last <- sum_df$elbo_last
  elbo_last[!is.finite(elbo_last)] <- -Inf
  if (all(elbo_last == -Inf)) stop("All initializations failed or produced empty ELBO traces.")

  best_idx <- which.max(elbo_last)
  best_method <- sum_df$method[best_idx]
  best <- fits[[best_method]]
  if (is.null(best)) stop("Best method returned NULL (unexpected).")

  # attach provenance
  if (is.null(best$meta)) best$meta <- list()
  best$meta$initial_search <- list(
    best_method = best_method,
    summary = sum_df,
    fits = fits
  )

  if (isTRUE(plot)) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    methods_use <- methods
    nM <- length(methods_use)

    # color palette (base R friendly)
    cols <- grDevices::hcl.colors(nM, palette = "Dark 3")
    names(cols) <- methods_use

    # helper: compute ylim with optional robust clipping + padding
    compute_ylim <- function(traces, robust = TRUE, q = c(0.01, 0.99), pad = 0.04) {
      y <- unlist(traces, use.names = FALSE)
      y <- y[is.finite(y)]
      if (!length(y)) return(NULL)

      if (robust && length(y) >= 20) {
        lo <- as.numeric(stats::quantile(y, q[1], na.rm = TRUE))
        hi <- as.numeric(stats::quantile(y, q[2], na.rm = TRUE))
      } else {
        lo <- min(y); hi <- max(y)
      }

      if (!is.finite(lo) || !is.finite(hi)) return(NULL)
      if (hi <= lo) {
        lo <- lo - 1; hi <- hi + 1
      }
      rng <- hi - lo
      c(lo - pad * rng, hi + pad * rng)
    }

    # helper: plot a family of traces
    plot_trace_family <- function(traces, ylab, main,
                                  best_method = NULL,
                                  robust_ylim = TRUE,
                                  ...) {
      # find one non-empty trace to initialize plot
      first_ok <- NULL
      for (mm in methods_use) {
        tr <- traces[[mm]] %||% numeric(0)
        if (length(tr) > 0) { first_ok <- mm; break }
      }
      if (is.null(first_ok)) {
        plot.new(); title(paste0(main, " (empty)"))
        return(invisible(NULL))
      }

      ylim <- compute_ylim(traces, robust = robust_ylim)

      tr0 <- traces[[first_ok]]
      plot(seq_along(tr0), tr0, type = "l",
           col = cols[first_ok], lwd = 2,
           xlab = "Iteration", ylab = ylab,
           main = main, ylim = ylim, ...)

      for (mm in methods_use) {
        tr <- traces[[mm]] %||% numeric(0)
        if (length(tr) > 0) {
          is_best <- (!is.null(best_method) && identical(mm, best_method))
          lines(seq_along(tr), tr,
                col = cols[mm],
                lwd = if (is_best) 4 else 2,
                lty = if (is_best) 1 else 2)
        }
      }

      legend("bottomright",
             legend = methods_use,
             col = cols[methods_use],
             lwd = ifelse(methods_use == best_method, 4, 2),
             bty = "n")
      invisible(NULL)
    }

    if (two_panel) {
      par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
    }

    # ---- ELBO panel ----
    elbo_traces <- lapply(methods_use, function(mm) fits[[mm]]$elbo_trace %||% numeric(0))
    names(elbo_traces) <- methods_use
    plot_trace_family(
      traces = elbo_traces,
      ylab = "ELBO",
      main = "ELBO traces across initializations",
      best_method = best_method,
      robust_ylim = TRUE
    )
    mtext(sprintf("Best: %s (last ELBO = %.6f)", best_method, sum_df$elbo_last[best_idx]),
          side = 3, line = 0.2)

    # ---- objective panel ----
    if (two_panel) {
      obj_traces <- lapply(methods_use, function(mm) fits[[mm]]$loglik_trace %||% numeric(0))
      names(obj_traces) <- methods_use
      plot_trace_family(
        traces = obj_traces,
        ylab = "Penalized objective",
        main = "Penalized observed-data objective traces",
        best_method = best_method,
        robust_ylim = TRUE
      )
    }
  }

  best
}



