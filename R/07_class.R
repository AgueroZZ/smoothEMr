#' Construct a smooth_em object from an EM_algorithm fit
#'
#' @param fit Output list from \code{EM_algorithm()}.
#' @param Q_prior Optional precision matrix used in fitting (legacy; discouraged for continuing).
#' @param Q_base Optional *base* precision matrix (without lambda). Recommended.
#' @param lambda Optional penalty strength; used with \code{Q_base}.
#' @param q Optional RW order, if relevant.
#' @param ridge Optional ridge used in building \code{Q_base}.
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
    Q_base  = NULL,
    lambda  = NULL,
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

  # ---- normalize prior storage ----
  # Prefer Q_base + lambda. Keep Q_prior only for backward compatibility.
  if (is.null(Q_base) && !is.null(Q_prior)) {
    Q_base <- Q_prior
    if (is.null(lambda)) lambda <- 1
  }
  if (is.null(lambda)) lambda <- 1
  lambda <- as.numeric(lambda)

  iter0 <- length(fit$elbo_trace %||% numeric(0))
  lambda_trace0 <- if (iter0 > 0L) rep(lambda, iter0) else numeric(0)

  obj <- list(
    params = params,
    gamma  = fit$gamma,
    data = fit$data %||% NULL,

    elbo_trace   = fit$elbo_trace   %||% numeric(0),
    loglik_trace = fit$loglik_trace %||% numeric(0),
    lambda_trace = fit$lambda_trace %||% lambda_trace0,

    iter = iter0,

    prior = list(
      # preferred fields
      Q_base = Q_base,
      lambda = lambda,

      # legacy (optional)
      Q_prior = Q_prior,

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
      mstep_iterate_once = FALSE,

      # adaptive-lambda knobs (default off)
      adapt_lambda = FALSE,
      lambda_min = 1e-8,
      lambda_max = 1e8
    ),

    meta = meta
  )

  class(obj) <- "smooth_em"
  obj
}

#' Run SmoothEM for a given number of iterations on a smooth_em object
#'
#' @param object A \code{smooth_em} object created by \code{as_smooth_em()}.
#' @param data Numeric matrix (n x d).
#' @param iter Integer >= 1; number of (E-step + M-step) iterations to run.
#' @param record Logical; whether to append objective values to traces.
#' @param check_decrease Logical; if TRUE, rollback if ELBO decreases materially.
#' @param tol_decrease Numeric; tolerance for considering ELBO decrease (default 1e-10).
#' @param adaptive Logical; if TRUE, update \code{lambda} each iteration (profile-style update).
#' @param lambda_min,lambda_max Bounds for adaptive lambda (ignored if adaptive=FALSE).
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
                        adaptive = TRUE,
                        lambda_min = NULL,
                        lambda_max = NULL,
                        verbose = FALSE) {

  if (!inherits(object, "smooth_em")) stop("object must be a 'smooth_em' object.")

  if (is.null(data)) {
    data <- object$data
    if (is.null(data)) stop("data must be provided either in the object or as an argument.")
  }
  data <- as.matrix(data)

  iter <- as.integer(iter)
  if (length(iter) != 1L || is.na(iter) || iter < 1L) stop("iter must be an integer >= 1.")

  # ---- pull settings ----
  Q_base          <- object$prior$Q_base %||% NULL
  lambda          <- as.numeric(object$prior$lambda %||% 1)
  rank_deficiency <- object$prior$rank_deficiency %||% 0

  modelName       <- object$control$modelName %||% "VVV"
  relative_lambda <- isTRUE(object$control$relative_lambda)
  eigen_tol       <- object$control$eigen_tol
  nugget          <- object$control$nugget %||% 0
  max_inner       <- object$control$max_inner %||% 10
  inner_tol       <- object$control$inner_tol %||% 1e-6

  # bounds: argument overrides object$control if provided
  if (is.null(lambda_min)) lambda_min <- object$control$lambda_min %||% 1e-8
  if (is.null(lambda_max)) lambda_max <- object$control$lambda_max %||% 1e8
  lambda_min <- as.numeric(lambda_min)
  lambda_max <- as.numeric(lambda_max)

  eps_quad <- 1e-12

  # ---- ensure cached params ----
  params <- object$params
  if (is.null(params$invSigma) || is.null(params$logdet)) {
    params <- init_cov_cache_fast(params)
  }

  # ---- helpers ----
  compute_Q_eval_base <- function(params, Q_base, relative_lambda) {
    Q_eval <- Q_base
    if (relative_lambda && !is.null(Q_eval)) {
      same <- all(vapply(params$sigma, function(S) isTRUE(all.equal(S, params$sigma[[1]])), logical(1)))
      if (!same) stop("relative_lambda requires identical covariances across clusters (typically EEI).")

      K <- length(params$pi)
      sigma_vec <- diag(params$sigma[[1]])
      scale_vec <- rep(1 / sqrt(pmax(sigma_vec, 1e-12)), times = K)
      Sscale <- Matrix::Diagonal(x = scale_vec)
      Q_eval <- Sscale %*% Q_eval %*% Sscale
    }
    Q_eval
  }

  compute_Q_eval <- function(params, Q_base, lambda, relative_lambda) {
    Qb <- compute_Q_eval_base(params, Q_base, relative_lambda)
    if (is.null(Qb)) return(NULL)
    lambda * Qb
  }

  stack_U_from_params <- function(params) {
    if (is.list(params$mu)) {
      unlist(lapply(params$mu, function(m) as.numeric(m)), use.names = FALSE)
    } else {
      as.numeric(params$mu)
    }
  }

  update_lambda_star <- function(params, Q_base, relative_lambda, rank_deficiency) {
    Qb <- compute_Q_eval_base(params, Q_base, relative_lambda)
    if (is.null(Qb)) return(NA_real_)

    U_vec <- stack_U_from_params(params)
    quad  <- as.numeric(crossprod(U_vec, Qb %*% U_vec))

    DK <- ncol(Qb)
    r  <- as.integer(DK - (rank_deficiency %||% 0L))
    r  <- max(r, 1L)

    r / pmax(quad, eps_quad)
  }

  # ---- main loop ----
  for (tt in seq_len(iter)) {

    gamma <- ESTEP(data, params)

    Q_prior_for_M <- if (is.null(Q_base)) NULL else (lambda * Q_base)

    new_params <- MSTEP(
      data = data,
      gamma = gamma,
      params = params,
      Q_prior = Q_prior_for_M,
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

    # ---- adaptive lambda (opt-in) ----
    if (isTRUE(adaptive) && !is.null(Q_base)) {
      lambda_star <- update_lambda_star(params, Q_base, relative_lambda, rank_deficiency)
      if (is.finite(lambda_star)) {
        lambda <- min(max(lambda_star, lambda_min), lambda_max)
        object$prior$lambda <- lambda
      }
    }

    if (record) {
      Q_eval <- compute_Q_eval(params, Q_base, lambda, relative_lambda)

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
            cat(sprintf("ELBO decreased at step %d: %.6f -> %.6f. Rolling back.\n",
                        tt, prev_elbo, elbo))
          }
          params <- last_params
          gamma  <- ESTEP(data, params)
          object$params <- params
          object$gamma  <- gamma
          return(object)
        }
      }

      object$loglik_trace <- c(object$loglik_trace %||% numeric(0), ll)
      object$elbo_trace   <- c(object$elbo_trace   %||% numeric(0), elbo)
      object$lambda_trace <- c(object$lambda_trace %||% numeric(0), lambda)
      object$iter <- length(object$elbo_trace)

      if (verbose) {
        cat(sprintf("do_smoothEM step %d/%d: penLogLik=%.6f, ELBO=%.6f, lambda=%.4g%s\n",
                    tt, iter, ll, elbo, lambda,
                    if (isTRUE(adaptive)) " (adaptive)" else ""))
      }
    }

    object$gamma <- gamma
  }

  object$params <- params
  object$prior$Q_base <- Q_base
  object$prior$lambda <- lambda
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
    adaptive = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)

  method <- match.arg(method)

  num_iter <- as.integer(num_iter)
  if (length(num_iter) != 1L || is.na(num_iter) || num_iter < 1L) {
    stop("num_iter must be a single integer >= 1.")
  }

  # rank deficiency for RW(q) separable prior: q per coordinate
  rank_def <- as.integer(rw_q * d)

  # ------------------------------------------------------------
  # helper: take an EM fit (1-iter) + build smooth_em + continue
  # ------------------------------------------------------------
  finalize_and_continue <- function(fit, Q_base, lambda, meta_init) {
    obj <- as_smooth_em(
      fit = fit,
      Q_base = Q_base,
      lambda = lambda,
      q = rw_q,
      ridge = ridge,
      relative_lambda = relative_lambda,
      modelName = modelName,
      eigen_tol = eigen_tol,
      rank_deficiency = rank_def
    )

    obj$meta <- obj$meta %||% list()
    obj$meta$init <- meta_init

    # record whether you *intend* to adapt later
    obj$control <- obj$control %||% list()
    obj$control$adapt_lambda <- isTRUE(adaptive)

    # initialize lambda_trace to align with current trace length (warm start)
    if (is.null(obj$lambda_trace)) obj$lambda_trace <- numeric(0)
    if (length(obj$lambda_trace) == 0L && length(obj$elbo_trace %||% numeric(0)) > 0L) {
      obj$lambda_trace <- rep(obj$prior$lambda %||% lambda, length(obj$elbo_trace))
    }

    # Continue with do_smoothEM for remaining iterations (this is where adaptive can be turned on)
    if (num_iter > 1L) {
      obj <- do_smoothEM(
        object = obj,
        data = X,
        iter = num_iter - 1L,
        adaptive = adaptive,
        record = TRUE,
        verbose = FALSE
      )
    }

    obj
  }

  # -------------------------
  # non-multi-scale branch
  # -------------------------
  if (method != "multi_scale") {

    if (is.null(K)) {
      K <- min(30L, floor(n / 5))
      K <- max(2L, as.integer(K))
    } else {
      K <- as.integer(K)
      if (length(K) != 1L || is.na(K) || K < 2L) stop("K must be a single integer >= 2.")
    }

    init_params <- initialize_ordering(X = X, K = K, method = method, ...)

    # store base precision + lambda (recommended)
    Q_base <- make_random_walk_precision(
      K = K, d = d, lambda = 1, q = rw_q, ridge = ridge
    )
    Q_prior <- lambda * Q_base

    # warm start: run exactly ONE EM outer iteration
    fit <- EM_algorithm(
      data = X,
      init_params = init_params,
      Q_prior = Q_prior,
      max_iter = 1L,
      tol = 0,
      modelName = modelName,
      eigen_tol = eigen_tol,
      rank_deficiency = rank_def,
      nugget = nugget,
      relative_lambda = relative_lambda,
      verbose = FALSE,
      include.data = include.data
    )

    return(finalize_and_continue(
      fit = fit,
      Q_base = Q_base,
      lambda = lambda,
      meta_init = list(
        method  = method,
        details = list(K = K)
      )
    ))
  }

  # -------------------------
  # multi-scale branch
  # -------------------------
  prog <- progressive_smoothEM(
    data = X,
    m_max = as.integer(m_max),
    lambda_final = lambda,
    q = rw_q,
    ridge = ridge,
    tol = 0,
    max_iter = 1L,              # <- IMPORTANT: only 1 iter per stage
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
  if (is.null(fit_final)) fit_final <- prog$fits[[length(prog$fits)]]
  if (is.null(fit_final)) stop("progressive_smoothEM did not return a final-stage fit in $fits.")

  Q_base_final <- make_random_walk_precision(
    K = K_final, d = d, lambda = 1, q = rw_q, ridge = ridge
  )

  obj <- finalize_and_continue(
    fit = fit_final,
    Q_base = Q_base_final,
    lambda = lambda,
    meta_init = list(
      method  = "multi_scale",
      details = list(m_max = m_max, K_final = K_final)
    )
  )

  if (isTRUE(keep_history)) {
    obj$meta$init$details$progressive <- prog
  }

  obj
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
#' @param methods Character vector of methods to try.
#' @param num_iter Integer >= 1. Number of EM iterations to run for each method.
#' @param num_cores Integer >= 1. Number of cores for parallel execution.
#' @param m_max Integer >= 1. Used by \code{multi_scale} and to set default \code{K}.
#' @param K Optional integer grid size for non-multi-scale methods. If NULL, uses \code{2^m_max+1}.
#' @param seed Optional base seed for reproducibility.
#' @param adaptive Logical; whether to enable adaptive lambda during the continuation phase
#'   (handled inside \code{initialize_smoothEM()} / \code{do_smoothEM()}).
#' @param quiet Logical; suppress messages from workers.
#' @param ... Extra args passed to \code{initialize_smoothEM()}.
#'
#' @return A named list of \code{smooth_em} objects (or \code{NULL} for failed fits),
#' with attribute \code{summary}.
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
    adaptive = TRUE,
    quiet = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required.")
  }

  X <- as.matrix(X)
  methods <- unique(as.character(methods))

  num_iter  <- as.integer(num_iter)
  num_cores <- as.integer(num_cores)
  m_max     <- as.integer(m_max)

  if (length(num_iter) != 1L || is.na(num_iter) || num_iter < 1L) stop("num_iter must be integer >= 1.")
  if (length(num_cores) != 1L || is.na(num_cores) || num_cores < 1L) stop("num_cores must be integer >= 1.")
  if (length(m_max) != 1L || is.na(m_max) || m_max < 1L) stop("m_max must be integer >= 1.")

  adaptive <- isTRUE(adaptive)

  # capture ...
  dots <- list(...)
  if ("adaptive" %in% names(dots)) {
    warning("Argument 'adaptive' was provided in ...; using the explicit adaptive= argument instead.")
    dots$adaptive <- NULL
  }

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
          m_max = m_max,
          num_iter = num_iter,
          adaptive = adaptive
        ), dots)
      )
    } else {
      do.call(
        initialize_smoothEM,
        c(list(
          X = X,
          method = "multi_scale",
          m_max = m_max,
          num_iter = num_iter,
          adaptive = adaptive
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

    out <- parallel::mclapply(methods, function(mm) {
      tryCatch(worker_one(mm), error = function(e) {
        if (!quiet) message(sprintf("[parallel_initial] %s failed: %s", mm, e$message))
        NULL
      })
    }, mc.cores = num_cores)

    names(out) <- methods
    results <- out

  } else {

    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Ensure functions are available on workers.
    # If you're running inside the package, workers need the package loaded.
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(smoothEMr))
      NULL
    })

    parallel::clusterExport(
      cl,
      varlist = c("X", "K_default", "m_max", "num_iter", "method_seeds",
                  "quiet", "dots", "adaptive", "worker_one"),
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
  get_lambda_last <- function(obj) {
    if (is.null(obj)) return(NA_real_)
    tr <- obj$lambda_trace %||% numeric(0)
    if (length(tr) > 0) return(tail(tr, 1))
    obj$prior$lambda %||% NA_real_
  }

  sum_df <- data.frame(
    method = methods,
    K = vapply(results, function(obj) {
      if (is.null(obj) || is.null(obj$params$pi)) NA_integer_ else length(obj$params$pi)
    }, integer(1)),
    iter = vapply(results, function(obj) {
      if (is.null(obj)) NA_integer_ else (obj$iter %||% length(obj$elbo_trace %||% numeric(0)))
    }, integer(1)),
    lambda_last = vapply(results, get_lambda_last, numeric(1)),
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
#' @param X Numeric matrix (n x d).
#' @param methods Methods to try.
#' @param num_iter Number of iterations for each fit.
#' @param num_cores Number of cores.
#' @param m_max Used for multi_scale and default K for others.
#' @param K Optional K for non-multi-scale methods.
#' @param adaptive Logical; whether to enable adaptive lambda during fitting.
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
    adaptive = TRUE,
    plot = FALSE,
    two_panel = FALSE,
    seed = NULL,
    quiet = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  adaptive <- isTRUE(adaptive)

  fits <- parallel_initial(
    X = X,
    methods = methods,
    num_iter = num_iter,
    num_cores = num_cores,
    m_max = m_max,
    K = K,
    seed = seed,
    adaptive = adaptive,
    quiet = quiet,
    ...
  )
  sum_df <- attr(fits, "summary")

  elbo_last <- sum_df$elbo_last
  elbo_last[!is.finite(elbo_last)] <- -Inf
  if (all(elbo_last == -Inf)) stop("All initializations failed or produced empty ELBO traces.")

  best_idx <- which.max(elbo_last)
  best_method <- sum_df$method[best_idx]
  best <- fits[[best_method]]
  if (is.null(best)) stop("Best method returned NULL (unexpected).")

  if (is.null(best$meta)) best$meta <- list()
  best$meta$initial_search <- list(
    best_method = best_method,
    adaptive = adaptive,
    summary = sum_df,
    fits = fits
  )

  if (isTRUE(plot)) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    methods_use <- methods
    nM <- length(methods_use)

    cols <- grDevices::hcl.colors(nM, palette = "Dark 3")
    names(cols) <- methods_use

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
      if (hi <= lo) { lo <- lo - 1; hi <- hi + 1 }
      rng <- hi - lo
      c(lo - pad * rng, hi + pad * rng)
    }

    plot_trace_family <- function(traces, ylab, main,
                                  best_method = NULL,
                                  robust_ylim = TRUE,
                                  ...) {
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

    if (two_panel) par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

    elbo_traces <- lapply(methods_use, function(mm) fits[[mm]]$elbo_trace %||% numeric(0))
    names(elbo_traces) <- methods_use
    plot_trace_family(
      traces = elbo_traces,
      ylab = "ELBO",
      main = "ELBO traces across initializations",
      best_method = best_method,
      robust_ylim = TRUE
    )
    mtext(sprintf("Best: %s (last ELBO = %.6f, adaptive=%s)",
                  best_method, sum_df$elbo_last[best_idx], as.character(adaptive)),
          side = 3, line = 0.2)

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
