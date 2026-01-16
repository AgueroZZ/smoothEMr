# =========================
# Helpers
# =========================

# helper: rescale to [0,1]
.scale01 <- function(x) {
  r <- range(x, finite = TRUE)
  if (!is.finite(r[1]) || r[1] == r[2]) return(rep(0.5, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

# helper: orient sign using a reference vector (e.g., PC1)
.orient_by_ref <- function(t, ref) {
  if (length(t) != length(ref)) return(t)
  if (stats::cor(t, ref, use = "pairwise.complete.obs") < 0) -t else t
}

# helper: get PC1 quickly (prcomp is fine for now)
.pc1 <- function(X, center = TRUE, scale. = FALSE) {
  stats::prcomp(X, center = center, scale. = scale.)$x[, 1]
}

# =========================
# Initializations
# =========================

#' Default random initialization for a Gaussian mixture model
#'
#' @param X numeric matrix (n x d).
#' @param K number of mixture components.
#' @param ordering if TRUE, reorder components by a pivot dimension of mu.
#' @return list(pi, mu, sigma)
#' @export
make_default_init <- function(X, K, ordering = TRUE) {
  X <- as.matrix(X)
  if (nrow(X) < 1) stop("X must have at least one row.")
  if (K < 1) stop("K must be >= 1.")

  d <- ncol(X)
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)

  out <- list(
    pi    = rep(1 / K, K),
    mu    = lapply(seq_len(K), function(k) stats::runif(d, mins, maxs)),
    sigma = lapply(seq_len(K), function(k) diag(d))
  )

  if (ordering) {
    mu_mat <- do.call(rbind, out$mu)
    pivot_index <- which.max(apply(mu_mat, 2, function(v) max(v) - min(v)))
    order_idx <- order(mu_mat[, pivot_index])
    out$pi    <- out$pi[order_idx]
    out$mu    <- out$mu[order_idx]
    out$sigma <- out$sigma[order_idx]
  }

  out
}


#' Initialization from an ordering vector
#'
#' @param X numeric matrix (n x d).
#' @param ordering_vec numeric vector of length n (can contain NA).
#' @param K number of mixture components.
#' @param assume_EEI logical; if TRUE, use shared diagonal covariance (EEI-like).
#' @param nugget nonnegative scalar added to diagonal variances/covariances.
#' @param discretization one of "equal", "quantile", "kmeans".
#' @param na_action how to handle NA in ordering_vec: "drop" or "error".
#' @return list(pi, mu, sigma, keep_idx, cluster_rank)
#' @export
make_init <- function(
    X, ordering_vec, K,
    assume_EEI = TRUE,
    nugget = 0,
    discretization = c("equal", "quantile", "kmeans"),
    na_action = c("drop", "error")
) {
  X <- as.matrix(X)
  discretization <- match.arg(discretization)
  na_action <- match.arg(na_action)

  if (K < 1) stop("K must be >= 1.")
  if (length(ordering_vec) != nrow(X)) stop("ordering_vec must have length nrow(X).")

  keep_idx <- seq_len(nrow(X))
  if (anyNA(ordering_vec)) {
    if (na_action == "error") stop("ordering_vec contains NA.")
    keep_idx <- which(!is.na(ordering_vec))
    ordering_vec <- ordering_vec[keep_idx]
    X <- X[keep_idx, , drop = FALSE]
  }

  n <- nrow(X)
  d <- ncol(X)
  if (n < K) stop("After NA handling, n must be >= K.")

  # ---- discretize ordering_vec into K bins / clusters ----
  if (discretization == "kmeans") {
    cl <- stats::kmeans(ordering_vec, centers = K)
    centers <- as.numeric(cl$centers)
    cluster_ids <- cl$cluster

  } else if (discretization == "quantile") {
    probs <- seq(0, 1, length.out = K + 1)
    brks <- as.numeric(stats::quantile(ordering_vec, probs = probs, na.rm = TRUE))
    brks <- unique(brks)
    if (length(brks) < 2) stop("ordering_vec has too few unique values for quantile discretization.")

    cluster_ids <- cut(ordering_vec, breaks = brks, include.lowest = TRUE, labels = FALSE)

    if (all(is.na(cluster_ids))) stop("Quantile discretization failed (all clusters NA).")
    if (anyNA(cluster_ids)) {
      warning("Some values were not assigned to any quantile bin; dropping them.")
      ok <- which(!is.na(cluster_ids))
      cluster_ids <- cluster_ids[ok]
      ordering_vec <- ordering_vec[ok]
      X <- X[ok, , drop = FALSE]
      n <- nrow(X)
    }

    K_eff <- max(cluster_ids, na.rm = TRUE)
    if (K_eff != K) {
      warning("Quantile discretization produced ", K_eff, " groups (requested K=", K, ").")
      K <- K_eff
    }

    centers <- tapply(ordering_vec, cluster_ids, mean)
    centers <- as.numeric(centers)

  } else { # "equal"
    brks <- seq(min(ordering_vec), max(ordering_vec), length.out = K + 1)
    cluster_ids <- cut(ordering_vec, breaks = brks, include.lowest = TRUE, labels = FALSE)
    if (all(is.na(cluster_ids))) stop("Equal discretization failed (all clusters NA).")

    centers <- tapply(ordering_vec, cluster_ids, mean)
    centers <- as.numeric(centers)
  }

  # reorder clusters from left to right by centers
  cluster_order <- order(centers)
  cluster_rank <- match(cluster_ids, cluster_order)

  # ---- mu ----
  mu_list <- lapply(seq_len(K), function(k) {
    idx <- which(cluster_rank == k)
    if (length(idx) == 0) rep(0, d) else colMeans(X[idx, , drop = FALSE])
  })

  # ---- sigma ----
  if (assume_EEI) {
    diffs <- do.call(rbind, lapply(seq_len(K), function(k) {
      idx <- which(cluster_rank == k)
      if (length(idx) == 0) return(NULL)
      sweep(X[idx, , drop = FALSE], 2, mu_list[[k]])
    }))

    if (is.null(diffs) || nrow(diffs) == 0) {
      diag_var <- rep(1, d)
    } else {
      diag_var <- colMeans(diffs^2) + nugget
    }

    Sigma_shared <- diag(diag_var, d)
    sigma_list <- replicate(K, Sigma_shared, simplify = FALSE)

  } else {
    sigma_list <- lapply(seq_len(K), function(k) {
      idx <- which(cluster_rank == k)
      if (length(idx) <= 1) diag(rep(1, d))
      else stats::cov(X[idx, , drop = FALSE]) + nugget * diag(d)
    })
  }

  Nk <- tabulate(cluster_rank, nbins = K)
  pi_vec <- Nk / n

  list(
    pi = pi_vec,
    mu = mu_list,
    sigma = sigma_list,
    keep_idx = keep_idx,
    cluster_rank = cluster_rank
  )
}

# =========================
# Ordering methods
# =========================

#' Fiedler ordering from a kNN graph
#'
#' Builds a kNN graph on rows of \code{X}, computes a 1D ordering coordinate using
#' the Fiedler direction of a normalized graph operator, and returns an ordering
#' \code{t} scaled to \eqn{[0,1]}.
#'
#' @param X Numeric matrix \eqn{n \times D} (rows are observations).
#' @param k Number of nearest neighbors.
#' @param weight Similarity type: \code{"rbf"}, \code{"inv"}, or \code{"binary"}.
#' @param sigma Bandwidth for \code{weight="rbf"}; if \code{NULL}, uses median kNN distance.
#' @param keep Component handling: \code{"giant"} keeps the largest connected component;
#'   \code{"all"} uses all nodes (may be unstable if disconnected).
#' @param return_full If \code{TRUE} (default), return \code{t} of length \eqn{n}
#'   with \code{NA} for nodes excluded by \code{keep="giant"}.
#'   If \code{FALSE}, return \code{t} only for the kept nodes.
#' @return A list with \code{t}, \code{keep_idx}, and \code{n_components}.
#' @export
fiedler_ordering <- function(
    X,
    k = 15,
    weight = c("rbf", "inv", "binary"),
    sigma = NULL,
    keep = c("giant", "all"),
    return_full = TRUE
) {
  X <- as.matrix(X)
  n0 <- nrow(X)
  if (n0 < 3) stop("Need n >= 3.")
  if (k < 2 || k >= n0) stop("Require 2 <= k < n.")

  weight <- match.arg(weight)
  keep <- match.arg(keep)

  nn <- RANN::nn2(data = X, query = X, k = k + 1, treetype = "kd")
  idx  <- nn$nn.idx[, -1, drop = FALSE]
  dist <- nn$nn.dists[, -1, drop = FALSE]

  i <- rep(seq_len(n0), each = k)
  j <- as.vector(t(idx))
  d <- as.vector(t(dist))

  w <- switch(
    weight,
    binary = rep(1, length(d)),
    inv    = 1 / pmax(d, 1e-12),
    rbf    = {
      if (is.null(sigma)) sigma <- stats::median(d)
      exp(-(d^2) / (2 * sigma^2))
    }
  )

  W <- Matrix::sparseMatrix(i = i, j = j, x = w, dims = c(n0, n0))
  W <- pmax(W, Matrix::t(W))  # symmetrize (max)

  g <- igraph::graph_from_adjacency_matrix(W > 0, mode = "undirected", diag = FALSE)
  comp <- igraph::components(g)
  keep_idx <- seq_len(n0)

  if (keep == "giant" && comp$no > 1) {
    giant <- which.max(comp$csize)
    keep_idx <- which(comp$membership == giant)
    Wk <- W[keep_idx, keep_idx, drop = FALSE]
  } else {
    Wk <- W
  }

  if (nrow(Wk) < 3) {
    if (return_full) {
      return(list(t = rep(NA_real_, n0), keep_idx = keep_idx, n_components = comp$no))
    }
    stop("Largest connected component has < 3 nodes; increase k or filter outliers.")
  }

  deg <- Matrix::rowSums(Wk)
  D_inv_sqrt <- Matrix::Diagonal(nrow(Wk), x = 1 / sqrt(pmax(deg, 1e-12)))
  S <- D_inv_sqrt %*% Wk %*% D_inv_sqrt
  S <- (S + Matrix::t(S)) / 2

  e <- RSpectra::eigs_sym(
    S, k = 2, which = "LA",
    opts = list(tol = 1e-10, maxitr = 20000, retvec = TRUE)
  )

  u2 <- e$vectors[, 2]
  t_kept <- .scale01(u2)

  if (return_full) {
    t_full <- rep(NA_real_, n0)
    t_full[keep_idx] <- t_kept
    t_out <- t_full
  } else {
    t_out <- t_kept
  }

  # if comp$no is larger than 1, warn the user
  if (comp$no > 1) {
    warning("The kNN graph has ", comp$no, " connected components; ",
            "consider increasing k or filtering outliers.")
  }

  list(
    t = t_out,
    keep_idx = keep_idx,
    n_components = comp$no
  )
}


#' t-SNE ordering
#'
#' @param X numeric matrix (n x D).
#' @param tSNE_dims integer; embedding dimension passed to Rtsne.
#' @param component which embedding coordinate to use as ordering (default 1).
#' @param perplexity t-SNE perplexity; must satisfy < (n-1)/3.
#' @param max_iter max iterations.
#' @param seed optional integer for reproducibility.
#' @param scale01 logical; if TRUE, rescale returned t to [0,1].
#' @param orient_by_pc1 logical; if TRUE, flip sign to align with PC1.
#' @return list(t, keep_idx)
#' @export
tSNE_ordering <- function(
    X,
    tSNE_dims = 1,
    component = 1,
    perplexity = 10,
    max_iter = 500,
    seed = NULL,
    scale01 = TRUE,
    orient_by_pc1 = TRUE
) {
  X <- as.matrix(X)
  n <- nrow(X)

  if (component < 1 || component > tSNE_dims) {
    stop("component must be between 1 and tSNE_dims.")
  }
  if (perplexity >= (n - 1) / 3) {
    stop("perplexity must be < (n-1)/3 for Rtsne (reduce perplexity or increase n).")
  }
  if (!is.null(seed)) set.seed(seed)

  tsne_result <- Rtsne::Rtsne(
    X = X, dims = tSNE_dims, perplexity = perplexity,
    verbose = FALSE, max_iter = max_iter
  )

  t <- tsne_result$Y[, component]

  if (orient_by_pc1) t <- .orient_by_ref(t, .pc1(X))
  if (scale01) t <- .scale01(t)

  list(t = t, keep_idx = seq_len(n))
}


#' PCA ordering (PC1)
#'
#' @param X numeric matrix (n x D).
#' @param center logical.
#' @param scale logical.
#' @param scale01 logical; if TRUE, rescale returned t to [0,1].
#' @return list(t, keep_idx)
#' @export
PCA_ordering <- function(
    X,
    center = TRUE,
    scale = FALSE,
    scale01 = TRUE
) {
  X <- as.matrix(X)
  t <- .pc1(X, center = center, scale. = scale)
  if (scale01) t <- .scale01(t)
  list(t = t, keep_idx = seq_len(nrow(X)))
}


#' Principal curve ordering
#'
#' @param X numeric matrix (n x D).
#' @param smoother a function (recommended) or a name like "smooth_spline".
#' @param thresh,maxit,stretch,approx_points passed to princurve::principal_curve.
#' @param scale01 logical; if TRUE, rescale returned t to [0,1].
#' @param orient_by_pc1 logical; if TRUE, flip sign to align with PC1.
#' @return list(t, keep_idx, fit)
#' @export
pcurve_ordering <- function(
    X,
    smoother = c("smooth_spline", "lowess", "periodic_lowess"),
    thresh = 0.001,
    maxit = 10,
    stretch = 2,
    approx_points = FALSE,
    scale01 = TRUE,
    orient_by_pc1 = TRUE
) {
  X <- as.matrix(X)
  smoother <- match.arg(smoother)
  pc_fit <- princurve::principal_curve(
    X,
    smoother = smoother,
    thresh = thresh,
    stretch = stretch,
    approx_points = approx_points,
    maxit = maxit
  )

  t <- pc_fit$lambda
  if (orient_by_pc1) t <- .orient_by_ref(t, .pc1(X))
  if (scale01) t <- .scale01(t)

  list(t = t, keep_idx = seq_len(nrow(X)), fit = pc_fit)
}

# =========================
# Wrapper
# =========================

#' Initialize GMM parameters using an ordering method
#'
#' @param X numeric matrix (n x d).
#' @param K number of mixture components.
#' @param method ordering method.
#' @param discretization discretization method passed to make_init().
#' @param ... forwarded to the ordering method function.
#' @return list(pi, mu, sigma, keep_idx, cluster_rank, ordering)
#' @export
initialize_ordering <- function(
    X, K,
    method = c("PCA", "fiedler", "pcurve", "tSNE", "random"),
    discretization = c("equal", "quantile", "kmeans"),
    ...
) {
  method <- match.arg(method)
  discretization <- match.arg(discretization)

  X <- as.matrix(X)

  if (method == "random") {
    init <- make_default_init(X = X, K = K, ordering = TRUE)
    init$keep_idx <- seq_len(nrow(X))
    init$ordering <- list(method = method, t = NA_real_, keep_idx = init$keep_idx)
    return(init)
  }

  ordering_result <- switch(
    method,
    PCA     = PCA_ordering(X, ...),
    tSNE    = tSNE_ordering(X, ...),
    pcurve  = pcurve_ordering(X, ...),
    fiedler = fiedler_ordering(X, ...)
  )

  init <- make_init(
    X = X,
    ordering_vec = ordering_result$t,
    K = K,
    discretization = discretization,
    na_action = "drop"
  )

  init$ordering <- c(ordering_result, list(method = method))
  init
}

#' Initialize SmoothEM (single- or multi-scale)
#'
#' @description
#' Creates a \code{smooth_em} object using a chosen initialization method and then runs
#' SmoothEM for a specified number of iterations.
#'
#' The initialization is always "warm-started" by running exactly one outer EM iteration
#' (via \code{EM_algorithm(max_iter = 1)} for non-multi-scale; or \code{progressive_smoothEM(max_iter = 1)}
#' per stage for multi-scale). If \code{num_iter > 1}, the remaining iterations are completed
#' by \code{do_smoothEM()}, which enables features such as adaptive penalty updates.
#'
#' When \code{adaptive = TRUE}, an initial \eqn{\lambda} is estimated ("profile-style") from
#' the current mean vector \eqn{u} (stacked \code{mu}) and the base precision \code{Q_base}:
#' \deqn{\lambda^\*(u) = r / (u^\top Q_\mathrm{base} u), \quad r = p - \mathrm{rank\_deficiency}}
#' where \eqn{p} is the dimension of \code{Q_base}. The estimated value is clipped to
#' \code{[lambda_min, lambda_max]} and used as the starting \eqn{\lambda}.
#'
#' @param X Numeric matrix \code{(n x d)} of observations.
#' @param method Initialization method. If not \code{"multi_scale"}, it is passed to
#'   \code{initialize_ordering()}. If \code{"multi_scale"}, a progressive coarse-to-fine
#'   initialization is used via \code{progressive_smoothEM()}.
#' @param rw_q Integer random-walk order \code{q} for the separable RW penalty (default 2).
#' @param lambda Nonnegative penalty strength. For \code{method="multi_scale"}, this is
#'   interpreted as \code{lambda_final}.
#' @param relative_lambda Logical; if TRUE, rescales the precision used for fitting/evaluation
#'   by the current marginal variances (intended for EEI-like shared covariance across clusters).
#' @param K Integer number of grid points (used only when \code{method != "multi_scale"}).
#'   If NULL, defaults to \code{min(50, floor(nrow(X)/5))} with minimum 2.
#' @param m_max Integer finest exponent for multi-scale grid; final grid size is
#'   \code{K_final = 2^m_max + 1}.
#' @param num_iter Integer >= 1; total number of SmoothEM outer iterations to run.
#'   Exactly one iteration is run in the warm-start step, and the remaining
#'   \code{num_iter - 1} iterations (if any) are run by \code{do_smoothEM()}.
#' @param modelName Covariance model passed to \code{EM_algorithm()} / \code{MSTEP()}.
#'   Common choices include \code{"VVV"}, \code{"VII"}, \code{"EII"}, \code{"EEI"}.
#' @param ridge Nonnegative ridge added in constructing the RW precision.
#' @param nugget Nonnegative diagonal jitter added to covariance estimates during fitting.
#' @param eigen_tol Optional tolerance passed to generalized log-determinant routines used
#'   in objective/ELBO evaluation.
#' @param keep_history Logical; if TRUE and \code{method=="multi_scale"}, attach the full
#'   progressive initialization result under \code{meta$init$details$progressive}.
#' @param include.data Logical; if TRUE, store the data matrix in the returned object.
#' @param adaptive Logical; if TRUE, estimate a starting \eqn{\lambda} from the initial
#'   mean vector and enable adaptive updates during continuation via \code{do_smoothEM()}.
#' @param lambda_min,lambda_max Positive bounds used to clip the estimated/updated \eqn{\lambda}
#'   when \code{adaptive=TRUE}.
#' @param ... Extra arguments passed to \code{initialize_ordering()} (when
#'   \code{method != "multi_scale"}), or to \code{progressive_smoothEM()} (when
#'   \code{method == "multi_scale"}).
#'
#' @return A \code{smooth_em} object with fitted parameters, responsibilities, and traces.
#'   Initialization provenance is stored in \code{meta$init}, including \code{lambda_init_est}
#'   and \code{lambda_init_source} when \code{adaptive=TRUE}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rnorm(200 * 2), 200, 2)
#'
#' # Fixed lambda
#' fit0 <- initialize_smoothEM(X, method = "PCA", num_iter = 10, adaptive = FALSE, lambda = 10)
#'
#' # Adaptive lambda (estimate start + adapt during continuation)
#' fit1 <- initialize_smoothEM(X, method = "tSNE", num_iter = 10, adaptive = TRUE, lambda = 10)
#'
#' # Multi-scale initialization
#' fit2 <- initialize_smoothEM(X, method = "multi_scale", m_max = 6, num_iter = 5, adaptive = TRUE)
#' }
initialize_smoothEM <- function(
    X,
    method = c("tSNE", "PCA", "random", "multi_scale", "fiedler"),
    rw_q = 2,
    lambda = 1,
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
    lambda_min = 1e-8,
    lambda_max = 1e8,
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

  lambda_min <- as.numeric(lambda_min)
  lambda_max <- as.numeric(lambda_max)
  if (!is.finite(lambda_min) || !is.finite(lambda_max) ||
      lambda_min <= 0 || lambda_max <= 0 || lambda_min > lambda_max) {
    stop("lambda_min/lambda_max must be positive finite numbers with lambda_min <= lambda_max.")
  }

  # rank deficiency for RW(q) separable prior: q per coordinate
  rank_def <- as.integer(rw_q * d)

  # -------------------------
  # helpers for lambda init
  # -------------------------
  stack_U_from_params <- function(params) {
    if (is.list(params$mu)) {
      unlist(lapply(params$mu, function(m) as.numeric(m)), use.names = FALSE)
    } else {
      as.numeric(params$mu)
    }
  }

  compute_Q_eval_base <- function(params, Q_base, relative_lambda) {
    # returns base precision for evaluation, WITHOUT multiplying lambda
    Q_eval <- Q_base
    if (is.null(Q_eval)) return(NULL)

    if (isTRUE(relative_lambda)) {
      # relative_lambda requires identical covariances across clusters (EEI-like)
      same <- TRUE
      if (!is.null(params$sigma) && length(params$sigma) > 1) {
        same <- all(vapply(params$sigma, function(S) isTRUE(all.equal(S, params$sigma[[1]])), logical(1)))
      }
      if (!same) {
        # For initialization, degrade gracefully rather than stop.
        # We'll estimate lambda without relative scaling.
        relative_lambda <- FALSE
      } else {
        K <- length(params$pi)
        sigma_vec <- diag(params$sigma[[1]])
        scale_vec <- rep(1 / sqrt(pmax(sigma_vec, 1e-12)), times = K)
        Sscale <- Matrix::Diagonal(x = scale_vec)
        Q_eval <- Sscale %*% Q_eval %*% Sscale
      }
    }

    Q_eval
  }

  estimate_lambda_star <- function(params, Q_base, relative_lambda, rank_deficiency, eps_quad = 1e-12) {
    Qb <- compute_Q_eval_base(params, Q_base, relative_lambda)
    if (is.null(Qb)) return(NA_real_)

    U <- stack_U_from_params(params)
    quad <- as.numeric(crossprod(U, Qb %*% U))

    p <- ncol(Qb)
    r <- as.integer(p - (rank_deficiency %||% 0L))
    r <- max(r, 1L)

    r / pmax(quad, eps_quad)
  }

  clamp_lambda <- function(lam) {
    lam <- as.numeric(lam)
    if (!is.finite(lam)) return(NA_real_)
    min(max(lam, lambda_min), lambda_max)
  }

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

    # remember intent
    obj$control <- obj$control %||% list()
    obj$control$adapt_lambda <- isTRUE(adaptive)
    obj$control$lambda_min <- lambda_min
    obj$control$lambda_max <- lambda_max

    # warm-start lambda_trace if your do_smoothEM supports it
    if (is.null(obj$lambda_trace)) obj$lambda_trace <- numeric(0)
    if (length(obj$lambda_trace) == 0L && length(obj$elbo_trace %||% numeric(0)) > 0L) {
      obj$lambda_trace <- rep(obj$prior$lambda %||% lambda, length(obj$elbo_trace))
    }

    # Continue with do_smoothEM for remaining iterations
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
      K <- min(50L, floor(n / 5))
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

    # ---- UPDATED: estimate lambda before EM warm-start (if adaptive) ----
    lambda_use <- as.numeric(lambda)
    lambda_init_est <- NA_real_
    lambda_init_source <- "user"

    if (isTRUE(adaptive)) {
      lambda_init_est <- estimate_lambda_star(
        params = init_params,
        Q_base = Q_base,
        relative_lambda = relative_lambda,
        rank_deficiency = rank_def
      )
      lambda_init_est <- clamp_lambda(lambda_init_est)

      if (is.finite(lambda_init_est)) {
        lambda_use <- lambda_init_est
        lambda_init_source <- "estimated"
      } else {
        warning("Estimated initial lambda is not finite; using provided lambda.")
      }
    }

    Q_prior <- lambda_use * Q_base

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
      lambda = lambda_use,
      meta_init = list(
        method  = method,
        details = list(
          K = K,
          lambda_init_est = lambda_init_est,
          lambda_init_source = lambda_init_source
        )
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
    max_iter = 1L,              # only 1 iter per stage
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

  # ---- UPDATED: for multi-scale, re-estimate lambda from fit_final$params as the continuation start ----
  lambda_use <- as.numeric(lambda)
  lambda_init_est <- NA_real_
  lambda_init_source <- "user"

  if (isTRUE(adaptive) && !is.null(fit_final$params)) {
    lambda_init_est <- estimate_lambda_star(
      params = fit_final$params,
      Q_base = Q_base_final,
      relative_lambda = relative_lambda,
      rank_deficiency = rank_def
    )
    lambda_init_est <- clamp_lambda(lambda_init_est)

    if (is.finite(lambda_init_est)) {
      lambda_use <- lambda_init_est
      lambda_init_source <- "estimated"
    } else {
      warning("Estimated initial lambda (multi_scale) is not finite; using provided lambda.")
    }
  }

  obj <- finalize_and_continue(
    fit = fit_final,
    Q_base = Q_base_final,
    lambda = lambda_use,
    meta_init = list(
      method  = "multi_scale",
      details = list(
        m_max = m_max,
        K_final = K_final,
        lambda_init_est = lambda_init_est,
        lambda_init_source = lambda_init_source
      )
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
#' Each method is fit for \code{num_iter} total SmoothEM iterations. Internally,
#' \code{initialize_smoothEM()} warm-starts with exactly one EM iteration, and if
#' \code{num_iter > 1}, continues via \code{do_smoothEM()} (where adaptive lambda updates
#' may be enabled).
#'
#' For compatibility with \code{method="multi_scale"}, if \code{K} is not provided,
#' then non-multi-scale methods default to \code{K = 2^m_max + 1}.
#'
#' @param X Numeric matrix (n x d).
#' @param methods Character vector of methods to try. Defaults to
#'   \code{c("PCA","tSNE","random","fiedler","multi_scale")}.
#' @param num_iter Integer >= 1. Total number of SmoothEM outer iterations to run for each method.
#' @param num_cores Integer >= 1. Number of cores for parallel execution.
#' @param m_max Integer >= 1. Used by \code{multi_scale} and to set default \code{K}.
#' @param K Optional integer grid size for non-multi-scale methods. If NULL, uses \code{2^m_max+1}.
#' @param seed Optional base seed for reproducibility. If provided, each method gets a deterministic
#'   derived seed.
#' @param adaptive Logical; if TRUE, enables adaptive lambda behavior in \code{initialize_smoothEM()}.
#' @param lambda_min,lambda_max Positive bounds used to clip lambda when \code{adaptive=TRUE}.
#' @param quiet Logical; suppress messages from workers.
#' @param ... Extra args passed to \code{initialize_smoothEM()} (e.g. \code{rw_q}, \code{lambda},
#'   \code{relative_lambda}, \code{modelName}, etc.).
#'
#' @return A named list of \code{smooth_em} objects (or \code{NULL} for failed fits),
#' with attributes:
#' \itemize{
#'   \item \code{summary}: a data.frame summarizing last ELBO / last objective / last lambda for each method.
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
    adaptive = TRUE,
    lambda_min = 1e-8,
    lambda_max = 1e8,
    quiet = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  X <- as.matrix(X)
  methods <- unique(as.character(methods))

  num_iter  <- as.integer(num_iter)
  num_cores <- as.integer(num_cores)
  m_max     <- as.integer(m_max)

  if (length(num_iter) != 1L || is.na(num_iter) || num_iter < 1L) stop("num_iter must be integer >= 1.")
  if (length(num_cores) != 1L || is.na(num_cores) || num_cores < 1L) stop("num_cores must be integer >= 1.")
  if (length(m_max) != 1L || is.na(m_max) || m_max < 1L) stop("m_max must be integer >= 1.")

  adaptive <- isTRUE(adaptive)
  lambda_min <- as.numeric(lambda_min)
  lambda_max <- as.numeric(lambda_max)
  if (!is.finite(lambda_min) || lambda_min <= 0) stop("lambda_min must be a positive finite number.")
  if (!is.finite(lambda_max) || lambda_max <= 0) stop("lambda_max must be a positive finite number.")
  if (lambda_max < lambda_min) stop("lambda_max must be >= lambda_min.")

  # capture ... and avoid duplicates
  dots <- list(...)
  dots[c("adaptive", "lambda_min", "lambda_max")] <- NULL

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
          num_iter = num_iter,
          adaptive = adaptive,
          lambda_min = lambda_min,
          lambda_max = lambda_max
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
          adaptive = adaptive,
          lambda_min = lambda_min,
          lambda_max = lambda_max
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

    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(smoothEMr))
      NULL
    })

    parallel::clusterExport(
      cl,
      varlist = c(
        "X", "K_default", "m_max", "num_iter",
        "method_seeds", "adaptive", "lambda_min", "lambda_max",
        "quiet", "dots", "worker_one"
      ),
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
    if (is.null(obj) || is.null(obj$prior) || is.null(obj$prior$lambda)) return(NA_real_)
    as.numeric(obj$prior$lambda)
  }
  get_lambda_trace_last <- function(obj) {
    if (is.null(obj) || is.null(obj$lambda_trace) || length(obj$lambda_trace) == 0) return(NA_real_)
    tail(obj$lambda_trace, 1)
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
    lambda_last = vapply(results, get_lambda_last, numeric(1)),
    lambda_trace_last = vapply(results, get_lambda_trace_last, numeric(1)),
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
#' If \code{plot=TRUE}, plots all ELBO traces (different colors) and optionally overlays
#' the penalized objective traces in a second panel.
#'
#' @param X Numeric matrix (n x d).
#' @param methods Methods to try (passed to \code{parallel_initial()}).
#' @param num_iter Total number of iterations for each fit.
#' @param num_cores Number of cores.
#' @param m_max Used for multi_scale and default K for others.
#' @param K Optional K for non-multi-scale methods.
#' @param adaptive Logical; passed to \code{parallel_initial()} / \code{initialize_smoothEM()}.
#' @param lambda_min,lambda_max Positive bounds for lambda when \code{adaptive=TRUE}.
#' @param plot Logical; if TRUE, plot traces.
#' @param two_panel Logical; if TRUE, show ELBO and objective in 2 panels.
#' @param seed Optional base seed.
#' @param quiet Logical.
#' @param ... Passed to \code{initialize_smoothEM()} via \code{parallel_initial()}.
#'
#' @return A \code{smooth_em} object (best by last ELBO). The returned object gains:
#' \itemize{
#'   \item \code{$meta$initial_search}: list with the full fits, summary table, and options.
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
    lambda_min = 1e-8,
    lambda_max = 1e8,
    plot = FALSE,
    two_panel = FALSE,
    seed = NULL,
    quiet = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  fits <- parallel_initial(
    X = X,
    methods = methods,
    num_iter = num_iter,
    num_cores = num_cores,
    m_max = m_max,
    K = K,
    seed = seed,
    adaptive = adaptive,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
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
    summary = sum_df,
    fits = fits,
    options = list(
      methods = methods,
      num_iter = num_iter,
      num_cores = num_cores,
      m_max = m_max,
      K = K,
      adaptive = adaptive,
      lambda_min = lambda_min,
      lambda_max = lambda_max,
      seed = seed
    )
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
             lty = ifelse(methods_use == best_method, 1, 2),
             bty = "n")
      invisible(NULL)
    }

    if (two_panel) {
      par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
    }

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



