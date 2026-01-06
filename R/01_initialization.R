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

