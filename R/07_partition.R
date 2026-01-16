#' Compute per-coordinate contributions to the collapsed objective C for csmooth_em
#'
#' @description
#' Returns a length-d vector C_j such that
#'   C_total = sum_j C_j + entropy(Gamma) + sum_k Nk * log(pi_k),
#' where C_total corresponds to the collapsed objective recorded as ml_trace
#' (up to constants consistent with include_constant).
#'
#' @param X n-by-d data matrix
#' @param Gamma n-by-K responsibilities
#' @param params list(pi, mu, sigma2)
#' @param Q_K base RW precision (K x K)
#' @param lambda_vec length-d vector
#' @param modelName "homoskedastic" or "heteroskedastic"
#' @param relative_lambda logical
#' @param eigen_tol numeric; passed to generalized_logdet
#' @param rw_q integer; rank deficiency of RW precision
#' @param include_constant logical; include +(K/2)log(2pi) per coordinate if TRUE
#'
#' @return a list with
#'   - C_coord: length-d vector of per-coordinate collapsed contributions
#'   - global_logpi: scalar sum_k Nk log(pi_k)
#'   - global_entropy: scalar entropy term -sum Gamma log Gamma
#'   - logdetH_coord: length-d vector log|A_j| (for diagnostics)
#' @keywords internal
compute_C_by_coord_csmooth <- function(X, Gamma, params, Q_K, lambda_vec,
                                       modelName, relative_lambda,
                                       eigen_tol = NULL, rw_q = 0L,
                                       include_constant = TRUE) {
  X <- as.matrix(X)
  Gamma <- as.matrix(Gamma)

  n <- nrow(X)
  d <- ncol(X)
  K <- ncol(Gamma)

  pi_k <- pmax(as.numeric(params$pi), .Machine$double.eps)
  Nk <- colSums(Gamma)
  Nk[Nk < 1e-8] <- 1e-8

  # sufficient stats
  GX  <- t(X)   %*% Gamma   # d x K
  GX2 <- t(X^2) %*% Gamma   # d x K

  # mu matrix: d x K
  mu_mat <- do.call(cbind, params$mu)
  if (!all(dim(mu_mat) == c(d, K))) stop("params$mu has incompatible dimensions.")

  # sigma2
  if (modelName == "homoskedastic") {
    sigma2 <- pmax(as.numeric(params$sigma2), .Machine$double.eps)  # length d
  } else {
    sigma2 <- pmax(as.matrix(params$sigma2), .Machine$double.eps)   # d x K
    if (!all(dim(sigma2) == c(d, K))) stop("params$sigma2 has incompatible dimensions.")
  }

  # global terms (not feature-specific)
  global_logpi <- sum(Nk * log(pi_k))
  global_entropy <- -sum(Gamma * log(Gamma + 1e-300))

  # outputs
  C_coord <- numeric(d)
  logdetH_coord <- numeric(d)

  const_j <- if (isTRUE(include_constant)) 0.5 * K * log(2 * base::pi) else 0

  for (j in seq_len(d)) {
    lam <- lambda_vec[j]
    if (!is.finite(lam) || lam < 0) lam <- 0

    # feature-specific sigma2 across K
    if (modelName == "homoskedastic") {
      s2jk <- rep(sigma2[j], K)
    } else {
      s2jk <- as.numeric(sigma2[j, ])
    }
    s2jk <- pmax(s2jk, .Machine$double.eps)

    # data term for feature j at current mu (MAP)
    # Q_j = sum_k [-0.5 Nk log(2pi*s2jk) - 0.5 * SSE_jk / s2jk]
    SSE_jk <- as.numeric(GX2[j, ]) - 2 * mu_mat[j, ] * as.numeric(GX[j, ]) + (mu_mat[j, ]^2) * Nk
    Q_data_j <- -0.5 * sum(Nk * log(2 * base::pi * s2jk)) - 0.5 * sum(SSE_jk / s2jk)

    # prior term for feature j at current mu
    prior_j <- 0
    if (lam > 0 && !is.null(Q_K)) {
      Qb <- .compute_Qbase_j(pi_k, sigma2, Q_K, relative_lambda, modelName, j)
      Qj <- lam * Qb

      mu_j <- as.numeric(mu_mat[j, ])
      prior_j <- -0.5 * as.numeric(crossprod(mu_j, Qj %*% mu_j)) +
        0.5 * generalized_logdet(Q = Qj, eigen_tol = eigen_tol, rank_deficiency = rw_q)
    }

    # Hessian block A_j = D_j + lam * Qb, where D_j = diag(Nk / s2jk)
    # If lam == 0, A_j = D_j (still SPD)
    if (lam > 0 && !is.null(Q_K)) {
      Qb <- .compute_Qbase_j(pi_k, sigma2, Q_K, relative_lambda, modelName, j)
      Dj <- Nk / s2jk
      A  <- diag(as.numeric(Dj), K, K) + lam * Qb
    } else {
      Dj <- Nk / s2jk
      A  <- diag(as.numeric(Dj), K, K)
    }
    A <- 0.5 * (A + t(A))
    L <- chol(A)
    logdetA <- 2 * sum(log(diag(L)))

    logdetH_coord[j] <- logdetA

    # collapsed contribution
    C_coord[j] <- Q_data_j + prior_j + const_j - 0.5 * logdetA
  }

  list(
    C_coord = C_coord,
    global_logpi = as.numeric(global_logpi),
    global_entropy = as.numeric(global_entropy),
    logdetH_coord = logdetH_coord
  )
}


# ============================================================
# Feature scoring + partition utilities (build on C_coord)
# ============================================================

#' Score features under a fitted csmooth_em model via per-coordinate collapsed contributions
#'
#' @param fit A csmooth_em object.
#' @param X Optional n-by-d data matrix. If NULL, uses fit$data.
#' @param include_constant Logical; passed to compute_C_by_coord_csmooth().
#'
#' @return A list with
#'   - C_coord: length-d vector
#'   - global_logpi, global_entropy: scalars
#'   - logdetH_coord: length-d vector
#' @export
score_features_onefit <- function(fit, X = NULL, include_constant = TRUE) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  if (!inherits(fit, "csmooth_em")) stop("fit must be a 'csmooth_em' object.")

  if (is.null(X)) {
    X <- fit$data
    if (is.null(X)) stop("X must be provided or stored in fit$data.")
  }
  X <- as.matrix(X)

  Gamma <- fit$gamma
  if (is.null(Gamma)) stop("fit$gamma is NULL. Run do_csmoothEM(..., record=TRUE) so gamma is stored.")

  compute_C_by_coord_csmooth(
    X = X,
    Gamma = Gamma,
    params = fit$params,
    Q_K = fit$prior$Q_K,
    lambda_vec = fit$prior$lambda_vec,
    modelName = fit$control$modelName %||% "homoskedastic",
    relative_lambda = isTRUE(fit$control$relative_lambda),
    eigen_tol = fit$control$eigen_tol,
    rw_q = as.integer(fit$prior$rw_q %||% 0L),
    include_constant = include_constant
  )
}

#' Partition features by comparing per-coordinate collapsed scores from two fits
#'
#' @param fitA,fitB csmooth_em objects (e.g., two different orderings / initializations).
#' @param X Optional data matrix. If NULL, uses fitA$data (and assumes same X for both fits).
#' @param delta Nonnegative margin. Assign to A only if (CA - CB) > delta; otherwise assign to B.
#' @param include_constant Logical; passed to score_features_onefit().
#'
#' @return A list with
#'   - assign: length-d character vector in {"A","B"}
#'   - score_diff: length-d vector (CA - CB)
#'   - CA, CB: length-d vectors of per-feature scores
#' @export
partition_features_twofits <- function(fitA, fitB, X = NULL, delta = 0, include_constant = TRUE) {
  if (!inherits(fitA, "csmooth_em") || !inherits(fitB, "csmooth_em")) {
    stop("fitA and fitB must both be 'csmooth_em' objects.")
  }
  if (!is.finite(delta) || delta < 0) stop("delta must be a nonnegative finite number.")

  if (is.null(X)) {
    X <- fitA$data
    if (is.null(X)) stop("X must be provided or stored in fitA$data.")
  }
  X <- as.matrix(X)

  resA <- score_features_onefit(fitA, X = X, include_constant = include_constant)
  resB <- score_features_onefit(fitB, X = X, include_constant = include_constant)

  CA <- as.numeric(resA$C_coord)
  CB <- as.numeric(resB$C_coord)
  if (length(CA) != length(CB)) stop("fitA and fitB produce different d; check inputs.")

  score_diff <- CA - CB
  assign <- ifelse(score_diff > delta, "A", "B")

  list(assign = assign, score_diff = score_diff, CA = CA, CB = CB)
}



#' Subset a csmooth_em object by features (columns)
#'
#' @description
#' Internal helper to subset a fitted \code{csmooth_em} object to a subset of features.
#' This keeps parameters aligned with a reduced feature set (mu, sigma2, lambda_vec, etc.).
#'
#' @param fit A \code{csmooth_em} object.
#' @param keep_cols Integer vector of feature indices to keep (1-based, in the current fit's feature space).
#'
#' @return A \code{csmooth_em} object restricted to the selected features.
#' @export
subset_csmooth_em_fit <- function(fit, keep_cols) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  if (!inherits(fit, "csmooth_em")) stop("fit must be a 'csmooth_em' object.")

  keep_cols <- as.integer(keep_cols)
  keep_cols <- keep_cols[!is.na(keep_cols)]
  if (length(keep_cols) == 0) stop("keep_cols is empty.")
  keep_cols <- sort(unique(keep_cols))

  # subset data
  if (!is.null(fit$data)) {
    fit$data <- as.matrix(fit$data)[, keep_cols, drop = FALSE]
  }

  # subset mu: list of length K, each element length d
  if (!is.null(fit$params$mu) && length(fit$params$mu) > 0) {
    fit$params$mu <- lapply(fit$params$mu, function(v) as.numeric(v)[keep_cols])
  }

  # subset sigma2
  if (!is.null(fit$params$sigma2)) {
    if (is.matrix(fit$params$sigma2)) {
      fit$params$sigma2 <- as.matrix(fit$params$sigma2)[keep_cols, , drop = FALSE]
    } else {
      fit$params$sigma2 <- as.numeric(fit$params$sigma2)[keep_cols]
    }
  }

  # subset lambda_vec
  if (!is.null(fit$prior$lambda_vec)) {
    fit$prior$lambda_vec <- as.numeric(fit$prior$lambda_vec)[keep_cols]
  }

  # subset lambda_trace if present
  if (!is.null(fit$lambda_trace) && length(fit$lambda_trace) > 0) {
    fit$lambda_trace <- lapply(fit$lambda_trace, function(lv) as.numeric(lv)[keep_cols])
  }

  # record which original columns kept (best-effort)
  fit$meta <- fit$meta %||% list()
  fit$meta$feature_subset <- list(keep_cols = keep_cols)

  fit
}


#' Greedy backward filtering of features for csmoothEM ordering inference
#'
#' @description
#' Implements a greedy "backward" feature filtering strategy to mitigate the case where
#' most features are noise and do not follow any latent ordering. The algorithm:
#' \enumerate{
#'   \item Fits csmoothEM on all features to obtain an initial ordering (via responsibilities \eqn{\Gamma}).
#'   \item Computes per-feature collapsed contributions \eqn{C_j} (via \code{compute_C_by_coord_csmooth}).
#'   \item Removes a batch of features with the smallest \eqn{C_j}.
#'   \item Refits csmoothEM on the remaining features for a few iterations.
#'   \item Repeats until a stopping rule is met.
#' }
#'
#' This procedure is intended as a preprocessing step before more ambitious tasks such as
#' feature partitioning across multiple orderings.
#'
#' @param X Numeric matrix \code{(n x d)}.
#' @param method Ordering method passed to \code{\link{initialize_csmoothEM}}. One of
#'   \code{"fiedler"}, \code{"PCA"}, \code{"tSNE"}, \code{"pcurve"}, \code{"random"}.
#' @param K Integer \eqn{\ge 2}. Number of mixture components.
#' @param modelName Either \code{"homoskedastic"} or \code{"heteroskedastic"}.
#' @param adaptive Adaptive mode passed to \code{\link{do_csmoothEM}} when refitting.
#'   Typically \code{"prior"} for speed (or \code{"ml"} if using collapsed-ML).
#' @param num_iter_init Integer \eqn{\ge 1}. Number of warm-start iterations for the initial fit.
#' @param num_iter_refit Integer \eqn{\ge 1}. Number of iterations for each refit after feature removal.
#' @param discretization Discretization method for initialization passed to \code{initialize_csmoothEM}.
#'   Recommended: \code{"quantile"} to avoid empty components.
#' @param batch Integer \eqn{\ge 1}. Number of lowest-scoring features (smallest \eqn{C_j}) removed per round.
#' @param min_keep Integer \eqn{\ge 1}. Minimum number of features to keep; stops if fewer would remain.
#' @param tau Optional numeric threshold. If provided, stops when \code{min(Cj) >= tau}.
#' @param max_rounds Integer \eqn{\ge 1}. Maximum number of greedy rounds.
#' @param verbose Logical; print a one-line summary each round.
#' @param ... Additional arguments passed to \code{initialize_csmoothEM} (e.g. ordering controls).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{keep_cols}: integer indices of retained features (w.r.t. the original X).
#'   \item \code{drop_cols}: integer indices of removed features (w.r.t. the original X).
#'   \item \code{fit}: final fitted \code{csmooth_em} object on the retained features.
#'   \item \code{history}: data.frame with per-round diagnostics (\code{C_total}, \code{min_Cj}, etc.).
#' }
#'
#' @export
greedy_backward_filter_csmooth <- function(
    X,
    method = c("fiedler", "PCA", "tSNE", "pcurve", "random"),
    K = 50,
    modelName = c("homoskedastic", "heteroskedastic"),
    adaptive = "prior",
    num_iter_init = 10,
    num_iter_refit = 5,
    discretization = c("equal", "quantile", "kmeans"),
    batch = 20,
    min_keep = 20,
    tau = NULL,
    max_rounds = 50,
    verbose = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  X <- as.matrix(X)
  d0 <- ncol(X)

  method <- match.arg(method)
  modelName <- match.arg(modelName)
  discretization <- match.arg(discretization)

  batch <- as.integer(batch)
  if (length(batch) != 1L || is.na(batch) || batch < 1L) stop("batch must be a single integer >= 1.")

  min_keep <- as.integer(min_keep)
  if (length(min_keep) != 1L || is.na(min_keep) || min_keep < 1L) stop("min_keep must be a single integer >= 1.")

  max_rounds <- as.integer(max_rounds)
  if (length(max_rounds) != 1L || is.na(max_rounds) || max_rounds < 1L) stop("max_rounds must be a single integer >= 1.")

  num_iter_init <- as.integer(num_iter_init)
  if (length(num_iter_init) != 1L || is.na(num_iter_init) || num_iter_init < 1L) stop("num_iter_init must be integer >= 1.")

  num_iter_refit <- as.integer(num_iter_refit)
  if (length(num_iter_refit) != 1L || is.na(num_iter_refit) || num_iter_refit < 1L) stop("num_iter_refit must be integer >= 1.")

  # initial fit on all features
  fit <- initialize_csmoothEM(
    X = X,
    method = method,
    K = K,
    modelName = modelName,
    adaptive = adaptive,
    num_iter = num_iter_init,
    discretization = discretization,
    ...
  )
  # ensure we have gamma + traces
  fit <- do_csmoothEM(fit, data = X, iter = 1, adaptive = NULL, record = TRUE, verbose = FALSE)

  keep <- seq_len(d0)
  history <- data.frame(
    round = integer(0),
    d = integer(0),
    C_total = numeric(0),
    min_Cj = numeric(0),
    q10_Cj = numeric(0),
    median_Cj = numeric(0),
    stringsAsFactors = FALSE
  )

  for (rr in seq_len(max_rounds)) {
    Gamma <- fit$gamma
    if (is.null(Gamma)) stop("fit$gamma is NULL; ensure do_csmoothEM stores gamma.")

    res <- compute_C_by_coord_csmooth(
      X = X[, keep, drop = FALSE],
      Gamma = Gamma,
      params = fit$params,
      Q_K = fit$prior$Q_K,
      lambda_vec = fit$prior$lambda_vec,
      modelName = fit$control$modelName,
      relative_lambda = isTRUE(fit$control$relative_lambda),
      eigen_tol = fit$control$eigen_tol,
      rw_q = as.integer(fit$prior$rw_q %||% 0L),
      include_constant = TRUE
    )

    Cj <- as.numeric(res$C_coord)
    C_total <- sum(Cj) + res$global_logpi + res$global_entropy
    minC <- min(Cj)
    q10 <- as.numeric(stats::quantile(Cj, 0.10, names = FALSE))
    med <- stats::median(Cj)

    history <- rbind(history, data.frame(
      round = rr, d = length(keep),
      C_total = C_total,
      min_Cj = minC,
      q10_Cj = q10,
      median_Cj = med
    ))

    if (verbose) {
      cat(sprintf(
        "[round %d] d=%d  C_total=%.3f  minC=%.3f  q10=%.3f  median=%.3f\n",
        rr, length(keep), C_total, minC, q10, med
      ))
    }

    # stopping rules
    if (!is.null(tau) && is.finite(tau) && minC >= tau) break
    if (length(keep) <= min_keep) break

    # drop batch smallest features
    ord <- order(Cj, decreasing = FALSE)
    drop_local <- ord[seq_len(min(batch, length(ord)))]
    keep_local <- setdiff(seq_along(keep), drop_local)

    if (length(keep_local) < min_keep) break

    # update keep indices in ORIGINAL feature space
    keep <- keep[keep_local]

    # subset fit to remaining features (in the CURRENT fit feature space)
    fit <- subset_csmooth_em_fit(fit, keep_cols = keep_local)

    # refit briefly on reduced feature set
    X_sub <- X[, keep, drop = FALSE]
    fit <- do_csmoothEM(fit, data = X_sub, iter = num_iter_refit,
                        adaptive = adaptive, record = TRUE, verbose = FALSE)
  }

  list(
    keep_cols = keep,
    drop_cols = setdiff(seq_len(d0), keep),
    fit = fit,
    history = history
  )
}




#' Append a coordinate to csmooth-style parameters (internal)
#'
#' @description
#' Appends a single coordinate (feature) to an existing csmooth-style parameter list.
#' This is used in greedy feature assignment to grow a partition without reinitializing
#' from scratch.
#'
#' @param params List with fields \code{pi}, \code{mu}, \code{sigma2} for the current partition.
#'   \code{mu} must be a list of length K, each element a numeric vector of length d_sub.
#'   \code{sigma2} must be a numeric vector of length d_sub (homoskedastic case).
#' @param one List describing the new coordinate. Must contain:
#'   \itemize{
#'     \item \code{mu_vec}: numeric vector of length K (component means for the new coordinate)
#'     \item \code{sigma2}: numeric scalar (variance for the new coordinate)
#'   }
#'
#' @return Updated \code{params} with the new coordinate appended.
#' @keywords internal
append_coord_to_params_csmooth <- function(params, one) {
  K <- length(params$pi)
  if (length(params$mu) != K) stop("params$mu must be a list of length K.")
  if (length(one$mu_vec) != K) stop("one$mu_vec must have length K.")

  params$mu <- lapply(seq_len(K), function(k) c(params$mu[[k]], one$mu_vec[k]))
  params$sigma2 <- c(as.numeric(params$sigma2), as.numeric(one$sigma2))
  params
}




#' Score a single feature given responsibilities (Gamma)
#'
#' @description
#' Computes an alignment score for a single feature \eqn{x_j} under a fixed responsibilities
#' matrix \eqn{\Gamma \in \mathbb{R}^{n\times K}}. The score is used for greedy feature
#' partitioning into multiple latent orderings.
#'
#' Two score modes are supported:
#' \describe{
#'   \item{\code{score_mode = "none"}}{
#'     Plug-in (ELBO/Q-like) score. The smoothing parameter \eqn{\lambda} is fixed at
#'     \code{lambda_init} (default 1). The score depends on the MAP component means
#'     \eqn{\hat\mu_{j\cdot}} via the weighted SSE.
#'   }
#'   \item{\code{score_mode = "ml"}}{
#'     Collapsed (marginal-like) score: plug-in score plus the Laplace curvature correction
#'     \eqn{+\frac{K}{2}\log(2\pi) - \frac12\log|A_j|}. In this mode, \eqn{\lambda} is optimized
#'     by 1D maximization over \eqn{\log \lambda \in [\log(\text{lambda_min}),\log(\text{lambda_max})]}.
#'   }
#' }
#'
#' The fitted 1D quantities \eqn{\hat\mu_{j\cdot}} and \eqn{\hat\sigma_j^2} are returned and can be
#' appended to a partition's csmooth-style parameters during greedy growth.
#'
#' @param xj Numeric vector of length \code{n} (one feature).
#' @param Gamma Numeric matrix \code{(n x K)} of responsibilities.
#' @param Q_K Numeric matrix \code{(K x K)}; base RW precision (lambda=1).
#' @param rw_q Integer \eqn{\ge 0}. Rank deficiency along K (RW order).
#' @param score_mode One of \code{"ml"} or \code{"none"}.
#' @param relative_lambda Logical; if TRUE use \eqn{Q_{base} = Q_K / \sigma_j^2} (homoskedastic scaling).
#' @param lambda_min,lambda_max Positive bounds for \eqn{\lambda} when \code{score_mode="ml"}.
#' @param nugget Nonnegative scalar added to \eqn{\sigma_j^2}.
#' @param optimize_lambda Logical or NULL. If NULL, defaults to TRUE for \code{"ml"} and FALSE for \code{"none"}.
#' @param lambda_init Positive scalar. Fixed \eqn{\lambda} when \code{score_mode="none"} (default 1),
#'   and initial value if \code{optimize_lambda=FALSE}.
#' @param max_sigma_iter Integer \eqn{\ge 1}. Number of (lambda -> mu -> sigma2) refresh steps.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{score}: scalar score for feature \eqn{j}.
#'   \item \code{lambda}: fitted (or fixed) \eqn{\lambda_j}.
#'   \item \code{sigma2}: fitted \eqn{\sigma_j^2}.
#'   \item \code{mu_vec}: numeric vector length K of \eqn{\hat\mu_{j\cdot}}.
#'   \item \code{logdetA}: scalar \eqn{\log|A_j|} (only meaningful for \code{"ml"}; still returned).
#' }
#'
#' @seealso \code{\link{score_one_coord_csmooth}}, \code{\link{forward_two_ordering_partition_csmooth}}
#' @export
score_feature_given_Gamma <- function(
    xj,
    Gamma,
    Q_K,
    rw_q = 0L,
    score_mode = c("ml", "none"),
    relative_lambda = TRUE,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    nugget = 0,
    optimize_lambda = NULL,
    lambda_init = 1,
    max_sigma_iter = 1L
) {
  score_mode <- match.arg(score_mode)

  # default semantics: ml -> optimize lambda; none -> fixed lambda_init (typically 1)
  if (is.null(optimize_lambda)) optimize_lambda <- (score_mode == "ml")

  xj <- as.numeric(xj)
  Gamma <- as.matrix(Gamma)

  n <- length(xj)
  if (nrow(Gamma) != n) stop("nrow(Gamma) must equal length(xj).")
  K <- ncol(Gamma)

  rw_q <- as.integer(rw_q)
  Nk <- colSums(Gamma)
  Nk[Nk < 1e-8] <- 1e-8

  # sufficient stats
  GX  <- as.numeric(crossprod(xj, Gamma))
  GX2 <- as.numeric(crossprod(xj^2, Gamma))

  # robust initial sigma2
  sigma2_j <- max(stats::var(xj) + nugget, .Machine$double.eps)

  # solve MAP mu for fixed (lambda, sigma2_j)
  solve_mu <- function(lambda, sigma2_j) {
    sigma2_j <- max(as.numeric(sigma2_j), .Machine$double.eps)

    Qb <- if (isTRUE(relative_lambda)) Q_K / sigma2_j else Q_K
    Qb <- 0.5 * (Qb + t(Qb))

    Dj <- Nk / sigma2_j
    if (any(!is.finite(Dj)) || any(Dj <= 0)) return(NULL)

    bj <- GX / sigma2_j
    if (any(!is.finite(bj))) return(NULL)

    A <- diag(as.numeric(Dj), K, K) + as.numeric(lambda) * Qb
    A <- 0.5 * (A + t(A))

    L <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(L)) {
      A2 <- A + 1e-10 * diag(K)
      L <- tryCatch(chol(A2), error = function(e) NULL)
      if (is.null(L)) return(NULL)
      A <- A2
    }

    y <- forwardsolve(t(L), bj)
    mu <- as.numeric(backsolve(L, y))

    logdetA <- 2 * sum(log(diag(L)))
    quad_bAinvb <- sum(bj * mu)

    list(mu = mu, logdetA = logdetA, quad = quad_bAinvb)
  }

  # likelihood contribution at MAP mu (depends on mu via SSE)
  lik_term <- function(mu, sigma2_j) {
    sigma2_j <- max(as.numeric(sigma2_j), .Machine$double.eps)
    SSE <- GX2 - 2 * mu * GX + (mu^2) * Nk
    -0.5 * sum(Nk) * log(2 * base::pi * sigma2_j) - 0.5 * sum(SSE / sigma2_j)
  }

  # prior term at MAP mu including pseudo-logdet normalization
  prior_at_mu <- function(mu, lambda, sigma2_j) {
    lambda <- as.numeric(lambda)
    if (!is.finite(lambda) || lambda <= 0) return(0)

    sigma2_j <- max(as.numeric(sigma2_j), .Machine$double.eps)
    Qb <- if (isTRUE(relative_lambda)) Q_K / sigma2_j else Q_K
    Qb <- 0.5 * (Qb + t(Qb))
    Qj <- lambda * Qb

    pen <- -0.5 * as.numeric(crossprod(mu, Qj %*% mu))
    ld  <- 0.5 * generalized_logdet(Q = Qj, eigen_tol = NULL, rank_deficiency = rw_q)
    as.numeric(pen + ld)
  }

  # objective as a function of log(lambda), for fixed sigma2_j
  obj_loglam <- function(loglam, sigma2_j) {
    lambda <- exp(loglam)
    s <- solve_mu(lambda, sigma2_j)
    if (is.null(s)) return(-Inf)

    mu <- s$mu
    base <- lik_term(mu, sigma2_j) + prior_at_mu(mu, lambda, sigma2_j)

    if (score_mode == "none") {
      as.numeric(base)
    } else {
      as.numeric(base + 0.5 * K * log(2 * base::pi) - 0.5 * s$logdetA)
    }
  }

  lambda_hat <- as.numeric(lambda_init)

  max_sigma_iter <- as.integer(max_sigma_iter)
  if (length(max_sigma_iter) != 1L || is.na(max_sigma_iter) || max_sigma_iter < 1L) {
    stop("max_sigma_iter must be a single integer >= 1.")
  }

  for (tt in seq_len(max_sigma_iter)) {
    if (isTRUE(optimize_lambda)) {
      opt <- optimize(
        f = function(loglam) obj_loglam(loglam, sigma2_j),
        interval = log(c(lambda_min, lambda_max)),
        maximum = TRUE
      )
      lambda_hat <- exp(opt$maximum)
    }

    s <- solve_mu(lambda_hat, sigma2_j)
    if (is.null(s)) {
      return(list(score = -Inf, lambda = lambda_hat, sigma2 = sigma2_j,
                  mu_vec = rep(NA_real_, K), logdetA = NA_real_))
    }

    mu_hat <- s$mu
    SSE <- GX2 - 2 * mu_hat * GX + (mu_hat^2) * Nk
    sigma2_j <- max(sum(SSE) / sum(Nk) + nugget, .Machine$double.eps)
  }

  s <- solve_mu(lambda_hat, sigma2_j)
  if (is.null(s)) {
    return(list(score = -Inf, lambda = lambda_hat, sigma2 = sigma2_j,
                mu_vec = rep(NA_real_, K), logdetA = NA_real_))
  }

  score <- obj_loglam(log(lambda_hat), sigma2_j)

  list(
    score = as.numeric(score),
    lambda = as.numeric(lambda_hat),
    sigma2 = as.numeric(sigma2_j),
    mu_vec = as.numeric(s$mu),
    logdetA = as.numeric(s$logdetA)
  )
}


#' Score one coordinate against fixed responsibilities (Gamma) for csmoothEM (internal)
#'
#' @description
#' Convenience wrapper used in greedy partitioning. Returns both the scalar score and a
#' 1D fitted object (\code{one}) that can be appended to csmooth-style parameters.
#'
#' @param X Numeric matrix \code{(n x d)}.
#' @param j Integer coordinate index in \code{1:d}.
#' @param Gamma Numeric matrix \code{(n x K)} of responsibilities.
#' @param Q_K Numeric matrix \code{(K x K)} base RW precision (lambda=1).
#' @param rw_q Integer \eqn{\ge 0}. Rank deficiency along K.
#' @param score_mode One of \code{"ml"} or \code{"none"}.
#' @param relative_lambda Logical.
#' @param lambda_min,lambda_max Bounds for lambda optimization (used when \code{score_mode="ml"}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{score}: scalar score for coordinate \code{j}.
#'   \item \code{one}: list with \code{mu_vec} (length K) and \code{sigma2} (scalar).
#' }
#' @keywords internal
#' @export
score_one_coord_csmooth <- function(
    X, j, Gamma, Q_K,
    rw_q = 2,
    score_mode = c("ml","none"),
    relative_lambda = TRUE,
    lambda_min = 1e-10,
    lambda_max = 1e10
) {
  score_mode <- match.arg(score_mode)
  j <- as.integer(j)
  if (length(j) != 1L || is.na(j) || j < 1L || j > ncol(X)) stop("j out of range.")

  out <- score_feature_given_Gamma(
    xj = X[, j],
    Gamma = Gamma,
    Q_K = Q_K,
    rw_q = rw_q,
    score_mode = score_mode,
    relative_lambda = relative_lambda,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    nugget = 0,
    optimize_lambda = NULL,   # default semantics: ml optimizes, none fixes
    lambda_init = 1,
    max_sigma_iter = 1L
  )

  list(
    score = out$score,
    one = list(mu_vec = out$mu_vec, sigma2 = out$sigma2)
  )
}




#' Forward greedy feature partition into two orderings (csmoothEM version)
#'
#' @description
#' Greedily partitions features (columns of \code{X}) into two groups, each associated with
#' its own latent ordering (represented by responsibilities \eqn{\Gamma} over \code{K} components).
#'
#' The algorithm mirrors the structure of \code{two_ordering_smoothEM_v2}:
#' \enumerate{
#'   \item Choose a seed feature \code{j1} (largest variance by default) and fit a 1D csmoothEM model
#'         on \code{X[, j1]} to obtain \eqn{\Gamma_1} and parameters \eqn{\theta_1}.
#'   \item Choose a second seed feature \code{j2} as the feature with the worst alignment score
#'         under \eqn{\Gamma_1}, then fit a 1D model on \code{X[, j2]} to obtain \eqn{\Gamma_2} and \eqn{\theta_2}.
#'   \item While unassigned features remain:
#'     \itemize{
#'       \item Score each remaining feature under \eqn{\Gamma_1} and \eqn{\Gamma_2}.
#'       \item Select the feature with the largest absolute score gap and assign it to the better ordering.
#'       \item Append the 1D fitted object to that ordering's parameters and update \eqn{\Gamma} by an E-step.
#'       \item Optionally run a short csmoothEM refinement on that ordering (controlled by \code{greedy_em_refine}).
#'     }
#' }
#'
#' Scoring:
#' \itemize{
#'   \item \code{score_mode="none"} uses a plug-in (ELBO/Q-like) score with fixed \eqn{\lambda=1}.
#'   \item \code{score_mode="ml"} uses a collapsed (marginal-like) score with \eqn{\lambda} optimized per feature.
#' }
#'
#' @param X Numeric matrix \code{(n x d)}.
#' @param K Integer \eqn{\ge 2}. Number of mixture components.
#' @param greedy_start_index Optional integer. First seed feature. If NULL, uses the largest-variance feature.
#' @param greedy_start_index2 Optional integer. Second seed feature. If NULL, chosen as worst under ordering 1.
#' @param score_mode One of \code{"ml"} or \code{"none"}.
#' @param rw_q Integer \eqn{\ge 0}. Rank deficiency along K (RW order).
#' @param relative_lambda Logical; whether relative-lambda scaling is used.
#' @param lambda_min,lambda_max Bounds for lambda optimization when \code{score_mode="ml"}.
#' @param greedy_em_refine Logical; if TRUE, run short csmoothEM refinement after each append.
#' @param greedy_em_max_iter Integer \eqn{\ge 0}. Number of refinement iterations per append.
#'   Set to 0 to disable refinement even if \code{greedy_em_refine=TRUE}.
#' @param discretization Discretization method used in seed initialization (recommended: \code{"quantile"}).
#' @param verbose Integer. \code{0}=silent, \code{1}=progress messages.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{coord_assign}: integer vector length d with values 1 or 2.
#'   \item \code{J1}, \code{J2}: feature indices in each ordering.
#'   \item \code{Gamma1}, \code{Gamma2}: responsibilities after greedy completion.
#'   \item \code{params1}, \code{params2}: csmooth-style params for each ordering.
#'   \item \code{fit1}, \code{fit2}: \code{csmooth_em} objects for each ordering (ready for further refinement).
#'   \item \code{seeds}: list with \code{j1}, \code{j2}.
#' }
#'
#' @seealso \code{\link{do_csmoothEM}}, \code{\link{score_feature_given_Gamma}}, \code{\link{compute_C_by_coord_csmooth}}
#' @export
forward_two_ordering_partition_csmooth <- function(
    X,
    K = 30,
    greedy_start_index = NULL,
    greedy_start_index2 = NULL,
    score_mode = c("ml","none"),
    rw_q = 2,
    relative_lambda = TRUE,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    greedy_em_refine = TRUE,
    greedy_em_max_iter = 10,
    discretization = c("quantile","equal","kmeans"),
    verbose = 1
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  score_mode <- match.arg(score_mode)
  discretization <- match.arg(discretization)

  v <- as.integer(verbose)
  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)

  rw_q <- as.integer(rw_q)
  if (rw_q < 0L) stop("rw_q must be >= 0.")

  Q1 <- make_random_walk_precision(K = K, d = 1, lambda = 1, q = rw_q, ridge = 0)
  Q2 <- make_random_walk_precision(K = K, d = 1, lambda = 1, q = rw_q, ridge = 0)

  seed_fit_1d <- function(xj, Q_K) {
    ord <- rank(xj, ties.method = "first")

    init <- make_init_csmooth(
      X = matrix(xj, ncol = 1),
      ordering_vec = ord,
      K = K,
      modelName = "homoskedastic",
      nugget = 0,
      discretization = discretization,
      na_action = "drop",
      eps = 1e-12
    )
    params <- list(pi = init$pi, mu = init$mu, sigma2 = init$sigma2)

    Gamma <- ESTEP_csmooth(matrix(xj, ncol = 1), params, "homoskedastic")
    if (any(!is.finite(Gamma))) stop("Seed ESTEP_csmooth produced non-finite Gamma. Check init/sigma2.")

    if (isTRUE(greedy_em_refine) && as.integer(greedy_em_max_iter) > 0L) {
      fit <- as_csmooth_em(
        params = params,
        gamma = Gamma,
        data = matrix(xj, ncol = 1),
        Q_K = Q_K,
        lambda_vec = 1,
        rw_q = rw_q,
        ridge = 0,
        modelName = "homoskedastic",
        relative_lambda = relative_lambda,
        nugget = 0,
        eigen_tol = NULL,
        meta = list(init = list(method = "seed_1d"))
      )
      fit$control <- list(
        adaptive = if (score_mode == "ml") "ml" else "none",
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        modelName = "homoskedastic",
        relative_lambda = relative_lambda
      )
      fit <- do_csmoothEM(
        fit, data = fit$data, iter = as.integer(greedy_em_max_iter),
        adaptive = fit$control$adaptive, record = TRUE, verbose = FALSE
      )
      Gamma <- fit$gamma
      params <- fit$params
    }

    # keep pi consistent with Gamma
    Nk <- pmax(colSums(Gamma), 1e-8)
    params$pi <- Nk / sum(Nk)

    list(Gamma = Gamma, params = params)
  }

  # ---- Seeds ----
  j1 <- if (is.null(greedy_start_index)) which.max(apply(X, 2, var)) else as.integer(greedy_start_index)
  if (length(j1) != 1L || is.na(j1) || j1 < 1L || j1 > d) stop("greedy_start_index out of range.")

  seed1 <- seed_fit_1d(X[, j1], Q1)
  Gamma1 <- seed1$Gamma
  params1 <- seed1$params

  if (!is.null(greedy_start_index2) && as.integer(greedy_start_index2) == j1) {
    warning("greedy_start_index2 equals greedy_start_index; selecting j2 automatically.")
    greedy_start_index2 <- NULL
  }

  if (!is.null(greedy_start_index2)) {
    j2 <- as.integer(greedy_start_index2)
    if (length(j2) != 1L || is.na(j2) || j2 < 1L || j2 > d) stop("greedy_start_index2 out of range.")
  } else {
    cand <- setdiff(seq_len(d), j1)
    scores_under_1 <- vapply(cand, function(j) {
      score_one_coord_csmooth(
        X = X, j = j, Gamma = Gamma1, Q_K = Q1,
        rw_q = rw_q, score_mode = score_mode,
        relative_lambda = relative_lambda,
        lambda_min = lambda_min, lambda_max = lambda_max
      )$score
    }, numeric(1))
    j2 <- cand[which.min(scores_under_1)]
  }

  seed2 <- seed_fit_1d(X[, j2], Q2)
  Gamma2 <- seed2$Gamma
  params2 <- seed2$params

  if (v >= 1) cat(sprintf("[Greedy-csmooth] seeds: j1=%d, j2=%d\n", j1, j2))

  J1 <- c(j1); J2 <- c(j2)
  remaining <- setdiff(seq_len(d), c(j1, j2))

  while (length(remaining) > 0) {
    idxs <- remaining

    svec1 <- vapply(idxs, function(j) {
      score_one_coord_csmooth(
        X = X, j = j, Gamma = Gamma1, Q_K = Q1,
        rw_q = rw_q, score_mode = score_mode,
        relative_lambda = relative_lambda,
        lambda_min = lambda_min, lambda_max = lambda_max
      )$score
    }, numeric(1))

    svec2 <- vapply(idxs, function(j) {
      score_one_coord_csmooth(
        X = X, j = j, Gamma = Gamma2, Q_K = Q2,
        rw_q = rw_q, score_mode = score_mode,
        relative_lambda = relative_lambda,
        lambda_min = lambda_min, lambda_max = lambda_max
      )$score
    }, numeric(1))

    gap_vec <- svec1 - svec2
    gap_vec[!is.finite(gap_vec)] <- 0

    pick_pos <- which.max(abs(gap_vec))
    j_star <- idxs[pick_pos]
    gap_star <- gap_vec[pick_pos]

    s1_star <- score_one_coord_csmooth(
      X = X, j = j_star, Gamma = Gamma1, Q_K = Q1,
      rw_q = rw_q, score_mode = score_mode,
      relative_lambda = relative_lambda,
      lambda_min = lambda_min, lambda_max = lambda_max
    )
    s2_star <- score_one_coord_csmooth(
      X = X, j = j_star, Gamma = Gamma2, Q_K = Q2,
      rw_q = rw_q, score_mode = score_mode,
      relative_lambda = relative_lambda,
      lambda_min = lambda_min, lambda_max = lambda_max
    )

    if (gap_star >= 0) {
      params1 <- append_coord_to_params_csmooth(params1, s1_star$one)
      J1 <- c(J1, j_star)

      X1 <- X[, J1, drop = FALSE]
      Gamma1 <- ESTEP_csmooth(X1, params1, "homoskedastic")
      if (any(!is.finite(Gamma1))) stop("Non-finite Gamma1 after append; check params/sigma2.")
      Nk <- pmax(colSums(Gamma1), 1e-8); params1$pi <- Nk / sum(Nk)

      if (isTRUE(greedy_em_refine) && as.integer(greedy_em_max_iter) > 0L) {
        fit1 <- as_csmooth_em(
          params = params1, gamma = Gamma1, data = X1,
          Q_K = Q1, lambda_vec = rep(1, ncol(X1)), rw_q = rw_q,
          ridge = 0, modelName = "homoskedastic",
          relative_lambda = relative_lambda,
          nugget = 0, eigen_tol = NULL, meta = NULL
        )
        fit1$control <- list(
          adaptive = if (score_mode == "ml") "ml" else "none",
          lambda_min = lambda_min, lambda_max = lambda_max,
          modelName = "homoskedastic", relative_lambda = relative_lambda
        )
        fit1 <- do_csmoothEM(
          fit1, data = X1, iter = as.integer(greedy_em_max_iter),
          adaptive = fit1$control$adaptive, record = TRUE, verbose = FALSE
        )
        Gamma1 <- fit1$gamma
        params1 <- fit1$params
        Nk <- pmax(colSums(Gamma1), 1e-8); params1$pi <- Nk / sum(Nk)
      }

      if (v >= 1) cat(sprintf("[Greedy-csmooth] pick %d -> ord1 (|Δ|=%.3f), remaining=%d\n",
                              j_star, abs(gap_star), length(remaining) - 1))
    } else {
      params2 <- append_coord_to_params_csmooth(params2, s2_star$one)
      J2 <- c(J2, j_star)

      X2 <- X[, J2, drop = FALSE]
      Gamma2 <- ESTEP_csmooth(X2, params2, "homoskedastic")
      if (any(!is.finite(Gamma2))) stop("Non-finite Gamma2 after append; check params/sigma2.")
      Nk <- pmax(colSums(Gamma2), 1e-8); params2$pi <- Nk / sum(Nk)

      if (isTRUE(greedy_em_refine) && as.integer(greedy_em_max_iter) > 0L) {
        fit2 <- as_csmooth_em(
          params = params2, gamma = Gamma2, data = X2,
          Q_K = Q2, lambda_vec = rep(1, ncol(X2)), rw_q = rw_q,
          ridge = 0, modelName = "homoskedastic",
          relative_lambda = relative_lambda,
          nugget = 0, eigen_tol = NULL, meta = NULL
        )
        fit2$control <- list(
          adaptive = if (score_mode == "ml") "ml" else "none",
          lambda_min = lambda_min, lambda_max = lambda_max,
          modelName = "homoskedastic", relative_lambda = relative_lambda
        )
        fit2 <- do_csmoothEM(
          fit2, data = X2, iter = as.integer(greedy_em_max_iter),
          adaptive = fit2$control$adaptive, record = TRUE, verbose = FALSE
        )
        Gamma2 <- fit2$gamma
        params2 <- fit2$params
        Nk <- pmax(colSums(Gamma2), 1e-8); params2$pi <- Nk / sum(Nk)
      }

      if (v >= 1) cat(sprintf("[Greedy-csmooth] pick %d -> ord2 (|Δ|=%.3f), remaining=%d\n",
                              j_star, abs(gap_star), length(remaining) - 1))
    }

    remaining <- remaining[-pick_pos]
  }

  coord_assign <- integer(d)
  coord_assign[J1] <- 1L
  coord_assign[J2] <- 2L

  if (v >= 1) cat(sprintf("[Greedy-csmooth] done: |J1|=%d, |J2|=%d\n", length(J1), length(J2)))

  fit1 <- as_csmooth_em(
    params = params1,
    gamma  = Gamma1,
    data   = X[, J1, drop = FALSE],
    Q_K    = Q1,
    lambda_vec = rep(1, length(J1)),
    rw_q   = rw_q,
    ridge  = 0,
    modelName = "homoskedastic",
    relative_lambda = relative_lambda,
    nugget = 0,
    eigen_tol = NULL,
    meta = list(init = list(method = "forward_greedy", seeds = list(j1 = j1, j2 = j2), group = 1))
  )
  fit1$control <- list(
    adaptive = if (score_mode == "ml") "ml" else "none",
    lambda_min = lambda_min, lambda_max = lambda_max,
    modelName = "homoskedastic", relative_lambda = relative_lambda
  )

  fit2 <- as_csmooth_em(
    params = params2,
    gamma  = Gamma2,
    data   = X[, J2, drop = FALSE],
    Q_K    = Q2,
    lambda_vec = rep(1, length(J2)),
    rw_q   = rw_q,
    ridge  = 0,
    modelName = "homoskedastic",
    relative_lambda = relative_lambda,
    nugget = 0,
    eigen_tol = NULL,
    meta = list(init = list(method = "forward_greedy", seeds = list(j1 = j1, j2 = j2), group = 2))
  )
  fit2$control <- fit1$control

  list(
    coord_assign = coord_assign,
    J1 = J1, J2 = J2,
    Gamma1 = Gamma1, Gamma2 = Gamma2,
    params1 = params1, params2 = params2,
    fit1 = fit1, fit2 = fit2,
    seeds = list(j1 = j1, j2 = j2)
  )
}



#' Drop one coordinate from csmooth-style parameters (internal)
#'
#' @description
#' Removes a single coordinate (feature) from a csmooth-style parameter list.
#' Used for warm-start backward greedy moves (remove from one partition).
#'
#' @param params List with fields \code{pi}, \code{mu}, \code{sigma2}.
#'   \code{mu} must be a list of length K, each element length d_sub.
#'   \code{sigma2} must be a numeric vector length d_sub (homoskedastic).
#' @param pos Integer in \code{1:d_sub}, position of the coordinate to remove
#'   (in the CURRENT partition's feature space).
#'
#' @return Updated \code{params}.
#' @keywords internal
drop_coord_from_params_csmooth <- function(params, pos) {
  pos <- as.integer(pos)
  if (length(pos) != 1L || is.na(pos) || pos < 1L) stop("pos must be a valid integer >= 1.")

  K <- length(params$pi)
  if (length(params$mu) != K) stop("params$mu must be a list of length K.")

  # drop from mu
  params$mu <- lapply(params$mu, function(v) {
    v <- as.numeric(v)
    if (pos > length(v)) stop("pos out of range for params$mu.")
    v[-pos]
  })

  # drop from sigma2
  s2 <- as.numeric(params$sigma2)
  if (pos > length(s2)) stop("pos out of range for params$sigma2.")
  params$sigma2 <- s2[-pos]

  params
}


#' Drop one coordinate from a csmooth_em fit (internal)
#'
#' @description
#' Removes a single feature column from a \code{csmooth_em} object, updating \code{data},
#' \code{params}, and \code{prior$lambda_vec}. This is a warm-start operation and does not
#' refit the model.
#'
#' @param fit A \code{csmooth_em} object.
#' @param pos Integer in \code{1:d_sub}, position of the coordinate to remove in the CURRENT fit.
#'
#' @return Updated \code{csmooth_em}.
#' @keywords internal
drop_coord_from_fit_csmooth <- function(fit, pos) {
  if (!inherits(fit, "csmooth_em")) stop("fit must be a csmooth_em object.")
  pos <- as.integer(pos)
  if (length(pos) != 1L || is.na(pos) || pos < 1L) stop("pos must be a valid integer >= 1.")

  if (!is.null(fit$data)) {
    if (pos > ncol(fit$data)) stop("pos out of range for fit$data.")
    fit$data <- as.matrix(fit$data)[, -pos, drop = FALSE]
  }

  fit$params <- drop_coord_from_params_csmooth(fit$params, pos)

  if (!is.null(fit$prior$lambda_vec)) {
    lv <- as.numeric(fit$prior$lambda_vec)
    if (pos > length(lv)) stop("pos out of range for fit$prior$lambda_vec.")
    fit$prior$lambda_vec <- lv[-pos]
  }

  # drop stored gamma unchanged (will be updated by do_csmoothEM)
  fit
}


#' Append one coordinate to a csmooth_em fit using a Gamma-based 1D initialization (internal)
#'
#' @description
#' Adds a new feature column to a \code{csmooth_em} object, initializing the new coordinate's
#' \eqn{\mu_{j\cdot}} and \eqn{\sigma_j^2} using \code{score_feature_given_Gamma} under the fit's
#' current responsibilities \code{fit$gamma}. Then appends the coordinate and sets its lambda
#' to 1 (it can be updated later by \code{do_csmoothEM} if adaptive is on).
#'
#' This is a warm-start operation: it does not reinitialize the whole model.
#'
#' @param fit A \code{csmooth_em} object (homoskedastic).
#' @param xj Numeric vector length n (the new feature column).
#' @param score_mode One of \code{"ml"} or \code{"none"}; controls how the 1D init is computed.
#' @param rw_q Integer RW rank deficiency along K.
#' @param relative_lambda Logical.
#' @param lambda_min,lambda_max Bounds for lambda optimization when \code{score_mode="ml"}.
#'
#' @return Updated \code{csmooth_em}.
#' @keywords internal
append_coord_to_fit_csmooth <- function(fit, xj,
                                        score_mode = c("ml","none"),
                                        rw_q = 2,
                                        relative_lambda = TRUE,
                                        lambda_min = 1e-10,
                                        lambda_max = 1e10) {
  if (!inherits(fit, "csmooth_em")) stop("fit must be a csmooth_em object.")
  if ((fit$control$modelName %||% "homoskedastic") != "homoskedastic") {
    stop("append_coord_to_fit_csmooth currently supports modelName='homoskedastic' only.")
  }

  score_mode <- match.arg(score_mode)

  xj <- as.numeric(xj)
  if (is.null(fit$data)) stop("fit$data is NULL; needed to append coordinate.")
  n <- nrow(fit$data)
  if (length(xj) != n) stop("length(xj) must match nrow(fit$data).")

  if (is.null(fit$gamma)) stop("fit$gamma is NULL; run do_csmoothEM(..., record=TRUE) first.")

  Q_K <- fit$prior$Q_K
  if (is.null(Q_K)) stop("fit$prior$Q_K is NULL.")

  # 1D init under current Gamma
  sc <- score_feature_given_Gamma(
    xj = xj,
    Gamma = fit$gamma,
    Q_K = Q_K,
    rw_q = rw_q,
    score_mode = score_mode,
    relative_lambda = relative_lambda,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    optimize_lambda = NULL,  # default: ml optimizes, none fixed lambda=1
    lambda_init = 1,
    max_sigma_iter = 1L
  )

  # append to data
  fit$data <- cbind(as.matrix(fit$data), xj)

  # append to params
  K <- length(fit$params$pi)
  fit$params$mu <- lapply(seq_len(K), function(k) c(fit$params$mu[[k]], sc$mu_vec[k]))
  fit$params$sigma2 <- c(as.numeric(fit$params$sigma2), as.numeric(sc$sigma2))

  # append lambda=1 (warm-start); adaptive updates can change it later
  fit$prior$lambda_vec <- c(as.numeric(fit$prior$lambda_vec), 1)

  fit
}





#' Backward greedy feature partition into two orderings (csmoothEM, warm-start only)
#'
#' @description
#' Warm-start backward greedy algorithm for partitioning features into two latent orderings.
#' The procedure initializes ordering 1 once on all features, initializes ordering 2 once
#' from a single seed feature via a rank-based 1D seed (no global ordering method), and then
#' performs greedy moves using only warm-start updates (drop/append + a few EM iterations).
#'
#' At each step, the algorithm selects the feature in ordering 1 with the largest positive
#' gain \eqn{\Delta_j = S_2(j \mid \Gamma_2) - S_1(j \mid \Gamma_1)} and moves it to ordering 2,
#' where \eqn{S_k} is computed by \code{\link{score_feature_given_Gamma}} under the current responsibilities.
#'
#' This function never re-calls ordering initializers (e.g. fiedler/pcurve/PCA) after the initial
#' construction. Ordering updates are done only via warm-start csmoothEM iterations.
#'
#' @param X Numeric matrix \code{(n x d)}.
#' @param K Integer \eqn{\ge 2}. Number of mixture components.
#' @param init_method1 Ordering initializer for ordering 1, passed to \code{\link{initialize_csmoothEM}}.
#' @param seed2 Optional integer. If NULL, chosen as the worst-aligned feature under ordering 1.
#' @param score_mode One of \code{"ml"} or \code{"none"}; passed to \code{\link{score_feature_given_Gamma}}.
#' @param adaptive Adaptive mode used in \code{\link{do_csmoothEM}} warm updates (\code{"ml"} or \code{"none"}).
#'   If NULL, defaults to \code{score_mode}.
#' @param rw_q Integer \eqn{\ge 0}. RW rank deficiency.
#' @param relative_lambda Logical.
#' @param lambda_min,lambda_max Bounds for lambda optimization when \code{score_mode="ml"}.
#' @param discretization Discretization used in initialization (recommended: \code{"quantile"}).
#' @param warm_iter_init Initial warm-start iterations for ordering 1.
#' @param warm_iter_refit Warm-start iterations after each move (applied to both fits).
#' @param max_steps Maximum number of greedy moves.
#' @param verbose Logical.
#' @param ... Passed to \code{\link{initialize_csmoothEM}} for ordering 1.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{assign}: character vector length d in \code{c("A","B")}.
#'   \item \code{J1}, \code{J2}: feature indices (original X column indices) in each ordering.
#'   \item \code{seedB}: seed feature index for ordering 2.
#'   \item \code{fit1}, \code{fit2}: \code{csmooth_em} objects for each ordering.
#'   \item \code{history}: data.frame of moves (feature, gain, sizes).
#' }
#'
#' @seealso \code{\link{initialize_csmoothEM}}, \code{\link{do_csmoothEM}}, \code{\link{score_feature_given_Gamma}},
#'   \code{\link{append_coord_to_fit_csmooth}}, \code{\link{drop_coord_from_fit_csmooth}}
#' @export
backward_two_ordering_partition_csmooth <- function(
    X,
    K = 50,
    init_method1 = c("fiedler", "PCA", "tSNE", "pcurve", "random"),
    seed2 = NULL,
    score_mode = c("ml", "none"),
    adaptive = NULL,
    rw_q = 2,
    relative_lambda = TRUE,
    lambda_min = 1e-10,
    lambda_max = 1e10,
    discretization = c("quantile", "equal", "kmeans"),
    warm_iter_init = 10,
    warm_iter_refit = 5,
    max_steps = 50,
    verbose = TRUE,
    ...
) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)
  fn <- colnames(X); if (is.null(fn)) fn <- paste0("V", seq_len(d))

  init_method1 <- match.arg(init_method1)
  score_mode <- match.arg(score_mode)
  discretization <- match.arg(discretization)

  if (is.null(adaptive)) adaptive <- if (score_mode == "ml") "ml" else "none"
  if (!adaptive %in% c("ml","none")) stop("adaptive must be 'ml' or 'none'.")

  warm_iter_init <- as.integer(warm_iter_init)
  warm_iter_refit <- as.integer(warm_iter_refit)
  max_steps <- as.integer(max_steps)
  if (warm_iter_init < 1L || warm_iter_refit < 1L || max_steps < 1L) stop("iters must be >= 1.")

  # ---------- (1) Fit ordering 1 ONCE on all features ----------
  fit1 <- initialize_csmoothEM(
    X = X,
    method = init_method1,
    K = K,
    modelName = "homoskedastic",
    adaptive = adaptive,
    discretization = discretization,
    num_iter = warm_iter_init,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    include.data = TRUE,
    ...
  )
  fit1 <- do_csmoothEM(fit1, data = X, iter = 1, adaptive = NULL, record = TRUE, verbose = FALSE)
  if (is.null(fit1$gamma)) stop("fit1$gamma is NULL; ensure do_csmoothEM stores gamma.")

  # ---------- (2) Choose seedB ----------
  if (is.null(seed2)) {
    scores1 <- vapply(seq_len(d), function(j) {
      score_feature_given_Gamma(
        xj = X[, j],
        Gamma = fit1$gamma,
        Q_K = fit1$prior$Q_K,
        rw_q = rw_q,
        score_mode = score_mode,
        relative_lambda = relative_lambda,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        optimize_lambda = NULL,
        lambda_init = 1,
        max_sigma_iter = 1L
      )$score
    }, numeric(1))
    seedB <- which.min(scores1)
  } else {
    seedB <- as.integer(seed2)
    if (seedB < 1L || seedB > d) stop("seed2 out of range.")
  }

  # ---------- (3) Initialize ordering 2 ONCE using a rank-based 1D seed ----------
  # This avoids calling fiedler_ordering/pcurve/PCA on a 1D matrix.
  x_seed <- X[, seedB]
  ord <- rank(x_seed, ties.method = "first")

  init2 <- make_init_csmooth(
    X = matrix(x_seed, ncol = 1),
    ordering_vec = ord,
    K = K,
    modelName = "homoskedastic",
    nugget = 0,
    discretization = discretization,
    na_action = "drop",
    eps = 1e-12
  )

  Q_K2 <- make_random_walk_precision(K = length(init2$pi), d = 1, lambda = 1, q = rw_q, ridge = 0)

  fit2 <- as_csmooth_em(
    params = list(pi = init2$pi, mu = init2$mu, sigma2 = init2$sigma2),
    gamma = NULL,
    data = matrix(x_seed, ncol = 1),
    Q_K = Q_K2,
    lambda_vec = 1,
    rw_q = rw_q,
    ridge = 0,
    modelName = "homoskedastic",
    relative_lambda = relative_lambda,
    nugget = 0,
    eigen_tol = NULL,
    meta = list(init = list(method = "seed_rank_1d", seedB = seedB))
  )
  fit2$control <- list(
    adaptive = adaptive,
    lambda_min = lambda_min,
    lambda_max = lambda_max,
    modelName = "homoskedastic",
    relative_lambda = relative_lambda
  )
  fit2 <- do_csmoothEM(fit2, data = fit2$data, iter = warm_iter_refit, adaptive = adaptive, record = TRUE, verbose = FALSE)
  if (is.null(fit2$gamma)) stop("fit2$gamma is NULL after seed initialization.")

  # ---------- (4) Make partitions disjoint immediately ----------
  # Global sets
  J2 <- seedB
  J1 <- setdiff(seq_len(d), seedB)

  # fit1 initially contains all d columns in original order -> drop seedB from fit1
  map1_global <- seq_len(d)
  pos_seed_in_fit1 <- match(seedB, map1_global)
  fit1 <- drop_coord_from_fit_csmooth(fit1, pos_seed_in_fit1)
  map1_global <- map1_global[-pos_seed_in_fit1]

  # fit2 currently contains only seedB
  map2_global <- seedB

  if (verbose) {
    cat(sprintf("[backward-init] seedB=%d (%s), |A|=%d |B|=%d\n",
                seedB, fn[seedB], length(J1), length(J2)))
  }

  history <- data.frame(
    step = integer(0),
    moved = integer(0),
    moved_name = character(0),
    gain = numeric(0),
    nA = integer(0),
    nB = integer(0),
    stringsAsFactors = FALSE
  )

  # ---------- (5) Greedy moves with WARM START ONLY ----------
  for (step in seq_len(max_steps)) {

    if (length(J1) < 1L) break

    gains <- vapply(J1, function(jg) {
      s1 <- score_feature_given_Gamma(
        xj = X[, jg],
        Gamma = fit1$gamma,
        Q_K = fit1$prior$Q_K,
        rw_q = rw_q,
        score_mode = score_mode,
        relative_lambda = relative_lambda,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        optimize_lambda = NULL,
        lambda_init = 1,
        max_sigma_iter = 1L
      )$score

      s2 <- score_feature_given_Gamma(
        xj = X[, jg],
        Gamma = fit2$gamma,
        Q_K = fit2$prior$Q_K,
        rw_q = rw_q,
        score_mode = score_mode,
        relative_lambda = relative_lambda,
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        optimize_lambda = NULL,
        lambda_init = 1,
        max_sigma_iter = 1L
      )$score

      s2 - s1
    }, numeric(1))

    best_pos <- which.max(gains)
    best_gain <- gains[best_pos]

    if (!is.finite(best_gain) || best_gain <= 0) {
      if (verbose) cat(sprintf("[backward %02d] no positive gain; stopping.\n", step))
      break
    }

    j_move <- J1[best_pos]

    # update global sets
    J1 <- setdiff(J1, j_move)
    J2 <- c(J2, j_move)

    # remove from fit1
    pos1 <- match(j_move, map1_global)
    if (is.na(pos1)) stop("internal error: j_move not found in map1_global.")
    fit1 <- drop_coord_from_fit_csmooth(fit1, pos1)
    map1_global <- map1_global[-pos1]

    # append to fit2 (Gamma-based 1D init under current Gamma2)
    fit2 <- append_coord_to_fit_csmooth(
      fit2,
      xj = X[, j_move],
      score_mode = score_mode,
      rw_q = rw_q,
      relative_lambda = relative_lambda,
      lambda_min = lambda_min,
      lambda_max = lambda_max
    )
    map2_global <- c(map2_global, j_move)

    # warm-start update gamma/params
    fit1 <- do_csmoothEM(fit1, data = fit1$data, iter = warm_iter_refit, adaptive = adaptive, record = TRUE, verbose = FALSE)
    fit2 <- do_csmoothEM(fit2, data = fit2$data, iter = warm_iter_refit, adaptive = adaptive, record = TRUE, verbose = FALSE)

    history <- rbind(history, data.frame(
      step = step,
      moved = j_move,
      moved_name = fn[j_move],
      gain = best_gain,
      nA = length(J1),
      nB = length(J2),
      stringsAsFactors = FALSE
    ))

    if (verbose) {
      cat(sprintf("[backward %02d] move %d (%s) gain=%.3f |A|=%d |B|=%d\n",
                  step, j_move, fn[j_move], best_gain, length(J1), length(J2)))
    }
  }

  assign <- rep("A", d)
  assign[J2] <- "B"

  list(
    assign = assign,
    J1 = J1,
    J2 = J2,
    seedB = seedB,
    fit1 = fit1,
    fit2 = fit2,
    history = history
  )
}
