# =========================
# VAR(1) precision
# =========================

#' Prior precision for VAR(1) on K time points (d-dimensional)
#'
#' Model: U_1 ~ N(0, P0^{-1}), and for k>=2: U_k | U_{k-1} ~ N(A U_{k-1}, Q)
#' Returns the joint precision matrix for vec(U_1,...,U_K).
#'
#' @param K integer number of time points.
#' @param d integer state dimension.
#' @param A d x d transition matrix (default 0.8 * I_d).
#' @param Q d x d innovation covariance (default 0.5 * I_d).
#' @param P0 d x d initial precision (default I_d).
#' @return (dK) x (dK) precision matrix (dense).
#' @export
make_VAR1_precision <- function(K, d, A = NULL, Q = NULL, P0 = NULL) {
  if (K < 1 || K != as.integer(K)) stop("K must be a positive integer.")
  if (d < 1 || d != as.integer(d)) stop("d must be a positive integer.")

  if (is.null(A))  A  <- diag(d) * 0.8
  if (is.null(Q))  Q  <- diag(d) * 0.5
  if (is.null(P0)) P0 <- diag(d)

  A  <- as.matrix(A)
  Q  <- as.matrix(Q)
  P0 <- as.matrix(P0)

  if (!all(dim(A)  == c(d, d))) stop("A must be d x d.")
  if (!all(dim(Q)  == c(d, d))) stop("Q must be d x d.")
  if (!all(dim(P0) == c(d, d))) stop("P0 must be d x d.")

  # Q^{-1} via Cholesky (more stable than solve() for SPD Q)
  Qinv <- tryCatch({
    chol2inv(chol(Q))
  }, error = function(e) stop("Q must be symmetric positive definite (chol failed)."))

  At_Qinv <- t(A) %*% Qinv
  At_Qinv_A <- At_Qinv %*% A

  H <- matrix(0, d * K, d * K)

  # prior on U1
  H[1:d, 1:d] <- P0

  if (K >= 2) {
    for (k in 2:K) {
      idx_k   <- ((k - 1) * d + 1):(k * d)
      idx_km1 <- ((k - 2) * d + 1):((k - 1) * d)

      # diagonal blocks
      H[idx_k,   idx_k]   <- H[idx_k,   idx_k]   + Qinv
      H[idx_km1, idx_km1] <- H[idx_km1, idx_km1] + At_Qinv_A

      # off-diagonal blocks
      H[idx_k,   idx_km1] <- H[idx_k,   idx_km1] - Qinv %*% A
      H[idx_km1, idx_k]   <- H[idx_km1, idx_k]   - At_Qinv
    }
  }

  # enforce symmetry (numerical)
  H <- (H + t(H)) / 2
  H
}

# =========================
# 1D random walk precision
# =========================

#' Prior precision for q-th order random walk on K time points (d-dimensional)
#'
#' Constructs Q = lambda * (D_q' D_q) \kron I_d, with optional ridge.
#'
#' @param K integer number of time points.
#' @param d integer dimension of the process at each time.
#' @param q integer order of the random walk (q>=1).
#' @param lambda nonnegative scalar precision multiplier.
#' @param ridge nonnegative scalar ridge added to the 1D precision.
#' @param sparse logical; return a sparse Matrix object if TRUE.
#' @return (dK) x (dK) precision matrix.
#' @export
make_random_walk_precision <- function(K, d, q = 1, lambda = 1, ridge = 0, sparse = FALSE) {
  if (K < 1 || K != as.integer(K)) stop("K must be a positive integer.")
  if (d < 1 || d != as.integer(d)) stop("d must be a positive integer.")
  if (q < 1 || q != as.integer(q) || q >= K) stop("q must be an integer in {1, ..., K-1}.")
  if (lambda < 0) stop("lambda must be >= 0.")
  if (ridge < 0) stop("ridge must be >= 0.")

  # Dense construction (fine for moderate K); for large K you may want a fully sparse builder later
  D_q <- diag(K)
  for (i in seq_len(q)) D_q <- diff(D_q)

  Q_scalar <- crossprod(D_q) + ridge * diag(K)
  Q_U <- lambda * kronecker(Q_scalar, diag(d))

  if (sparse) {
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("Please install.packages('Matrix')")
    Q_U <- Matrix::Matrix(Q_U, sparse = TRUE)
  }

  Q_U
}

# =========================
# 2D lattice random walk precision
# =========================

#' Prior precision for q-th order random walk on a K x K lattice (d-dimensional)
#'
#' Builds a 2D intrinsic RW(q) precision via Kronecker sum:
#' Q2D = kron(Q1D, I) + kron(I, Q1D), then replicates across d independent processes.
#'
#' @param K integer lattice side length (total sites = K^2).
#' @param d integer number of independent processes.
#' @param q integer order of the random walk (q>=1).
#' @param lambda nonnegative scalar precision multiplier.
#' @param ridge nonnegative ridge added to improve conditioning.
#' @return sparse (d*K^2) x (d*K^2) precision matrix (Matrix).
#' @export
make_lattice_rwq_precision <- function(K, d, q = 1, lambda = 1, ridge = 0) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Please install.packages('Matrix')")
  }
  if (K < 1 || K != as.integer(K)) stop("K must be a positive integer.")
  if (d < 1 || d != as.integer(d)) stop("d must be a positive integer.")
  if (q < 1 || q != as.integer(q) || q >= K) stop("q must be an integer in {1, ..., K-1}.")
  if (lambda < 0) stop("lambda must be >= 0.")
  if (ridge < 0) stop("ridge must be >= 0.")

  # 1D q-th order difference operator (dense, but small-ish unless K huge)
  D <- diag(K)
  for (i in seq_len(q)) D <- diff(D)

  Q1D <- Matrix::Matrix(crossprod(D), sparse = TRUE)

  I_K <- Matrix::Diagonal(K)
  Q2D <- Matrix::kronecker(Q1D, I_K) + Matrix::kronecker(I_K, Q1D)

  if (ridge > 0) {
    Q2D <- Q2D + ridge * Matrix::Diagonal(K * K)
  }

  Q_U <- lambda * Matrix::kronecker(Q2D, Matrix::Diagonal(d))
  Q_U
}

#' Sparse prior precision for q-th order random walk on K time points (d-dimensional)
#'
#' Builds the q-th forward-difference matrix D (size (K-q) x K) with coefficients
#' (-1)^(q-j) * choose(q, j), then returns:
#'   Q_U = lambda * (t(D) %*% D + ridge * I_K) \kron I_d
#'
#' @param K integer number of time points.
#' @param d integer dimension per time point.
#' @param q integer RW order (1 <= q < K).
#' @param lambda nonnegative scalar multiplier.
#' @param ridge nonnegative scalar ridge added to Q1D.
#' @return sparse Matrix of size (d*K) x (d*K).
#' @export
make_random_walk_precision_sparse <- function(K, d, q = 1, lambda = 1, ridge = 0) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Please install.packages('Matrix')")
  }
  if (K < 1 || K != as.integer(K)) stop("K must be a positive integer.")
  if (d < 1 || d != as.integer(d)) stop("d must be a positive integer.")
  if (q < 1 || q != as.integer(q) || q >= K) stop("q must be an integer in {1, ..., K-1}.")
  if (lambda < 0) stop("lambda must be >= 0.")
  if (ridge < 0) stop("ridge must be >= 0.")

  nrows <- K - q

  # q-th forward difference coefficients:
  # Δ^q x_i = sum_{j=0}^q (-1)^(q-j) * choose(q, j) * x_{i+j}
  coeff <- stats::choose(q, 0:q) * (-1)^(q - (0:q))

  # Build sparse D: rows i=1..K-q, cols i+j
  # Each offset j contributes one diagonal band
  i_list <- vector("list", q + 1)
  j_list <- vector("list", q + 1)
  x_list <- vector("list", q + 1)

  for (jj in 0:q) {
    ii <- seq_len(nrows)
    cc <- ii + jj
    i_list[[jj + 1]] <- ii
    j_list[[jj + 1]] <- cc
    x_list[[jj + 1]] <- rep(coeff[jj + 1], nrows)
  }

  Dq <- Matrix::sparseMatrix(
    i = unlist(i_list),
    j = unlist(j_list),
    x = unlist(x_list),
    dims = c(nrows, K)
  )

  Q1D <- Matrix::crossprod(Dq) # K x K sparse band matrix
  if (ridge > 0) {
    Q1D <- Q1D + ridge * Matrix::Diagonal(K)
  }

  Q_U <- lambda * Matrix::kronecker(Q1D, Matrix::Diagonal(d))
  Q_U
}


#' Sparse prior precision for q-th order random walk on a K x K lattice (d-dimensional)
#'
#' Builds sparse 1D RW(q) precision Q1D = t(Dq)%*%Dq (optionally + ridge*I),
#' then constructs 2D lattice precision via Kronecker sum:
#'   Q2D = kron(Q1D, I_K) + kron(I_K, Q1D),
#' and replicates across d independent processes:
#'   Q_U = lambda * kron(Q2D, I_d).
#'
#' @param K integer lattice side length (total sites = K^2).
#' @param d integer number of independent processes.
#' @param q integer RW order (1 <= q < K).
#' @param lambda nonnegative scalar multiplier.
#' @param ridge nonnegative ridge added to Q1D (and thus Q2D) for conditioning.
#' @return sparse Matrix of size (d*K^2) x (d*K^2).
#' @export
make_lattice_rwq_precision_sparse <- function(K, d, q = 1, lambda = 1, ridge = 0) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Please install.packages('Matrix')")
  }
  if (K < 1 || K != as.integer(K)) stop("K must be a positive integer.")
  if (d < 1 || d != as.integer(d)) stop("d must be a positive integer.")
  if (q < 1 || q != as.integer(q) || q >= K) stop("q must be an integer in {1, ..., K-1}.")
  if (lambda < 0) stop("lambda must be >= 0.")
  if (ridge < 0) stop("ridge must be >= 0.")

  # ---- 1D RW(q): build sparse forward-difference Dq of size (K-q) x K ----
  nrows <- K - q
  coeff <- stats::choose(q, 0:q) * (-1)^(q - (0:q))  # length q+1

  i_list <- vector("list", q + 1)
  j_list <- vector("list", q + 1)
  x_list <- vector("list", q + 1)

  for (jj in 0:q) {
    ii <- seq_len(nrows)
    cc <- ii + jj
    i_list[[jj + 1]] <- ii
    j_list[[jj + 1]] <- cc
    x_list[[jj + 1]] <- rep(coeff[jj + 1], nrows)
  }

  Dq <- Matrix::sparseMatrix(
    i = unlist(i_list),
    j = unlist(j_list),
    x = unlist(x_list),
    dims = c(nrows, K)
  )

  Q1D <- Matrix::crossprod(Dq)  # K x K sparse band matrix
  if (ridge > 0) Q1D <- Q1D + ridge * Matrix::Diagonal(K)

  # ---- 2D lattice via Kronecker sum ----
  I_K <- Matrix::Diagonal(K)
  Q2D <- Matrix::kronecker(Q1D, I_K) + Matrix::kronecker(I_K, Q1D)

  # ---- replicate across d independent processes ----
  Q_U <- lambda * Matrix::kronecker(Q2D, Matrix::Diagonal(d))
  Q_U
}
