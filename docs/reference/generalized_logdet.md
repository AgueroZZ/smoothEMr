# Generalized log-determinant (sum of log positive eigenvalues)

Computes \\\sum_i \log(\lambda_i)\\ over eigenvalues \\\lambda_i\\ of a
symmetric matrix `Q` that exceed `eigen_tol`.

## Usage

``` r
generalized_logdet(Q, eigen_tol = NULL, rank_deficiency = 0)
```

## Arguments

- Q:

  Symmetric matrix.

- eigen_tol:

  Nonnegative threshold; eigenvalues \<= eigen_tol are ignored. If
  `NULL`, uses a heuristic based on machine precision.

- rank_deficiency:

  Optional integer. If \> 0, drops the smallest `rank_deficiency`
  eigenvalues (useful when you know the null-space dimension a priori).

## Value

Scalar generalized log-determinant.
