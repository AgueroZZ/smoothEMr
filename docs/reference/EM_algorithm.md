# Penalized Expectation-Maximization (EM) algorithm for SmoothEM

Fits a Gaussian mixture model with optional quadratic prior penalty on
the stacked mean vector `mu`. The algorithm alternates E-steps and
penalized M-steps and records both a penalized observed-data objective
and a penalized ELBO trace.

## Usage

``` r
EM_algorithm(
  data,
  init_params,
  Q_prior = NULL,
  iterate_once = TRUE,
  max_inner = 10,
  modelName = "VVV",
  max_iter = 100,
  tol = 1e-04,
  inner_tol = 1e-06,
  eigen_tol = NULL,
  rank_deficiency = 0,
  nugget = 0,
  relative_lambda = FALSE,
  verbose = TRUE,
  include.data = TRUE
)
```

## Arguments

- data:

  Numeric matrix `(n x d)` of observations.

- init_params:

  Initial parameters as a list with `pi`, `mu`, `sigma`.

- Q_prior:

  Optional precision matrix on stacked means (`d*K x d*K`).

- iterate_once:

  Logical; if TRUE, stop after one EM iteration (often used inside
  wrappers).

- max_inner:

  Maximum number of inner iterations inside the M-step.

- modelName:

  Covariance model: one of `"VVV"`, `"VII"`, `"EII"`, `"EEI"`.

- max_iter:

  Maximum number of EM iterations.

- tol:

  Convergence tolerance on ELBO changes.

- inner_tol:

  Tolerance for inner iterations inside the M-step.

- eigen_tol:

  Optional tolerance passed to
  [`generalized_logdet()`](https://aguerozz.github.io/MPCurver/reference/generalized_logdet.md).

- rank_deficiency:

  Rank deficiency for pseudo-determinant correction (commonly `q*d`).

- nugget:

  Nonnegative diagonal jitter added to covariance estimates.

- relative_lambda:

  Logical; if TRUE, rescales `Q_prior` by current marginal variances.

- verbose:

  Logical; print progress.

- include.data:

  Logical; if TRUE, include `data` in the returned list.

## Value

A list with:

- `params`: final parameters (including cached `invSigma` and `logdet`).

- `gamma`: responsibility matrix `(n x K)`.

- `elbo_trace`: numeric vector.

- `loglik_trace`: numeric vector.
