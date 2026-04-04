# Initialize SmoothEM (single- or multi-scale)

Creates a `smooth_em` object using a chosen initialization method and
then runs SmoothEM for a specified number of iterations.

The initialization is always "warm-started" by running exactly one outer
EM iteration (via `EM_algorithm(max_iter = 1)` for non-multi-scale; or
`progressive_smoothEM(max_iter = 1)` per stage for multi-scale). If
`num_iter > 1`, the remaining iterations are completed by
[`do_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_smoothEM.md),
which enables features such as adaptive penalty updates.

When `adaptive = TRUE`, an initial \\\lambda\\ is estimated
("profile-style") from the current mean vector \\u\\ (stacked `mu`) and
the base precision `Q_base`: \$\$\lambda^\\(u) = r / (u^\top
Q\_\mathrm{base} u), \quad r = p - \mathrm{rank\\deficiency}\$\$ where
\\p\\ is the dimension of `Q_base`. The estimated value is clipped to
`[lambda_min, lambda_max]` and used as the starting \\\lambda\\.

## Usage

``` r
initialize_smoothEM(
  X,
  method = c("tSNE", "PCA", "random", "multi_scale", "fiedler", "pcurve", "isomap"),
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
  lambda_min = 1e-08,
  lambda_max = 1e+08,
  ...
)
```

## Arguments

- X:

  Numeric matrix `(n x d)` of observations.

- method:

  Initialization method. If not `"multi_scale"`, it is passed to
  [`initialize_ordering()`](https://aguerozz.github.io/MPCurver/reference/initialize_ordering.md).
  If `"multi_scale"`, a progressive coarse-to-fine initialization is
  used via
  [`progressive_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/progressive_smoothEM.md).

- rw_q:

  Integer random-walk order `q` for the separable RW penalty (default
  2).

- lambda:

  Nonnegative penalty strength. For `method="multi_scale"`, this is
  interpreted as `lambda_final`.

- relative_lambda:

  Logical; if TRUE, rescales the precision used for fitting/evaluation
  by the current marginal variances (intended for EEI-like shared
  covariance across clusters).

- K:

  Integer number of grid points (used only when
  `method != "multi_scale"`). If NULL, defaults to
  `min(50, floor(nrow(X)/5))` with minimum 2.

- m_max:

  Integer finest exponent for multi-scale grid; final grid size is
  `K_final = 2^m_max + 1`.

- num_iter:

  Integer \>= 1; total number of SmoothEM outer iterations to run.
  Exactly one iteration is run in the warm-start step, and the remaining
  `num_iter - 1` iterations (if any) are run by
  [`do_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_smoothEM.md).

- modelName:

  Covariance model passed to
  [`EM_algorithm()`](https://aguerozz.github.io/MPCurver/reference/EM_algorithm.md)
  / `MSTEP()`. Common choices include `"VVV"`, `"VII"`, `"EII"`,
  `"EEI"`.

- ridge:

  Nonnegative ridge added in constructing the RW precision.

- nugget:

  Nonnegative diagonal jitter added to covariance estimates during
  fitting.

- eigen_tol:

  Optional tolerance passed to generalized log-determinant routines used
  in objective/ELBO evaluation.

- keep_history:

  Logical; if TRUE and `method=="multi_scale"`, attach the full
  progressive initialization result under
  `meta$init$details$progressive`.

- include.data:

  Logical; if TRUE, store the data matrix in the returned object.

- adaptive:

  Logical; if TRUE, estimate a starting \\\lambda\\ from the initial
  mean vector and enable adaptive updates during continuation via
  [`do_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_smoothEM.md).

- lambda_min, lambda_max:

  Positive bounds used to clip the estimated/updated \\\lambda\\ when
  `adaptive=TRUE`.

- ...:

  Extra arguments passed to
  [`initialize_ordering()`](https://aguerozz.github.io/MPCurver/reference/initialize_ordering.md)
  (when `method != "multi_scale"`), or to
  [`progressive_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/progressive_smoothEM.md)
  (when `method == "multi_scale"`).

## Value

A `smooth_em` object with fitted parameters, responsibilities, and
traces. Initialization provenance is stored in `meta$init`, including
`lambda_init_est` and `lambda_init_source` when `adaptive=TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(rnorm(200 * 2), 200, 2)

# Fixed lambda
fit0 <- initialize_smoothEM(X, method = "PCA", num_iter = 10, adaptive = FALSE, lambda = 10)

# Adaptive lambda (estimate start + adapt during continuation)
fit1 <- initialize_smoothEM(X, method = "tSNE", num_iter = 10, adaptive = TRUE, lambda = 10)

# Multi-scale initialization
fit2 <- initialize_smoothEM(X, method = "multi_scale", m_max = 6, num_iter = 5, adaptive = TRUE)
} # }
```
