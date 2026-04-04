# Initialize csmoothEM (coordinate-specific SmoothEM)

Create a `csmooth_em` object using the same ordering-initialization
ideas as
[`initialize_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/initialize_smoothEM.md),
but with coordinate-specific penalties `lambda_vec` and diagonal
covariance structures only.

Workflow:

1.  Initialize an ordering and discretize into `K` components using
    [`initialize_ordering_csmooth()`](https://aguerozz.github.io/MPCurver/reference/initialize_ordering_csmooth.md).

2.  Build the base random-walk precision `Q_K` along components (with
    `lambda=1`).

3.  Initialize `lambda_vec`. If `adaptive` is not `"none"`, estimate an
    initial `lambda_vec` from the initialized parameters (and clamp to
    `[lambda_min, lambda_max]`).

4.  Construct a `csmooth_em` object via
    [`as_csmooth_em()`](https://aguerozz.github.io/MPCurver/reference/as_csmooth_em.md),
    then run a warm start via
    [`do_csmoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM.md)
    for `num_iter` iterations.

## Usage

``` r
initialize_csmoothEM(
  X,
  method = c("tSNE", "PCA", "pcurve", "random", "fiedler", "multi_scale", "isomap"),
  rw_q = 2,
  lambda = 1,
  relative_lambda = TRUE,
  K = NULL,
  num_iter = 1,
  modelName = c("homoskedastic", "heteroskedastic"),
  ridge = 0,
  nugget = 0,
  eigen_tol = NULL,
  include.data = TRUE,
  adaptive = "ml",
  lambda_min = 1e-08,
  lambda_max = 1e+08,
  sigma_update = c("ml", "mstep"),
  sigma_min = 1e-10,
  sigma_max = 1e+10,
  discretization = c("equal", "quantile", "kmeans"),
  ...
)
```

## Arguments

- X:

  n-by-d numeric matrix.

- method:

  One of
  `"tSNE","PCA","random","fiedler","multi_scale", "pcurve","isomap"`.
  Currently `"multi_scale"` is not implemented for csmoothEM and will
  error.

- rw_q:

  Integer RW order along K for `Q_K`.

- lambda:

  Scalar or length-d vector. If scalar, recycled to length d. Used as an
  initial value when `adaptive="none"`.

- relative_lambda:

  Logical; if TRUE, scale the base prior for coordinate `j` by a
  coordinate-specific variance scale (see Details).

- K:

  Number of mixture components. If NULL, defaults to
  `min(50, floor(n/5))` (at least 2).

- num_iter:

  Integer \>= 1; number of warm-start iterations to run immediately.

- modelName:

  Either `"homoskedastic"` or `"heteroskedastic"`.

- ridge:

  Nonnegative ridge added when building RW precision `Q_K`.

- nugget:

  Nonnegative nugget used in M-step variance updates.

- eigen_tol:

  Optional eigen tolerance used for generalized log-determinants (if
  needed downstream).

- include.data:

  Logical; store the (kept) data matrix in the returned object.

- adaptive:

  Character specifying how to initialize/update `lambda_vec`:

  `"none"`

  :   Do not estimate `lambda_vec` from the initialization; use
      `lambda`.

  `"prior"`

  :   Obsolete. Retained only for backward compatibility; uses a
      prior/profile estimate from the initialized component means
      \\\mu\\.

  `"ml"`

  :   Collapsed-ML warm start (dispatches to
      [`do_csmoothEM_ml_collapsed()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_ml_collapsed.md)
      inside
      [`do_csmoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM.md)),
      allowing optional ML updates of \\\sigma^2\\.

  Logical values are accepted for backward compatibility: `TRUE` is
  interpreted as `"ml"` and `FALSE` as `"none"`.

- lambda_min, lambda_max:

  Positive bounds for `lambda_vec` (used when `adaptive!="none"`).

- sigma_update:

  Character. Only used when `adaptive="ml"` during the warm start.
  Defaults to `"ml"`; legacy `"mstep"` is retained only for backward
  compatibility. See
  [`do_csmoothEM_ml_collapsed`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_ml_collapsed.md).

- sigma_min, sigma_max:

  Positive bounds for `sigma2` when `adaptive="ml"` and
  `sigma_update="ml"`.

- discretization:

  Discretization method passed to
  [`initialize_ordering_csmooth()`](https://aguerozz.github.io/MPCurver/reference/initialize_ordering_csmooth.md).

- ...:

  Passed to the ordering method (e.g. PCA/tSNE/pcurve/fiedler/isomap).

## Value

A `csmooth_em` object.

## Details

Let \\Q_K\\ denote the RW(q) precision along components. For a separable
prior \\Q\_{\mathrm{full}} = I_d \otimes Q_K\\, the rank deficiency is
\\rw_q\\ per coordinate, so a convenient degrees-of-freedom term is \\r
= K - rw_q\\.

When `adaptive = "prior"` (obsolete), the initializer estimates each
coordinate penalty by: \$\$\lambda_j = r \big/ \left( \mu\_{j\cdot}^\top
Q\_{j,\mathrm{base}} \mu\_{j\cdot} \right),\$\$ where \\\mu\_{j\cdot}\\
is the length-K vector of component means for coordinate \\j\\, and
\\Q\_{j,\mathrm{base}}\\ is the coordinate-specific base precision
returned by
[`.compute_Qbase_j()`](https://aguerozz.github.io/MPCurver/reference/dot-compute_Qbase_j.md)
(equivalently, \\Q_K\\ possibly rescaled when `relative_lambda = TRUE`).
The estimate is clamped to \\\[\text{lambda_min},\text{lambda_max}\]\\.
This mode is heuristic and should not be treated as a
convergence-guaranteed update.

When `adaptive = "ml"`, the warm start uses the collapsed-ML routine. In
this mode,
[`do_csmoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM.md)
dispatches to
[`do_csmoothEM_ml_collapsed()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_ml_collapsed.md),
which optimizes a collapsed (Laplace-exact, Gaussian) objective and
records `ml_trace` as the collapsed objective \\\mathcal{C}\\.

When `relative_lambda = TRUE`, the base precision
\\Q\_{j,\mathrm{base}}\\ is scaled by a coordinate variance proxy:

- `modelName = "homoskedastic"`: scale by `1 / sigma2[j]`.

- `modelName = "heteroskedastic"`: scale by
  `1 / sum_k pi_k sigma2[j,k]`.
