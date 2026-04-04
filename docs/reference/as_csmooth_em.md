# Construct a csmooth_em object

Create a coordinate-specific SmoothEM object (`csmooth_em`) with
diagonal covariance and a separable RW prior along the K dimension. Each
coordinate \\j\\ has its own penalty strength \\\lambda_j\\.

Supported covariance models (diagonal only):

- `"homoskedastic"`: \\\sigma^2_j\\ shared across clusters (but varies
  across coordinates).

- `"heteroskedastic"`: \\\sigma^2\_{j,k}\\ varies across both
  coordinates and clusters.

## Usage

``` r
as_csmooth_em(
  params,
  gamma = NULL,
  data = NULL,
  Q_K,
  lambda_vec,
  rw_q = 2,
  ridge = 0,
  modelName = c("homoskedastic", "heteroskedastic"),
  relative_lambda = TRUE,
  nugget = 0,
  eigen_tol = NULL,
  meta = NULL
)
```

## Arguments

- params:

  List with fields:

  - `pi`: length-K mixing proportions.

  - `mu`: list of length K; each element is a length-d mean vector.

  - `sigma2`: either

    - length-d numeric vector (homoskedastic), or

    - d-by-K numeric matrix (heteroskedastic).

- gamma:

  (optional) n-by-K responsibility matrix.

- data:

  (optional) n-by-d data matrix.

- Q_K:

  K-by-K base precision matrix along components (RW prior along K).
  Should be built with `lambda = 1`. Coordinate penalties are in
  `lambda_vec`.

- lambda_vec:

  length-d nonnegative vector of per-coordinate lambdas.

- rw_q:

  RW order along K (e.g. 1 or 2). Used as rank deficiency in generalized
  logdet.

- ridge:

  Ridge used in building `Q_K` (stored for provenance).

- modelName:

  One of `"homoskedastic"` or `"heteroskedastic"`.

- relative_lambda:

  Logical; if TRUE, scale the prior for coordinate j by \\1/\sigma_j^2\\
  (homoskedastic) or \\1/\bar\sigma_j^2\\ (heteroskedastic; see
  details).

- nugget:

  Nonnegative jitter added to variances after updates.

- eigen_tol:

  Optional tolerance for generalized logdet.

- meta:

  Optional list of metadata.

## Value

An object of class `csmooth_em`.

## Details

For `modelName="heteroskedastic"` and `relative_lambda=TRUE`, we use
\\\bar\sigma_j^2 = \sum_k \pi_k \sigma^2\_{j,k}\\ as a reference scale
for the prior scaling. This is a pragmatic analogue of the EEI-style
scaling.
