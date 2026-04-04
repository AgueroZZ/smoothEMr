# Initialization from an ordering vector for csmoothEM

Initializes csmoothEM parameters from a 1D ordering score by:

1.  Discretizing `ordering_vec` into `K` ordered groups
    (equal/quantile/kmeans).

2.  Computing component means `mu`.

3.  Estimating diagonal variances `sigma2` under either:

    - `modelName="homoskedastic"`: one variance per coordinate shared
      across clusters (length `d`).

    - `modelName="heteroskedastic"`: one variance per coordinate per
      cluster (`d x K`).

This function mirrors the discretization and mean construction logic of
[`make_init()`](https://aguerozz.github.io/MPCurver/reference/make_init.md),
but returns `sigma2` (diagonal variances) instead of a list of full
covariance matrices.

## Usage

``` r
make_init_csmooth(
  X,
  ordering_vec,
  K,
  modelName = c("homoskedastic", "heteroskedastic"),
  nugget = 0,
  discretization = c("equal", "quantile", "kmeans"),
  na_action = c("drop", "error"),
  eps = 1e-12
)
```

## Arguments

- X:

  Numeric matrix `(n x d)`.

- ordering_vec:

  Numeric vector of length `n` (can contain NA).

- K:

  Integer \>= 2; number of mixture components.

- modelName:

  Either `"homoskedastic"` or `"heteroskedastic"`.

- nugget:

  Nonnegative scalar added to variance estimates.

- discretization:

  One of `"equal"`, `"quantile"`, `"kmeans"`.

- na_action:

  How to handle NA in `ordering_vec`: `"drop"` or `"error"`.

- eps:

  Small positive floor for `pi` and `sigma2`.

## Value

A list with:

- `pi`: length-`K` mixing proportions.

- `mu`: list of length `K`, each a length-`d` mean vector.

- `sigma2`: diagonal variances; length-`d` vector (homoskedastic) or
  `d x K` matrix (heteroskedastic).

- `keep_idx`: row indices kept after NA handling.

- `cluster_rank`: integer vector in `1..K` for each kept row.
