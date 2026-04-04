# Initialize ordering for csmoothEM (diagonal-variance version)

Same ordering logic as
[`initialize_ordering()`](https://aguerozz.github.io/MPCurver/reference/initialize_ordering.md),
but returns csmooth-style parameters with diagonal variances `sigma2`.

## Usage

``` r
initialize_ordering_csmooth(
  X,
  K,
  method = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
  discretization = c("equal", "quantile", "kmeans"),
  modelName = c("homoskedastic", "heteroskedastic"),
  nugget = 0,
  eps = 1e-12,
  ...
)
```

## Arguments

- X:

  Numeric matrix (n x d).

- K:

  Integer \>= 2; number of mixture components.

- method:

  One of `"PCA","fiedler","pcurve","tSNE","random", "isomap"`.

- discretization:

  One of `"equal","quantile","kmeans"`.

- modelName:

  Either `"homoskedastic"` or `"heteroskedastic"`.

- nugget:

  Nonnegative scalar added to variance estimates.

- eps:

  Small positive floor for `pi` and `sigma2`.

- ...:

  Extra arguments passed to the ordering method
  (PCA/tSNE/pcurve/fiedler).

## Value

A list with fields:

- `params`: list(pi, mu, sigma2)

- `keep_idx`: indices kept after NA handling

- `ordering`: ordering metadata including score `t` and method name
