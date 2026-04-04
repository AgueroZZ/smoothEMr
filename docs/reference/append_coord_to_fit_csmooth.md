# Append one coordinate to a csmooth_em fit using a Gamma-based 1D initialization (internal)

Adds a new feature column to a `csmooth_em` object, initializing the new
coordinate's \\\mu\_{j\cdot}\\ and \\\sigma_j^2\\ using
`score_feature_given_Gamma` under the fit's current responsibilities
`fit$gamma`. Then appends the coordinate and sets its lambda to 1 (it
can be updated later by `do_csmoothEM` if adaptive is on).

This is a warm-start operation: it does not reinitialize the whole
model.

## Usage

``` r
append_coord_to_fit_csmooth(
  fit,
  xj,
  score_mode = c("ml", "none"),
  rw_q = 2,
  relative_lambda = TRUE,
  lambda_min = 1e-10,
  lambda_max = 1e+10
)
```

## Arguments

- fit:

  A `csmooth_em` object (homoskedastic).

- xj:

  Numeric vector length n (the new feature column).

- score_mode:

  One of `"ml"` or `"none"`; controls how the 1D init is computed.

- rw_q:

  Integer RW rank deficiency along K.

- relative_lambda:

  Logical.

- lambda_min, lambda_max:

  Bounds for lambda optimization when `score_mode="ml"`.

## Value

Updated `csmooth_em`.
