# Subset a csmooth_em object by features (columns)

Internal helper to subset a fitted `csmooth_em` object to a subset of
features. This keeps parameters aligned with a reduced feature set (mu,
sigma2, lambda_vec, etc.).

## Usage

``` r
subset_csmooth_em_fit(fit, keep_cols)
```

## Arguments

- fit:

  A `csmooth_em` object.

- keep_cols:

  Integer vector of feature indices to keep (1-based, in the current
  fit's feature space).

## Value

A `csmooth_em` object restricted to the selected features.
