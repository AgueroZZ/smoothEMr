# Krige/interpolate from observed nodes using a Gaussian precision matrix

Conditional mean under a Gaussian model specified by a precision matrix
`Q_new`.

## Usage

``` r
kriging_from_precision(f_obs, idx_obs, Q_new, nugget = 0, keep_obs = TRUE)
```

## Arguments

- f_obs:

  Numeric vector length \|O\| or matrix \|O\| x p.

- idx_obs:

  Integer vector of observed indices (subset of 1:K_new).

- Q_new:

  K_new x K_new precision matrix (dense or sparse).

- nugget:

  Nonnegative scalar: add nugget\*I to Q_UU for stability.

- keep_obs:

  If TRUE, returns full length K_new (or K_new x p) with obs filled in.

## Value

Full vector/matrix if keep_obs=TRUE; else list(pred_idx, f_pred).
