# Drop one coordinate from a csmooth_em fit (internal)

Removes a single feature column from a `csmooth_em` object, updating
`data`, `params`, and `prior$lambda_vec`. This is a warm-start operation
and does not refit the model.

## Usage

``` r
drop_coord_from_fit_csmooth(fit, pos)
```

## Arguments

- fit:

  A `csmooth_em` object.

- pos:

  Integer in `1:d_sub`, position of the coordinate to remove in the
  CURRENT fit.

## Value

Updated `csmooth_em`.
