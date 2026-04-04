# Convert an algorithm fit object to an `mpcurve` model object

Convert an algorithm fit object to an `mpcurve` model object

## Usage

``` r
as_mpcurve(x, ...)
```

## Arguments

- x:

  A `cavi`, `smooth_em`, or `csmooth_em` object.

- ...:

  Ignored.

## Value

An object of class `"mpcurve"`.

For single-ordering fits, the returned object stores a unified `$params`
block, the responsibility matrix `$gamma`, and a `$locations` field
containing inferred cell locations derived from `$gamma`:

- `$locations$mean$index`: posterior-mean component index in `[1, K]`

- `$locations$mean$pseudotime`: responsibility-weighted pseudotime in
  `[0, 1]`

- `$locations$map$index`: MAP component index from `max.col(gamma)`

- `$locations$map$pseudotime`: pseudotime corresponding to the MAP
  component

For partition fits (`intrinsic_dim >= 2`), `$locations` is a named list
with one such location object per ordering.
