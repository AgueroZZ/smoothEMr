# Summarise an `mpcurve` model fit

Delegates to the underlying `summary.cavi`, `summary.csmooth_em`, or
`summary.smooth_em`, wrapping the result in a `summary.mpcurve` object
so users can access all fields (`$ml_last`, `$elbo_last`, convergence
diagnostics, etc.) that the underlying class exposes.

## Usage

``` r
# S3 method for class 'mpcurve'
summary(object, ...)
```

## Arguments

- object:

  An `mpcurve` object.

- ...:

  Passed to the underlying summary method.

## Value

An object of class `summary.mpcurve` with fields: `$algorithm`,
`$modelName`, `$K`, `$n`, `$d`, and `$underlying` (the full underlying
summary).
