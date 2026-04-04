# Cache inverse-sigma2 and log-determinant for csmooth_em (diagonal covariance)

For `csmooth_em` the covariance is diagonal, so no Cholesky is needed.
Stores `invsig2` and `logdet` into `params` so that
[`ESTEP_csmooth()`](https://aguerozz.github.io/MPCurver/reference/ESTEP_csmooth.md)
can skip recomputing them on every call.

## Usage

``` r
cache_csmooth_params(params, modelName = c("homoskedastic", "heteroskedastic"))
```

## Arguments

- params:

  A `csmooth_em` parameter list with a `sigma2` field.

- modelName:

  One of `"homoskedastic"` or `"heteroskedastic"`.

## Value

Updated `params` with `$invsig2` and `$logdet`.
