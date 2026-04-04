# Cache covariance inverses and log-determinants

Adds `invSigma` and `logdet` to a parameter list.

## Usage

``` r
init_cov_cache_fast(params, jitter = 0)
```

## Arguments

- params:

  A list with at least `pi` and `sigma`, where `sigma` is a list of
  covariance matrices.

- jitter:

  Nonnegative diagonal jitter added if Cholesky fails.

## Value

Updated params list with `invSigma` and `logdet`.
