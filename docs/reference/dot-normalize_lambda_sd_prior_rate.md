# Normalize the optional exponential rate for the lambda SD prior

Public APIs use `NULL` to mean "no lambda prior penalty". For backward
compatibility, an explicit numeric `0` is treated the same way.

## Usage

``` r
.normalize_lambda_sd_prior_rate(rate, arg_name = "lambda_sd_prior_rate")
```

## Arguments

- rate:

  Optional non-negative scalar.

- arg_name:

  Argument name used in error messages.

## Value

Either `NULL` (penalty off) or a strictly positive numeric scalar.
