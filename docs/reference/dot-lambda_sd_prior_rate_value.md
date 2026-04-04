# Numeric value of the optional lambda SD prior rate

Converts the public `NULL`-means-off convention into a numeric value
used internally by ELBO bookkeeping helpers.

## Usage

``` r
.lambda_sd_prior_rate_value(rate)
```

## Arguments

- rate:

  Optional rate.

## Value

Zero when the penalty is off, otherwise the positive numeric rate.
