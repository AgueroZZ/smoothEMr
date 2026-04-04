# Partition features by comparing CAVI feature scores from two fits

Partition features by comparing CAVI feature scores from two fits

## Usage

``` r
partition_features_twofits_cavi(
  fitA,
  fitB,
  X = NULL,
  delta = 0,
  include_prior = TRUE
)
```

## Arguments

- fitA, fitB:

  `cavi` objects.

- X:

  Optional data matrix.

- delta:

  Nonnegative assignment margin.

- include_prior:

  Logical.

## Value

A list with `assign`, `score_diff`, `CA`, and `CB`.
