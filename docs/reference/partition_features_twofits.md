# Partition features by comparing per-coordinate collapsed scores from two fits

Partition features by comparing per-coordinate collapsed scores from two
fits

## Usage

``` r
partition_features_twofits(
  fitA,
  fitB,
  X = NULL,
  delta = 0,
  include_constant = TRUE
)
```

## Arguments

- fitA, fitB:

  csmooth_em objects (e.g., two different orderings / initializations).

- X:

  Optional data matrix. If NULL, uses fitA\$data (and assumes same X for
  both fits).

- delta:

  Nonnegative margin. Assign to A only if (CA - CB) \> delta; otherwise
  assign to B.

- include_constant:

  Logical; passed to score_features_onefit().

## Value

A list with - assign: length-d character vector taking values `"A"` or
`"B"` - score_diff: length-d vector (CA - CB) - CA, CB: length-d vectors
of per-feature scores
