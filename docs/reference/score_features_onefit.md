# Score features under a fitted csmooth_em model via per-coordinate collapsed contributions

Score features under a fitted csmooth_em model via per-coordinate
collapsed contributions

## Usage

``` r
score_features_onefit(fit, X = NULL, include_constant = TRUE)
```

## Arguments

- fit:

  A csmooth_em object.

- X:

  Optional n-by-d data matrix. If NULL, uses fit\$data.

- include_constant:

  Logical; passed to compute_C_by_coord_csmooth().

## Value

A list with - C_coord: length-d vector - global_logpi, global_entropy:
scalars - logdetH_coord: length-d vector
