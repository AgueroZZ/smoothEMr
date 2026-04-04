# Score features under a fitted `cavi` model

Score features under a fitted `cavi` model

## Usage

``` r
score_features_onefit_cavi(fit, X = NULL, include_prior = TRUE)
```

## Arguments

- fit:

  A `cavi` object.

- X:

  Optional data matrix. If NULL, uses `fit$data`.

- include_prior:

  Logical; include the GMRF prior contribution.

## Value

A list with `feature_score`, `global_terms`, and `entropy_u`.
