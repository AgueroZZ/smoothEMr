# Weighted E-step for csmooth_em

Identical to ESTEP_csmooth() except per-feature log-likelihood
contributions are multiplied by feature_weights\[j\] before summing over
features. In the adaptive = "ml" path, this E-step is paired with a
weighted collapsed hyper-step for lambda / sigma2.

## Usage

``` r
ESTEP_csmooth_weighted(X, params, modelName, feature_weights)
```

## Arguments

- X:

  n-by-d data matrix.

- params:

  csmooth params list (pi, mu, sigma2).

- modelName:

  "homoskedastic" or "heteroskedastic".

- feature_weights:

  length-d non-negative weight vector.

## Value

n-by-K responsibility matrix.
