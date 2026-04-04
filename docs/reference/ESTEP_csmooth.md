# Internal E-step for csmooth_em

Internal E-step for csmooth_em

## Usage

``` r
ESTEP_csmooth(X, params, modelName)
```

## Arguments

- X:

  n-by-d data matrix

- params:

  csmooth params list (pi, mu, sigma2)

- modelName:

  "homoskedastic" or "heteroskedastic"

## Value

n-by-K responsibility matrix
