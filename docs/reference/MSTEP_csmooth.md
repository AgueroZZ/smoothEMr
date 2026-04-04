# Internal M-step for csmooth_em (coordinate-wise)

Internal M-step for csmooth_em (coordinate-wise)

## Usage

``` r
MSTEP_csmooth(
  X,
  Gamma,
  params,
  Q_K,
  lambda_vec,
  modelName,
  relative_lambda,
  nugget = 0,
  rw_q = 0L,
  iterate_once = TRUE
)
```

## Arguments

- X:

  n-by-d data

- Gamma:

  n-by-K responsibilities

- params:

  current params (pi, mu, sigma2)

- Q_K:

  K-by-K base precision

- lambda_vec:

  length-d vector

- modelName:

  "homoskedastic" or "heteroskedastic"

- relative_lambda:

  logical

- nugget:

  nonnegative

- rw_q:

  rank deficiency along K (RW order)

- iterate_once:

  logical; if FALSE, you could add inner loops later (kept simple here)

## Value

updated params list (pi, mu, sigma2)
