# Compute per-coordinate contributions to the collapsed objective C for csmooth_em

Returns a length-d vector C_j such that C_total = sum_j C_j +
entropy(Gamma) + sum_k Nk \* log(pi_k), where C_total corresponds to the
collapsed objective recorded as ml_trace (up to constants consistent
with include_constant).

## Usage

``` r
compute_C_by_coord_csmooth(
  X,
  Gamma,
  params,
  Q_K,
  lambda_vec,
  modelName,
  relative_lambda,
  eigen_tol = NULL,
  rw_q = 0L,
  include_constant = TRUE
)
```

## Arguments

- X:

  n-by-d data matrix

- Gamma:

  n-by-K responsibilities

- params:

  list(pi, mu, sigma2)

- Q_K:

  base RW precision (K x K)

- lambda_vec:

  length-d vector

- modelName:

  "homoskedastic" or "heteroskedastic"

- relative_lambda:

  logical

- eigen_tol:

  numeric; passed to generalized_logdet

- rw_q:

  integer; rank deficiency of RW precision

- include_constant:

  logical; include +(K/2)log(2pi) per coordinate if TRUE

## Value

a list with - C_coord: length-d vector of per-coordinate collapsed
contributions - global_logpi: scalar sum_k Nk log(pi_k) -
global_entropy: scalar entropy term -sum Gamma log Gamma -
logdetH_coord: length-d vector log\|A_j\| (for diagnostics)
