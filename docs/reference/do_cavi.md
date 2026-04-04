# Continue CAVI sweeps on an existing `cavi` fit

Continue CAVI sweeps on an existing `cavi` fit

## Usage

``` r
do_cavi(
  object,
  iter = 1L,
  tol = NULL,
  lambda = NULL,
  lambda_sd_prior_rate = NULL,
  lambda_min = NULL,
  lambda_max = NULL,
  sigma_min = NULL,
  sigma_max = NULL,
  verbose = FALSE
)
```

## Arguments

- object:

  A `cavi` object.

- iter:

  Integer maximum number of additional CAVI sweeps.

- tol:

  Optional convergence tolerance. If `NULL`, reuse the original fit's
  tolerance.

- lambda:

  Optional scalar or d-vector. If provided, overrides `lambda_vec` in
  the fit and fixes lambda (skips variational lambda updates).

- lambda_sd_prior_rate:

  Optional positive rate for the induced exponential prior on
  `1 / sqrt(lambda_j)`. If `NULL`, reuse the value stored in the
  original fit's `$control`. An explicit `0` is treated as "no penalty"
  for backward compatibility and does not represent a literal
  exponential prior with rate zero.

- lambda_min, lambda_max:

  Optional bounds for `lambda_j`. If `NULL`, reuse the values stored in
  the original fit's `$control`. Ignored when `lambda` is provided.

- sigma_min, sigma_max:

  Optional bounds for `sigma_j^2`. If `NULL`, reuse the values stored in
  the original fit's `$control`.

- verbose:

  Logical.

## Value

An updated `cavi` object with extended traces. The returned object keeps
the full accumulated `$elbo_trace`, `$loglik_trace`, `$lambda_trace`,
`$sigma2_trace`, and `$pi_trace`.
