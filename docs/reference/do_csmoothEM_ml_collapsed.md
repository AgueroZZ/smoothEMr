# Run collapsed-ML csmoothEM iterations (homoskedastic only)

Iteration order (collapsed-ML): 1) E-step: update responsibilities using
current (MAP) \\\mu\\. 2) Hyper-step: update \\\pi\\, \\\lambda\\, and
optionally \\\sigma^2\\. 3) M-step: update MAP \\\mu\\ given current
responsibilities and hyperparameters. 4) Record: penalized observed
objective (`loglik_trace`), penalized ELBO (`elbo_trace`), collapsed
objective \\\mathcal{C}\\ (`ml_trace`), and `logdetH_trace`.

## Usage

``` r
do_csmoothEM_ml_collapsed(
  object,
  data = NULL,
  iter = 1,
  record = TRUE,
  lambda_min = NULL,
  lambda_max = NULL,
  sigma_update = c("ml", "mstep"),
  sigma_min = 1e-10,
  sigma_max = 1e+10,
  verbose = FALSE
)
```

## Arguments

- object:

  A `csmooth_em` object.

- data:

  Optional n-by-d data matrix. If NULL, uses `object$data`.

- iter:

  Integer \>= 1 number of iterations.

- record:

  Logical; if TRUE, record traces.

- lambda_min, lambda_max:

  Bounds for lambda optimization.

- sigma_update:

  Either `"ml"` (default; optimize the collapsed objective jointly with
  \\\lambda\\) or legacy `"mstep"` (standard variance M-step).

- sigma_min, sigma_max:

  Bounds for sigma2 when `sigma_update="ml"`.

- verbose:

  Logical.

## Value

Updated `csmooth_em` object.

## Details

The recorded `ml_trace` is the collapsed objective \$\$\mathcal{C} \\=\\
\text{penELBO} \\+\\ \frac{dK}{2}\log(2\pi) \\-\\ \frac12 \log\|H\|,\$\$
where \\H\\ is the Hessian/precision w.r.t. \\\mu\\ (block-diagonal over
coordinates).
