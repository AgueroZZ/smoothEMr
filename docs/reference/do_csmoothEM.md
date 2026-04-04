# Run csmoothEM iterations

Runs EM iterations for coordinate-specific SmoothEM (`csmooth_em`).
Supports adaptive updates of the coordinate-wise smoothing parameters
\\\lambda_j\\.

When `adaptive = "ml"`, this function dispatches to
[`do_csmoothEM_ml_collapsed`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_ml_collapsed.md)
(collapsed-ML variant). In that case, `sigma_update` controls whether
\\\sigma^2\\ is updated via the standard M-step or via collapsed-ML
optimization.

## Usage

``` r
do_csmoothEM(
  object,
  data = NULL,
  iter = 1,
  record = TRUE,
  adaptive = NULL,
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

  Integer \\\ge 1\\; number of EM iterations to run.

- record:

  Logical; if TRUE, record traces.

- adaptive:

  Adaptive option controlling \\\lambda\\ updates:

  `NULL`

  :   Inherit `object$control$adaptive`.

  `"none"`

  :   No adaptive update (fixed `lambda_vec`).

  `"prior"`

  :   Obsolete. Retained only for backward compatibility; updates
      \\\lambda_j\\ from the current M-step \\\mu\\ via a prior/profile
      rule and should not be used for convergence claims.

  `"ml"`

  :   Collapsed-ML update for \\\lambda_j\\ (dispatches to
      [`do_csmoothEM_ml_collapsed`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_ml_collapsed.md)).

  Logical values are accepted for backward compatibility: `TRUE` is
  interpreted as `"ml"` and `FALSE` as `"none"`.

- lambda_min, lambda_max:

  Positive bounds used for adaptive \\\lambda\\ updates.

- sigma_update:

  Character. Only used when `adaptive="ml"` (i.e., in the collapsed-ML
  path). Defaults to `"ml"`; `"mstep"` is retained only for backward
  compatibility; see
  [`do_csmoothEM_ml_collapsed`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_ml_collapsed.md).

- sigma_min, sigma_max:

  Positive bounds for `sigma2` when `adaptive="ml"` and
  `sigma_update="ml"`.

- verbose:

  Logical; if TRUE, print a short progress line each iteration.

## Value

Updated `csmooth_em` object.
