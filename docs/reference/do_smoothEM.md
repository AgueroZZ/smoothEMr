# Run SmoothEM for a given number of iterations on a smooth_em object

Run SmoothEM for a given number of iterations on a smooth_em object

## Usage

``` r
do_smoothEM(
  object,
  data = NULL,
  iter = 1,
  record = TRUE,
  check_decrease = TRUE,
  tol_decrease = 1e-10,
  adaptive = TRUE,
  lambda_min = NULL,
  lambda_max = NULL,
  verbose = FALSE
)
```

## Arguments

- object:

  A `smooth_em` object created by
  [`as_smooth_em()`](https://aguerozz.github.io/MPCurver/reference/as_smooth_em.md).

- data:

  Numeric matrix (n x d).

- iter:

  Integer \>= 1; number of (E-step + M-step) iterations to run.

- record:

  Logical; whether to append objective values to traces.

- check_decrease:

  Logical; if TRUE, rollback if ELBO decreases materially.

- tol_decrease:

  Numeric; tolerance for considering ELBO decrease (default 1e-10).

- adaptive:

  Logical; if TRUE, update `lambda` each iteration (profile-style
  update).

- lambda_min, lambda_max:

  Bounds for adaptive lambda (ignored if adaptive=FALSE).

- verbose:

  Logical.

## Value

Updated `smooth_em` object.
