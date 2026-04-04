# Continue CAVI iterations on an existing `mpcurve` fit

Public continuation wrapper for MPCurve's recommended CAVI paths. For
single-ordering fits, this delegates to
[`do_cavi`](https://aguerozz.github.io/MPCurver/reference/do_cavi.md).
For partition fits, it continues the exact \\T = 1\\ phase of
[`soft_partition_cavi`](https://aguerozz.github.io/MPCurver/reference/soft_partition_cavi.md)
using the stored frozen-ordering and assignment-prior controls.

## Usage

``` r
do_mpcurve(
  object,
  iter = 1,
  lambda = NULL,
  tol = NULL,
  lambda_sd_prior_rate = NULL,
  lambda_min = NULL,
  lambda_max = NULL,
  sigma_min = NULL,
  sigma_max = NULL,
  tol_outer = NULL,
  freeze_unused_ordering = NULL,
  freeze_unused_ordering_threshold = NULL,
  freeze_feature = NULL,
  freeze_feature_weight_threshold = NULL,
  drop_unused_ordering = NULL,
  assignment_prior = NULL,
  ordering_alpha = NULL,
  verbose = FALSE
)
```

## Arguments

- object:

  An `mpcurve` object wrapping a `cavi` or `soft_partition_cavi` fit.

- iter:

  Integer \>= 1. Maximum number of additional CAVI sweeps (single
  ordering) or exact `T = 1` partition steps.

- lambda:

  Optional scalar or d-vector. For single-ordering CAVI only: fix
  `lambda_j` at this value for the continuation run.

- tol:

  Optional relative ELBO tolerance for single-ordering CAVI. If `NULL`,
  reuse the stored value.

- lambda_sd_prior_rate:

  Optional positive rate for the induced exponential prior on
  `1 / sqrt(lambda_j)`. If `NULL`, the stored control value is reused.
  An explicit `0` is treated as "no penalty" for backward compatibility
  and does not represent a literal exponential prior with rate zero.

- lambda_min, lambda_max:

  Optional positive bounds for `lambda_j`. If `NULL`, reuse the stored
  values.

- sigma_min, sigma_max:

  Optional positive bounds for `sigma_j^2`. If `NULL`, reuse the stored
  values.

- tol_outer:

  For partition fits only: relative objective tolerance used by the
  phase-2 convergence rule. If `NULL`, reuse the stored value.

- freeze_unused_ordering:

  For partition fits only: if `TRUE`, inactive orderings remain frozen
  and skipped in subsequent weighted updates. If `NULL`, reuse the
  stored value.

- freeze_unused_ordering_threshold:

  For partition fits only: non-negative feature-mass threshold used to
  decide when an ordering is effectively unused. If `NULL`, reuse the
  stored value.

- freeze_feature:

  For partition fits only: if `TRUE`, feature-ordering pairs with
  sufficiently small posterior weight are frozen and their trajectory
  contribution is neutralized. If `NULL`, reuse the stored value.

- freeze_feature_weight_threshold:

  For partition fits only: non-negative threshold on `w_{jm}` used when
  `freeze_feature = TRUE`. If `NULL`, reuse the stored value.

- drop_unused_ordering:

  For partition fits only: if `TRUE`, the returned `mpcurve` view may
  compact away frozen orderings. If `NULL`, reuse the stored value. When
  `FALSE`, the reported partition `$objective_history` retains
  fixed-requested-`M` comparison semantics; when `TRUE`, it becomes the
  post-drop fitting objective of the active model.

- assignment_prior:

  For partition fits only: either `"uniform"` or `"dirichlet"`. If
  `NULL`, reuse the stored value.

- ordering_alpha:

  For partition fits only: positive scalar concentration used when
  `assignment_prior = "dirichlet"`. If `NULL`, reuse the stored value.

- verbose:

  Logical. Print per-iteration progress?

## Value

An updated `mpcurve` object.

## Details

Legacy smoothEM/csmoothEM continuation remains available through the
lower level
[`do_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_smoothEM.md)
and
[`do_csmoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM.md)
functions, not through this high-level wrapper.
