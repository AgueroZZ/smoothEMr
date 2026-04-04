# Soft dual-trajectory partition using explicit mean-field CAVI

Wrapper around
[`soft_partition_cavi`](https://aguerozz.github.io/MPCurver/reference/soft_partition_cavi.md)
with `M = 2`. Provides the original dual-trajectory interface for
backward compatibility.

## Usage

``` r
soft_two_trajectory_cavi(
  X,
  fit1_init = NULL,
  fit2_init = NULL,
  init_method = c("score", "mincor", "random_split"),
  init_method1 = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
  init_seed = 42L,
  K = NULL,
  discretization = c("quantile", "equal", "kmeans"),
  T_start = 5,
  T_end = 1,
  n_outer = 25L,
  inner_iter = 1L,
  max_converge_iter = 100L,
  tol_outer = 1e-05,
  rw_q = 2L,
  ridge = 0,
  lambda_sd_prior_rate = NULL,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  sigma_min = 1e-10,
  sigma_max = 1e+10,
  assignment_prior = c("uniform", "dirichlet"),
  ordering_alpha = 0.5,
  hard_assign_final = FALSE,
  freeze_unused_ordering = TRUE,
  freeze_unused_ordering_threshold = 0.5,
  freeze_feature = TRUE,
  freeze_feature_weight_threshold = 0.1,
  drop_unused_ordering = FALSE,
  verbose = TRUE
)
```

## Arguments

- X:

  Numeric matrix (`n x d`).

- fit1_init, fit2_init:

  Optional initial `cavi` fits. If NULL, they are initialized via
  [`init_two_trajectories_cavi()`](https://aguerozz.github.io/MPCurver/reference/init_two_trajectories_cavi.md).

- init_method, init_method1:

  Initialization methods.

- init_seed:

  Random seed for `"random_split"`.

- K:

  Number of pseudotime bins.

- discretization:

  Initial discretization scheme used when converting ordering scores
  into `K` bins. Dual-ordering initialisation enforces a common `K`
  across the two orderings.

- T_start, T_end:

  Annealing temperatures.

- n_outer:

  Number of annealing steps.

- inner_iter:

  Number of weighted CAVI sweeps per outer step.

- max_converge_iter:

  Maximum number of exact `T=1` iterations.

- tol_outer:

  Relative objective tolerance for the stricter `T=1` convergence rule
  used after annealing.

- rw_q:

  Random-walk order.

- ridge:

  Optional nugget added to the RW precision for both orderings.

- lambda_sd_prior_rate:

  Optional positive rate for an exponential prior on
  `1 / sqrt(lambda_j^{(m)})`. The default `NULL` means no lambda prior
  penalty. For backward compatibility, an explicit `0` is treated the
  same way; it is only an alias for "no penalty" and does not correspond
  to a literal exponential prior with rate zero.

- lambda_min, lambda_max:

  Bounds passed through to the weighted single-ordering `cavi` updates.

- sigma_min, sigma_max:

  Bounds for the feature-specific noise variance `sigma2_j`. Forwarded
  to internal
  [`cavi()`](https://aguerozz.github.io/MPCurver/reference/cavi.md)
  calls.

- assignment_prior:

  Either `"uniform"` or `"dirichlet"`.

- ordering_alpha:

  Positive scalar used when `assignment_prior = "dirichlet"`.

- hard_assign_final:

  Logical.

- freeze_unused_ordering:

  Logical; if `TRUE`, unused orderings may be frozen and skipped in
  subsequent updates.

- freeze_unused_ordering_threshold:

  Non-negative feature-mass threshold used when
  `freeze_unused_ordering = TRUE`.

- freeze_feature:

  Logical; if `TRUE`, feature-ordering pairs with sufficiently small
  posterior weight are frozen and their trajectory contribution is
  neutralized.

- freeze_feature_weight_threshold:

  Non-negative threshold on `w_{jm}` used when `freeze_feature = TRUE`.

- drop_unused_ordering:

  Logical; if `TRUE`, user-facing wrappers may compact the returned
  partition view to the active orderings after fitting. With `FALSE`,
  the exposed `$objective_history` keeps the fixed-requested-`M`
  comparison semantics of
  [`soft_partition_cavi`](https://aguerozz.github.io/MPCurver/reference/soft_partition_cavi.md);
  with `TRUE`, it becomes the post-drop fitting objective of the active
  model.

- verbose:

  Logical.

## Value

An object of class `"soft_partition_cavi"` with fields: `$fits` (list of
2 cavi fits), `$pi_weights` (d x 2), `$assign`, `$M`,
`$objective_history`, etc. The `$objective_history` semantics match
[`soft_partition_cavi`](https://aguerozz.github.io/MPCurver/reference/soft_partition_cavi.md).
