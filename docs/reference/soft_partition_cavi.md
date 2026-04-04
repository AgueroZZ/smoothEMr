# Soft M-ordering partition using mean-field CAVI

Generalised soft partition routine supporting \\M \ge 2\\ orderings.
Features are assigned probabilistic weights across orderings via
temperature-annealed softmax on per-feature ELBO scores.

## Usage

``` r
soft_partition_cavi(
  X,
  M = 2L,
  fits_init = NULL,
  init_methods = NULL,
  pca_components = NULL,
  partition_init = c("similarity", "ordering_methods"),
  similarity_metric = c("spearman", "pearson", "smooth_fit"),
  smooth_fit_lambda_mode = c("optimize", "fixed"),
  smooth_fit_lambda_value = 1,
  cluster_linkage = "single",
  similarity_min_feature_sd = 1e-08,
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

- M:

  Integer \>= 2. Number of orderings.

- fits_init:

  Optional list of M `cavi` fits. If NULL, initialized via
  [`init_m_trajectories_cavi()`](https://aguerozz.github.io/MPCurver/reference/init_m_trajectories_cavi.md).

- init_methods:

  Character vector of length M. Ordering initialisation methods (e.g.,
  `c("fiedler", "PCA", "PCA")`). If NULL, defaults to `rep("PCA", M)`.

- pca_components:

  Integer vector of length M. Which PCA component to use for PCA-based
  methods. If NULL, auto-assigned sequentially.

- partition_init:

  Partition initialisation strategy. Use `"similarity"` to cluster
  features into `M` blocks via `1 - S(X)` before building
  ordering-specific warm starts, or `"ordering_methods"` for the
  existing per-ordering method starts. The default is `"similarity"`.

- similarity_metric:

  For `partition_init = "similarity"` only: pairwise feature similarity
  metric used to form `S(X)`. One of `"spearman"`, `"pearson"`, or
  `"smooth_fit"`. `"smooth_fit"` can be used with `ridge = 0`, but for
  `M > 1` this uses an intrinsic-RW pseudo-evidence rather than a fully
  proper marginal likelihood; use a small positive `ridge` if you want
  the smoother evidence to be theoretically proper.

- smooth_fit_lambda_mode:

  For `similarity_metric = "smooth_fit"` only: whether the directional
  smoother optimizes `lambda` or keeps it fixed at
  `smooth_fit_lambda_value`.

- smooth_fit_lambda_value:

  For `similarity_metric = "smooth_fit"` only: fixed `lambda` value used
  when `smooth_fit_lambda_mode = "fixed"`, and the starting value when
  `smooth_fit_lambda_mode = "optimize"`.

- cluster_linkage:

  For `partition_init = "similarity"` only: hierarchical-clustering
  linkage passed to [`hclust`](https://rdrr.io/r/stats/hclust.html). The
  default is `"single"`.

- similarity_min_feature_sd:

  For `partition_init = "similarity"` only: features with standard
  deviation below this threshold are treated as low-information and
  assigned zero off-diagonal similarity.

- K:

  Number of pseudotime bins.

- discretization:

  Initial discretization scheme used when converting ordering scores
  into `K` bins. Partition fits require a common `K` across orderings,
  so if quantile discretization collapses the number of groups, MPCurver
  falls back to equal-width bins.

- T_start, T_end:

  Annealing temperatures.

- n_outer:

  Number of annealing steps.

- inner_iter:

  Number of weighted CAVI sweeps per outer step.

- max_converge_iter:

  Maximum number of post-annealing iterations at `T = T_end`.

- tol_outer:

  Relative objective tolerance for the phase-2 convergence rule.
  Convergence is declared only after at least three exact `T = T_end`
  iterations and three consecutive nonnegative relative objective
  improvements below this threshold.

- rw_q:

  Random-walk order.

- ridge:

  Optional nugget added to the RW precision for all orderings.

- lambda_sd_prior_rate:

  Optional positive rate for an exponential prior on
  `1 / sqrt(lambda_j^{(m)})`. The default `NULL` means no lambda prior
  penalty. For backward compatibility, an explicit `0` is treated the
  same way; it is only an alias for "no penalty" and does not correspond
  to a literal exponential prior with rate zero.

- lambda_min, lambda_max:

  Lambda bounds.

- sigma_min, sigma_max:

  Bounds for `sigma_j^2`. These are used both in the weighted CAVI
  updates and by the `"smooth_fit"` similarity metric when
  `partition_init = "similarity"`.

- assignment_prior:

  Either `"uniform"` or `"dirichlet"` for the partition assignment
  prior.

- ordering_alpha:

  Positive scalar concentration parameter used when
  `assignment_prior = "dirichlet"`.

- hard_assign_final:

  Logical; make hard assignments at the end?

- freeze_unused_ordering:

  Logical; if `TRUE`, orderings whose feature-posterior mass falls below
  `freeze_unused_ordering_threshold` are frozen and no longer updated.

- freeze_unused_ordering_threshold:

  Non-negative threshold on the summed feature-posterior mass
  `sum(pi_weights[, m])` used to decide whether an ordering is
  effectively unused.

- freeze_feature:

  Logical; if `TRUE`, feature-ordering pairs with posterior weight below
  `freeze_feature_weight_threshold` are frozen, their trajectory
  `U`-block is no longer updated, and their feature-trajectory
  contribution is neutralized to zero.

- freeze_feature_weight_threshold:

  Non-negative threshold on `w_{jm}` used to decide whether a
  feature-ordering pair is effectively unused.

- drop_unused_ordering:

  Logical; if `TRUE`, the user-facing `mpcurve` wrapper returned by
  [`fit_mpcurve()`](https://aguerozz.github.io/MPCurver/reference/fit_mpcurve.md)
  drops frozen orderings from its displayed partition view. Internally,
  fitting still proceeds by freezing rather than deleting orderings.
  When `drop_unused_ordering = FALSE`, the exposed `$objective_history`
  is a fixed-requested-`M` comparison objective: frozen orderings keep
  their preserved feature-assignment weights and cell-ordering block,
  while their trajectory `U`-block contribution is neutralized to zero.
  When `drop_unused_ordering = TRUE`, the exposed `$objective_history`
  is the post-drop fitting objective of the active model and should not
  be compared across requested values of `M`.

- verbose:

  Logical.

## Value

An object of class `"soft_partition_cavi"` with fields: `$fits` (list of
M cavi fits), `$pi_weights` (d x M matrix), `$assign` (character
vector), `$M`, `$objective_history`, `$score_history`,
`$weight_history`, `$T_schedule`, `$n_anneal`, `$converged`,
`$convergence_info`, `$active_orderings`, `$active_feature_pairs`,
`$frozen_orderings`, `$ordering_events`, `$feature_events`, and
`$control`, plus optional initialisation diagnostics `$init_info` and
`$ordering_similarity`, and optional similarity metadata
`$similarity_init`. When similarity initialisation is used,
`$similarity_init` stores the full symmetric similarity matrix `S`,
`distance = 1 - S`, and for `similarity_metric = "smooth_fit"` the raw
directional score diagnostics used to construct `S`, including
directional `lambda` and `sigma^2` estimates. The exposed
`$objective_history` has two user-facing semantics. If
`drop_unused_ordering = FALSE`, it is the fixed-requested-`M` comparison
objective intended for comparing fits across requested `M`. If
`drop_unused_ordering = TRUE`, it is the post-drop fitting objective of
the active model and is not intended for cross-`M` comparison.
