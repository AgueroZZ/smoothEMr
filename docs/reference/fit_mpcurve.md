# Fit an MPCurve model with the public CAVI interface

High-level user-facing wrapper around MPCurve's recommended CAVI fitting
paths. Use `intrinsic_dim = 1` for a standard single-ordering fit and
`intrinsic_dim >= 2` for the partition-CAVI model with one ordering per
latent dimension.

## Usage

``` r
fit_mpcurve(
  X,
  method = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
  K = NULL,
  rw_q = 2,
  lambda = 1,
  fix_lambda = FALSE,
  iter = 100,
  tol = 1e-06,
  num_cores = 1L,
  intrinsic_dim = 1L,
  greedy = c("none", "forward", "backward"),
  partition_init = c("similarity", "ordering_methods"),
  discretization = NULL,
  ridge = 0,
  lambda_sd_prior_rate = NULL,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  sigma_min = 1e-10,
  sigma_max = 1e+10,
  assignment_prior = c("uniform", "dirichlet"),
  ordering_alpha = 0.5,
  similarity_metric = c("spearman", "pearson", "smooth_fit"),
  smooth_fit_lambda_mode = c("optimize", "fixed"),
  smooth_fit_lambda_value = 1,
  cluster_linkage = "single",
  similarity_min_feature_sd = 1e-08,
  T_start = 5,
  T_end = 1,
  n_outer = 25L,
  inner_iter = 1L,
  max_converge_iter = 100L,
  tol_outer = 1e-05,
  freeze_unused_ordering = TRUE,
  freeze_unused_ordering_threshold = 0.5,
  freeze_feature = TRUE,
  freeze_feature_weight_threshold = 0.1,
  drop_unused_ordering = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- X:

  Numeric matrix (n x d) of observations.

- method:

  Initialisation method(s) for the trajectory ordering. If
  `num_cores > 1` and `length(method) > 1` in the single-ordering case,
  all methods are run in parallel and a named list of fits is returned.
  For partition fits, `length(method)` must be either 1 or
  `intrinsic_dim`.

- K:

  Integer number of mixture components (grid knots).

- rw_q:

  Integer random-walk order for the GMRF prior.

- lambda:

  Positive scalar initial value for `lambda_j` in the single-ordering
  CAVI fit.

- fix_lambda:

  Logical; if `TRUE`, keep `lambda_j` fixed at the supplied initial
  value in the single-ordering CAVI fit.

- iter:

  Maximum number of CAVI sweeps for the single-ordering fit.

- tol:

  Relative ELBO tolerance for the single-ordering CAVI fit.

- num_cores:

  Integer \>= 1. Workers for parallel multi-method single-ordering runs.

- intrinsic_dim:

  Integer intrinsic dimensionality of the latent ordering system.
  `intrinsic_dim = 1` fits the standard single-ordering model. Values
  `>= 2` fit the partition-CAVI model.

- greedy:

  Dimension-selection mode. `"none"` keeps the requested
  `intrinsic_dim`. `"forward"` treats `intrinsic_dim` as an upper bound
  and compares \\M\\ versus \\M+1\\ sequentially starting from \\M=1\\.
  `"backward"` fits the upper bound first, then compares \\M\\ versus
  \\M-1\\ while greedily removing the most correlated active ordering.

- partition_init:

  For partition fits only: either `"similarity"` for
  feature-similarity-driven block initialization or `"ordering_methods"`
  for the existing per-ordering warm starts. The default is
  `"similarity"`.

- discretization:

  Optional discretization method passed to ordering-based
  initialization. For partition fits, MPCurver enforces a common `K`
  across orderings; if quantile cuts collapse, it falls back to
  equal-width bins.

- ridge:

  Optional nugget added to the RW precision.

- lambda_sd_prior_rate:

  Optional positive rate for an exponential prior on
  `1 / sqrt(lambda_j)`. Passed to the CAVI backend and the partition
  CAVI path. The default `NULL` means no lambda prior penalty. For
  backward compatibility, an explicit `0` is treated the same way; it is
  only an alias for "no penalty" and does not correspond to a literal
  exponential prior with rate zero.

- lambda_min, lambda_max:

  Positive bounds for `lambda_j`.

- sigma_min, sigma_max:

  Positive bounds for `sigma_j^2`. These apply to the single-ordering
  CAVI path and to the `"smooth_fit"` similarity metric when
  `partition_init = "similarity"`.

- assignment_prior:

  For partition fits only: either `"uniform"` or `"dirichlet"`.

- ordering_alpha:

  For partition fits only: positive scalar concentration used when
  `assignment_prior = "dirichlet"`.

- similarity_metric:

  For `partition_init = "similarity"` only: feature-similarity metric
  used to construct feature blocks. One of `"spearman"`, `"pearson"`, or
  `"smooth_fit"`. The `"smooth_fit"` metric can be used with
  `ridge = 0`, but for `intrinsic_dim > 1` this uses an intrinsic-RW
  pseudo-evidence rather than a fully proper marginal likelihood; use a
  small positive `ridge` if you want the smoother evidence to be
  theoretically proper.

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
  linkage applied to `1 - S(X)`. Defaults to `"single"`.

- similarity_min_feature_sd:

  For `partition_init = "similarity"` only: low-variance feature
  threshold used when building `S(X)`.

- T_start, T_end:

  Annealing temperatures for partition CAVI.

- n_outer:

  Number of annealing steps for partition CAVI.

- inner_iter:

  Number of weighted CAVI sweeps per annealing step.

- max_converge_iter:

  Maximum number of exact `T = 1` partition iterations after annealing.

- tol_outer:

  Relative objective tolerance for partition phase-2 convergence.

- freeze_unused_ordering:

  For partition fits only: if `TRUE`, inactive orderings are frozen and
  skipped in subsequent weighted updates.

- freeze_unused_ordering_threshold:

  For partition fits only: non-negative feature-mass threshold used to
  decide when an ordering is effectively unused.

- freeze_feature:

  For partition fits only: if `TRUE`, feature-ordering pairs with
  sufficiently small posterior weight are frozen and their trajectory
  contribution is neutralized.

- freeze_feature_weight_threshold:

  For partition fits only: non-negative threshold on `w_{jm}` used when
  `freeze_feature = TRUE`.

- drop_unused_ordering:

  For partition fits only: if `TRUE`, the returned `mpcurve` view may
  compact away frozen orderings. With `FALSE`, the reported partition
  `$objective_history` is the fixed-requested-`M` comparison objective;
  with `TRUE`, it is the post-drop fitting objective of the active
  model.

- verbose:

  Logical; print per-iteration progress?

- ...:

  Advanced CAVI / partition-CAVI options forwarded to the underlying
  backend. This is intended for advanced initialization controls such as
  `responsibilities_init`, `pi_init`, `sigma2_init`, `fits_init`,
  `init_methods`, `pca_components`, or `hard_assign_final`. Legacy
  smoothEM/csmoothEM controls are rejected.

## Value

An [`mpcurve`](https://aguerozz.github.io/MPCurver/reference/mpcurve.md)
object, or (for parallel multi-method single-ordering runs) a named list
of `mpcurve` objects with a `summary` attribute. Every returned fit
stores inferred cell locations in `$locations`; for single-ordering fits
this includes both posterior-mean and MAP locations derived from
`$gamma`. Greedy runs additionally attach a `$greedy_selection` metadata
block recording the stepwise search history and stop reason.

## Details

Legacy `smooth_em` and `csmooth_em` algorithms remain available as
lower-level compatibility functions, but they are no longer part of the
public `fit_mpcurve()` interface and are not the active development path
for the package.
