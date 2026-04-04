# Initialize M CAVI trajectories for multi-ordering partitioning

Initialize M CAVI trajectories for multi-ordering partitioning

## Usage

``` r
init_m_trajectories_cavi(
  X,
  M = 2L,
  methods = NULL,
  pca_components = NULL,
  K = NULL,
  rw_q = 2L,
  ridge = 0,
  lambda_sd_prior_rate = NULL,
  smooth_fit_lambda_mode = c("optimize", "fixed"),
  smooth_fit_lambda_value = 1,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  sigma_min = 1e-10,
  sigma_max = 1e+10,
  discretization = c("quantile", "equal", "kmeans"),
  partition_init = c("similarity", "ordering_methods"),
  similarity_metric = c("spearman", "pearson", "smooth_fit"),
  cluster_linkage = "single",
  similarity_min_feature_sd = 1e-08,
  num_iter = 5L,
  verbose = FALSE
)
```

## Arguments

- X:

  Numeric matrix (`n x d`).

- M:

  Integer \>= 2. Number of orderings.

- methods:

  Character vector of length M. Each element is an ordering method
  (e.g., `"PCA"`, `"fiedler"`, `"random"`).

- pca_components:

  Integer vector of length M. For PCA-based methods, which PC component
  to use. Ignored for non-PCA methods.

- K:

  Number of pseudotime bins.

- rw_q:

  Random-walk order.

- ridge:

  Optional nugget added to the RW precision. Mainly useful for internal
  diagnostics that compare intrinsic and properized priors.

- lambda_sd_prior_rate:

  Optional positive rate for an exponential prior on
  `1 / sqrt(lambda_j^{(m)})`. The default `NULL` means no lambda prior
  penalty. For backward compatibility, an explicit `0` is treated the
  same way; it is only an alias for "no penalty" and does not correspond
  to a literal exponential prior with rate zero.

- smooth_fit_lambda_mode:

  For `similarity_metric = "smooth_fit"` only: whether the directional
  smoother optimizes `lambda` or keeps it fixed at
  `smooth_fit_lambda_value`.

- smooth_fit_lambda_value:

  For `similarity_metric = "smooth_fit"` only: fixed `lambda` value used
  when `smooth_fit_lambda_mode = "fixed"`, and the starting value when
  `smooth_fit_lambda_mode = "optimize"`.

- lambda_min, lambda_max:

  Lambda bounds.

- sigma_min, sigma_max:

  Bounds used by the `"smooth_fit"` similarity metric and the warm-start
  single-ordering CAVI fits.

- discretization:

  Initial discretization scheme used when converting ordering scores
  into `K` bins. For partition fits, MPCurver enforces a common `K`
  across orderings and may fall back to equal-width bins if quantile
  discretization collapses.

- partition_init:

  Either `"similarity"` for feature-similarity clustering followed by
  block-specific ordering initialization or `"ordering_methods"` for the
  existing per-ordering warm starts. The default is `"similarity"`.

- similarity_metric:

  For `partition_init = "similarity"` only: feature-similarity metric
  used to construct the clustering. One of `"spearman"`, `"pearson"`, or
  `"smooth_fit"`. `"smooth_fit"` can be used with `ridge = 0`, but for
  `M > 1` this uses an intrinsic-RW pseudo-evidence rather than a fully
  proper marginal likelihood; use a small positive `ridge` if you want
  the smoother evidence to be theoretically proper.

- cluster_linkage:

  For `partition_init = "similarity"` only: hierarchical-clustering
  linkage applied to `1 - S(X)`.

- similarity_min_feature_sd:

  For `partition_init = "similarity"` only: features below this
  standard-deviation threshold are treated as low-information and get
  zero off-diagonal similarity.

- num_iter:

  Warm-start CAVI sweeps per fit.

- verbose:

  Logical.

## Value

A list with `$fits` (list of M cavi objects), `$init_info` (per-ordering
method metadata), and `$ordering_similarity` (pairwise ordering
correlations for `partition_init = "ordering_methods"`) plus optional
`$similarity_init` metadata for `partition_init = "similarity"`. When
similarity initialization is used, `$similarity_init` stores the full
symmetric similarity matrix `S`, the corresponding `distance = 1 - S`,
and for `similarity_metric = "smooth_fit"` the raw directional score
diagnostics used to build `S`, including directional `lambda` and
`sigma^2` estimates.
