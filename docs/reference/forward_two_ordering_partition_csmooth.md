# Forward greedy feature partition into two orderings (csmoothEM version)

Greedily partitions features (columns of `X`) into two groups, each
associated with its own latent ordering (represented by responsibilities
\\\Gamma\\ over `K` components).

The algorithm mirrors the structure of `two_ordering_smoothEM_v2`:

1.  Choose a seed feature `j1` and fit a 1D csmoothEM model on `X[, j1]`
    to obtain \\\Gamma_1\\ and parameters \\\theta_1\\.

2.  Choose a second seed feature `j2` and fit a 1D csmoothEM model on
    `X[, j2]` to obtain \\\Gamma_2\\ and parameters \\\theta_2\\.

3.  While unassigned features remain:

    - Score each remaining feature under \\\Gamma_1\\ and \\\Gamma_2\\.

    - Select the feature with the largest absolute score gap and assign
      it to the better ordering.

    - Append the 1D fitted object to that ordering's parameters and
      update \\\Gamma\\ by an E-step.

    - Optionally run a short csmoothEM refinement on that ordering
      (controlled by `greedy_em_refine`).

Seeding:

- If `greedy_start_index` is provided, it is used as `j1`.

- Otherwise, `seeding="variance"` chooses `j1` as the maximum-variance
  feature.

- Otherwise, `seeding="correlation"` chooses `j1` as the feature most
  correlated (by default, absolute correlation) with an ordering vector
  `seed_ordering_vec`. If `seed_ordering_vec` is NULL, the ordering
  vector defaults to the first PC score of `X`.

For `j2`:

- If `greedy_start_index2` is provided, it is used as `j2`.

- Else if `seeding="correlation"` and `greedy_start_index` is NULL,
  choose `j2` as the feature least correlated with the same ordering
  vector used for `j1`.

- Else (default), choose `j2` as the feature with the worst alignment
  score under \\\Gamma_1\\.

Scoring:

- `score_mode="none"` uses a plug-in (ELBO/Q-like) score with fixed
  \\\lambda=1\\.

- `score_mode="ml"` uses a collapsed (marginal-like) score with
  \\\lambda\\ optimized per feature.

## Usage

``` r
forward_two_ordering_partition_csmooth(
  X,
  K = 30,
  greedy_start_index = NULL,
  greedy_start_index2 = NULL,
  seeding = c("variance", "correlation"),
  seed_ordering_vec = NULL,
  seed_pc_scale = TRUE,
  seed_cor_abs = TRUE,
  score_mode = c("ml", "none"),
  rw_q = 2,
  relative_lambda = TRUE,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  greedy_em_refine = TRUE,
  greedy_em_max_iter = 10,
  discretization = c("quantile", "equal", "kmeans"),
  verbose = 1
)
```

## Arguments

- X:

  Numeric matrix `(n x d)`.

- K:

  Integer \\\ge 2\\. Number of mixture components.

- greedy_start_index:

  Optional integer. First seed feature. If NULL, chosen by `seeding`.

- greedy_start_index2:

  Optional integer. Second seed feature. If NULL, chosen by `seeding`
  rule above.

- seeding:

  Seeding strategy when `greedy_start_index` is NULL.

  `"variance"`

  :   Use the maximum-variance feature as the first seed (default).

  `"correlation"`

  :   Use correlation with an ordering vector to choose seeds.

- seed_ordering_vec:

  Optional numeric vector of length `n`. Only used when
  `seeding="correlation"` and `greedy_start_index` is NULL. If NULL,
  defaults to the first PC score of `X`.

- seed_pc_scale:

  Logical; only used when `seed_ordering_vec` is NULL and
  `seeding="correlation"`. Passed as `scale.` to `prcomp`.

- seed_cor_abs:

  Logical; if TRUE (default) seed using absolute correlations. If FALSE,
  use signed correlations (max/min correlation).

- score_mode:

  One of `"ml"` or `"none"`.

- rw_q:

  Integer \\\ge 0\\. Rank deficiency along K (RW order).

- relative_lambda:

  Logical; whether relative-lambda scaling is used.

- lambda_min, lambda_max:

  Bounds for lambda optimization when `score_mode="ml"`.

- greedy_em_refine:

  Logical; if TRUE, run short csmoothEM refinement after each append.

- greedy_em_max_iter:

  Integer \\\ge 0\\. Number of refinement iterations per append. Set to
  0 to disable refinement even if `greedy_em_refine=TRUE`.

- discretization:

  Discretization method used in seed initialization (recommended:
  `"quantile"`).

- verbose:

  Integer. `0`=silent, `1`=progress messages.

## Value

A list with:

- `coord_assign`: integer vector length d with values 1 or 2.

- `J1`, `J2`: feature indices in each ordering.

- `Gamma1`, `Gamma2`: responsibilities after greedy completion.

- `params1`, `params2`: csmooth-style params for each ordering.

- `fit1`, `fit2`: `csmooth_em` objects for each ordering (ready for
  further refinement).

- `seeds`: list with `j1`, `j2`.

## Details

This is the legacy `csmooth_em`-based hard-partition routine. It is
retained for benchmark comparison and internal regression testing. For
new work, prefer
[`fit_mpcurve`](https://aguerozz.github.io/MPCurver/reference/fit_mpcurve.md)
with the soft-CAVI partition path and optional greedy dimension
selection.

## See also

[`fit_mpcurve`](https://aguerozz.github.io/MPCurver/reference/fit_mpcurve.md),
[`do_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM.md),
[`score_feature_given_Gamma`](https://aguerozz.github.io/MPCurver/reference/score_feature_given_Gamma.md),
[`compute_C_by_coord_csmooth`](https://aguerozz.github.io/MPCurver/reference/compute_C_by_coord_csmooth.md)
