# Backward greedy feature partition into two orderings (csmoothEM, warm-start only)

Warm-start backward greedy algorithm for partitioning features into two
latent orderings. The procedure initializes ordering 1 once on all
features, initializes ordering 2 once from a single seed feature via a
rank-based 1D seed (no global ordering method), and then performs greedy
moves using only warm-start updates (drop/append + a few EM iterations).

At each step, the algorithm selects the feature in ordering 1 with the
largest positive gain \\\Delta_j = S_2(j \mid \Gamma_2) - S_1(j \mid
\Gamma_1)\\ and moves it to ordering 2, where \\S_k\\ is computed by
[`score_feature_given_Gamma`](https://aguerozz.github.io/MPCurver/reference/score_feature_given_Gamma.md)
under the current responsibilities.

This function never re-calls ordering initializers (e.g.
fiedler/pcurve/PCA) after the initial construction. Ordering updates are
done only via warm-start csmoothEM iterations.

## Usage

``` r
backward_two_ordering_partition_csmooth(
  X,
  K = 50,
  init_method1 = c("fiedler", "PCA", "tSNE", "pcurve", "random"),
  seed2 = NULL,
  score_mode = c("ml", "none"),
  adaptive = NULL,
  rw_q = 2,
  relative_lambda = TRUE,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  discretization = c("quantile", "equal", "kmeans"),
  warm_iter_init = 10,
  warm_iter_refit = 5,
  max_steps = 50,
  verbose = TRUE,
  ...
)
```

## Arguments

- X:

  Numeric matrix `(n x d)`.

- K:

  Integer \\\ge 2\\. Number of mixture components.

- init_method1:

  Ordering initializer for ordering 1, passed to
  [`initialize_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/initialize_csmoothEM.md).

- seed2:

  Optional integer. If NULL, chosen as the worst-aligned feature under
  ordering 1.

- score_mode:

  One of `"ml"` or `"none"`; passed to
  [`score_feature_given_Gamma`](https://aguerozz.github.io/MPCurver/reference/score_feature_given_Gamma.md).

- adaptive:

  Adaptive mode used in
  [`do_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM.md)
  warm updates (`"ml"` or `"none"`). If NULL, defaults to `score_mode`.

- rw_q:

  Integer \\\ge 0\\. RW rank deficiency.

- relative_lambda:

  Logical.

- lambda_min, lambda_max:

  Bounds for lambda optimization when `score_mode="ml"`.

- discretization:

  Discretization used in initialization (recommended: `"quantile"`).

- warm_iter_init:

  Initial warm-start iterations for ordering 1.

- warm_iter_refit:

  Warm-start iterations after each move (applied to both fits).

- max_steps:

  Maximum number of greedy moves.

- verbose:

  Logical.

- ...:

  Passed to
  [`initialize_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/initialize_csmoothEM.md)
  for ordering 1.

## Value

A list with:

- `assign`: character vector length d in `c("A","B")`.

- `J1`, `J2`: feature indices (original X column indices) in each
  ordering.

- `seedB`: seed feature index for ordering 2.

- `fit1`, `fit2`: `csmooth_em` objects for each ordering.

- `history`: data.frame of moves (feature, gain, sizes).

## Details

This is the legacy `csmooth_em`-based backward partition routine. It is
retained mainly for benchmark comparison and internal testing against
the newer variational pipeline. For user-facing analyses, prefer
[`fit_mpcurve`](https://aguerozz.github.io/MPCurver/reference/fit_mpcurve.md)
with the soft-CAVI partition path and optional greedy dimension
selection.

## See also

[`fit_mpcurve`](https://aguerozz.github.io/MPCurver/reference/fit_mpcurve.md),
[`initialize_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/initialize_csmoothEM.md),
[`do_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM.md),
[`score_feature_given_Gamma`](https://aguerozz.github.io/MPCurver/reference/score_feature_given_Gamma.md),
[`append_coord_to_fit_csmooth`](https://aguerozz.github.io/MPCurver/reference/append_coord_to_fit_csmooth.md),
[`drop_coord_from_fit_csmooth`](https://aguerozz.github.io/MPCurver/reference/drop_coord_from_fit_csmooth.md)
