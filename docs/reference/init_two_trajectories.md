# Initialise two trajectory fits guaranteed to face different directions

Wrapper that produces a pair of `csmooth_em` objects suitable as
`fit1_init` / `fit2_init` for
[`soft_two_trajectory_EM()`](https://aguerozz.github.io/MPCurver/reference/soft_two_trajectory_EM.md).
Three strategies are available:

- **"score"**:

  fit1 is initialised with PCA on all features. fit2 is seeded from the
  single feature that fit1 explains least (lowest per-feature marginal
  score).

- **"mincor"**:

  fit1 is initialised via `method1` (default PCA). fit2 is seeded from
  the feature with the smallest absolute correlation with the feature
  that drives fit1 most (top PC1 loading). Ensures fit2 starts from a
  direction maximally orthogonal to fit1.

- **"random_split"**:

  Features are randomly split 50/50. fit1 is trained on the first half,
  fit2 on the second half; both are then re-expanded to the full feature
  set. Stochastic: set `seed` for reproducibility.

## Usage

``` r
init_two_trajectories(
  X,
  method = c("score", "mincor", "random_split"),
  method1 = c("PCA", "fiedler", "tSNE", "pcurve", "random"),
  modelName = c("homoskedastic", "heteroskedastic"),
  K = NULL,
  rw_q = 2L,
  adaptive = "ml",
  relative_lambda = TRUE,
  num_iter = 5L,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  seed = 42L,
  verbose = FALSE
)
```

## Arguments

- X:

  Numeric matrix (n x d).

- method:

  Initialisation strategy: "score" (default), "mincor", or
  "random_split".

- method1:

  Method for fit1 ordering: "PCA" (default), "fiedler", "tSNE",
  "pcurve", "random".

- modelName:

  Variance structure: "homoskedastic" (default) or "heteroskedastic".

- K:

  Number of pseudotime bins. If NULL (default), auto-selected as
  `max(2, min(50, floor(n/5)))`.

- rw_q:

  Random-walk order for the smoothness prior (default 2).

- adaptive:

  Adaptive hyperparameter mode passed to inner EM: `"ml"` (default) or
  `"none"`. `"prior"` is obsolete and retained only for backward
  compatibility.

- relative_lambda:

  Scale smoothness prior by feature variance (default TRUE).

- num_iter:

  Warm-start EM iterations for each fit (default 5).

- lambda_min, lambda_max:

  Lambda search bounds.

- seed:

  Integer seed (used by "random_split"; ignored otherwise).

- verbose:

  Logical.

## Value

List with `fit1`, `fit2` (both `csmooth_em`), and `seed2` (the column
index used to seed fit2, or NA for "random_split").
