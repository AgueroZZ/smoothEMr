# Greedy backward filtering of features for csmoothEM ordering inference

Implements a greedy "backward" feature filtering strategy to mitigate
the case where most features are noise and do not follow any latent
ordering. The algorithm:

1.  Fits csmoothEM on all features to obtain an initial ordering (via
    responsibilities \\\Gamma\\).

2.  Computes per-feature collapsed contributions \\C_j\\ (via
    `compute_C_by_coord_csmooth`).

3.  Removes a batch of features with the smallest \\C_j\\.

4.  Refits csmoothEM on the remaining features for a few iterations.

5.  Repeats until a stopping rule is met.

This procedure is intended as a preprocessing step before more ambitious
tasks such as feature partitioning across multiple orderings.

## Usage

``` r
greedy_backward_filter_csmooth(
  X,
  method = c("fiedler", "PCA", "tSNE", "pcurve", "random"),
  K = 50,
  modelName = c("homoskedastic", "heteroskedastic"),
  adaptive = "ml",
  num_iter_init = 10,
  num_iter_refit = 5,
  discretization = c("equal", "quantile", "kmeans"),
  batch = 20,
  min_keep = 20,
  tau = NULL,
  max_rounds = 50,
  verbose = TRUE,
  ...
)
```

## Arguments

- X:

  Numeric matrix `(n x d)`.

- method:

  Ordering method passed to
  [`initialize_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/initialize_csmoothEM.md).
  One of `"fiedler"`, `"PCA"`, `"tSNE"`, `"pcurve"`, `"random"`.

- K:

  Integer \\\ge 2\\. Number of mixture components.

- modelName:

  Either `"homoskedastic"` or `"heteroskedastic"`.

- adaptive:

  Adaptive mode passed to
  [`do_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM.md)
  when refitting. `"ml"` is the recommended default. `"prior"` is
  obsolete and retained only for backward compatibility.

- num_iter_init:

  Integer \\\ge 1\\. Number of warm-start iterations for the initial
  fit.

- num_iter_refit:

  Integer \\\ge 1\\. Number of iterations for each refit after feature
  removal.

- discretization:

  Discretization method for initialization passed to
  `initialize_csmoothEM`. Recommended: `"quantile"` to avoid empty
  components.

- batch:

  Integer \\\ge 1\\. Number of lowest-scoring features (smallest
  \\C_j\\) removed per round.

- min_keep:

  Integer \\\ge 1\\. Minimum number of features to keep; stops if fewer
  would remain.

- tau:

  Optional numeric threshold. If provided, stops when `min(Cj) >= tau`.

- max_rounds:

  Integer \\\ge 1\\. Maximum number of greedy rounds.

- verbose:

  Logical; print a one-line summary each round.

- ...:

  Additional arguments passed to `initialize_csmoothEM` (e.g. ordering
  controls).

## Value

A list with components:

- `keep_cols`: integer indices of retained features (w.r.t. the original
  X).

- `drop_cols`: integer indices of removed features (w.r.t. the original
  X).

- `fit`: final fitted `csmooth_em` object on the retained features.

- `history`: data.frame with per-round diagnostics (`C_total`, `min_Cj`,
  etc.).
