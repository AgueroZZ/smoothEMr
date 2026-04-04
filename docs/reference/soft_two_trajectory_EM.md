# Soft dual-trajectory partition via annealed CAVI

Jointly infers two latent 1-D orderings and a soft feature-level
partition via coordinate-ascent variational inference (CAVI). The
algorithm alternates between:

Outer E-step: score each feature under current Gamma_1, Gamma_2 via
[`compute_C_by_coord_csmooth()`](https://aguerozz.github.io/MPCurver/reference/compute_C_by_coord_csmooth.md);
compute `Softmax(S/T)` to get soft feature weights `w_jm`.

Inner EM: refit each trajectory via
[`do_csmoothEM_weighted()`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_weighted.md).
The E-step uses weighted log-likelihood contributions, and when
`adaptive = "ml"` the collapsed `lambda`/`sigma2` hyper-step is
reweighted as well. The final MAP refresh of `mu` is unchanged.

Phase 1 anneals the temperature T from `T_start` to `T_end` to escape
local optima. Each outer step then performs `inner_iter` weighted
trajectory updates before recomputing feature weights. Phase 2 runs the
same untempered updates at `T = T_end` until the recorded objective
stabilises.

This is the legacy `csmooth_em`-based soft partition routine. It is
retained for benchmark comparison and internal regression testing. For
new analyses, prefer
[`soft_two_trajectory_cavi`](https://aguerozz.github.io/MPCurver/reference/soft_two_trajectory_cavi.md),
which uses the explicit variational `cavi` model throughout and is the
recommended dual-ordering soft-partition backend.

## Usage

``` r
soft_two_trajectory_EM(
  X,
  fit1_init = NULL,
  fit2_init = NULL,
  init_method = c("score", "mincor", "random_split"),
  init_method1 = c("PCA", "fiedler", "tSNE", "pcurve", "random"),
  init_seed = 42L,
  modelName = c("homoskedastic", "heteroskedastic"),
  K = NULL,
  T_start = 5,
  T_end = 1,
  n_outer = 25,
  inner_iter = 1L,
  max_converge_iter = 100L,
  tol_outer = 1e-04,
  score_mode = c("ml", "none"),
  rw_q = 2L,
  relative_lambda = TRUE,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  adaptive = "ml",
  sigma_update = c("ml", "mstep"),
  hard_assign_final = FALSE,
  verbose = TRUE
)
```

## Arguments

- X:

  Numeric matrix (n x d).

- fit1_init:

  Pre-built `csmooth_em` for trajectory A. If NULL, auto-initialised via
  `init_method`.

- fit2_init:

  Pre-built `csmooth_em` for trajectory B. If NULL, auto-initialised via
  `init_method`.

- init_method:

  Initialisation strategy when fits are NULL: "score" (default),
  "mincor", or "random_split".

- init_method1:

  Ordering method for fit1: "PCA" (default), "fiedler", "tSNE",
  "pcurve", or "random".

- init_seed:

  Integer seed for "random_split" (default 42).

- modelName:

  Variance structure: "homoskedastic" (default) or "heteroskedastic".

- K:

  Number of pseudotime bins. If NULL (default), auto-selected as
  `max(2, min(50, floor(n/5)))`. Ignored when `fit1_init`/`fit2_init`
  are supplied.

- T_start:

  Starting annealing temperature (default 5).

- T_end:

  Final temperature; also the Phase 2 temperature (default 1.0).

- n_outer:

  Number of Phase 1 annealing steps (default 25).

- inner_iter:

  Weighted inner EM sweeps per outer step (default 1).

- max_converge_iter:

  Maximum Phase 2 iterations (default 100).

- tol_outer:

  Relative ELBO convergence tolerance for Phase 2 (default 1e-4).
  Checked as \|Delta F\| / (\|F\| + 1).

- score_mode:

  Outer scoring mode: "ml" (collapsed marginal, default) or "none"
  (penalised ELBO only).

- rw_q:

  Random-walk order for the smoothness prior (default 2).

- relative_lambda:

  Scale smoothness prior by feature variance (default TRUE).

- lambda_min, lambda_max:

  Lambda search bounds.

- adaptive:

  Adaptive hyperparameter mode for the inner EM: `"ml"` (collapsed
  marginal, default) or `"none"`. `"prior"` is obsolete and retained
  only for backward compatibility.

- sigma_update:

  Sigma estimation method inside inner EM when `adaptive = "ml"`: `"ml"`
  (collapsed marginal, default) or legacy `"mstep"` (weighted SSE). Use
  `"ml"` to jointly optimize \\\sigma^2\\ and \\\lambda\\; `"mstep"` is
  retained only for backward compatibility.

- hard_assign_final:

  Snap final soft weights to hard 0/1 (default FALSE).

- verbose:

  Print iteration log (default TRUE).

## Value

A list with:

- pi_weights:

  d x 2 matrix of soft feature assignments for trajectories A and B.

- assign:

  Character vector of hard assignments ("A" or "B").

- fit1, fit2:

  Final `csmooth_em` objects for each trajectory.

- ll_history:

  Full ELBO F at each post-M-step iteration.

- score_history, weight_history:

  Per-iteration scores and weights.

- T_schedule:

  Annealing temperature schedule.

- converged:

  Logical: did Phase 2 converge?

- n_anneal:

  Number of Phase 1 iterations.

## See also

[`soft_two_trajectory_cavi`](https://aguerozz.github.io/MPCurver/reference/soft_two_trajectory_cavi.md)
