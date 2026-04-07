# Fix: Known-Noise Partition CAVI Initialization

- Agent: `claude-sonnet-4-6`
- Update title: `fix-known-noise-partition-init`
- Update date: `2026-04-06`

## Problem

Partition CAVI with known measurement SD (`S` non-NULL) gave ~50% partition
accuracy (random chance) on dual-trajectory data, while estimated-σ² mode
gave 100%. Since known S carries strictly more information than estimated σ²,
this was a bug.

## Root Cause

In `.cavi_init_m_trajectories_similarity()`, after similarity-based clustering
correctly assigned features to two groups, each group's subset CAVI fit was
expanded to full-d by calling:

```r
cavi(X = X, S = S, max_iter = num_iter, responsibilities_init = subset_R)
```

with `num_iter = 5`. With known S, **all d features contribute equally** to
the R (responsibilities) update in every CAVI iteration. After 5 iterations
on all features without feature weights, both orderings converged to the same
compromise R, destroying the subset-specialized cell orderings before the
partition loop even started.

With estimated σ²: mismatched features self-downweighted via 1/σ²_j in the R
update ("σ² specialization feedback loop"), so 5 warm-up iterations were
harmless.

## Fix

Set `max_iter = 0L` unconditionally in the full-d expansion `cavi()` call
(lines 1014–1030 of `R/10_partition_cavi.R`). With zero iterations:

1. The subset-specialized R is preserved.
2. `q_u`, `lambda`, and (for estimated-σ²) `sigma2` are initialized once for
   all d features from the preserved R.
3. The partition CAVI loop (with feature weights `w_jm`) then handles all
   refinement — it already does this correctly for both known and estimated
   noise models.

## Files Changed

- `R/10_partition_cavi.R`: `max_iter = 0L` in `.cavi_init_m_trajectories_similarity()`
- `tests/testthat/test-partition-cavi.R`: restrict the objective_history
  monotonicity check to the convergence phase only (post-annealing); during
  the annealing phase the objective function itself changes with T, so
  monotonicity is not guaranteed by design
- `important_derivations/known_noise_cavi_derivation.Rmd`: replaced Section
  4.4 (which described a reverted estimated-σ² band-aid) with the correct
  `max_iter = 0` explanation

## Validation

All 432 tests pass (`FAIL 0 | WARN 0 | SKIP 0 | PASS 432`).

Partition accuracy with known S (K=30, seeds 1–10):

| Seed | Accuracy |
|------|----------|
| 1–10 | 100%     |

Also validated:
- K=50, seeds 1–5: 100%
- Heteroscedastic S (n×d matrix), seeds 1–5: 100%

## Notes for Future Agents

- The similarity-based clustering (`partition_init = "similarity"`) correctly
  separates features for independent orderings via |Spearman| correlation.
  The bug was entirely in the subsequent warm-up step, not in the clustering.
- The partition CAVI weighted updates (`.do_cavi_weighted()`) are
  mathematically exact: every micro-step (R → π → λ → q_u) monotonically
  increases the ELBO. Confirmed by the micro-step ELBO diagnostic in
  `scripts/elbo_microstep_diagnostic.R`.
- Orderings are always assumed independent (never "crossing"). Crossing
  trajectories are by definition the same trajectory in different directions
  and need not be handled specially.
- `objective_history` covers both the annealing phase (T changes, objective
  not monotone) and the convergence phase (T fixed, objective monotone). Tests
  should only assert monotonicity over the convergence phase
  (`seq(from = res$n_anneal, to = length(res$objective_history))`).
