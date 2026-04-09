# Add Known-Noise `S` Support to CAVI and Public Wrappers

## Summary
Add a new public argument `S = NULL` to the CAVI path and all CAVI-facing public wrappers.

Default behavior stays exactly the same when `S = NULL`.

When `S` is supplied, interpret it as known observation standard deviation:
- `length(S) == d`: feature-level known SD, shared across observations for each feature
- `dim(S) == c(n, d)`: observation-level known SD for each cell-feature pair

This support applies to:
- single-ordering `cavi()`
- `do_cavi()`
- `fit_mpcurve()`
- `do_mpcurve()`
- partition-CAVI internals as well, including weighted updates and `smooth_fit` similarity initialization

## Key Changes
### Public API and object behavior
- Add `S = NULL` to `cavi()`, `do_cavi()`, `fit_mpcurve()`, and `do_mpcurve()`.
- Treat `S` as standard deviation, not variance; internally convert once to variance as needed.
- Validate `S` strictly:
  - `NULL`, length-`d` numeric vector, or `n x d` numeric matrix only
  - finite and nonnegative
  - no recycling beyond the explicit length-`d` vector case
- Preserve `S` in the fit object as a separate field, not in `params$sigma2`.
- Add a control marker such as `control$noise_model` with values like:
  - `"estimated_feature_variance"` when `S = NULL`
  - `"known_feature_sd"` when `S` is a length-`d` vector
  - `"known_observation_sd"` when `S` is an `n x d` matrix
- Keep `params$sigma2` unchanged only for the current estimated-noise path.
- For known-noise fits, set `params$sigma2` to `NULL` and update print/summary/as_mpcurve code to report the known-noise mode without assuming a feature-level `sigma2` vector exists.

### Core CAVI updates
- Introduce an internal normalization helper that turns user input `S` into a standard internal representation, ideally both:
  - `sd_mat` as `n x d`
  - `var_mat` as `n x d`
  Even for length-`d` input, expand to matrix form once so downstream math is uniform.
- Split CAVI into two noise paths:
  - existing path when `S = NULL`, with no behavioral changes
  - known-noise path when `S` is supplied
- In the known-noise path:
  - remove the `sigma2` update step entirely
  - keep `lambda_j` variational updates unchanged in spirit
  - replace all likelihood/posterior calculations that assume feature-level `sigma2[j]` with observation-aware weights `1 / var_mat[i, j]`
- Update the posterior `q(U_j)` calculation so each feature uses:
  - per-component weighted counts from `colSums(R / var_mat[, j])`
  - per-component weighted responses from `crossprod(R / var_mat[, j], X[, j])`
- Update responsibility, ELBO, and plug-in log-likelihood computations to use the known variances row-by-row and feature-by-feature.
- Keep ELBO as the monotonic objective in the known-noise path as well.

### Continuation and wrappers
- `do_cavi()` should reuse the stored `S` by default; allow an explicit `S` argument only if it exactly replaces the stored noise specification.
- `do_mpcurve()` should forward `S` into `do_cavi()` for single-ordering fits and reuse stored `S` by default.
- `fit_mpcurve()` should forward `S` through:
  - single-ordering `cavi()` calls
  - partition initialization calls
  - partition weighted-CAVI updates
- `as_mpcurve.cavi()` should preserve the known-noise metadata/field on the wrapped `mpcurve` object.
- For partition fits, all orderings must use the same supplied `S`; this is global measurement noise, not ordering-specific noise.

### Partition-CAVI support
- Thread `S` through the internal builders and weighted update functions in `R/10_partition_cavi.R`.
- Add known-noise variants or generalized implementations for:
  - weighted posterior update
  - weighted responsibility update
  - weighted ELBO/objective scoring
  - single-feature local fits used for feature scoring
- Update `smooth_fit` similarity initialization to use the supplied known-noise model instead of estimating a fresh residual variance when `S` is present.
- Keep `spearman` and `pearson` similarity behavior unchanged.

## Test Plan
- Add CAVI tests for `S` as a length-`d` vector:
  - fit runs without error
  - ELBO is nondecreasing within tolerance
  - no sigma-update path is used
  - returned object records known-noise mode and stored `S`
- Add CAVI tests for `S` as an `n x d` matrix:
  - fit runs without error
  - ELBO is nondecreasing within tolerance
  - responsibilities remain row-stochastic
  - continuation via `do_cavi()` preserves and reuses `S`
- Add regression tests that `S = NULL` reproduces current behavior.
- Add equivalence tests:
  - if `S` is a constant-by-row expansion of a length-`d` vector, matrix and vector inputs produce the same fit up to tolerance
- Add wrapper tests:
  - `fit_mpcurve(..., S=...)` works for intrinsic_dim `= 1`
  - `do_mpcurve()` continues a known-noise fit correctly
- Add partition tests:
  - `soft_partition_cavi` / `fit_mpcurve(intrinsic_dim >= 2, S=...)` run without error
  - weighted objective remains nondecreasing where currently claimed
  - `partition_init = "similarity", similarity_metric = "smooth_fit", S=...` uses the known-noise path successfully

## Assumptions and defaults
- `S` is always standard deviation, never variance.
- No residual extra noise term is added; the observation noise is exactly the supplied `S`.
- Only `NULL`, length-`d`, and `n x d` inputs are supported for `S`.
- Known-noise fits store the noise specification in a dedicated field plus `control$noise_model`; they do not overload `params$sigma2`.
- Existing behavior, output, and tests for `S = NULL` should remain unchanged except for harmless new metadata fields.
