# Known-Noise CAVI Support

- Agent: `codex`
- Update title: `known-noise-cavi-support`
- Update date: `2026-04-05`

## Key updates

- Added public `S = NULL` support to `cavi()`, `do_cavi()`, `fit_mpcurve()`, and `do_mpcurve()`.
- Defined `S` as known measurement standard deviation:
  - length-`d` vector: one known SD per feature, shared across observations
  - `n x d` matrix: one known SD per observation-feature pair
- Preserved the old behavior exactly when `S = NULL`.
- Implemented a known-noise CAVI path that uses the supplied measurement SD directly and skips `sigma2` inference.
- Added `measurement_sd` and `control$noise_model` to CAVI-facing objects so downstream code can tell whether the fit uses inferred feature variance or known measurement noise.
- Updated partition-CAVI weighted updates and `smooth_fit` similarity initialization so they also respect the supplied `S`.
- Extended `as_mpcurve()`, summary, and print methods to handle known-noise fits cleanly.

## Key files touched

- `R/09_cavi.R`
- `R/10_partition_cavi.R`
- `R/07_mpcurve.R`
- `tests/testthat/test-cavi.R`
- `tests/testthat/test-mpcurve-cavi.R`
- `tests/testthat/test-partition-cavi.R`

## Validation

- Ran the full `tests/testthat` suite after implementation.
- Result: `PASS 432`, `FAIL 0`, `WARN 0`, `SKIP 0`.

## Notes for future agents

- For known-noise fits, `params$sigma2` is intentionally `NULL`; use `measurement_sd` plus `control$noise_model` instead.
- Partition fits should treat `S` as a global observation model shared across all orderings, not as an ordering-specific parameter.
