Agent: codex
Update title: cavi-error-diagnostics
Update date: 2026-04-08

Key updates:
- Added detailed CAVI diagnostics when the posterior precision update fails a Cholesky factorization, including iteration, feature index, effective component count, lambda scale, and a stabilization hint.
- Updated `fit_mpcurve()` to surface the underlying `cavi()` error message instead of replacing it with a generic "failed to produce a valid cavi fit" stop.
- Added a regression test covering the known-measurement-SD failure mode with collapsed component mass.

Key files touched:
- `R/09_cavi.R`
- `R/07_mpcurve.R`
- `tests/testthat/test-cavi.R`

Follow-up notes:
- In the extreme known-noise eLife case, the failure is caused by a numerically non-positive-definite posterior precision matrix after responsibilities collapse onto very few components under an intrinsic RW prior with `ridge = 0`.
- A small positive ridge such as `1e-6` stabilizes that edge case.
