Agent: codex
Update title: unified-prior-api
Update date: 2026-04-08

Key updates:
- Added the unified prior API across the public CAVI wrappers:
  `position_prior`, `position_prior_init`, `partition_prior`, and
  `partition_prior_init`.
- Updated single-ordering `cavi()` and `do_cavi()` so `position_prior = "fixed"`
  keeps `pi` constant, while `position_prior = "adaptive"` keeps the existing
  empirical-Bayes update with a positivity floor.
- Updated partition CAVI to support fixed and adaptive ordering priors:
  adaptive mode now uses empirical-Bayes `omega_hat ∝ colSums(weights)` over
  the active ordering set, and exact zero-mass active orderings are dropped
  before the next objective/prior evaluation.
- Kept the old `assignment_prior` / `ordering_alpha` interface as a deprecated
  compatibility path. `assignment_prior = "uniform"` maps to the new fixed
  partition prior; `assignment_prior = "dirichlet"` remains available only via
  the deprecated legacy path.
- Added fixed-prior propagation through `fit_mpcurve()`, `soft_partition_cavi()`,
  `soft_two_trajectory_cavi()`, and `do_mpcurve()`, including continuation on
  existing fits.
- Extended tests to cover fixed `position_prior`, fixed `partition_prior`, the
  new public formals, and the updated adaptive-partition objective semantics.
- Regenerated `man/` for the updated public API.

Key files touched:
- `R/07_mpcurve.R`
- `R/09_cavi.R`
- `R/10_partition_cavi.R`
- `tests/testthat/test-cavi.R`
- `tests/testthat/test-mpcurve-cavi.R`
- `tests/testthat/test-partition-m-cavi.R`
- `man/cavi.Rd`
- `man/do_mpcurve.Rd`
- `man/fit_mpcurve.Rd`
- `man/soft_partition_cavi.Rd`
- `man/soft_two_trajectory_cavi.Rd`

Tests run:
- `Rscript -e "testthat::test_local(filter='mpcurve-cavi')"`
- `Rscript -e "testthat::test_local(filter='partition-m-cavi')"`
- `Rscript -e "testthat::test_local(filter='cavi')"`

Follow-up notes:
- The deprecated `assignment_prior = "dirichlet"` path is still exercised by
  tests but now emits deprecation warnings by design.
- The package still accepts deprecated `pi_init` for `cavi()`, but new code and
  examples should use `position_prior_init`.
