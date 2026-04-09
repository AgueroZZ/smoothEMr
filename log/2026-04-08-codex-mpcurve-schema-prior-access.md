Agent: codex
Update title: mpcurve-schema-prior-access
Update date: 2026-04-08

Key updates:
- Unified `mpcurve` wrappers so `intrinsic_dim` now represents the requested/model dimension, with new `active_intrinsic_dim` and `displayed_intrinsic_dim` fields.
- Added top-level `data`, `measurement_sd`, `control`, `converged`, and `priors` fields to wrapped fits, plus aggregated per-ordering state for partition wrappers.
- Added exported `fitted_prior()` methods for `cavi`, `soft_partition_cavi`, and `mpcurve`.
- Moved effective partition prior storage to `priors$partition`; `assignment_posterior` now keeps only assignment-side / legacy-Dirichlet metadata.
- Updated `print.mpcurve`, `summary.mpcurve`, and `plot.mpcurve` to dispatch on actual partition state instead of `intrinsic_dim >= 2`, fixing compacted partition plotting.

Key files touched:
- `R/07_mpcurve.R`
- `R/09_cavi.R`
- `R/10_partition_cavi.R`
- `tests/testthat/test-mpcurve-cavi.R`
- `tests/testthat/test-partition-m-cavi.R`
- `tests/testthat/test-partition-cavi.R`

Tests run:
- `Rscript -e "pkgload::load_all('.', quiet=TRUE); testthat::test_file('tests/testthat/test-mpcurve-cavi.R')"`
- `Rscript -e "pkgload::load_all('.', quiet=TRUE); testthat::test_file('tests/testthat/test-partition-m-cavi.R')"`
- `Rscript -e "pkgload::load_all('.', quiet=TRUE); testthat::test_file('tests/testthat/test-cavi.R')"`
- `Rscript -e "pkgload::load_all('.', quiet=TRUE); testthat::test_file('tests/testthat/test-partition-cavi.R')"`
- `Rscript -e "roxygen2::roxygenise()"`

Follow-up notes:
- The deprecated `assignment_prior = "dirichlet"` compatibility tests still warn, which is expected.
- Existing dirty worktree items outside the files above were left untouched.
