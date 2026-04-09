Agent: codex
Update title: greedy-similarity-provenance
Update date: 2026-04-07

Key updates:
- Preserved full similarity initialization metadata across greedy compaction/refit paths, including backward greedy runs that temporarily switch to `fits_init`.
- Normalized greedy model comparison to the active intrinsic dimension after compaction, so collapsed candidates no longer count as evidence for a larger model.
- Added greedy comparison diagnostics for requested vs active dimensions and collapse status.
- Preserved greedy provenance through `do_mpcurve()` for both partition fits and greedy-selected single-ordering fits.
- Added regression tests for backward similarity provenance retention, collapsed-candidate rejection, and `do_mpcurve()` provenance preservation on greedy-selected 1D fits.

Key files touched:
- `R/07_mpcurve.R`
- `R/10_partition_cavi.R`
- `tests/testthat/test-mpcurve-cavi.R`
- `man/fit_mpcurve.Rd`
- `man/soft_partition_cavi.Rd`

Tests run:
- `Rscript -e 'suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE)); testthat::test_file("tests/testthat/test-mpcurve-cavi.R")'`
- `Rscript -e 'suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE)); testthat::test_file("tests/testthat/test-partition-m-cavi.R")'`
