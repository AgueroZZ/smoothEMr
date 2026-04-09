Agent: codex
Update title: public-api-documentation-tidy-up
Update date: 2026-04-09

Key updates:
- Shrunk the exported namespace to the curated `mpcurve`-first public API and removed legacy/internal topics from generated `man/` pages.
- Tightened `fit_mpcurve()` to a CAVI-only public wrapper, kept the `algorithm` formal for compatibility, and restricted `do_mpcurve()` to CAVI-backed `mpcurve` fits.
- Simplified `print.mpcurve()` so summary-only diagnostics stay in `summary.mpcurve()` / `print.summary.mpcurve()`.
- Reworked `README`, `_pkgdown.yml`, and the main public vignettes around `fit_mpcurve()`, `do_mpcurve()`, and `intrinsic_dim`.
- Rebuilt `README`, `NAMESPACE`, `man/`, and `docs/` using an allowlisted pkgdown reference build so stale legacy/internal reference pages were removed from `docs/reference/`.

Key files touched:
- `R/00_utils.R`
- `R/07_mpcurve.R`
- `README.Rmd`
- `_pkgdown.yml`
- `DESCRIPTION`
- `tests/testthat/test-public-api.R`

Verification:
- `devtools::test(filter = "public-api|mpcurve-cavi", stop_on_failure = TRUE)`
- `R CMD build . --no-build-vignettes --no-manual`

Notes:
- The public site keeps `Intro` only as an archived article section without a navbar entry.
- `pkgdown` reported missing alt text in `vignettes/fitness.rmd`; this did not block the site build, but the article still has accessibility follow-up work.
