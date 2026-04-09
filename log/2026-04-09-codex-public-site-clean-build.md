Agent: codex
Update title: public-site-clean-build
Update date: 2026-04-09

Key updates:
- Added `scripts/build_public_site.R` to build only the intended public pkgdown homepage, reference pages, and articles.
- Removed legacy and internal materials from the generated website output, including archived `Intro`, internal reference pages, and root markdown pages such as `AGENTS.md` / `CLAUDE.md`.
- Simplified public website copy so the homepage no longer links to internal derivation documents.
- Kept the public reference surface aligned with the curated `mpcurve`-first API.
- Filtered generated `docs/pkgdown.yml` article metadata so the committed site index, search metadata, and sitemap expose only the public article whitelist.

Key files touched:
- `scripts/build_public_site.R`
- `_pkgdown.yml`
- `README.Rmd`
- `DESCRIPTION`
- `vignettes/partition.Rmd`

Verification:
- `Rscript scripts/build_public_site.R`
- `devtools::test(filter = "public-api|mpcurve-cavi", stop_on_failure = TRUE)`
- `R CMD build . --no-build-vignettes --no-manual`

Notes:
- The build script intentionally avoids `pkgdown::build_site()` because that default path picks up top-level markdown files that are not part of the public website.
