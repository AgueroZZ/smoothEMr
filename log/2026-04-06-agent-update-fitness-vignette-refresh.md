Agent: codex
Update title: fitness-vignette-refresh
Update date: 2026-04-06

Key updates:
- Rewrote `vignettes/fitness.rmd` around a clearer end-to-end user workflow built on `fit_mpcurve()`.
- Reorganized the vignette into four questions: single ordering, smoothing, multiple-ordering search, and known-measurement-SD refinement.
- Clarified that `S` in the public API represents known measurement standard deviations and removed the earlier inconsistent `sqrt(S)` usage.
- Tightened the prose so the vignette reads as a guided analysis rather than a rough analysis notebook.

Key files touched:
- `vignettes/fitness.rmd`

Checks:
- Extracted and parsed all vignette R chunks with `knitr::purl()` and `parse()`.

Follow-up notes:
- Full rendering of this vignette is computationally heavier than a syntax check because it fits several real-data models.
