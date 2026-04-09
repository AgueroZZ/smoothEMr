Agent: codex
Update title: package-structure-cleanup
Update date: 2026-04-09

Key updates:
- Excluded `internal/`, `log/`, `plan/`, `scripts/`, generated vignette outputs, and agent-only metadata from package builds via `.Rbuildignore`.
- Restored normal version control behavior for `log/` and `plan/` by removing them from `.gitignore`.
- Fixed the `fit_mpcurve` documentation text so `R CMD build` no longer emits the unknown `\propto` Rd macro warning.

Key files touched:
- `.gitignore`
- `.Rbuildignore`
- `R/07_mpcurve.R`
- `man/fit_mpcurve.Rd`
- `docs/reference/fit_mpcurve.html`

Verification:
- Rebuilt the package with `R CMD build . --no-build-vignettes --no-manual`.
- Checked the generated tarball contents for leaked helper directories and generated vignette artifacts.

Follow-up notes:
- The repo still contains generated pkgdown output changes under `docs/`; those look like a normal site rebuild rather than structural corruption.
- The top-level `elife-61271-fig2-data1-v2.csv` is still shipped in the source tarball; that is separate from this cleanup and can be revisited later if you want a leaner release bundle.
