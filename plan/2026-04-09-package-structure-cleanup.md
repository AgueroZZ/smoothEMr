Agent: codex
Plan title: package-structure-cleanup
Plan date: 2026-04-09

Scope:
- Keep package-helper directories out of the source tarball.
- Preserve versioned coordination files in `log/` and `plan/`.
- Remove the current Rd build warning in `fit_mpcurve`.

Planned steps:
1. Audit the current `R CMD build` output and identify leaked directories/files.
2. Update `.Rbuildignore` so package builds exclude helper artifacts and generated vignette outputs.
3. Adjust `.gitignore` so future `log/` and `plan/` files can be tracked normally.
4. Fix the `fit_mpcurve` documentation text that currently emits an unknown Rd macro warning.
5. Rebuild the package and verify the tarball contents.
