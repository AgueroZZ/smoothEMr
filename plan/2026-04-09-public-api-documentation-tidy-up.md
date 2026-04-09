Agent: codex
Plan title: public-api-documentation-tidy-up
Plan date: 2026-04-09

Scope:
- Shrink the package to a curated `mpcurve`-first public API.
- Make the public wrapper CAVI-only while retaining legacy implementations internally.
- Rework documentation, tests, and pkgdown reference pages around the curated public surface.

Planned steps:
1. Remove exports and public Rd topics for legacy backends, converters, and internal helpers.
2. Tighten `fit_mpcurve()` / `do_mpcurve()` semantics to the CAVI-only public workflow.
3. Simplify `print.mpcurve()` while keeping detailed diagnostics in `summary.mpcurve()`.
4. Update `README`, public vignettes, and `_pkgdown.yml` to document only the curated public interface.
5. Regenerate `NAMESPACE`, `man/`, and `docs/`, remove stale reference pages, and verify the exported API with tests.
