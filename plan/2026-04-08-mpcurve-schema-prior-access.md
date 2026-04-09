# Plan: Unify `mpcurve` Schema and Prior Access

## Summary
- Make `mpcurve` a stable user-facing wrapper across intrinsic dimensions.
- Keep `intrinsic_dim` as the requested/model dimension and add explicit active/displayed dimension fields.
- Expose fitted priors through a unified top-level `priors` block and a public `fitted_prior()` accessor.
- Move effective partition prior state out of `assignment_posterior`.

## Implementation
- Refactor `as_mpcurve.*` so all wrapped fits expose top-level `data`, `measurement_sd`, `priors`, dimension metadata, and consistent high-level fields.
- For partition fits, keep `fits` as the displayed-ordering view while storing requested/active/displayed dimensions separately.
- Update `print.mpcurve`, `summary.mpcurve`, and `plot.mpcurve` to branch on partition-vs-single fit state instead of `intrinsic_dim >= 2`.
- Add `fitted_prior()` methods for `cavi`, `soft_partition_cavi`, and `mpcurve`.
- Store `priors` on raw `cavi` and raw `soft_partition_cavi` objects.
- Restrict `assignment_posterior` to assignment-side / legacy-Dirichlet metadata only.

## Tests
- Cover single-ordering `priors` and `fitted_prior()` access.
- Cover compacted partition views with requested/active/displayed dimensions.
- Reproduce and prevent the compacted-partition `plot(..., dims = 1)` regression.
- Check that partition effective prior is available through `priors` / `fitted_prior()` and removed from `assignment_posterior`.

## Assumptions
- `requested_intrinsic_dim` remains as a compatibility alias for now.
- `fitted_prior(type = "partition")` returns the full partition-prior record, while specifying `ordering` returns a single effective `omega_m`.
- Fixed partition priors are reported after active-set renormalization in `priors$partition$omega`.
