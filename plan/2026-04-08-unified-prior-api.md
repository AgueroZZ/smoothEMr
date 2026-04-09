# Unified Prior API for CAVI and Partition-CAVI

## Summary
Replace the prior-facing public API with:

- `position_prior = c("adaptive", "fixed")`
- `position_prior_init = NULL`
- `partition_prior = c("adaptive", "fixed")`
- `partition_prior_init = NULL`

Keep backward compatibility through deprecated aliases:

- `pi_init` -> `position_prior_init`
- `assignment_prior` / `ordering_alpha` -> deprecated partition controls

## Decisions
- `position_prior = "adaptive"` keeps the current empirical-Bayes `pi` update with a strict positivity floor.
- `position_prior = "fixed"` holds `pi` fixed at `position_prior_init`, defaulting to uniform.
- `partition_prior = "adaptive"` uses empirical-Bayes `omega <- normalize(colSums(weights))`.
- Adaptive partition priors may hit exact-zero active ordering mass; those orderings are dropped from the active set before the next prior/objective/softmax evaluation.
- `partition_prior = "fixed"` holds a fixed `omega`, defaulting to uniform over the current active orderings.
- For `intrinsic_dim > 1`, `position_prior_init` is a single length-`K` vector broadcast to every ordering. Per-ordering differences still require `fits_init`.
- Continuation inherits stored prior modes and current prior values; no continuation override API is introduced.

## Implementation Areas
- Update `R/09_cavi.R` to use `position_prior` / `position_prior_init`.
- Update `R/10_partition_cavi.R` to add the new assignment-block abstraction and deprecated legacy mapping.
- Update `R/07_mpcurve.R` to expose the new public arguments and preserve compatibility warnings.
- Add test coverage for fixed/adaptive position priors, fixed/adaptive partition priors, active-set dropping, and deprecated alias routing.
- Add a short `log/` entry after implementation.
