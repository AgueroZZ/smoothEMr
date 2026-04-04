# Simulate a dual-trajectory dataset with configurable trajectory families

Backward-compatible convenience wrapper around
[`simulate_intrinsic_trajectories`](https://aguerozz.github.io/MPCurver/reference/simulate_intrinsic_trajectories.md)
for the two-ordering case. Preserves the old `t1`/`t2` return values and
the `crossing` shortcut.

## Usage

``` r
simulate_dual_trajectory(
  n = 200,
  d1 = 30,
  d2 = 30,
  d_noise = 20,
  sigma = 0.3,
  crossing = FALSE,
  seed = 42L,
  trajectory_family = c("sinusoidal", "sinusoidal"),
  signal_range = c(0.5, 2),
  linear_slope_range = c(0.8, 2),
  monotone_power_range = c(0.7, 2.5),
  quadratic_curvature_range = c(0.8, 2),
  quadratic_center_range = c(0.25, 0.75),
  intercept_sd = 0,
  sinusoid_freq = 1:4
)
```

## Arguments

- n:

  Integer number of samples.

- d1, d2:

  Integer numbers of signal features assigned to trajectories A and B.

- d_noise:

  Integer number of pure-noise features.

- sigma:

  Observation noise standard deviation added independently to each
  simulated feature trajectory.

- crossing:

  Logical; if `TRUE`, the second latent ordering is generated as a noisy
  reversed version of the first, making the two orderings harder to
  disentangle.

- seed:

  Optional integer seed.

- trajectory_family:

  Character scalar or length-`M` vector specifying the feature family
  for each ordering block. Each entry must be one of `"sinusoidal"`,
  `"linear"`, `"monotone"`, or `"quadratic"`.

- signal_range:

  Length-2 positive range controlling the per-feature signal amplitude.

- linear_slope_range:

  Length-2 positive range for linear slopes.

- monotone_power_range:

  Length-2 positive range controlling the power-like basis shapes used
  in the monotone family. The monotone generator mixes several
  increasing bases internally to create more varied monotone curves.

- quadratic_curvature_range:

  Length-2 positive range for quadratic curvature magnitudes.

- quadratic_center_range:

  Length-2 range inside `[0, 1]` for the quadratic vertex.

- intercept_sd:

  Standard deviation of feature-specific intercept shifts. Defaults to
  `0` to preserve the historical simulation baseline.

- sinusoid_freq:

  Integer vector of admissible sinusoid frequencies.

## Value

A list with the same fields as
[`simulate_intrinsic_trajectories`](https://aguerozz.github.io/MPCurver/reference/simulate_intrinsic_trajectories.md),
plus backward-compatible aliases `t1` and `t2`.
