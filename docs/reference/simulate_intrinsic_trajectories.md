# Simulate a multi-ordering intrinsic-trajectory dataset

Simulates a feature-partition dataset with an arbitrary number of
intrinsic orderings. Each ordering contributes its own signal block, and
optional noise features are pure Gaussian noise.

## Usage

``` r
simulate_intrinsic_trajectories(
  n = 200,
  d_signal = c(30, 30),
  d_noise = 20,
  sigma = 0.3,
  seed = 42L,
  trajectory_family = "sinusoidal",
  signal_range = c(0.5, 2),
  linear_slope_range = c(0.8, 2),
  monotone_power_range = c(0.7, 2.5),
  quadratic_curvature_range = c(0.8, 2),
  quadratic_center_range = c(0.25, 0.75),
  intercept_sd = 0,
  sinusoid_freq = 1:4,
  latent_positions = NULL
)
```

## Arguments

- n:

  Integer number of samples.

- d_signal:

  Integer vector giving the number of signal features assigned to each
  intrinsic ordering.

- d_noise:

  Integer number of pure-noise features.

- sigma:

  Observation noise standard deviation added independently to each
  simulated feature trajectory.

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

- latent_positions:

  Optional `n x M` matrix of true latent sample positions. If `NULL`,
  each intrinsic ordering is generated independently from
  `Uniform(0, 1)`.

## Value

A list with components:

- `X`: shuffled simulated data matrix

- `true_assign`: feature labels such as `"A"`, `"B"`, ..., or `"noise"`

- `latent_positions`: `n x M` matrix of true latent sample locations

- `ordering_labels`: ordering labels used in `true_assign`

- `original_order`: column permutation applied to the feature blocks

- `trajectory_family`: the realized family labels for each ordering
