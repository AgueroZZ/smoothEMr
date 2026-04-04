# Simulate a two-ordering GP dataset (Matern) for feature partitioning

Simulates an \\N \times D\\ dataset with two latent sample orderings
\\t_1\\ and \\t_2\\. The first \\D/2\\ features are generated from a
Gaussian process over \\t_1\\, and the remaining \\D/2\\ features are
generated from an independent Gaussian process over \\t_2\\.

Each feature is a GP draw evaluated at the \\N\\ sample locations.
Optionally:

- permute feature columns (`permute_cols`),

- permute the sample order within block 2 (`permute_rows_block2`),

- shift the data to be positive (`shift_positive`),

- add i.i.d. Gaussian noise (`noise_sd`).

Column names are assigned \*before\* permutation and then permuted
consistently with the columns, so that names remain aligned with the
returned `true_group`.

## Usage

``` r
simulate_two_order_gp_dataset(
  N = 1000,
  D = 16,
  t_range = c(0, 10),
  range = 5,
  smoothness = 2.5,
  variance = 3,
  noise_sd = 0.05,
  shift_positive = TRUE,
  permute_cols = FALSE,
  permute_rows_block2 = TRUE,
  seed = NULL
)
```

## Arguments

- N:

  Integer \\\ge 2\\. Number of samples (rows).

- D:

  Integer \\\ge 2\\ and even. Number of features (columns).

- t_range:

  Numeric length-2 vector. Range for sampling `t1` and `t2`.

- range, smoothness, variance:

  Matern GP hyperparameters.

- noise_sd:

  Nonnegative numeric. Standard deviation of i.i.d. Gaussian noise added
  to `X`.

- shift_positive:

  Logical; if TRUE, shift each GP block so its minimum is 1.

- permute_cols:

  Logical; if TRUE, permute feature columns.

- permute_rows_block2:

  Logical; if TRUE, permute rows of the second GP block before
  combining.

- seed:

  Optional integer seed. If not NULL, sets `set.seed(seed)`.

## Value

A list with components:

- `X`: numeric matrix \\N \times D\\.

- `t1`, `t2`: numeric vectors length N (latent orderings).

- `permut_cols`: integer vector length D (the applied column
  permutation; identity if `permute_cols=FALSE`).

- `true_group`: integer vector length D in `1,2` indicating which
  ordering generated each feature (after permutation).

- `row_perm_block2`: integer vector length N giving the row permutation
  applied to block 2 (or `NULL`).

## Details

Requires [`fields::rdist`](https://rdrr.io/pkg/fields/man/rdist.html),
[`fields::Matern`](https://rdrr.io/pkg/fields/man/Exponential.html), and
[`MASS::mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html).
