# Simulate a toy dataset from the current single-ordering CAVI model

Simulate a toy dataset from the current single-ordering CAVI model

## Usage

``` r
simulate_cavi_toy(
  n = 150,
  d = 20,
  K = 8,
  rw_q = 2L,
  lambda_range = c(0.5, 3),
  sigma_range = c(0.08, 0.18),
  pi = NULL,
  ridge = 0.001,
  seed = NULL
)
```

## Arguments

- n:

  Number of observations.

- d:

  Number of features.

- K:

  Number of latent positions/components.

- rw_q:

  Random-walk order for the GMRF prior.

- lambda_range:

  Length-2 positive range for feature-specific `lambda_j`.

- sigma_range:

  Length-2 positive range for feature-specific `sigma_j`.

- pi:

  Optional component probabilities; defaults to uniform.

- ridge:

  Small ridge used only for simulation to make the prior proper.

- seed:

  Optional integer seed.

## Value

A list with `X`, `z`, `pi`, `mu`, `sigma2`, `lambda_vec`, and `Q_K`.
