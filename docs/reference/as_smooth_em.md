# Construct a smooth_em object from an EM_algorithm fit

Construct a smooth_em object from an EM_algorithm fit

## Usage

``` r
as_smooth_em(
  fit,
  Q_prior = NULL,
  Q_base = NULL,
  lambda = NULL,
  q = NULL,
  ridge = NULL,
  modelName = NULL,
  relative_lambda = FALSE,
  rank_deficiency = 0,
  eigen_tol = NULL,
  nugget = 0,
  max_inner = 10,
  inner_tol = 1e-06,
  meta = NULL
)
```

## Arguments

- fit:

  Output list from
  [`EM_algorithm()`](https://aguerozz.github.io/MPCurver/reference/EM_algorithm.md).

- Q_prior:

  Optional precision matrix used in fitting (legacy; discouraged for
  continuing).

- Q_base:

  Optional \*base\* precision matrix (without lambda). Recommended.

- lambda:

  Optional penalty strength; used with `Q_base`.

- q:

  Optional RW order, if relevant.

- ridge:

  Optional ridge used in building `Q_base`.

- modelName:

  Covariance model used in M-step (e.g. "EEI").

- relative_lambda:

  Logical; whether relative-lambda scaling is used.

- rank_deficiency:

  Rank deficiency used in generalized logdet / EEI update.

- eigen_tol:

  Optional tolerance for generalized logdet.

- nugget:

  Diagonal jitter used in covariance updates.

- max_inner, inner_tol:

  M-step inner loop controls.

- meta:

  Optional list of extra metadata to store.

## Value

An object of class `smooth_em`.
