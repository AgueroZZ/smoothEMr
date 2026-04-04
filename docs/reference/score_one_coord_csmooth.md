# Score one coordinate against fixed responsibilities (Gamma) for csmoothEM (internal)

Convenience wrapper used in greedy partitioning. Returns both the scalar
score and a 1D fitted object (`one`) that can be appended to
csmooth-style parameters.

## Usage

``` r
score_one_coord_csmooth(
  X,
  j,
  Gamma,
  Q_K,
  rw_q = 2,
  score_mode = c("ml", "none"),
  relative_lambda = TRUE,
  lambda_min = 1e-10,
  lambda_max = 1e+10
)
```

## Arguments

- X:

  Numeric matrix `(n x d)`.

- j:

  Integer coordinate index in `1:d`.

- Gamma:

  Numeric matrix `(n x K)` of responsibilities.

- Q_K:

  Numeric matrix `(K x K)` base RW precision (lambda=1).

- rw_q:

  Integer \\\ge 0\\. Rank deficiency along K.

- score_mode:

  One of `"ml"` or `"none"`.

- relative_lambda:

  Logical.

- lambda_min, lambda_max:

  Bounds for lambda optimization (used when `score_mode="ml"`).

## Value

A list with:

- `score`: scalar score for coordinate `j`.

- `one`: list with `mu_vec` (length K) and `sigma2` (scalar).
