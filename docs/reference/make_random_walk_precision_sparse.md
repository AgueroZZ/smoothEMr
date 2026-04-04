# Sparse prior precision for q-th order random walk on K time points (d-dimensional)

Builds the q-th forward-difference matrix D (size (K-q) x K) with
coefficients (-1)^(q-j) \* choose(q, j), then returns: Q_U = lambda \*
(t(D)

## Usage

``` r
make_random_walk_precision_sparse(K, d, q = 1, lambda = 1, ridge = 0)
```

## Arguments

- K:

  integer number of time points.

- d:

  integer dimension per time point.

- q:

  integer RW order (1 \<= q \< K).

- lambda:

  nonnegative scalar multiplier.

- ridge:

  nonnegative scalar ridge added to Q1D.

## Value

sparse Matrix of size (d\*K) x (d\*K).
