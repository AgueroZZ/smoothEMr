# Prior precision for VAR(1) on K time points (d-dimensional)

Model: \\U_1 \sim N(0, P0^{-1})\\, and for \\k \>= 2\\, \\U_k \mid
U\_{k-1} \sim N(A U\_{k-1}, Q)\\. Returns the joint precision matrix for
\\\mathrm{vec}(U_1, \dots, U_K)\\.

## Usage

``` r
make_VAR1_precision(K, d, A = NULL, Q = NULL, P0 = NULL)
```

## Arguments

- K:

  integer number of time points.

- d:

  integer state dimension.

- A:

  d x d transition matrix (default 0.8 \* I_d).

- Q:

  d x d innovation covariance (default 0.5 \* I_d).

- P0:

  d x d initial precision (default I_d).

## Value

(dK) x (dK) precision matrix (dense).
