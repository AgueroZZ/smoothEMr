# Prior precision for q-th order random walk on K time points (d-dimensional)

Constructs Q = lambda \* (D_q' D_q) kronecker I_d, with optional ridge.

## Usage

``` r
make_random_walk_precision(K, d, q = 1, lambda = 1, ridge = 0, sparse = FALSE)
```

## Arguments

- K:

  integer number of time points.

- d:

  integer dimension of the process at each time.

- q:

  integer order of the random walk (q\>=1).

- lambda:

  nonnegative scalar precision multiplier.

- ridge:

  nonnegative scalar ridge added to the 1D precision.

- sparse:

  logical; return a sparse Matrix object if TRUE.

## Value

(dK) x (dK) precision matrix.
