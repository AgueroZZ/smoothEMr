# Prior precision for q-th order random walk on a K x K lattice (d-dimensional)

Builds a 2D intrinsic RW(q) precision via Kronecker sum: Q2D = kron(Q1D,
I) + kron(I, Q1D), then replicates across d independent processes.

## Usage

``` r
make_lattice_rwq_precision(K, d, q = 1, lambda = 1, ridge = 0)
```

## Arguments

- K:

  integer lattice side length (total sites = K^2).

- d:

  integer number of independent processes.

- q:

  integer order of the random walk (q\>=1).

- lambda:

  nonnegative scalar precision multiplier.

- ridge:

  nonnegative ridge added to improve conditioning.

## Value

sparse (d\*K^2) x (d\*K^2) precision matrix (Matrix).
