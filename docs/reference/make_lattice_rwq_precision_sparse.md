# Sparse prior precision for q-th order random walk on a K x K lattice (d-dimensional)

Builds sparse 1D RW(q) precision Q1D = t(Dq) then constructs 2D lattice
precision via Kronecker sum: Q2D = kron(Q1D, I_K) + kron(I_K, Q1D), and
replicates across d independent processes: Q_U = lambda \* kron(Q2D,
I_d).

## Usage

``` r
make_lattice_rwq_precision_sparse(K, d, q = 1, lambda = 1, ridge = 0)
```

## Arguments

- K:

  integer lattice side length (total sites = K^2).

- d:

  integer number of independent processes.

- q:

  integer RW order (1 \<= q \< K).

- lambda:

  nonnegative scalar multiplier.

- ridge:

  nonnegative ridge added to Q1D (and thus Q2D) for conditioning.

## Value

sparse Matrix of size (d\*K^2) x (d\*K^2).
