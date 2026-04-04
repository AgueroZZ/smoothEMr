# Scale random-walk penalty strength to account for grid spacing

Rescale `lambda` across nested equally-spaced grids on \[0,1\].

Using the common heuristic for RW(q): \$\$\lambda(h) \propto
h^{-(2q-1)}\$\$ and `h = 1/(K-1)`, we use: \$\$\lambda_level =
\lambda_final \* ((K_level-1)/(K_final-1))^{(2q-1)}\$\$.

## Usage

``` r
lambda_scale_for_spacing(lambda_final, K_final, K_level, q = 2)
```

## Arguments

- lambda_final:

  Scalar lambda used on final grid.

- K_final:

  Number of points on the final grid.

- K_level:

  Number of points on the current level grid.

- q:

  RW order (e.g. 2 for RW2).

## Value

Scalar `lambda_level`.
