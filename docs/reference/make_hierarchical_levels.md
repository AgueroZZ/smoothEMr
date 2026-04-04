# Construct a hierarchy of nested grids inside a final grid

Build nested grid levels indexed into a final grid. Final grid has
`K_final = 2^m_max + 1` equally-spaced points on \[0,1\]. Level `m` has
`K_m = 2^m + 1` points, appearing as a subset of final indices.

## Usage

``` r
make_hierarchical_levels(m_max = 6)
```

## Arguments

- m_max:

  Nonnegative integer; the finest level exponent.

## Value

List with `m_max`, `K_final`, `u_final`, `idx_levels`, `K_levels`.
