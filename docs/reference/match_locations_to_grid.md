# Match locations on an old grid to indices on a new grid

Match locations on an old grid to indices on a new grid

## Usage

``` r
match_locations_to_grid(loc_old, loc_new, tol = 1e-08)
```

## Arguments

- loc_old:

  Numeric vector of locations to be matched (e.g. u_obs).

- loc_new:

  Numeric vector of candidate locations (e.g. u_final).

- tol:

  Nonnegative numeric tolerance for matching.

## Value

Integer vector of indices into loc_new.
