# Krige SmoothEM mean functions from an observed grid to the final grid

Krige SmoothEM mean functions from an observed grid to the final grid

## Usage

``` r
krige_mu_list_to_full_grid(mu_list_obs, u_obs, u_final, Q_final_1d, nugget = 0)
```

## Arguments

- mu_list_obs:

  List length K_obs; each element a numeric vector length d.

- u_obs:

  Numeric length K_obs.

- u_final:

  Numeric length K_final.

- Q_final_1d:

  K_final x K_final precision acting along the grid.

- nugget:

  Nonnegative scalar passed to kriging_from_precision().

## Value

List(Mu_full, mu_full_list, idx_obs).
