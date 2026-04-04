# Plot kriged mean curves from mu_full_history

Plot kriged mean curves from mu_full_history

## Usage

``` r
plot_mu_history(
  mu_full_history,
  u_final,
  history_i = 1,
  coords = c(1, 2, 3),
  main = NULL,
  xlab = "u (grid position)",
  ylab = expression(mu[j](u)),
  lty = 1,
  add_points = FALSE,
  u_obs = NULL,
  mu_obs_list = NULL,
  res = NULL,
  pch = 16,
  cex = 0.8,
  legend_loc = "topright",
  legend_prefix = "coord ",
  bty = "n"
)
```

## Arguments

- mu_full_history:

  List of matrices (K_final x d).

- u_final:

  Numeric length K_final.

- history_i:

  Which history element to plot.

- coords:

  Which coordinates (columns) to plot.

- main:

  Plot title.

- xlab:

  X-axis label.

- ylab:

  Y-axis label.

- lty:

  Line type for the kriged mean curves.

- add_points:

  Overlay observed means at design points.

- u_obs, mu_obs_list:

  Optional overrides for points.

- res:

  Optional progressive_smoothEM() result used to auto-fill points/title.

- pch:

  Point character for observed means.

- cex:

  Point expansion for observed means.

- legend_loc:

  Legend location.

- legend_prefix:

  Prefix for coordinate labels in the legend.

- bty:

  Legend box type.

## Value

Invisibly returns a list with plot inputs.
