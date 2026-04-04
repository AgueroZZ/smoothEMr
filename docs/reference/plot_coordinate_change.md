# Plot change in a single coordinate between two histories

Plot change in a single coordinate between two histories

## Usage

``` r
plot_coordinate_change(
  mu_full_history,
  u_final,
  coord,
  history_i,
  history_j,
  type = c("both", "overlay", "diff"),
  highlight = c("both", "i", "j", "none"),
  highlight_u = NULL,
  res = NULL,
  main = NULL,
  xlab = "u (grid position)",
  ylab_overlay = expression(mu[j](u)),
  ylab_diff = expression(Delta * mu[j](u)),
  lty_i = 1,
  lty_j = 2,
  lty_diff = 1,
  pch_i = 16,
  pch_j = 17,
  pch_diff = 16,
  cex_pts = 0.9,
  add_zero_line = TRUE,
  legend_loc = "topright",
  bty = "n",
  ylim_overlay = NULL,
  pad_ylim = 0.04,
  include_highlight_in_ylim = TRUE,
  label_points = FALSE
)
```

## Arguments

- mu_full_history:

  List of kriged mean matrices.

- u_final:

  Final grid locations.

- coord:

  Coordinate index to plot.

- history_i:

  First history index.

- history_j:

  Second history index.

- type:

  One of `"both"`, `"overlay"`, or `"diff"`.

- highlight:

  Which history's observed locations to highlight.

- highlight_u:

  Optional explicit locations to highlight.

- res:

  Optional
  [`progressive_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/progressive_smoothEM.md)
  result used to recover observed design locations.

- main:

  Plot title.

- xlab:

  X-axis label.

- ylab_overlay:

  Y-axis label for the overlay plot.

- ylab_diff:

  Y-axis label for the difference plot.

- lty_i, lty_j, lty_diff:

  Line types for history i, history j, and the difference curve.

- pch_i, pch_j, pch_diff:

  Point characters for highlighted locations.

- cex_pts:

  Point expansion for highlighted locations.

- add_zero_line:

  Logical; add a horizontal zero line on the difference plot?

- legend_loc:

  Legend location.

- bty:

  Legend box type.

- ylim_overlay:

  Optional y-limits for the overlay panel.

- pad_ylim:

  Fractional padding applied when computing overlay y-limits.

- include_highlight_in_ylim:

  Logical; should highlighted points affect the computed overlay
  y-limits?

- label_points:

  Logical; label highlighted points by grid index?
