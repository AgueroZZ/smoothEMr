# Plot soft feature weights (4-panel diagnostic)

Plot soft feature weights (4-panel diagnostic)

## Usage

``` r
plot_soft_weights(
  result,
  feature_names = NULL,
  show_convergence = TRUE,
  sort_by_weight = TRUE,
  title = "Soft feature partition"
)
```

## Arguments

- result:

  Output of soft_two_trajectory_EM().

- feature_names:

  Optional character vector.

- show_convergence:

  Show convergence panels 2, 3 & 4.

- sort_by_weight:

  Sort features by w_A descending in bar chart.

- title:

  Plot title.
