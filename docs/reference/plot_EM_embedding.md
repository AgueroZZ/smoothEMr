# Plot data embedding with SmoothEM component means

Visualize an EM fit (from
[`EM_algorithm()`](https://aguerozz.github.io/MPCurver/reference/EM_algorithm.md))
on 1D or 2D selected coordinates. For 2D, draws arrows connecting
component means in order. For 1D, draws a scatter of posterior position
vs the selected coordinate and overlays the mean curve.

## Usage

``` r
plot_EM_embedding(
  fit,
  X,
  dims = c(1, 2),
  position = NULL,
  use_posterior_mean = TRUE,
  pch = 19,
  col = "grey60",
  cex = 0.7,
  mu_pch = 8,
  mu_col = "orange",
  mu_cex = 1,
  arrow_col = "orange",
  arrow_lwd = 3,
  arrow_len = 0.08,
  line_col = "orange",
  line_lwd = 2,
  add = FALSE,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted object returned by
  [`EM_algorithm()`](https://aguerozz.github.io/MPCurver/reference/EM_algorithm.md).

- X:

  Numeric matrix (n x d). The data used to fit the model.

- dims:

  Integer vector of length 1 or 2 indicating which coordinates (columns
  of X) to plot.

- position:

  Optional numeric vector of length n giving x-axis "position". If NULL,
  uses `fit$position` if present; otherwise uses posterior mean position
  computed as `fit$gamma %*% seq_len(K)`.

- use_posterior_mean:

  Logical; if TRUE (default), use posterior mean position. If FALSE, use
  MAP position via `max.col(fit$gamma)`.

- pch, col, cex:

  Point style for data scatter.

- mu_pch, mu_col, mu_cex:

  Mean point style.

- arrow_col, arrow_lwd, arrow_len:

  Arrow style (2D only).

- line_col, line_lwd:

  Mean curve style (1D only).

- add:

  Logical; if TRUE, add to existing plot.

- xlab, ylab, main:

  Labels; if NULL, auto-generated.

- ...:

  Passed to [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  the scatter.

## Value

Invisibly returns a list with `dims`, `pos`, `mu_mat`.
