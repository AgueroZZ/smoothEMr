# Plot a 2D embedding colored by an ordering vector with EM means overlaid

Given an \\n \times 2\\ embedding (`coords`) and a continuous ordering
vector (`order_vec`), this function bins observations by `order_vec`
(either by quantiles or equal-width bins) and colors points by the
resulting bins. It then overlays SmoothEM component means (`mu_list`) as
points connected by arrows.

## Usage

``` r
plot_order_EM_overlay2D(
  coords,
  order_vec,
  mu_list,
  method = c("auto", "quantiles", "seq"),
  probs = seq(0, 1, by = 0.1),
  nbins = 30,
  include_lowest = TRUE,
  right = TRUE,
  only = NULL,
  fixed_limits = TRUE,
  pal = NULL,
  pch = 16,
  cex = 0.7,
  add_centroid = FALSE,
  centroid_pch = 4,
  centroid_cex = 1.2,
  centroid_lwd = 2,
  centroid_labels = FALSE,
  label_cex = 0.8,
  label_pos = 3,
  mu_order = NULL,
  mu_pch = 8,
  mu_col = "orange",
  mu_cex = 1,
  arrow_col = "orange",
  arrow_lwd = 2.5,
  arrow_len = 0.08,
  add_legend = FALSE,
  legend_loc = "none",
  legend_cex = 0.7,
  legend_title = NULL,
  asp = 1,
  xlab = "dim1",
  ylab = "dim2",
  main = "Ordering + EM means",
  ...
)
```

## Arguments

- coords:

  Numeric matrix of dimension \\n \times 2\\.

- order_vec:

  Numeric vector of length \\n\\. Any continuous ordering/score vector.

- mu_list:

  A length-\\K\\ list of mean vectors; each element must have at least 2
  entries. The first two entries are used as coordinates in the same 2D
  space as `coords`.

- method:

  Binning method for `order_vec`: `"quantiles"` or `"seq"`. If `"auto"`,
  defaults to `"seq"`.

- probs:

  Quantile breakpoints (used when `method = "quantiles"`).

- nbins:

  Number of equal-width bins (used when `method = "seq"`).

- include_lowest, right:

  Passed to [`cut`](https://rdrr.io/r/base/cut.html).

- only:

  Optional subset of bins to plot; either integer bin indices or bin
  labels.

- fixed_limits:

  Logical; if `TRUE`, keep `xlim`/`ylim` fixed to the full `coords`
  range even when `only` is used.

- pal:

  Optional vector of colors (length = number of plotted bins). If
  `NULL`, a blue-white-red palette is used.

- pch, cex:

  Point style for observations.

- add_centroid:

  Logical; if `TRUE`, add per-bin centroids (computed on plotted
  subset).

- centroid_pch, centroid_cex, centroid_lwd, centroid_labels, label_cex,
  label_pos:

  Centroid styling.

- mu_order:

  Optional permutation of `1:K` controlling the overlay order of
  `mu_list`.

- mu_pch, mu_col, mu_cex:

  Styling for mean points.

- arrow_col, arrow_lwd, arrow_len:

  Styling for arrows.

- add_legend, legend_loc, legend_cex, legend_title:

  Legend controls. Set `legend_loc = "none"` (or `NULL`) to suppress
  legend.

- asp, xlab, ylab, main:

  Plot controls passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

- ...:

  Additional arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns a list with kept indices, bins, colors, optional bin
centroids, and the plotted mean coordinates.
