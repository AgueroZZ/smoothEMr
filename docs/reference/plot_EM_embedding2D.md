# Plot a 2D embedding with Smooth-EM component means overlaid

Plot an \\n \times 2\\ embedding (e.g., two selected dimensions,
PCA/t-SNE/UMAP) and overlay the Smooth-EM component means (`mu_list`) as
points connected by arrows. Optionally color observations by a
continuous vector (e.g., a Fiedler vector).

## Usage

``` r
plot_EM_embedding2D(
  mu_list,
  X2,
  t_vec = NULL,
  mu_order = NULL,
  pch = 19,
  col = "grey60",
  cex = 0.7,
  mu_pch = 8,
  mu_col = "orange",
  mu_cex = 1,
  arrow_col = "orange",
  arrow_lwd = 2.5,
  arrow_len = 0.08,
  add = FALSE,
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  add_legend = !is.null(t_vec),
  legend_title = "t",
  pal = (grDevices::colorRampPalette(c("#2b8cbe", "white", "#de2d26")))(256),
  ...
)
```

## Arguments

- mu_list:

  A length-\\K\\ list of mean vectors. Each element should be a numeric
  vector of length at least 2, where the first two entries correspond to
  the coordinates in `X2`.

- X2:

  Numeric matrix of dimension \\n \times 2\\ giving the embedding
  coordinates for observations.

- t_vec:

  Optional numeric vector of length \\n\\. If provided, points are
  colored by `t_vec` using `pal`.

- mu_order:

  Optional integer permutation of `1:K` specifying the plotting order of
  components. Defaults to `1:K`.

- pch:

  Point character for observations (passed to
  [`points`](https://rdrr.io/r/graphics/points.html) /
  [`plot`](https://rdrr.io/r/graphics/plot.default.html)).

- col:

  Observation color used when `t_vec` is `NULL`.

- cex:

  Point size for observations.

- mu_pch:

  Point character for component means.

- mu_col:

  Color for component means.

- mu_cex:

  Point size for component means.

- arrow_col:

  Color for arrows connecting consecutive component means.

- arrow_lwd:

  Line width for arrows.

- arrow_len:

  Arrow head length (see
  [`arrows`](https://rdrr.io/r/graphics/arrows.html)).

- add:

  Logical; if `TRUE`, add to an existing plot. If `FALSE`, create a new
  plot.

- xlab, ylab:

  Axis labels. If `NULL`, defaults to `"dim 1"` and `"dim 2"`.

- main:

  Plot title. If `NULL`, uses a default title.

- add_legend:

  Logical; if `TRUE` and `t_vec` is provided, add a small legend
  (quantile ticks).

- legend_title:

  Title for the legend when `t_vec` is provided.

- pal:

  A vector of colors used to map `t_vec` to point colors. Defaults to a
  blue-white-red palette.

- ...:

  Additional arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) when
  `add = FALSE`.

## Value

An invisible list with elements:

- `X2`: the embedding coordinates used for plotting.

- `mu`: the \\K \times 2\\ matrix of plotted component means (after
  `mu_order`).

- `mu_order`: the component order used.

## Examples

``` r
set.seed(1)
n <- 50
X2 <- cbind(rnorm(n), rnorm(n))
mu_list <- list(c(-1, -1), c(0, 0), c(1, 1))
t_vec <- rnorm(n)

plot_EM_embedding2D(mu_list, X2)

plot_EM_embedding2D(mu_list, X2, t_vec = t_vec, legend_title = "Fiedler")

```
