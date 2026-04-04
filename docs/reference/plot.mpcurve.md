# Plot an `mpcurve` model fit

Visualise an `mpcurve` fit. The default (`plot_type="scatterplot"`)
colors each observation by its estimated pseudotime position along the
inferred trajectory, mapped to a \\\[0, 1\]\\ scale where component
\\k\\ occupies position \\(k-1)/(K-1)\\. The pseudotime for observation
\\i\\ is the responsibility-weighted average \$\$t_i = \sum\_{k=1}^{K}
\gamma\_{ik} \cdot \frac{k-1}{K-1}.\$\$ Component means are overlaid as
orange stars connected by arrows. A vertical gradient color bar is drawn
in the top-right corner of the plot.

- `"scatterplot"`:

  Scatter of two chosen dimensions, colored by pseudotime (default).

- `"elbo"`:

  ELBO / log-likelihood / collapsed-ML traces. Delegates to the
  underlying `plot.cavi`, `plot.csmooth_em`, or `plot.smooth_em`.

- `"mu"`:

  Component-means-only plot. Delegates to the underlying fit's plot
  method.

## Usage

``` r
# S3 method for class 'mpcurve'
plot(
  x,
  plot_type = c("scatterplot", "elbo", "mu"),
  dims = c(1L, 2L),
  data = NULL,
  pal = (grDevices::colorRampPalette(c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00",
    "#FF0000")))(256L),
  add_legend = TRUE,
  ...
)
```

## Arguments

- x:

  An `mpcurve` object.

- plot_type:

  One of `"scatterplot"` (default), `"elbo"`, `"mu"`.

- dims:

  Integer vector of length 2 (or 1). Columns of the data matrix to use
  for the scatter plot. Defaults to `c(1, 2)`.

- data:

  Optional numeric matrix (n x d). If `NULL` (default), uses
  `x$fit$data`. Must be supplied if the data were not stored at fit
  time.

- pal:

  Colour palette (length-256 character vector) used to map pseudotime to
  point colors. Defaults to the classic pseudotime rainbow (navy \\\to\\
  cyan \\\to\\ green \\\to\\ gold \\\to\\ red), which avoids white and
  is easy to read against a white background.

- add_legend:

  Logical; draw a gradient color-bar legend? Default `TRUE`.

- ...:

  Further arguments forwarded to `plot_EM_embedding2D` (for
  `"scatterplot"`) or the underlying algorithm's `plot` method (for
  `"elbo"` / `"mu"`).

## Value

Invisibly returns `x`.
