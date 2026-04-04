# Plot a `cavi` fit

Plot a `cavi` fit

## Usage

``` r
# S3 method for class 'cavi'
plot(
  x,
  plot_type = c("scatterplot", "elbo", "mu"),
  dims = c(1L, 2L),
  data = NULL,
  two_panel = FALSE,
  pal = (grDevices::colorRampPalette(c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00",
    "#FF0000")))(256L),
  add_legend = TRUE,
  ...
)
```

## Arguments

- x:

  A `cavi` object.

- plot_type:

  One of `"scatterplot"`, `"elbo"`, or `"mu"`.

- dims:

  Integer vector of length 1 or 2 used for plotting.

- data:

  Optional data matrix. Defaults to `x$data`.

- two_panel:

  Logical; for `plot_type = "elbo"`, show ELBO and plug-in
  log-likelihood in separate panels?

- pal:

  Colour palette for `"scatterplot"`.

- add_legend:

  Logical; draw the pseudotime legend on scatterplots?

- ...:

  Passed through to low-level plotting functions.

## Value

Invisibly returns `x`.
