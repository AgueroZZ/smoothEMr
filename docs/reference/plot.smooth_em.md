# Plot a smooth_em object

Visualization for a `smooth_em` object.

- `plot_type="scatterplot"`: calls
  [`plot_EM_embedding()`](https://aguerozz.github.io/MPCurver/reference/plot_EM_embedding.md).

- `plot_type="elbo"`: plots ELBO and penalized observed-data objective
  traces.

## Usage

``` r
# S3 method for class 'smooth_em'
plot(
  x,
  data = NULL,
  plot_type = c("scatterplot", "elbo"),
  dims = c(1, 2),
  two_panel = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A `smooth_em` object.

- data:

  Numeric matrix (n x d). Required when `plot_type="scatterplot"`.

- plot_type:

  One of `"scatterplot"`, `"elbo"`.

- dims:

  Integer vector of length 1 or 2 (only used for `"scatterplot"`).

- two_panel:

  Logical; if TRUE and `plot_type="elbo"`, draw ELBO and objective in
  two panels.

- verbose:

  Logical; not passed to graphics (prevents "not a graphical parameter"
  warnings).

- ...:

  Passed to the underlying plotting functions.

## Value

Invisibly returns `x`.
