# Plot a csmooth_em object

Visualization for a `csmooth_em` object.

- `plot_type="scatterplot"`: calls
  [`plot_EM_embedding()`](https://aguerozz.github.io/MPCurver/reference/plot_EM_embedding.md)
  (same behavior as `smooth_em`).

- `plot_type="elbo"`: plots ELBO and penalized observed-data objective
  traces.

- `plot_type="mu"`: plots only the component means (with arrows for 2D).

## Usage

``` r
# S3 method for class 'csmooth_em'
plot(
  x,
  data = NULL,
  plot_type = c("scatterplot", "elbo", "mu"),
  dims = c(1, 2),
  two_panel = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A `csmooth_em` object.

- data:

  Numeric matrix (n x d). Required when `plot_type="scatterplot"` unless
  stored in `x$data`.

- plot_type:

  One of `"scatterplot"`, `"elbo"`, `"mu"`.

- dims:

  Integer vector of length 1 or 2. Used for `"scatterplot"` and `"mu"`.

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
