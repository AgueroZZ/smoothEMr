# Choose the best SmoothEM initialization by ELBO

Runs
[`parallel_initial()`](https://aguerozz.github.io/MPCurver/reference/parallel_initial.md)
and returns the `smooth_em` object with the largest last-iteration ELBO.

If `plot=TRUE`, plots all ELBO traces (different colors) and optionally
overlays the penalized objective traces in a second panel.

## Usage

``` r
optimize_initial(
  X,
  methods = c("PCA", "tSNE", "random", "fiedler", "multi_scale", "pcurve"),
  num_iter = 1,
  num_cores = 2,
  m_max = 6,
  K = NULL,
  adaptive = TRUE,
  lambda_min = 1e-08,
  lambda_max = 1e+08,
  plot = FALSE,
  two_panel = FALSE,
  seed = NULL,
  quiet = TRUE,
  ...
)
```

## Arguments

- X:

  Numeric matrix (n x d).

- methods:

  Methods to try (passed to
  [`parallel_initial()`](https://aguerozz.github.io/MPCurver/reference/parallel_initial.md)).

- num_iter:

  Total number of iterations for each fit.

- num_cores:

  Number of cores.

- m_max:

  Used for multi_scale and default K for others.

- K:

  Optional K for non-multi-scale methods.

- adaptive:

  Logical; passed to
  [`parallel_initial()`](https://aguerozz.github.io/MPCurver/reference/parallel_initial.md)
  /
  [`initialize_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/initialize_smoothEM.md).

- lambda_min, lambda_max:

  Positive bounds for lambda when `adaptive=TRUE`.

- plot:

  Logical; if TRUE, plot traces.

- two_panel:

  Logical; if TRUE, show ELBO and objective in 2 panels.

- seed:

  Optional base seed.

- quiet:

  Logical.

- ...:

  Passed to
  [`initialize_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/initialize_smoothEM.md)
  via
  [`parallel_initial()`](https://aguerozz.github.io/MPCurver/reference/parallel_initial.md).

## Value

A `smooth_em` object (best by last ELBO). The returned object gains:

- `$meta$initial_search`: list with the full fits, summary table, and
  options.
