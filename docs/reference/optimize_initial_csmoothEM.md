# Optimize csmoothEM initialization by comparing multiple methods after a warm start

Runs several initialization methods (via
[`parallel_initial_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/parallel_initial_csmoothEM.md)),
each for `num_iter` warm-start iterations, then selects the best fit.

Selection criterion:

- If `adaptive == "ml"`, selects the fit with the largest `ml_last`
  (collapsed objective \\\mathcal{C}\\).

- Otherwise, selects the fit with the largest `elbo_last`.

If `plot=TRUE`, plots traces across methods and highlights the selected
best method. When `adaptive=="ml"`, the primary plot shows `ml_trace`;
otherwise it shows `elbo_trace`. If `two_panel=TRUE`, a second panel
shows `loglik_trace`.

## Usage

``` r
optimize_initial_csmoothEM(
  X,
  methods = c("PCA", "tSNE", "random", "fiedler", "pcurve"),
  num_iter = 5,
  num_cores = 2,
  K = NULL,
  adaptive = "ml",
  lambda_min = 1e-08,
  lambda_max = 1e+08,
  sigma_update = c("ml", "mstep"),
  sigma_min = 1e-10,
  sigma_max = 1e+10,
  plot = FALSE,
  two_panel = FALSE,
  seed = NULL,
  quiet = TRUE,
  ...
)
```

## Arguments

- X:

  Numeric matrix `(n x d)`.

- methods:

  Character vector of initialization methods. Default includes `"PCA"`,
  `"tSNE"`, `"random"`, `"fiedler"`, `"pcurve"`, `"isomap"`.

- num_iter:

  Integer \\\ge 1\\. Number of warm-start iterations to run per method.

- num_cores:

  Integer \\\ge 1\\. Number of cores to use (see
  [`parallel_initial_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/parallel_initial_csmoothEM.md)).

- K:

  Optional integer \\\ge 2\\. Number of mixture components. If `NULL`,
  the default logic in `initialize_csmoothEM` is used.

- adaptive:

  Character. One of `"none"`, `"prior"`, `"ml"`. `"prior"` is obsolete
  and retained only for backward compatibility. Logical values are
  accepted for backward compatibility: `TRUE` is interpreted as `"ml"`
  and `FALSE` as `"none"`.

- lambda_min, lambda_max:

  Positive bounds for `lambda_vec` when `adaptive!="none"`.

- sigma_update:

  Character. Only used when `adaptive="ml"` (passed through). Defaults
  to `"ml"`; legacy `"mstep"` is retained only for backward
  compatibility.

- sigma_min, sigma_max:

  Positive bounds for `sigma2` when `adaptive="ml"` and
  `sigma_update="ml"`.

- plot:

  Logical; if `TRUE`, plot traces across methods.

- two_panel:

  Logical; if `TRUE`, also plot penalized objective traces in a second
  panel.

- seed:

  Optional integer seed. If provided, a different derived seed is used
  per method.

- quiet:

  Logical; reserved for future use.

- ...:

  Additional arguments passed to `initialize_csmoothEM` (and downstream
  ordering routines).

## Value

A `csmooth_em` object corresponding to the selected best method. The
returned object has an additional `meta$initial_search` list containing:

- `best_method`: selected method name

- `criterion`: selection criterion (`"ml_last"` or `"elbo_last"`)

- `summary`: the per-method summary data.frame

- `fits`: the full named list of fits

- `options`: options used for the search

## See also

[`parallel_initial_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/parallel_initial_csmoothEM.md),
[`initialize_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/initialize_csmoothEM.md)
