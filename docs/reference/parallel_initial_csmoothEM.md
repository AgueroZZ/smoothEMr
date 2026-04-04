# Run multiple csmoothEM initializations in parallel and summarize results

Runs
[`initialize_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/initialize_csmoothEM.md)
with multiple initialization `methods`, each for `num_iter` warm-start
iterations, and returns a named list of fits. A per-method summary table
is attached as an attribute `"summary"`.

The summary table reports the last values of `elbo_trace`,
`loglik_trace`, and `ml_trace`. In the current codebase, `ml_trace` is
the collapsed objective \\\mathcal{C}\\ (Laplace-exact in the Gaussian
csmooth setting).

## Usage

``` r
parallel_initial_csmoothEM(
  X,
  methods = c("PCA", "tSNE", "random", "fiedler", "pcurve"),
  num_iter = 1,
  num_cores = 2,
  K = NULL,
  adaptive = "ml",
  lambda_min = 1e-08,
  lambda_max = 1e+08,
  sigma_update = c("ml", "mstep"),
  sigma_min = 1e-10,
  sigma_max = 1e+10,
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
  `"tSNE"`, `"random"`, `"fiedler"`, `"pcurve"`, `"isomap"`. See
  [`initialize_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/initialize_csmoothEM.md)
  for details.

- num_iter:

  Integer \\\ge 1\\. Number of warm-start iterations to run inside
  `initialize_csmoothEM` for each method.

- num_cores:

  Integer \\\ge 1\\. Number of cores to use. On non-Windows systems,
  uses [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html);
  on Windows, falls back to a PSOCK cluster.

- K:

  Optional integer \\\ge 2\\. Number of mixture components. If `NULL`,
  the default logic in `initialize_csmoothEM` is used.

- adaptive:

  Character. One of `"none"`, `"prior"`, `"ml"`. `"prior"` is obsolete
  and retained only for backward compatibility. Logical values are
  accepted for backward compatibility: `TRUE` is interpreted as `"ml"`
  and `FALSE` as `"none"`.

- lambda_min, lambda_max:

  Positive bounds for `lambda_vec` (passed to `initialize_csmoothEM`).

- sigma_update:

  Character. Only used when `adaptive="ml"` (passed to
  `initialize_csmoothEM`). Defaults to `"ml"`; legacy `"mstep"` is
  retained only for backward compatibility. See
  [`do_csmoothEM_ml_collapsed`](https://aguerozz.github.io/MPCurver/reference/do_csmoothEM_ml_collapsed.md).

- sigma_min, sigma_max:

  Positive bounds for `sigma2` when `adaptive="ml"` and
  `sigma_update="ml"`.

- seed:

  Optional integer seed. If provided, a different derived seed is used
  per method.

- quiet:

  Logical; reserved for future use.

- ...:

  Additional arguments passed to `initialize_csmoothEM` (and downstream
  ordering routines).

## Value

A named list of fits, with names equal to `methods`. Each entry is
either a `csmooth_em` object (on success) or `NULL` (on failure). The
list has an attribute `"summary"` which is a data.frame with one row per
method:

- `method`: method name

- `success`: logical

- `n,d,K`: scalar integers from `summary(csmooth_em)` when available

- `elbo_last`: last penalized ELBO (or `NA`)

- `obj_last`: last penalized observed objective (or `NA`)

- `ml_last`: last collapsed objective \\\mathcal{C}\\ (or `NA`)

- `error`: error message if failed

## See also

[`optimize_initial_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/optimize_initial_csmoothEM.md),
[`initialize_csmoothEM`](https://aguerozz.github.io/MPCurver/reference/initialize_csmoothEM.md)
