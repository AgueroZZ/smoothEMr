# Run multiple SmoothEM initializations in parallel

Runs
[`initialize_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/initialize_smoothEM.md)
for a set of initialization methods in parallel. Each method is fit for
`num_iter` total SmoothEM iterations. Internally,
[`initialize_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/initialize_smoothEM.md)
warm-starts with exactly one EM iteration, and if `num_iter > 1`,
continues via
[`do_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/do_smoothEM.md)
(where adaptive lambda updates may be enabled).

For compatibility with `method="multi_scale"`, if `K` is not provided,
then non-multi-scale methods default to `K = 2^m_max + 1`.

## Usage

``` r
parallel_initial(
  X,
  methods = c("PCA", "tSNE", "random", "fiedler", "multi_scale", "pcurve"),
  num_iter = 1,
  num_cores = 2,
  m_max = 6,
  K = NULL,
  seed = NULL,
  adaptive = TRUE,
  lambda_min = 1e-08,
  lambda_max = 1e+08,
  quiet = TRUE,
  ...
)
```

## Arguments

- X:

  Numeric matrix (n x d).

- methods:

  Character vector of methods to try. Defaults to
  `c("PCA","tSNE","random","fiedler","multi_scale")`.

- num_iter:

  Integer \>= 1. Total number of SmoothEM outer iterations to run for
  each method.

- num_cores:

  Integer \>= 1. Number of cores for parallel execution.

- m_max:

  Integer \>= 1. Used by `multi_scale` and to set default `K`.

- K:

  Optional integer grid size for non-multi-scale methods. If NULL, uses
  `2^m_max+1`.

- seed:

  Optional base seed for reproducibility. If provided, each method gets
  a deterministic derived seed.

- adaptive:

  Logical; if TRUE, enables adaptive lambda behavior in
  [`initialize_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/initialize_smoothEM.md).

- lambda_min, lambda_max:

  Positive bounds used to clip lambda when `adaptive=TRUE`.

- quiet:

  Logical; suppress messages from workers.

- ...:

  Extra args passed to
  [`initialize_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/initialize_smoothEM.md)
  (e.g. `rw_q`, `lambda`, `relative_lambda`, `modelName`, etc.).

## Value

A named list of `smooth_em` objects (or `NULL` for failed fits), with
attributes:

- `summary`: a data.frame summarizing last ELBO / last objective / last
  lambda for each method.
