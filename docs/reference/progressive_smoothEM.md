# Progressive-resolution initialization for SmoothEM via kriging on a final grid

Coarse-to-fine continuation scheme for initializing SmoothEM means over
nested grids.

## Usage

``` r
progressive_smoothEM(
  data,
  m_max = 6,
  lambda_final = 500,
  q = 2,
  ridge = 0,
  nugget_kriging = 0,
  tol = 0.001,
  max_iter = 1000,
  relative_lambda = TRUE,
  modelName = "EEI",
  coords_show = c(1, 2, 3),
  plot_each_stage = TRUE,
  verbose = TRUE,
  include.data = TRUE
)
```

## Arguments

- data:

  Numeric matrix n x d.

- m_max:

  Finest grid exponent; final grid size is 2^m_max + 1.

- lambda_final:

  Penalty strength on final grid.

- q:

  RW order (e.g. 2 for RW2).

- ridge:

  Ridge added in RW precision construction.

- nugget_kriging:

  Nugget added to Q_UU during kriging solve.

- tol:

  Convergence tolerance on ELBO changes.

- max_iter:

  Maximum number of EM iterations.

- relative_lambda:

  Logical; if TRUE, rescales `Q_prior` by current marginal variances.

- modelName:

  Covariance model: one of `"VVV"`, `"VII"`, `"EII"`, `"EEI"`.

- coords_show:

  Coords to visualize when plot_each_stage=TRUE.

- plot_each_stage:

  If TRUE, plot kriged mean curves each stage.

- verbose:

  If TRUE, print stage summaries and forward verbose to EM_algorithm().

- include.data:

  If TRUE, include data in returned fits.

## Value

List(grid, Q_final_1d, fits, mu_full_history, mu_full_list_final,
meta_history, mclust).
