# Initialize two CAVI trajectories for dual-ordering partitioning

Initialize two CAVI trajectories for dual-ordering partitioning

## Usage

``` r
init_two_trajectories_cavi(
  X,
  method = c("score", "mincor", "random_split"),
  method1 = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
  K = NULL,
  rw_q = 2L,
  ridge = 0,
  lambda_sd_prior_rate = NULL,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  sigma_min = 1e-10,
  sigma_max = 1e+10,
  discretization = c("quantile", "equal", "kmeans"),
  num_iter = 5L,
  seed = 42L,
  verbose = FALSE
)
```

## Arguments

- X:

  Numeric matrix (`n x d`).

- method:

  Initialisation strategy for the second trajectory.

- method1:

  Initialisation method for the first trajectory.

- K:

  Number of pseudotime bins.

- rw_q:

  Random-walk order.

- ridge:

  Optional nugget added to the RW precision. Mainly useful for internal
  diagnostics that compare intrinsic and properized priors.

- lambda_sd_prior_rate:

  Optional positive rate for an exponential prior on
  `1 / sqrt(lambda_j^{(m)})`. The default `NULL` means no lambda prior
  penalty. For backward compatibility, an explicit `0` is treated the
  same way; it is only an alias for "no penalty" and does not correspond
  to a literal exponential prior with rate zero.

- lambda_min, lambda_max:

  Bounds used when clipping the initial feature-specific `lambda_j`
  values for each warm-start fit.

- sigma_min, sigma_max:

  Bounds for the feature-specific noise variance `sigma2_j`. Forwarded
  to internal
  [`cavi()`](https://aguerozz.github.io/MPCurver/reference/cavi.md)
  calls.

- discretization:

  Initial discretization scheme used when converting ordering scores
  into `K` bins. For dual-ordering initialisation the requested `K` is
  enforced strictly; if quantile cuts collapse, the code falls back to
  equal-width bins so both orderings remain comparable.

- num_iter:

  Warm-start CAVI sweeps per fit.

- seed:

  Random seed used by `"random_split"`.

- verbose:

  Logical.

## Value

A list with `fit1`, `fit2`, and `seed2`.
