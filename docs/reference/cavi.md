# Fit the recommended CAVI model for a single MPCurver ordering

Recommended variational pipeline for a single-ordering MPCurver model
with feature-specific noise levels \\\sigma_j^2\\ and smoothness
parameters \\\lambda_j\\. The latent trajectory for feature \\j\\ is a
length-\\K\\ vector with a GMRF prior along the component axis.

This function is now the main single-ordering MPCurver backend. In
contrast to the legacy `csmooth_em` path, it keeps an explicit Gaussian
posterior over feature trajectories, with posterior means \\m_j\\ and
covariances \\S_j\\.

## Usage

``` r
cavi(
  X,
  K = NULL,
  method = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
  responsibilities_init = NULL,
  pi_init = NULL,
  sigma2_init = NULL,
  lambda_init = NULL,
  rw_q = 2L,
  ridge = 0,
  discretization = c("quantile", "equal", "kmeans"),
  fix_lambda = FALSE,
  lambda_sd_prior_rate = NULL,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  sigma_min = 1e-10,
  sigma_max = 1e+10,
  max_iter = 100L,
  tol = 1e-06,
  verbose = FALSE
)
```

## Arguments

- X:

  Numeric matrix (\\n \times d\\).

- K:

  Number of components. If NULL, defaults to
  `max(2, min(50, floor(n/5)))`.

- method:

  Initialization ordering method when `responsibilities_init` is NULL.

- responsibilities_init:

  Optional initial responsibility matrix (\\n \times K\\).

- pi_init:

  Optional initial mixture proportions.

- sigma2_init:

  Optional initial feature-specific variances.

- lambda_init:

  Optional initial feature-specific smoothness vector.

- rw_q:

  Random-walk order for the prior precision along components.

- ridge:

  Ridge added to `Q_K`.

- discretization:

  Discretization method for ordering-based initialization.

- fix_lambda:

  Logical. If `TRUE`, skip variational lambda updates and keep
  `lambda_j` fixed at the initial values throughout fitting.

- lambda_sd_prior_rate:

  Optional positive rate for an exponential prior on
  `1 / sqrt(lambda_j)`. The default `NULL` means no lambda prior
  penalty. For backward compatibility, an explicit `0` is treated the
  same way; it is only an alias for "no penalty" and does not correspond
  to a literal exponential prior with rate zero.

- lambda_min, lambda_max:

  Bounds for `lambda_j`.

- sigma_min, sigma_max:

  Bounds for `sigma_j^2`.

- max_iter:

  Maximum number of CAVI sweeps. `0` returns the initialization-only
  state.

- tol:

  Relative ELBO tolerance.

- verbose:

  Logical.

## Value

An object of class `"cavi"` with components including:

- `$params`: current `pi`, `mu`, and `sigma2`

- `$posterior`: Gaussian posterior summaries for feature trajectories

- `$gamma`: cell-by-component responsibilities

- `$lambda_vec`: feature-specific smoothness parameters

- `$elbo_trace`: the monotone variational objective used for diagnostics

- `$loglik_trace`: a plug-in observed log-likelihood diagnostic
