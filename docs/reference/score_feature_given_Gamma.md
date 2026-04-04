# Score a single feature given responsibilities (Gamma)

Computes an alignment score for a single feature \\x_j\\ under a fixed
responsibilities matrix \\\Gamma \in \mathbb{R}^{n\times K}\\. The score
is used for greedy feature partitioning into multiple latent orderings.

Two score modes are supported:

- `score_mode = "none"`:

  Plug-in (ELBO/Q-like) score. The smoothing parameter \\\lambda\\ is
  fixed at `lambda_init` (default 1). The score depends on the MAP
  component means \\\hat\mu\_{j\cdot}\\ via the weighted SSE.

- `score_mode = "ml"`:

  Collapsed (marginal-like) score: plug-in score plus the Laplace
  curvature correction \\+\frac{K}{2}\log(2\pi) - \frac12\log\|A_j\|\\.
  In this mode, \\\lambda\\ is optimized by 1D maximization over \\\log
  \lambda \in \[\log(\text{lambda_min}),\log(\text{lambda_max})\]\\.

The fitted 1D quantities \\\hat\mu\_{j\cdot}\\ and \\\hat\sigma_j^2\\
are returned and can be appended to a partition's csmooth-style
parameters during greedy growth.

## Usage

``` r
score_feature_given_Gamma(
  xj,
  Gamma,
  Q_K,
  rw_q = 0L,
  score_mode = c("ml", "none"),
  relative_lambda = TRUE,
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  nugget = 0,
  optimize_lambda = NULL,
  lambda_init = 1,
  max_sigma_iter = 1L
)
```

## Arguments

- xj:

  Numeric vector of length `n` (one feature).

- Gamma:

  Numeric matrix `(n x K)` of responsibilities.

- Q_K:

  Numeric matrix `(K x K)`; base RW precision (lambda=1).

- rw_q:

  Integer \\\ge 0\\. Rank deficiency along K (RW order).

- score_mode:

  One of `"ml"` or `"none"`.

- relative_lambda:

  Logical; if TRUE use \\Q\_{base} = Q_K / \sigma_j^2\\ (homoskedastic
  scaling).

- lambda_min, lambda_max:

  Positive bounds for \\\lambda\\ when `score_mode="ml"`.

- nugget:

  Nonnegative scalar added to \\\sigma_j^2\\.

- optimize_lambda:

  Logical or NULL. If NULL, defaults to TRUE for `"ml"` and FALSE for
  `"none"`.

- lambda_init:

  Positive scalar. Fixed \\\lambda\\ when `score_mode="none"` (default
  1), and initial value if `optimize_lambda=FALSE`.

- max_sigma_iter:

  Integer \\\ge 1\\. Number of (lambda -\> mu -\> sigma2) refresh steps.

## Value

A list with components:

- `score`: scalar score for feature \\j\\.

- `lambda`: fitted (or fixed) \\\lambda_j\\.

- `sigma2`: fitted \\\sigma_j^2\\.

- `mu_vec`: numeric vector length K of \\\hat\mu\_{j\cdot}\\.

- `logdetA`: scalar \\\log\|A_j\|\\ (only meaningful for `"ml"`; still
  returned).

## See also

[`score_one_coord_csmooth`](https://aguerozz.github.io/MPCurver/reference/score_one_coord_csmooth.md),
[`forward_two_ordering_partition_csmooth`](https://aguerozz.github.io/MPCurver/reference/forward_two_ordering_partition_csmooth.md)
