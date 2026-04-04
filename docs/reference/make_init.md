# Initialization from an ordering vector

Initialization from an ordering vector

## Usage

``` r
make_init(
  X,
  ordering_vec,
  K,
  assume_EEI = TRUE,
  nugget = 0,
  discretization = c("equal", "quantile", "kmeans"),
  na_action = c("drop", "error")
)
```

## Arguments

- X:

  numeric matrix (n x d).

- ordering_vec:

  numeric vector of length n (can contain NA).

- K:

  number of mixture components.

- assume_EEI:

  logical; if TRUE, use shared diagonal covariance (EEI-like).

- nugget:

  nonnegative scalar added to diagonal variances/covariances.

- discretization:

  one of "equal", "quantile", "kmeans".

- na_action:

  how to handle NA in ordering_vec: "drop" or "error".

## Value

list(pi, mu, sigma, keep_idx, cluster_rank)
