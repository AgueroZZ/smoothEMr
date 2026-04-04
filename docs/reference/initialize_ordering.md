# Initialize GMM parameters using an ordering method

Initialize GMM parameters using an ordering method

## Usage

``` r
initialize_ordering(
  X,
  K,
  method = c("PCA", "fiedler", "pcurve", "tSNE", "random", "isomap"),
  discretization = c("equal", "quantile", "kmeans"),
  ...
)
```

## Arguments

- X:

  numeric matrix (n x d).

- K:

  number of mixture components.

- method:

  ordering method.

- discretization:

  discretization method passed to make_init().

- ...:

  forwarded to the ordering method function.

## Value

list(pi, mu, sigma, keep_idx, cluster_rank, ordering)
