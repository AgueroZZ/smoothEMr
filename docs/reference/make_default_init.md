# Default random initialization for a Gaussian mixture model

Default random initialization for a Gaussian mixture model

## Usage

``` r
make_default_init(X, K, ordering = TRUE)
```

## Arguments

- X:

  numeric matrix (n x d).

- K:

  number of mixture components.

- ordering:

  if TRUE, reorder components by a pivot dimension of mu.

## Value

list(pi, mu, sigma)
