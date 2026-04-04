# PCA ordering (PC1)

PCA ordering (PC1)

## Usage

``` r
PCA_ordering(X, center = TRUE, scale = FALSE, scale01 = TRUE, component = 1L)
```

## Arguments

- X:

  numeric matrix (n x D).

- center:

  logical.

- scale:

  logical.

- scale01:

  logical; if TRUE, rescale returned t to \[0,1\].

- component:

  Integer principal-component index to use as the ordering score.
  Defaults to 1.

## Value

list(t, keep_idx)
