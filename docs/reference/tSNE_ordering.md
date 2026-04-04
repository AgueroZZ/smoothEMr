# t-SNE ordering

t-SNE ordering

## Usage

``` r
tSNE_ordering(
  X,
  tSNE_dims = 1,
  component = 1,
  perplexity = 10,
  max_iter = 500,
  seed = NULL,
  scale01 = TRUE,
  orient_by_pc1 = TRUE
)
```

## Arguments

- X:

  numeric matrix (n x D).

- tSNE_dims:

  integer; embedding dimension passed to Rtsne.

- component:

  which embedding coordinate to use as ordering (default 1).

- perplexity:

  t-SNE perplexity; must satisfy \< (n-1)/3.

- max_iter:

  max iterations.

- seed:

  optional integer for reproducibility.

- scale01:

  logical; if TRUE, rescale returned t to \[0,1\].

- orient_by_pc1:

  logical; if TRUE, flip sign to align with PC1.

## Value

list(t, keep_idx)
