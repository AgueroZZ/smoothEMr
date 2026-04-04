# Landmark Isomap ordering (no vegan; avoids O(n^2) dist matrix)

Landmark approximation to Isomap for fast ordering. Builds a weighted
kNN graph on rows of `X`, computes geodesic distances to a subset of
landmark points, embeds landmarks by classical MDS, then projects all
points by inverse-distance weighted averaging.

## Usage

``` r
isomap_ordering(
  X,
  k = 15,
  ndim = 1,
  component = 1,
  landmark = NULL,
  landmark_method = c("random", "kmeans"),
  seed = NULL,
  keep = c("giant", "all"),
  return_full = TRUE,
  orient_by_pc1 = TRUE,
  scale01 = TRUE,
  eps = 1e-08
)
```

## Arguments

- X:

  numeric matrix (n x D).

- k:

  number of nearest neighbors in the kNN graph (2 \<= k \< n).

- ndim:

  embedding dimension (\>= component).

- component:

  which embedding coordinate to use as ordering.

- landmark:

  integer; number of landmark points. If `NULL`, uses `min(1000, n)`.

- landmark_method:

  how to choose landmarks: `"random"` (default) or `"kmeans"`.

- seed:

  optional integer for reproducibility.

- keep:

  component handling: `"giant"` keeps the largest connected component;
  `"all"` uses all nodes (may yield Inf distances if disconnected).

- return_full:

  if TRUE, return length-n t with NA for excluded nodes (keep="giant").

- orient_by_pc1:

  logical; if TRUE, flip sign to align with PC1.

- scale01:

  logical; if TRUE, rescale returned t to \[0,1\].

- eps:

  small positive number to stabilize inverse-distance weights.

## Value

list(t, keep_idx, n_components, embed, landmark_idx,
geodesic_to_landmark)
