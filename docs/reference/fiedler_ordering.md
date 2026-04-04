# Fiedler ordering from a kNN graph

Builds a kNN graph on rows of `X`, computes a 1D ordering coordinate
using the Fiedler direction of a normalized graph operator, and returns
an ordering `t` scaled to \\\[0,1\]\\.

## Usage

``` r
fiedler_ordering(
  X,
  k = 15,
  weight = c("binary", "rbf", "inv"),
  sigma = NULL,
  keep = c("giant", "all"),
  return_full = TRUE
)
```

## Arguments

- X:

  Numeric matrix \\n \times D\\ (rows are observations).

- k:

  Number of nearest neighbors. If left at the default and the resulting
  kNN graph is disconnected, the function automatically increases `k`
  until the graph becomes connected.

- weight:

  Similarity type: `"rbf"`, `"inv"`, or `"binary"`.

- sigma:

  Bandwidth for `weight="rbf"`; if `NULL`, uses median kNN distance.

- keep:

  Component handling: `"giant"` keeps the largest connected component;
  `"all"` uses all nodes (may be unstable if disconnected).

- return_full:

  If `TRUE` (default), return `t` of length \\n\\ with `NA` for nodes
  excluded by `keep="giant"`. If `FALSE`, return `t` only for the kept
  nodes.

## Value

A list with `t`, `keep_idx`, `n_components`, and `k_used`.
