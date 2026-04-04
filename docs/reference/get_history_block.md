# Extract a (u, gamma) block for a given progressive history index

Extract a (u, gamma) block for a given progressive history index

## Usage

``` r
get_history_block(res, h, normalize_gamma = TRUE)
```

## Arguments

- res:

  Result returned by
  [`progressive_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/progressive_smoothEM.md).

- h:

  History index to extract.

- normalize_gamma:

  Logical; renormalize rows of the returned gamma matrix?
