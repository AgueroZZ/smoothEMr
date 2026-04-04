# Compare posterior position summaries across two histories

Compare posterior position summaries across two histories

## Usage

``` r
compare_positions_by_history(
  res,
  h1,
  h2,
  type = c("mean", "max"),
  add_identity = TRUE,
  main = NULL
)
```

## Arguments

- res:

  Result returned by
  [`progressive_smoothEM()`](https://aguerozz.github.io/MPCurver/reference/progressive_smoothEM.md).

- h1, h2:

  History indices to compare.

- type:

  Position summary to compare: posterior mean or MAP.

- add_identity:

  Logical; draw the identity line?

- main:

  Optional plot title.
