# Run csmoothEM with weighted inner updates

Identical interface to do_csmoothEM(), with one extra argument:
`feature_weights`. The inner loop reuses the same csmoothEM updates, but
replaces the responsibility update with
[`ESTEP_csmooth_weighted()`](https://aguerozz.github.io/MPCurver/reference/ESTEP_csmooth_weighted.md)
so all inner iterations use the correct weighted log-likelihood. When
`adaptive = "ml"`, the collapsed `lambda`/`sigma2` hyper-step is also
reweighted; the final MAP refresh of `mu` is unchanged.

## Usage

``` r
do_csmoothEM_weighted(
  object,
  data,
  feature_weights,
  iter = 3L,
  adaptive = "ml",
  sigma_update = c("ml", "mstep"),
  lambda_min = 1e-10,
  lambda_max = 1e+10,
  record = TRUE,
  verbose = FALSE
)
```

## Arguments

- object:

  A csmooth_em object.

- data:

  n-by-d data matrix.

- feature_weights:

  length-d weight vector in \[0,1\].

- iter:

  Number of inner EM iterations.

- adaptive:

  Adaptive lambda mode (`"ml"` or NULL). `"prior"` is obsolete and
  retained only for backward compatibility.

- sigma_update:

  Sigma update mode for `adaptive = "ml"`. Defaults to `"ml"`; legacy
  `"mstep"` is retained only for backward compatibility.

- lambda_min, lambda_max:

  Lambda bounds.

- record:

  Logical; record traces.

- verbose:

  Logical.

## Value

Updated csmooth_em object.
