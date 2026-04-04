# Internal metadata for RW precision normalizing constants

Distinguishes the intrinsic RW(q) case from a properized/full-rank
precision (for example after adding a ridge nugget). In the proper case
the Gaussian prior normalizer uses the full matrix rank; in the
intrinsic case we retain the known RW(q) null-space dimension.

## Usage

``` r
.rw_precision_metadata(Q, rw_q = 0L, eigen_tol = NULL)
```

## Arguments

- Q:

  Symmetric precision matrix.

- rw_q:

  Expected RW(q) nullity in the intrinsic case.

- eigen_tol:

  Optional eigenvalue threshold used to decide whether `Q` is
  numerically full-rank.

## Value

A list with `rank`, `logdet`, `proper`, `nullity`, and `eigen_tol`.
