# Principal curve ordering

Principal curve ordering

## Usage

``` r
pcurve_ordering(
  X,
  smoother = c("smooth_spline", "lowess", "periodic_lowess"),
  thresh = 0.001,
  maxit = 10,
  stretch = 2,
  approx_points = FALSE,
  scale01 = TRUE,
  orient_by_pc1 = TRUE
)
```

## Arguments

- X:

  numeric matrix (n x D).

- smoother:

  a function (recommended) or a name like "smooth_spline".

- thresh, maxit, stretch, approx_points:

  passed to princurve::principal_curve.

- scale01:

  logical; if TRUE, rescale returned t to \[0,1\].

- orient_by_pc1:

  logical; if TRUE, flip sign to align with PC1.

## Value

list(t, keep_idx, fit)
