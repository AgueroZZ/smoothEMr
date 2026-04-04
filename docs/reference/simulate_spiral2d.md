# Simulate a 2D Archimedean spiral (helix-like) with noise

Simulates a 2D Archimedean spiral curve parameterized by a latent
ordering variable \\t \in \[0,1\]\\. The radius increases linearly with
\\t\\ and the angle makes `turns` revolutions. Independent Gaussian
noise is added to each coordinate.

## Usage

``` r
simulate_spiral2d(n = 500, turns = 3, noise = 0.05, r0 = 0.2, r1 = 1, seed = 1)
```

## Arguments

- n:

  Integer \\\ge 1\\. Number of samples.

- turns:

  Positive numeric. Number of revolutions over \\t \in \[0,1\]\\.

- noise:

  Nonnegative numeric. Standard deviation of Gaussian noise added to
  both coordinates.

- r0, r1:

  Nonnegative numerics with `r1 >= r0`. Start/end radius.

- seed:

  Optional integer. If not NULL, sets the RNG seed via
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

A list with components:

- `obs`: numeric matrix `(n x 2)` with column names `c("x1","x2")`.

- `t`: numeric vector of length `n`, the latent ordering in \\\[0,1\]\\.

## Examples

``` r
sim <- simulate_spiral2d(n = 200, turns = 4, noise = 0.05, seed = 1)
head(sim$obs)
#>               x1         x2
#> [1,]  0.16817785 0.11261864
#> [2,]  0.20099800 0.01721633
#> [3,]  0.13659109 0.21956446
#> [4,]  0.15110509 0.15879060
#> [5,] -0.01061102 0.32886305
#> [6,]  0.09283838 0.32499952
head(sim$t)
#> [1] 0.01307758 0.01339033 0.02333120 0.03554058 0.05893438 0.06178627
```
