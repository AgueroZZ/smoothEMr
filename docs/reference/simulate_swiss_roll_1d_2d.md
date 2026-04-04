# Simulate a 2D "swiss roll" spiral with a 1D latent parameter

Simulates a 2D spiral (often used as a 1D analogue of a swiss roll)
parameterized by \\t \in \[t\_{\min}, t\_{\max}\]\\. The noiseless curve
is \\(x(t),y(t)) = (t\cos t, t\sin t)\\. Independent Gaussian noise is
added to each coordinate.

## Usage

``` r
simulate_swiss_roll_1d_2d(
  n = 800,
  t_range = c(1.5 * base::pi, 6 * base::pi),
  sigma = 0.15,
  seed = 1
)
```

## Arguments

- n:

  Integer \\\ge 1\\. Number of samples.

- t_range:

  Numeric vector of length 2 specifying \\\[t\_{\min}, t\_{\max}\]\\.

- sigma:

  Nonnegative numeric. Standard deviation of Gaussian noise. May be a
  scalar (applied to both coordinates) or a length-2 vector for
  coordinate-specific noise.

- seed:

  Optional integer. If not NULL, sets the RNG seed via
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

A list with components:

- `t`: numeric vector of length `n`, the latent parameter.

- `truth`: data.frame with columns `x`, `y` for the noiseless curve.

- `obs`: data.frame with columns `x`, `y` for the noisy observations.

## Examples

``` r
sim <- simulate_swiss_roll_1d_2d(n = 300, sigma = 0.1, seed = 1)
head(sim$obs)
#>           x          y
#> 1 -4.818342   7.025793
#> 2 -8.512613  -5.259583
#> 3 12.397960   3.026388
#> 4  4.641363 -17.056883
#> 5  2.016847   7.101558
#> 6  2.224280 -17.250608
```
