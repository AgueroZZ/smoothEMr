# Log-prior contribution induced by an exponential prior on 1/sqrt(lambda)

For \\\tau = 1 / \sqrt{\lambda}\\ with \\\tau \sim \mathrm{Exp}(\rho)\\,
the induced prior on \\\lambda \> 0\\ contributes \$\$ \log p(\lambda) =
\log(\rho / 2) - \frac{3}{2}\log \lambda - \rho \lambda^{-1/2}. \$\$

## Usage

``` r
.lambda_sd_prior_terms(lambda_vec, rate, include_constant = TRUE)
```

## Arguments

- lambda_vec:

  Positive numeric vector.

- rate:

  Non-negative scalar rate parameter \\\rho\\.

- include_constant:

  Logical; include \\\log(\rho / 2)\\ when `rate > 0`?

## Value

Numeric vector of log-prior contributions.
