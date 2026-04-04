# Optimize lambda under the induced exponential prior on 1/sqrt(lambda)

Solves the penalized one-dimensional problem \$\$ \max\_{\lambda \in
\[\lambda\_{\min}, \lambda\_{\max}\]} \frac{r_Q - 3}{2}\log \lambda -
\frac{q}{2}\lambda - \rho \lambda^{-1/2}, \$\$ where \\q\\ is the
current posterior second moment \\m^\top Q m + \mathrm{tr}(Q S)\\.

## Usage

``` r
.optimize_lambda_induced_exp(eq_quad, r_rank, rate, lambda_min, lambda_max)
```

## Arguments

- eq_quad:

  Positive numeric vector of posterior second moments.

- r_rank:

  Positive scalar RW precision rank.

- rate:

  Positive scalar rate parameter \\\rho\\.

- lambda_min, lambda_max:

  Positive bounds with `lambda_min <= lambda_max`.

## Value

Numeric vector of optimized `lambda`.

## Details

The optimization is carried out in \\\eta = \log \lambda\\, for which
the objective is strictly concave and the derivative is monotone
decreasing.
