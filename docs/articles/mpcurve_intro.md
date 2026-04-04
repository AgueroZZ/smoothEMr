# Getting Started with MPCurver CAVI

``` r
library(MPCurver)
set.seed(42)
```

`MPCurver` now uses `cavi` as the recommended default fitting backend.
This vignette introduces the statistical model, the variational
algorithm, and then walks through a larger but still easy spiral
example.

The legacy `csmooth_em` pipeline remains available only through
lower-level compatibility functions and is no longer part of the public
wrapper workflow.

## 1 Statistical model

We observe a data matrix
``` math

X = (x_{ij}) \in \mathbb{R}^{n \times d},
```
where row $`i`$ is one cell and column $`j`$ is one feature.

The single-ordering MPCurve model assumes:

1.  each cell has a latent position
    ``` math

    c_i \in \{1,\dots,K\},
    ```
2.  each feature has a smooth trajectory over the $`K`$ ordered
    components, written as
    ``` math

    u_j = (u_{j1},\dots,u_{jK})^\top \in \mathbb{R}^K,
    ```
3.  observations follow a Gaussian likelihood
    ``` math

    x_{ij} \mid c_i = k, u_j, \sigma_j^2 \sim N(u_{jk}, \sigma_j^2),
    ```
4.  and each trajectory is regularized by a Gaussian Markov random field
    prior along the component axis:
    ``` math

    u_j \sim N\!\left(0, (\lambda_j Q_K)^{-1}\right),
    ```
    where $`Q_K`$ is the random-walk precision matrix on the $`K`$ bins.

In the current `cavi` implementation:

- `sigma_j^2` is feature-specific,
- `lambda_j` is feature-specific,
- the covariance model is homoskedastic across components,
- and the posterior is approximated with a mean-field family
  ``` math

  q(C,U) = \prod_{i=1}^n q_i(c_i)\prod_{j=1}^d q_j(u_j).
  ```

## 2 Variational algorithm

The [`cavi()`](https://aguerozz.github.io/MPCurver/reference/cavi.md)
routine performs coordinate-ascent variational inference on an explicit
ELBO. Write
``` math

r_{ik} := q_i(c_i = k),
\qquad
N_k := \sum_{i=1}^n r_{ik},
```
and for each feature $`j`$,
``` math

q_j(u_j) = N(m_j, S_j).
```

Under the mean-field family
``` math

q(C,U)=\prod_{i=1}^n q_i(c_i)\prod_{j=1}^d q_j(u_j),
```
the ELBO can be written, up to additive constants, as
``` math

\mathcal L(q)
 =
\sum_{i=1}^n\sum_{k=1}^K r_{ik}\log \pi_k
- \sum_{i=1}^n\sum_{k=1}^K r_{ik}\log r_{ik}
```
``` math

- \frac12 \sum_{i=1}^n\sum_{k=1}^K\sum_{j=1}^d
r_{ik}
\left[
\log(2\pi \sigma_j^2)
+
\frac{(x_{ij}-m_{jk})^2 + S_{j,kk}}{\sigma_j^2}
\right]
```
``` math

+ \frac12 \sum_{j=1}^d
\left[
\operatorname{rank}(Q_K)\log \lambda_j
- \lambda_j\Big(m_j^\top Q_K m_j + \operatorname{tr}(Q_K S_j)\Big)
+ \log |S_j|
\right].
```

This is the quantity recorded in `elbo_trace`. The updates used in
[`cavi()`](https://aguerozz.github.io/MPCurver/reference/cavi.md) are
obtained by maximizing this ELBO one block at a time.

### 2.1 Update for $`q(u_j)`$

Fix $`r_{ik}`$, $`\sigma_j^2`$, and $`\lambda_j`$. The ELBO terms
involving $`u_j`$ are quadratic, so the optimal variational factor is
Gaussian:
``` math

\log q_j^\star(u_j)
=
-\frac12 u_j^\top A_j u_j + u_j^\top b_j + \text{const},
```
where
``` math

A_j
=
\lambda_j Q_K +
\operatorname{diag}\!\left(N_1/\sigma_j^2,\dots,N_K/\sigma_j^2\right),
```
and
``` math

b_{jk}
=
\frac{1}{\sigma_j^2}\sum_{i=1}^n r_{ik} x_{ij}.
```
Therefore,
``` math

S_j = A_j^{-1},
\qquad
m_j = S_j b_j.
```

This is exactly the linear-algebra step implemented in
[`cavi()`](https://aguerozz.github.io/MPCurver/reference/cavi.md).

### 2.2 Update for $`q(c_i)`$

Fix $`q(U)`$, $`\pi`$, and $`\sigma_j^2`$. Collecting the ELBO terms
involving the row $`r_{i\cdot}`$ gives
``` math

\log r_{ik}
\propto
\log \pi_k
- \frac12 \sum_{j=1}^d
\left[
\log(2\pi \sigma_j^2)
+
\frac{(x_{ij}-m_{jk})^2 + S_{j,kk}}{\sigma_j^2}
\right].
```
So the cell-position posterior is updated by a row-wise softmax:
``` math

r_{i\cdot} \leftarrow \operatorname{softmax}\!\big(\ell_{i1},\dots,\ell_{iK}\big),
```
with $`\ell_{ik}`$ equal to the right-hand side above.

The extra $`S_{j,kk}`$ term is the uncertainty correction that
distinguishes this update from a pure plug-in EM step.

### 2.3 Update for $`\pi_k`$

Fixing everything else, we maximize
``` math

\sum_{k=1}^K N_k \log \pi_k
```
subject to $`\sum_k \pi_k = 1`$. A Lagrange multiplier gives
``` math

\pi_k = \frac{N_k}{n}.
```

### 2.4 Update for $`\sigma_j^2`$

Fixing $`q(U)`$ and $`q(C)`$, differentiate the ELBO with respect to
$`\sigma_j^2`$. The stationary point is
``` math

\sigma_j^2
=
\frac{1}{n}
\sum_{i=1}^n \sum_{k=1}^K
r_{ik}\Big[(x_{ij}-m_{jk})^2 + S_{j,kk}\Big].
```

So each feature variance is updated from the posterior expected residual
sum of squares.

### 2.5 Update for $`\lambda_j`$

Fixing $`q(U)`$, only the prior part of the ELBO depends on
$`\lambda_j`$:
``` math

\frac12 \operatorname{rank}(Q_K)\log \lambda_j
- \frac12 \lambda_j
\Big(m_j^\top Q_K m_j + \operatorname{tr}(Q_K S_j)\Big).
```
Setting the derivative to zero yields
``` math

\lambda_j
=
\frac{\operatorname{rank}(Q_K)}
{m_j^\top Q_K m_j + \operatorname{tr}(Q_K S_j)}.
```

### 2.6 Why the ELBO is monotone here

Each of the five updates above is the exact optimizer of one ELBO block
while holding the others fixed. That is why the current `cavi`
implementation has a clean monotonicity story: the recorded `elbo_trace`
should be nondecreasing up to numerical tolerance.

## 3 Simulated spiral example

We use a larger and more stable spiral than the older vignette examples:

- `n = 1000` cells,
- a shorter latent range `t in [1.5*pi, 4.5*pi]`,
- lower noise,
- and a finer mixture grid `K = 60`.

``` r
sim <- simulate_swiss_roll_1d_2d(
  n = 1000,
  t_range = c(1.5 * pi, 4.5 * pi),
  sigma = 0.12,
  seed = 1
)

X <- as.matrix(sim$obs)
t_true <- sim$t

cat("Data dimensions:", nrow(X), "x", ncol(X), "\n")
#> Data dimensions: 1000 x 2
cat("t range:", round(range(t_true), 3), "\n")
#> t range: 4.725 14.137
```

``` r
pal <- grDevices::colorRampPalette(
  c("#1F3A5F", "#2E8A99", "#F2C14E", "#E76F51")
)(256)
col_true <- pal[cut(t_true, breaks = 256, labels = FALSE)]

plot(X, col = col_true, pch = 19, cex = 0.55,
     xlab = "x", ylab = "y",
     main = "Observed spiral, coloured by true pseudotime")
```

![Shorter spiral used in the introductory CAVI
example.](mpcurve_intro_files/figure-html/plot-data-1.png)

Shorter spiral used in the introductory CAVI example.

## 4 Fit the default `cavi` backend

``` r
fit <- fit_mpcurve(
  X,
  method = "fiedler",
  K = 60,
  iter = 120
)

stopifnot(inherits(fit, "mpcurve"))
stopifnot(identical(fit$algorithm, "cavi"))

fit
#> MPCurve model fit
#>   Algorithm : cavi
#>   Model     : homoskedastic
#>   n / d / K : 1000 / 2 / 60
#>   Iterations: 46
#>   ELBO (last)  : -4787.638294
#>   logLik (last): -4474.988366
```

``` r
summary(fit)
#> MPCurve Model Summary
#> Algorithm : cavi  |  Model : homoskedastic  |  n=1000  d=2  K=60
#> Underlying fit summary
#> <summary.cavi>
#>   n = 1000, d = 2, K = 60
#>   model = homoskedastic, RW(q) = 2
#>   init = fiedler, discretization = quantile, adaptive = variational
#>   iter = 46, converged = yes
#>   last ELBO = -4787.6383 (delta last = 0.004311)
#>   last logLik = -4474.9884 (delta last = 0.004553)
#>   pi range = [5.265e-08, 0.02495]
#>   sigma2 range = [0.1072, 0.1197]
#>   lambda range = [0.7521, 2.967]
#>   relative change (last step):
#>     ELBO: 9e-07
#>     logLik: 1.02e-06
#>     lambda_vec (L2 rel): 4.74e-05
```

## 5 Check convergence and recovery

The default `cavi` fit should have a nondecreasing ELBO trace.

``` r
stopifnot(all(diff(fit$elbo_trace) >= -1e-8))
cat("ELBO monotone:", TRUE, "\n")
#> ELBO monotone: TRUE
cat("ELBO range:", round(range(fit$elbo_trace), 3), "\n")
#> ELBO range: -7595.206 -4787.638
```

We can also compare the inferred ordering to the true latent parameter
using the posterior-mean component index
``` math

\hat t_i = \sum_{k=1}^K k \, r_{ik}.
```

``` r
pseudotime_hat <- as.numeric(fit$gamma %*% seq_len(ncol(fit$gamma)))
rho_t <- abs(cor(pseudotime_hat, t_true, method = "spearman"))

cat("Absolute Spearman correlation with true t:", round(rho_t, 4), "\n")
#> Absolute Spearman correlation with true t: 0.9998
```

## 6 Visualize the fitted ordering

``` r
col_fit <- pal[cut(pseudotime_hat, breaks = 256, labels = FALSE)]

op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
plot(X, col = col_true, pch = 19, cex = 0.5,
     xlab = "x", ylab = "y", main = "Truth")
plot(X, col = col_fit, pch = 19, cex = 0.5,
     xlab = "x", ylab = "y", main = "CAVI fit")
```

![True pseudotime (left) and inferred CAVI pseudotime
(right).](mpcurve_intro_files/figure-html/plot-fit-1.png)

True pseudotime (left) and inferred CAVI pseudotime (right).

``` r
par(op)
```

The built-in `mpcurve` plotting methods still work on top of the default
`cavi` backend.

``` r
plot(fit, plot_type = "elbo")
```

![CAVI trace plot from the fitted mpcurve
object.](mpcurve_intro_files/figure-html/plot-trace-1.png)

CAVI trace plot from the fitted mpcurve object.

``` r
plot(fit, plot_type = "mu")
```

![Posterior mean trajectories across the K
bins.](mpcurve_intro_files/figure-html/plot-mu-1.png)

Posterior mean trajectories across the K bins.

## 7 Takeaway

For this easier spiral regime, `cavi` is:

- fast enough to use as the default frontend algorithm,
- explicit about its optimization target through the ELBO,
- and able to recover the ordering very accurately on a moderately large
  example with a fairly fine grid.
