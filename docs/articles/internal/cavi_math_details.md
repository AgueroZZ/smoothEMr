# Mathematical Details of the MPCurver CAVI Algorithm

``` r
library(MPCurver)
```

This note is written in the style of supplementary methods. Its purpose
is to state the model, define the variational objective, and then derive
the explicit update equations implemented by
[`cavi()`](https://aguerozz.github.io/MPCurver/reference/cavi.md).

Throughout, we focus on the current recommended single-ordering CAVI
model used in `MPCurver`.

We use the following notation throughout:

``` math

C := (c_1,\dots,c_n),
\qquad
U := (u_1,\dots,u_d),
\qquad
\theta := (\pi,\sigma^2,\lambda),
```
where
``` math

\sigma^2 := (\sigma_1^2,\dots,\sigma_d^2),
\qquad
\lambda := (\lambda_1,\dots,\lambda_d).
```

We write $`q`$ for the full variational distribution $`q(C,U)`$,
$`r_{ik}`$ for the variational responsibility of cell $`i`$ for
component $`k`$, and $`r_Q := \operatorname{rank}(Q_K)`$ for the rank of
the random-walk precision.

## 1 Model setup

We observe a matrix
``` math

X = (x_{ij}) \in \mathbb{R}^{n \times d},
```
with cells indexed by $`i = 1,\dots,n`$ and features indexed by
$`j = 1,\dots,d`$.

### 1.1 Latent positions

Each cell is assigned to one of $`K`$ ordered mixture components:
``` math

c_i \in \{1,\dots,K\},
\qquad
\Pr(c_i = k) = \pi_k,
\qquad
\sum_{k=1}^K \pi_k = 1.
```

### 1.2 Feature trajectories

For each feature $`j`$, we introduce a latent trajectory
``` math

u_j = (u_{j1},\dots,u_{jK})^\top \in \mathbb{R}^K,
```
where $`u_{jk}`$ is the mean of feature $`j`$ at component $`k`$.

Conditionally on $`c_i = k`$, the observation model is
``` math

x_{ij} \mid c_i = k, u_j, \sigma_j^2 \sim N(u_{jk}, \sigma_j^2).
```

The current CAVI implementation uses feature-specific variances
$`\sigma_j^2`$, shared across components $`k`$.

### 1.3 GMRF prior

Each trajectory is regularized by an intrinsic Gaussian Markov random
field prior:
``` math

u_j \sim N\!\left(0, (\lambda_j Q_K)^{-}\right),
```
where:

- $`Q_K`$ is the $`K \times K`$ random-walk precision matrix,
- $`\lambda_j > 0`$ is the feature-specific smoothness parameter,
- $`(\cdot)^{-}`$ denotes the generalized inverse because $`Q_K`$ is
  rank-deficient for intrinsic random-walk priors.

Let
``` math

r_Q := \operatorname{rank}(Q_K) = K - q,
```
where $`q`$ is the random-walk order.

Up to additive constants, the prior density contributes
``` math

\log p(u_j \mid \lambda_j)
=
\frac{r_Q}{2}\log \lambda_j
+ \frac{1}{2}\log |Q_K|_+
- \frac{\lambda_j}{2} u_j^\top Q_K u_j,
```
where $`|Q_K|_+`$ is the generalized determinant over the nonzero
eigenvalues.

The current package also supports an optional penalty on $`\lambda_j`$
through the induced prior
``` math

\tau_j := \lambda_j^{-1/2} \sim \operatorname{Exp}(\rho),
\qquad \rho \ge 0.
```
Interpreted as a prior on $`\lambda_j`$, this contributes
``` math

\log p(\lambda_j)
=
\log(\rho/2)
- \frac{3}{2}\log \lambda_j
- \rho \lambda_j^{-1/2},
```
when $`\rho > 0`$. The default setting $`\rho = 0`$ recovers the
original unpenalized CAVI objective.

### 1.4 Joint distribution

The complete-data model factorizes as
``` math

p(X,C,U \mid \theta)
=
\left[\prod_{i=1}^n p(c_i \mid \pi)\right]
\left[\prod_{j=1}^d p(u_j \mid \lambda_j)\right]
\left[\prod_{i=1}^n \prod_{j=1}^d p(x_{ij} \mid c_i, u_j, \sigma_j^2)\right].
```

More explicitly,
``` math

\log p(X,C,U \mid \theta)
=
\sum_{i=1}^n \log \pi_{c_i}
- \frac12 \sum_{i=1}^n \sum_{j=1}^d
\left[
\log(2\pi \sigma_j^2)
+
\frac{(x_{ij} - u_{j,c_i})^2}{\sigma_j^2}
\right]
+ \sum_{j=1}^d
\left[
\frac{r_Q}{2}\log \lambda_j
+ \frac{1}{2}\log |Q_K|_+
- \frac{\lambda_j}{2} u_j^\top Q_K u_j
\right]
+ \text{const}.
```

## 2 General CAVI idea

We approximate the posterior with the mean-field family
``` math

q(C,U)
=
\prod_{i=1}^n q_i(c_i)\prod_{j=1}^d q_j(u_j).
```

Write
``` math

r_{ik} := q_i(c_i = k),
\qquad
\sum_{k=1}^K r_{ik} = 1.
```

For the trajectory factors we use Gaussians:
``` math

q_j(u_j) = N(m_j, S_j),
```
with posterior mean $`m_j \in \mathbb{R}^K`$ and covariance
$`S_j \in \mathbb{R}^{K \times K}`$.

The parameters $`\pi`$, $`\sigma_j^2`$, and $`\lambda_j`$ are treated as
point-estimated blocks rather than random variational factors.

### 2.1 ELBO in expectation form

The variational objective is the evidence lower bound (ELBO),
``` math

\operatorname{ELBO}(q;\theta)
=
\mathbb E_q \big[\log p(X,C,U \mid \theta)\big]
- \mathbb E_q \big[\log q(C,U)\big].
```

For fixed $`\theta`$, this is a lower bound on the marginal log-evidence
``` math

\log p(X \mid \theta)
=
\log \int \sum_C p(X,C,U \mid \theta)\, dU.
```
More precisely,
``` math

\log p(X \mid \theta)
=
\operatorname{ELBO}(q;\theta)
+
\mathrm{KL}\!\left(
q(C,U)\,\|\,p(C,U \mid X,\theta)
\right).
```
So maximizing the ELBO over $`q`$ for fixed hyperparameters is
equivalent to driving the variational approximation toward the exact
posterior.

Once $`q`$ has been updated, we then maximize the same ELBO with respect
to the point-estimated block $`\theta`$. This is a standard variational
EM / variational empirical-Bayes viewpoint:

- first optimize the variational approximation $`q`$,
- then optimize the ELBO as a surrogate for $`\log p(X \mid \theta)`$,
- hoping that when the KL gap is small, the optimized ELBO is close to
  the true marginal log-likelihood.

Equivalently, if we define the profiled variational objective
``` math

F(\theta)
:=
\max_q \operatorname{ELBO}(q;\theta),
```
then
``` math

F(\theta)
\le
\log p(X \mid \theta).
```
So at a high level the algorithm alternates between:

1.  updating the variational approximation $`q`$ to improve the bound,
2.  then updating $`\theta`$ to increase the optimized bound $`F`$.

This is why the procedure can be viewed as “approximately maximizing the
marginal log-likelihood”: we optimize a lower bound to it, rather than
the intractable quantity directly.

Substituting the model terms into the ELBO definition gives
``` math

\operatorname{ELBO}(q;\theta)
=
\mathbb E_q[\log p(C \mid \pi)]
+ \mathbb E_q[\log p(X \mid C,U,\sigma^2)]
+ \mathbb E_q[\log p(U \mid \lambda)]
+ H(q_C)
+ H(q_U),
```
where $`H(\cdot)`$ denotes entropy.

The major coordinate-ascent idea is:

1.  fix all blocks except one,
2.  maximize the ELBO with respect to that block,
3.  cycle until convergence.

For the variational factors, the generic mean-field update is
``` math

\log q_\ell^\star(\xi_\ell)
=
\mathbb E_{q_{-\ell}}[\log p(X,C,U \mid \theta)]
+ \text{const},
```
where $`\xi_\ell`$ denotes the argument of the current variational
block. In the present model this means:

- $`\xi_\ell = c_i`$ when the current block is the categorical factor
  $`q_i(c_i)`$,
- $`\xi_\ell = u_j`$ when the current block is the Gaussian factor
  $`q_j(u_j)`$.

Here $`q_{-\ell}`$ means “all variational factors except the current
one”.

For the point-estimated parameter block $`\theta`$, we maximize the same
ELBO directly with respect to that parameter block.

## 3 Expanded ELBO

To derive the code formulas, define the effective component sizes
``` math

N_k := \sum_{i=1}^n r_{ik}.
```

Then the mixture-weight term is
``` math

\mathbb E_q[\log p(C \mid \pi)]
=
\sum_{i=1}^n \sum_{k=1}^K r_{ik}\log \pi_k.
```

The entropy of the categorical factors is
``` math

H(q_C)
=
- \sum_{i=1}^n \sum_{k=1}^K r_{ik}\log r_{ik}.
```

For the Gaussian likelihood,
``` math

\mathbb E_q[\log p(X \mid C,U,\sigma^2)]
=
-\frac12 \sum_{i=1}^n\sum_{k=1}^K\sum_{j=1}^d
r_{ik}
\left[
\log(2\pi \sigma_j^2)
+
\frac{\mathbb E_q[(x_{ij} - u_{jk})^2]}{\sigma_j^2}
\right].
```

Since $`q_j(u_j)=N(m_j,S_j)`$,
``` math

\mathbb E_q[(x_{ij} - u_{jk})^2]
=
(x_{ij} - m_{jk})^2 + S_{j,kk}.
```

So the likelihood contribution becomes
``` math

-\frac12 \sum_{i=1}^n\sum_{k=1}^K\sum_{j=1}^d
r_{ik}
\left[
\log(2\pi \sigma_j^2)
+
\frac{(x_{ij} - m_{jk})^2 + S_{j,kk}}{\sigma_j^2}
\right].
```

For the prior,
``` math

\mathbb E_q[\log p(U \mid \lambda)]
=
\frac12 \sum_{j=1}^d
\left[
r_Q \log \lambda_j
+ \log |Q_K|_+
- \lambda_j \mathbb E_q(u_j^\top Q_K u_j)
\right].
```

Using Gaussian second moments,
``` math

\mathbb E_q(u_j^\top Q_K u_j)
=
m_j^\top Q_K m_j + \operatorname{tr}(Q_K S_j).
```

Finally, the entropy of $`q_j(u_j)=N(m_j,S_j)`$ is
``` math

H(q_j)
=
\frac12 \left[K(1+\log 2\pi) + \log |S_j|\right],
```
so
``` math

H(q_U)
=
\frac12 \sum_{j=1}^d \left[K(1+\log 2\pi) + \log |S_j|\right].
```

Putting everything together,
``` math

\operatorname{ELBO}(q;\theta)
=
\sum_{i=1}^n \sum_{k=1}^K r_{ik}\log \pi_k
- \sum_{i=1}^n \sum_{k=1}^K r_{ik}\log r_{ik}
```
``` math

-\frac12 \sum_{i=1}^n\sum_{k=1}^K\sum_{j=1}^d
r_{ik}
\left[
\log(2\pi \sigma_j^2)
+
\frac{(x_{ij} - m_{jk})^2 + S_{j,kk}}{\sigma_j^2}
\right]
```
``` math

+ \frac12 \sum_{j=1}^d
\left[
r_Q \log \lambda_j
+ \log |Q_K|_+
- \lambda_j\left(m_j^\top Q_K m_j + \operatorname{tr}(Q_K S_j)\right)
+ \log |S_j|
+ K(1+\log 2\pi)
\right].
```

This expanded form is the starting point for every concrete update in
the code.

## 4 Specific coordinate updates

## 4.1 Update for $`q_j(u_j)`$

Fix $`r_{ik}`$, $`\sigma_j^2`$, and $`\lambda_j`$, and keep only the
ELBO terms that depend on $`u_j`$.

Ignoring constants, we get
``` math

\log q_j^\star(u_j)
=
-\frac{1}{2\sigma_j^2}
\sum_{i=1}^n \sum_{k=1}^K r_{ik}(x_{ij} - u_{jk})^2
- \frac{\lambda_j}{2} u_j^\top Q_K u_j
+ \text{const}.
```

Expanding the quadratic term in $`u_j`$,
``` math

\log q_j^\star(u_j)
=
-\frac12 u_j^\top
\left[
\lambda_j Q_K + \operatorname{diag}(N_1/\sigma_j^2,\dots,N_K/\sigma_j^2)
\right]
u_j
+ u_j^\top b_j
+ \text{const},
```
where
``` math

b_{jk}
=
\frac{1}{\sigma_j^2}\sum_{i=1}^n r_{ik} x_{ij}.
```

Therefore $`q_j^\star(u_j)`$ is Gaussian with
``` math

A_j
:=
\lambda_j Q_K + \operatorname{diag}(N_1/\sigma_j^2,\dots,N_K/\sigma_j^2),
```
``` math

S_j = A_j^{-1},
\qquad
m_j = S_j b_j.
```

This is exactly the matrix solve used in the implementation.

## 4.2 Update for $`q_i(c_i)`$

Fix $`q(U)`$, $`\pi`$, and $`\sigma_j^2`$, and focus on the terms
involving the single latent variable $`c_i`$.

The general mean-field rule gives
``` math

\log q_i^\star(c_i = k)
\propto
\log \pi_k
+ \mathbb E_q[\log p(x_i \mid c_i = k, U, \sigma^2)].
```

Using the Gaussian likelihood,
``` math

\log q_i^\star(c_i = k)
\propto
\log \pi_k
- \frac12 \sum_{j=1}^d
\left[
\log(2\pi \sigma_j^2)
+
\frac{(x_{ij} - m_{jk})^2 + S_{j,kk}}{\sigma_j^2}
\right].
```

So the update is a row-wise softmax:
``` math

r_{ik}
=
\frac{\exp(\ell_{ik})}{\sum_{k'=1}^K \exp(\ell_{ik'})},
```
with
``` math

\ell_{ik}
:=
\log \pi_k
- \frac12 \sum_{j=1}^d
\left[
\log(2\pi \sigma_j^2)
+
\frac{(x_{ij} - m_{jk})^2 + S_{j,kk}}{\sigma_j^2}
\right].
```

The posterior variance term $`S_{j,kk}`$ is the uncertainty correction
that distinguishes this update from a pure plug-in EM step.

## 4.3 Update for $`\pi_k`$

Holding everything else fixed, maximize
``` math

\sum_{i=1}^n \sum_{k=1}^K r_{ik}\log \pi_k
```
subject to $`\sum_k \pi_k = 1`$.

Using a Lagrange multiplier,
``` math

\pi_k = \frac{N_k}{n}.
```

## 4.4 Update for $`\sigma_j^2`$

Fix $`q(U)`$ and $`q(C)`$, and keep only the ELBO terms involving
$`\sigma_j^2`$:
``` math

-\frac12 \sum_{i=1}^n\sum_{k=1}^K r_{ik}
\left[
\log \sigma_j^2
+
\frac{(x_{ij} - m_{jk})^2 + S_{j,kk}}{\sigma_j^2}
\right].
```

Differentiating with respect to $`\sigma_j^2`$ and setting the
derivative to zero yields
``` math

\sigma_j^2
=
\frac{1}{n}
\sum_{i=1}^n \sum_{k=1}^K
r_{ik}\left[(x_{ij} - m_{jk})^2 + S_{j,kk}\right].
```

In the code this is followed by projection onto the allowed range
`[sigma_min, sigma_max]`.

## 4.5 Update for $`\lambda_j`$

Fix $`q(U)`$, and write
``` math

q_j
:=
m_j^\top Q_K m_j + \operatorname{tr}(Q_K S_j).
```

Without the optional $`\tau_j \sim \operatorname{Exp}(\rho)`$ penalty,
the only ELBO terms involving $`\lambda_j`$ are
``` math

\frac{r_Q}{2}\log \lambda_j - \frac{\lambda_j}{2} q_j.
```

Differentiating and setting to zero gives the familiar closed-form
update
``` math

\lambda_j = \frac{r_Q}{q_j},
```
followed by projection onto `[lambda_min, lambda_max]`.

With the induced prior active, the $`\lambda_j`$ block becomes
``` math

\frac{r_Q-3}{2}\log \lambda_j
- \frac{\lambda_j}{2} q_j
- \rho \lambda_j^{-1/2}
+ \log(\rho/2).
```

This no longer has a closed-form maximizer in $`\lambda_j`$. The
implementation therefore optimizes the same scalar objective in
``` math

\eta_j := \log \lambda_j,
```
which yields
``` math

\frac{r_Q-3}{2}\eta_j
- \frac{q_j}{2} e^{\eta_j}
- \rho e^{-\eta_j/2}
+ \log(\rho/2).
```

Its derivative is
``` math

\frac{r_Q-3}{2}
- \frac{q_j}{2} e^{\eta_j}
+ \frac{\rho}{2} e^{-\eta_j/2},
```
which is strictly decreasing in $`\eta_j`$, so the optimum is unique. In
the code this is solved by a bounded one-dimensional root search on
`log(lambda_j)`, and the resulting maximizer is then mapped back to
`[lambda_min, lambda_max]`.

## 5 What the code actually stores

The main quantities in the derivation correspond to the following
objects in
[`cavi()`](https://aguerozz.github.io/MPCurver/reference/cavi.md):

| Mathematical quantity | Code object |
|----|----|
| $`r_{ik}`$ | `R[i, k]` or final `fit$gamma[i, k]` |
| $`N_k = \sum_i r_{ik}`$ | `Nk <- colSums(R)` |
| $`m_{jk}`$ | `q_u$m_mat[j, k]` or `fit$posterior$mean[j, k]` |
| $`S_j`$ | `q_u$S_list[[j]]` or `fit$posterior$cov[[j]]` |
| $`S_{j,kk}`$ | `q_u$sdiag_mat[j, k]` or `fit$posterior$var[j, k]` |
| $`\sum_i r_{ik} x_{ij}`$ | `GX[j, k]` |
| $`m_j^\top Q_K m_j + \operatorname{tr}(Q_K S_j)`$ | `q_u$eq_quad[j]` |

So the derivation above is not just conceptual: it is exactly the
algebra that drives the implementation.

## 6 Why this gives a clean monotonicity story

The key point is that all five updates above are exact optimizers of
their own ELBO blocks:

- Gaussian update for $`q_j(u_j)`$,
- softmax update for $`q_i(c_i)`$,
- normalized counts for $`\pi_k`$,
- closed-form posterior residual update for $`\sigma_j^2`$,
- closed-form or one-dimensional concave update for $`\lambda_j`$,
  depending on whether the induced $`\lambda`$ prior is active.

Because
[`cavi()`](https://aguerozz.github.io/MPCurver/reference/cavi.md) cycles
through exact coordinate maximizers of the same explicit ELBO, the
recorded `elbo_trace` should be nondecreasing up to numerical tolerance.
Importantly, when `lambda_sd_prior_rate > 0`, the stored `elbo_trace`
includes the full induced-prior contribution as well. That is the core
optimization guarantee that makes the current `cavi` pipeline cleaner
than the older hybrid collapsed-EM path.
