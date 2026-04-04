# Conceptual and Mathematical Comparison of CAVI and csmooth_em

``` r
library(MPCurver)
```

This note compares the two main single-ordering fitting strategies
currently available in `MPCurver`:

- `cavi`, which is the recommended default,
- and the legacy `csmooth_em` pipeline.

The emphasis is conceptual and mathematical rather than empirical. The
goal is to place the two procedures under a common notation and then
identify which parts of the inferential target are shared and which are
genuinely different.

We use the same notation as in the CAVI theory note:
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

## 1 Shared statistical skeleton

At the model level, `cavi` and `csmooth_em` are not solving unrelated
problems. They start from the same basic MPCurve ingredients.

We observe
``` math

X = (x_{ij}) \in \mathbb{R}^{n \times d},
```
with latent cell positions
``` math

c_i \in \{1,\dots,K\},
\qquad
\Pr(c_i = k) = \pi_k,
```
and feature-specific trajectories
``` math

u_j = (u_{j1},\dots,u_{jK})^\top \in \mathbb{R}^K.
```

The observation model is
``` math

x_{ij} \mid c_i = k, u_j, \sigma_j^2 \sim N(u_{jk}, \sigma_j^2),
```
and the smoothness prior is
``` math

u_j \sim N\!\left(0, (\lambda_j Q_K)^{-}\right),
```
where $`Q_K`$ is the random-walk precision matrix along the ordered
component axis.

So at the level of likelihood and prior, the two pipelines are closely
related. The main differences arise in how inference is carried out.

## 2 The key conceptual split

The principal inferential distinction is:

- `csmooth_em` treats the trajectories $`U`$ as deterministic quantities
  to be estimated, typically by MAP / penalized optimization;
- `cavi` treats the trajectories as latent random quantities and carries
  an explicit approximate posterior for them.

More concretely:

### `csmooth_em`

- keeps a soft posterior over cell positions $`C`$ through the
  responsibility matrix $`\Gamma`$,
- but uses point estimates for the trajectories,
- and in the adaptive-ML path augments this with a collapsed /
  Laplace-style objective for $`\lambda_j`$ and $`\sigma_j^2`$.

### `cavi`

- keeps a soft posterior over cell positions $`C`$,
- keeps a Gaussian variational posterior over each trajectory $`u_j`$,
- and updates $`\pi_k`$, $`\sigma_j^2`$, and $`\lambda_j`$ by maximizing
  the same explicit ELBO used for the variational factors.

That is the main philosophical distinction: plug-in MAP versus explicit
mean-field uncertainty quantification.

## 3 What objective does each algorithm optimize?

## 3.1 `cavi`

For `cavi`, the central object is the variational lower bound
``` math

\operatorname{ELBO}(q;\theta)
=
\mathbb E_q[\log p(X,C,U \mid \theta)]
- \mathbb E_q[\log q(C,U)].
```

For fixed $`\theta`$, this is a lower bound on the marginal log-evidence
``` math

\log p(X \mid \theta).
```
Indeed,
``` math

\log p(X \mid \theta)
=
\operatorname{ELBO}(q;\theta)
+
\mathrm{KL}\!\left(
q(C,U)\,\|\,p(C,U \mid X,\theta)
\right).
```

So `cavi` can be understood as variational empirical Bayes:

1.  update the variational approximation $`q`$,
2.  then optimize the same bound with respect to $`\theta`$.

If we define the profiled variational objective
``` math

F(\theta)
:=
\max_q \operatorname{ELBO}(q;\theta),
```
then the algorithm is effectively trying to increase $`F`$, which is
still a lower bound on the true marginal log-likelihood.

## 3.2 `csmooth_em`

For `csmooth_em`, there are really two regimes.

### Fixed-$`\lambda`$ regime

When $`\lambda`$ is fixed, `csmooth_em` behaves like a penalized EM /
MAP algorithm. A representative objective is
``` math

J_{\mathrm{MAP}}(\pi,U,\sigma^2;\lambda)
=
\log p(X \mid \pi,U,\sigma^2)
- \frac12 \sum_{j=1}^d \lambda_j u_j^\top Q_K u_j
+ \text{const},
```
or equivalently an EM lower bound where $`U`$ is optimized as a point
estimate.

In this regime, the algorithm is close to standard penalized EM: update
responsibilities, then update $`(\pi,U,\sigma^2)`$.

### Adaptive-ML regime

When `adaptive = "ml"`, the code follows a different logic. It no longer
treats the plain ELBO as the main optimization certificate. Instead it
tracks a collapsed objective of the form
``` math

J_{\mathrm{coll}} = \mathrm{penELBO} + \mathrm{const} - \frac12 \log |H|,
```
where $`H`$ is the local Hessian around the MAP solution for the
trajectory means.

Conceptually, this is trying to move closer to a marginal-likelihood
objective by integrating out or locally collapsing the trajectory block
$`U`$, but it is not the same object as the explicit `cavi` ELBO.

So the two algorithms are genuinely optimizing different things:

- `cavi`: an explicit variational lower bound for a model with latent
  $`C`$ and $`U`$,
- `csmooth_em`: a plug-in MAP/EM objective, or in adaptive-ML mode a
  collapsed / Laplace-style objective for the point-estimated trajectory
  model.

## 4 Side-by-side update structure

The following table is a useful high-level summary.

| Block | `cavi` | `csmooth_em` |
|----|----|----|
| Cell positions | variational posterior $`q(C)`$ | EM responsibilities $`\Gamma`$ |
| Trajectories | variational posterior $`q(U)`$ | deterministic MAP estimates $`U`$ |
| Mixture weights | ELBO maximization | EM / MAP update |
| Variances | ELBO maximization using posterior moments | plug-in or collapsed hyper-step |
| Smoothness $`\lambda_j`$ | ELBO maximization using posterior moments | profile / collapsed hyper-step |
| Main recorded objective | `elbo_trace` | fixed-$`\lambda`$: ELBO-like; adaptive-ML: `ml_trace` |

## 5 Specific mathematical differences

## 5.1 Treatment of $`U`$

This is the most important difference.

### In `cavi`

Each feature trajectory has a Gaussian posterior
``` math

q_j(u_j) = N(m_j, S_j).
```

So all later updates depend on both:

- posterior means $`m_j`$,
- posterior variances / covariances $`S_j`$.

### In `csmooth_em`

The trajectory is reduced to a point estimate
``` math

\hat u_j = \arg\max_{u_j} \Big\{\text{penalized complete-data objective}\Big\}.
```

So downstream updates are based on $`\hat u_j`$ only, not on an
uncertainty matrix $`S_j`$.

This is why `cavi` has an uncertainty correction in its responsibility
update, while `csmooth_em` does not.

## 5.2 E-step / responsibility update

### `csmooth_em`

A plug-in responsibility update has the form
``` math

\gamma_{ik}
\propto
\pi_k
\exp\!\left[
-\frac12 \sum_{j=1}^d
\left(
\log(2\pi \sigma_j^2)
+
\frac{(x_{ij} - \hat u_{jk})^2}{\sigma_j^2}
\right)
\right].
```

### `cavi`

The variational update for $`q(c_i)`$ is
``` math

r_{ik}
\propto
\pi_k
\exp\!\left[
-\frac12 \sum_{j=1}^d
\left(
\log(2\pi \sigma_j^2)
+
\frac{(x_{ij} - m_{jk})^2 + S_{j,kk}}{\sigma_j^2}
\right)
\right].
```

The extra term $`S_{j,kk}`$ is the posterior uncertainty correction.

This immediately explains one visible empirical difference: early in
optimization, `cavi` can move more aggressively because it is not
treating the current trajectory means as certainty.

## 5.3 Update for the trajectories

### `cavi`

The update is a full Gaussian variational block:
``` math

q_j(u_j) = N(m_j,S_j),
```
with
``` math

S_j^{-1}
=
\lambda_j Q_K + \operatorname{diag}(N_1/\sigma_j^2,\dots,N_K/\sigma_j^2),
\qquad
m_j = S_j b_j.
```

### `csmooth_em`

The update is a penalized least-squares / MAP step:
``` math

\hat u_j
=
\arg\min_{u_j}
\left\{
\frac{1}{2\sigma_j^2}\sum_{i=1}^n\sum_{k=1}^K \gamma_{ik}(x_{ij}-u_{jk})^2
+ \frac{\lambda_j}{2} u_j^\top Q_K u_j
\right\}.
```

So one can think of the `csmooth_em` update as the mode of the
corresponding Gaussian block, whereas `cavi` keeps both the mode and the
covariance.

## 5.4 Hyperparameter updates

### `cavi`

The updates for $`\sigma_j^2`$ and $`\lambda_j`$ are derived from the
same ELBO and use posterior second moments:
``` math

\sigma_j^2
=
\frac{1}{n}\sum_{i=1}^n\sum_{k=1}^K
r_{ik}\Big[(x_{ij}-m_{jk})^2 + S_{j,kk}\Big],
```
``` math

\lambda_j
=
\frac{\operatorname{rank}(Q_K)}
{m_j^\top Q_K m_j + \operatorname{tr}(Q_K S_j)}.
```

### `csmooth_em`

In the legacy adaptive-ML path, the hyperparameter updates are instead
tied to the collapsed objective $`C`$. So even when the two algorithms
share the same likelihood and prior, their $`\lambda_j`$ and
$`\sigma_j^2`$ updates are conceptually different.

This is one reason why `cavi` and `csmooth_em` may agree closely in some
problems but still not trace the same optimization path.

## 6 A useful limiting interpretation

There is a practical sense in which `csmooth_em` can be viewed as a
limiting or plug-in version of `cavi`.

If the posterior $`q_j(u_j)`$ in `cavi` becomes very concentrated, then:

- $`S_j`$ becomes small,
- $`m_j`$ behaves like a point estimate,
- the `cavi` responsibility update approaches the plug-in `csmooth_em`
  responsibility update.

So in a regime where posterior uncertainty around the trajectories is
tiny, the two algorithms may produce very similar final fits.

But this should not hide the conceptual difference:

- `cavi` is still optimizing a variational lower bound for a
  latent-$`U`$ model,
- `csmooth_em` is still optimizing a point-estimation objective for
  $`U`$, with an optional collapsed hyperparameter correction.

## 7 Why the convergence stories differ

## 7.1 `cavi`

`cavi` has the cleaner monotonicity story because each block update is
an exact coordinate maximizer of the same explicit objective:
``` math

\operatorname{ELBO}(q;\theta).
```

So the recorded `elbo_trace` should be nondecreasing up to numerical
tolerance.

## 7.2 `csmooth_em`

The convergence story depends on the regime:

- with fixed $`\lambda`$, the usual EM / penalized-EM intuition applies;
- with adaptive `ml`, the main quantity to watch is the collapsed
  objective `ml_trace`, not the plain ELBO.

So the legacy method can still have a coherent optimization story, but
it is a different story from `cavi`.

## 8 Bottom line

The most honest high-level summary is:

- `cavi` and `csmooth_em` share the same model skeleton,
- but they do not share the same inferential target,
- and they do not share the same primary optimization objective.

`cavi` is the more principled choice if you want:

- an explicit approximate posterior over trajectories,
- a single clean ELBO to optimize,
- and a transparent monotonicity guarantee.

`csmooth_em` remains useful as:

- a strong legacy MAP/EM baseline,
- a benchmark path,
- and a collapsed-ML alternative for internal comparison.
