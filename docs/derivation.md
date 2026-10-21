# Derivation of the possible-range prior and its Bayes factor

This document derives the two formulas that `R/possible_range.R` uses:

1. the prior density on the regression coefficients (used for estimation in brms), and
2. the closed-form Bayes factor (used by `bayes_factor_possible_range()`).

Each section ends with the function that implements it.

---

## 1. Model and notation

Linear regression with $n$ observations and $k$ predictors:

$$
y = \alpha \mathbf{1} + X_c \beta + \varepsilon, \qquad \varepsilon \sim N(0, \sigma^2 I_n),
$$

where $X_c$ holds the **centered** predictor columns (each column has mean zero).
Centering makes the intercept $\alpha$ the mean of $y$ and separates it from the slopes.

Notation:

| Symbol | Meaning |
|---|---|
| $\Sigma = X_c^\top X_c / n$ | sample covariance matrix of the predictors |
| $\Sigma = L L^\top$ | Cholesky factorization ($L$ lower triangular) |
| $\mathrm{SST} = \sum_i (y_i - \bar y)^2$ | total sum of squares |
| $R^2$ | sample explained share of variance of the full regression |
| $\rho^2$ | **population** explained share (a parameter, defined below) |

The prior is defined conditional on the observed predictors, as in Zellner's
g-prior and the JZS default: the predictors are treated as fixed.

---

## 2. The prior in three statements

**(a) Neutral priors on intercept and noise level** (the same as the JZS default):

$$
p(\alpha) \propto 1, \qquad p(\sigma) \propto 1/\sigma .
$$

These are improper, but they are shared by every model being compared, so
their arbitrary constants cancel in Bayes factors.

**(b) The explained share gets the same prior for any number of predictors.**
The signal variance is $\beta^\top \Sigma \beta$, so the population explained share is

$$
\rho^2 = \frac{\beta^\top \Sigma \beta}{\beta^\top \Sigma \beta + \sigma^2},
\qquad \rho^2 \sim \mathrm{Beta}(\tfrac12, 1).
$$

**(c) Every direction is equally plausible.** Write the coefficients in
"whitened" form: in units of the noise and with the predictors decorrelated,

$$
z = L^\top \beta / \sigma, \qquad s = \lVert z \rVert, \qquad u = z / s .
$$

Then $s^2 = \beta^\top\Sigma\beta/\sigma^2$ is the signal-to-noise ratio, so
$\rho^2 = s^2 / (1 + s^2)$. The direction $u$ is uniform on the unit sphere and
independent of $s$.

"Direction" here means how the effect is divided over the predictors, measured
in the whitened coordinates. That is how correlated predictors are handled
automatically.

### Why Beta(½, 1)? It is the flat correlation prior for one predictor

With one predictor, the population correlation is $r = \pm\sqrt{\rho^2}$ with
a random sign. If $\rho^2 \sim \mathrm{Beta}(\tfrac12, 1)$, then for $0 \le a \le 1$:

$$
P(|r| \le a) = P(\rho^2 \le a^2) = (a^2)^{1/2} = a ,
$$

so $|r|$ is uniform on $(0, 1)$, and $r$ is uniform on $(-1, 1)$. That is
Jeffreys' default prior for a correlation: every possible value equally
plausible. Beta(½, 1) is the only Beta prior on $\rho^2$ with this property,
and (b) keeps it the same for every $k$, so that the implied $\rho^2$ does not drift
toward 1 as predictors are added (see the tutorial, Step 2).

---

## 3. The prior density on the coefficients (estimation)

brms puts a prior directly on $\beta$, so we need its density.

**Step 1: the length $s$.** Transform $\rho^2 \sim \mathrm{Beta}(\tfrac12, 1)$,
whose density is $\tfrac12 (\rho^2)^{-1/2}$, to $s$ via
$\rho^2 = s^2/(1+s^2)$, $\; d\rho^2/ds = 2s/(1+s^2)^2$:

$$
p(s) = \tfrac12 \left(\frac{s^2}{1+s^2}\right)^{-1/2} \frac{2s}{(1+s^2)^2}
     = (1 + s^2)^{-3/2}, \qquad s > 0 .
$$

(Check: $\int_0^\infty (1+s^2)^{-3/2} ds = [\,s/\sqrt{1+s^2}\,]_0^\infty = 1$.)

**Step 2: the whitened vector $z$.** A uniform direction spreads the density of
$s$ evenly over the sphere of radius $s$, whose surface area is
$A_k s^{k-1}$ with $A_k = 2\pi^{k/2}/\Gamma(k/2)$:

$$
p(z) = \frac{(1 + \lVert z\rVert^2)^{-3/2}}{A_k \, \lVert z\rVert^{k-1}} .
$$

For $k = 1$, $A_1 = 2$ and $p(z) = \tfrac12 (1+z^2)^{-3/2}$: the flat-correlation
prior on the effect-size scale (the tutorial's Step 1 curve).

**Step 3: back to $\beta$.** Since $z = L^\top \beta / \sigma$, the Jacobian is
$\lvert\det L^\top\rvert / \sigma^k = \prod_j L_{jj} / \sigma^k$:

$$
\log p(\beta \mid \sigma)
= -\tfrac32 \log(1 + s^2) - (k-1)\log s - \log A_k
  + \sum_j \log L_{jj} - k \log \sigma,
\qquad s = \lVert L^\top \beta \rVert / \sigma .
$$

Note the factor $s^{-(k-1)}$: for $k \ge 2$ the density has an infinitely high
(but integrable) peak at $\beta = 0$. That peak is why brms can show
occasional divergent transitions; `archive/stan_unit_vector_model.R`
samples $(\rho^2, u)$ directly and avoids it.

**In the code:** the Stan function `possible_range_lpdf` in `R/possible_range.R`
(`length_squared` is $s^2$, `covariance_root` is $L$), plus
`target += -log(sigma)` for the noise prior.

---

## 4. The closed-form Bayes factor

We want the Bayes factor of the model with $k$ predictors against the
intercept-only model:

$$
\mathrm{BF}_k = \frac{p(y \mid \text{model with } k \text{ predictors})}{p(y \mid \text{intercept only})} .
$$

The idea: **fix a direction $u$ first.** Given $u$, the effect lies along one
line, so the model becomes a regression on a single constructed predictor.
For one predictor we know the answer, and then we average over directions.

### Step 1: one predictor along direction $u$

Let $Q = X_c (\sqrt{n} L^\top)^{-1}$. Its columns are orthonormal and
$X_c\beta = Q w$ with $w = \sqrt n\, L^\top\beta = \sigma\sqrt n\, s\,u$.
For a fixed direction $u$, the signal is $\sigma \sqrt n\, s\,(Qu)$: one
predictor $v = Qu$ with $\lVert v \rVert = 1$, and coefficient
$\gamma = \sigma\sqrt n\, s$.

Allowing $s$ to carry a sign (equivalent, since $u$ and $-u$ are equally
likely), the coefficient in noise units, $\tilde s = \gamma/(\sigma\sqrt n)$, has density
$\tfrac12(1 + \tilde s^2)^{-3/2}$: the one-predictor flat-correlation prior.

The explained share of this one-predictor regression is

$$
R_u^2 = \frac{(v^\top y_c)^2}{\mathrm{SST}} = R^2\, t^2 ,
$$

where $t$ is the cosine of the angle between $u$ and the direction in which
the data point (the least-squares fit in whitened coordinates). So given $u$,
the data enter only through $R^2$ and $t$.

### Step 2: the one-predictor Bayes factor as a mixture of g-priors

The density $\tfrac12(1+\tilde s^2)^{-3/2}$ is a scale mixture of normals:

$$
\tilde s \mid g \sim N(0, g/n), \qquad g \sim \text{Inverse-Gamma}(1, \; n/2),
$$

since $\int (g/n)^{-1/2} e^{-n\tilde s^2/(2g)}\, g^{-2} e^{-n/(2g)}\, dg \propto (1+\tilde s^2)^{-3/2}$.
Given $g$, the coefficient on $v$ is $\gamma \sim N(0, g\sigma^2)$, which is
exactly Zellner's g-prior for one predictor with $v^\top v = 1$. Its Bayes
factor against the intercept-only model is the standard g-prior result
(Liang et al., 2008), with $R_u^2$ in place of $R^2$:

$$
\mathrm{BF}(g, u) = (1+g)^{(n-2)/2}\,\bigl(1 + g(1 - R_u^2)\bigr)^{-(n-1)/2}
= (1+g)^{-1/2} \left(1 - R^2 t^2 \frac{g}{1+g}\right)^{-(n-1)/2} .
$$

(Integrating $\alpha$ with the flat prior and $\sigma$ with $1/\sigma$ is already
contained in this formula.)

Compare JZS: Zellner–Siow uses $g \sim$ Inverse-Gamma(½, $n r^2/2$), a Cauchy
on the standardized effect. Here the shape is 1 instead of ½, which is a $t_2$
instead of a Cauchy. That is the difference between "flat on the correlation"
and "Cauchy on the effect size" (tutorial, Step 1).

### Step 3: average over $g$ and over directions

For a uniform direction $u$ on the sphere in $k$ dimensions, the cosine $t$
with any fixed direction has density

$$
p_k(t) \propto (1 - t^2)^{(k-3)/2}, \qquad -1 < t < 1 .
$$

Averaging gives the closed form:

$$
\boxed{\;
\mathrm{BF}_k(R^2, n) = \int_0^\infty \int_{-1}^{1}
\mathrm{IG}(g; 1, n/2)\,(1+g)^{-1/2}
\left(1 - R^2 t^2 \frac{g}{1+g}\right)^{-(n-1)/2} p_k(t)\, dt\, dg
\;}
$$

In words: **the Bayes factor with $k$ predictors is the one-predictor
flat-correlation Bayes factor, averaged over every direction the effect could
take.** The data enter only through $R^2$, $n$ and $k$.

For $k = 1$ there are only two directions, both with $t^2 = 1$, and the
formula reduces to the one-predictor flat-correlation Bayes factor.

### Step 4: numerical evaluation

The double integral is evaluated on grids:

- $g$ on a log scale from $10^{-6}$ to $10^{10}$ (800 points); the extra $+\log g$
  in the code is the Jacobian of the log scale;
- $t = \cos\varphi$ with $\varphi$ on a grid over $(0, \pi/2)$ (200 points; $\varphi$ and
  $\pi - \varphi$ give the same $t^2$), so that $p_k(t)\,dt \propto \sin^{k-2}\varphi\, d\varphi$,
  normalized with its exact total $\int_0^{\pi/2} \sin^{k-2}\varphi\,d\varphi$.
  The grid is **packed toward $\varphi = 0$**: $\varphi = \tfrac{\pi}{2} u^3$ with $u$ evenly
  spaced. When the evidence is strong the integrand is a narrow peak at $\varphi = 0$
  (width about $\sqrt{(1 - z)/(2 m z)}$, with $z = R^2 g/(1+g)$ and $m = (n-1)/2$), which an
  even grid steps over;
- everything on the log scale with log-sum-exp, because the terms can be as
  small as $10^{-300}$.

Against a ground truth computed by adaptive integration (`integrate()`), the
largest error over 162 settings ($R^2$ from 0.01 to 0.99, $n$ from 10 to
10,000, $k$ from 1 to 30) is $4 \times 10^{-5}$ on the log scale
(`validation/01`, check A1). An even grid was off by up to 28 on the log scale
at $n = 10{,}000$, $R^2 = 0.99$.

**In the code:** `log_bayes_factor_possible_range_vs_intercept_only()` in `R/possible_range.R`:

| Code | Formula |
|---|---|
| `scale_values` | $g$ |
| `log_scale_weights` | $\log \mathrm{IG}(g; 1, n/2) + \log g$ (density on the log-$g$ grid) |
| `shrink_factor` | $g/(1+g)$ |
| `alignment_squared` | $t^2$ |
| `log_direction_weights` | $\log p_k(t)\,dt$ on the packed angle grid |
| `-0.5 * log1p(scale_values)` | $(1+g)^{-1/2}$ |
| `log_fit_matrix` | $\log\bigl(1 - R^2 t^2 g/(1+g)\bigr)^{-(n-1)/2}$ plus the direction weight |

---

## 5. Comparing any two models

Every model's Bayes factor is taken against the same intercept-only model,
with the same priors on $\alpha$ and $\sigma$. So for any two models A and B
fitted to the same outcome and rows, nested or not:

$$
\mathrm{BF}_{AB} = \frac{\mathrm{BF}_{k_A}(R^2_A, n)}{\mathrm{BF}_{k_B}(R^2_B, n)} .
$$

The evidence for one coefficient (or a block of columns, such as a factor) is
the ratio for the model with and without it. Each model uses the prior with
its own predictors' covariance and its own $k$.

**In the code:** `bayes_factor_possible_range()` and `bayes_factor_possible_range_per_term()`.

---

## 6. Why no MCMC is needed for the Bayes factor

The Bayes factor is a ratio of averaged likelihoods. Everything except
$g$ and $t$ integrates out exactly, so a two-dimensional grid suffices. Posterior
estimates do need MCMC: the posterior of $(\alpha, \sigma, \beta)$ under this
prior has no closed form, which is why estimation goes through brms.

The two routes describe the same prior: brms bridge sampling gives 18.76
against 18.70 from the closed form for the vignette's x1, and 2.833 against
2.821 for a three-level factor tested as a block (`validation/01`, check A5).

---

## 7. Properties

These are checked numerically in `validation/01_checks.R`, part B (invariance: part A), not proved:

| Property | Meaning | Check |
|---|---|---|
| Predictive matching | with $n = k + 1$ observations the data can't separate the models: BF = 1 | 01: B2a |
| Information consistency | BF → ∞ as $R^2 \to 1$ | 01: B2b |
| Model-selection consistency | under the null, $\log \mathrm{BF}_{01}$ grows like $\tfrac12 \log n$ | 01: B2c |
| Invariance | BF depends only on $R^2$, $n$, $k$: unaffected by units or linear recoding of predictors | by construction (§4); 01: A3 |
| Correct reduction | $k = 1$ equals the flat-correlation prior | §2; 01: B1 |

Proofs (for example along the lines of Bayarri et al., 2012) are still open.

---

## 8. Assumptions and limits

- Normal errors, fixed effects, an intercept.
- The prior uses the **sample** covariance of the predictors (it conditions
  on $X$), like the g-prior and JZS.
- The improper priors on $\alpha$ and $\sigma$ are only valid for comparing models
  that share them, which all models here do.
- Not covered: group-level effects, other likelihoods, distributional models.
