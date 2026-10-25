# Derivation of the flat-R prior and its Bayes factor

The whole prior in one sentence:

$$
p(\beta \mid \sigma) \;=\; \underbrace{P(\text{total})}_{1}\;\times\;
\underbrace{P(\text{this combination} \mid \text{total})}_{2}\;\times\;
\underbrace{\text{unit conversion}}_{3}
$$

1. **P(total):** how likely is the total effect? Flat on the multiple correlation $R$.
2. **P(this combination | total):** given the total, every way of dividing it
   over the predictors is equally likely.
3. **Unit conversion:** from "noise units" to the units of your predictors and outcome.

Only factor 1 carries a judgement. Factors 2 and 3 turn "per slice of $R$" into
"per point of $\beta$", which is what estimation software asks for.

This document derives the two formulas that `R/flat_r.R` uses:

1. the prior density on the coefficients, factor by factor (used for estimation in brms), and
2. the closed-form Bayes factor (used by `bayes_factor_flat_r()`).

Each section ends with the code that implements it.

---

## 1. Model and notation

Linear regression with $n$ observations and $k$ predictors:

$$
y = \alpha \mathbf{1} + X_c \beta + \varepsilon, \qquad \varepsilon \sim N(0, \sigma^2 I_n),
$$

where $X_c$ holds the **centered** predictor columns (each column has mean zero).
Centering makes the intercept $\alpha$ the mean of $y$ and separates it from the slopes.

| Symbol | Meaning |
|---|---|
| $\Sigma = X_c^\top X_c / n$ | sample covariance matrix of the predictors |
| $\beta^\top \Sigma \beta$ | variance explained by the predictors (handles correlated predictors by itself) |
| $\sigma^2$ | unexplained variance (the noise) |
| $f^2 = \beta^\top \Sigma \beta / \sigma^2$ | explained / unexplained variance (**Cohen's $f^2$**) |
| $f$ | the **length of $\beta$ in noise units**: one number for the whole vector |
| $\rho^2 = f^2 / (1 + f^2)$ | **population** explained share of variance; $R = \sqrt{\rho^2}$ |
| $R^2$ (sample) | explained share of the fitted regression |
| $\mathrm{SST} = \sum_i (y_i - \bar y)^2$ | total sum of squares |

**Noise units:** an effect measured with the unexplained spread $\sigma$ as the
ruler, like Cohen's $d$ divides a difference by an SD. $f = 1$ means the
predictors explain as much variance as the noise leaves, so $\rho^2 = 0.5$.

The prior is defined conditional on the observed predictors, as in Zellner's
g-prior and the JZS default: the predictors are treated as fixed.

Intercept and noise get the usual neutral priors (as in JZS):
$p(\alpha) \propto 1$, $p(\sigma) \propto 1/\sigma$. They are improper, but
shared by every model being compared, so their constants cancel in Bayes factors.

---

## 2. Factor 1: P(total), flat on R

**The judgement.** The multiple correlation $R$ is flat from 0 to 1: every
overall strength of relationship is equally plausible, for any number of
predictors. Equivalently, $\rho^2 = R^2 \sim \mathrm{Beta}(\tfrac12, 1)$.

**Why flat on $R$.** With one predictor, $R = |r|$, so flat on $R$ is flat on
$r$ from $-1$ to $1$: Jeffreys' default prior for a correlation. For one
predictor this is the stretched Beta$(b, b)$ family of Ly, Verhagen &
Wagenmakers with $b = 1$; $b = \tfrac12$ is what a Cauchy(0, 1) on $\beta/\sigma$
implies (JZS). Putting the same prior on the **total** keeps it fixed as
predictors are added, so the implied $\rho^2$ does not drift towards 1
(tutorial, Step 2).

**Why Beta(½, 1) is the same statement.** If $R$ is uniform on $(0,1)$, then for $0 \le a \le 1$:

$$
P(\rho^2 \le a^2) = P(R \le a) = a = (a^2)^{1/2},
$$

which is the Beta(½, 1) distribution function. Its density is $\tfrac12 (\rho^2)^{-1/2}$.

**Converted to the length $f$.** Estimation needs a density per unit of $f$
(the length of $\beta$), not per unit of $R$. With $R = f / \sqrt{1 + f^2}$:

$$
\frac{dR}{df} = (1 + f^2)^{-3/2}
\qquad\Longrightarrow\qquad
p(f) = \underbrace{1}_{\text{flat on } R} \times (1 + f^2)^{-3/2}, \qquad f > 0 .
$$

Equal slices of $R$ cover longer and longer stretches of $f$ (because $R$
saturates near 1), so the same probability is spread thinner per unit of $f$.
(Check: $\int_0^\infty (1+f^2)^{-3/2} df = [\,f/\sqrt{1+f^2}\,]_0^\infty = 1$.)

**In the code:** `log_p_total = -1.5 * log1p(f_squared)`.

---

## 3. Factor 2: P(this combination | total), every split equally likely

Given the total $f$, many coefficient vectors have that length: with two
predictors a whole circle, with $k$ predictors the surface of a sphere of radius $f$
(measured in noise units, with the overlap between predictors accounted for).

Every combination on that sphere is equally likely, so the probability of one
combination is 1 divided by the size of all combinations:

$$
P(\text{this combination} \mid f) = \frac{1}{A_k \, f^{k-1}},
\qquad A_k = \frac{2\pi^{k/2}}{\Gamma(k/2)} .
$$

| $k$ | All combinations with length $f$ | Size |
|---|---|---|
| 1 | two points, $+f$ and $-f$ | 2 |
| 2 | a circle | $2\pi f$ |
| 3 | a sphere surface | $4\pi f^2$ |

For $k = 1$ factors 1 and 2 give $\tfrac12 (1+f^2)^{-3/2}$: the flat-correlation
prior on the effect-size scale (tutorial, Step 1 curve). This equals a $t$
distribution with 2 degrees of freedom and scale $1/\sqrt2$ on $\beta/\sigma$.

**Consequence: the peak at zero.** For $k \ge 2$ the size $A_k f^{k-1}$ goes to
0 as $f \to 0$, so the density per point has an infinitely high (but integrable)
peak at $\beta = 0$. That is where the shrinkage at small $n$ comes from, and why
brms can show occasional divergent transitions; `archive/stan_unit_vector_model.R`
samples $(\rho^2, u)$ directly and avoids the peak.

**Why only this split works with correlated predictors.** A uniform split looks
the same in every orientation, so the prior does not depend on the order of
the predictors. Other splits (for example a flat Dirichlet on the shares) would
make the prior depend on that order when predictors correlate.

**In the code:**
`log_p_combination_given_total = -(log(2) + 0.5 * k * log(pi()) - lgamma(0.5 * k) + (k - 1) * log(f))`.

---

## 4. Factor 3: unit conversion

Factors 1 and 2 give a density per point in noise units. The coefficients
$\beta$ are in the units of the data (outcome units per predictor unit).
Converting the axes stretches or squeezes space, and the density changes by
the inverse of that stretch:

$$
\text{unit conversion} = \frac{\sqrt{\det \Sigma}}{\sigma^k} .
$$

This factor is the same for every $\beta$: it carries no judgement, only units.
In real units the circles of factor 2 become tilted ellipses (different spreads
and correlations between predictors).

**In the code:** `log_conversion = 0.5 * log_determinant(predictor_covariance) - k * log(sigma)`.

---

## 5. The whole density

$$
\log p(\beta \mid \sigma)
= \underbrace{-\tfrac32 \log(1 + f^2)}_{1}
\;\underbrace{-\,\log A_k - (k-1)\log f}_{2}
\;+\underbrace{\tfrac12 \log\det\Sigma - k \log \sigma}_{3},
\qquad f^2 = \frac{\beta^\top\Sigma\beta}{\sigma^2} .
$$

**In the code:** the Stan function `flat_r_lpdf` in `R/flat_r.R` computes the
three factors as `log_p_total`, `log_p_combination_given_total` and
`log_conversion` and returns their sum; `target += -log(sigma)` adds the noise
prior. Checked against simulation (counting prior draws in small squares) in
`tutorial/`.

**Drawing from the prior** follows the same sentence in reverse:
`draw_coefficients_from_flat_r_prior()` draws $R$ flat, turns it into $f$,
draws a direction with no preference in terms of explained variance (normal
draws with covariance $\Sigma^{-1}$), and scales it to length $f$.

---

## 6. The closed-form Bayes factor

We want the Bayes factor of the model with $k$ predictors against the
intercept-only model:

$$
\mathrm{BF}_k = \frac{p(y \mid \text{model with } k \text{ predictors})}{p(y \mid \text{intercept only})} .
$$

The idea: **fix a direction $u$ first.** Given $u$, the effect lies along one
line, so the model becomes a regression on a single constructed predictor.
For one predictor we know the answer, and then we average over directions.

### Step 1: one predictor along direction $u$

Work with the **fitted values** $X_c\beta$ instead of the coefficients. Take any
orthonormal basis $Q$ ($n \times k$) of the columns of $X_c$, for example
`qr.Q(qr(X_c))`, the same QR decomposition `lm()` uses. Every vector of fitted
values is $X_c\beta = Q w$ for some $w$ with $k$ entries, and

$$
\lVert w \rVert^2 = \lVert X_c \beta \rVert^2 = n\, \beta^\top \Sigma \beta .
$$

So in $w$ the explained variance is ordinary squared length, and the flat-R
prior (§2–3) reads: length $\lVert w \rVert = \sigma\sqrt n\, f$, direction
$u = w / \lVert w \rVert$ uniform on the sphere. Which orthonormal basis you
take does not matter: a uniform direction looks the same in every orientation.

For a fixed direction $u$, the signal is $\sigma \sqrt n\, f\,(Qu)$: one
predictor $v = Qu$ with $\lVert v \rVert = 1$, and coefficient
$\gamma = \sigma\sqrt n\, f$.

Allowing $f$ to carry a sign (equivalent, since $u$ and $-u$ are equally
likely), the coefficient in noise units, $\tilde f = \gamma/(\sigma\sqrt n)$, has density
$\tfrac12(1 + \tilde f^2)^{-3/2}$: the one-predictor flat-correlation prior.

The explained share of this one-predictor regression is

$$
R_u^2 = \frac{(v^\top y_c)^2}{\mathrm{SST}} = R^2\, t^2 ,
$$

where $t$ is the cosine of the angle between $u$ and the direction in which
the data point (the least-squares fit, $Q^\top y_c$). So given $u$,
the data enter only through $R^2$ and $t$.

### Step 2: the one-predictor Bayes factor as a mixture of g-priors

The density $\tfrac12(1+\tilde f^2)^{-3/2}$ is a scale mixture of normals:

$$
\tilde f \mid g \sim N(0, g/n), \qquad g \sim \text{Inverse-Gamma}(1, \; n/2),
$$

since $\int (g/n)^{-1/2} e^{-n\tilde f^2/(2g)}\, g^{-2} e^{-n/(2g)}\, dg \propto (1+\tilde f^2)^{-3/2}$.
Given $g$, the coefficient on $v$ is $\gamma \sim N(0, g\sigma^2)$, which is
exactly Zellner's g-prior for one predictor with $v^\top v = 1$. Its Bayes
factor against the intercept-only model is the standard g-prior result
(Liang et al., 2008), with $R_u^2$ in place of $R^2$:

$$
\mathrm{BF}(g, u) = (1+g)^{(n-2)/2}\,\bigl(1 + g(1 - R_u^2)\bigr)^{-(n-1)/2}
= (1+g)^{-1/2} \left(1 - R^2 t^2 \frac{g}{1+g}\right)^{-(n-1)/2} .
$$

(Integrating $\alpha$ with the flat prior and $\sigma$ with $1/\sigma$ is already
contained in this formula.) The mixing density on $g$ is factors 1 and 2 for one predictor, written as a mixture of normals.

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

**In the code:** `log_bayes_factor_flat_r_vs_intercept_only()` in `R/flat_r.R`:

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

## 7. Comparing any two models

Every model's Bayes factor is taken against the same intercept-only model,
with the same priors on $\alpha$ and $\sigma$. So for any two models A and B
fitted to the same outcome and rows, nested or not:

$$
\mathrm{BF}_{AB} = \frac{\mathrm{BF}_{k_A}(R^2_A, n)}{\mathrm{BF}_{k_B}(R^2_B, n)} .
$$

The evidence for one coefficient (or a block of columns, such as a factor) is
the ratio for the model with and without it. Each model uses the prior with
its own predictors' covariance and its own $k$.

**In the code:** `bayes_factor_flat_r()` and `bayes_factor_flat_r_per_term()`.

---

## 8. Why no MCMC is needed for the Bayes factor

The Bayes factor is a ratio of averaged likelihoods. Everything except
$g$ and $t$ integrates out exactly, so a two-dimensional grid suffices. Posterior
estimates do need MCMC: the posterior of $(\alpha, \sigma, \beta)$ under this
prior has no closed form, which is why estimation goes through brms.

The two routes describe the same prior: brms bridge sampling gives 18.76
against 18.70 from the closed form for the vignette's x1, and 2.833 against
2.821 for a three-level factor tested as a block (`validation/01`, check A5).

---

## 9. Properties

These are checked numerically in `validation/01_checks.R`, part B (invariance: part A), not proved:

| Property | Meaning | Check |
|---|---|---|
| Predictive matching | with $n = k + 1$ observations the data can't separate the models: BF = 1 | 01: B2a |
| Information consistency | BF → ∞ as $R^2 \to 1$ | 01: B2b |
| Model-selection consistency | under the null, $\log \mathrm{BF}_{01}$ grows like $\tfrac12 \log n$ | 01: B2c |
| Invariance | BF depends only on $R^2$, $n$, $k$: unaffected by units or linear recoding of predictors | by construction (§6); 01: A3 |
| Correct reduction | $k = 1$ equals the flat-correlation prior | §2–3; 01: B1 |

Proofs (for example along the lines of Bayarri et al., 2012) are still open.

---

## 10. Assumptions and limits

- Normal errors, fixed effects, an intercept.
- The prior uses the **sample** covariance of the predictors (it conditions
  on $X$), like the g-prior and JZS.
- The improper priors on $\alpha$ and $\sigma$ are only valid for comparing models
  that share them, which all models here do.
- Not covered: group-level effects, other likelihoods, distributional models.
