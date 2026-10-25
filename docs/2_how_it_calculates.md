# 2. How the prior calculates

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

Factors 1 and 2 carry the judgements (why: [`1_from_behaviour.md`](1_from_behaviour.md),
steps 2 and 3). Factor 3 carries none: it only converts units, so that the
prior is a density per point of $\beta$, which is what estimation software asks for.

And the Bayes factor in one sentence:

> **The Bayes factor with $k$ predictors is the one-predictor Bayes factor,
> averaged over every direction the effect could take.**

Sections 1–5 build the density (used by brms for estimation), sections 6–7 the
Bayes factor (used by `bayes_factor_flat_r()`). Each section ends with the code
that implements it.

---

## 1. Model and notation

Linear regression with $n$ observations and $k$ predictors:

$$
y = \alpha \mathbf{1} + X_c \beta + \varepsilon, \qquad \varepsilon \sim N(0, \sigma^2 I_n),
$$

where $X_c$ holds the **centered** predictor columns (each column has mean zero),
so the intercept $\alpha$ is the mean of $y$ and separate from the slopes.

| Symbol | Meaning |
|---|---|
| $\Sigma = X_c^\top X_c / n$ | covariance matrix of the predictors |
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

The prior conditions on the observed predictors (they are treated as fixed),
as in Zellner's g-prior and the JZS default. Intercept and noise get the usual
neutral priors (as in JZS): $p(\alpha) \propto 1$, $p(\sigma) \propto 1/\sigma$.
They are improper, but shared by every model being compared, so their constants
cancel in Bayes factors.

---

## 2. Factor 1: P(total), flat on R

$R$ is flat from 0 to 1, for any number of predictors. Equivalently
$\rho^2 = R^2 \sim \mathrm{Beta}(\tfrac12, 1)$, because for $0 \le a \le 1$:

$$
P(\rho^2 \le a^2) = P(R \le a) = a = (a^2)^{1/2},
$$

which is the Beta(½, 1) distribution function.

**Converted to the length $f$.** Estimation needs a density per unit of $f$,
not per unit of $R$. With $R = f / \sqrt{1 + f^2}$:

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
predictors a whole circle, with $k$ predictors the surface of a sphere of
radius $f$ (in noise units, with the overlap between predictors accounted for).
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
prior on the effect-size scale (figure 1 in `1_from_behaviour.md`), a $t$
distribution with 2 degrees of freedom and scale $1/\sqrt2$ on $\beta/\sigma$.

![](figures/4_three_factors.png)

**The peak at zero.** For $k \ge 2$ the size $A_k f^{k-1}$ goes to 0 as
$f \to 0$, so the density per point has an infinitely high (but integrable)
peak at $\beta = 0$. With many predictors brms under-samples that peak
without warnings (`validation/README.md`, open problem);
`archive/stan_unit_vector_model.R` samples $(\rho^2,$ direction$)$ directly
and avoids the peak.

**In the code:**
`log_p_combination_given_total = -(log(2) + 0.5 * k * log(pi()) - lgamma(0.5 * k) + (k - 1) * log(f))`.

---

## 4. Factor 3: unit conversion

Factors 1 and 2 give a density per point in noise units. The coefficients
$\beta$ are in the units of the data. Converting the axes stretches or squeezes
space, and the density changes by the inverse of that stretch:

$$
\text{unit conversion} = \frac{\sqrt{\det \Sigma}}{\sigma^k} .
$$

It is the same for every $\beta$. In real units the circles of factor 2 become
tilted ellipses (different spreads and correlations between predictors).

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

**In the code:** the Stan function `flat_r_lpdf` in `R/flat_r.R` returns the
sum of `log_p_total`, `log_p_combination_given_total` and `log_conversion`;
`target += -log(sigma)` adds the noise prior.

**Drawing from the prior** is the same sentence read forwards:
`draw_coefficients_from_flat_r_prior()` draws $R$ flat (factor 1), a direction
with no preference in explained variance (factor 2), and scales it to length $f$.

**Checks.** The formula and the drawing recipe describe the same prior:
counting draws in small squares gives the formula's density (figure below,
correlated predictors, so factor 3 is tested too). And brms, sampling from
`flat_r_lpdf` itself, gives back the flat $R$ and the per-predictor shares
(`validation/01_prior_in_brms.R`).

![](figures/5_formula_vs_draws.png)

---

## 6. The closed-form Bayes factor

We want the Bayes factor of the model with $k$ predictors against the
intercept-only model,
$\mathrm{BF}_k = p(y \mid k \text{ predictors}) / p(y \mid \text{intercept only})$.

**The one-sentence version:** fix a direction first. Along one direction the
model is a regression on a single constructed predictor, whose Bayes factor we
know. Then average over directions.

### Step 1: one predictor along a fixed direction

Work with the **fitted values** $X_c\beta$ instead of the coefficients. Take any
orthonormal basis $Q$ ($n \times k$) of the columns of $X_c$, for example
`qr.Q(qr(X_c))`, the QR decomposition `lm()` uses. Every vector of fitted
values is $X_c\beta = Q w$ for some $w$ with $k$ entries, and

$$
\lVert w \rVert^2 = \lVert X_c \beta \rVert^2 = n\, \beta^\top \Sigma \beta .
$$

So in $w$ the explained variance is ordinary squared length, and the flat-R
prior reads: length $\lVert w \rVert = \sigma\sqrt n\, f$, direction
$w / \lVert w \rVert$ uniform on the sphere. Which orthonormal basis does not
matter: a uniform direction looks the same in every orientation.

For a fixed direction $d$ (a unit vector), the signal is $\sigma \sqrt n\, f\,(Qd)$:
one predictor $v = Qd$ with $\lVert v \rVert = 1$. Letting $f$ carry a sign
($d$ and $-d$ are equally likely), the coefficient in noise units has density
$\tfrac12(1 + f^2)^{-3/2}$: the one-predictor flat-correlation prior.

The explained share of this one-predictor regression is

$$
R_d^2 = \frac{(v^\top y_c)^2}{\mathrm{SST}} = R^2\, c^2 ,
$$

where $c$ is the **alignment**: the cosine of the angle between $d$ and the
direction the data point in (the least-squares fit, $Q^\top y_c$). Given $d$,
the data enter only through $R^2$ and $c$.

### Step 2: the one-predictor Bayes factor, via g

The one-predictor Bayes factor needs an integral over $\sigma$ (and $\alpha$).
For a *normal* prior on the coefficient that integral has a known closed form:
Zellner's g-prior, where $g$ is the prior variance of the coefficient relative
to the noise. So the flat-correlation prior is written as a mixture of normals:

$$
f \mid g \sim N(0, g/n), \qquad g \sim \text{Inverse-Gamma}(1, \; n/2),
$$

since $\int (g/n)^{-1/2} e^{-nf^2/(2g)}\, g^{-2} e^{-n/(2g)}\, dg \propto (1+f^2)^{-3/2}$.
For each $g$ the g-prior Bayes factor is known (Liang et al., 2008), with
$R_d^2$ in place of $R^2$:

$$
\mathrm{BF}(g, d) = (1+g)^{(n-2)/2}\,\bigl(1 + g(1 - R_d^2)\bigr)^{-(n-1)/2}
= (1+g)^{-1/2} \left(1 - R^2 c^2 \frac{g}{1+g}\right)^{-(n-1)/2} .
$$

The Zellner–Siow (JZS) default uses Inverse-Gamma(½, $n r^2/2$) instead: a
Cauchy on the effect. Shape 1 instead of ½ is the difference between "flat on
the correlation" and "Cauchy on the effect size".

### Step 3: average over g and over directions

For a uniform direction in $k$ dimensions, the alignment $c$ with any fixed
direction has density

$$
p_k(c) \propto (1 - c^2)^{(k-3)/2}, \qquad -1 < c < 1 .
$$

Averaging gives the closed form:

$$
\boxed{\;
\mathrm{BF}_k(R^2, n) = \int_0^\infty \int_{-1}^{1}
\mathrm{IG}(g; 1, n/2)\,(1+g)^{-1/2}
\left(1 - R^2 c^2 \frac{g}{1+g}\right)^{-(n-1)/2} p_k(c)\, dc\, dg
\;}
$$

The data enter only through $R^2$, $n$ and $k$. For $k = 1$ there are only two
directions, both with $c^2 = 1$, and the formula is the one-predictor
flat-correlation Bayes factor.

### Step 4: numerical evaluation

The double integral is evaluated on grids:

- $g$ on a log scale from $10^{-6}$ to $10^{10}$ (800 points); the extra $+\log g$
  in the weights is because the step is in $\log g$;
- $c = \cos\varphi$ with the angle $\varphi$ on a grid over $(0, \pi/2)$ (200 points;
  $\varphi$ and $\pi - \varphi$ give the same $c^2$), so that
  $p_k(c)\,dc \propto \sin^{k-2}\varphi\, d\varphi$, normalized with its exact total.
  The grid is **packed toward $\varphi = 0$** ($\varphi = \tfrac{\pi}{2} q^3$ with $q$
  evenly spaced): when the evidence is strong the integrand is a narrow peak
  there, which an even grid steps over;
- everything on the log scale with log-sum-exp, because terms can be as small
  as $10^{-300}$.

Against adaptive integration (`integrate()`), the largest error over 162
settings ($R^2$ 0.01–0.99, $n$ 10–10,000, $k$ 1–30) is $4 \times 10^{-5}$ on
the log scale (`validation/02_bayes_factor.R`, check A1).

**In the code:** `log_bayes_factor_flat_r_vs_intercept_only()` in `R/flat_r.R`:

| Code | Formula |
|---|---|
| `g_values` | $g$ |
| `log_g_weights` | $\log \mathrm{IG}(g; 1, n/2) + \log g$ |
| `shrink_factor` | $g/(1+g)$ |
| `alignment_squared` | $c^2$ |
| `log_direction_weights` | $\log p_k(c)\,dc$ on the packed angle grid |
| `-0.5 * log1p(g_values)` | $(1+g)^{-1/2}$ |
| `log_fit_matrix` | $\log\bigl(1 - R^2 c^2 g/(1+g)\bigr)^{-(n-1)/2}$ plus the direction weight |

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

Everything except $g$ and $c$ integrates out exactly, so the Bayes factor needs
no MCMC. Posterior estimates do, which is why estimation goes through brms.
Both routes describe the same prior: brms bridge sampling agrees with the
closed form (`validation/02_bayes_factor.R`, check A3).

**In the code:** `bayes_factor_flat_r()` and `bayes_factor_flat_r_per_term()`.

---

## 8. Properties

Checked numerically in `validation/02_bayes_factor.R`, not proved:

| Property | Meaning | Check |
|---|---|---|
| Invariance | BF depends only on $R^2$, $n$, $k$: unaffected by units or linear recoding of predictors | by construction (§6); A2 |
| Predictive matching | with $n = k + 1$ observations the data can't separate the models: BF = 1 | B1 |
| Information consistency | BF → ∞ as $R^2 \to 1$ | B2 |
| Model-selection consistency | under the null, $\log \mathrm{BF}_{01}$ grows like $\tfrac12 \log n$ | B3 |
| Correct reduction | $k = 1$ is the flat-correlation prior | §3 |

Proofs (for example along the lines of Bayarri et al., 2012) are still open.

---

## 9. Assumptions and limits

- Normal errors, fixed effects, an intercept.
- The prior uses the **sample** covariance of the predictors (it conditions
  on $X$), like the g-prior and JZS.
- The improper priors on $\alpha$ and $\sigma$ are only valid for comparing models
  that share them, which all models here do.
- Not covered: group-level effects, other likelihoods, distributional models.
