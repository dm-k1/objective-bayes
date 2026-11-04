# 2. How the prior calculates

> **Which priors this document is about.** All three: flat partial r
> (`prior = "flat_partial_r"`), flat R (`"flat_r"`) and flat R-squared
> (`"flat_r_squared"`). They are one family, $R^2 \sim \mathrm{Beta}(a, 1)$ with
> a uniform direction, $a = k/2$, ½ and 1, and differ only in factor 1. Where a
> derivation is written out for one shape (flat R, $a = ½$, in sections 7–8),
> the general case follows it. Where they come from:
> [`1_from_behaviour.md`](1_from_behaviour.md); how they compare:
> [`4_comparison.md`](4_comparison.md).

The whole prior in one sentence:

$$
p(\beta \mid \sigma) \;=\; \underbrace{P(\text{total})}_{1}\;\times\;
\underbrace{P(\text{this combination} \mid \text{total})}_{2}\;\times\;
\underbrace{\text{unit conversion}}_{3}
$$

1. **P(total):** how likely is the total effect? $R^2 \sim \mathrm{Beta}(a, 1)$:
   flat on the correlations ($a = k/2$), flat on $R$ ($a = ½$) or flat on $R^2$ ($a = 1$).
2. **P(this combination | total):** given the total, every way of dividing it
   over the predictors is equally likely.
3. **Unit conversion:** from "noise units" to the units of your predictors and outcome.

Factors 1 and 2 carry the judgements (why: [`1_from_behaviour.md`](1_from_behaviour.md),
steps 2 and 3). Factor 3 carries none: it only converts units, so that the
prior is a density per point of $\beta$, which is what estimation software asks for.

And the Bayes factor in one sentence:

> **The Bayes factor with $k$ predictors is the one-predictor Bayes factor,
> averaged over every direction the effect could take.** For any shape $a$ the
> same number comes from testing the model's $k$ columns as one block against
> none (section 7, the block route).

| Sections | What | Used by |
|---|---|---|
| 1–5 | the density, factor by factor | `draw_coefficients_from_prior()` |
| 6 | how estimation samples it | `fit_regression()` |
| 7, 9 | the Bayes factor for whole models | `bayes_factor_models()`, `log_bayes_factor_vs_intercept_only()` |
| 8 | the same prior as a mixture of g-priors | (theory, and check A4) |
| 10 | the test for one coefficient (the conditional test) | `bayes_factor_per_term()`, `log_bayes_factor_conditional()` |
| 11–12 | the default-prior criteria; assumptions and limits | |

Where this fits in the whole project: [`0_overview.md`](0_overview.md).

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

## 2. Factor 1: P(total), R-squared ~ Beta(a, 1)

The total gets $\rho^2 = R^2 \sim \mathrm{Beta}(a, 1)$: $P(R^2 \le x) = x^a$, so
$P(R \le u) = u^{2a}$ and the density on $R$ is $2a\,R^{2a-1}$.

| Prior | $a$ | Density on $R$ | In words |
|---|---|---|---|
| flat partial r | $k/2$ | $k R^{k-1}$ | every correlation matrix of $(y, x_1 \ldots x_k)$ equally likely |
| flat R | ½ | 1 | $R$ flat for every $k$: $P(R \le u) = u$ |
| flat R-squared | 1 | $2R$ | $R^2$ flat for every $k$ |

With one predictor flat partial r and flat R are the same ($a = ½$): Jeffreys'
flat correlation.

**Converted to the length $f$.** Estimation needs a density per unit of $f$,
not per unit of $R$. With $R = f / \sqrt{1 + f^2}$:

$$
\frac{dR}{df} = (1 + f^2)^{-3/2}
\qquad\Longrightarrow\qquad
p(f) = 2a\left(\frac{f^2}{1+f^2}\right)^{a - 1/2} (1 + f^2)^{-3/2}, \qquad f > 0,
$$

and $P(|f| < f_0) = \bigl(f_0^2/(1+f_0^2)\bigr)^a$. For flat R ($a = ½$) this is
$(1+f^2)^{-3/2}$: equal slices of $R$ cover longer and longer stretches of $f$
(because $R$ saturates near 1), so the same probability is spread thinner per
unit of $f$. (Check: $\int_0^\infty (1+f^2)^{-3/2} df = [\,f/\sqrt{1+f^2}\,]_0^\infty = 1$.)
For flat R-squared ($a = 1$) it is $2f(1+f^2)^{-2}$, which is 0 at $f = 0$.

All three have the same tail, $p(f) \approx 2a f^{-3}$ for large $f$: that is
why they share the information-consistency boundary (section 11).

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

For $k = 1$ factors 1 and 2 give $\tfrac12 (1+f^2)^{-3/2}$ for flat R and flat
partial r: the flat-correlation prior on the effect-size scale (figure 1 in
`1_from_behaviour.md`), a $t$ distribution with 2 degrees of freedom and scale
$1/\sqrt2$ on $\beta/\sigma$. For flat R-squared they give $|f|(1+f^2)^{-2}$.

![](figures/4_three_factors.png)

**The peak at zero.** For $k \ge 2$ the size $A_k f^{k-1}$ goes to 0 as
$f \to 0$. Factor 1 goes to 0 like $f^{2a-1}$, so the density per point behaves
like $f^{2a-k}$ near $\beta = 0$: an infinitely high (but integrable) peak for
flat R ($f^{1-k}$) and for flat R-squared with $k \ge 3$ ($f^{2-k}$), none for
flat partial r ($f^0$). A sampler that moves through $\beta$ directly misses
that peak for larger $k$, without warnings; that is why estimation samples
$R$ and the direction instead (section 6).

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

---

## 5. The whole density

$$
\log p(\beta \mid \sigma)
= \underbrace{\log(2a) + \left(a - \tfrac12\right)\log\frac{f^2}{1+f^2} - \tfrac32 \log(1 + f^2)}_{1}
\;\underbrace{-\,\log A_k - (k-1)\log f}_{2}
\;+\underbrace{\tfrac12 \log\det\Sigma - k \log \sigma}_{3},
\qquad f^2 = \frac{\beta^\top\Sigma\beta}{\sigma^2} .
$$

**Drawing from the prior** is the same sentence read forwards:
`draw_coefficients_from_prior()` draws $R = U^{1/(2a)}$ with $U$ flat (factor 1), a direction
with no preference in explained variance (factor 2), and scales it to length $f$
(factor 3). Counting those draws in small squares gives the formula's density
(figure below; correlated predictors, so factor 3 is tested too).

![](figures/5_formula_vs_draws.png)

---

## 6. Estimation: sampling in the prior's own coordinates

The density above is what the prior *says*; it is not how estimation should
*sample* it. Near $\beta = 0$ the density has the peak of section 3, and a
sampler that moves through $\beta$ steps over it without warning (with 8
predictors it draws far too few small effects).

So `inst/stan/flatcor.stan` samples the parameters the prior is built from:

| Parameter | Prior | Factor |
|---|---|---|
| `R` | density $2a\,\lvert R\rvert^{2a-1}$ on (0, 1) (on (−1, 1) when $k = 1$: the sign) | 1 |
| `direction` $= z/\lVert z\rVert$ | uniform on the sphere | 2 |
| `sigma`, `alpha` | $1/\sigma$, flat | |

and computes $\beta = \sigma f\, L^{-\top} d$, with $f = R/\sqrt{1-R^2}$, $d$
the direction and $L$ a square root of $\Sigma$ (Cholesky). Then
$\beta^\top \Sigma \beta = \sigma^2 f^2$, so the length is right (factor 3).
Here a matrix root *is* needed: it turns a uniform direction in explained
variance into predictor units. Which root does not matter.

The direction is $z/\lVert z\rVert$ for a free vector $z$. Any density for $z$
that depends only on its length makes the direction uniform; the model uses one
that keeps $\lVert z\rVert$ near 1. Stan's built-in `unit_vector` does the same
with a standard normal $z$, which lets $\lVert z\rVert$ come near 0, where the
direction is undefined; with 2 predictors that gave 1–3 divergent transitions
per 4,000 draws, and none with the length kept near 1.

**Three priors, one model.** The prior enters only through factor 1: the
data value `r_squared_shape` is $a$, and one line adds
`target += (2 * r_squared_shape - 1) * log(abs(R))` (nothing for $a = ½$).
Factors 2 and 3 are the same for all three.

In these coordinates there is no peak. `validation/01_prior_and_posterior.R`
checks that the prior comes out right for $k$ = 1, 3 and 8 for each prior
(P1–P10), and that the posterior is calibrated (simulation-based calibration,
$k = 8$, all three).

---

## 7. The Bayes factor for whole models

We want the Bayes factor of the model with $k$ predictors against the
intercept-only model,
$\mathrm{BF}_k = p(y \mid k \text{ predictors}) / p(y \mid \text{intercept only})$.
Steps 1–4 build it for flat R ($a = ½$); step 5 gives the route that works for
every $a$, which checks the closed forms and computes flat R-squared with one
predictor.

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

So in $w$ the explained variance is ordinary squared length, and the flat R
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
direction has density (its square $c^2$ is the **share** of
[`1_from_behaviour.md`](1_from_behaviour.md), step 4: Beta(½, (k−1)/2))

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

**In the code:** `log_bayes_factor_flat_r_vs_intercept_only()` in `R/bayes_factor.R`
(internal; `log_bayes_factor_vs_intercept_only(..., prior = "flat_r")` calls it):

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

### Step 5: any shape a, the block route

Testing the model's $k$ columns as one block against none is the same question
as the whole model against the intercept-only model. Section 10 gives the
block likelihood ratio for a length $f$ in noise units with a uniform
direction; with $q = k$ it is

$$
\frac{p(y \mid f)}{p(y \mid 0)} = e^{-nf^2/2}\, {}_1F_1\!\left(\tfrac{n-1}{2}; \tfrac k2; \tfrac{n R^2 f^2}{2}\right),
$$

(there $m = n - k + q - 2 = n - 2$), and the Bayes factor is its average over
factor 1,

$$
\mathrm{BF}_k = \int_0^\infty e^{-nf^2/2}\, {}_1F_1\!\left(\tfrac{n-1}{2}; \tfrac k2; \tfrac{n R^2 f^2}{2}\right)
\, dP(|f| < f), \qquad P(|f| < f) = \left(\frac{f^2}{1+f^2}\right)^{a}.
$$

The code takes the average on a grid in $f$ around the peak of the likelihood
ratio, each cell weighted by its exact prior probability (a difference of the
CDF), as the conditional test does.

- **It reproduces the closed forms:** with $a = ½$ it gives flat R's Bayes
  factor of steps 1–4, with $a = k/2$ flat partial r's (section 8), within
  $2\times10^{-4}$ on the log scale (`tests/testthat/test-three-priors.R`;
  `validation/02`, A4). So it is a check on the closed forms, and the route for
  flat R-squared with one predictor, where no g-mixture exists (section 8); the
  closed forms are used where they exist because they are faster.
- **At $n = k + 1$** every model fits the data exactly ($R^2 = 1$) and the data
  say nothing about the effect; the Bayes factor is 1 for any proper prior on
  it (predictive matching, Bayarri et al., 2012). The closed forms give 1 there;
  the block route returns 1 directly, because its grid needs $n \ge k + 2$.
- **Very strong evidence.** The average over directions is a Bessel function,
  $I_\nu(t)$ with $\nu = k/2 - 1$. Above $t = 10^4$ R's `besselI()` loses
  accuracy and returns 0 beyond about $10^5$, so there the large-argument
  expansion $I_\nu(t)\,e^{-t} = (2\pi t)^{-1/2}\bigl(1 - \tfrac{4\nu^2-1}{8t} +
  \tfrac{(4\nu^2-1)(4\nu^2-9)}{2(8t)^2} - \dots\bigr)$ is used (error on the log scale at $t = 10^4$:
  below $10^{-12}$ up to 15 columns, $4\times10^{-10}$ at 30, $10^{-7}$ at 60). This also keeps the block conditional test exact at
  extreme evidence.

**In the code:** `log_bayes_factor_vs_intercept_only(r_squared, n, k, prior)`:
flat R and flat partial r through their closed forms; flat R-squared through
its g-mixture for $k \ge 3$ (section 8; the hyperprior grid is computed once per
$k$ and kept), the flat partial r closed form at $k = 2$ (the same prior there),
and `log_bayes_factor_shape_vs_intercept_only()` (the block route, internal) at
$k = 1$.

---

## 8. The same prior as a mixture of g-priors

**In one sentence:** the family $R^2 \sim$ Beta($a$, 1) with a uniform
direction is a family of mixtures of Zellner g-priors (for $a \le k/2$), and
flat R is the only member whose implied $R^2$ is the same for every $k$ *and*
equals Jeffreys' flat correlation at $k = 1$.

Zellner's g-prior is $\beta \mid g, \sigma \sim N\bigl(0,\ g\,\sigma^2 (X_c^\top X_c)^{-1}\bigr)$:
a uniform direction (in explained variance) with $f^2 = (g/n)\,\chi^2_k$. Putting a
prior on $g$ gives a *mixture of g-priors* (Liang et al., 2008); Zellner–Siow
(JZS), hyper-g, hyper-g/n, beta-prime, robust and intrinsic priors are all
members. The flat R prior is the member with

$$
\frac{g}{n} = B \cdot W, \qquad B \sim \mathrm{Beta}\!\left(\tfrac12, \tfrac{k-1}{2}\right),
\qquad W \sim \text{Inverse-Gamma}\!\left(1, \tfrac12\right).
$$

- $W$ is the one-predictor part: for $k = 1$, $B = 1$ and this is the
  Inverse-Gamma(1, $n/2$) of section 7, step 2.
- $B$ is the **share** again: the squared alignment $c^2$ of section 7 and the
  share of [`1_from_behaviour.md`](1_from_behaviour.md), step 4. Averaging over
  directions there is the same as averaging over $B$ here.
- Why it works: $f^2 = B W \chi^2_k$, and $B\,\chi^2_k$ has the distribution of
  $\chi^2_1$ (a Beta(½, (k−1)/2) share of a $\chi^2_k$ is a $\chi^2_1$). So
  $f^2 = W\chi^2_1$ for every $k$: exactly the one-predictor prior on the total.
- Constant in $k$ alone does not single it out. Any $W$ gives an implied $R^2$
  that does not depend on $k$; with $W \sim$ Inverse-Gamma($b$, ½) it is
  $R^2 \sim$ Beta(½, $b$) for every $k$ ($W\chi^2_1$ is then a ratio of a
  Gamma(½) and a Gamma($b$) variable; checked in `validation/09`, V5). What fixes $b = 1$ is the second condition:
  at $k = 1$ the prior must be flat on $r$, i.e. $r^2 \sim$ Beta(½, 1). More
  generally, for a g-mixture $f^2 = (g/n)\,\chi^2_k$, so the distribution of
  $f^2$ at each $k$ determines that of $g/n$ (Mellin transforms factor); fixing
  the $k = 1$ distribution to Jeffreys' and asking it to stay the same for every
  $k$ leaves exactly one hyperprior per $k$, this one.

The hyperprior density is $p(g/n = t) = \dfrac{1}{2kt^2}\, {}_1F_1\!\left(\tfrac32;
\tfrac k2 + 1; -\tfrac{1}{2t}\right)$ (Kummer's confluent hypergeometric
function). That factor is why it is not a special case of the tCCH family of
Li & Clyde (2018), which covers hyper-g, beta-prime, robust and intrinsic: their
densities are products of powers and an exponential.

**The Bayes factor** is then the standard g-prior Bayes factor averaged over
this hyperprior, a single integral over $g$:

$$
\mathrm{BF}_k = \int_0^\infty (1+g)^{(n-1-k)/2}\,\bigl(1+g(1-R^2)\bigr)^{-(n-1)/2}\, p(g)\, dg ,
$$

which agrees with section 7 to $3.5\times10^{-4}$ on the log scale
(`validation/02_bayes_factor.R`, check A4, flat R).

**Where it sits in the family.** Median implied $R^2$ before any data ($n = 100$):

| Prior | $k = 1$ | $k = 3$ | $k = 10$ | $k = 30$ |
|---|---|---|---|---|
| Zellner–Siow (scale 1) | 0.50 | 0.84 | 0.95 | 0.98 |
| hyper-g/n ($a = 3$) | 0.54 | 0.87 | 0.96 | 0.99 |
| robust (Bayarri et al., 2012) | 0.52 | 0.72 | 0.77 | 0.78 |
| intrinsic | 0.39 | 0.59 | 0.65 | 0.66 |
| **flat R** | 0.25 | 0.25 | 0.25 | 0.25 |
| **flat R-squared** ($a = 1$) | 0.50 | 0.50 | 0.50 | 0.50 |
| **flat partial r** ($a = k/2$) | 0.25 | 0.63 | 0.87 | 0.95 |

The robust and intrinsic priors also scale $g$ with $k$, and drift less, but
still drift.

**The whole family.** The same construction with $R^2 \sim$ Beta($a$, 1)
and a uniform direction is the g-mixture with $B \sim$ Beta($a$, $k/2 - a$) and
the same $W$, for $a \le k/2$ (checked in `validation/02`, A4, and
`validation/09`, V2):

| Prior | $B$ | A g-mixture? |
|---|---|---|
| flat R ($a = ½$) | Beta(½, (k−1)/2) | for every $k$ |
| flat R-squared ($a = 1$) | Beta(1, $k/2 - 1$) | for $k \ge 3$ (the package computes it this way); at $k = 2$, $B \equiv 1$ and it is flat partial r; at $k = 1$ no ($a > k/2$: it would need $B > 1$), and the block route of section 7 is used |
| flat partial r ($a = k/2$) | $B \equiv 1$: $g/n \sim$ Inverse-Gamma(1, ½) for every $k$ | for every $k$ |

At $a = k/2$ the prior is the one you get when every correlation matrix of
$(y, x_1 \ldots x_k)$ is equally likely, and rstanarm `stan_lm`'s $R^2$ prior
with η = 1. Zellner–Siow (scale 1) is $R^2 \sim$ Beta($k/2$, ½) exactly (its
row in the table). How the three compare:
[`relation_to_others.md`](relation_to_others.md) and [`4_comparison.md`](4_comparison.md).

![](figures/6_g_prior_family.png)

---

## 9. Comparing any two models

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
no MCMC. It agrees with a brute-force Bayes factor that draws from the prior and
uses no formula (`validation/02_bayes_factor.R`, check A3).

**In the code:** `bayes_factor_models()`.

---

## 10. Testing one coefficient: the conditional test

**In one sentence:** to test one coefficient, look only at its own contribution
(its partial correlation), with the prior that the joint prior implies for one
coefficient. A term with several columns is tested the same way, as a
block (end of this section).

Comparing the models with and without x1 (the *joint* test) judges x1 by how it
moves the total, so it depends on how much the other predictors explain.
Example: x1 has partial correlation 0.2 ($n = 100$, $k = 3$) in two datasets.
When the others explain 5%, the joint flat R Bayes factor for x1 is 2.25; when they
explain 60%, it is 0.70, with the same evidence about x1. This is the
**conditional Lindley paradox** of Som, Hans and MacEachern (2016): every prior
with one scale shared by all coefficients has it, including the JZS default
(weaker; [`5_validation.md`](5_validation.md), study 04). Its cause is the shared scale.

The conditional test gives the tested coefficient its own scale: its share of
the total ([`1_from_behaviour.md`](1_from_behaviour.md), step 4). With $f$ its effect in noise units,

$$
f^2 = \text{share}\cdot\frac{R^2}{1-R^2}, \qquad R^2 \sim \mathrm{Beta}(a, 1), \qquad
\text{share} \sim \mathrm{Beta}\!\left(\tfrac12, \tfrac{k-1}{2}\right),
$$

with a random sign. Given the share $s$, $|f| = \sqrt s\, R/\sqrt{1-R^2}$, so

$$
P(|f| < f_0 \mid s) = P\!\left(R^2 < \frac{f_0^2}{s + f_0^2}\right) = \left(\frac{f_0^2}{s + f_0^2}\right)^{a},
$$

and the prior's CDF is this averaged over the share (a grid even in
logit(share), each cell with its exact Beta probability). This is exactly what the joint prior implies for one
coefficient, whatever the correlations between predictors. The other
coefficients and the intercept get flat priors and the noise $1/\sigma$, the
same under both hypotheses. Integrating them out:

$$
\frac{p(y \mid f)}{p(y \mid 0)} = e^{-nf^2/2}\,\frac{I(\sqrt n\, r f)}{I(0)},
\qquad I(b) = \int_0^\infty v^{m} e^{-v^2/2 + bv}\,dv, \qquad m = n - k - 1,
$$

with $r$ the sample partial correlation. The Bayes factor averages this over the
prior for $f$, and depends only on $r$, $n$ and $k$. With one predictor it equals
the joint test.

**Under flat partial r** ($a = k/2$) the mixture is exactly simpler: one
coefficient's partial correlation is flat on (−1, 1) for every $k$
([`relation_to_others.md`](relation_to_others.md), "Why η = 1"), so
$P(|f| < f) = f/\sqrt{1+f^2}$, Jeffreys' one-predictor prior, whatever $k$ is.
The code uses that closed form; the share mixture with $a = k/2$ reproduces it
within $10^{-4}$ (`validation/05`, C4b). Under flat R and flat R-squared the
prior narrows as $k$ grows (study 00, versions B and C): the prior expects each
coefficient to be a small share of the total, on average $1/(3k)$ and $1/(2k)$
of the variance.

How the three tests behave across designs (separating nulls from small
effects, confirming nulls, calling small effects null) is in studies 04 and 06.
The three conditional tests rank datasets identically (their AUC for null vs
0.1 is the same, 0.001 to 0.045 above JZS), but they set the threshold
differently: a true null gets BF₀₁ > 3 in 86–88% of datasets under flat
partial r for every $k$, against 86% → 26% (flat R-squared) and 55% → 6%
(flat R) from $k = 2$ to 20. Summary in [`4_comparison.md`](4_comparison.md).

**In the code:** `log_bayes_factor_conditional(partial_r, n, k, prior)`, used by
`bayes_factor_per_term()` for every one-column term.

### Terms with several columns: the block test

A factor with $q + 1$ levels, or an interaction with $q$ columns, is one
question: does the term matter? It is tested as a block of $q$ columns, in the
same way. Now $f$ is the length of the block's effect in noise units (its
$f^2$ is the variance it adds beyond the other predictors, over $\sigma^2$),
and its direction within the block has no preference. The implied prior,
checked by simulation in `validation/05`, check C9:

$$
\text{flat R, flat R-squared:}\quad f^2 = \text{share}\cdot\frac{R^2}{1-R^2},\quad
\text{share} \sim \mathrm{Beta}\!\left(\tfrac q2, \tfrac{k-q}{2}\right);
\qquad
\text{flat partial r:}\quad \frac{f^2}{1+f^2} \sim \mathrm{Beta}\!\left(\tfrac q2, 1\right)
\text{ for every } k.
$$

The second is the block's population partial $R^2$: under flat partial r
it is Beta($q/2$, 1) whatever $k$ is, as one coefficient's partial correlation
is flat whatever $k$ is (the case $q = 1$). Integrating out the other
coefficients, the intercept and $\sigma$ as before:

$$
\frac{p(y \mid f)}{p(y \mid 0)} = e^{-nf^2/2}\,\frac{\mathrm E_c\, I(\sqrt n\, r f c)}{I(0)}
= e^{-nf^2/2}\; {}_1F_1\!\left(\tfrac{m+1}{2}; \tfrac q2; \tfrac{n r^2 f^2}{2}\right),
\qquad m = n - k + q - 2,
$$

with $r^2$ the block's sample partial $R^2$ (from the models with and without
it) and $c$ the cosine between the effect's direction and the data's, averaged
over directions. For $q = 1$, $c = \pm 1$ and this is the formula above. So the
evidence depends only on $r^2$, $n$, $k$ and $q$, and not on what the other
predictors explain. The code evaluates the integral over $v$ with the average
over directions in closed form (a Bessel function, with its large-argument
expansion for very strong evidence: section 7, step 5), which stays accurate
for large $n$; it agrees with `stats::integrate()` within $1\cdot10^{-4}$ on
the log scale for $q$ = 2 to 5, $n$ = 50 to $10^6$ (check C10, every prior),
and the likelihood ratio agrees with a direct integration over $\sigma$ and
the direction from simulated data (check C11). With $q = k$ the block test is
the whole-model Bayes factor (section 7, step 5).

**In the code:** `log_bayes_factor_conditional()` with `number_of_columns = q`, used by
`bayes_factor_per_term()` for every term with several columns (shown as
"conditional (block)").

**The support curve** is this likelihood ratio itself, as a function of the
effect size: `support_per_partial_correlation()` for one column (signed, over
the partial correlation) and `support_per_partial_r_squared()` for a term with
any number of columns (over its partial $R^2$). It shows which effect sizes the
data can and cannot tell apart from zero; the conditional Bayes factor is this
curve averaged over the prior (`tests/testthat/test-support-curve.R`).

---

**One coherent model space for flat partial r; not for flat R-squared or flat R.**
The conditional test puts flat priors on the other coefficients, so it is not a
ratio of marginal likelihoods under the joint prior, and a coefficient's
conditional Bayes factor differs from `bayes_factor_models(full, without it)`
(the joint test). How far is it from "one prior for everything"? Under a joint
model, the Bayes factor for a zero partial correlation of x1 is the
Savage–Dickey ratio: the prior density of that partial correlation at 0 over
its posterior density at 0, estimated from `fit_regression()` draws
(`validation/10`, part 2; run of 4 November 2026). On 20 datasets per prior
(k = 3, 8; n = 50, 200):

- **Flat partial r** (prior density 1/2 at 0, flat on (−1, 1)): the conditional
  log Bayes factor is within 0.09 of Savage–Dickey wherever the density at 0
  can be estimated precisely (16 datasets; mean absolute difference 0.03, a few
  percent on the Bayes factor), and within Monte Carlo error for the other four
  (log BF₁₀ 2 to 9). The joint test is further away (mean absolute 0.23, up to
  0.42), because the smaller model's own prior is not the joint prior given a
  zero partial correlation. So the conditional test *is* the Bayes factor of
  one model with one prior, for the hypothesis "this partial correlation is
  0", to within a few percent.
- **Flat R-squared** (prior density at 0 from the share mixture, π/4 at k = 3):
  close at k = 3 (mean difference −0.04), but not at k = 8 (mean −0.48 at
  n = 50, −0.34 at n = 200; largest 0.87 on the log scale, a factor 2.4; 19
  estimable datasets, mean absolute 0.28). This is check P4, the one failure
  of the run, and it is a property of the prior, not an error: Savage–Dickey
  equals the conditional test only when the null's prior for the other
  parameters is the joint prior given a zero partial correlation. Under flat
  partial r (the LKJ vine) the partial correlations are independent, so that
  holds; under flat R-squared they are tied together through the total
  $R^2$, so what the other coefficients explain moves the posterior of x1's
  partial correlation, and the gap grows with $k$. The joint test is further
  still (mean absolute 1.0, up to 2.4).
- **Flat R**: not defined; its implied prior for one partial correlation has
  an infinite density at 0 for $k \ge 2$.

Each prior is also calibrated on itself: with data drawn from the joint prior
(half with x1 = 0), the share of nulls among datasets with Bayes factor BF₁₀ is
1/(1 + BF₁₀) within Monte Carlo error, at k = 3 and 8, n = 50 and 200, for all
three (`validation/10`, P1). This justifies calling a slope of 0.1 null in
about half the datasets under flat partial r (study 06): under a flat partial
correlation an effect of 0.1 is rare, and at these n it predicts data much like
no effect, so among datasets with BF₀₁ > 3, 83–89% are true nulls (78–86% for
flat R-squared, 78% for flat R). Calibration is relative to the prior: with
data from the JZS prior, which expects many tiny effects, the same Bayes
factors are miscalibrated (P2).

Use the conditional test to ask "is there evidence for this coefficient?", and
the joint Bayes factors to compare or average over whole models.

---

## 11. The criteria for a default prior (Bayarri et al., 2012)

Bayarri, Berger, Forte and García-Donato (2012) list seven criteria a default
prior for model choice should meet. Status for each of the three priors;
"checked" means a numerical check in `validation/02_bayes_factor.R` (which
runs every check for each prior), "proof" and "sketch" refer to
[`3_consistency.md`](3_consistency.md).

| # | Criterion | Meaning | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|---|---|
| 1 | Basic | the prior on the tested coefficients is proper, and not arbitrarily vague | **met**, by construction | **met** | **met** |
| 2 | Model selection consistency | the true model's posterior probability goes to 1 as $n$ grows | **sketch**, checked (B3, B3b, B4): null rate $\tfrac q2 \log n$ in every case | **sketch**, checked: ½ log n against the intercept-only model, $\tfrac q2 \log n$ otherwise | **sketch**: $\log n$ against the intercept-only model, $\tfrac q2 \log n$ otherwise; checked (B3, B3b, B4): slopes 0.99, 1.00, 1.01 and 0.53, 1.64, 2.88 |
| 3 | Information consistency | BF → ∞ as $R^2 \to 1$ | **met for $n \ge k + 3$**: proof, checked (B2, B2b) | **met for $n \ge k + 3$**: proof, checked | **met for $n \ge k + 3$**: proof (through $f$ for $k = 1$), checked (B2, B2b) |
| 4 | Intrinsic consistency | the prior converges to a fixed proper prior as $n$ grows | **met**: proof | **met**: proof | **met**: proof |
| 5 | Predictive matching | on minimal samples ($n = k + 1$) BF = 1 with the null and between models of equal size | **met**: checked (B1) | **met**: checked (B1) | **met**: by the general argument, checked against the ground truth for $k \ge 2$ (B1) |
| 6 | Measurement invariance | units of the data do not change the answer | **met**, by construction, checked (A2) | **met**, checked (A2) | **met**, checked (A2) |
| 7 | Group invariance | the common parameters get the right-Haar prior, and the prior respects the invariance | **met** | **met** | **met** |

Dimensional matching follows from predictive matching for all three: the Bayes
factor depends only on $R^2$, $n$ and $k$, and $R^2 = 1$ for every model at
$n = k + 1$. Group invariance: $p(\alpha, \sigma) \propto 1/\sigma$ and $\beta$ is
scaled by $\sigma$ in all three.

**Information consistency, the argument.** As $R^2 \to 1$ the g-prior
Bayes factor of section 8 tends to $\int (1+g)^{(n-1-k)/2}\,p(g)\,dg$, so the
Bayes factor grows without limit exactly when that integral diverges (Liang et
al., 2008). For large $g$ every member of the family has
$p(g/n = t) \asymp t^{-2}$ (for flat R the confluent hypergeometric factor
tends to 1, $p(t) \approx 1/(2kt^2)$; for $B \sim$ Beta($a$, $k/2-a$) the bounds
$0 < B \le 1$ give the same order). So the integrand behaves like
$g^{(n-1-k)/2 - 2}$, which diverges if and only if $(n-1-k)/2 - 2 \ge -1$, that
is $n \ge k + 3$. The same tail on the effect scale: $p(f) \approx 2a f^{-3}$ for
every $a$ (section 2). Check B2b confirms the boundary: at $n = k + 2$ the
evidence levels off (for flat R, $\log\mathrm{BF} \to 0.91$ for $k = 1$), at
$n = k + 3$ it keeps growing. For comparison, the robust prior of Bayarri et
al. needs $n \ge k + 2$.

What is still open for criterion 2: [`0_overview.md`](0_overview.md), "Open".

---

## 12. Assumptions and limits

- Normal errors, fixed effects, an intercept.
- The prior uses the **sample** covariance of the predictors (it conditions
  on $X$), like the g-prior and JZS.
- The improper priors on $\alpha$ and $\sigma$ are only valid for comparing models
  that share them, which all models here do.
- Not covered: group-level effects, other likelihoods, distributional models.
