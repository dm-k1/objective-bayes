# 3. Consistency: arguments behind the default-prior criteria

Section 10 of [`2_how_it_calculates.md`](2_how_it_calculates.md) lists the seven
criteria of Bayarri et al. (2012). This document gives the arguments for the
three that need one. Each is labelled **proof** (complete, short) or **sketch**
(the steps are right, but a careful write-up with all conditions is still to do).
Numerical checks back each one up.

Throughout: the flat-R prior is the mixture of Zellner g-priors with

$$
t = \frac{g}{n} = B \cdot W, \qquad B \sim \mathrm{Beta}\!\left(\tfrac12, \tfrac{k-1}{2}\right),
\qquad W \sim \text{Inverse-Gamma}\!\left(1, \tfrac12\right)
$$

(section 8 of document 2), so the density of $t$ is

$$
p(t) = \frac{1}{2t^2}\, E\!\left[B\, e^{-B/(2t)}\right]
$$

(take the expectation over $B$ of the Inverse-Gamma density of $W = t/B$, times
the change of variables $1/B$). Two bounds on it do all the work below:

- **large $t$:** since $0 < B \le 1$ and $E[B] = 1/k$,
  $\;\dfrac{e^{-1/(2t)}}{2k\,t^2} \;\le\; p(t) \;\le\; \dfrac{1}{2k\,t^2}$;
- **small $t$:** only small $B$ contribute, where the Beta density behaves like
  $b^{-1/2}/\mathrm{B}(\tfrac12, \tfrac{k-1}{2})$, so
  $p(t) \approx c_k\, t^{-1/2}$ as $t \to 0$ for $k \ge 2$, with
  $c_k = \sqrt2\,\Gamma(\tfrac32) / \mathrm{B}(\tfrac12, \tfrac{k-1}{2})$ (checked numerically: for $k = 4$, $\sqrt t\,p(t) = 0.7978$ at $t = 10^{-4}$, against $c_4 = 0.7979$).
  (For $k = 1$, $B = 1$ and $p(t)$ vanishes faster than any power as $t \to 0$.)

The Bayes factor against the intercept-only model is (Liang et al., 2008)

$$
\mathrm{BF}_k = \int_0^\infty (1+nt)^{(n-1-k)/2}\,\bigl(1 + nt(1-R^2)\bigr)^{-(n-1)/2}\, p(t)\, dt .
$$

---

## Information consistency: proof

**Claim.** $\mathrm{BF}_k \to \infty$ as $R^2 \to 1$ if and only if $n \ge k + 3$.

As $R^2 \to 1$ the integrand increases to $(1+nt)^{(n-1-k)/2} p(t)$ (monotone
convergence), so $\mathrm{BF}_k$ tends to $\int_0^\infty (1+nt)^{(n-1-k)/2} p(t)\,dt$.
Near 0 this integral is finite. For large $t$, by the bounds above, the
integrand lies between two constants times $t^{(n-1-k)/2 - 2}$, so the integral
diverges if and only if $(n-1-k)/2 - 2 \ge -1$, that is $n \ge k + 3$. ∎

Checked in `validation/02_bayes_factor.R` (B2b): at $n = k + 2$ the log Bayes
factor levels off (0.91 for $k = 1$, 0.63 for $k = 3$), at $n = k + 3$ it keeps rising.

---

## Intrinsic consistency: proof

**Claim.** As $n \to \infty$, the prior for the coefficients converges to a fixed
proper prior.

The prior for $\beta$ given $\sigma$ depends on the data only through
$\Sigma_n = X_c^\top X_c / n$ (the hyperprior on $t = g/n$ does not involve $n$).
If $\Sigma_n \to \Sigma$, positive definite, as for predictors sampled from a
population, the flat-R density (section 5 of document 2) is continuous in
$\Sigma$ for every $\beta \ne 0$, so it converges pointwise to the flat-R density
with $\Sigma$, which is proper. By Scheffé's lemma the priors converge in total
variation. ∎

---

## Model selection consistency, fixed k: sketch

**Claim.** For nested models with a fixed number of predictors, the smaller
model $M_0$ with $k_0$ predictors inside the larger $M_1$ with $k = k_0 + q$:

- if $M_0$ is true and is the **intercept-only model** ($k_0 = 0$),
  $\log \mathrm{BF}_{01}$ grows like $\tfrac12 \log n$, **the same rate for every $q$**;
- if $M_0$ is true and **has a real effect** ($k_0 \ge 1$, $\rho_0^2 > 0$),
  $\log \mathrm{BF}_{01}$ grows like $\tfrac q2 \log n$, the usual rate for $q$
  extra parameters;
- if $M_1$ is true, $\log \mathrm{BF}_{10}$ grows linearly in $n$.

In both null cases the true model wins, so the prior is consistent; the two
rates differ, and that difference is a property of flat-R worth knowing.

**Why the rates differ, in one picture.** The evidence for a model is roughly
its best fit times its *Occam factor*: the prior mass within the posterior's
reach, a neighbourhood of size $n^{-1/2}$ around the estimate.

- When the estimate is near $\beta = 0$ (the null is true), that neighbourhood
  is a ball of radius $\propto n^{-1/2}$ around 0. Flat on $R$ gives it prior
  mass $P(R \lesssim n^{-1/2}) \propto n^{-1/2}$, for every $k$: the prior density
  of $\beta$ peaks at 0 (like $|\beta|^{-(k-1)}$), so leaving the null model
  costs only $n^{-1/2}$.
- When the estimate is away from 0 (some effect is real), the neighbourhood is a
  box of volume $\propto n^{-k/2}$ where the prior density is finite and
  positive (flat in $R$, uniform in direction). The Occam factor is
  $n^{-k/2}$, and adding $q$ null predictors costs $n^{-q/2}$.

So flat-R has no extra penalty for leaving the null model, but the ordinary
penalty for each extra dimension once some effect is in the model.

*The smaller model is the intercept-only model and is true.* Then $nR^2$
stays bounded in probability ($\approx \chi^2_k$). Put $t = s/n$:

$$
\mathrm{BF}_k = \int_0^\infty (1+s)^{-k/2}
\underbrace{\left(\frac{1+s}{1+s(1-R^2)}\right)^{(n-1)/2}}_{\to\ \exp\left(\frac{nR^2}{2}\frac{s}{1+s}\right),\ \text{bounded}}
\; p\!\left(\tfrac{s}{n}\right) \frac{ds}{n} .
$$

With $p(s/n) \approx c_k (s/n)^{-1/2}$ for small $s/n$:

$$
\mathrm{BF}_k \approx n^{-1/2}\, c_k \int_0^\infty (1+s)^{-k/2} s^{-1/2}\,
e^{\frac{nR^2}{2}\frac{s}{1+s}}\, ds ,
$$

and the integral is finite and bounded in probability for $k \ge 2$
(for $k = 1$ the same rate follows from $E[t^{-1/2}] < \infty$). So
$\mathrm{BF}_{10} = O_p(n^{-1/2})$: the evidence for the true null grows like
$\tfrac12 \log n$ for every fixed $k$. The integral concentrates at
$t \approx 1/n$, where $p(t) \propto t^{-1/2}$: that is the $n^{-1/2}$ of the
picture above.

*The smaller model has a real effect and is true.* Now
$R^2 \to \rho_0^2 > 0$ under both models, and the integrand of $\mathrm{BF}_k$
sits at $t$ of order 1 instead of $1/n$ (a Laplace argument). For fixed $t > 0$,

$$
(1+nt)^{(n-1-k)/2}\bigl(1+nt(1-R^2)\bigr)^{-(n-1)/2}
\approx (nt)^{-k/2}\,(1-R^2)^{-(n-1)/2}\,\exp\!\left(-\frac{R^2}{2(1-R^2)\,t}\right),
$$

so

$$
\log \mathrm{BF}_k = -\tfrac{n-1}{2}\log(1-R^2) - \tfrac k2 \log n
+ \log \int_0^\infty t^{-k/2} e^{-R^2/(2(1-R^2)t)}\, p(t)\, dt + o_p(1),
$$

and the last integral is a finite positive constant (the exponential removes
small $t$, and $p(t) \propto t^{-2}$ for large $t$). Taking the difference between
the two models,

$$
\log \mathrm{BF}_{01} = \tfrac q2 \log n - \tfrac{n-1}{2}\log\frac{1-R^2_0}{1-R^2_1} + O_p(1),
$$

and the middle term tends to $\tfrac12\chi^2_q$ (the likelihood-ratio statistic),
so it is bounded in probability. The $(nt)^{-k/2}$ is the $n^{-k/2}$ box of the
picture: it is the average over a uniform direction, which spreads the prior
over all $k$ dimensions. Under the null this factor stays of order 1
(because $s = nt$ stays of order 1), which is why the rate there does not depend
on $k$.

Checked numerically (`validation/02`, B3b; 200 datasets per $n$ at
$n$ = 50, 200, 1000, 5000, slope of the median $\log \mathrm{BF}_{01}$ against
$\log n$):

| $q$ null predictors added to | $q = 1$ | $q = 3$ | $q = 6$ |
|---|---|---|---|
| intercept only (expected $\tfrac12$) | 0.51 | 0.51 | 0.51 |
| $x_1$, with $y = 0.5\,x_1 + \varepsilon$ (expected $\tfrac q2$) | 0.55 | 1.59 | 2.91 |

Check B3 ($x_1, x_2$ real, $x_3$ null: $q = 1$) gives slope 0.53.

*The larger model is true.* Then $R^2_{\text{large}} \to \rho_1^2 > \rho_0^2 \leftarrow R^2_{\text{small}}$,
and for $t$ bounded away from 0 the ratio of the two integrands behaves like
$\bigl((1-R^2_{\text{small}})/(1-R^2_{\text{large}})\bigr)^{(n-1)/2}$, which grows
exponentially; the prior contributes only polynomial factors in $n$. So
$\log \mathrm{BF}_{10}$ grows linearly in $n$ (`validation/02`, B4: 0.9, 9.3, 54.1
at $n$ = 200, 1000, 5000).

**What a full proof still needs:** the dominated-convergence steps written out
with explicit bounds, and the non-nested case.

---

## When k grows with n: numerical only

The sketch above also shows why this is harder. Under the null,
$nR^2 \approx \chi^2_k$ grows with $k$, so the factor
$\exp\bigl(\tfrac{nR^2}{2}\tfrac{s}{1+s}\bigr)$ is no longer bounded: the evidence
for the null still grows, but more slowly. `validation/08_growing_k.R` checks
this for $k = n/4$ and $k = \sqrt n$ (200 null datasets per $n$, median
$\log\mathrm{BF}_{01}$ of the null model against the full model):

| | $n$ = 50 | 200 | 1000 | 2000 |
|---|---|---|---|---|
| flat-R, $k = n/4$ | 1.0 | 1.3 | 1.7 | 1.9 |
| flat-R, $k = \sqrt n$ | 1.1 | 1.7 | 2.3 | 2.5 |
| uniform correlations, $k = n/4$ | 8.8 | 39.2 | 202.8 | 407.7 |
| uniform correlations, $k = \sqrt n$ | 6.3 | 18.1 | 52.6 | 82.6 |

Flat-R stays consistent, but its evidence for the null barely moves: with 500
null predictors at $n = 2000$, $\mathrm{BF}_{01} \approx 7$. Its alternative stays at
"$R$ flat", which predicts data close to the null's; uniform correlations and
Zellner–Siow expect $R^2$ near 1 when $k$ is large, which the data easily rule
out. That is the shared-total property, seen from the null side.

---|---|---|---|---|---|
| $k$ fixed at 5 | 1.19 | 1.86 | 2.59 | 3.43 | 4.12 |
| $k = \sqrt n$ | 1.18 | 1.64 | 2.25 | 2.98 | 3.49 |
| $k = n/4$ | 0.99 | 1.34 | 1.80 | 2.15 | 2.49 |

Consistent, but with $k = n/4$ at about half the rate of fixed $k$. Zellner–Siow
accumulates evidence for the null much faster here, because its prior expects
$R^2$ near 1 when $k$ is large, which the data easily rule out; flat-R's
alternative stays at "$R$ flat", which is harder to tell apart from the null.
That is the same shared-total property, seen from the null side.

---

## The uniform-correlations prior: the same criteria

**In one sentence:** uniform correlations meets the same criteria by the same
arguments, with one difference: leaving the null model costs the usual
$\tfrac k2 \log n$, because its prior has no peak at $\beta = 0$.

As a g-mixture its hyperprior is the same for every $k$ (document 2, section 8,
$B \equiv 1$):

$$
t = \frac gn \sim \text{Inverse-Gamma}\!\left(1, \tfrac12\right), \qquad
p(t) = \frac{1}{2t^2}\, e^{-1/(2t)} .
$$

| Criterion | Flat-R | Uniform correlations | Why |
|---|---|---|---|
| Information consistency | from $n \ge k + 3$ | from $n \ge k + 3$ | the same $t^{-2}$ tail; the proof above with $p(t) \le 1/(2t^2)$ (`validation/02`, B2b and CB2) |
| Intrinsic consistency | yes | yes | the hyperprior does not involve $n$; the prior density is continuous in $\Sigma$ |
| Predictive matching | yes | yes | every g-mixture gives BF = 1 at $n = k + 1$ (B1, CB1) |
| Null rate, smaller model = intercept only | $\tfrac12 \log n$ | $\tfrac q2 \log n$ | the prior mass near $\beta = 0$, below |
| Null rate, smaller model with an effect | $\tfrac q2 \log n$ | $\tfrac q2 \log n$ | the Laplace argument above holds for any finite, positive prior density |
| Under an effect | linear in $n$ | linear in $n$ | (B4, CB4) |

**Why the null rate differs.** The Occam factor of the picture above is the
prior mass within $n^{-1/2}$ of $\beta = 0$, which is $P(R \lesssim n^{-1/2})$:

| | prior on $R$ | $P(R \le a)$ | Occam factor at $a = n^{-1/2}$ | $\log \mathrm{BF}_{01}$ grows like |
|---|---|---|---|---|
| flat-R | flat | $a$ | $n^{-1/2}$ | $\tfrac12 \log n$ |
| uniform correlations | density $k R^{k-1}$ ($R^2 \sim$ Beta($k/2$, 1)) | $a^k$ | $n^{-k/2}$ | $\tfrac k2 \log n$ |

So under uniform correlations the null model gets the ordinary penalty per
dimension, like Zellner–Siow and BIC. In the g-mixture form: $p(t)$ vanishes
faster than any power as $t \to 0$, so the integral of the null case above has
no $t^{-1/2}$ region near $t \approx 1/n$, and the $(1+s)^{-k/2}$ factor sets the rate.

Checked (`validation/02`, CB3; slope of the median $\log \mathrm{BF}_{01}$
against $\log n$, expected $\tfrac q2$ in every cell):

| $q$ null predictors added to | $q = 1$ | $q = 3$ | $q = 6$ |
|---|---|---|---|
| intercept only | 0.50 | 1.46 | 2.87 |
| $x_1$, with $y = 0.5\,x_1 + \varepsilon$ | 0.50 | 1.47 | 2.89 |

With $k$ growing with $n$ it accumulates evidence for the null at the
Zellner–Siow rate (the table in "When k grows with n" above).
