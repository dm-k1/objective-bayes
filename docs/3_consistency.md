# 3. Consistency: arguments behind the default-prior criteria

Section 11 of [`2_how_it_calculates.md`](2_how_it_calculates.md) lists the seven
criteria of Bayarri et al. (2012). This document gives the arguments for the
three that need one, for the whole family $R^2 \sim$ Beta($a$, 1) with a
uniform direction: flat partial r ($a = k/2$), flat R ($a = ½$) and flat
R-squared ($a = 1$). Each argument is labelled **proof** (complete, short) or
**sketch** (the steps are right, but a careful write-up with all conditions is
still to do). Numerical checks back each one up (`validation/02`, which runs
every check for each prior).

**In one sentence:** all three meet the same criteria by the same arguments;
they differ only in how fast evidence for the null model grows, $a \log n$
against the intercept-only model, because that rate is set by the prior's mass
near zero, $P(R^2 < 1/n) \approx n^{-a}$.

Throughout, for $a < k/2$ the prior is the mixture of Zellner g-priors with

$$
t = \frac{g}{n} = B \cdot W, \qquad B \sim \mathrm{Beta}\!\left(a, \tfrac k2 - a\right),
\qquad W \sim \text{Inverse-Gamma}\!\left(1, \tfrac12\right)
$$

(section 8 of document 2), and for $a = k/2$ it is the same with $B \equiv 1$.
That covers flat R (every $k$), flat partial r (every $k$) and flat R-squared
for $k \ge 3$; at $k = 2$ flat R-squared *is* flat partial r ($a = k/2 = 1$),
and at $k = 1$ it is not a g-mixture, so the arguments there go through the
effect $f$ directly. The density of $t$ is

$$
p(t) = \frac{1}{2t^2}\, E\!\left[B\, e^{-B/(2t)}\right]
$$

(take the expectation over $B$ of the Inverse-Gamma density of $W = t/B$, times
the change of variables $1/B$). Two bounds on it do the work below:

- **large $t$:** since $0 < B \le 1$ and $E[B] = 2a/k$,
  $\;\dfrac{e^{-1/(2t)}\,E[B]}{2t^2} \;\le\; p(t) \;\le\; \dfrac{E[B]}{2t^2}$, the
  same order $t^{-2}$ for every member;
- **small $t$:** only small $B$ contribute, where the Beta($a$, $k/2-a$) density
  behaves like $b^{a-1}/\mathrm{B}(a, \tfrac k2 - a)$, so
  $p(t) \approx c_{a,k}\, t^{a-1}$ as $t \to 0$ for $a < k/2$, with
  $c_{a,k} = 2^{a}\,\Gamma(a+1) / \mathrm{B}(a, \tfrac k2 - a)$. For flat R that is
  $t^{-1/2}$ (checked numerically: for $k = 4$, $\sqrt t\,p(t) = 0.7978$ at
  $t = 10^{-4}$, against $c_{1/2,4} = 0.7979$); for flat R-squared, $k \ge 3$, it
  is a constant. For $a = k/2$ ($B \equiv 1$) $p(t) = e^{-1/(2t)}/(2t^2)$ vanishes
  faster than any power as $t \to 0$.

The Bayes factor against the intercept-only model is (Liang et al., 2008)

$$
\mathrm{BF}_k = \int_0^\infty (1+nt)^{(n-1-k)/2}\,\bigl(1 + nt(1-R^2)\bigr)^{-(n-1)/2}\, p(t)\, dt .
$$

---

## Information consistency: proof

**Claim.** For every member of the family, $\mathrm{BF}_k \to \infty$ as
$R^2 \to 1$ if and only if $n \ge k + 3$.

*g-mixture members.* As $R^2 \to 1$ the integrand increases to
$(1+nt)^{(n-1-k)/2} p(t)$ (monotone convergence), so $\mathrm{BF}_k$ tends to
$\int_0^\infty (1+nt)^{(n-1-k)/2} p(t)\,dt$. Near 0 this integral is finite.
For large $t$, by the bounds above, the integrand lies between two constants
times $t^{(n-1-k)/2 - 2}$, so the integral diverges if and only if
$(n-1-k)/2 - 2 \ge -1$, that is $n \ge k + 3$.

*Every member, through $f$.* Write $\mathrm{BF}_k = \int_0^\infty L(f)\, p_a(f)\, df$
with $L$ the likelihood ratio of document 2, section 7 (step 5), and $p_a$ the
prior density of the length $f$ (section 2 there). For $f \ge 1$,

$$
\frac{p_a(f)}{p_{1/2}(f)} = 2a\left(\frac{f^2}{1+f^2}\right)^{a-1/2}
$$

lies between two positive constants ($2a\,2^{1/2-a}$ and $2a$ for $a \ge ½$). On
$f < 1$ the integral stays bounded as $R^2 \to 1$ ($L$ is continuous in $R^2$ and
$f$ on that compact range). So $\mathrm{BF}_k$ diverges for shape $a$ exactly
when it diverges for flat R, whose boundary is $n \ge k + 3$ by the g-mixture
argument. This covers flat R-squared at $k = 1$, the one case that is not a
g-mixture. ∎

The same tail in one line: $p_a(f) \approx 2a f^{-3}$ for large $f$, for every $a$.

Checked in `validation/02_bayes_factor.R` (B2b): at $n = k + 2$ the log Bayes
factor levels off (flat R: 0.91 for $k = 1$, 0.63 for $k = 3$), at $n = k + 3$
it keeps rising. Flat partial r and flat R-squared pass the same check
(run of 4 November 2026).

---

## Intrinsic consistency: proof

**Claim.** As $n \to \infty$, the prior for the coefficients converges to a fixed
proper prior, for every member.

The prior for $\beta$ given $\sigma$ depends on the data only through
$\Sigma_n = X_c^\top X_c / n$ (the shape $a$ and the hyperprior on $t = g/n$ do
not involve $n$). If $\Sigma_n \to \Sigma$, positive definite, as for predictors
sampled from a population, the density (section 5 of document 2) is continuous
in $\Sigma$ for every $\beta \ne 0$, so it converges pointwise to the density with
$\Sigma$, which is proper. By Scheffé's lemma the priors converge in total
variation. ∎

---

## Model selection consistency, fixed k: sketch

**Claim.** For nested models with a fixed number of predictors, the smaller
model $M_0$ with $k_0$ predictors inside the larger $M_1$ with $k = k_0 + q$:

- if $M_0$ is true and is the **intercept-only model** ($k_0 = 0$),
  $\log \mathrm{BF}_{01}$ grows like $a \log n$, with $a$ the shape of the larger
  model's prior: ½ log n for flat R and $\log n$ for flat R-squared, **the same
  for every $q$**; $\tfrac q2 \log n$ for flat partial r;
- if $M_0$ is true and **has a real effect** ($k_0 \ge 1$, $\rho_0^2 > 0$),
  $\log \mathrm{BF}_{01}$ grows like $\tfrac q2 \log n$, the usual rate for $q$
  extra parameters, for all three;
- if $M_1$ is true, $\log \mathrm{BF}_{10}$ grows linearly in $n$.

In every case the true model wins, so all three are consistent; the rates under
the intercept-only null differ, and that difference is worth knowing.

**Why the rates differ, in one picture.** The evidence for a model is roughly
its best fit times its *Occam factor*: the prior mass within the posterior's
reach, a neighbourhood of size $n^{-1/2}$ around the estimate.

- When the estimate is near $\beta = 0$ (the null is true), that neighbourhood
  is a ball of radius $\propto n^{-1/2}$ in $f$ around 0. Its prior mass is
  $P(f \lesssim n^{-1/2}) \approx P(R^2 \lesssim 1/n) = n^{-a}$:

  | | prior on $R^2$ | $P(R^2 \le x)$ | Occam factor at $x = 1/n$ | $\log \mathrm{BF}_{01}$ grows like |
  |---|---|---|---|---|
  | flat R | Beta(½, 1) | $x^{1/2}$ | $n^{-1/2}$ | $\tfrac12 \log n$ |
  | flat R-squared | Beta(1, 1) | $x$ | $n^{-1}$ | $\log n$ |
  | flat partial r | Beta($k/2$, 1) | $x^{k/2}$ | $n^{-k/2}$ | $\tfrac k2 \log n$ |

  For flat R the prior density of $\beta$ peaks at 0 (like $|\beta|^{1-k}$), so
  leaving the null model costs only $n^{-1/2}$ whatever $k$ is; flat R-squared's
  peak is lower ($|\beta|^{2-k}$), so it costs $n^{-1}$; flat partial r has no
  peak and pays the ordinary $n^{-k/2}$.
- When the estimate is away from 0 (some effect is real), the neighbourhood is
  a box of volume $\propto n^{-k/2}$ where the prior density is finite and
  positive, for all three. The Occam factor is $n^{-k/2}$, and adding $q$ null
  predictors costs $n^{-q/2}$.

*The smaller model is the intercept-only model and is true.* Then $nR^2$
stays bounded in probability ($\approx \chi^2_k$). Put $t = s/n$:

$$
\mathrm{BF}_k = \int_0^\infty (1+s)^{-k/2}
\underbrace{\left(\frac{1+s}{1+s(1-R^2)}\right)^{(n-1)/2}}_{\to\ \exp\left(\frac{nR^2}{2}\frac{s}{1+s}\right),\ \text{bounded}}
\; p\!\left(\tfrac{s}{n}\right) \frac{ds}{n} .
$$

For $a < k/2$, with $p(s/n) \approx c_{a,k} (s/n)^{a-1}$ for small $s/n$:

$$
\mathrm{BF}_k \approx n^{-a}\, c_{a,k} \int_0^\infty (1+s)^{-k/2} s^{a-1}\,
e^{\frac{nR^2}{2}\frac{s}{1+s}}\, ds ,
$$

and the integral is finite and bounded in probability because $a - 1 - k/2 < -1$.
So $\mathrm{BF}_{10} = O_p(n^{-a})$: the evidence for the true null grows like
$a \log n$ for every fixed $k$. The integral concentrates at $t \approx 1/n$,
where $p(t) \propto t^{a-1}$: that is the $n^{-a}$ of the picture above. For
$a = k/2$, $p(t)$ vanishes faster than any power near 0, the integral has no
such region, and the $(1+s)^{-k/2}$ factor sets the rate $\tfrac k2 \log n$.
For flat R-squared at $k = 1$ ($a > k/2$, not a g-mixture) the same rate,
$\log n$, follows from the Occam picture with $P(f^2 < 1/n) \approx n^{-1}$; it is
checked numerically, not proven here (study 02, B3b, $q = 1$).

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
over all $k$ dimensions. The argument needs only a prior density that is finite
and positive away from $\beta = 0$, which all three have (for flat R-squared at
$k = 1$ the same Laplace step applies to the integral over $f$).

Checked numerically (`validation/02`, B3b; 200 datasets per $n$ at
$n$ = 50, 200, 1000, 5000, slope of the median $\log \mathrm{BF}_{01}$ against
$\log n$; criterion: within 20% of $\tfrac q2$ when the smaller model is $x_1$,
within 30% of $a$ when it is the intercept-only model):

| $q$ null predictors added to | prior | expected | $q = 1$ | $q = 3$ | $q = 6$ |
|---|---|---|---|---|---|
| intercept only | flat R | ½, ½, ½ | 0.47 | 0.50 | 0.48 |
| intercept only | flat partial r | ½, 3/2, 3 | 0.51 | 1.46 | 2.86 |
| intercept only | flat R-squared | 1, 1, 1 | 0.99 | 1.00 | 1.01 |
| $x_1$, with $y = 0.5\,x_1 + \varepsilon$ | flat R | $\tfrac q2$ | 0.54 | 1.56 | 2.98 |
| $x_1$ | flat partial r | $\tfrac q2$ | 0.47 | 1.48 | 2.91 |
| $x_1$ | flat R-squared | $\tfrac q2$ | 0.53 | 1.64 | 2.88 |

Check B3 ($x_1, x_2$ real, $x_3$ null: $q = 1$) gives slope 0.53 for flat R and
0.49 for flat partial r and 0.51 for flat R-squared.

*The larger model is true.* Then $R^2_{\text{large}} \to \rho_1^2 > \rho_0^2 \leftarrow R^2_{\text{small}}$,
and for $t$ bounded away from 0 the ratio of the two integrands behaves like
$\bigl((1-R^2_{\text{small}})/(1-R^2_{\text{large}})\bigr)^{(n-1)/2}$, which grows
exponentially; the prior contributes only polynomial factors in $n$. So
$\log \mathrm{BF}_{10}$ grows linearly in $n$ (`validation/02`, B4, median
$\log\mathrm{BF}_{10}$ at $n$ = 200, 1000, 5000: flat R 0.9, 9.3, 54.1; flat
partial r −0.0, 8.5, 53.3; flat R-squared 1.1, 9.5, 54.3).

**What a full proof still needs:** the dominated-convergence steps written out
with explicit bounds, the non-nested case, and the null rate of flat R-squared
at $k = 1$.

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
| flat R, $k = n/4$ | 1.0 | 1.3 | 1.7 | 1.9 |
| flat R, $k = \sqrt n$ | 1.1 | 1.7 | 2.3 | 2.5 |
| flat R-squared, $k = n/4$ | 1.9 | 2.6 | 3.4 | 3.7 |
| flat R-squared, $k = \sqrt n$ | 2.1 | 3.3 | 4.4 | 5.0 |
| flat partial r, $k = n/4$ | 8.7 | 39.2 | 202.8 | 407.7 |
| flat partial r, $k = \sqrt n$ | 6.3 | 18.1 | 52.6 | 82.6 |

(Run of 4 November 2026; for comparison, JZS gives 3.8, 15.1, 74.8, 149.3 and
Zellner–Siow 10.1, 41.2, 205.6, 410.9 at $k = n/4$.) All three pass G1 and G2
(evidence keeps rising). Flat R stays consistent, but its evidence for the null barely moves:
with 500 null predictors at $n = 2000$, $\mathrm{BF}_{01} \approx 7$. Its
alternative stays at "$R$ flat", which predicts data close to the null's; flat
partial r and Zellner–Siow expect $R^2$ near 1 when $k$ is large, which the
data easily rule out. Flat R-squared's alternative stays at "$R^2$ flat", and
it sits between the two, much nearer flat R: 3.7 at $n = 2000$, $k = n/4$
($\mathrm{BF}_{01} \approx 40$). With a small true model (G2) flat R and flat
R-squared give nearly the same evidence for it (175.1 and 175.1 at
$n = 2000$), flat partial r more (443.6).

---

## The three priors side by side

| Criterion | Flat partial r | Flat R | Flat R-squared | Why |
|---|---|---|---|---|
| Information consistency | from $n \ge k + 3$ | from $n \ge k + 3$ | from $n \ge k + 3$ | the same $f^{-3}$ tail, $p(t) \asymp t^{-2}$ (proof above; `validation/02`, B2b) |
| Intrinsic consistency | yes | yes | yes | neither $a$ nor the hyperprior involves $n$; the density is continuous in $\Sigma$ |
| Predictive matching | yes | yes | yes | at $n = k + 1$ the data cannot tell the models apart, for any proper prior on the effect (B1) |
| Null rate, smaller model = intercept only | $\tfrac q2 \log n$ | $\tfrac12 \log n$ | $\log n$ | the prior mass near $\beta = 0$, $n^{-a}$ |
| Null rate, smaller model with an effect | $\tfrac q2 \log n$ | $\tfrac q2 \log n$ | $\tfrac q2 \log n$ | the Laplace argument above holds for any finite, positive prior density |
| Under an effect | linear in $n$ | linear in $n$ | linear in $n$ | (B4) |

So flat partial r gives the null model the ordinary penalty per dimension, like
Zellner–Siow and BIC; flat R gives the smallest penalty for leaving the null
model, whatever the number of predictors; flat R-squared sits between, with
one full $\log n$ whatever the number of predictors.
