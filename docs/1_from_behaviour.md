# 1. Where the priors come from

**In one sentence:** with one predictor a flat correlation is the only
objective choice; with several, the valid correlation matrices force a choice
of what stays flat, and the three answers are one family, $R^2 \sim$
Beta($a$, 1) with every split of the total equally plausible.

This document builds that sentence from what effects in a regression can
actually do. It starts from the first simulation of the project (study 00),
then goes step by step. Each step has a figure that is a simulation, made by
`docs/figures/make_figures.R` or by the validation scripts. How the priors
become a density and a Bayes factor is in
[`2_how_it_calculates.md`](2_how_it_calculates.md).

| Step | What effects can do | What the prior says | A judgement? |
|---|---|---|---|
| 0 | valid correlation matrices can't keep everything flat as $k$ grows | a choice of what stays flat | yes: the fork, three branches |
| 1 | one correlation lies between −1 and 1 | flat on the correlation | yes: which ruler is flat |
| 2 | together, predictors explain at most 100% | a prior on the total, $R^2 \sim$ Beta($a$, 1) | yes: $a = k/2$, ½ or 1 |
| 3 | the total can be divided in many ways | every division equally plausible | yes: symmetry |
| 4 | — | each predictor's share is Beta(½, (k−1)/2) | no: follows from 2 and 3 |

---

## Step 0. The first result: not everything can stay flat

Study 00 (`validation/00_initial_simulation.R`) draws valid correlation
matrices of $(y, x_1 \ldots x_k)$ for $k = 1 \ldots 10$, 200,000 per $k$, and
looks at three things in each: the partial correlation of $y$ with $x_1$ given
the other predictors, $R$, and $R^2$.

1. **Simulate.** Draw every valid correlation matrix with equal probability
   (the vine method of Lewandowski, Kurowicka & Joe, 2009, with η = 1). That
   is version A. For versions B and C, keep each matrix's predictor
   correlations and the direction of $y$'s correlations with them, and only
   rescale how large $y$'s correlations are, so that $R$ is flat (B) or $R^2$
   is flat (C). Rescaling keeps every matrix valid as long as $R < 1$.
2. **Look.** In A each partial correlation is flat on (−1, 1) for every $k$, and
   $R^2$ moves toward 1. In B and C the total stays where it was put, and the
   partial correlation narrows toward 0 as $k$ grows:

   | $k$ | 1 | 2 | 3 | 5 | 10 |
   |---|---|---|---|---|---|
   | A: P(\|partial r\| > 0.3) | 0.70 | 0.70 | 0.70 | 0.70 | 0.70 |
   | A: P($R^2$ > 0.5) | 0.29 | 0.50 | 0.65 | 0.82 | 0.97 |
   | B: P(\|partial r\| > 0.3) | 0.70 | 0.50 | 0.41 | 0.32 | 0.22 |
   | C: P(\|partial r\| > 0.3) | 0.91 | 0.70 | 0.60 | 0.49 | 0.36 |

   (B keeps P($R^2$ > 0.5) at 0.29 and C at 0.50, for every $k$.)
3. **Check.** Each version is flat where it should be: A's partial correlation,
   B's $R$ and C's $R^2$, against the flat distribution for every $k$ (checks
   F1–F3, KS distance 0.002–0.003).

The figure is `validation/results/figures/00_partial_r_R_R2_by_k.png` (made
by running the study). At $k = 1$, A and B are the same prior: Jeffreys' flat
correlation. From $k = 2$ on they part. So the question for every later
document is not "which prior is right" but **"which of the three branches"**;
the studies compare them.

---

## Step 1. One predictor: the correlation is bounded

With one predictor, the whole effect is the correlation $r$ between $x$ and
$y$, and a correlation can only lie between −1 and 1. Saying nothing more
than that, every value equally plausible, is a flat prior on $r$: Jeffreys'
default for a correlation.

On the effect-size scale (the effect in noise units, Cohen's $f$ with a sign,
$f = r/\sqrt{1-r^2}$) the flat correlation becomes a heavy-tailed curve,
$\tfrac12(1+f^2)^{-3/2}$. It resembles the JZS default's Cauchy, but its width
was not chosen: it follows from the bounds of $r$. The JZS default puts more
weight on small effects (grey): that is its built-in expectation.

![](figures/1_one_predictor.png)

**The judgement here** is which ruler is flat. Flat on $r$ and flat on $r^2$
are different priors: flat on $r^2$ (version C at $k = 1$) gives $r$ the
density $|r|$, so correlations near zero get little weight. The correlation is
the natural scale of an effect: bounded, signed, and comparable across
studies. That is why flat partial r and flat R start from Jeffreys, and flat
R-squared does not. (In the stretched-beta family of Ly, Verhagen &
Wagenmakers, flat on $r$ is shape $a = 1$; a Cauchy on $f$ with scale 1
corresponds to $a = \tfrac12$. The JZS default's scale 0.354 is not a member of
that family: it is narrower.)

---

## Step 2. Several predictors: the bound is on the total

With several predictors, no single coefficient has a bound, but together they
can explain at most all of $y$'s variance: the multiple correlation $R$ lies
between 0 and 1. Every valid correlation matrix respects that bound; a prior
per coefficient does not.

A prior per predictor ignores the shared bound. Each predictor brings its own
expected share, so the explained variance the prior claims before seeing data
grows with every predictor added:

![](figures/2_drift.png)

| Prior | median implied $R^2$, $k = 1$ | $k = 10$ | $k = 30$ |
|---|---|---|---|
| JZS default (Cauchy, scale 0.354) | 0.11 | 0.72 | 0.89 |
| independent flat partial correlation per predictor | 0.25 | 0.997 | 1.00 |
| **flat partial r** (every correlation matrix equally likely) | 0.25 | 0.87 | 0.95 |
| **flat R** | 0.25 | 0.25 | 0.25 |
| **flat R-squared** | 0.50 | 0.50 | 0.50 |

The three branches of step 0 are the three ways to keep the bound and stay
objective:

- **Flat on the correlations** (flat partial r): every valid correlation
  matrix of $(y, x_1 \ldots x_k)$ equally likely. It implies
  $R^2 \sim$ Beta($k/2$, 1) (how that was found:
  [`relation_to_others.md`](relation_to_others.md)), so the claimed explained
  variance grows with $k$, but much less than independent partial correlations
  would (the second row of the table: those are not one prior on the matrix).
  It is rstanarm `stan_lm`'s $R^2$ prior with η = 1.
- **Flat on the total $R$** (flat R): $R$ flat from 0 to 1 for any number of
  predictors; $R^2 \sim$ Beta(½, 1). With one predictor $R = |r|$, so this is
  exactly step 1.
- **Flat on the total $R^2$** (flat R-squared): $R^2$ flat from 0 to 1 for any
  number of predictors; $R^2 \sim$ Beta(1, 1).

**The judgement here** is which of the three is flat. The first says each
coefficient keeps its Jeffreys prior, whatever else is in the model; the other
two say the total is unknown whatever the number of predictors, on the $R$ or
the $R^2$ ruler. Which describes real studies better is partly an empirical
question: how explained variance grows with the number of predictors in
practice.

---

## Step 3. Dividing the total: every split equally plausible

A given total can be divided over the predictors in many ways: with two
predictors every point on a circle, with three every point on a sphere. Without
knowing anything about the predictors, every division is equally plausible.
All three branches do this; they differ only in the total.

"Equally plausible" is measured in explained variance, so correlated
predictors are handled by the prior itself: two predictors that overlap can't
both claim the shared part.

**The judgement here** is symmetry. It is also the only split that does not
depend on the order or coding of the predictors: a flat Dirichlet on the
shares, for example, gives a different prior when correlated predictors are
listed in another order. For flat partial r it is not a separate choice:
"every correlation matrix equally likely" already implies a uniform direction.

---

## Step 4. What follows: the share per predictor

Steps 2 and 3 fix the total and spread it evenly. What one predictor gets is
then a consequence, not an extra assumption: alone it gets everything; with
others it gets a random slice, on average $1/k$ of the total.

![](figures/3_share_per_predictor.png)

**How large is the slice? Beta(½, (k−1)/2), for every branch, found by
simulating, then deriving, then checking.**

1. **Simulate.** A direction with no preference is $k$ independent normal draws
   $z_1 \ldots z_k$, divided by their length (in explained-variance units, so
   correlated predictors are already accounted for). Predictor 1's share of the
   total is $z_1^2 / (z_1^2 + \dots + z_k^2)$.
2. **Derive.** The numerator is one squared normal ($\chi^2_1$) and the
   remainder is $k - 1$ of them ($\chi^2_{k-1}$). A $\chi^2_a$ divided by itself
   plus an independent $\chi^2_b$ is Beta($a/2$, $b/2$). So

   $$
   \text{share} \sim \mathrm{Beta}\!\left(\tfrac12, \tfrac{k-1}{2}\right),
   \qquad E[\text{share}] = \frac{1/2}{1/2 + (k-1)/2} = \frac1k .
   $$

3. **Check.** 40,000 directions per $k$: KS distance to that Beta 0.003–0.006,
   mean share 0.501, 0.333, 0.201, 0.100 for $k$ = 2, 3, 5, 10. For $k = 1$ the
   share is 1: the predictor gets the whole total.

Predictor 1's own explained/unexplained variance (its Cohen's $f^2$) is then
share × $R^2/(1-R^2)$ with $R^2 \sim$ Beta($a$, 1): the prior of the
conditional test ([`2_how_it_calculates.md`](2_how_it_calculates.md),
section 10).

- **Flat partial r** ($a = k/2$): the mixture is exactly a flat partial
  correlation, for every $k$. The total grows with $k$ just fast enough to keep
  each slice flat (validation 05, C4).
- **Flat R and flat R-squared** ($a$ fixed): the total stays put and the slice
  shrinks, so each coefficient's prior narrows as $k$ grows (study 00, versions
  B and C).

---

## What each prior does expect

The motto "define what's possible, not what you expect" holds for whatever is
flat. Each branch does carry expectations, and they are stated here rather
than left for a reader to find:

| | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|
| $E[R^2]$ | $k/(k+2)$: 1/3, 3/5, 5/6 at $k$ = 1, 3, 10 | 1/3 | 1/2 |
| median $R^2$ | $0.5^{2/k}$: 0.25, 0.63, 0.87 | 0.25 | 0.50 |
| expected share per predictor | $1/(k+2)$ | $1/(3k)$ | $1/(2k)$ |
| one coefficient's prior | flat partial correlation, any $k$ | narrows with $k$ | narrows with $k$ |
| evidence for the null against the intercept-only model grows like | $\tfrac k2 \log n$ | $\tfrac12 \log n$ | $\log n$ |

- For flat R and flat R-squared, adding predictors does not add explained
  variance to the prior: a built-in correction for testing many predictors, and
  shrinkage toward zero at small $n$.
- For flat partial r, each coefficient keeps the same prior whatever the model
  around it, and a large model claims a large total.
- How much prior mass sits near zero sets the null rate in the last row
  ([`3_consistency.md`](3_consistency.md)): $P(R^2 < 1/n) \approx n^{-a}$.

No branch gets every property for free; the trade-offs are measured in
[`4_comparison.md`](4_comparison.md), which ends with the recommendation (flat
partial r) and what it costs.
