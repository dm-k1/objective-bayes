# 1. Where the priors come from

With one predictor, a flat prior on the correlation is the natural objective
choice. With several, the set of valid correlation matrices forces a choice of
what stays flat, and the three possible answers form one family:
$R^2 \sim$ Beta($a$, 1), with every split of the total over the predictors
equally plausible.

This document builds that statement from what effects in a regression can
actually do. It starts from the project's first simulation (study 00) and
proceeds in steps. Steps 1, 2 and 4 each have a figure based on a simulation,
made by `docs/figures/make_figures.R`; step 0 uses the figure of study 00, and
step 3 is a definition without a figure. How the priors
become a density and a Bayes factor is described in
[`2_how_it_calculates.md`](2_how_it_calculates.md).

| Step | What effects can do                                                  | What the prior says                               | A judgement?                  |
| ---- | -------------------------------------------------------------------- | ------------------------------------------------- | ----------------------------- |
| 0    | valid correlation matrices can't keep everything flat as $k$ grows | a choice of what stays flat                       | yes: the fork, three branches |
| 1    | one correlation lies between −1 and 1                               | flat on the correlation                           | yes: which ruler is flat      |
| 2    | together, predictors explain at most 100%                            | a prior on the total, $R^2 \sim$ Beta($a$, 1) | yes: $a = k/2$, ½ or 1     |
| 3    | the total can be divided in many ways                                | every division equally plausible                  | yes: symmetry                 |
| 4    | —                                                                   | each predictor's share is Beta(½, (k−1)/2)      | no: follows from 2 and 3      |

---

## Step 0. The first result: not everything can stay flat

Study 00 (`motivating_simulation.R` in this folder) draws valid correlation
matrices of $(y, x_1 \ldots x_k)$ for $k = 1 \ldots 10$, 200,000 per $k$. In
each matrix it records three quantities: the **partial correlation** (= the
correlation of $y$ with $x_1$ after removing the other predictors from both), $R$, and
$R^2$.

1. **Simulate.** Version A draws every valid correlation matrix with equal
   probability (the vine method of Lewandowski, Kurowicka & Joe, 2009, with
   η = 1). Versions B and C keep each matrix's predictor correlations and the
   direction of $y$'s correlations with them, and only rescale the size of
   $y$'s correlations, so that $R$ is flat (B) or $R^2$ is flat (C).
   Rescaling keeps every matrix valid as long as $R < 1$.
2. **Look.** In version A each partial correlation is flat on (−1, 1) for
   every $k$, and $R^2$ moves towards 1. In versions B and C the total stays
   where it was placed, and the partial correlation narrows towards 0 as $k$
   grows:

   | $k$                     | 1    | 2    | 3    | 5    | 10   |
   | ------------------------- | ---- | ---- | ---- | ---- | ---- |
   | A: P(\|partial r\| > .30) | .70 | .70 | .70 | .70 | .70 |
   | A: P($R^2$ > .50)       | .29 | .50 | .65 | .82 | .97 |
   | B: P(\|partial r\| > .30) | .70 | .50 | .41 | .32 | .22 |
   | C: P(\|partial r\| > .30) | .91 | .70 | .60 | .50 | .36 |

   Version B keeps P($R^2$ > .50) at .29 and version C keeps it at .50, for
   every $k$.
3. **Check.** Each version is flat on the scale it should be: A's partial
   correlation, B's $R$ and C's $R^2$, compared with the flat distribution for
   every $k$ (checks F1–F3, largest Kolmogorov–Smirnov (KS) distance
   .0024–.0038).

The figure is `figures/motivation_partial_r_R_R2_by_k.png` (made
by running the study). At $k = 1$, versions A and B are the same prior,
Jeffreys' flat correlation; from $k = 2$ on they differ. Every later document
therefore asks which of the three branches to use, and the studies compare
them on that question.

---

## Step 1. One predictor: the correlation is bounded

With one predictor, the whole effect is the correlation $r$ between $x$ and
$y$, which can only lie between −1 and 1. Saying nothing beyond that, with
every value equally plausible, gives a flat prior on $r$: Jeffreys' default
for a correlation.

On the effect-size scale, the flat correlation becomes a heavy-tailed curve,
$\tfrac12(1+f^2)^{-3/2}$. Here $f = r/\sqrt{1-r^2}$ is the effect in noise
units (Cohen's $f$ with a sign). The curve resembles the Cauchy of the JZS
default, but its width was not chosen: it follows from the bounds of $r$. The
JZS default (grey in the figure) puts more weight on small effects, which is
its built-in expectation.

![](figures/1_one_predictor.png)

**The judgement in this step** is which ruler is flat. A flat prior on $r$ and
a flat prior on $r^2$ are different priors: flat on $r^2$ (version C at
$k = 1$) gives $r$ the density $|r|$, so correlations near zero get little
weight. The correlation is the natural scale of an effect, because it is
bounded, signed and comparable across studies. For that reason flat partial r
and flat R start from Jeffreys' prior and flat R-squared does not. In the
stretched-beta family of Ly, Verhagen and Wagenmakers, the flat prior on $r$
is shape $a = 1$, and a Cauchy on $f$ with scale 1 corresponds to
$a = \tfrac12$. The JZS default's scale of 0.354 is narrower and is not a
member of that family.

---

## Step 2. Several predictors: the bound is on the total

With several predictors, no single coefficient is bounded, but together the
predictors can explain at most all of $y$'s variance: the **multiple
correlation** $R$ (= the correlation between $y$ and its best linear
prediction from all predictors) lies between 0 and 1. Every valid
correlation matrix respects that bound; a separate prior for each coefficient
does not.

A prior per predictor ignores the shared bound. Each predictor brings its own
expected contribution, so the explained variance that the prior claims before
seeing data grows with every predictor added:

![](figures/2_drift.png)

| Prior                                                              | median implied $R^2$, $k = 1$ | $k = 10$ | $k = 30$ |
| ------------------------------------------------------------------ | --------------------------------- | ---------- | ---------- |
| JZS default (Cauchy, scale 0.354)                                  | .11                              | .72       | .89       |
| independent flat partial correlation per predictor                 | .25                              | .997      | 1.00       |
| **flat partial r** (every correlation matrix equally likely) | .25                              | .87       | .95       |
| **flat R**                                                   | .25                              | .25       | .25       |
| **flat R-squared**                                           | .50                              | .50       | .50       |

The three branches of step 0 are the three ways to respect the bound and stay
objective:

- **Flat on the correlations** (flat partial r): every valid correlation
  matrix of $(y, x_1 \ldots x_k)$ is equally likely. This implies
  $R^2 \sim$ Beta($k/2$, 1) (how that was found is described in
  [`relation_to_others.md`](relation_to_others.md)), so the claimed explained
  variance grows with $k$. It grows much less than under independent partial
  correlations (the second row of the table), which do not form one prior on
  the correlation matrix. It is the $R^2$ prior of rstanarm's `stan_lm` with
  η = 1.
- **Flat on the total $R$** (flat R): $R$ is flat from 0 to 1 for any number
  of predictors, so $R^2 \sim$ Beta(½, 1). With one predictor $R = |r|$, so
  this is exactly step 1.
- **Flat on the total $R^2$** (flat R-squared): $R^2$ is flat from 0 to 1 for
  any number of predictors, so $R^2 \sim$ Beta(1, 1).

**The judgement in this step** is which of the three is flat. The first says
that each coefficient keeps its Jeffreys prior, whatever else is in the model.
The other two say that the total is unknown whatever the number of
predictors, on the $R$ ruler or the $R^2$ ruler. Which of these describes real
studies better is partly an empirical question: how explained variance grows
with the number of predictors in practice.

---

## Step 3. Dividing the total: every split equally plausible

A given total can be divided over the predictors in many ways: with two
predictors, any point on a circle; with three, any point on a sphere. Without
knowledge about the predictors, every division is equally plausible. All
three branches make this assumption and differ only in the total.

"Equally plausible" is measured in explained variance, so the prior itself
handles correlated predictors: two predictors that overlap cannot both claim
the shared part.

**The judgement in this step** is symmetry. It is also the only split that
does not depend on the order or coding of the predictors. A flat Dirichlet on
the shares, for example, gives a different prior when correlated predictors
are listed in another order. For flat partial r it is not a separate choice,
because "every correlation matrix equally likely" already implies a uniform
direction.

---

## Step 4. What follows: the share per predictor

Steps 2 and 3 fix the total and spread it evenly. What one predictor gets is
then a consequence, not an extra assumption: a single predictor gets the
whole total, and with others it gets a random part of it, on average $1/k$.

![](figures/3_share_per_predictor.png)

The size of that part follows a Beta(½, (k−1)/2) distribution for every
branch. This was found by simulating, then deriving, then checking:

1. **Simulate.** A direction with no preference consists of $k$ independent
   normal draws $z_1 \ldots z_k$, divided by their length (in
   explained-variance units, so correlated predictors are already accounted
   for). Predictor 1's share of the total is
   $z_1^2 / (z_1^2 + \dots + z_k^2)$.
2. **Derive.** The numerator is one squared normal ($\chi^2_1$) and the rest
   of the denominator is the sum of $k - 1$ of them ($\chi^2_{k-1}$). A
   $\chi^2_a$ divided by itself plus an independent $\chi^2_b$ follows a
   Beta($a/2$, $b/2$) distribution. So

   $$
   \text{share} \sim \mathrm{Beta}\!\left(\tfrac12, \tfrac{k-1}{2}\right),
   \qquad E[\text{share}] = \frac{1/2}{1/2 + (k-1)/2} = \frac1k .
   $$
3. **Check** (printed by `docs/figures/make_figures.R`). With 40,000
   directions per $k$, the KS distance to that Beta distribution is
   .003–.006, and the mean share is .499, .330, .200 and .101 for $k$ = 2, 3,
   5 and 10. For $k = 1$ the share is 1: the predictor gets the whole total.

Predictor 1's own ratio of explained to unexplained variance (its Cohen's
$f^2$) is then share × $R^2/(1-R^2)$, with $R^2 \sim$ Beta($a$, 1). This is
the prior of the conditional test
([`2_how_it_calculates.md`](2_how_it_calculates.md), section 10).

- **Flat partial r** ($a = k/2$): the mixture is exactly a flat partial
  correlation, for every $k$. The total grows with $k$ just fast enough to
  keep each coefficient's partial correlation flat (validation 05, C4).
- **Flat R and flat R-squared** ($a$ fixed): the total stays where it is and
  the share shrinks, so each coefficient's prior narrows as $k$ grows (study
  00, versions B and C).

---

## What each prior does expect

The motto "define what's possible, not what you expect" applies to whatever
is kept flat. Each branch still carries expectations, and they are stated here
rather than left for the reader to find:

|                                                                   | Flat partial r                                 | Flat R              | Flat R-squared     |
| ----------------------------------------------------------------- | ---------------------------------------------- | ------------------- | ------------------ |
| $E[R^2]$                                                        | $k/(k+2)$: 1/3, 3/5, 5/6 at $k$ = 1, 3, 10 | 1/3                 | 1/2                |
| median $R^2$                                                    | $0.5^{2/k}$: .25, .63, .87                | .25                | .50               |
| expected explained variance per predictor, $E[\text{share} \times R^2]$ | $1/(k+2)$                                    | $1/(3k)$          | $1/(2k)$         |
| one coefficient's prior                                           | flat partial correlation, any $k$            | narrows with $k$  | narrows with $k$ |
| evidence for the null against the intercept-only model grows like | $\tfrac k2 \log n$                           | $\tfrac12 \log n$ | $\log n$         |

- Under flat R and flat R-squared, adding predictors does not add explained
  variance to the prior. This acts as a built-in correction for testing many
  predictors, and as shrinkage towards zero at small $n$.
- Under flat partial r, each coefficient keeps the same prior whatever the
  model around it, and a large model claims a large total.
- The amount of prior mass near zero sets the rate in the last row:
  $P(R^2 < 1/n) \approx n^{-a}$ for $R^2 \sim$ Beta($a$, 1). Study 02 measures
  this rate as the slope of the evidence against $\log n$ (check B3b).

No branch has every property, and the trade-offs are measured in
[`3_comparison.md`](3_comparison.md), which ends with the recommendation (flat
partial r) and its costs.
