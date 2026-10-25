# 1. Where the prior comes from

**In one sentence:** the total strength of all predictors together gets one
flat prior, the same for any number of predictors, and every way of dividing
that total over the predictors is equally plausible.

This document builds that sentence from what effects in a regression can
actually do, in four steps. Each step has a figure that is a simulation, made
by `docs/figures/make_figures.R`. How the sentence becomes a density and a
Bayes factor is in [`2_how_it_calculates.md`](2_how_it_calculates.md).

| Step | What effects can do | What the prior says | A judgement? |
|---|---|---|---|
| 1 | one correlation lies between −1 and 1 | flat on the correlation | yes: which ruler is flat |
| 2 | together, predictors explain at most 100% | the prior goes on the total, $R$ | yes: the total is what gets the prior |
| 3 | the total can be divided in many ways | every division equally plausible | yes: symmetry |
| 4 | — | each predictor's share shrinks as $k$ grows | no: follows from 2 and 3 |

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
weight on small effects (orange): that is its built-in expectation.

![](figures/1_one_predictor.png)

**The judgement here** is which ruler is flat. Flat on $r$ and flat on $r^2$
are different priors. The correlation is chosen because it is the natural
scale of an effect: bounded, signed, and comparable across studies. (In the
stretched-beta family of Ly, Verhagen & Wagenmakers, flat on $r$ is shape
$a = 1$; the JZS Cauchy corresponds to $a = \tfrac12$.)

---

## Step 2. Several predictors: the bound is on the total

With several predictors, no single coefficient has a bound, but together they
can explain at most all of $y$'s variance: the multiple correlation $R$ lies
between 0 and 1.

A prior per predictor ignores that shared bound. Each predictor brings its own
expected share, so the explained variance the prior claims before seeing data
grows with every predictor added:

![](figures/2_drift.png)

| Prior | median implied $R^2$, $k = 1$ | $k = 10$ | $k = 30$ |
|---|---|---|---|
| JZS default (Cauchy, scale 0.354) | 0.11 | 0.72 | 0.89 |
| flat partial correlation per predictor | 0.25 | 0.997 | 1.00 |
| **flat-R** (flat on the total) | 0.25 | 0.25 | 0.25 |

So the flat prior of step 1 goes on the total: **$R$ is flat from 0 to 1, for
any number of predictors.** With one predictor $R = |r|$, so this is exactly
step 1. On the $R^2$ scale it is $\mathrm{Beta}(\tfrac12, 1)$: the same
statement on another ruler.

---

## Step 3. Dividing the total: every split equally plausible

A given total can be divided over the predictors in many ways: with two
predictors every point on a circle, with three every point on a sphere. Without
knowing anything about the predictors, every division is equally plausible.

"Equally plausible" is measured in explained variance, so correlated
predictors are handled by the prior itself: two predictors that overlap can't
both claim the shared part.

**The judgement here** is symmetry. It is also the only split that does not
depend on the order or coding of the predictors: a flat Dirichlet on the
shares, for example, gives a different prior when correlated predictors are
listed in another order.

---

## Step 4. What follows: the share per predictor

Steps 2 and 3 fix the total and spread it evenly. What one predictor gets is
then a consequence, not an extra assumption: alone it gets everything (flat);
with others it gets a random slice, on average $1/k$ of the total.

![](figures/3_share_per_predictor.png)

- **Expected share per predictor:** $E[R^2]/k = \tfrac{1}{3k}$. Adding
  predictors does not add explained variance to the prior: a built-in
  correction for testing many predictors.
- **Shrinkage:** with several predictors at a small sample, estimates are
  pulled toward zero, because each predictor's share peaks at zero. The pull
  fades as $n$ grows.

---

## What the prior does expect

The motto "define what's possible, not what you expect" holds for the total:
$R$ is flat. It does carry expectations, and they are stated here rather than
left for a reader to find:

- flat on $R$ means $E[R^2] = \tfrac13$ and median $R^2 = 0.25$;
- the per-predictor expectation of step 4 follows from the budget.

Other rulers make other trade-offs. Flat on $R^2$ (Beta(1, 1)) shrinks less
but favours the null more for small effects; the Beta(½, b) family on $R^2$
runs between them. No member gets both properties for free: both come from
how much prior mass sits on small effects.
