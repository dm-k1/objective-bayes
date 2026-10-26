# 4. Why uniform correlations is the default

**In one sentence:** the weights of the prior follow the space of what is
possible, every correlation structure equally likely, instead of a shape and
a scale chosen by convention; and in the validation studies that choice is the
only one whose test of a coefficient does not change with the number of
covariates.

This document gives the stance behind the default, how it compares with the
standard defaults (JZS, Zellner–Siow, BIC) and with flat-R, and the one cost
that follows from the stance. The numbers come from the validation reruns of
26 October 2026 ([`../validation/README.md`](../validation/README.md)).

---

## 1. The stance: a prior as a shared standard

Probability can be read as a degree of belief, private to each researcher
(de Finetti). But when researchers have to agree on what the data show, a
private prior is not enough. They agree on a **convention**, and a convention
can be defended only if it comes from a principle anyone can check, not from
someone's taste. This is the objective-Bayes position (Jeffreys; Bayarri,
Berger, Forte & García-Donato, 2012): a default prior is a public standard.

The principle here: **the weights follow the space of what is possible**
(the motto "define what's possible, not what you expect").

| Step | What it says |
|---|---|
| What is possible | every valid correlation matrix of the outcome and the predictors; nothing else (a set of correlations that would explain more than 100% of the variance is not a valid matrix) |
| How the weights are spread | evenly over that space: every possible correlation structure equally likely (LKJ with η = 1) |
| What follows for one coefficient | its partial correlation is flat on (−1, 1), Jeffreys' prior, whatever the number of other predictors ([`relation_to_others.md`](relation_to_others.md), "Why η = 1") |

**What the other defaults do differently.** JZS and Zellner–Siow do not leave
out anything that can happen: they give some probability to every possible
effect, as this prior does. The difference is **where the weights come from**:

| | Where the prior lives | Where its weights come from |
|---|---|---|
| uniform correlations | the correlations: a bounded space | an even spread over everything possible; no number to choose |
| JZS default (BayesFactor, JASP) | the coefficients: an unbounded space | a Cauchy shape with scale 0.354, chosen by convention |
| Zellner–Siow | the coefficients | a Cauchy with scale 1 |

So the objection to JZS is not that it ignores possibilities, but that its
weights rest on a chosen number with no principle behind it.

**Still a choice, stated openly.** "Evenly" is evenly *on the correlation
scale*, as Jeffreys' prior is flat on $r$ rather than on $r^2$
([`1_from_behaviour.md`](1_from_behaviour.md), step 1: "which ruler is flat").
Correlations are the bounded, comparable scale researchers reason in, which
makes this the natural ruler; it is a judgement, not a theorem.

---

## 2. The comparison with the standard defaults

### Testing one coefficient (study 06, averaged over designs; k = 2 → 20)

| | True null found (BF₀₁ > 3) | Small effect (0.1) called null | Effect 0.2 found |
|---|---|---|---|
| **uniform correlations** (conditional test, the default) | **0.86 → 0.88**, stable in k | 0.52 → 0.60 | 0.51 → 0.44 |
| Zellner–Siow (scale 1, joint test) | 0.87 → 0.52, falls with k | 0.53 → 0.22 | 0.51 → 0.52 |
| JZS default (scale 0.354, joint test) | 0.54 → 0.30, falls with k | 0.21 → 0.11 | 0.56 → 0.52 |

False effects for a true null: at most 3.3% for every method.

Separating nulls from small effects (AUC, null against 0.1, minus JZS):

| | uncorrelated predictors | correlation 0.7 |
|---|---|---|
| uniform correlations (conditional) | +0.001 to +0.007 | **+0.031 to +0.045** |
| Zellner–Siow | −0.000 to +0.003 | +0.012 to +0.027 |

### Whole models and estimation (studies 07, 08, 09)

| | Variable selection (09) | Null evidence, k = n/4 at n = 2000 (08, median log BF₀₁) | Prediction error vs OLS, k = 8, n = 30 (07) | Lowest coverage of 95% intervals, real effects (07) |
|---|---|---|---|---|
| uniform correlations | strict end; ranks at least as well as JZS in every cell | 408 | 0.74 | **0.92** |
| Zellner–Siow | strict end, next to uniform correlations | 411 | 0.77 | 0.91 |
| JZS default | in between | 149 | **0.56** | 0.87 |
| BIC | strictest at the 0.5 cutoff, yet ranks best | — | — | — |

### What the comparison says

- **Against Zellner–Siow:** for whole models and estimation nearly the same
  prior (η = 1 against η = ½ in the same family), slightly better in ranking and
  interval coverage: *Zellner–Siow with a reason*, the same behaviour derived
  from "every correlation structure equally possible" instead of a Cauchy chosen
  by convention. For one coefficient it differs: Zellner–Siow's test drifts with
  the number of covariates (null found 0.87 → 0.52), this one does not.
- **Against JZS default:** JZS is more sensitive to very small effects and
  shrinks harder, which gives better predictions at small n. It pays with
  finding nulls far less often (0.30 at k = 20), intervals that undercover more
  (down to 0.87), and a scale of 0.354 that is a convention.
- **Against flat-R** (the other extension of Jeffreys' flat correlation, built
  here too): flat-R's R² prior stays the same for every k, and that costs the
  evidence for a null once covariates are in the model (null found 0.55 at k = 2,
  0.06 at k = 20) and interval coverage (down to 0.76). Details:
  [`0_overview.md`](0_overview.md) and [`relation_to_others.md`](relation_to_others.md).

---

## 3. The one cost, and why it follows from the stance

The default calls an effect of 0.1 (a partial correlation of about 0.1, below
Cohen's "small") null in about half the datasets at moderate n (study 06;
34% at n = 500 in study 03, and falling as n grows).

This follows from the principle, not from a flaw in it. A flat partial
correlation spreads its probability evenly over (−1, 1), so only 10% of it lies
below |ρ| = 0.1. Weighting every possibility equally therefore treats an effect
of 0.1 as close to "no effect": the data have to be strong before a correlation
that small is told apart from zero.

So the rule for users:

> **A researcher who expects effects this small has to say so in the prior; the
> default will not assume it for them.**

That is the same rule the stance applies everywhere: what the prior assumes is
stated, never hidden in a chosen scale.

---

## 4. In one paragraph, for a paper

> We use a prior that spreads probability evenly over every correlation
> structure that is possible for the outcome and the predictors, and gives none
> to what is impossible. It needs no scale, and each coefficient's partial
> correlation receives Jeffreys' flat prior whatever the number of covariates.
> In simulations it is the only default we compared (JZS, Zellner–Siow, flat-R)
> whose test of a single coefficient behaves the same with 2 or 20 covariates;
> it confirms true nulls more often than JZS and keeps interval coverage near
> nominal, at the cost of treating partial correlations of about 0.1 as close to
> zero at moderate sample sizes. Its Bayes factors are calibrated under the
> prior, and the per-coefficient test, the fixed-predictor form of the exact
> flat partial-correlation test, agrees with the joint prior's Savage–Dickey
> ratio to within a few percent.

---

## 5. Checks of the justification

Four checks on whether the case above holds (26 October 2026; all pass).

| Question | Result | Where |
|---|---|---|
| **Is the per-coefficient test new?** | Not as a test. It is the **fixed-predictor version** of the exact partial-correlation Bayes factor of Kucharský et al. (2026) with a flat prior. Giving ours their one extra assumption (random predictors, so the tested predictor's residual SD is unknown) makes the two agree to 5e-7. Without it, ours is slightly lower: at most −0.26 on the log scale (n = 20, k = 8, r = 0), shrinking like 1/n | validation 05, C8; [`relation_to_others.md`](relation_to_others.md) |
| **Do its Bayes factors mean what they say?** | Yes, under its own prior: datasets drawn from the joint uniform-correlations prior, half of them null, give the share of nulls per bin of BF predicted by 1/(1 + BF₁₀) (χ² p = 0.40–0.79, calibration slope 0.98–1.03). Data drawn from the JZS prior instead are clearly miscalibrated (χ² 525–593): calibration is always relative to a prior | validation 10, P1–P2 |
| **Is "one prior for everything" true in practice?** | Closely, for the conditional test: it is within 0.092 on the log scale of the Savage–Dickey ratio under the joint prior (mean absolute difference 0.028; 16 of 20 datasets where the ratio can be estimated, the other 4 within their Monte Carlo error). The joint test is further off (mean 0.23, up to 0.42) | validation 10, P3–P5; [`2_how_it_calculates.md`](2_how_it_calculates.md), section 10 |
| **Does it cover factors?** | Yes: a block version for terms with q columns, with the block's partial R² ~ Beta(q/2, 1) under uniform correlations for every k; `bayes_factor_per_term()` uses it for multi-column terms | validation 05, C9–C11; doc 2, section 10 |

**What this means for section 3's cost.** Among datasets with BF₀₁ > 3 under
this prior, 83–89% are true nulls, and the real effects among them are tiny:
97–98% have a partial correlation below 0.2 at n = 200 (78–80% at n = 50). So
calling an effect of 0.1 null is the correct posterior verdict *under this
prior*, not a failure of the test. This confirms the computation and shows the cost is the
prior's, not the test's; whether the prior fits a field is the question of
section 1 (every correctly computed Bayes factor is calibrated under its own
prior, so calibration alone cannot justify one).

**What this means for the claim of novelty.** The test itself exists. What is
new is that it follows from one joint prior that also gives the whole-model
Bayes factors and the estimates, that this joint prior is singled out by η = 1
(the only shape with a flat partial correlation for every k), and the fixed-X
form that conditions on the observed predictors as regression does.

**Still open:** a simulation study of the block test against JZS and the joint
test (started 26 October 2026), and the empirical question of how R² grows with
the number of predictors in published regressions ([`0_overview.md`](0_overview.md)).
