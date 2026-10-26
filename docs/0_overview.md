# 0. Overview: what this project claims, and how far each claim is shown

**In one sentence:** with one predictor, a flat prior on the correlation; with
$k$ predictors, a flat prior on the total ($R$, the multiple correlation), and
every way of dividing that total over the predictors equally plausible.

This page is the map. Every other document goes into one part of it.

---

## The two paths from one predictor to many

Start: one predictor, a flat prior on its correlation $r$ (Jeffreys). With $k$
predictors you must choose *what* stays flat. There are two natural choices:

```
                     k = 1: flat on r
                            |
          ------------------+-------------------
          |                                    |
  A. FLAT-R (this package)          B. UNIFORM CORRELATIONS
  flat on the total R               every correlation matrix of
                                    (y, x1 ... xk) equally likely
          |                                    |
  total:  R² ~ Beta(½, 1)           total:  R² ~ Beta(k/2, 1)
          the same for every k              grows with k
          median 0.25                       median 0.25, 0.63, 0.87 (k = 1, 3, 10)
          |
  per predictor:
  share ~ Beta(½, (k−1)/2)
  on average 1/k of the total
```

| Beta | Path | What it is about | Where it is derived |
|---|---|---|---|
| Beta(½, 1) | A | $R^2$ of the whole model | [`2_how_it_calculates.md`](2_how_it_calculates.md), section 2 |
| Beta(½, (k−1)/2) | A | one predictor's **share** of that $R^2$ | [`1_from_behaviour.md`](1_from_behaviour.md), step 4 (simulate → derive → check) |
| Beta(k/2, 1) | B | $R^2$ of the whole model | [`relation_to_others.md`](relation_to_others.md), "How Beta(k/2, 1) was found" |

Path B is rstanarm `stan_lm`'s $R^2$ prior with η = 1. Both paths are mixtures of
Zellner g-priors, the family of the JZS default
([`2_how_it_calculates.md`](2_how_it_calculates.md), section 8). Flat-R is the
one member whose $R^2$ prior is the same for every $k$ *and* is Jeffreys' flat
correlation at $k = 1$.

---

## The three questions the package answers

| Question | Function | How it computes | Explained in |
|---|---|---|---|
| What are the coefficients? | `fit_flat_r_model()` | Stan, sampling $R$ and a direction (not the coefficients) | doc 2, section 6 |
| Which of two models? (nested or not) | `bayes_factor_flat_r()` | closed form: needs only $R^2$, $n$, $k$ | doc 2, sections 7–9 |
| Is there evidence for this one coefficient? | `bayes_factor_flat_r_per_term()` | the conditional test: its partial correlation, with its share of the flat total as prior | doc 2, section 10 |

Two ways to test one coefficient, and why the package uses the second:

| | Joint test (`test = "joint"`) | Conditional test (default) |
|---|---|---|
| What it compares | the model with against without the term | the term's partial correlation against zero |
| Depends on what the other predictors explain? | yes (the conditional Lindley paradox) | no |
| Separates nulls from small effects, against JZS | worse | equal, better with correlated predictors |
| Used for | factors with more than two levels (a block version is not built yet) | every one-column term |

---

## Which document for what

| Document | Read it for |
|---|---|
| [`1_from_behaviour.md`](1_from_behaviour.md) | *Why* this prior: from what effects can do to the one sentence, in four steps with simulated figures |
| [`2_how_it_calculates.md`](2_how_it_calculates.md) | *How* it computes: the density, the Stan model, the Bayes factor, the conditional test, the default-prior criteria |
| [`3_consistency.md`](3_consistency.md) | The proofs and sketches behind the consistency criteria |
| [`relation_to_others.md`](relation_to_others.md) | How it relates to JZS, `stan_lm` (path B), R2-D2 and other priors |
| [`../validation/README.md`](../validation/README.md) | Every check and simulation study, with its numbers |
| [`../vignette.R`](../vignette.R) | Use it on one dataset, start to finish |

---

## One name per idea

| Name used | Means | Also appears as |
|---|---|---|
| **total** | the multiple correlation $R$ of the whole model; flat from 0 to 1 | "shared total" |
| **share** | the fraction of the total one predictor gets; Beta(½, (k−1)/2) | $c^2$ (squared alignment, doc 2 section 7), $B$ (doc 2 section 8, doc 3) |
| **noise units**, $f$ | an effect divided by the noise SD $\sigma$; $f^2$ is Cohen's $f^2$ | |
| **joint test** | the Bayes factor of the model with against without a term | |
| **conditional test** | the test on one coefficient's partial correlation | "per-coefficient test" |
| **uniform correlations** | path B | `stan_lm` with η = 1, LKJ(1), $a = k/2$ |

---

## Status of each claim

**proven** = a complete proof in doc 3; **checked** = a numerical check that
passes on the current code; **stale** = results from a run before the current
method set, to be rerun; **not run** = written, not yet run in full.

| Claim | Status | Where |
|---|---|---|
| The Stan model samples the prior and the posterior correctly | checked | validation 01 (P1–P4, S) |
| The closed-form Bayes factor is exact (against `integrate()`, brute force, the g-mixture) | checked | validation 02 (A1–A4) |
| Units and recoding don't change the answer | checked | validation 02 (A2) |
| Predictive matching (BF = 1 at $n = k + 1$) | checked | validation 02 (B1) |
| Information consistency, for $n \ge k + 3$ | proven, checked | doc 3; validation 02 (B2, B2b) |
| Intrinsic consistency | proven | doc 3 |
| Model-selection consistency, fixed $k$: evidence for the null grows like ½ log n against the intercept-only model, like q/2 log n between models that both contain an effect | sketch, checked | doc 3; validation 02 (B3, B3b, B4) |
| Model-selection consistency, $k$ growing with $n$ | development run only | validation 08 |
| Constant in $k$ alone does not single out flat-R; with Jeffreys at $k = 1$ it does | checked | doc 2 section 8; validation 09 (V5) |
| The conditional test is computed correctly (equals the joint test at $k = 1$; implied prior; accuracy up to $n = 10^6$) | checked | validation 05 |
| The conditional test separates nulls from small effects at least as well as JZS | stale | validation 05, 06 |
| Variable selection: flat-R is the most inclusive; the trade-off curves nearly coincide | stale | validation 09 |
| Estimation quality (bias, error, coverage) against `lm()`, JZS, `stan_lm` | not run | validation 07 |

---

## What to watch out for when using it

- **Shrinkage at small $n$.** With several predictors at small $n$ the prior pulls
  estimates toward zero (each share peaks at zero). The vignette prints the pull
  next to the ordinary estimates.
- **The conditional test is cautious with many predictors.** It rarely calls a
  small effect null, but it also rarely confirms a null: the prior expects each
  coefficient to be small (on average $1/(3k)$ of the variance), so "no effect"
  and "a tiny effect" predict nearly the same data.
- **Leaving the null model is cheap.** Against the intercept-only model, evidence
  for the null grows only like ½ log n, whatever the number of predictors. So in
  variable selection flat-R is the most inclusive of the priors compared. How
  strict a selection is then rests on the prior over models: choose it deliberately.
- **Information consistency needs $n \ge k + 3$.** With fewer observations the
  evidence for a perfect fit stays bounded.
- **Main effects with interactions.** `bayes_factor_flat_r_per_term()` removes
  one term at a time, so testing a main effect while its interaction stays in
  the model follows the usual marginality caveats.
- **Stan messages.** At the start of warmup Stan may print a few "Informational
  Message ... proposal is about to be rejected" lines (the sampler tried $R$ at
  its boundary). A few are harmless; check the divergent transitions that
  `print(fit)` reports.
- **Why not brms.** brms samples the coefficients directly, and for $k \ge 2$ the
  prior has a sharp peak at all-zero coefficients that such a sampler misses,
  without warnings (doc 2, section 6).

---

## Open

- **Rerun** studies 03–06 and 09 with the common method set; run 07 (estimation)
  and 08 (growing $k$) in full; record the results in `validation/README.md` and
  update the status table above.
- **A block version of the conditional test** for multi-column terms (factors),
  with share ~ Beta(q/2, (k − q)/2) for $q$ columns.
- **Model-selection consistency, in full:** the dominated-convergence steps of the
  sketch, non-nested models, and a proof for $k$ growing with $n$ (doc 3).
- **The empirical question between paths A and B:** how explained variance grows
  with the number of predictors in published regressions. Path B predicts a
  median $R^2$ of 0.76 with 5 predictors; flat-R predicts 0.25 at any size.
- **Benchmarks in variable selection** against the robust, hyper-g/n and intrinsic
  priors (BAS package), and a real-data reanalysis.
