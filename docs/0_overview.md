# 0. Overview: what this project claims, and how far each claim is shown

**In one sentence:** two objective priors for linear regression, both
Jeffreys' flat correlation with one predictor, built and validated side by
side: every correlation matrix equally likely (uniform correlations, the
default of the package `flatcor`), or flat on the total $R$ (flat-R).

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
  A. FLAT-R                         B. UNIFORM CORRELATIONS (default)
  flat on the total R               every correlation matrix of
                                    (y, x1 ... xk) equally likely
          |                                    |
  total:  R² ~ Beta(½, 1)           total:  R² ~ Beta(k/2, 1)
          the same for every k              grows with k
          median 0.25                       median 0.25, 0.63, 0.87 (k = 1, 3, 10)
          |                                    |
  per predictor:                    per predictor:
  share ~ Beta(½, (k−1)/2)          partial correlation flat on (−1, 1)
  on average 1/k of the total       for every k (Jeffreys, per coefficient)
```

| Beta | Path | What it is about | Where it is derived |
|---|---|---|---|
| Beta(½, 1) | A | $R^2$ of the whole model | [`2_how_it_calculates.md`](2_how_it_calculates.md), section 2 |
| Beta(½, (k−1)/2) | A | one predictor's **share** of that $R^2$ | [`1_from_behaviour.md`](1_from_behaviour.md), step 4 (simulate → derive → check) |
| Beta(k/2, 1) | B | $R^2$ of the whole model | [`relation_to_others.md`](relation_to_others.md), "How Beta(k/2, 1) was found" |

Path B is rstanarm `stan_lm`'s $R^2$ prior with η = 1. Both paths are mixtures
of Zellner g-priors, the family of the JZS default
([`2_how_it_calculates.md`](2_how_it_calculates.md), section 8), and each is
singled out by one property:

| | Unique property | Shown in |
|---|---|---|
| flat-R | the only g-mixture whose $R^2$ prior is the same for every $k$ *and* is Jeffreys at $k = 1$ | doc 2, section 8; validation 09, V5 |
| uniform correlations | η = 1 is the only `stan_lm` shape where every coefficient's partial correlation is Jeffreys' flat correlation, whatever $k$ is | [`relation_to_others.md`](relation_to_others.md), "Why η = 1"; validation 05, C7 |

Every function takes `prior = "uniform_correlations"` (the default, chosen on
the validation reruns of 26 October 2026) or `prior = "flat_r"`.

---

## The three questions the package answers

| Question | Function | How it computes | Explained in |
|---|---|---|---|
| What are the coefficients? | `fit_regression()` | Stan, sampling $R$ and a direction (not the coefficients) | doc 2, section 6 |
| Which of two models? (nested or not) | `bayes_factor_models()` | closed form: needs only $R^2$, $n$, $k$ | doc 2, sections 7–9 |
| Is there evidence for this one coefficient? | `bayes_factor_per_term()` | the conditional test: its partial correlation, with the prior the joint prior implies for it (flat-R: its share of the flat total; uniform correlations: a flat partial correlation) | doc 2, section 10 |

Two ways to test one coefficient, and why the package uses the second:

| | Joint test (`test = "joint"`) | Conditional test (default) |
|---|---|---|
| What it compares | the model with against without the term | the term's partial correlation against zero |
| Depends on what the other predictors explain? | yes (the conditional Lindley paradox) | no |
| Depends on how many other predictors there are? | yes | flat-R: yes (the share shrinks with $k$); uniform correlations: no |
| Separates nulls from small effects, against JZS | worse | equal, better with correlated predictors |
| Used for | only when asked (`test = "joint"`), or to compare whole models | every term; a term with several columns (a factor with more than two levels) as a block, by its partial $R^2$ |

---

## Which document for what

| Document | Read it for |
|---|---|
| [`1_from_behaviour.md`](1_from_behaviour.md) | *Why* this prior: from what effects can do to the one sentence, in four steps with simulated figures |
| [`2_how_it_calculates.md`](2_how_it_calculates.md) | *How* it computes: the density, the Stan model, the Bayes factor, the conditional test, the default-prior criteria |
| [`3_consistency.md`](3_consistency.md) | The proofs and sketches behind the consistency criteria |
| [`relation_to_others.md`](relation_to_others.md) | How it relates to JZS, `stan_lm` (path B), R2-D2 and other priors |
| [`4_why_this_default.md`](4_why_this_default.md) | *Why* uniform correlations is the default: the stance (weights follow what is possible), the comparison with JZS, Zellner–Siow, BIC and flat-R, and its one cost |
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
passes on the current code; **sketch** = the steps of a proof, with a check;
**run** = a simulation study run on the current code (26 October 2026), with
its result; details in [`../validation/README.md`](../validation/README.md).

| Claim | Flat-R | Uniform correlations | Where |
|---|---|---|---|
| The Stan model samples the prior correctly | checked (P1–P4) | checked (P5–P7) | validation 01 |
| The Stan posterior is calibrated ($k = 8$) | checked | checked | validation 01, S |
| The closed-form Bayes factor is exact | checked (A1–A4) | checked (CA1, CA3) | validation 02 |
| Units and recoding don't change the answer | checked (A2) | checked (CA2) | validation 02 |
| Predictive matching (BF = 1 at $n = k + 1$) | checked (B1) | checked (CB1) | validation 02 |
| Information consistency, for $n \ge k + 3$ | proven, checked | proven, checked | doc 3; validation 02 (B2b, CB2) |
| Intrinsic consistency | proven | proven | doc 3 |
| Model-selection consistency, fixed $k$ | sketch, checked: ½ log n against the intercept-only model, q/2 log n otherwise | sketch, checked: q/2 log n in every case | doc 3; validation 02 (B3, B3b, B4; CB3, CB4) |
| Model-selection consistency, $k$ growing with $n$ | run: consistent, but slow (log BF₀₁ 1.9 at n = 2000, k = n/4) | run: consistent (408) | validation 08 |
| What singles the prior out | checked (V5) | checked (C7) | validation 09, 05 |
| The conditional test is computed correctly | checked (C1, C2, C5) | checked (C4, C6) | validation 05 |
| The block conditional test (terms with $q$ columns) is computed correctly | checked (C9–C11) | checked (C9–C11) | validation 05 |
| The conditional test separates nulls from small effects at least as well as JZS | run: yes (the two tests rank identically) | run: yes | validation 04, 06 |
| Confirming a true null, one coefficient | run: falls with k (0.55 at k = 2, 0.06 at k = 20) | run: about 0.87 at every k | validation 06 |
| Calling a small effect (0.1) null | run: rare (0.19 to 0.005) | run: about half (0.52–0.60) | validation 06 |
| The conditional test is calibrated on its own prior (share of nulls = 1/(1 + BF₁₀)) | – | checked (P1); miscalibrated on JZS data, as expected (P2) | validation 10 |
| The conditional test is close to one prior for everything (Savage–Dickey under the joint model) | – | checked: within 0.09 on the log scale (mean 0.03) where estimable; the joint test is further away (up to 0.42) | validation 10 |
| Variable selection: the trade-off curves nearly coincide | run: yes; the most inclusive | run: yes; next to Zellner–Siow | validation 09 |
| Estimation: prediction error, interval coverage | run: largest gains at small n; coverage down to 0.76 | run: moderate gains; coverage 0.92 or higher | validation 07 |

---

## What to watch out for when using it

- **Shrinkage at small $n$.** With several predictors at small $n$ the prior pulls
  estimates toward zero (each share peaks at zero). The vignette prints the pull
  next to the ordinary estimates.
- **Flat-R: the conditional test is cautious with many predictors.** It rarely
  calls a small effect null, but it also rarely confirms a null: the prior
  expects each coefficient to be small (on average $1/(3k)$ of the variance), so
  "no effect" and "a tiny effect" predict nearly the same data. Under uniform
  correlations each coefficient keeps a flat partial correlation, whatever $k$ is.
- **Flat-R: leaving the null model is cheap.** Against the intercept-only model,
  evidence for the null grows only like ½ log n, whatever the number of
  predictors. So in variable selection flat-R is the most inclusive of the priors
  compared. How strict a selection is then rests on the prior over models: choose
  it deliberately. Uniform correlations has the usual k/2 log n.
- **Uniform correlations: the claimed $R^2$ grows with $k$.** Median 0.63 with 3
  predictors, 0.87 with 10. A large model whose fit is modest is then penalised.
- **Information consistency needs $n \ge k + 3$.** With fewer observations the
  evidence for a perfect fit stays bounded.
- **Main effects with interactions.** `bayes_factor_per_term()` removes
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

- **The block conditional test in simulation studies.** It is built and checked
  for correctness (validation 05, C9–C11; doc 2, section 10), but not yet
  compared with JZS and the joint test on factors in a study like 04 or 06.
- **Model-selection consistency, in full:** the dominated-convergence steps of the
  sketch, non-nested models, and a proof for $k$ growing with $n$ (doc 3).
- **Decided: uniform correlations is the default** (26 October 2026), on the
  reruns (`validation/README.md`, "In short"): for testing one coefficient it is
  stable in $k$, confirms nulls and keeps its interval coverage, at the cost of
  calling effects of about 0.1 null at moderate $n$. Flat-R stays available as
  `prior = "flat_r"`. Still worth checking: the empirical question below.
- **The empirical question between paths A and B:** how explained variance grows
  with the number of predictors in published regressions. Path B predicts a
  median $R^2$ of 0.76 with 5 predictors; flat-R predicts 0.25 at any size.
- **Benchmarks in variable selection** against the robust, hyper-g/n and intrinsic
  priors (BAS package), and a real-data reanalysis.
