# 0. Overview: what this project claims, and how far each claim is shown

**In one sentence:** with one predictor, Jeffreys' flat correlation is the
only objective choice; with several, not everything can stay flat, and the
three priors that remain are built and validated side by side, so the choice
between them follows from the results.

This page is the map. Every other document goes into one part of it.

---

## Where it started: study 00

The first result, the one that started the comparison, is a simulation
(`validation/00_initial_simulation.R`). Draw valid correlation matrices of
$(y, x_1 \ldots x_k)$, 200,000 per $k$, and see what happens to one
coefficient's partial correlation $r$ and to the total $R^2$ as $k$ grows. Three
versions, all valid correlation matrices, differing only in what is kept flat:

- **A. partial r flat:** every valid correlation matrix equally likely; each
  partial correlation is then flat on (−1, 1);
- **B. R flat:** the same matrices, with $R$ rescaled to be flat on (0, 1);
- **C. R² flat:** the same matrices, with $R^2$ rescaled to be flat on (0, 1).

| $k$ | A: P(\|partial r\| > 0.3) | A: P($R^2$ > 0.5) | B: P(\|partial r\| > 0.3) | B: P($R^2$ > 0.5) | C: P(\|partial r\| > 0.3) | C: P($R^2$ > 0.5) |
|---|---|---|---|---|---|---|
| 1 | 0.70 | 0.29 | 0.70 | 0.29 | 0.91 | 0.50 |
| 2 | 0.70 | 0.50 | 0.50 | 0.29 | 0.70 | 0.50 |
| 3 | 0.70 | 0.64 | 0.41 | 0.29 | 0.60 | 0.50 |
| 5 | 0.70 | 0.82 | 0.32 | 0.29 | 0.49 | 0.50 |
| 10 | 0.70 | 0.97 | 0.22 | 0.29 | 0.36 | 0.50 |

Keep each partial correlation flat and the total grows with $k$; keep the total
flat and each partial correlation shrinks with $k$. With one predictor A and B
are the same (Jeffreys); C is not flat on $r$ even then. The figure
(`validation/results/figures/00_partial_r_R_R2_by_k.png`, made by running
study 00) shows the three distributions per $k$; the table is also written to
`validation/results/00_first_result.md`.

Studies 01–10 test the three options this leaves.

---

## The three branches

All three are one family: **$R^2 \sim \mathrm{Beta}(a, 1)$ with a uniform
direction** (every split of the total over the predictors equally likely, in
explained variance). They differ only in the shape $a$:

```
                         k = 1: flat on r (Jeffreys)
                                  |
        ---------------------------+---------------------------
        |                          |                          |
  FLAT PARTIAL R              FLAT R                     FLAT R-SQUARED
  every correlation matrix    the total R flat           the total R² flat
  equally likely              for every k                for every k
  a = k/2                     a = 1/2                    a = 1
        |                          |                          |
  R² ~ Beta(k/2, 1)           R² ~ Beta(1/2, 1)          R² ~ Beta(1, 1)
  grows with k                the same for every k       the same for every k
  median 0.25, 0.63, 0.87     median 0.25                median 0.50
  (k = 1, 3, 10)                                         
        |                          |                          |
  per coefficient:            per coefficient:           per coefficient:
  partial r flat on (−1, 1)   its share of the total,    its share of the total,
  for every k (Jeffreys)      narrows as k grows         narrows as k grows
```

| Prior | `prior =` | $a$ | At $k = 1$ | Per coefficient | Where it is derived |
|---|---|---|---|---|---|
| flat partial r | `"flat_partial_r"` | $k/2$ | Jeffreys | flat partial correlation, any $k$ | [`relation_to_others.md`](relation_to_others.md), "How Beta(k/2, 1) was found" |
| flat R | `"flat_r"` | ½ | Jeffreys | share ~ Beta(½, (k−1)/2) of the total | [`1_from_behaviour.md`](1_from_behaviour.md), step 4 |
| flat R-squared | `"flat_r_squared"` | 1 | flat on $r^2$, not on $r$ | share ~ Beta(½, (k−1)/2) of the total | [`2_how_it_calculates.md`](2_how_it_calculates.md), section 2 |

**There is no default.** Every function asks for `prior =`. The choice is to be
made from the comparison in [`4_comparison.md`](4_comparison.md), which is
filled from the validation reruns.

Flat partial r is rstanarm `stan_lm`'s $R^2$ prior with η = 1. All three are
mixtures of Zellner g-priors, the family of the JZS default, for the shapes
where such a mixture exists ([`2_how_it_calculates.md`](2_how_it_calculates.md),
section 8); each Bayes factor is also computed by one route that works for
every $a$ (section 7). What singles each one out:

| | Unique property | Shown in |
|---|---|---|
| flat partial r | the only `stan_lm` shape (η = 1) where every coefficient's partial correlation is Jeffreys' flat correlation, whatever $k$ is | [`relation_to_others.md`](relation_to_others.md), "Why η = 1"; validation 05, C7 |
| flat R | the only g-mixture whose $R^2$ prior is the same for every $k$ *and* is Jeffreys at $k = 1$ | doc 2, section 8; validation 09, V5 |
| flat R-squared | the total explained share is flat on its own scale, for every $k$ | by construction; validation 00, F3 |

---

## The three questions the package answers

| Question | Function | How it computes | Explained in |
|---|---|---|---|
| What are the coefficients? | `fit_regression()` | Stan, sampling $R$ and a direction (not the coefficients) | doc 2, section 6 |
| Which of two models? (nested or not) | `bayes_factor_models()`, `log_bayes_factor_vs_intercept_only()` | exact: needs only $R^2$, $n$, $k$ | doc 2, sections 7–9 |
| Is there evidence for this one coefficient? | `bayes_factor_per_term()`, `log_bayes_factor_conditional()` | the conditional test: its partial correlation, with the prior the joint prior implies for it | doc 2, section 10 |

Two ways to test one coefficient, and why the package uses the second:

| | Joint test (`test = "joint"`) | Conditional test (default) |
|---|---|---|
| What it compares | the model with against without the term | the term's partial correlation against zero |
| Depends on what the other predictors explain? | yes (the conditional Lindley paradox) | no |
| Depends on how many other predictors there are? | yes | flat partial r: no; flat R and flat R-squared: yes (the share shrinks with $k$) |
| Used for | only when asked (`test = "joint"`), or to compare whole models | every term; a term with several columns (a factor with more than two levels) as a block, by its partial $R^2$ |

---

## Which document for what

| Document | Read it for |
|---|---|
| [`1_from_behaviour.md`](1_from_behaviour.md) | *Why* these priors: from study 00 and what effects can do, in steps with simulated figures |
| [`2_how_it_calculates.md`](2_how_it_calculates.md) | *How* it computes: the density, the Stan model, the Bayes factor, the conditional test, the default-prior criteria |
| [`3_consistency.md`](3_consistency.md) | The proofs and sketches behind the consistency criteria, for the whole family |
| [`relation_to_others.md`](relation_to_others.md) | How the three relate to JZS, `stan_lm`, R2-D2 and other priors |
| [`4_comparison.md`](4_comparison.md) | The three compared per criterion, from the validation results; the recommendation (flat partial r) |
| [`5_validation.md`](5_validation.md) | Every check and simulation study, with its numbers |
| [`../vignette.R`](../vignette.R) | Use it on one dataset, start to finish |

---

## One name per idea

| Name used | Means | Also appears as |
|---|---|---|
| **total** | the multiple correlation $R$ of the whole model | "shared total" |
| **shape** $a$ | the first parameter of $R^2 \sim$ Beta($a$, 1): $k/2$, ½ or 1 | `r_squared_shape` (Stan), `prior_shape()` |
| **share** | the fraction of the total one predictor gets; Beta(½, (k−1)/2) for every $a$ | $c^2$ (squared alignment, doc 2 section 7), $B$ (doc 2 section 8, doc 3) |
| **noise units**, $f$ | an effect divided by the noise SD $\sigma$; $f^2$ is Cohen's $f^2$ | |
| **joint test** | the Bayes factor of the model with against without a term | |
| **conditional test** | the test on one coefficient's partial correlation | "per-coefficient test" |
| **block route** | a model's $k$ columns tested as one block against none: the whole-model Bayes factor for any $a$ | doc 2, section 7 |
| **flat partial r** | $a = k/2$ | "uniform correlations" (earlier versions), `stan_lm` with η = 1, LKJ(1) |
| **flat R** | $a = ½$ | "flat-R" (earlier versions) |
| **flat R-squared** | $a = 1$ | "flat R²" |
| **references** | JZS default (Cauchy, scale 0.354) and Zellner–Siow (scale 1) | |

---

## Status of each claim

**proven** = a complete proof in doc 3; **checked** = a numerical check that
passes on the current code; **sketch** = the steps of a proof, with a check;
**run** = a simulation study run on the current code, with its result; details
in [`5_validation.md`](5_validation.md). All from the full run of
4 November 2026 (commit e221561): 175 of 176 checks pass; the one failure is a
finding about flat R-squared (validation 10, P4; see
[`4_comparison.md`](4_comparison.md), section 5).

| Claim | Flat partial r | Flat R | Flat R-squared | Where |
|---|---|---|---|---|
| Valid correlation matrices give the three versions | checked (F1) | checked (F2) | checked (F3) | validation 00 |
| The Stan model samples the prior correctly | checked (P5–P7, KS ≤ 0.012) | checked (P1–P4, KS ≤ 0.009) | checked (P8–P10, KS ≤ 0.009) | validation 01, P1–P10 |
| The Stan posterior is calibrated ($k = 8$) | checked (SBC p ≥ 0.43) | checked (SBC p ≥ 0.61) | checked (SBC p ≥ 0.45) | validation 01, S |
| The whole-model Bayes factor is exact | checked (A1–A4) | checked (A1–A4) | checked (A1–A4) | validation 02 |
| Units and recoding don't change the answer | checked (A2) | checked (A2) | checked (A2) | validation 02 |
| Predictive matching (BF = 1 at $n = k + 1$) | checked (B1) | checked (B1) | checked (B1) | validation 02 |
| Information consistency, for $n \ge k + 3$ | proven, checked | proven, checked | proven, checked | doc 3; validation 02, B2b |
| Intrinsic consistency | proven | proven | proven | doc 3 |
| Model-selection consistency, fixed $k$ | sketch, checked: $\tfrac q2 \log n$ in every case | sketch, checked: ½ log n against the intercept-only model, $\tfrac q2 \log n$ otherwise | sketch: $\log n$ against the intercept-only model, $\tfrac q2 \log n$ otherwise (slopes 0.99, 1.00, 1.01 and 0.53, 1.64, 2.88) | doc 3; validation 02, B3, B3b, B4 |
| Model-selection consistency, $k$ growing with $n$ | checked (G1, G2): evidence grows fast (log BF₀₁ 408 at $n = 2000$, $k = n/4$) | checked, but slowly under the null (log BF₀₁ 1.9) | checked, but slowly under the null (log BF₀₁ 3.7) | validation 08 |
| What singles the prior out | checked (C7) | checked (V5) | – | validation 05, 09 |
| The conditional test is computed correctly | checked (C1, C3–C5) | checked (C1, C2, C5) | checked (C1, C2, C5) | validation 05, C1–C5 |
| The block conditional test is computed correctly | checked (C9–C11) | checked (C9–C11) | checked (C9–C11) | validation 05, C9–C11 |
| Separating nulls from small effects, against JZS | run: equal or better (AUC +0.001 to +0.045) | run: the same as flat partial r (conditional) | run: the same as flat partial r (conditional) | validation 04, 06 |
| Confirming a true null, one coefficient | run: 86–88%, for any $k$ | run: 55% at $k = 2$, 6% at $k = 20$ | run: 86% at $k = 2$, 26% at $k = 20$ | validation 06 |
| Calling a small effect (0.1) null | run: 52–60% | run: 19% at $k = 2$, 0.5% at $k = 20$ | run: 52% at $k = 2$, 4% at $k = 20$ | validation 06 |
| Calibrated on its own prior | checked (P1) | checked (P1) | checked (P1) | validation 10 |
| Variable selection | run: fewest false inclusions of the three, misses more | run: most inclusions, fewest misses | run: close to flat R | validation 09 |
| Estimation: prediction error, interval coverage | run: better than OLS everywhere, coverage ≥ 0.916 | run: lowest error at small $n$, coverage down to 0.76 | run: in between, coverage down to 0.82 | validation 07 |

---

## What to watch out for when using it

- **Choose the prior.** There is no default; `prior =` must be given. What each
  one assumes is in the branches above and in [`4_comparison.md`](4_comparison.md).
- **Shrinkage at small $n$.** With several predictors at small $n$ the prior
  pulls estimates toward zero. The vignette prints the pull next to the
  ordinary estimates.
- **Flat R and flat R-squared: the conditional test is cautious with many
  predictors.** The prior expects each coefficient to be a small share of the
  total, so "no effect" and "a tiny effect" predict nearly the same data. Under
  flat partial r each coefficient keeps a flat partial correlation, whatever
  $k$ is.
- **Leaving the null model.** Against the intercept-only model, evidence for
  the null grows like $a \log n$: ½ log n for flat R and log n for flat
  R-squared, whatever the number of predictors; $\tfrac k2 \log n$ for flat
  partial r. How strict a variable selection is then rests partly on the prior
  over models: choose it deliberately.
- **Flat partial r: the claimed $R^2$ grows with $k$.** Median 0.63 with 3
  predictors, 0.87 with 10. A large model whose fit is modest is then penalised.
- **Flat R-squared with one predictor is not Jeffreys.** Its prior on $r$ has
  density $|r|$, which gives little weight to correlations near zero.
- **Information consistency needs $n \ge k + 3$**, for all three.
- **Main effects with interactions.** `bayes_factor_per_term()` removes one
  term at a time, so testing a main effect while its interaction stays in the
  model follows the usual marginality caveats.
- **Stan messages.** At the start of warmup Stan may print a few "Informational
  Message ... proposal is about to be rejected" lines (the sampler tried $R$ at
  its boundary). A few are harmless; check the divergent transitions that
  `print(fit)` reports.
- **Why not brms.** brms samples the coefficients directly, and for $k \ge 2$ the
  prior has a sharp peak at all-zero coefficients that such a sampler misses,
  without warnings (doc 2, section 6).

---

## Open

- **Model-selection consistency, in full:** the dominated-convergence steps of
  the sketch, non-nested models, and a proof for $k$ growing with $n$ (doc 3).
- **The empirical question between the branches:** how explained variance
  grows with the number of predictors in published regressions. Flat partial r
  predicts a median $R^2$ of 0.76 with 5 predictors; flat R predicts 0.25 and
  flat R-squared 0.5 at any size.
- **Benchmarks in variable selection** against the robust, hyper-g/n and
  intrinsic priors (BAS package), and a real-data reanalysis.
