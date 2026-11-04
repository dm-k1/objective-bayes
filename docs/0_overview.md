# 0. Overview: what this project claims, and how far each claim is shown

With one predictor, Jeffreys' flat prior on the correlation is the natural
objective choice. With several predictors, not every correlation-based
quantity can stay flat, which leaves three priors. This project builds and
validates all three side by side and lets the results decide between them;
they recommend flat partial r ([`3_comparison.md`](3_comparison.md),
section 8).

This page is the map of the project, and every other document covers one part
of it. The finding itself is stated on one page in [`findings.md`](findings.md).

---

## Where it started: study 00

The comparison started from a simulation (`validation/00_initial_simulation.R`).
It draws valid correlation matrices of $(y, x_1 \ldots x_k)$, 200,000 per $k$,
and follows one coefficient's partial correlation $r$ and the total $R^2$ as
$k$ grows. The three versions start from the same matrices and differ only in
what is kept flat:

- **A. partial r flat:** every valid correlation matrix is equally likely, so
  each partial correlation is flat on (−1, 1);
- **B. R flat:** the same matrices, with $R$ rescaled to be flat on (0, 1);
- **C. R² flat:** the same matrices, with $R^2$ rescaled to be flat on (0, 1).

| $k$ | A: P(\|partial r\| > .30) | A: P($R^2$ > .50) | B: P(\|partial r\| > .30) | B: P($R^2$ > .50) | C: P(\|partial r\| > .30) | C: P($R^2$ > .50) |
|---|---|---|---|---|---|---|
| 1 | .70 | .29 | .70 | .29 | .91 | .50 |
| 2 | .70 | .50 | .50 | .29 | .70 | .50 |
| 3 | .70 | .65 | .41 | .29 | .60 | .50 |
| 5 | .70 | .82 | .32 | .29 | .50 | .50 |
| 10 | .70 | .97 | .22 | .29 | .36 | .50 |

Keeping each partial correlation flat makes the total grow with $k$; keeping
the total flat makes each partial correlation shrink with $k$. With one
predictor, A and B coincide (Jeffreys' prior), while C is not flat on $r$
even then. The figure `validation/results/figures/00_partial_r_R_R2_by_k.png`
(made by running study 00) shows the three distributions for each $k$, and
the table is also written to `validation/results/00_first_result.md`.

Studies 01–11 test the three options this leaves.

---

## The three branches

All three priors belong to one family: **$R^2 \sim \mathrm{Beta}(a, 1)$ with a
uniform direction**, where the **direction** (= how the total is split over
the predictors, in explained variance) is equally likely to be any split.
They differ only in the shape $a$:

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

The default prior is **flat partial r**; a function that uses it because no
prior was given says so once per session, and `prior =` selects the others.
The comparison in [`3_comparison.md`](3_comparison.md), based on the validation
run of 4 November 2026, recommends flat partial r (section 8) and lists its
costs; flat R and flat R-squared remain available.

Flat partial r is the $R^2$ prior of rstanarm's `stan_lm` with η = 1. All
three are mixtures of Zellner g-priors, the family of the JZS default, for the
shapes where such a mixture exists ([`2_how_it_calculates.md`](2_how_it_calculates.md),
section 8). Each Bayes factor is also computed by one route that works for
every $a$ (section 7). Each prior has one property that singles it out:

| | Unique property | Shown in |
|---|---|---|
| flat partial r | the only `stan_lm` shape (η = 1) where every coefficient's partial correlation is Jeffreys' flat correlation, whatever $k$ is | [`relation_to_others.md`](relation_to_others.md), "Why η = 1"; validation 05, C7 |
| flat R | the only g-mixture whose $R^2$ prior is the same for every $k$ *and* is Jeffreys' prior at $k = 1$ | doc 2, section 8; validation 09, V5 |
| flat R-squared | the total explained share is flat on its own scale, for every $k$ | by construction; validation 00, F3 |

---

## The three questions the package answers

| Question | Function | How it computes | Explained in |
|---|---|---|---|
| What are the coefficients? | `fit_regression()` | Stan, sampling $R$ and a direction (not the coefficients) | doc 2, section 6 |
| Which of two models? (nested or not) | `bayes_factor_models()`, `log_bayes_factor_vs_intercept_only()` | exact: needs only $R^2$, $n$, $k$ | doc 2, sections 7–9 |
| Is there evidence for this one coefficient? | `bayes_factor_per_term()`, `log_bayes_factor_conditional()` | the conditional test: its partial correlation, with the prior the joint prior implies for it | doc 2, section 10 |

There are two ways to test one coefficient, and the package uses the second
by default:

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
| [`findings.md`](findings.md) | What the project shows: the claim and the evidence for it, on one page |
| [`1_from_behaviour.md`](1_from_behaviour.md) | *Why* these priors: from study 00 and what effects can do, in steps with simulated figures |
| [`2_how_it_calculates.md`](2_how_it_calculates.md) | *How* it computes: the density, the Stan model, the Bayes factor, the conditional test, the default-prior criteria |
| [`relation_to_others.md`](relation_to_others.md) | How the three relate to JZS, `stan_lm`, R2-D2 and other priors; what is new here and what is not |
| [`3_comparison.md`](3_comparison.md) | The three compared per criterion, from the validation results; the recommendation (flat partial r) |
| [`4_validation.md`](4_validation.md) | Every check and simulation study, with its numbers |
| [`example.R`](example.R) | The package used on one dataset, start to finish |

---

## One name per idea

| Name used | Means | Also appears as |
|---|---|---|
| **total** | the explained variance $R^2$ of the whole model (flat R keeps its square root $R$ flat) | "shared total" |
| **shape** $a$ | the first parameter of $R^2 \sim$ Beta($a$, 1): $k/2$, ½ or 1 | `r_squared_shape` (Stan), `prior_shape()` |
| **share** | the fraction of the total that one predictor gets; Beta(½, (k−1)/2) for every $a$ | $c^2$ (squared alignment, doc 2 section 7), $B$ (doc 2 section 8) |
| **noise units**, $f$ | an effect divided by the noise SD $\sigma$; $f^2$ is Cohen's $f^2$ | |
| **joint test** | the Bayes factor of the model with against without a term | |
| **conditional test** | the test of one coefficient's partial correlation | "per-coefficient test" |
| **block route** | a model's $k$ columns tested as one block against none: the whole-model Bayes factor for any $a$ | doc 2, section 7 |
| **flat partial r** | $a = k/2$ | "uniform correlations" (earlier versions), `stan_lm` with η = 1, LKJ(1) |
| **flat R** | $a = ½$ | "flat-R" (earlier versions) |
| **flat R-squared** | $a = 1$ | "flat R²" |
| **references** | JZS default (Cauchy, scale 0.354) and Zellner–Siow (scale 1) | |

---

## Status of each claim

Two labels are used. **Checked** means a numerical check passes on the current
code, in the settings the study covers; no claim is proven analytically here.
**Run** means a simulation study was run on the current code, and its result
is given. Details are in [`4_validation.md`](4_validation.md). All results are
from the full run of 4 November 2026 (commit 21d42dc), in which 137 of 138
checks pass. The one failure is a finding about flat R-squared (validation 10,
P4; see [`3_comparison.md`](3_comparison.md), section 5).

| Claim | Flat partial r | Flat R | Flat R-squared | Where |
|---|---|---|---|---|
| Valid correlation matrices give the three versions | checked (F1) | checked (F2) | checked (F3) | validation 00 |
| The Stan model samples the prior correctly | checked (01-P5–P7, KS ≤ .012) | checked (01-P1–P4, KS ≤ .009) | checked (01-P8–P10, KS ≤ .009) | validation 01, P1–P10 |
| The Stan posterior is calibrated ($k = 8$) | checked (SBC *p* ≥ .02) | checked (SBC *p* ≥ .38) | checked (SBC *p* ≥ .22) | validation 01, S |
| The whole-model Bayes factor is exact | checked (A1–A4) | checked (A1–A4) | checked (A1–A4) | validation 02 |
| Units and recoding don't change the answer | checked (A2) | checked (A2) | checked (A2) | validation 02 |
| Predictive matching (BF = 1 at $n = k + 1$) | checked (B1) | checked (B1) | checked (B1) | validation 02 |
| Information consistency, for $n \ge k + 3$ | checked (B2, B2b) | checked (B2, B2b) | checked (B2, B2b) | validation 02 |
| Intrinsic consistency | checked (B5) | checked (B5) | checked (B5) | validation 02 |
| Model-selection consistency, fixed $k$ | checked: $\tfrac q2 \log n$ in every case | checked: ½ log n against the intercept-only model, $\tfrac q2 \log n$ otherwise | checked: $\log n$ against the intercept-only model, $\tfrac q2 \log n$ otherwise (slopes 1.01, 0.96, 0.99 and 0.55, 1.55, 2.95) | validation 02, B3, B3b, B4 |
| Model-selection consistency, $k$ growing with $n$ | checked (G1, G2): evidence grows fast (log BF₀₁ 408 at $n = 2000$, $k = n/4$) | checked, but slow under the null (log BF₀₁ 1.9) | checked, but slow under the null (log BF₀₁ 3.7) | validation 08 |
| What singles the prior out | checked (C7) | checked (V5) | – | validation 05, 09 |
| The conditional test is computed correctly | checked (C1, C3–C5) | checked (C1, C2, C5) | checked (C1, C2, C5) | validation 05, C1–C5 |
| The block conditional test is computed correctly | checked (C9–C11) | checked (C9–C11) | checked (C9–C11) | validation 05, C9–C11 |
| Separating nulls from small effects, against JZS | run: equal or better (AUC +.001 to +.045, mean over $n$ and pattern per $k$ and correlation) | run: the same as flat partial r (conditional) | run: the same as flat partial r (conditional) | validation 04, 06 |
| Confirming a true null, one coefficient | run: 86–88%, for any $k$ | run: 55% at $k = 2$, 6% at $k = 20$ | run: 86% at $k = 2$, 26% at $k = 20$ | validation 06 |
| Calling a small effect (0.1) null | run: 52–60% | run: 19% at $k = 2$, 0.5% at $k = 20$ | run: 52% at $k = 2$, 4% at $k = 20$ | validation 06 |
| Calibrated on its own prior | checked (P1) | checked (P1) | checked (P1) | validation 10 |
| Variable selection | run: fewest false inclusions of the three, misses more (mean 0.36 of 8 nulls included, 2.09 of 4 real missed) | run: many more inclusions, fewer misses (1.14, 1.55) | run: nearly the same as flat R (1.13, 1.56) | validation 09 |
| Estimation: prediction error, interval coverage | run: better than OLS everywhere, coverage ≥ .916 | run: lowest error at small $n$, coverage down to .76 | run: in between, coverage down to .82 | validation 07 |

---

## What to watch out for when using the package

- **Choosing the prior.** Without `prior =`, the functions use flat partial
  r and say so once per session. What each prior assumes is described in "The three branches"; the
  recommendation (flat partial r) and its costs are in
  [`3_comparison.md`](3_comparison.md), section 8.
- **Shrinkage at small $n$.** With several predictors and a small sample, the
  prior pulls estimates towards zero. The example (`example.R`) prints this
  pull next to the ordinary estimates.
- **Flat R and flat R-squared make the conditional test cautious with many
  predictors.** These priors expect each coefficient to be a small share of
  the total, so "no effect" and "a tiny effect" predict nearly the same data.
  Under flat partial r each coefficient keeps a flat partial correlation,
  whatever $k$ is.
- **Leaving the null model.** Against the intercept-only model, evidence for
  the null grows like $a \log n$: ½ log n for flat R and log n for flat
  R-squared, whatever the number of predictors, and $\tfrac k2 \log n$ for
  flat partial r. How strict a variable selection is therefore depends partly
  on the prior over models, which should be chosen deliberately.
- **Under flat partial r the claimed $R^2$ grows with $k$.** The median is .63
  with 3 predictors and .87 with 10, so a large model with a modest fit is
  penalised.
- **Flat R-squared with one predictor is not Jeffreys' prior.** Its prior on
  $r$ has density $|r|$, which gives little weight to correlations near zero.
- **Information consistency needs $n \ge k + 3$**, for all three priors.
- **Main effects with interactions.** `bayes_factor_per_term()` removes one
  term at a time, so testing a main effect while its interaction stays in the
  model is subject to the usual marginality caveats.
- **Stan messages.** At the start of warmup, Stan may print a few
  "Informational Message ... proposal is about to be rejected" lines, because
  the sampler tried $R$ at its boundary. A few of these are harmless; check
  the divergent transitions that `print(fit)` reports.
- **Why the package does not use brms.** brms samples the coefficients
  directly. For $k \ge 2$ the prior has a sharp peak at all-zero coefficients
  that such a sampler misses without warning (doc 2, section 6).

---

## Open

- **Proofs.** The consistency criteria are checked numerically (studies 02
  and 08), in the settings those studies cover. None is proven analytically
  here, and non-nested models are not covered.
- **The empirical question between the branches:** how explained variance
  grows with the number of predictors in published regressions. Flat partial r
  predicts a median $R^2$ of .76 with 5 predictors; flat R predicts .25 and
  flat R-squared .50 at any size.
- **Benchmarks in variable selection** against the robust, hyper-g/n and
  intrinsic priors (BAS package), and a reanalysis of real data.
