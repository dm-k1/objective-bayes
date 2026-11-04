# 4. The three priors compared

**In one sentence:** the three priors that study 00 leaves (flat partial r,
flat R, flat R-squared) are compared on the criteria from the literature, each
with a column, and the choice between them is made from this table, not
before it.

The numbers come from the full validation run of 4 November 2026 (commit
e221561; [`5_validation.md`](5_validation.md)): 175 of 176
checks pass, and the one failure (10, P4, flat R-squared) is a finding,
explained in section 5. JZS (Cauchy, scale
0.354, the BayesFactor/JASP default) and Zellner–Siow (scale 1) are reference
lines throughout, not contenders.

---

## 1. What each prior stands for

All three start from the same principle: the weights follow what is possible
(every valid correlation matrix), with no scale chosen by convention. They
differ in which scale they keep flat as predictors are added.

| | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|
| `prior =` | `"flat_partial_r"` | `"flat_r"` | `"flat_r_squared"` |
| $R^2$ prior | Beta($k/2$, 1) | Beta(½, 1) | Beta(1, 1) |
| What stays flat as $k$ grows | each partial correlation | the total $R$ | the total $R^2$ |
| What moves as $k$ grows | the total grows (median 0.25, 0.63, 0.87 at $k$ = 1, 3, 10) | each partial correlation narrows | each partial correlation narrows, less than under flat R |
| At $k = 1$ | Jeffreys' flat correlation | Jeffreys' flat correlation | flat on $r^2$: little weight near $r = 0$ |
| Study 00, P(\|partial r\| > 0.3) at $k$ = 1 / 10 | 0.70 / 0.70 | 0.70 / 0.22 | 0.91 / 0.36 |
| Study 00, P($R^2$ > 0.5) at $k$ = 1 / 10 | 0.29 / 0.97 | 0.29 / 0.29 | 0.50 / 0.50 |

"Evenly" is always evenly *on some ruler*; which ruler is the judgement
([`1_from_behaviour.md`](1_from_behaviour.md), steps 1 and 2).

---

## 2. Correctness (does the code compute what the prior says?)

| Check | Flat partial r | Flat R | Flat R-squared | Study |
|---|---|---|---|---|
| Valid correlation matrices give the version | pass (F1, KS 0.003) | pass (F2, KS 0.002) | pass (F3, KS 0.003) | 00 |
| Stan samples the prior (KS < 0.03) | pass (P5–P7, KS ≤ 0.012) | pass (P1–P4, KS ≤ 0.009) | pass (P8–P10, KS ≤ 0.009) | 01, P1–P10 |
| Simulation-based calibration (Talts et al., 2018) | pass ($k = 8$, $n = 60$; χ² p ≥ 0.43) | pass (p ≥ 0.61) | pass (p ≥ 0.45) | 01, S |
| Whole-model BF against `integrate()` (162 settings) | pass, largest error 2.6e-07 | pass, largest error 4.1e-05 | pass, largest error 8.4e-05 | 02, A1 |
| Against a brute-force BF (no formula) | pass | pass | pass | 02, A3 |
| A second derivation | pass (block route, 5.1e-05) | pass (g-mixture and block route, 3.5e-04) | pass (block route and g-mixture by `integrate()`, 2.4e-05) | 02, A4 |
| Conditional test against `integrate()` | pass, 3.5e-06 | pass, 6.7e-05 | pass, 5.0e-05 | 05, C5 |
| Block conditional test against `integrate()` | pass (largest over the three 1.0e-04) | pass | pass | 05, C10 |

---

## 3. The criteria of Bayarri, Berger, Forte & García-Donato (2012)

| Criterion | Flat partial r | Flat R | Flat R-squared | Where |
|---|---|---|---|---|
| Basic (proper, no arbitrary scale) | met | met | met | doc 2, section 11 |
| Predictive matching ($n = k + 1$) | met, checked | met, checked | met, checked | 02, B1 |
| Measurement invariance | met, checked | met, checked | met, checked | 02, A2 |
| Group invariance | met | met | met | doc 2, section 11 |
| Information consistency ($n \ge k + 3$) | proof, checked | proof, checked | proof, checked | doc 3; 02, B2, B2b |
| Intrinsic consistency | proof | proof | proof | doc 3 |
| Model selection consistency, fixed $k$ | sketch, checked | sketch, checked | sketch, checked | doc 3; 02, B3, B3b, B4 |
| Null rate against the intercept-only model | $\tfrac q2 \log n$ (slopes 0.51, 1.46, 2.86 for $q$ = 1, 3, 6) | $\tfrac12 \log n$ (0.47, 0.50, 0.48) | $\log n$ (0.99, 1.00, 1.01) | 02, B3b |
| Null rate between models with an effect | $\tfrac q2 \log n$ (0.47, 1.48, 2.91) | $\tfrac q2 \log n$ (0.54, 1.56, 2.98) | $\tfrac q2 \log n$ (0.53, 1.64, 2.88) | 02, B3b |
| Model selection consistency, $k$ growing with $n$ (median log BF₀₁, $k = n/4$, $n = 2000$) | pass, 407.7 | pass, but slow: 1.9 | pass, but slow: 3.7 | 08 (JZS 149.3, Zellner–Siow 410.9) |

---

## 4. Testing one coefficient (operating characteristics)

Conditional tests unless stated; JZS and Zellner–Siow test jointly.

| | Flat partial r | Flat R | Flat R-squared | JZS | Zellner–Siow | Study |
|---|---|---|---|---|---|---|
| Separating nulls from small effects (AUC null vs 0.1, minus JZS), uncorrelated | +0.001 to +0.007 | the same | the same | 0 | −0.000 to +0.003 | 04, 06 |
| The same, correlation 0.7 | +0.031 to +0.045 | the same | the same | 0 | +0.012 to +0.027 | 06 |
| True null found (BF₀₁ > 3), $k = 2$ / $k = 20$ | 0.86 / 0.88 | 0.55 / 0.06 | 0.86 / 0.26 | 0.54 / 0.30 | 0.87 / 0.52 | 06 |
| Small effect (0.1) called null, $k = 2$ / $k = 20$ | 0.52 / 0.60 | 0.19 / 0.005 | 0.52 / 0.04 | 0.21 / 0.11 | 0.53 / 0.22 | 06 |
| Effect 0.2 found, $k = 2$ / $k = 20$ | 0.51 / 0.44 | 0.53 / 0.48 | 0.51 / 0.49 | 0.56 / 0.52 | 0.51 / 0.52 | 06 |
| False effects for a true null (BF₁₀ > 3), largest | 0.010 | 0.017 | 0.018 | 0.029 | 0.033 | 06 |
| Joint test: dependence on what the other predictors explain (M4; log BF per unit $R^2$ of the others) | −0.51 (conditional 0.002) | −5.19 (conditional −0.011) | −5.10 (conditional −0.001) | −1.93 | −0.52 | 04 |


What the rows above leave out is the third answer, "undecided". At $k = 20$
(conditional tests, share of datasets; full table in
[`5_validation.md`](5_validation.md), study 06):

| Truth | | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|---|
| null | says null / undecided / says effect | **0.88** / 0.11 / 0.01 | **0.06** / 0.92 / 0.02 | **0.26** / 0.72 / 0.02 |
| slope 0.1 | says null / undecided / says effect | 0.60 / 0.26 / **0.15** | 0.005 / 0.78 / **0.21** | 0.04 / 0.75 / **0.21** |
| slope 0.2 | says null / undecided / says effect | 0.31 / 0.26 / **0.44** | 0.00 / 0.52 / **0.48** | 0.00 / 0.52 / **0.49** |
| AUC null vs 0.1 (ranking) | | 0.701 | 0.701 | 0.701 |

The ranking is identical, so no prior is more accurate here. Flat R's low
rate of calling small effects null comes from rarely deciding at all (92% of
true nulls undecided), not from seeing effects better. Flat partial r decides
on most datasets and leans to the null when an effect is small.

---

## 5. Calibration of the Bayes factor

Is the share of true nulls among datasets with BF₁₀ equal to 1/(1 + BF₁₀),
with data drawn from the prior itself (and miscalibrated with data from
another prior, as it should be)?

| | Flat partial r | Flat R | Flat R-squared | Study |
|---|---|---|---|---|
| Calibrated on its own prior (χ² p, logistic slope) | pass: p 0.40–0.79, slope 0.98–1.03 | pass: p 0.52–0.98, slope 0.99–1.02 | pass: p 0.37–0.94, slope 0.98–1.02 | 10, P1 |
| Miscalibration shows with JZS data | yes: χ² 503–556, slope 0.89–0.95 | yes: χ² 233–1198, slope 1.16–1.41 | yes, least: χ² 50–265, slope 0.95–1.12 | 10, P2 |
| Conditional test against Savage–Dickey under the joint model (largest log difference) | pass: 0.092 (factor 1.10) | not defined (infinite prior density at 0) | **fail (P4)**: 0.874 (factor 2.40), at $k = 8$ | 10, P3–P5 |

---

## 6. Variable selection and estimation

| | Flat partial r | Flat R | Flat R-squared | JZS | Zellner–Siow | Study |
|---|---|---|---|---|---|---|
| Real predictors missed (of 4), $n = 500$, uncorrelated (beta-binomial / uniform model prior) | 0.84 / 0.64 | 0.45 / 0.43 | 0.44 / 0.43 | 0.58 / 0.47 | 0.89 / 0.67 | 09 |
| Null predictors included (of 8), the same | 0.13 / 0.20 | 0.91 / 0.60 | 0.92 / 0.61 | 0.42 / 0.47 | 0.12 / 0.17 | 09 |
| Ranking (AUC real vs null), difference from flat R, mean (range) over cells | +0.003 (−0.021 to +0.022) | 0 | −0.002 (−0.018 to +0.002) | −0.002 (−0.029 to +0.007) | +0.002 (−0.030 to +0.023) | 09 |
| Prediction error relative to OLS, $k = 8$, $n = 30$ (range over 4 cells) | 0.71–0.81 | 0.37–0.77 | 0.42–0.75 | 0.48–0.74 | 0.74–0.82 | 07 |
| 95% interval coverage, real effects, lowest cell (OLS 0.920) | 0.916 | 0.760 | 0.815 | 0.870 | 0.910 | 07 |
| Mean bias of real effects, toward zero (share of the effect) | −8% | −29% | −21% | −19% | −7% | 07 |

---

## 7. What the results show, per criterion

- **Correctness and the Bayarri et al. criteria:** no difference. All three
  pass every correctness check (00, 01, 02, 05) and meet predictive matching,
  invariance, information consistency ($n \ge k + 3$) and model selection
  consistency for fixed $k$. They differ in the null rate against the
  intercept-only model: $\tfrac k2 \log n$ for flat partial r, $\log n$ for flat
  R-squared, $\tfrac12 \log n$ for flat R.
- **$k$ growing with $n$ (08):** all three are consistent, but under the null
  flat R and flat R-squared gain evidence very slowly (log BF₀₁ 1.9 and 3.7 with
  500 null predictors at $n = 2000$), flat partial r fast (407.7).
- **One coefficient (04, 06):** the conditional tests rank datasets the same
  way for all three priors; they differ in the threshold. Flat partial r's
  answer does not depend on $k$ (true null found 0.86–0.88; a 0.1 effect
  called null 0.52–0.60). Flat R and flat R-squared confirm nulls less often as
  $k$ grows (0.55 → 0.06 and 0.86 → 0.26 from $k = 2$ to 20) and call a 0.1
  effect null less often too (0.19 → 0.005 and 0.52 → 0.04). False effects stay below 0.02 for all three. Their joint tests depend
  strongly on what the other predictors explain (M4: −5.2 and −5.1 per unit
  $R^2$, against −0.5 for flat partial r).
- **Calibration (10):** each is calibrated on its own prior. Coherence with
  the joint model holds for flat partial r (within a factor 1.10 of
  Savage–Dickey), not for flat R-squared at $k = 8$ (a factor 2.40; P4, the
  one failed check: its partial correlations are dependent through $R^2$, so
  the conditional test is not the Bayes factor of the joint model), and is
  undefined for flat R.
- **Variable selection (09):** ranking is nearly identical. Flat partial r
  includes few null predictors and misses more real ones, like Zellner–Siow;
  flat R and flat R-squared include more and miss fewer.
- **Estimation (07):** flat R and flat R-squared shrink more, so they have the
  lowest prediction error at small $n$ with many predictors, at the cost of
  bias (−29% and −21% of a real effect) and lower interval coverage (down to
  0.76 and 0.82). Flat partial r shrinks less (−8%), is better than OLS in
  every cell, and keeps coverage at 0.916 or above.

---

## 8. Recommendation

> **Left to the author.** This document reports the results; it does not pick
> a prior. Which one to use depends on which of the trade-offs above matters
> for the intended use (testing one coefficient, comparing whole models,
> estimation), and that judgement is the author's.

The package has no default on purpose: every function asks for `prior =`, so
whichever is chosen is stated in every analysis.
