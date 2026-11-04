# 4. The three priors compared

**In one sentence:** the three priors that study 00 leaves (flat partial r,
flat R, flat R-squared) are compared on the criteria from the literature, each
with a column, and the choice between them is made from this table, not
before it.

The numbers come from the validation reruns for three priors
([`../validation/README.md`](../validation/README.md)). Cells marked
"(rerun pending)" are filled when that study's rerun is in. JZS (Cauchy, scale
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
| Stan samples the prior (KS < 0.03) | (rerun pending) | (rerun pending) | (rerun pending) | 01, P1–P10 |
| Simulation-based calibration (Talts et al., 2018) | (rerun pending) | (rerun pending) | (rerun pending) | 01, S |
| Whole-model BF against `integrate()` (162 settings) | pass, largest error 2.6e-07 | pass, largest error 4.1e-05 | (rerun pending) | 02, A1 |
| Against a brute-force BF (no formula) | pass | pass | (rerun pending) | 02, A3 |
| A second derivation | pass (block route, 5.1e-05) | pass (g-mixture and block route, 3.5e-04) | (rerun pending) | 02, A4 |
| Conditional test against `integrate()` | (rerun pending) | (rerun pending) | (rerun pending) | 05, C5 |
| Block conditional test against `integrate()` | (rerun pending) | (rerun pending) | (rerun pending) | 05, C10 |

---

## 3. The criteria of Bayarri, Berger, Forte & García-Donato (2012)

| Criterion | Flat partial r | Flat R | Flat R-squared | Where |
|---|---|---|---|---|
| Basic (proper, no arbitrary scale) | met | met | met | doc 2, section 11 |
| Predictive matching ($n = k + 1$) | met, checked | met, checked | met, (rerun pending) | 02, B1 |
| Measurement invariance | met, checked | met, checked | met, (rerun pending) | 02, A2 |
| Group invariance | met | met | met | doc 2, section 11 |
| Information consistency ($n \ge k + 3$) | proof, checked | proof, checked | proof, (rerun pending) | doc 3; 02, B2, B2b |
| Intrinsic consistency | proof | proof | proof | doc 3 |
| Model selection consistency, fixed $k$ | sketch, checked | sketch, checked | sketch, (rerun pending) | doc 3; 02, B3, B3b, B4 |
| Null rate against the intercept-only model | $\tfrac q2 \log n$ (slopes 0.50, 1.47, 2.81 for $q$ = 1, 3, 6) | $\tfrac12 \log n$ (0.50, 0.49, 0.51) | $\log n$ expected; (rerun pending) | 02, B3b |
| Null rate between models with an effect | $\tfrac q2 \log n$ (0.47, 1.49, 2.93) | $\tfrac q2 \log n$ (0.54, 1.55, 2.91) | $\tfrac q2 \log n$ expected; (rerun pending) | 02, B3b |
| Model selection consistency, $k$ growing with $n$ (median log BF₀₁, $k = n/4$, $n = 2000$) | (rerun pending) | (rerun pending) | (rerun pending) | 08 |

---

## 4. Testing one coefficient (operating characteristics)

Conditional tests unless stated; JZS and Zellner–Siow test jointly.

| | Flat partial r | Flat R | Flat R-squared | JZS | Zellner–Siow | Study |
|---|---|---|---|---|---|---|
| Separating nulls from small effects (AUC null vs 0.1, minus JZS), uncorrelated | (rerun pending) | (rerun pending) | (rerun pending) | 0 | (rerun pending) | 04, 06 |
| The same, correlation 0.7 | (rerun pending) | (rerun pending) | (rerun pending) | 0 | (rerun pending) | 06 |
| True null found (BF₀₁ > 3), $k = 2$ / $k = 20$ | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 06 |
| Small effect (0.1) called null, $k = 2$ / $k = 20$ | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 06 |
| Effect 0.2 found, $k = 2$ / $k = 20$ | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 06 |
| False effects for a true null (BF₁₀ > 3), largest | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 06 |
| Joint test: dependence on what the other predictors explain (M4) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 04 |

---

## 5. Calibration of the Bayes factor

Is the share of true nulls among datasets with BF₁₀ equal to 1/(1 + BF₁₀),
with data drawn from the prior itself (and miscalibrated with data from
another prior, as it should be)?

| | Flat partial r | Flat R | Flat R-squared | Study |
|---|---|---|---|---|
| Calibrated on its own prior (χ² p, logistic slope) | (rerun pending) | (rerun pending) | (rerun pending) | 10, P1 |
| Miscalibration shows with JZS data | (rerun pending) | (rerun pending) | (rerun pending) | 10, P2 |
| Conditional test against Savage–Dickey under the joint model (largest log difference) | (rerun pending) | (rerun pending) | (rerun pending) | 10, P3–P5 |

---

## 6. Variable selection and estimation

| | Flat partial r | Flat R | Flat R-squared | JZS | Zellner–Siow | Study |
|---|---|---|---|---|---|---|
| Real predictors missed (of 4), $n = 500$, uncorrelated | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 09 |
| Null predictors included (of 8), the same | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 09 |
| Ranking (AUC real vs null), difference from flat R | (rerun pending) | 0 | (rerun pending) | (rerun pending) | (rerun pending) | 09 |
| Prediction error relative to OLS, $k = 8$, $n = 30$ | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 07 |
| 95% interval coverage, real effects, lowest cell | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 07 |
| Mean bias of real effects | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | 07 |

---

## 7. Recommendation

To be written from the results above, once every "(rerun pending)" cell is
filled. It will say which prior the results favour for which use (testing one
coefficient, comparing whole models, estimation), and what that choice costs.

The choice is the author's. The package has no default on purpose: every
function asks for `prior =`, so whichever is chosen is stated in every analysis.
