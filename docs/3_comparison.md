# 3. The three priors compared

This document compares the three priors that study 00 leaves (flat partial
r, flat R and flat R-squared) on the criteria from the literature, with one
column per prior. The choice between them is made from these tables, not
before them; it is given in section 8.

The numbers come from the full validation run of 4 November 2026 (started
12:24, commit 21d42dc; [`4_validation.md`](4_validation.md)). Of its 138
checks, 137 pass. The one failure (study 10, check P4, flat R-squared) is a
finding about that prior and is explained in section 5. JZS (a Cauchy prior
with scale 0.354, the BayesFactor and JASP default) and Zellner–Siow (scale 1)
are shown throughout as references; they are not candidates.

---

## 1. What each prior keeps flat

All three priors start from the same principle: the prior weight follows
what is possible (every valid correlation matrix), with no scale chosen by
convention. They differ in which quantity they keep flat as predictors are
added.

| | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|
| `prior =` | `"flat_partial_r"` | `"flat_r"` | `"flat_r_squared"` |
| $R^2$ prior | Beta($k/2$, 1) | Beta(½, 1) | Beta(1, 1) |
| What stays flat as $k$ grows | each partial correlation | the total $R$ | the total $R^2$ |
| What moves as $k$ grows | the total grows (median .25, .63, .87 at $k$ = 1, 3, 10) | each partial correlation narrows | each partial correlation narrows, less than under flat R |
| At $k = 1$ | Jeffreys' flat correlation | Jeffreys' flat correlation | flat on $r^2$: little weight near $r = 0$ |
| Study 00, P(\|partial r\| > .30) at $k$ = 1 / 10 | .70 / .70 | .70 / .22 | .91 / .36 |
| Study 00, P($R^2$ > .50) at $k$ = 1 / 10 | .29 / .97 | .29 / .29 | .50 / .50 |

A prior that spreads its weight evenly always does so on some scale, and
which scale is the judgement ([`1_from_behaviour.md`](1_from_behaviour.md),
steps 1 and 2).

---

## 2. Correctness: does the code compute what the prior says?

| Check | Flat partial r | Flat R | Flat R-squared | Study |
|---|---|---|---|---|
| Valid correlation matrices give the version | pass (F1, KS .003) | pass (F2, KS .002) | pass (F3, KS .004) | 00 |
| Stan samples the prior (KS < .03) | pass (P5–P7, KS ≤ .012) | pass (P1–P4, KS ≤ .009) | pass (P8–P10, KS ≤ .009) | 01, P1–P10 |
| Simulation-based calibration (Talts et al., 2018) | pass ($k = 8$, $n = 60$; χ² *p* ≥ .02) | pass (*p* ≥ .38) | pass (*p* ≥ .22) | 01, S |
| Whole-model BF against `integrate()` (162 settings) | pass, largest error 2.6e-07 | pass, largest error 4.1e-05 | pass, largest error 8.4e-05 | 02, A1 |
| Against a brute-force BF (no formula) | pass | pass | pass | 02, A3 |
| A second derivation | pass (block route, 5.1e-05) | pass (g-mixture and block route, 3.5e-04) | pass (block route and g-mixture by `integrate()`, 2.4e-05) | 02, A4 |
| Conditional test against `integrate()` | pass, 3.5e-06 | pass, 6.7e-05 | pass, 5.0e-05 | 05, C5 |
| Block conditional test against `integrate()` | pass (largest over the three 1.0e-04) | pass | pass | 05, C10 |

KS is the Kolmogorov–Smirnov distance between the simulated and the stated
distribution. The closest call is flat partial r's calibration of the noise
SD (χ² *p* = .02, against a threshold of .01).

---

## 3. The criteria of Bayarri, Berger, Forte and García-Donato (2012)

| Criterion | Flat partial r | Flat R | Flat R-squared | Where |
|---|---|---|---|---|
| Basic (proper, no arbitrary scale) | met | met | met | doc 2, section 11 |
| Predictive matching ($n = k + 1$) | met, checked ($k$ = 1–6) | met, checked ($k$ = 1–6) | met, checked ($k$ = 2–6) | 02, B1 |
| Measurement invariance | met, checked | met, checked | met, checked | 02, A2 |
| Group invariance | met | met | met | doc 2, section 11 |
| Information consistency ($n \ge k + 3$) | checked | checked | checked | 02, B2, B2b |
| Intrinsic consistency | checked | checked | checked | 02, B5 |
| Model selection consistency, fixed $k$ | checked | checked | checked | 02, B3, B3b, B4 |
| Null rate against the intercept-only model | $\tfrac k2 \log n$ (slopes 0.51, 1.44, 2.90 for $k$ = 1, 3, 6) | $\tfrac12 \log n$ (0.51, 0.48, 0.50) | $\log n$ (1.01, 0.96, 0.99) | 02, B3b |
| Null rate between models with an effect | $\tfrac q2 \log n$ (0.50, 1.47, 2.95) | $\tfrac q2 \log n$ (0.56, 1.58, 2.97) | $\tfrac q2 \log n$ (0.55, 1.55, 2.95) | 02, B3b |
| Model selection consistency, $k$ growing with $n$ (median log BF₀₁, $k = n/4$, $n = 2000$) | pass, 407.6 | pass, but slow: 1.9 | pass, but slow: 3.7 | 08 (JZS 149.3, Zellner–Siow 410.8) |

"Checked" means a numerical check in the settings the study covers; none of
these criteria is proven analytically here. The **null rate** (= how fast the
evidence for a true null grows with the sample size) is the slope of the
median log BF₀₁ against $\log n$. Against the intercept-only model it is
$a \log n$, which is $\tfrac k2 \log n$ for flat partial r; between two models
that differ by $q$ predictors it is $\tfrac q2 \log n$ for all three.

---

## 4. Testing one coefficient: operating characteristics

The rows use the **conditional tests** (= the test of one coefficient's
partial correlation against zero, with the prior the joint prior implies for
it) unless stated otherwise; JZS and Zellner–Siow test jointly. Slopes are in
units of the noise SD, with standardised predictors. The study 06 rows come
from 72 designs ($k$ = 2, 5, 10, 20; $n$ = 50, 200, 1,000; predictor
correlation 0, .30, .70; sparse or dense other coefficients), each at three
slopes of $x_1$ (0, 0.1, 0.2) with 200 datasets each; the rates are averaged
over $n$, correlation and pattern.

| | Flat partial r | Flat R | Flat R-squared | JZS | Zellner–Siow | Study |
|---|---|---|---|---|---|---|
| Separating nulls from small effects (AUC null vs 0.1, minus JZS), uncorrelated | +.001 to +.007 | the same | the same | 0 | −.000 to +.003 | 04, 06 |
| The same, correlation .70 | +.031 to +.045 | the same | the same | 0 | +.012 to +.027 | 06 |
| True null found (BF₀₁ > 3), $k = 2$ / $k = 20$ | .86 / .88 | .55 / .06 | .86 / .26 | .54 / .30 | .87 / .52 | 06 |
| Small effect (0.1) called null, $k = 2$ / $k = 20$ | .52 / .60 | .19 / .005 | .52 / .04 | .21 / .11 | .53 / .22 | 06 |
| Effect 0.2 found, $k = 2$ / $k = 20$ | .51 / .44 | .53 / .48 | .51 / .49 | .56 / .52 | .51 / .52 | 06 |
| False effects for a true null (BF₁₀ > 3), largest | .010 | .017 | .018 | .029 | .033 | 06 |
| Joint test: dependence on what the other predictors explain (M4; log BF per unit $R^2$ of the others) | −0.51 (conditional 0.002) | −5.19 (conditional −0.011) | −5.10 (conditional −0.001) | −1.93 | −0.52 | 04 |

These rows count only the correct and the wrong answers. A test can also stay
undecided, and that third answer is where the priors differ most. At $k = 20$
(conditional tests, share of datasets; full table in
[`4_validation.md`](4_validation.md), study 06):

| Truth | | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|---|
| null | says null / undecided / says effect | **.88** / .11 / .01 | **.06** / .92 / .01 | **.26** / .72 / .02 |
| slope 0.1 | says null / undecided / says effect | .60 / .26 / **.15** | .005 / .78 / **.21** | .04 / .75 / **.21** |
| slope 0.2 | says null / undecided / says effect | .31 / .26 / **.44** | .00 / .52 / **.48** | .00 / .52 / **.49** |
| AUC null vs 0.1 (ranking) | | .701 | .701 | .701 |

Correct answers are in bold. The three conditional tests rank the datasets
identically, so none of them separates nulls from small effects better.
Flat R rarely calls a small effect null because it rarely calls anything
null: it leaves 92% of true nulls undecided. It finds a 0.1 effect about as
often as flat R-squared (.21) and JZS (.24), and more often than flat
partial r (.15). Flat partial r decides on most datasets and leans towards
the null when an effect is small.

---

## 5. Calibration and coherence of the Bayes factor

A Bayes factor is calibrated when, with data drawn from the prior itself, the
share of true nulls among datasets with a given BF₁₀ equals 1/(1 + BF₁₀). With
data from another prior it should be miscalibrated.

| | Flat partial r | Flat R | Flat R-squared | Study |
|---|---|---|---|---|
| Calibrated on its own prior (χ² *p*, logistic slope) | pass: *p* = .40–.79, slope 0.98–1.03 | pass: *p* = .52–.98, slope 0.99–1.02 | pass: *p* = .37–.94, slope 0.98–1.02 | 10, P1 |
| "*x*₁'s partial correlation = 0" is a proper restriction of the joint prior (density of ρ₁ near 0, ε = 0.1 → 0.005, $k$ = 3 / 8) | yes (P6): 0.50 → 0.49 / 0.50 → 0.49 | no (P7): 1.51 → 2.99 / 2.23 → 5.20, others' median $R^2$ → .001 / .003 | yes (P6): 0.74 → 0.78 / 1.34 → 1.56 | 10, P6–P7 |
| Miscalibration shows with JZS data | yes: χ² = 503–556, slope 0.89–0.95 | yes: χ² = 233–1198, slope 1.16–1.41 | yes, least: χ² = 50–265, slope 0.95–1.12 | 10, P2 |
| Conditional test against Savage–Dickey under the joint model (largest log difference, on the datasets where Savage–Dickey is estimable) | pass: 0.092 (factor 1.10; 16 of 20 datasets, the other 4 within 3 Monte Carlo SE) | no joint counterpart (P7): the prior density of one partial correlation at 0 is infinite, and fixing it at 0 drags the other coefficients to 0 | **fail (P4)**: 0.874 (factor 2.40), at $k = 8$; 19 of 20 datasets | 10, P3–P5, P7 |

The **Savage–Dickey ratio** (= the posterior density of a parameter at zero
divided by its prior density there) gives the Bayes factor of the full model
against the same model with that parameter at zero. The row on the
restriction shows whether that restricted model is a proper alternative.
Under flat partial r and flat R-squared the density of ρ₁ near zero settles
as the interval shrinks. Under flat R it keeps growing, and conditioning on
ρ₁ near zero drags the other coefficients to zero, so the restricted model
collapses to the intercept-only model (study 10, P7).

---

## 6. Variable selection and estimation

| | Flat partial r | Flat R | Flat R-squared | JZS | Zellner–Siow | Study |
|---|---|---|---|---|---|---|
| Real predictors missed (of 4), $n = 500$, uncorrelated (beta-binomial / uniform model prior) | 0.84 / 0.64 | 0.45 / 0.43 | 0.44 / 0.43 | 0.58 / 0.47 | 0.89 / 0.67 | 09 |
| Null predictors included (of 8), the same | 0.13 / 0.20 | 0.91 / 0.60 | 0.92 / 0.61 | 0.42 / 0.47 | 0.12 / 0.17 | 09 |
| Ranking (AUC real vs null), difference from flat R, mean (range) over cells | +.003 (−.021 to +.022) | 0 | −.002 (−.018 to +.002) | −.002 (−.029 to +.007) | +.002 (−.030 to +.023) | 09 |
| Prediction error relative to OLS, $k = 8$, $n = 30$ (range over 4 cells) | 0.72–0.81 | 0.37–0.77 | 0.42–0.75 | 0.48–0.74 | 0.74–0.82 | 07 |
| 95% interval coverage of real effects, each method's lowest cell (OLS's own lowest cell: .920) | .916 | .760 | .815 | .870 | .910 | 07 |
| Mean bias of real effects towards zero (share of the effect) | 8% | 28% | 21% | 19% | 7% | 07 |

---

## 7. What the results show, per criterion

- **Correctness and the Bayarri et al. criteria.** The priors do not differ
  here. All three pass every correctness check (studies 00, 01, 02 and 05) and
  meet predictive matching, invariance, information consistency
  ($n \ge k + 3$) and model selection consistency for fixed $k$. They differ
  only in the null rate against the intercept-only model: $\tfrac k2 \log n$
  for flat partial r, $\log n$ for flat R-squared and $\tfrac12 \log n$ for
  flat R.
- **$k$ growing with $n$ (study 08).** All three are consistent. Under a true
  null, flat R and flat R-squared gain evidence very slowly (log BF₀₁ of 1.9
  and 3.7 with 500 null predictors at $n = 2000$); flat partial r gains it
  quickly (407.6).
- **One coefficient (studies 04 and 06).** The conditional tests rank
  datasets the same way for all three priors and differ only in where they
  decide. Flat partial r confirms true nulls at the same rate at every $k$
  (.86 at $k = 2$, .88 at $k = 20$); its rates for small effects shift
  somewhat (0.1 effects called null from .52 to .60, 0.2 effects found from
  .51 to .44). As $k$ grows from 2
  to 20, flat R and flat R-squared confirm fewer nulls (from .55 to .06, and
  from .86 to .26) and call fewer 0.1 effects null (from .19 to .005, and from
  .52 to .04), because they stay undecided more often. False effects stay
  below .02 for all three. The joint tests of flat R and flat R-squared also
  separate nulls from small effects worse than JZS when predictors are
  correlated, and depend strongly on what the other predictors explain (M4: −5.2 and −5.1 per unit $R^2$, against −0.5 for flat
  partial r).
- **Calibration and coherence (study 10).** Each prior is calibrated on data
  from its own prior. The conditional test agrees with the joint model under
  flat partial r (within a factor 1.10 of Savage–Dickey on the 16 of 20
  datasets where it is estimable). Under flat R-squared
  it does not at $k = 8$ (a factor 2.40; P4, the one failed check), because the
  partial correlations are dependent through $R^2$, so the conditional test
  is not the Bayes factor of the joint model. Under flat R the restricted
  model is not a proper alternative: the prior density of one partial
  correlation is infinite at zero, and conditioning on it being near zero
  drags the other coefficients to zero (study 10, P7: the others' median
  $R^2$ falls to .001 at $k = 3$), so the restricted model collapses to the
  intercept-only model.
- **Variable selection (study 09).** The ranking of predictors is nearly
  identical. Flat partial r, like Zellner–Siow, includes few null predictors
  and misses more real ones; flat R and flat R-squared include more and miss
  fewer.
- **Estimation (study 07).** Flat R and flat R-squared shrink more, which
  gives them the lowest prediction error at small $n$ with many predictors.
  The cost is bias (28% and 21% of a real effect, towards zero) and lower
  interval coverage (down to .76 and .82). Flat partial r shrinks less (8%),
  is better than OLS in every cell, and keeps coverage at .916 or above.

---

## 8. Recommendation

> **Recommendation (the author's judgement, from the results above):**
> **flat partial r.** It is the only one of the three whose test of one
> coefficient confirms true nulls at the same rate at every number of
> predictors, and the only one whose conditional test agrees with one joint
> model. Sections 2–6 are the evidence; this section is the conclusion drawn
> from them.

### Why flat partial r

- **Its test of one coefficient confirms true nulls at the same rate at
  every $k$.** It finds .86 of true nulls at $k = 2$ and .88 at $k = 20$
  (study 06); its rates for small effects shift somewhat (0.1 effects called
  null from .52 to .60). Under flat R (from .55 to .06) and flat R-squared
  (from .86 to .26), adding predictors changes the verdict about $x_1$ while
  the data on $x_1$ stay the same.
- **It is coherent.** Its conditional test agrees with the Bayes factor of the
  joint model, within a factor 1.10 of Savage–Dickey on the 16 of 20 datasets
  where Savage–Dickey is estimable (study 10). Flat R-squared fails this check
  at $k = 8$ (a factor 2.40). Under flat R the restricted model is not a
  proper alternative: fixing one coefficient's partial correlation near zero
  drags the other coefficients to zero as well (study 10, P7).
- **Its joint test barely depends on what the other predictors explain**
  (−0.5 per unit $R^2$, against −5.2 and −5.1; study 04).
- **Evidence for a true null keeps growing in large models** (log BF₀₁ of
  407.6 with 500 null predictors at $n = 2000$, against 1.9 and 3.7; study 08).
- **Its estimates stay close to the truth**: 8% bias towards zero, 95%
  interval coverage of .916 or above, and better than OLS in every cell
  (study 07). Flat R has 28% bias and coverage as low as .760.
- **The other two priors are not more accurate.** All three conditional tests
  rank the datasets identically (AUC .701 at $k = 20$). Flat R's low rate of
  calling small effects null comes from rarely deciding (92% of true nulls
  undecided at $k = 20$; section 4).

### Why a stable rate of confirming nulls is preferable (the author's argument)

This is an argument, not a result of the studies.

- A default test should let the conclusion about $x_1$ depend on the data
  about $x_1$. Under flat R and flat R-squared, the same data about $x_1$ lead
  to a different verdict when covariates are added.
- Their extra caution is not learned from the data. The three conditional
  tests rank datasets identically, so the larger undecided share comes from
  the prior, which shrinks each partial correlation as $k$ grows (study 00),
  and not from weaker evidence.
- A researcher who adds covariates would, under flat R, nearly lose the
  ability to confirm a null (.06 at $k = 20$), with nothing in the data to
  show for it.
- Flat partial r's own cost, its lean towards the null for small effects, is
  visible and can be reported with the support curve, and only its
  conditional test agrees with one joint model (by conditioning, the notion
  the Savage–Dickey ratio assumes; see `relation_to_others.md`, "What is new
  here, and what is not", for the other notions of compatible priors).

### What it costs (to be stated when it is used)

- **It leans towards the null for small effects.** At $k = 20$ it calls 60%
  of 0.1 effects and 31% of 0.2 effects null (study 06, averaged over
  $n$ = 50 to 1,000).
  This follows from the prior, which keeps real weight near zero at every $k$,
  so a weak signal reads as "probably nothing". A BF₀₁ above 3 from this prior
  therefore does not rule out a small effect; report it together with the
  support curve ([`R/support_curve.R`](https://github.com/dm-k1/flatcor/blob/main/R/support_curve.R) in flatcor), which shows the effect sizes the data
  cannot distinguish from zero.
- **It shrinks less, so it predicts somewhat worse** with many predictors and
  a small sample (prediction error 0.72–0.81 times that of OLS at $k = 8$ and
  $n = 30$, against 0.37–0.77 for flat R; study 07). For prediction alone,
  flat R does better.
- **It is conservative in variable selection**: it includes few null
  predictors and misses more real ones (0.84 of 4 missed at $n = 500$ with the
  beta-binomial model prior, against 0.45 for flat R; study 09), much like
  Zellner–Siow.
- **Against the intercept-only model its evidence depends on $k$.** For a
  true null it grows like $\tfrac k2 \log n$ (section 3), because with many
  predictors the prior expects a large total (median $R^2$ of .87 at
  $k = 10$). A model whose predictors together explain little is judged
  against that expectation.

### What stays the same

Following this recommendation, flat partial r is the package's default; a
function that uses it because no prior was given says so once per session,
so the choice stays visible. Flat R and flat R-squared remain
available for uses where their trade-off fits, such as prediction with many
predictors and a small sample.
