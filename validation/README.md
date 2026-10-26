# Validation

Three questions, one script each. Run from the project folder:

```bash
Rscript validation/01_prior_and_posterior.R   # needs cmdstanr (or rstan)
Rscript validation/02_bayes_factor.R          # no Stan needed
Rscript validation/03_simulations.R           # a few minutes
Rscript validation/04_matched_comparison.R    # needs the results of 03
```

The scripts load the package from this folder with `pkgload::load_all()`.

| Script | Question | How | Run time |
|---|---|---|---|
| `01_prior_and_posterior.R` | Does the Stan model (`inst/stan/flat_r.stan`) get the prior and the posterior right? | P1–P4: sample the prior alone; it must give a flat correlation (k = 1), R² ~ Beta(½, 1) (k = 3 correlated, k = 8), and the same per-predictor shares as the drawing recipe. S: simulation-based calibration of the posterior, k = 8, n = 60, 400 datasets. Contrast: the same prior sampled in coefficient coordinates | about 5 minutes |
| `02_bayes_factor.R` | Is the Bayes factor right? | A1: the grid against `integrate()` over a sweep of R², n, k. A2: invariance to units and recoding. A3: against a brute-force Bayes factor that draws from the prior and uses no formula. A4: against the same prior written as a mixture of g-priors. B1–B3: predictive matching, information and null consistency | about a minute |
| `03_simulations.R` | How does it behave across many datasets, compared with JZS, per-coefficient flat partial correlation and R2-D2? | 4 correlated predictors, 200 replications per cell, closed-form Bayes factors | long |
| `04_matched_comparison.R` | Is the difference in 03 discrimination or calibration? | AUC (threshold-free) of BF01 for null vs small effect, paired bootstrap flat-R vs JZS, matched thresholds, and why: dependence of the evidence on what the other predictors explain | about 2 minutes |
| `05_conditional_test.R` | Experiment: a per-coefficient test that uses only that coefficient's partial correlation, with the prior the joint flat-R prior implies for it | checks C1–C3, then AUC and calibration on the design of 03 | about 3 minutes |
| `06_broad_simulation.R` | Does the picture hold beyond one design? | k = 2–20, n = 50–1000, correlation 0–0.7, sparse and dense other coefficients; JZS, joint, conditional, flat partial r | about 15 minutes |
| `07_estimation.R` | How good are the estimates? | bias, error relative to ordinary regression, coverage of 95% intervals; flat-R vs JZS vs `lm()`; 2,400 Stan fits, resumes if interrupted | 1.5–3 hours |
| `08_growing_k.R` | Consistency when k grows with n (k = n/4, k = √n) | null true and small model true; flat-R vs Zellner–Siow vs JZS | 10–20 minutes |
| `09_variable_selection.R` | Choosing among 4,096 subsets of 12 candidates (4 real): does the prior's claim about model size matter? | flat-R, JZS, Zellner–Siow, BIC; uniform and beta-binomial model priors; true model found, false inclusions, missed | about 10 minutes |
| `conditional_flat_r.R` | helper: the per-coefficient test, used by 06 | – | – |
| `comparison_methods.R` | helper: the comparison methods | – | – |
| `settings.R` | helper: parallel cores, caching, `record_check()`, `save_figure()` | – | – |

Package tests (`devtools::test()`): 22 pass on cmdstanr.

Each check prints PASS/FAIL with its criterion; 01 and 02 stop with an error if
anything failed. Simulation results and figures are cached in
`validation/results/` (not tracked by git).

**Simulation-based calibration** (01, S): draw the true parameters from the
prior, simulate a dataset, fit it, and record where the true value falls among
the posterior draws. If the sampler is right, every position is equally likely,
so the counts per bin should be flat.

## Last results

| Check | Criterion | Last result |
|---|---|---|
| 01 P1–P4 | KS distance < 0.03 | pass (cmdstanr, current model): k = 1 0.006, k = 3 0.006, k = 8 0.008, shares 0.004–0.009 |
| 01 S | rank bins, χ² p > 0.01 | pass (cmdstanr, current model): R² p = 0.65, β₁ 0.73, σ 0.52. Contrast, coefficient coordinates: all p < 0.001; R² bins 38 15 24 23 23 18 24 19 26 **190**, σ bins **191** 15 17 … (40 expected per bin) |
| 02 A1 | log BF within 0.001 of `integrate()`, 162 settings | pass, largest error 4.1e-05 |
| 02 A2 | identical up to 1e-8 | pass |
| 02 A3 | within 0.05 of brute force (400,000 draws) | pass: 1.644 vs 1.638 (k = 3), −0.692 vs −0.691 (k = 1), 2.507 vs 2.518 (k = 5) |
| 02 A4 | within 0.002 of the Bayes factor computed as a mixture of g-priors | pass, largest difference 3.5e-04 |
| 02 B1 | BF within 0.001 of 1, k = 1–6 | pass |
| 02 B2 | increasing; log BF > 50 at R² = 0.9999 | pass (12.6, 28.6, 44.7, 60.9) |
| 02 B2b | bounded at n = k + 2, growing at n = k + 3 (the analytic boundary) | pass |
| 02 B3 | median log BF₀₁ increasing, slope vs log n in 0.35–0.65 | pass (slope 0.53) |
| 02 B4 | effect x1 = 0.15: median log BF₁₀ increasing, > 20 at n = 5000 | pass (0.9, 9.3, 54.1 at n = 200, 1000, 5000) |

**Why the Stan model samples R and a direction.** The prior as a density on the
coefficients has a sharp peak at zero (`docs/2_how_it_calculates.md`, section
3). A sampler moving through the coefficients, as brms does with a custom prior,
misses it for larger k, without divergent transitions: with brms at k = 8 the
10% quantile of the implied R² was 0.027 where the prior gives 0.010 (KS
distance 0.074; `archive/check_brms_prior.R`). Check 01 shows the contrast.

## 03: how the methods compare (last run: October 2026, 200 datasets per cell)

The null predictor x3, n = 200, P(BF₀₁ > 3), i.e. the true null found:

| Method | predictors uncorrelated | correlated (0.5) |
|---|---|---|
| flat partial correlation, per coefficient | 0.91 | 0.90 |
| JZS default | 0.61 | 0.65 |
| **joint flat-R** | 0.40 | 0.54 |
| R2-D2 matched | 0.16 | 0.12 |

A small real effect (x1 = 0.1) called null, P(BF₀₁ > 3):

| Method | n = 200, uncorrelated | n = 200, correlated | n = 500, uncorrelated | n = 500, correlated |
|---|---|---|---|---|
| flat partial correlation | 0.62 | 0.64 | 0.34 | 0.52 |
| JZS default | 0.22 | 0.26 | 0.13 | 0.27 |
| **joint flat-R** | 0.08 | 0.15 | 0.08 | 0.22 |
| R2-D2 matched | 0.00 | 0.00 | 0.00 | 0.02 |

- For real effects of 0.2 or more, flat-R and JZS are practically identical.
- False effects for a true null (x1 = 0) are rare for every method (≤ 5%).
- **Trade-off figure** (`docs/figures/7_trade_off.png`): plotted as "nulls found"
  against "small effects called null", the flat-R curve lies left of JZS's. But
  check 04 shows this is **calibration, not discrimination**: where each
  method's Bayes factors sit relative to 3, not how well they separate nulls
  from small effects (see below).
- Per-coefficient flat partial correlations find nulls fastest but call small
  real effects null in 60–70% of datasets: a flat prior on each coefficient
  puts too much weight on large effects.

## 04: compared at equal error rates

**Discrimination** (AUC: how often BF₀₁ is larger for a null dataset than for a
small effect of 0.1; threshold-free), flat-R minus JZS, paired bootstrap 95% CI:

| | n = 100 | n = 200 | n = 500 |
|---|---|---|---|
| predictors uncorrelated | −0.010 [−0.034, 0.012] | −0.006 [−0.019, 0.007] | −0.004 [−0.009, 0.001] |
| correlated (0.5) | **−0.039** [−0.058, −0.020] | **−0.034** [−0.049, −0.019] | **−0.028** [−0.036, −0.021] |

- Uncorrelated predictors: all four methods discriminate equally well. The
  difference in 03 is only where the Bayes factors sit relative to 3.
- Correlated predictors: flat-R discriminates slightly but clearly **worse** than
  JZS (and than R2-D2 and per-coefficient flat partial correlations).
- **Why** (M4, correlated, n = 200): at a fixed t-statistic for x1, flat-R's log
  BF for x1 drops by 5.2 per unit of R² explained by the other predictors,
  against 1.9 for JZS. Its evidence for x1 depends 2.7 times as strongly on
  what the rest of the model explains: the shared budget at work. With
  correlated predictors, that R² partly contains x1's own effect, so the
  dependence adds noise (rank correlation of log BF with t²: 0.78 vs 0.93).
- This is the cost of the budget for testing **one coefficient** in a model with
  correlated predictors. It does not affect comparing whole models.

## 05: a per-coefficient test (experiment)

**In one sentence:** test one coefficient's partial correlation only, with the
prior the joint flat-R prior implies for it, so the multiplicity correction
stays but the evidence no longer depends on what the other predictors explain.

The implied prior per coefficient: $f_j^2 = \text{share} \cdot R^2/(1-R^2)$,
$R$ flat, share ~ Beta(½, (k−1)/2). It depends only on k, not on how the
predictors correlate (C2, KS 0.002). The other coefficients get flat priors, as
in Kucharský et al. (2026), so the Bayes factor depends only on the sample
partial correlation, n and k.

Checks: C1 equals the joint test for k = 1 (difference 3e-12); C2 as above;
C3 with a flat prior it matches the approximate flat partial r (within 0.012).

Discrimination (AUC, null vs 0.1), paired bootstrap 95% CI, design of 03:

| | n = 100 | n = 200 | n = 500 |
|---|---|---|---|
| conditional − JZS, correlated (0.5) | **+0.052** [0.030, 0.074] | **+0.034** [0.022, 0.046] | **+0.025** [0.018, 0.034] |
| conditional − joint, correlated (0.5) | **+0.113** [0.069, 0.150] | **+0.078** [0.054, 0.102] | **+0.053** [0.039, 0.068] |
| conditional − JZS, uncorrelated | −0.004 [−0.022, 0.012] | +0.004 [−0.008, 0.016] | −0.002 [−0.007, 0.002] |

Calibration at 3 (correlated, n = 200): nulls found 0.43 (JZS 0.55, joint 0.22);
slope 0.1 called null 0.21 (JZS 0.27); slope 0.3 found 0.94 (JZS 0.94).

So in this design the conditional test separates nulls from small effects
**better than JZS** with correlated predictors and equally well without, while
keeping the flat-R expectation per coefficient. Study 06 checks whether this
holds across designs. Not yet covered: factors and other multi-column terms
(the block version would use a share ~ Beta(q/2, (k−q)/2) for q columns).

## 06: across designs (k = 2–20, n = 50–1000, correlation 0–0.7, 72 designs)

**Discrimination** (AUC null vs 0.1, minus JZS, mean over n and pattern of the
other coefficients):

| k | correlation 0 | 0.3 | 0.7 |
|---|---|---|---|
| conditional flat-R: 2 / 5 / 10 / 20 | +0.001 / +0.004 / +0.002 / +0.007 | +0.009 / +0.011 / +0.013 / +0.012 | **+0.035 / +0.045 / +0.041 / +0.031** |
| joint flat-R: 2 / 5 / 10 / 20 | −0.015 / −0.011 / −0.007 / −0.005 | −0.022 / −0.020 / −0.009 / −0.008 | −0.064 / −0.035 / −0.022 / −0.013 |

- The conditional test is at least as good as JZS in every group (worst single
  cell −0.013, noise level), and better with correlated predictors.
- The conditional test and the per-coefficient flat partial r have **identical**
  AUCs: both depend only on the coefficient's partial correlation, so they rank
  datasets the same way. That is the best ranking any prior can give; the prior
  only changes calibration.
- The joint test ranks worse than JZS everywhere.

**Calibration at 3** (figure `06_calibration.png`): with many predictors the
conditional test almost never calls a small real effect null (k = 10–20: 0–7%,
JZS 11–16%, flat partial r about 70%), but it also rarely confirms a true null
(k = 20, n = 1000: 18%, JZS 65%). With k predictors the prior expects each
coefficient to be small (on average 1/(3k) of the variance), so "no effect" and
"a tiny effect" predict nearly the same data: the Bayes factor stays near 1,
undecided rather than wrong. That is the multiplicity correction, seen from the
null side.

## Not yet done

- Whether the conditional test (05) should replace the joint one in `bayes_factor_flat_r_per_term()`: depends on study 06. A block version for factors.
- An estimation simulation (coverage, bias, RMSE with `fit_flat_r_model()`),
  and a mean-squared-error comparison to test whether the shrinkage pays off.
- Consistency when the number of predictors grows with n.
- Formal proofs of B1–B3; these are numerical checks.
