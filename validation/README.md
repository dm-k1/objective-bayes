# Validation

Three questions, one script each. Run from the project folder:

```bash
Rscript validation/01_prior_and_posterior.R   # needs cmdstanr (or rstan)
Rscript validation/02_bayes_factor.R          # no Stan needed
Rscript validation/03_simulations.R           # long
```

| Script | Question | How | Run time |
|---|---|---|---|
| `01_prior_and_posterior.R` | Does the Stan model (`R/flat_r.stan`) get the prior and the posterior right? | P1–P4: sample the prior alone; it must give a flat correlation (k = 1), R² ~ Beta(½, 1) (k = 3 correlated, k = 8), and the same per-predictor shares as the drawing recipe. S: simulation-based calibration of the posterior, k = 8, n = 60, 400 datasets. Contrast: the same prior sampled in coefficient coordinates | about 5 minutes |
| `02_bayes_factor.R` | Is the Bayes factor right? | A1: the grid against `integrate()` over a sweep of R², n, k. A2: invariance to units and recoding. A3: against a brute-force Bayes factor that draws from the prior and uses no formula. A4: against the same prior written as a mixture of g-priors. B1–B3: predictive matching, information and null consistency | about a minute |
| `03_simulations.R` | How does it behave across many datasets, compared with JZS, per-coefficient flat partial correlation and R2-D2? | 4 correlated predictors, 200 replications per cell, closed-form Bayes factors | long |
| `comparison_methods.R` | helper: the comparison methods | – | – |
| `settings.R` | helper: parallel cores, caching, `record_check()`, `save_figure()` | – | – |

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
| 01 P1–P4 | KS distance < 0.03 | pass (cmdstanr): k = 1 0.005, k = 3 0.007, k = 8 0.007, shares 0.005–0.007 |
| 01 S | rank bins, χ² p > 0.01 | pass (cmdstanr): R² p = 0.19, β₁ 0.45, σ 0.76. Contrast, coefficient coordinates: all p < 0.001; R² bins 38 15 24 23 23 18 24 19 26 **190**, σ bins **191** 15 17 … (40 expected per bin) |
| 02 A1 | log BF within 0.001 of `integrate()`, 162 settings | pass, largest error 4.1e-05 |
| 02 A2 | identical up to 1e-8 | pass |
| 02 A3 | within 0.05 of brute force (400,000 draws) | pass: 1.644 vs 1.638 (k = 3), −0.692 vs −0.691 (k = 1), 2.507 vs 2.518 (k = 5) |
| 02 A4 | within 0.002 of the Bayes factor computed as a mixture of g-priors | pass, largest difference 3.5e-04 |
| 02 B1 | BF within 0.001 of 1, k = 1–6 | pass |
| 02 B2 | increasing; log BF > 50 at R² = 0.9999 | pass (12.6, 28.6, 44.7, 60.9) |
| 02 B3 | median log BF₀₁ increasing, slope vs log n in 0.35–0.65 | pass (slope 0.53) |

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
  against "small effects called null" over n, the flat-R curve lies to the left
  of the JZS curve in both conditions. At the same rate of finding nulls, flat-R
  called the small effect null less often (e.g. at about 40% nulls found,
  uncorrelated: about 8% vs 13%, interpolated). That suggests a better
  trade-off, not just a different point on the same one, in this design; a
  precision-matched comparison (equal error rates by adjusting thresholds)
  is needed before claiming it.
- Per-coefficient flat partial correlations find nulls fastest but call small
  real effects null in 60–70% of datasets: a flat prior on each coefficient
  puts too much weight on large effects.

## Not yet done

- Precision-matched null thresholds for all four methods in 03, to test whether flat-R's better trade-off holds at equal error rates.
- An estimation simulation (coverage, bias, RMSE with `fit_flat_r_model()`),
  and a mean-squared-error comparison to test whether the shrinkage pays off.
- Consistency when the number of predictors grows with n.
- Formal proofs of B1–B3; these are numerical checks.
