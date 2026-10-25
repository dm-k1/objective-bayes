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
| `02_bayes_factor.R` | Is the Bayes factor right? | A1: the grid against `integrate()` over a sweep of R², n, k. A2: invariance to units and recoding. A3: against a brute-force Bayes factor that draws from the prior and uses no formula. B1–B3: predictive matching, information and null consistency | about a minute |
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
| 01 P1–P4 | KS distance < 0.03 | pass: k = 1 0.006, k = 3 0.006, k = 8 0.007, shares 0.006–0.012 |
| 01 S | rank bins, χ² p > 0.01 | pass: R² p = 0.29, β₁ 0.50, σ 0.24. Contrast, coefficient coordinates: R² p = 0.03, lowest bin 61 where 40 is expected |
| 02 A1 | log BF within 0.001 of `integrate()`, 162 settings | pass, largest error 4.1e-05 |
| 02 A2 | identical up to 1e-8 | pass |
| 02 A3 | within 0.05 of brute force (400,000 draws) | pass: 1.644 vs 1.638 (k = 3), −0.692 vs −0.691 (k = 1), 2.507 vs 2.518 (k = 5) |
| 02 B1 | BF within 0.001 of 1, k = 1–6 | pass |
| 02 B2 | increasing; log BF > 50 at R² = 0.9999 | pass (12.6, 28.6, 44.7, 60.9) |
| 02 B3 | median log BF₀₁ increasing, slope vs log n in 0.35–0.65 | pass (slope 0.53) |

**Why the Stan model samples R and a direction.** The prior as a density on the
coefficients has a sharp peak at zero (`docs/2_how_it_calculates.md`, section
3). A sampler moving through the coefficients, as brms does with a custom prior,
misses it for larger k, without divergent transitions: with brms at k = 8 the
10% quantile of the implied R² was 0.027 where the prior gives 0.010 (KS
distance 0.074; `archive/check_brms_prior.R`). Check 01 shows the contrast.

**03 has not been rerun since the scripts changed.** The earlier run found, for
the null predictor x3 at n = 200, P(BF₀₁ > 3) = 0.61 (JZS), 0.41 (flat-R),
0.16 (R2-D2), and that the flat-R prior rarely called a small real effect
(slope 0.1) null (0.03, against 0.19 for JZS). Rerun before quoting.

## Not yet done

- Precision-matched null thresholds for all four methods in 03.
- An estimation simulation (coverage, bias, RMSE with `fit_flat_r_model()`),
  and a mean-squared-error comparison to test whether the shrinkage pays off.
- Consistency when the number of predictors grows with n.
- Formal proofs of B1–B3; these are numerical checks.
