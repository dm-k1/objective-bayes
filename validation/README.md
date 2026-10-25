# Validation

Three questions, one script each. Run from the project folder on your own
computer (brms with CmdStan needed for 01 and for 02's check A3):

```bash
Rscript validation/01_prior_in_brms.R
Rscript validation/02_bayes_factor.R
Rscript validation/03_simulations.R
```

| Script | Question | How | Run time |
|---|---|---|---|
| `01_prior_in_brms.R` | Does the brms model carry the prior the docs describe? | brms samples the Stan prior with no data; the draws must give a flat correlation (k = 1), R² ~ Beta(½, 1) (k = 3 correlated, k = 8), and the same per-predictor shares as the drawing recipe | a few minutes |
| `02_bayes_factor.R` | Is the Bayes factor right? | A1: the grid against `integrate()` over a sweep of R², n, k. A2: invariance to units and recoding. A3: closed form against brms bridge sampling with the Stan prior. B1–B3: predictive matching, information and null consistency | a few minutes |
| `03_simulations.R` | How does it behave across many datasets, compared with JZS, per-coefficient flat partial correlation and R2-D2? | 4 correlated predictors, 200 replications per cell, closed-form Bayes factors | long |
| `comparison_methods.R` | helper: the comparison methods | – | – |
| `settings.R` | helper: parallel cores, caching, `record_check()`, `save_figure()` | – | – |

Each check prints PASS/FAIL with its criterion; 01 and 02 stop with an error if
anything failed. Results, brms fits and figures are cached in
`validation/results/` (not tracked by git) and reused when nothing changed.

## Last results

| Check | Criterion | Last result |
|---|---|---|
| 01 P1 | k = 1, KS distance < 0.03 | pass (0.007) |
| 01 P2 | k = 3 correlated | pass, narrowly (0.025) |
| 01 P3 | k = 8 | **FAIL (0.074)**: see below |
| 01 P4 | per-predictor shares, k = 3 | pass (0.013–0.017) |
| 02 A1 | log BF within 0.001 of `integrate()`, 162 settings | pass, largest error 4.1e-05 |
| 02 A2 | identical up to 1e-8 | pass |
| 02 A3 | bridge sampling within 10% | pass in the previous version (18.76 vs 18.70; 2.833 vs 2.821); rerun |
| 02 B1 | BF within 0.001 of 1, k = 1–6 | pass |
| 02 B2 | increasing; log BF > 50 at R² = 0.9999 | pass (12.6, 28.6, 44.7, 60.9) |
| 02 B3 | median log BF₀₁ increasing, slope vs log n in 0.35–0.65 | pass (slope 0.53) |

**03 has not been rerun since the scripts changed.** The earlier run found, for
the null predictor x3 at n = 200, P(BF₀₁ > 3) = 0.61 (JZS), 0.41 (flat-R),
0.16 (R2-D2), and that the flat-R prior rarely called a small real effect
(slope 0.1) null (0.03, against 0.19 for JZS). Rerun before quoting.

## Open problem: brms under-samples near zero for larger k

Check 01 P3 fails, and the failure is in the sampler, not the formula. With 8
predictors, brms draws too few small effects: the 10% quantile of the implied
R² is 0.027 where Beta(½, 1) gives 0.010, and the 1% quantile 0.002 where it
gives 0.0001. The density has a sharp peak at all-zero coefficients for k ≥ 2
(docs/2_how_it_calculates.md, section 3), and HMC does not get into it, **without
reporting divergent transitions** (0 out of 20,000). The formula itself is
right: it integrates to the flat R for every k, matches counted draws (docs
figure 5), and k = 3 passes.

What it means: estimates for null predictors in models with many predictors
can be shrunk less than the prior implies. The Bayes factor is unaffected (it
is a closed form, not sampled).

Possible fix: estimate in coordinates without the peak, the explained share
and a direction (`archive/stan_unit_vector_model.R`), instead of the
coefficients directly. That needs Stan code outside brms' standard `b` vector.

## Not yet done

- Precision-matched null thresholds for all four methods in 03.
- An estimation simulation (coverage, bias, RMSE with `fit_flat_r_model()`),
  and a mean-squared-error comparison to test whether the shrinkage pays off.
- Consistency when the number of predictors grows with n.
- Formal proofs of B1–B3; these are numerical checks.
