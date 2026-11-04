# Validation run

- **started:** 2026-11-04 04:23 CET
- **commit:** e221561
- **R:** R version 4.6.1 (2026-06-24)
- **CmdStan:** 2.40.0
- **machine:** Darwin arm64
- **cores:** 7
- **studies:** 00 01 02 03 04 05 06 07 08 09 10 11

| Study | Status | Minutes | Checks passed | Checks failed |
|---|---|---|---|---|
| 00_initial_simulation | pass | 1.2 | 3 | 0 |
| 01_prior_and_posterior | pass | 7.4 | 19 | 0 |
| 02_bayes_factor | pass | 1.3 | 45 | 0 |
| 03_simulations | pass | 11.5 | 2 | 0 |
| 04_matched_comparison | pass | 2.5 | 0 | 0 |
| 05_conditional_test | pass | 0.4 | 13 | 0 |
| 06_broad_simulation | pass | 21.3 | 0 | 0 |
| 07_estimation | pass | 83.7 | 0 | 0 |
| 08_growing_k | pass | 1.4 | 12 | 0 |
| 09_variable_selection | pass | 3.0 | 5 | 0 |
| 10_calibration | FAIL | 68.3 | 29 | 1 |
| 11_factor_terms | pass | 63.6 | 3 | 0 |

## 00: the first result

200,000 valid correlation matrices of (y, x1 ... xk) per k. Keep each partial r flat
and R-squared grows with k; keep R or R-squared flat and each partial r shrinks with k.

| k | P(\|partial r\| > 0.3) | P(R² > 0.5) | P(\|partial r\| > 0.3) | P(R² > 0.5) | P(\|partial r\| > 0.3) | P(R² > 0.5) |
| | partial r flat | partial r flat | R flat | R flat | R-squared flat | R-squared flat |
|---|---|---|---|---|---|---|
| 1 | 0.70 | 0.29 | 0.70 | 0.29 | 0.91 | 0.50 |
| 2 | 0.70 | 0.50 | 0.50 | 0.29 | 0.70 | 0.50 |
| 3 | 0.70 | 0.65 | 0.41 | 0.29 | 0.60 | 0.50 |
| 4 | 0.70 | 0.75 | 0.36 | 0.29 | 0.54 | 0.50 |
| 5 | 0.70 | 0.82 | 0.32 | 0.29 | 0.49 | 0.50 |
| 6 | 0.70 | 0.87 | 0.29 | 0.29 | 0.46 | 0.50 |
| 7 | 0.70 | 0.91 | 0.27 | 0.29 | 0.43 | 0.50 |
| 8 | 0.70 | 0.94 | 0.25 | 0.29 | 0.40 | 0.50 |
| 9 | 0.70 | 0.96 | 0.23 | 0.29 | 0.38 | 0.50 |
| 10 | 0.70 | 0.97 | 0.22 | 0.29 | 0.36 | 0.50 |

Studies 03, 04, 06, 07 and 09 are comparisons (results in their logs and
figures); they have few or no PASS/FAIL checks. Full output: `cache/logs/`.

## Failed

- 10_calibration: [FAIL] P4. flat R-squared: conditional within a factor 2 of Savage-Dickey (estimable) largest |difference| 0.874 on the log scale (factor 2.40); mean |difference| 0.283
