# Validation

All scripts run from the project folder, on your own computer:

```bash
cd ~/github/objective-bayes
Rscript validation/01_checks.R
```

**Speed.** `settings.R` (sourced by every script) uses all cores but one:
Stan chains and simulation cells run in parallel.

**Caching.** Results are saved in `validation/results/` (not tracked by git)
and reused on the next run:

- simulations (02, 03) rerun only when their settings or design change, or
  when `rerun_simulation <- TRUE`;
- brms fits (01, check A5) are saved in `validation/results/fits/` and refitted
  only when the model, prior or data change.

**Figures.** 02 and 03 save their figures as PNG files in
`validation/results/figures/` (not tracked by git); in an interactive session
they are shown as well.

**Note:** the numbers recorded below come from runs of earlier versions of the
scripts (same methods, different random-number streams). Rerun 02 and 03 and
update them before quoting.

| Script | Question | Run time |
|---|---|---|
| `01_checks.R` | Is everything correct? Part A: calculator (ground truth by `integrate()` over a sweep of R², n and k; brute force; invariance; tutorial values; brms). Part B: default properties (k = 1, predictive matching, information and null consistency). One pass/fail summary | a few minutes (brms part cached) |
| `02_sim_one_coefficient.R` | One coefficient at a time: flat prior on the partial correlation vs JZS | long |
| `03_sim_joint_prior.R` | Several predictors: the joint possible-range prior vs JZS, per-coefficient flat r and R2-D2 | long |
| `comparison_methods.R` | Helper: the comparison methods (JZS, flat partial r, R2-D2) | – |
| `settings.R` | Helper: parallel cores, cache folders, cached simulation runner, `record_check()`, `save_figure()` | – |

## 01. Checks: is everything correct?

### Part A. Calculator

Does the code compute the formula correctly? The Bayes factor depends only on
R², n and k, so these inputs are swept directly, extremes included, instead of
simulating datasets.

| Check | Criterion | Last result |
|---|---|---|
| A1a. Ground truth, ordinary evidence (\|log BF\| < 50) | within 0.001 on the log scale | pass (largest error 4e-05) |
| A1b. Ground truth, all settings | within 0.001 | pass (largest error 4e-05) after packing the direction grid; the even grid was off by up to 28 at n = 10,000, R² = 0.99 |
| A2. Brute force (k = 1, 3, 5) | within 0.05 on the log scale | pass (e.g. 1.644 vs 1.652) |
| A3. Invariance (rescaling predictors or outcome, recoding) | identical up to 1e-8 | pass |
| A4. Reproduces the tutorial | within 0.1% of 18.696, 0.3827, 1.5872 | pass |
| A5. brms bridge sampling vs closed form (continuous x1; factor as a block) | within 10% | pass: 18.76 vs 18.70 (continuous), 2.833 vs 2.821 (factor) |

The ground truth computes the same Bayes factor with R's `integrate()`:
adaptive quadrature that adds points where the curve changes quickly and
reports its own error estimate.

### Part B. Default properties

| Check | Criterion | Last result |
|---|---|---|
| B1. k = 1 vs per-coefficient flat r | within 3% (flat r plugs in σ̂) | pass |
| B2a. Predictive matching (n = k + 1) | BF within 0.001 of 1, k = 1–6 | pass |
| B2b. Information consistency (R² → 1) | increasing; log BF > 50 at R² = 0.9999 | pass (12.6, 28.6, 44.9, 59.9) |
| B2c. Null consistency, simulated | median log BF₀₁ increasing; slope vs log n in 0.35–0.65 | pass (slope 0.53) |

## 02. One coefficient: flat partial r vs JZS

Two predictors (x and covariate z). True effect 0–0.5, n = 20–500,
correlation between x and z 0 or 0.5, 200 replications per cell.

Main findings:

- **Testing:** at n ≥ 100 the two are equivalent once the flat-r null threshold
  is raised to about BF₀₁ > 4–4.5 (instead of 3); this is the precision-matching
  analysis in section D.
- **Small samples:** the flat-r prior gives much more evidence for the null,
  including for small real effects. It finds nulls faster but calls small
  effects null more often. This follows from its prior putting less weight near
  zero than JZS (tutorial Step 1).
- **Estimation:** about half the shrinkage bias of JZS. RMSE is worse when
  effects are null or small at small n.
- The partial-correlation interval (`r_coverage`) ignores uncertainty in the
  conversion factor and is unreliable; use a Fisher z interval for ρ.

## 03. Joint prior vs alternatives

Four predictors (equicorrelation 0 or 0.5). x1's effect varies, x2 = 0.3, x3
and x4 are null. Tested: x1 and x3. 200 replications per cell.

| Method | Finds true nulls (x3, n = 200, BF₀₁ > 3) | Small effect called null (x1 = 0.1, n = 200) |
|---|---|---|
| per-coefficient flat r | fastest | most often |
| JZS default | 0.61 | 0.19 |
| joint possible-range | 0.41 | 0.03 |
| R2-D2 (same R² prior) | 0.16 | almost never |

- For real effects the joint prior is practically identical to JZS.
- The four methods sit on a trade-off between finding nulls quickly and not
  calling small real effects null. Ordered by how readily they support the
  null: per-coefficient flat r > JZS > joint > R2-D2.
- R2-D2 has slightly more power with correlated predictors.

## Not yet done

- Precision-matched thresholds for all four methods in study 03 (as done for
  two methods in 02, section D).
- An estimation simulation for the joint prior (coverage, bias, RMSE with
  `fit_possible_range_model()`).
- Consistency when the number of predictors grows with n (checked
  numerically outside this repository: the joint prior stays consistent, while
  Zellner–Siow with scale 1 fails at k/n = 0.25). Should be added as a script.
- Bridge-sampling agreement (check 6) for factors, interactions and
  polynomials, not only continuous predictors.
- Formal proofs of the criteria in check 4; these are numerical checks.
