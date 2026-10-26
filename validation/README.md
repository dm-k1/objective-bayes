# Validation

One script per question. Run from the project folder:

```bash
Rscript validation/01_prior_and_posterior.R   # needs cmdstanr (or rstan)
Rscript validation/02_bayes_factor.R          # no Stan needed
Rscript validation/03_simulations.R           # a few minutes
Rscript validation/04_matched_comparison.R    # needs the results of 03
Rscript validation/05_conditional_test.R
Rscript validation/06_broad_simulation.R
Rscript validation/07_estimation.R            # hours; resumes if interrupted
Rscript validation/08_growing_k.R
Rscript validation/09_variable_selection.R
Rscript validation/10_calibration.R           # needs cmdstanr (or rstan)
```

The scripts load the package from this folder with `pkgload::load_all()`.

**All at once:** `Rscript validation/run_all.R` runs every study, each in its
own R process, and writes `results/run_summary.md` (the commit, R and CmdStan
versions, cores, and every PASS/FAIL check) plus a log per study in
`results/logs/`. `run_all.R quick` runs the correctness checks 01, 02 and 05;
`run_all.R 03 04 09` runs a selection. `VALIDATION_CORES=32` sets the cores.

### Running it on another machine

The environment is pinned: `renv.lock` (the R packages) and `Dockerfile` (R,
those packages and CmdStan). renv is not activated in this folder; the
lockfile is the record, and the image installs from it.

| Where | How |
|---|---|
| Any machine with Docker (a cloud VM, a SURF Research Cloud workspace) | `git clone` the repository, then `docker build -t flatcor-validation .` (once, 20–40 minutes) and `docker run --rm -v "$PWD":/project -e VALIDATION_CORES=32 flatcor-validation`. Results land in `validation/results/` of the clone |
| A machine with R but no Docker | `Rscript -e 'install.packages("renv"); renv::restore(library = .libPaths()[1], prompt = FALSE)'` (into your normal library), then `Rscript -e 'cmdstanr::install_cmdstan(version = "2.40.0")'`, then `Rscript validation/run_all.R` |
| A cluster without Docker (Snellius) | the same as the previous row, or the image through Apptainer: `apptainer build flatcor.sif docker-daemon://flatcor-validation:latest` |

To update the lockfile after installing or upgrading packages:
`Rscript -e 'renv::snapshot(type = "implicit", exclude = c("flatcor", "rstan"), prompt = FALSE)'`
(rstan is only the fallback when cmdstanr is missing, and is left out to keep
the image small).

| Script | Question | How | Run time |
|---|---|---|---|
| `01_prior_and_posterior.R` | Does the Stan model (`inst/stan/flatcor.stan`) get the prior and the posterior right, for both priors? | P1–P4 (flat-R): sample the prior alone; it must give a flat correlation (k = 1), R² ~ Beta(½, 1) (k = 3 correlated, k = 8), and the same per-predictor shares as the drawing recipe. P5–P7 (uniform correlations): R² ~ Beta(k/2, 1), and every partial correlation flat. S: simulation-based calibration of the posterior, k = 8, n = 60, 400 datasets, both priors. Contrast: flat-R sampled in coefficient coordinates | about 10 minutes |
| `02_bayes_factor.R` | Is the Bayes factor right, for both priors? | Flat-R: A1 the grid against `integrate()` over a sweep of R², n, k; A2 invariance to units and recoding; A3 against a brute-force Bayes factor that draws from the prior and uses no formula; A4 against the same prior written as a mixture of g-priors; B1–B4 predictive matching, information and model-selection consistency (B3b: the rate with and without an effect in the smaller model). Uniform correlations: the same as CA1–CA3 and CB1–CB4 | a few minutes |
| `03_simulations.R` | Testing one coefficient across many datasets | 4 correlated predictors, x1's slope varies, x3 null; all one-coefficient methods (below); checks U1–U2 | 5–10 minutes |
| `04_matched_comparison.R` | Is a difference in 03 discrimination or calibration? | AUC (threshold-free), each method minus JZS with a paired bootstrap, matched thresholds, and why: how much the evidence depends on what the other predictors explain | about 3 minutes |
| `05_conditional_test.R` | Are the conditional tests right, and why η = 1? | C1–C4: equals the joint test at k = 1; the implied priors for one coefficient (flat-R: depends only on k; uniform correlations: a flat partial correlation for any k). C5–C6: both grids against `integrate()`, n up to 10⁶. C7: under LKJ(η) one partial correlation is Beta(η, η) for every k, so flat only at η = 1. C8: against the exact flat-prior partial-correlation Bayes factor of Kucharský et al. (2026). C9–C11: the block test for a term with q columns (implied prior by simulation; against `integrate()`; the likelihood ratio from the data) | a few minutes |
| `06_broad_simulation.R` | Does the picture of 03–04 hold beyond one design? | k = 2–20, n = 50–1000, correlation 0–0.7, sparse and dense other coefficients; all one-coefficient methods | 20–30 minutes |
| `07_estimation.R` | How good are the estimates? | bias, coefficient error and prediction error relative to `lm()`, coverage of 95% intervals; flat-R, uniform correlations (`stan_lm`, η = 1), JZS, Zellner–Siow, `lm()`; needs rstanarm; resumes if interrupted | 2.5–4 hours |
| `08_growing_k.R` | Consistency when k grows with n (k = n/4, k = √n) | null true and small model true; all whole-model methods | under a minute |
| `09_variable_selection.R` | Choosing among 4,096 subsets of 12 candidates (4 real) | all whole-model methods and BIC; uniform and beta-binomial model priors; true model found, false inclusions, missed; whole trade-off curves and AUC; checks V1–V5 | 3–5 minutes |
| `10_calibration.R` | Is the conditional test calibrated, and how far is it from one prior for everything? | P1–P2: 160,000 datasets from the joint uniform-correlations prior (half with x1 = 0) and from the JZS prior; share of nulls per bin of log BF₁₀ against 1/(1 + BF₁₀). P3–P5: 20 datasets, the conditional Bayes factor against the Savage–Dickey ratio for x1's partial correlation under the joint model (`fit_regression()` draws, logspline) | about 10 minutes |
| `comparison_methods.R` | helper: the common method set (names, colours, Bayes factors) | – | – |
| `settings.R` | helper: parallel cores, caching, `record_check()`, `save_figure()` | – | – |

## The methods compared (03 on)

One set, defined once in `comparison_methods.R` and used by every study.

| Prior | What it says | Whole models | One coefficient |
|---|---|---|---|
| **flat-R** | R² ~ Beta(½, 1) for every k | ✓ | joint and conditional |
| **uniform correlations** | R² ~ Beta(k/2, 1): every correlation matrix equally likely; `stan_lm`'s prior with η = 1 | ✓ | joint and conditional |
| JZS default | Cauchy, scale 0.354 (BayesFactor, JASP) | ✓ | joint |
| Zellner–Siow | Cauchy, scale 1; R² ~ Beta(k/2, ½) | ✓ | joint |

- **Joint:** the model with against without the coefficient.
- **Conditional:** only the coefficient's partial correlation, with the prior the
  joint prior implies for one coefficient. For uniform correlations that is a
  flat partial correlation whatever k is, so this is the "flat partial r per
  coefficient" of earlier versions, now computed exactly.
- Estimation (07): flat-R and uniform correlations (both `fit_regression()`,
  the same Stan model), JZS, Zellner–Siow and `lm()`; the conditional tests have
  no estimation version.
- Dropped: R2-D2 and the approximate flat partial correlation from 03
  (`archive/03_simulations_with_r2d2.R`), and the experiment R² ~ Beta(k/(k+1), 1)
  from 09 (its result is kept below).

Package tests (`devtools::test()`): 56 pass on cmdstanr.

Each check prints PASS/FAIL with its criterion; 01 and 02 stop with an error if
anything failed. Simulation results and figures are cached in
`validation/results/` (not tracked by git).

**Simulation-based calibration** (01, S): draw the true parameters from the
prior, simulate a dataset, fit it, and record where the true value falls among
the posterior draws. If the sampler is right, every position is equally likely,
so the counts per bin should be flat.

## Last results

What each result means for the project's claims, with its status, is summarised
in [`docs/0_overview.md`](../docs/0_overview.md). Every study below was run on
26 October 2026 on the current code (commit 604db78 plus the portability
changes), with the common method set, on an 8-core Mac: 03–06, 08 and 09 in
about 20 minutes, 07 in 56 minutes. The earlier results are kept outside git in
`validation/results_before_2026-10-26/`.

| Check | Criterion | Last result |
|---|---|---|
| 01 P1–P4 | KS distance < 0.03 | pass (cmdstanr, current model): k = 1 0.006, k = 3 0.006, k = 8 0.008, shares 0.004–0.009 |
| 01 S | rank bins, χ² p > 0.01 | pass (cmdstanr, current model). Flat-R: R² p = 0.61, β₁ 0.86, σ 0.67. Uniform correlations: R² 0.43, β₁ 0.94, σ 0.76. Contrast, flat-R in coefficient coordinates: all p < 0.001; R² bins 36 18 24 24 26 18 19 19 25 **191**, σ bins **189** 14 25 … (40 expected per bin) |
| 02 A1 | log BF within 0.001 of `integrate()`, 162 settings | pass, largest error 4.1e-05 |
| 02 A2 | identical up to 1e-8 | pass |
| 02 A3 | within 0.05 of brute force (400,000 draws) | pass: 1.644 vs 1.629 (k = 3), −0.692 vs −0.689 (k = 1), 2.507 vs 2.524 (k = 5) |
| 02 A4 | within 0.002 of the Bayes factor computed as a mixture of g-priors | pass, largest difference 3.5e-04 |
| 02 B1 | BF within 0.001 of 1, k = 1–6 | pass |
| 02 B2 | increasing; log BF > 50 at R² = 0.9999 | pass (12.6, 28.6, 44.7, 60.9) |
| 02 B2b | bounded at n = k + 2, growing at n = k + 3 (the analytic boundary) | pass |
| 02 B3 | median log BF₀₁ increasing, slope vs log n in 0.35–0.65 | pass (slope 0.53) |
| 02 B3b | q = 1, 3, 6 null predictors added: slope within 20% of q/2 when the smaller model is x1 (y = 0.5 x1 + e); in 0.35–0.65 for every q when it is the intercept-only model | pass (0.55, 1.55, 2.89; 0.51, 0.49, 0.49) |
| 02 B4 | effect x1 = 0.15: median log BF₁₀ increasing, > 20 at n = 5000 | pass (0.9, 9.3, 54.1 at n = 200, 1000, 5000) |
| 01 P5–P7 | uniform correlations, KS distance < 0.03 | pass: R² k = 3 0.003, k = 8 0.004; partial correlations flat 0.012, 0.012, 0.006 |
| 02 CA1 | uniform correlations, log BF within 0.001 of `integrate()`, 162 settings | pass, largest error 2.6e-07 |
| 02 CA2 | identical up to 1e-8 | pass |
| 02 CA3 | within 0.05 of brute force (400,000 flat-R draws weighted by k R^(k−1): plain draws from this prior rarely reach small R and come out too low, 0.161 for k = 5) | pass: 0.817 vs 0.815 (k = 3), −0.692 vs −0.692 (k = 1), 0.284 vs 0.319 (k = 5) |
| 02 CB1–CB2 | as B1, B2b | pass; information consistency from n = k + 3, as flat-R |
| 02 CB3 | slope within 20% of q/2 for both smaller models | pass (x1: 0.50, 1.47, 2.89; intercept only: 0.50, 1.46, 2.87) |
| 02 CB4 | as B4 | pass (−0.0, 8.5, 53.3) |
| 05 C1–C7 | see the script | pass; C5 largest error 1.7e-04, C6 3.5e-06; C7 KS to Beta(η, η) ≤ 0.018, to flat at η = 2 ≥ 0.099 |
| 05 C8 | Kucharský et al. (2026, eq. 2, flat prior): their example within 0.01 of 45.88; ours with the predictor's residual SD unknown within 1e-5 of theirs (48 cases) | pass: 45.877; largest 4.7e-07. Ours (fixed X) − theirs: largest −0.262 (n = 20, k = 8, r = 0), ≤ 0.005 at n = 1000 |
| 10 P1 | data from the joint uniform-correlations prior, half with x1 = 0: share of nulls per bin of log BF₁₀ equals 1/(1 + BF₁₀), χ² over the bins p > 0.01 and logistic slope within 3 SE of 1 | pass in all four cells (k = 3, 8; n = 50, 200; 20,000 datasets each): χ² p = 0.40–0.79, slopes 0.98–1.03 (SE 0.02), intercepts within 2 SE of 0 |
| 10 P2 | the same with data from the JZS prior: miscalibration shows (χ² p < 0.001) | pass: χ² 525–593 on 10–11 bins; intercept −0.24 to −0.49 (fewer nulls than 1/(1 + BF₁₀) says) |
| 10 P3 | Savage–Dickey estimate: logspline and kernel within 0.25, where the chain SE is below 0.1 (16 of 20 datasets) | pass, largest difference 0.10 |
| 10 P4 | conditional log BF₁₀ within log 2 of Savage–Dickey under the joint model (the 16) | pass: mean difference −0.005, mean absolute 0.028, largest 0.092 (factor 1.10) |
| 10 P5 | the other 4 (log BF₁₀ 2–9, too few draws near 0): within 3 Monte Carlo SE | pass: 4.98 vs 4.39 (SE 0.27), 2.98 vs 2.82 (0.15), 8.97 vs 8.33 (0.70), 1.96 vs 2.17 (0.12) |

The criteria of P3–P5 were set after a first run: P3 originally covered all 20
datasets with a tolerance of 0.1 and failed on the 4 with strong evidence, where
the density at 0 cannot be estimated from draws. P3 and P4 were then restricted
to datasets with a chain SE below 0.1 (tolerance 0.25 for P3), and P5 added for
the other 4.
| 05 C9–C11 | see the script | pass; C9 KS ≤ 0.020 (LKJ), ≤ 0.005 (draws); C10 largest error 1.2e-04; C11 1.6e-14 |

**Why the Stan model samples R and a direction.** The prior as a density on the
coefficients has a sharp peak at zero (`docs/2_how_it_calculates.md`, section
3). A sampler moving through the coefficients, as brms does with a custom prior,
misses it for larger k, without divergent transitions: with brms at k = 8 the
10% quantile of the implied R² was 0.027 where the prior gives 0.010 (KS
distance 0.074; `archive/check_brms_prior.R`). Check 01 shows the contrast.

## In short: what the studies say about the two priors

| Question | Flat-R | Uniform correlations | Study |
|---|---|---|---|
| Separating nulls from small effects, one coefficient (AUC) | conditional test = uniform's, ≥ JZS | the same | 04, 06 |
| Confirming a true null, one coefficient (BF₀₁ > 3) | falls with the number of covariates: 0.55 (k = 2) to 0.06 (k = 20) | about 0.87 for every k | 06 |
| Calling a small effect (0.1) null | rare: 0.19 (k = 2) to 0.005 (k = 20) | about half: 0.52–0.60 | 06 |
| Evidence for the null when k grows with n (k = n/4, n = 2000) | log BF₀₁ 1.9 | 408 | 08 |
| Variable selection | the most inclusive | next to Zellner–Siow, the strict end | 09 |
| Estimation: prediction error against OLS | the largest gains at small n | moderate gains | 07 |
| Estimation: 95% interval coverage, real effects | as low as 0.76 in sparse designs | 0.92 or higher in every cell | 07 |

Both conditional tests rank datasets by the same partial correlation, so they
separate nulls from effects equally well; the priors differ only in where the
Bayes factors land relative to 3 (calibration).

## 03: one coefficient, one design (200 datasets per cell)

Four predictors, x1's slope varies, x3 is null; n = 20 to 500; predictors
uncorrelated or correlated 0.5. The null predictor found, P(BF₀₁ > 3):

| Method | n = 200, uncorrelated | n = 200, correlated | n = 500, uncorrelated | n = 500, correlated |
|---|---|---|---|---|
| JZS default | 0.55 | 0.49 | 0.81 | 0.76 |
| Zellner–Siow | 0.85 | 0.81 | 0.96 | 0.94 |
| flat-R joint | 0.18 | 0.17 | 0.67 | 0.62 |
| **flat-R conditional** | 0.41 | 0.39 | 0.77 | 0.72 |
| uniform correlations joint | 0.82 | 0.79 | 0.94 | 0.92 |
| **uniform correlations conditional** | 0.90 | 0.87 | 0.97 | 0.96 |

A small real effect (x1 = 0.1) called null, P(BF₀₁ > 3):

| Method | n = 200, uncorrelated | n = 200, correlated | n = 500, uncorrelated | n = 500, correlated |
|---|---|---|---|---|
| JZS default | 0.22 | 0.26 | 0.13 | 0.27 |
| Zellner–Siow | 0.52 | 0.55 | 0.28 | 0.45 |
| flat-R joint | 0.08 | 0.15 | 0.08 | 0.22 |
| **flat-R conditional** | 0.16 | 0.20 | 0.12 | 0.22 |
| uniform correlations joint | 0.47 | 0.51 | 0.25 | 0.40 |
| **uniform correlations conditional** | 0.62 | 0.64 | 0.34 | 0.52 |

- An effect of 0.3 is found (BF₁₀ > 3) in 0.69–0.80 of datasets at n = 100 and
  in 0.975 or more from n = 200 on, by every method; false effects for a true
  null stay at 5.5% or below.
- Checks U1 (at k = 1 all four flat-R and uniform-correlations versions agree,
  4.9e-06) and U2 (the uniform-correlations grid against Monte Carlo, 0.0015) pass.

## 04: discrimination, threshold-free

AUC (how often BF₀₁ is larger for a null dataset than for a small effect of
0.1), each method minus JZS, paired bootstrap 95% CI:

| | n = 100 | n = 200 | n = 500 |
|---|---|---|---|
| both conditional tests, uncorrelated | +0.003 [−0.018, 0.023] | +0.003 [−0.010, 0.015] | −0.000 [−0.004, 0.003] |
| both conditional tests, correlated (0.5) | **+0.047** [0.027, 0.067] | **+0.027** [0.014, 0.040] | **+0.028** [0.018, 0.039] |
| uniform correlations joint, correlated | **+0.029** [0.015, 0.042] | **+0.018** [0.009, 0.028] | **+0.022** [0.014, 0.031] |
| flat-R joint, correlated | **−0.039** [−0.058, −0.020] | **−0.034** [−0.049, −0.019] | **−0.028** [−0.036, −0.021] |

- The two conditional tests have identical AUCs: both rank datasets by the
  coefficient's partial correlation, the best ranking any prior can give.
- **Why the joint tests differ** (M4, correlated, n = 200): at a fixed t-statistic
  for x1, the log BF for x1 changes by −5.19 per unit of R² explained by the
  other predictors under flat-R joint, −1.93 under JZS, −0.52 under Zellner–Siow
  and uniform correlations joint, and about 0 under both conditional tests. Rank
  correlation of the log BF with t²: 0.78 (flat-R joint), 0.93 (JZS), 0.99
  (Zellner–Siow, uniform correlations joint), 1.00 (conditional tests). The
  joint flat-R test is the one most affected by the conditional Lindley paradox.

## 05: the conditional tests are computed correctly

Checks C1–C7 (the table under "Last results"): both conditional tests against
`integrate()` up to n = 10⁶, their implied priors, and C7, why η = 1.

**C8: the uniform-correlations conditional test is not the Kucharský et al.
(2026) test, but differs from it in one assumption only.** Their exact Bayes
factor for a partial correlation with a flat prior (eq. 2, α = 1) is

$$\mathrm{BF}_{10} = \frac{B(\tfrac12, 1 + h)}{B(\tfrac12, 1)}\,
{}_2F_1\!\left(h, h;\ \tfrac32 + h;\ r^2\right),\qquad h = \frac{n - k}{2},$$

with $k$ our number of predictors (their $k - 1$ controlling variables). It
reproduces their worked example (BF₁₀ = 45.88). Ours is lower in 47 of 48 cases
(the 48th: +1e-5, r = 0.5, n = 1000, k = 1) on the grid r = 0, 0.1, 0.3, 0.5; n = 20, 50, 200, 1000; k = 1, 3, 8 (log BF):

| | n = 20 | 50 | 200 | 1000 |
|---|---|---|---|---|
| k = 1, r = 0 | −0.054 | −0.023 | −0.006 | −0.001 |
| k = 3, r = 0 | −0.105 | −0.044 | −0.011 | −0.002 |
| k = 8, r = 0 | −0.262 | −0.098 | −0.024 | −0.005 |
| k = 8, r = 0.5 | −0.128 | −0.029 | −0.005 | −0.001 |

The difference is largest at r = 0, shrinks as 1/n and grows with k. Its cause is pinned
down exactly: their model is multivariate normal, so the predictors are random
and the tested predictor's residual SD (given the others) is an unknown
parameter, with a 1/σ prior. Ours conditions on X, so the effect size f is in
units of the observed residual SD. Put their uncertainty into ours (f in sample
units = f in population units × √(W/n), W ~ χ²(n − k)) and the two agree to
5e-7 in all 48 cases. The priors on the noise and the other coefficients play
no part: they are the same in both. One visible consequence: at r = 0 our
Bayes factor does not depend on k, theirs moves toward 1 as k grows. How they
compare with the other methods is in 03, 04 and 06.

## 06: across designs (k = 2–20, n = 50–1000, correlation 0–0.7, 216 cells)

**Discrimination** (AUC null vs 0.1, minus JZS, mean over n and the pattern of
the other coefficients):

| k | 2 | 5 | 10 | 20 |
|---|---|---|---|---|
| both conditional tests, correlation 0 | +0.001 | +0.004 | +0.002 | +0.007 |
| both conditional tests, correlation 0.7 | **+0.035** | **+0.045** | **+0.041** | **+0.031** |
| uniform correlations joint, correlation 0.7 | +0.022 | +0.027 | +0.020 | +0.012 |
| flat-R joint, correlation 0.7 | −0.064 | −0.035 | −0.022 | −0.013 |

**Calibration at 3**, averaged over everything except the method and k:

| | k = 2 | k = 5 | k = 10 | k = 20 |
|---|---|---|---|---|
| True null found: flat-R conditional | 0.55 | 0.29 | 0.20 | **0.06** |
| True null found: uniform correlations conditional | 0.86 | 0.88 | 0.87 | 0.88 |
| True null found: JZS | 0.54 | 0.38 | 0.30 | 0.30 |
| Small effect (0.1) called null: flat-R conditional | 0.19 | 0.05 | 0.02 | 0.005 |
| Small effect (0.1) called null: uniform correlations conditional | 0.52 | 0.56 | 0.57 | **0.60** |
| Small effect (0.1) called null: JZS | 0.21 | 0.12 | 0.08 | 0.11 |
| Effect 0.2 found: flat-R / uniform / JZS | .53 / .51 / .56 | .52 / .47 / .55 | .50 / .45 / .53 | .48 / .44 / .52 |

False effects for a true null: at most 3.3% for every method.

- **Uniform correlations behaves the same whatever the number of covariates**:
  its prior for one coefficient is a flat partial correlation for every k. It
  confirms true nulls, and it calls an effect of 0.1 (a partial correlation of
  about 0.1, below Cohen's "small") null in about half the datasets at these n.
- **Flat-R's conditional test loses the ability to confirm a null as covariates
  are added**: its prior expects each coefficient to be small (on average
  1/(3k) of the variance), so "no effect" and "a tiny effect" predict nearly the
  same data.

## 08: k growing with n

Median log BF₀₁ for the true null model against the full model (200 datasets
per n):

| | n = 200 | n = 2000 |
|---|---|---|
| k = n/4: flat-R / uniform correlations / JZS / Zellner–Siow | 1.3 / 39.2 / 15.1 / 41.2 | **1.9** / 407.7 / 149.3 / 410.9 |
| k = √n: flat-R / uniform correlations / JZS / Zellner–Siow | 1.7 / 18.1 / 8.8 / 19.6 | **2.5** / 82.6 / 44.2 / 84.7 |

All four keep growing (checks G1–G2 pass for flat-R), but flat-R's evidence for
the null barely moves: with 500 null predictors at n = 2000, BF₀₁ ≈ 7. When the
small model is true, every method finds it (flat-R: log BF 175 at n = 2000,
k = n/4; uniform correlations 444).

## 07: estimation (24 designs, 100 datasets each; 0 divergent transitions in 4,800 Stan fits)

Excess prediction error for a new observation, relative to OLS (below 1 =
better), mean over correlation and pattern:

| | k = 3: n = 30 / 60 / 200 | k = 8: n = 30 / 60 / 200 |
|---|---|---|
| flat-R | 0.70 / 0.85 / 1.02 | **0.51** / 0.73 / 0.94 |
| uniform correlations | 0.83 / 0.89 / 0.97 | 0.74 / 0.84 / 0.95 |
| JZS default | 0.68 / 0.80 / 0.95 | 0.56 / 0.72 / 0.91 |
| Zellner–Siow | 0.86 / 0.91 / 0.97 | 0.77 / 0.85 / 0.95 |

Coverage of the 95% intervals (OLS: 0.945–0.953):

| | real effects, k = 3 / 8 | true zeros, k = 3 / 8 | lowest cell, real effects |
|---|---|---|---|
| flat-R | 0.947 / 0.928 | 0.977 / 0.982 | **0.76** (k = 8, n = 60, sparse) |
| uniform correlations | 0.957 / 0.940 | 0.952 / 0.948 | 0.92 |
| JZS default | 0.959 / 0.951 | 0.965 / 0.971 | 0.87 |
| Zellner–Siow | 0.952 / 0.939 | 0.953 / 0.946 | 0.91 |

Mean bias of real effects: flat-R −0.054, JZS −0.037, uniform correlations
−0.016, Zellner–Siow −0.014, OLS 0.000.

- Flat-R shrinks the most: the best predictions at small n with many predictors,
  but intervals that miss real effects in sparse designs (0.76).
- Uniform correlations shrinks less, with intervals close to nominal in every
  cell. It is the same Stan model as flat-R with one line changed
  (`fit_regression(prior = "uniform_correlations")`).

## 09: variable selection (12 candidates, 4 real: 0.3, 0.2, 0.15, 0.1; 100 datasets per cell)

Finding the exact true model is hard for every method; flat-R and JZS tie for
best (n = 500: 0.37 vs 0.38 uncorrelated, 0.15 vs 0.15 correlated, uniform
model prior). The methods differ in where they sit on a trade-off
(n = 500, uncorrelated, beta-binomial model prior):

| | real predictors missed (of 4) | null predictors included (of 8) | inclusion of the 0.1 effect |
|---|---|---|---|
| flat-R | **0.45** | **0.91** | **0.64** |
| JZS default | 0.58 | 0.42 | 0.58 |
| Zellner–Siow (1) | 0.89 | 0.12 | 0.43 |
| BIC | 1.01 | 0.05 | 0.34 |

- Flat-R is the most inclusive: it misses the fewest real predictors and admits
  the most null ones. At n = 50 its inclusion probabilities for null predictors
  are near 0.5: close to undecided.
- Why: the other priors claim more explained variance for bigger models, so
  bigger models that don't show it are penalised automatically. Flat-R claims
  the same R² for every model size, so the extra penalty for leaving the null
  model is gone: against the intercept-only model, evidence for the null grows
  like ½ log n whatever the size of the alternative. Between models that both
  contain an effect, the usual penalty of ½ log n per extra predictor remains
  (`docs/3_consistency.md`; check 02 B3b). Control over how many predictors
  enter then has to come from the prior over models.
- With correlated predictors, flat-R and JZS are nearly identical.
- No method dominates; flat-R sits at the sensitive end, BIC at the strict end.

### Experiment: growing from flat-R toward flat-R² (negative result; earlier run, method since dropped)

Idea: let the prior on R² open up with model size, R² ~ Beta(a, 1) with
a = k / (k + 1): flat on R for one predictor (a = ½), flat on R² as k grows
(a → 1). Same uniform direction; as a g-mixture, g/n = B·W with
B ~ Beta(a, k/2 − a), W ~ Inverse-Gamma(1, ½). Checks: V2 (a = ½ reproduces
flat-R, largest difference 7e-5) and V3 (implied R² is Beta(k/(k+1), 1), KS ≤ 0.004) pass.

Result: it barely differs from flat-R. Only with uncorrelated predictors and
the beta-binomial model prior does it move, about a quarter of the way toward
JZS, and it pays with more missed predictors:

| n | false inclusions: flat-R / growing / JZS | missed: flat-R / growing / JZS |
|---|---|---|
| 50 | 2.33 / 1.81 / 1.00 | 1.98 / 2.22 / 2.80 |
| 100 | 2.38 / 1.78 / 0.83 | 1.45 / 1.69 / 2.22 |
| 200 | 2.02 / 1.79 / 0.79 | 0.98 / 1.11 / 1.58 |
| 500 | 0.91 / 0.89 / 0.42 | 0.45 / 0.49 / 0.58 |

With correlation 0.5 or a uniform model prior the two are indistinguishable;
finding the true model is identical (0.31 vs 0.31, n = 500). It slides along
the same trade-off rather than improving it.

Why: the prior median R² each prior claims as the model grows

| k | flat-R | growing | JZS |
|---|---|---|---|
| 1 | 0.25 | 0.25 | 0.11 |
| 4 | 0.25 | 0.42 | 0.48 |
| 12 | 0.25 | 0.47 | 0.76 |

JZS puts its scale on each coefficient, so every added predictor adds claimed
variance; a 12-predictor model promises R² ≈ 0.76 and is punished when the data
don't show it. Even flat on R² tops out at a median of 0.5, too close to flat-R
to matter. The gap to JZS is a difference in stance (the total explained
variance is unknown, whatever the size, versus each predictor brings its own
effect), not something a better a(k) schedule closes. Conclusion: keep flat-R
as is; put any wanted strictness in the prior over models, where it is visible
and adjustable. The method stays in 09 as a documented experiment.

### Uniform correlations in variable selection

If every valid correlation matrix of (y, x₁ … x_k) is equally likely (LKJ,
η = 1), the population R² is exactly Beta(k/2, η): flat-R at k = 1, flat R² at
k = 2, then growing (1 − R² is a product of independent partial-correlation
terms that telescopes; check V4, KS ≤ 0.013). This is a known result, and the
construction is rstanarm `stan_lm`'s prior (R² ~ Beta(k/2, η) with a uniform
direction), where η is chosen by the user and η = 1 is not singled out. The
same growth law a = k/2 is forced by requiring models of different sizes to
agree (a fixed g-mixture, B ≡ 1), and keeping every single correlation flat
(Barnard et al.'s marginally uniform prior) implies it too
(`docs/relation_to_others.md`).

Result: it lands on Zellner–Siow, the strict end, marginally ahead of it
(slightly fewer misses at the same false inclusions):

| n = 500, uncorrelated | false inclusions (beta-binomial / uniform) | missed (beta-binomial / uniform) | true model found (uniform) |
|---|---|---|---|
| flat-R | 0.91 / 0.60 | 0.45 / 0.43 | 0.37 |
| JZS default | 0.42 / 0.47 | 0.58 / 0.47 | 0.38 |
| uniform correlations | 0.13 / 0.20 | 0.84 / 0.64 | 0.36 |
| Zellner–Siow (1) | 0.12 / 0.17 | 0.89 / 0.67 | 0.33 |

At the 0.5 cutoff each principled stance sits at its own place: one claim
about the total (flat-R, sensitive), every correlation structure equally likely
(a = k/2, strict), and JZS in between (its median prior R² is close to LKJ with
η ≈ 2: 0.11 vs 0.12 at k = 1, 0.67 vs 0.69 at k = 8; medians only). Zellner–Siow
(scale 1) is exactly R² ~ Beta(k/2, ½), the LKJ η = ½ member, which is why it
lands next to uniform correlations.

### Equal or worse: whole trade-off curves

The 0.5 cutoff gives each method one point. Sliding the cutoff on the inclusion
probabilities traces a curve of missed against false inclusions per method
(figure `09_trade_off_curves.png`); the AUC for "real vs null predictor"
summarises it. Difference from flat-R, beta-binomial model prior, 95% bootstrap
over datasets (* = interval excludes 0):

| | n | JZS | uniform correlations | Zellner–Siow | BIC |
|---|---|---|---|---|---|
| uncorrelated | 50 | −0.025* | −0.019* | −0.021* | +0.005 |
| | 100 | −0.029* | −0.021* | −0.030* | +0.012 |
| | 200 | −0.006 | +0.006 | +0.001 | +0.028* |
| | 500 | +0.002 | +0.003 | +0.002 | +0.009* |
| correlation 0.5 | 50 | +0.001 | +0.011 | +0.011 | +0.030* |
| | 100 | +0.006 | +0.020* | +0.022* | +0.039* |
| | 200 | +0.007* | +0.022* | +0.023* | +0.035* |
| | 500 | +0.003* | +0.011* | +0.011* | +0.019* |

With the uniform model prior the differences are within about ±0.01 (BIC up
to +0.017 with correlated predictors).

- The curves nearly coincide; the priors mainly choose where on the curve the
  0.5 cutoff lands.
- Flat-R ranks predictors best in small samples with uncorrelated predictors;
  with correlated predictors it is 0.01–0.02 behind uniform correlations and
  Zellner–Siow (a null predictor that correlates with a real one can ride along
  when leaving the null model carries no extra penalty; cf. study 06).
- Uniform correlations is at least as good as JZS in every cell (by AUC,
  without a bootstrap for that pair).
- BIC ranks best or equal-best almost everywhere, while missing the most at its
  0.5 cutoff: ranking and calibration are separate properties.
- One design (12 candidates, 4 real, 100 datasets per cell): indicative.
- The rerun of 26 October 2026 reproduces every number in this section.

## 10: is the conditional test calibrated, and how far is it from one prior for everything?

**Calibration (P1–P2).** Study 06 found that the uniform-correlations
conditional test calls a slope of 0.1 null in about half the datasets. Is
that a fault? A Bayes factor is right if, among datasets that give it BF₁₀,
the share of true nulls is 1/(1 + BF₁₀) (prior odds 1:1), *when the data come
from the prior it assumes*. Datasets drawn from the joint uniform-correlations
prior, with x1's coefficient set to 0 in half of them, k = 8, n = 50 (the other
three cells are as close):

| log BF₁₀ | −2 to −1 | −1 to 0 | 0 to 1 | 1 to 2 | 2 to 4 | above 4 |
|---|---|---|---|---|---|---|
| datasets | 9,517 | 2,492 | 996 | 643 | 859 | 5,493 |
| share of nulls, observed | 0.823 | 0.638 | 0.409 | 0.191 | 0.054 | 0.000 |
| predicted, mean of 1/(1 + BF₁₀) | 0.823 | 0.643 | 0.394 | 0.191 | 0.058 | 0.001 |
| the same, data from JZS: observed | 0.743 | 0.567 | 0.352 | 0.195 | 0.052 | 0.002 |
| predicted | 0.821 | 0.640 | 0.393 | 0.193 | 0.061 | 0.001 |

The test is calibrated on its own prior, so its Bayes factors are real
posterior odds there, even though it is computed with flat priors on the other
coefficients. With data from the JZS prior (Cauchy, scale 0.354) it says
"null" too readily: JZS expects many tiny effects, which look like nulls. So
calibration holds relative to the prior, and the choice of prior is what
decides whether an effect of 0.1 should be called null.

*Why a slope of 0.1 gets called null.* Among datasets with BF₀₁ > 3 drawn from
the prior, 83–89% are true nulls. The real effects among them are tiny: at
n = 200, 68–70% have a partial correlation below 0.1 and 97–98% below 0.2;
at n = 50, 45–47% and 78–80%. Under a flat partial correlation an effect of
0.1 is rare (10% of the prior mass is below 0.1), and at these n it predicts
data much like no effect at all. BF₀₁ > 3 in such a dataset is the correct
posterior, given that prior.

**Coherence (P3–P5).** The conditional test puts flat priors on the other
coefficients, so it is not a ratio of marginal likelihoods under the joint
prior. The Savage–Dickey ratio is: the joint model's prior density of x1's
partial correlation at 0 (1/2) over its posterior density at 0 (`fit_regression()`,
20,000 draws, logspline), the Bayes factor against the joint model restricted
to a zero partial correlation. On 20 datasets (k = 3, 8; n = 50, 200; x1's t
about 0 to 2.5; the other coefficients drawn from the prior):

| | mean difference | mean absolute | largest absolute |
|---|---|---|---|
| conditional − Savage–Dickey (16 estimable) | −0.005 | 0.028 | 0.092 |
| joint test − Savage–Dickey (the same 16) | −0.008 | 0.226 | 0.422 |

For the other four (log BF₁₀ 2 to 9) there are too few draws near 0 for a
precise density; there the conditional test is within 0.6 of the Savage–Dickey
estimate, within 2.2 Monte Carlo SEs. So in practice the conditional
test is one prior for everything, to within a few percent of the Bayes factor;
the joint test (`bayes_factor_models()`, with against without x1) is further
away (up to a factor 1.5), because the reduced model's own uniform-correlations
prior is not the joint prior given a zero partial correlation.
Figures: `results/figures/10_calibration.png`, `10_coherence.png`.

## Not yet done

The list of open work, for the whole project, is in
[`docs/0_overview.md`](../docs/0_overview.md) ("Open").
