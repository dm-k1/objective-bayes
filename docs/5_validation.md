# Validation

One script per question, for the three priors side by side: flat partial r
(`flat_partial_r`, $R^2 \sim$ Beta(k/2, 1)), flat R (`flat_r`, Beta(½, 1)) and
flat R-squared (`flat_r_squared`, Beta(1, 1)). Study 00 comes first: it is the
result that started the comparison. Run from the project folder:

```bash
Rscript validation/00_initial_simulation.R    # the first result; about a minute, no Stan
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
Rscript validation/11_factor_terms.R
```

The scripts (in `validation/`) load the package from the project folder with `pkgload::load_all()`.

**All at once:** `Rscript validation/run_all.R` runs every study, 00 to 11,
each in its own R process, and writes `validation/results/run_summary.md` (the commit, R
and CmdStan versions, cores, study 00's table, and every PASS/FAIL check) plus
a log per study in `cache/logs/`. `run_all.R quick` runs 00 and the
correctness checks 01, 02 and 05; `run_all.R 03 04 09` runs a selection.
`VALIDATION_CORES=32` sets the cores.

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
| `00_initial_simulation.R` | The first result: what happens to partial r, R and R² as k grows, in valid correlation matrices? | 200,000 matrices per k = 1–10; three versions (partial r flat, R flat, R² flat); F1–F3: each version flat where it should be; writes `validation/results/00_first_result.md` | about a minute |
| `01_prior_and_posterior.R` | Does the Stan model (`inst/stan/flatcor.stan`) get the prior and the posterior right, for each prior? | P1–P4 (flat R), P5–P7 (flat partial r), P8–P10 (flat R-squared): sample the prior alone and compare with its distribution. S: simulation-based calibration, k = 8, n = 60, 400 datasets, all three. Contrast: flat R sampled in coefficient coordinates | about 10 minutes |
| `02_bayes_factor.R` | Is the whole-model Bayes factor right, for each prior? | One loop over the three priors: A1 against `integrate()` by a route that is not the package's; A2 invariance; A3 against a brute-force Bayes factor; A4 a second derivation; B1–B4 predictive matching, information and model-selection consistency | a few minutes |
| `03_simulations.R` | Testing one coefficient across many datasets | 4 correlated predictors, x1's slope varies, x3 null; all one-coefficient methods (below) | 5–10 minutes |
| `04_matched_comparison.R` | Is a difference in 03 discrimination or calibration? | AUC (threshold-free), each method minus JZS with a paired bootstrap, matched thresholds, and why (M1–M4) | about 3 minutes |
| `05_conditional_test.R` | Are the conditional tests right, and why η = 1? | C1–C5 for every prior; C7 why η = 1; C8 against Kucharský et al. (2026); C9–C11 the block test | a few minutes |
| `06_broad_simulation.R` | Does the picture of 03–04 hold beyond one design? | k = 2–20, n = 50–1000, correlation 0–0.7, sparse and dense other coefficients; all one-coefficient methods | 20–30 minutes |
| `07_estimation.R` | How good are the estimates? | bias, coefficient error and prediction error relative to `lm()`, coverage of 95% intervals; the three priors, JZS, Zellner–Siow, `lm()`; resumes if interrupted | 3–5 hours |
| `08_growing_k.R` | Consistency when k grows with n (k = n/4, k = √n) | null true and small model true; all whole-model methods (G1–G2) | under a minute |
| `09_variable_selection.R` | Choosing among 4,096 subsets of 12 candidates (4 real) | all whole-model methods and BIC; uniform and beta-binomial model priors; trade-off curves and AUC; V-checks | 3–5 minutes |
| `10_calibration.R` | Is the conditional test calibrated on its own prior, and how far is it from one prior for everything? | P1–P2 calibration per prior, with data drawn from that prior and from JZS; P3–P5 against the Savage–Dickey ratio under the joint model | about 15 minutes |
| `comparison_methods.R` | helper: the common method set (names, colours, Bayes factors), the LKJ sampler, the g-mixture check | – | – |
| `settings.R` | helper: parallel cores, caching, `record_check()`, `save_figure()` | – | – |

## The criteria from the literature, and where each is checked

| Criterion | Source | Study and checks |
|---|---|---|
| Basic: proper prior on the tested coefficients, no arbitrary scale | Bayarri, Berger, Forte & García-Donato (2012) | by construction (`docs/2_how_it_calculates.md`, section 11); 01, P-checks (the prior integrates to its stated distribution) |
| Model selection consistency, fixed k | Bayarri et al. (2012) | 02, B3, B3b (the rate, with and without an effect in the smaller model), B4; proofs and sketches in `docs/3_consistency.md` |
| Model selection consistency, k growing with n | Bayarri et al. (2012) | 08, G1–G2 |
| Information consistency | Bayarri et al. (2012) | 02, B2, B2b (the boundary n = k + 3) |
| Intrinsic consistency | Bayarri et al. (2012) | 02, B5 (the prior with the sample covariance tends to the prior with the population covariance); proof in `docs/3_consistency.md` |
| Predictive matching (null and dimensional) | Bayarri et al. (2012) | 02, B1 (null: n = k + 1 gives BF = 1). Dimensional matching (two models of the same size at their minimal sample give BF = 1) holds by construction: the Bayes factor depends only on R², n and k, and at n = k + 1 both fit exactly |
| Measurement and group invariance | Bayarri et al. (2012) | 02, A2; by construction |
| Correctness of the computation | the g-prior closed form: Liang, Paulo, Molina, Clyde & Berger (2008) | 02, A1, A3, A4; 05, C1–C5, C9–C11 |
| Simulation-based calibration of the posterior | Talts, Betancourt, Simpson, Vehtari & Gelman (2018) | 01, S |
| Frequentist operating characteristics: nulls found, small effects called null, effects found, false effects, discrimination (AUC) | Schönbrodt & Wagenmakers (2018), Bayes factor design analysis | 03, 04, 06; variable selection 09; estimation 07 |
| Calibration of the Bayes factor (share of nulls = 1/(1 + BF₁₀) under the prior) | Schad, Nicenboim, Bürkner, Betancourt & Vasishth (2023) | 10, P1–P2 |
| Coherence with one prior for everything (Savage–Dickey) | Wagenmakers, Lodewyckx, Kuriyal & Grasman (2010); the conditional Lindley paradox it avoids: Som, Hans & MacEachern (2016) | 10, P3–P5 |

## The methods compared (03 on)

One set, defined once in `comparison_methods.R` and used by every study. The
three contenders keep the colours of study 00 (darker: conditional, lighter:
joint); the references are grey.

| Method | What it says | Whole models | One coefficient |
|---|---|---|---|
| **flat partial r** | R² ~ Beta(k/2, 1): every correlation matrix equally likely; `stan_lm`'s prior with η = 1 | ✓ | joint and conditional |
| **flat R** | R² ~ Beta(½, 1) for every k | ✓ | joint and conditional |
| **flat R-squared** | R² ~ Beta(1, 1) for every k | ✓ | joint and conditional |
| JZS default (reference) | Cauchy, scale 0.354 (BayesFactor, JASP) | ✓ | joint |
| Zellner–Siow (reference) | Cauchy, scale 1; R² ~ Beta(k/2, ½) | ✓ | joint |

- **Joint:** the model with against without the coefficient.
- **Conditional:** only the coefficient's partial correlation, with the prior
  the joint prior implies for one coefficient. For flat partial r that is a
  flat partial correlation whatever k is.
- Estimation (07): the three priors (all `fit_regression()`, the same Stan
  model with the shape a as data), JZS, Zellner–Siow and `lm()`; the
  conditional tests have no estimation version.
- Dropped: R2-D2 and the approximate flat partial correlation from 03, and the
  experiment R² ~ Beta(k/(k+1), 1) from 09, which flat R-squared now replaces
  as a full member (the code is in git history).

Each check prints PASS/FAIL with its criterion; 00, 01, 02 and 05 stop with an
error if anything failed. Simulation results, Stan fits and logs are cached in
`cache/` at the project root (not tracked by git); `validation/results/` holds
only what is tracked: the run summary, the first result and the figures.

**Simulation-based calibration** (01, S): draw the true parameters from the
prior, simulate a dataset, fit it, and record where the true value falls among
the posterior draws. If the sampler is right, every position is equally likely,
so the counts per bin should be flat (Talts et al., 2018).

**Why the Stan model samples R and a direction.** The prior as a density on the
coefficients has a sharp peak at zero (`docs/2_how_it_calculates.md`, section
3). A sampler moving through the coefficients, as brms does with a custom prior,
misses it for larger k, without divergent transitions: with brms at k = 8 the
10% quantile of the implied R² was 0.027 where the prior gives 0.010 (KS
distance 0.074; the earlier brms check, now in git history). Check 01 shows
the contrast.

## Last results

The run of 4 November 2026 at commit e221561, all studies (00–11), 7 cores
(`validation/results/run_summary.md`, logs in `cache/logs/`): **175 of 176 checks
pass**. The one failure is 10 P4 for flat R-squared, a finding rather than a
bug (below and in section 10).

| Check | Criterion | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|---|
| 00 F1–F3 | KS distance < 0.01 | pass (F1, 0.0028) | pass (F2, 0.0021) | pass (F3, 0.0027) |
| 01 P-checks | KS distance < 0.03 | pass (P5–P7, 0.003–0.012) | pass (P1–P4, 0.004–0.009) | pass (P8–P10, 0.006–0.009) |
| 01 S | rank bins, χ² p > 0.01 (R², b1, σ) | pass (p = 0.43, 0.94, 0.76) | pass (0.61, 0.86, 0.67) | pass (0.89, 0.91, 0.45) |
| 02 A1 | log BF within 0.001 of `integrate()`, 162 settings | pass, largest error 2.6e-07 | pass, 4.1e-05 | pass, 8.4e-05 |
| 02 A2 | identical up to 1e-8 | pass, 2.4e-14 | pass, 1.9e-14 | pass, 2.2e-14 |
| 02 A3 | within 0.05 of brute force (400,000 draws) | pass: 0.817 vs 0.797 (k = 3), −0.692 vs −0.689 (k = 1), 0.284 vs 0.308 (k = 5) | pass: 1.644 vs 1.629, −0.692 vs −0.689, 2.507 vs 2.524 | pass: 1.325 vs 1.308, −1.730 vs −1.727, 2.155 vs 2.175 |
| 02 A4 | within 0.002 of a second derivation | pass, 5.1e-05 | pass, 3.5e-04 | pass, 2.4e-05 |
| 02 B1 | BF within 0.001 of 1 at n = k + 1 | pass (k = 1–6) | pass (k = 1–6) | pass (k = 2–6) |
| 02 B2, B2b | increasing, > 50 at R² = 0.9999; boundary n = k + 3 | pass (13.5 29.7 45.8 61.9) | pass (12.6 28.6 44.7 60.9) | pass (13.2 29.3 45.4 61.5) |
| 02 B3 | slope of median log BF₀₁ in 0.35–0.65 | pass (0.49) | pass (0.53) | pass (0.51) |
| 02 B3b | slopes within 20% of q/2 (smaller model x1), within 30% of a (intercept only) | pass: 0.47 1.48 2.91; 0.51 1.46 2.86 | pass: 0.54 1.56 2.98; 0.47 0.50 0.48 | pass: 0.53 1.64 2.88; 0.99 1.00 1.01 |
| 02 B4 | median log BF₁₀ increasing, > 20 at n = 5000 | pass (−0.0, 8.5, 53.3) | pass (0.9, 9.3, 54.1) | pass (1.1, 9.5, 54.3) |
| 02 B5 | KS distance at n = 30, 300, 10⁵: < 0.01 at 10⁵ and falling | pass (0.034, 0.005, 0.005) | pass (0.026, 0.005, 0.005) | pass (0.032, 0.005, 0.005) |
| 05 C1–C5 | see the script | pass (C1 1.8e-06 for all; C4a KS ≤ 0.0025; C4b 1.8e-05; C5 3.5e-06) | pass (C2 KS 0.0024; C5 6.7e-05) | pass (C2 KS 0.0021; C5 5.0e-05) |
| 05 C7–C11 | see the script; C10 (1.0e-04) and C11 (1.6e-14) cover all three | pass (C7, C8; C9 KS ≤ 0.016) | pass (C9 KS 0.0033) | pass (C9 KS 0.0038) |
| 08 G1–G2 | median log BF for the truth rises with n (slope > 0.1) | pass (G1 slopes 98.8, 19.8; G2 108.4, 19.4) | pass (G1 0.24, 0.39; G2 43.1, 12.6) | pass (G1 0.49, 0.78; G2 43.0, 12.6) |
| 09 V1–V5 | tables and draws (shared by the methods) | pass (V3, V4: KS ≤ 0.013) | pass (V2: 7.7e-06) | pass (V1 tables: 2.8e-06) |
| 10 P1 | calibrated on its own prior: χ² p > 0.01, slope near 1 | pass (p 0.40–0.79; slopes 0.98–1.03) | pass (0.52–0.98; 0.99–1.02) | pass (0.37–0.94; 0.98–1.02) |
| 10 P2 | miscalibration shows with JZS data | pass (χ² 503–556) | pass (233–1198) | pass (50–265) |
| 10 P3–P5 | Savage–Dickey estimable and agreeing; conditional within a factor 2 | pass (P4 largest 0.092, factor 1.10) | not defined (see 10) | P3, P5 pass; **P4 fail** (largest 0.874, factor 2.40) |

**The one failure, 10 P4 for flat R-squared, is a finding, not a bug.** The
conditional test and the Savage–Dickey ratio under the joint model differ by
0.283 on the log scale on average (largest 0.874, a factor 2.40) on the 19
datasets where Savage–Dickey can be estimated. The difference is small at
k = 3 (mean −0.04 for n = 50 and n = 200) and grows with k (mean −0.48 at
k = 8, n = 50; −0.34 at k = 8, n = 200). For flat partial r it is at most
0.092 (factor 1.10). The Savage–Dickey ratio equals the conditional test only
when the prior of the other parameters given a zero partial correlation is
the joint prior's conditional. Under flat partial r (the LKJ vine) the partial
correlations are independent, so that holds. Under flat R-squared they are
dependent through R², so the joint posterior of x1's partial correlation also
draws information from the other coefficients; the more of them, the more it
does. The check stays as written: it records that flat R-squared's conditional
test is not "one prior for everything".

## In short: what the studies say about the three priors

The per-criterion comparison is [`docs/4_comparison.md`](4_comparison.md).

- **Correctness and the Bayarri et al. (2012) criteria** (01, 02, 05): all
  three pass every check: the sampler, SBC, the Bayes factors against
  independent routes, predictive matching, invariance, information
  consistency from n = k + 3, model-selection consistency at fixed k, and
  intrinsic consistency. Against the intercept-only model the evidence for a
  true null grows like a·log n: ½ log n for flat R, log n for flat R-squared,
  k/2 log n for flat partial r.
- **k growing with n** (08): all three are consistent, but for a true null at
  k = n/4 the evidence grows slowly under flat R (median log BF₀₁ 1.9 at
  n = 2000) and flat R-squared (3.7), quickly under flat partial r (408) and
  Zellner–Siow (411); JZS lies between (149).
- **One coefficient** (03, 04, 06): the three conditional tests rank datasets
  identically (same AUC); they differ in where the Bayes factor lands. Flat
  partial r confirms true nulls most often and does not depend on k; flat R
  and flat R-squared confirm fewer as k grows (at k = 20: 0.88, 0.06, 0.26)
  and call a small effect null less often. Their joint tests depend on what
  the other predictors explain (M4: −5.2 and −5.1 log BF per unit R², against
  −0.5 for flat partial r) and separate nulls from small effects worse than
  JZS under correlated predictors.
- **Calibration and coherence** (10): each is calibrated on data from its own
  prior. Only flat partial r's conditional test matches Savage–Dickey under
  its joint model; flat R-squared's drifts with k, and for flat R the ratio is
  not defined.
- **Variable selection and estimation** (07, 09): flat R and flat R-squared
  include more predictors (fewer real ones missed, more nulls included), and
  shrink more: lower prediction error at small n (k = 8, n = 30: 0.51 and
  0.53 of OLS, against 0.74 for flat partial r), but more bias in real effects
  and lower interval coverage in sparse designs (lowest cell 0.76 and 0.82,
  against 0.92). Flat partial r behaves much like Zellner–Siow.

| Question | Flat partial r | Flat R | Flat R-squared | Study |
|---|---|---|---|---|
| Separating nulls from small effects, one coefficient (AUC minus JZS, conditional, correlation 0.7, mean over k) | +0.038 | +0.038 | +0.038 | 04, 06 |
| Confirming a true null, one coefficient (BF₀₁ > 3), conditional, k = 2 / 20 | 0.86 / 0.88 | 0.55 / 0.06 | 0.86 / 0.26 | 06 |
| Calling a small effect (0.1) null, conditional, k = 2 / 20 | 0.52 / 0.60 | 0.19 / 0.005 | 0.52 / 0.04 | 06 |
| Evidence for the null when k grows with n (median log BF₀₁, k = n/4, n = 2000) | 407.7 | 1.9 | 3.7 | 08 |
| Variable selection (n = 500, uncorrelated): real missed of 4 / nulls included of 8 | 0.84 / 0.13 | 0.45 / 0.91 | 0.44 / 0.92 | 09 |
| Estimation: prediction error against OLS, k = 8, n = 30 / 200 | 0.74 / 0.95 | 0.51 / 0.94 | 0.53 / 0.92 | 07 |
| Estimation: 95% interval coverage, real effects, k = 3 / 8 (lowest cell) | 0.957 / 0.940 (0.92) | 0.947 / 0.928 (0.76) | 0.957 / 0.942 (0.82) | 07 |

## 00: the first result

200,000 valid correlation matrices of (y, x1 … xk) per k. Keep each partial r
flat and R² grows with k; keep R or R² flat and each partial r shrinks with k.

| k | P(\|partial r\| > 0.3) | P(R² > 0.5) | P(\|partial r\| > 0.3) | P(R² > 0.5) | P(\|partial r\| > 0.3) | P(R² > 0.5) |
|---|---|---|---|---|---|---|
| | partial r flat | partial r flat | R flat | R flat | R² flat | R² flat |
| 1 | 0.70 | 0.29 | 0.70 | 0.29 | 0.91 | 0.50 |
| 2 | 0.70 | 0.50 | 0.50 | 0.29 | 0.70 | 0.50 |
| 3 | 0.70 | 0.64 | 0.41 | 0.29 | 0.60 | 0.50 |
| 4 | 0.70 | 0.75 | 0.36 | 0.29 | 0.54 | 0.50 |
| 5 | 0.70 | 0.82 | 0.32 | 0.29 | 0.49 | 0.50 |
| 6 | 0.70 | 0.87 | 0.29 | 0.29 | 0.45 | 0.50 |
| 7 | 0.70 | 0.91 | 0.27 | 0.29 | 0.43 | 0.50 |
| 8 | 0.70 | 0.94 | 0.25 | 0.29 | 0.40 | 0.50 |
| 9 | 0.70 | 0.96 | 0.23 | 0.29 | 0.38 | 0.50 |
| 10 | 0.70 | 0.97 | 0.22 | 0.29 | 0.36 | 0.50 |

- At k = 1 partial r flat and R flat are the same prior (Jeffreys); R² flat is
  not flat on r even then.
- Figures: `validation/results/figures/00_partial_r_R_R2_by_k.png` (the three
  distributions overlaid per k) and `00_three_versions_against_k.png`.
- This is why every later study compares three priors.

## 03: one coefficient, one design (200 datasets per cell)

Four predictors, x1's slope varies, x3 is null; n = 20 to 500; predictors
uncorrelated or correlated 0.5. The null predictor found, P(BF₀₁ > 3), and a
small real effect (x1 = 0.1) called null, P(BF₀₁ > 3):

| Method | null found, n = 200 unc. / corr. | null found, n = 500 unc. / corr. | 0.1 called null, n = 200 unc. / corr. | 0.1 called null, n = 500 unc. / corr. |
|---|---|---|---|---|
| flat partial r conditional | 0.91 / 0.90 | 0.94 / 0.94 | 0.61 / 0.64 | 0.34 / 0.52 |
| flat R conditional | 0.41 / 0.40 | 0.70 / 0.72 | 0.16 / 0.20 | 0.12 / 0.21 |
| flat R-squared conditional | 0.80 / 0.79 | 0.89 / 0.88 | 0.46 / 0.50 | 0.23 / 0.38 |
| flat partial r joint | 0.82 / 0.83 | 0.91 / 0.90 | 0.47 / 0.51 | 0.25 / 0.40 |
| flat R joint | 0.40 / 0.54 | 0.73 / 0.78 | 0.07 / 0.15 | 0.08 / 0.21 |
| flat R-squared joint | 0.42 / 0.54 | 0.73 / 0.78 | 0.08 / 0.17 | 0.09 / 0.21 |
| JZS default | 0.61 / 0.65 | 0.80 / 0.82 | 0.22 / 0.26 | 0.12 / 0.27 |
| Zellner–Siow | 0.85 / 0.85 | 0.92 / 0.91 | 0.52 / 0.55 | 0.28 / 0.45 |

## 04: discrimination, threshold-free

AUC (how often BF₀₁ is larger for a null dataset than for a small effect of
0.1), each method minus JZS, paired bootstrap 95% CI:

| | n = 100 | n = 200 | n = 500 |
|---|---|---|---|
| conditional tests (all three), uncorrelated | +0.003 (−0.018, 0.023) | +0.003 (−0.010, 0.015) | −0.000 (−0.004, 0.003) |
| conditional tests (all three), correlated (0.5) | +0.047 (0.027, 0.067) | +0.027 (0.014, 0.040) | +0.028 (0.018, 0.039) |
| flat partial r joint, correlated | +0.029 (0.015, 0.042) | +0.018 (0.009, 0.028) | +0.022 (0.014, 0.031) |
| flat R joint, correlated | −0.039 (−0.058, −0.020) | −0.034 (−0.049, −0.019) | −0.028 (−0.036, −0.021) |
| flat R-squared joint, correlated | −0.032 (−0.047, −0.016) | −0.029 (−0.042, −0.017) | −0.027 (−0.035, −0.020) |

- The conditional tests are expected to have identical AUCs: all three rank
  datasets by the coefficient's partial correlation, so they differ only in
  where their Bayes factors land relative to 3.
- M4 (why the joint tests differ): the change in log BF for x1 per unit of R²
  explained by the other predictors, at a fixed t-statistic: about 0 for the
  three conditional tests (0.002, −0.011, −0.001); −0.51 flat partial r joint,
  −5.19 flat R joint, −5.10 flat R-squared joint; −1.93 JZS, −0.52
  Zellner–Siow.

## 05: the conditional tests are computed correctly

Checks C1–C5 for every prior (equal to the joint test at k = 1; the implied
prior for one coefficient from correlated joint draws; the share mixture with
a = k/2 is a flat partial correlation; the grids against `integrate()` up to
n = 10⁶), C7 (why η = 1), C9–C11 (the block test): all pass. C1 largest
difference 1.8e-06; C2 KS 0.0024 (flat R), 0.0021 (flat R-squared); C4a KS
≤ 0.0025; C4b 1.8e-05; C5 largest difference 3.5e-06 (flat partial r),
6.7e-05 (flat R), 5.0e-05 (flat R-squared); C9 KS ≤ 0.016; C10 1.0e-04;
C11 1.6e-14.

**C8: the flat partial r conditional test is not the Kucharský et al. (2026)
test, but differs from it in one assumption only.** (Unchanged code; numbers
from the run of 26 October 2026.) Their exact Bayes factor for a partial
correlation with a flat prior (eq. 2, α = 1) is

$$\mathrm{BF}_{10} = \frac{B(\tfrac12, 1 + h)}{B(\tfrac12, 1)}\,
{}_2F_1\!\left(h, h;\ \tfrac32 + h;\ r^2\right),\qquad h = \frac{n - k}{2},$$

with $k$ our number of predictors (their $k - 1$ controlling variables). It
reproduces their worked example (BF₁₀ = 45.88). Ours is lower in 47 of 48 cases
on the grid r = 0, 0.1, 0.3, 0.5; n = 20, 50, 200, 1000; k = 1, 3, 8 (log BF):

| | n = 20 | 50 | 200 | 1000 |
|---|---|---|---|---|
| k = 1, r = 0 | −0.054 | −0.023 | −0.006 | −0.001 |
| k = 3, r = 0 | −0.105 | −0.044 | −0.011 | −0.002 |
| k = 8, r = 0 | −0.262 | −0.098 | −0.024 | −0.005 |
| k = 8, r = 0.5 | −0.128 | −0.029 | −0.005 | −0.001 |

The difference is largest at r = 0, shrinks as 1/n and grows with k. Their
model is multivariate normal, so the predictors are random and the tested
predictor's residual SD (given the others) is an unknown parameter, with a
1/σ prior. Ours conditions on X, so the effect size f is in units of the
observed residual SD. Put their uncertainty into ours (f in sample units = f in
population units × √(W/n), W ~ χ²(n − k)) and the two agree to 5e-7 in all 48
cases.

## 06: across designs (k = 2–20, n = 50–1000, correlation 0–0.7, 216 cells)

**Discrimination** (AUC null vs 0.1, minus JZS, mean over n and the pattern of
the other coefficients):

| k | 2 | 5 | 10 | 20 |
|---|---|---|---|---|
| conditional tests, correlation 0 | +0.001 | +0.004 | +0.002 | +0.007 |
| conditional tests, correlation 0.7 | +0.035 | +0.045 | +0.041 | +0.031 |
| flat partial r joint, correlation 0.7 | +0.022 | +0.027 | +0.020 | +0.012 |
| flat R joint, correlation 0.7 | −0.064 | −0.035 | −0.022 | −0.013 |
| flat R-squared joint, correlation 0.7 | −0.059 | −0.033 | −0.020 | −0.012 |

**Calibration at 3**, averaged over everything except the method and k:

| | k = 2 | k = 5 | k = 10 | k = 20 |
|---|---|---|---|---|
| True null found: flat partial r conditional | 0.86 | 0.88 | 0.87 | 0.88 |
| True null found: flat R conditional | 0.55 | 0.29 | 0.20 | 0.06 |
| True null found: flat R-squared conditional | 0.86 | 0.65 | 0.47 | 0.26 |
| True null found: JZS | 0.54 | 0.38 | 0.30 | 0.30 |
| Small effect (0.1) called null: flat partial r conditional | 0.52 | 0.56 | 0.57 | 0.60 |
| Small effect (0.1) called null: flat R conditional | 0.19 | 0.05 | 0.02 | 0.005 |
| Small effect (0.1) called null: flat R-squared conditional | 0.52 | 0.30 | 0.16 | 0.04 |
| Small effect (0.1) called null: JZS | 0.21 | 0.12 | 0.08 | 0.11 |
| Effect 0.2 found: flat partial r / flat R / flat R-squared / JZS | 0.51 / 0.53 / 0.51 / 0.56 | 0.47 / 0.52 / 0.51 / 0.55 | 0.45 / 0.50 / 0.50 / 0.53 | 0.44 / 0.48 / 0.49 / 0.52 |

False effects for a true null (BF₁₀ > 3), largest over methods: 0.047 (flat
R-squared joint, k = 2); the conditional tests at most 0.018.

**What each test says** at 3 (share of datasets; conditional tests and JZS,
averaged over n, correlation and pattern). Correct answers in bold.

| Truth | k | Test | says null | undecided | says effect |
|---|---|---|---|---|---|
| null | 2 | flat partial r / flat R-squared | **0.86** | 0.13 | 0.01 |
| null | 2 | flat R | **0.55** | 0.44 | 0.02 |
| null | 2 | JZS | **0.54** | 0.44 | 0.02 |
| null | 20 | flat partial r | **0.88** | 0.11 | 0.01 |
| null | 20 | flat R-squared | **0.26** | 0.72 | 0.02 |
| null | 20 | flat R | **0.06** | 0.92 | 0.02 |
| null | 20 | JZS | **0.30** | 0.68 | 0.03 |
| slope 0.1 | 2 | flat partial r / flat R-squared | 0.52 | 0.29 | **0.19** |
| slope 0.1 | 2 | flat R | 0.19 | 0.58 | **0.23** |
| slope 0.1 | 20 | flat partial r | 0.60 | 0.26 | **0.15** |
| slope 0.1 | 20 | flat R-squared | 0.04 | 0.75 | **0.21** |
| slope 0.1 | 20 | flat R | 0.005 | 0.78 | **0.21** |
| slope 0.2 | 20 | flat partial r | 0.31 | 0.26 | **0.44** |
| slope 0.2 | 20 | flat R-squared | 0.00 | 0.52 | **0.49** |
| slope 0.2 | 20 | flat R | 0.00 | 0.52 | **0.48** |

The three conditional tests rank the datasets identically (AUC null vs 0.1:
0.728 at k = 2 and 0.701 at k = 20 for all three; JZS 0.713 and 0.684), so
none of them is more accurate. They differ only in where they decide. At
k = 20 flat R almost never calls a small effect null (0.005) because it
almost never calls anything null: 92% of true nulls stay undecided, and it
finds effects no more often than the others (0.21 for slope 0.1). Flat
partial r decides on most datasets, which finds 88% of true nulls but also
calls 60% of 0.1 effects and 31% of 0.2 effects null: with a small effect and
little data it leans to the null. Neither is a gain in accuracy; one is
silent, the other leans.

## 07: estimation (24 designs, 100 datasets each)

Excess prediction error for a new observation, relative to OLS (below 1 =
better), mean over correlation and pattern; coverage of the 95% intervals;
mean bias of real effects:

| | prediction error, k = 3: n = 30 / 60 / 200 | k = 8: n = 30 / 60 / 200 | coverage, real effects, k = 3 / 8 | lowest cell, real effects | mean bias, real effects |
|---|---|---|---|---|---|
| flat partial r | 0.83 / 0.89 / 0.97 | 0.74 / 0.84 / 0.95 | 0.957 / 0.940 | 0.916 | −0.016 (−8%) |
| flat R | 0.72 / 0.85 / 1.02 | 0.51 / 0.73 / 0.94 | 0.947 / 0.928 | 0.760 | −0.054 (−28%) |
| flat R-squared | 0.74 / 0.84 / 0.97 | 0.53 / 0.72 / 0.92 | 0.957 / 0.942 | 0.815 | −0.040 (−21%) |
| JZS default | 0.70 / 0.80 / 0.95 | 0.56 / 0.72 / 0.91 | 0.959 / 0.951 | 0.870 | −0.037 (−19%) |
| Zellner–Siow | 0.86 / 0.91 / 0.97 | 0.77 / 0.85 / 0.95 | 0.952 / 0.939 | 0.910 | −0.014 (−7%) |
| OLS | 1 | 1 | 0.953 / 0.945 | 0.920 | −0.000 (0%) |

The lowest cell for every method except flat partial r (k = 8, n = 30, dense)
and Zellner–Siow (k = 8, n = 30, sparse) is k = 8, n = 60, uncorrelated,
sparse. Mean bias of real effects: mean error, in brackets relative to the
true value.

Divergent transitions over all Stan fits: none (0 of 2400 fits for each of
the three priors).

## 08: k growing with n

Median log BF₀₁ for the true null model against the full model (200 datasets
per n):

| | n = 200 | n = 2000 |
|---|---|---|
| k = n/4: flat partial r / flat R / flat R-squared / JZS / Zellner–Siow | 39.2 / 1.3 / 2.6 / 15.1 / 41.2 | 407.7 / 1.9 / 3.7 / 149.3 / 410.9 |
| k = √n: flat partial r / flat R / flat R-squared / JZS / Zellner–Siow | 18.1 / 1.7 / 3.3 / 8.8 / 19.6 | 82.6 / 2.5 / 5.0 / 44.2 / 84.7 |

When the small model is true, the evidence for it at n = 2000, k = n/4
(median log BF, true vs full): 443.6 flat partial r, 175.1 flat R, 175.1 flat
R-squared, 245.5 JZS, 445.9 Zellner–Siow. Every method's evidence keeps
rising (G1, G2 pass for all three priors), but under flat R and flat
R-squared the evidence for a true null grows slowly: their prior puts R² of
the full model near 0 often enough that an empty model of 500 predictors is
hardly penalised.

## 09: variable selection (12 candidates, 4 real: 0.3, 0.2, 0.15, 0.1; 100 datasets per cell)

At the 0.5 cutoff (n = 500, uncorrelated, beta-binomial model prior):

| | real predictors missed (of 4) | null predictors included (of 8) | inclusion of the 0.1 effect | true model found (uniform model prior) |
|---|---|---|---|---|
| flat partial r | 0.84 | 0.13 | 0.45 | 0.37 |
| flat R | 0.45 | 0.91 | 0.64 | 0.36 |
| flat R-squared | 0.44 | 0.92 | 0.65 | 0.35 |
| JZS default | 0.58 | 0.42 | 0.58 | 0.37 |
| Zellner–Siow (1) | 0.89 | 0.12 | 0.43 | 0.35 |
| BIC | 1.01 | 0.05 | 0.34 | 0.27 |

**Whole trade-off curves.** Sliding the cutoff on the inclusion probabilities
traces a curve of missed against false inclusions per method (figure
`09_trade_off_curves.png`); the AUC for "real vs null predictor" summarises it.
Difference from flat R, beta-binomial model prior, 95% bootstrap over datasets:

| | n | JZS | flat partial r | flat R-squared | Zellner–Siow | BIC |
|---|---|---|---|---|---|---|
| uncorrelated | 50 / 100 / 200 / 500 | −0.025* / −0.029* / −0.006 / +0.002 | −0.019* / −0.021* / +0.006 / +0.003 | −0.017* / −0.018* / −0.004 / +0.000* | −0.021* / −0.030* / +0.001 / +0.002 | +0.005 / +0.012 / +0.028* / +0.009* |
| correlation 0.5 | 50 / 100 / 200 / 500 | +0.001 / +0.006 / +0.007* / +0.003* | +0.011 / +0.020* / +0.022* / +0.011* | +0.000 / −0.000 / −0.000 / +0.000 | +0.011 / +0.022* / +0.023* / +0.011* | +0.030* / +0.039* / +0.035* / +0.019* |

\* the 95% interval excludes 0.

- The earlier experiment R² ~ Beta(k/(k+1), 1), growing from flat R toward
  flat R-squared, is superseded by flat R-squared as a full member.
- Against the intercept-only model the evidence for the null grows like
  a log n (½ log n for flat R, log n for flat R-squared, k/2 log n for flat
  partial r; `docs/3_consistency.md`), so how inclusive each prior is follows
  from its shape a; the prior over models sets the rest.

## 10: is the conditional test calibrated, and how far is it from one prior for everything?

**Calibration (P1–P2).** A Bayes factor is right if, among datasets that give
it BF₁₀, the share of true nulls is 1/(1 + BF₁₀) (prior odds 1:1), *when the
data come from the prior it assumes*. Datasets are drawn from each prior's
joint version, with x1's coefficient set to 0 in half of them; and from the
JZS prior, where miscalibration should show:

| | χ² p over bins (own prior) | logistic slope (own prior) | JZS data: χ² |
|---|---|---|---|
| flat partial r | 0.40–0.79 | 0.98–1.03 | 503–556 (slopes 0.89–0.95) |
| flat R | 0.52–0.98 | 0.99–1.02 | 233–1198 (slopes 1.16–1.41) |
| flat R-squared | 0.37–0.94 | 0.98–1.02 | 50–265 (slopes 0.95–1.12) |

**Coherence (P3–P5).** The conditional test puts flat priors on the other
coefficients, so it is not a ratio of marginal likelihoods under the joint
prior. The Savage–Dickey ratio is: the joint model's prior density of x1's
partial correlation at 0 over its posterior density at 0 (`fit_regression()`
draws, logspline). Conditional and joint tests minus Savage–Dickey, on 20
datasets (k = 3, 8; n = 50, 200):

| | conditional: mean absolute / largest | joint: mean absolute / largest |
|---|---|---|
| flat partial r | 0.028 / 0.092 (16 of 20 estimable) | 0.226 / 0.422 |
| flat R | not defined: infinite prior density at 0 | not defined |
| flat R-squared | 0.283 / 0.874 (19 of 20 estimable) | 1.007 / 2.365 |

P1 ranges are over k = 3, 8 and n = 50, 200 (χ² over 8–11 bins); all four
cells pass for every prior. The coherence numbers are for the datasets where
Savage–Dickey can be estimated (chain SE < 0.1); the others agree within 3
Monte Carlo SE (P5). For flat R-squared the conditional test sits below
Savage–Dickey and more so at larger k (mean −0.04 at k = 3; −0.48 and −0.34
at k = 8, n = 50 and 200), so P4 (within a factor 2) fails: its partial
correlations are dependent through R², and the joint posterior of x1's
partial correlation also uses the other coefficients (section Last results).

Figures: `validation/results/figures/10_calibration.png`, `10_coherence.png`.

## Not yet done

The list of open work, for the whole project, is in
[`docs/0_overview.md`](0_overview.md) ("Open").

## References

- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012). Criteria
  for Bayesian model choice with application to variable selection. *The Annals
  of Statistics, 40*(3), 1550–1577.
- Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
  correlation matrices based on vines and extended onion method. *Journal of
  Multivariate Analysis, 100*(9), 1989–2001.
- Liang, F., Paulo, R., Molina, G., Clyde, M. A., & Berger, J. O. (2008).
  Mixtures of g priors for Bayesian variable selection. *Journal of the American
  Statistical Association, 103*(481), 410–423.
- Schad, D. J., Nicenboim, B., Bürkner, P.-C., Betancourt, M., & Vasishth, S.
  (2023). Workflow techniques for the robust use of Bayes factors.
  *Psychological Methods, 28*(6), 1404–1426.
- Schönbrodt, F. D., & Wagenmakers, E.-J. (2018). Bayes factor design analysis:
  Planning for compelling evidence. *Psychonomic Bulletin & Review, 25*(1),
  128–142.
- Som, A., Hans, C. M., & MacEachern, S. N. (2016). A conditional Lindley
  paradox in Bayesian linear models. *Biometrika, 103*(4), 993–999.
- Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018).
  Validating Bayesian inference algorithms with simulation-based calibration.
  arXiv:1804.06788.
- Wagenmakers, E.-J., Lodewyckx, T., Kuriyal, H., & Grasman, R. (2010). Bayesian
  hypothesis testing for psychologists: A tutorial on the Savage–Dickey method.
  *Cognitive Psychology, 60*(3), 158–189.
