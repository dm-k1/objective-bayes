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
```

The scripts load the package from this folder with `pkgload::load_all()`.

**All at once:** `Rscript validation/run_all.R` runs every study, 00 to 10,
each in its own R process, and writes `results/run_summary.md` (the commit, R
and CmdStan versions, cores, study 00's table, and every PASS/FAIL check) plus
a log per study in `results/logs/`. `run_all.R quick` runs 00 and the
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
| `00_initial_simulation.R` | The first result: what happens to partial r, R and R² as k grows, in valid correlation matrices? | 200,000 matrices per k = 1–10; three versions (partial r flat, R flat, R² flat); F1–F3: each version flat where it should be; writes `results/00_first_result.md` | about a minute |
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
error if anything failed. Simulation results and figures are cached in
`validation/results/` (not tracked by git).

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

To be refilled from the rerun for three priors (`results/run_summary.md`).
Known so far, on the current code:

| Check | Criterion | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|---|
| 00 F1–F3 | KS distance < 0.01 | pass (F1, 0.0028) | pass (F2, 0.0021) | pass (F3, 0.0027) |
| 01 P-checks | KS distance < 0.03 | (rerun pending) | (rerun pending) | (rerun pending) |
| 01 S | rank bins, χ² p > 0.01 | (rerun pending) | (rerun pending) | (rerun pending) |
| 02 A1 | log BF within 0.001 of `integrate()`, 162 settings | pass, largest error 2.6e-07 | pass, 4.1e-05 | (rerun pending) |
| 02 A2 | identical up to 1e-8 | pass, 2.4e-14 | pass, 1.9e-14 | (rerun pending) |
| 02 A3 | within 0.05 of brute force (400,000 draws) | pass: 0.817 vs 0.797 (k = 3), −0.692 vs −0.689 (k = 1), 0.284 vs 0.308 (k = 5) | pass: 1.644 vs 1.629, −0.692 vs −0.689, 2.507 vs 2.524 | (rerun pending) |
| 02 A4 | within 0.002 of a second derivation | pass, 5.1e-05 | pass, 3.5e-04 | (rerun pending) |
| 02 B1 | BF within 0.001 of 1 at n = k + 1 | pass (k = 1–6) | pass (k = 1–6) | (rerun pending) |
| 02 B2, B2b | increasing, > 50 at R² = 0.9999; boundary n = k + 3 | pass (13.5 29.7 45.8 61.9) | pass (12.6 28.6 44.7 60.9) | (rerun pending) |
| 02 B3 | slope of median log BF₀₁ in 0.35–0.65 | pass (0.49) | pass (0.53) | (rerun pending) |
| 02 B3b | slopes within 20% of q/2 (smaller model x1), within 30% of a (intercept only) | pass: 0.47 1.49 2.93; 0.50 1.47 2.81 | pass: 0.54 1.55 2.91; 0.50 0.49 0.51 | (rerun pending) |
| 02 B4 | median log BF₁₀ increasing, > 20 at n = 5000 | pass (−0.0, 8.5, 53.3) | pass (0.9, 9.3, 54.1) | (rerun pending) |
| 05 C1–C5 | see the script | (rerun pending) | (rerun pending) | (rerun pending) |
| 05 C7–C11 | see the script | (rerun pending) | (rerun pending) | (rerun pending) |
| 08 G1–G2 | see the script | (rerun pending) | (rerun pending) | (rerun pending) |
| 09 V-checks | see the script | (rerun pending) | (rerun pending) | (rerun pending) |
| 10 P1–P5 | see the script | (rerun pending) | (rerun pending) | (rerun pending) |

## In short: what the studies say about the three priors

To be written from the rerun; the per-criterion comparison is
[`docs/4_comparison.md`](../docs/4_comparison.md).

| Question | Flat partial r | Flat R | Flat R-squared | Study |
|---|---|---|---|---|
| Separating nulls from small effects, one coefficient (AUC) | (rerun pending) | (rerun pending) | (rerun pending) | 04, 06 |
| Confirming a true null, one coefficient (BF₀₁ > 3) | (rerun pending) | (rerun pending) | (rerun pending) | 06 |
| Calling a small effect (0.1) null | (rerun pending) | (rerun pending) | (rerun pending) | 06 |
| Evidence for the null when k grows with n (k = n/4, n = 2000) | (rerun pending) | (rerun pending) | (rerun pending) | 08 |
| Variable selection | (rerun pending) | (rerun pending) | (rerun pending) | 09 |
| Estimation: prediction error against OLS | (rerun pending) | (rerun pending) | (rerun pending) | 07 |
| Estimation: 95% interval coverage, real effects | (rerun pending) | (rerun pending) | (rerun pending) | 07 |

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
- Figures: `results/figures/00_partial_r_R_R2_by_k.png` (the three
  distributions overlaid per k) and `00_three_versions_against_k.png`.
- This is why every later study compares three priors.

## 03: one coefficient, one design (200 datasets per cell)

Four predictors, x1's slope varies, x3 is null; n = 20 to 500; predictors
uncorrelated or correlated 0.5. The null predictor found, P(BF₀₁ > 3), and a
small real effect (x1 = 0.1) called null, P(BF₀₁ > 3):

| Method | null found, n = 200 unc. / corr. | null found, n = 500 unc. / corr. | 0.1 called null, n = 200 unc. / corr. | 0.1 called null, n = 500 unc. / corr. |
|---|---|---|---|---|
| flat partial r conditional | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R conditional | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R-squared conditional | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat partial r joint | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R joint | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R-squared joint | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| JZS default | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| Zellner–Siow | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |

## 04: discrimination, threshold-free

AUC (how often BF₀₁ is larger for a null dataset than for a small effect of
0.1), each method minus JZS, paired bootstrap 95% CI:

| | n = 100 | n = 200 | n = 500 |
|---|---|---|---|
| conditional tests (all three), uncorrelated | (rerun pending) | (rerun pending) | (rerun pending) |
| conditional tests (all three), correlated (0.5) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat partial r joint, correlated | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R joint, correlated | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R-squared joint, correlated | (rerun pending) | (rerun pending) | (rerun pending) |

- The conditional tests are expected to have identical AUCs: all three rank
  datasets by the coefficient's partial correlation, so they differ only in
  where their Bayes factors land relative to 3.
- M4 (why the joint tests differ): the change in log BF for x1 per unit of R²
  explained by the other predictors, at a fixed t-statistic: (rerun pending).

## 05: the conditional tests are computed correctly

Checks C1–C5 for every prior (equal to the joint test at k = 1; the implied
prior for one coefficient from correlated joint draws; the share mixture with
a = k/2 is a flat partial correlation; the grids against `integrate()` up to
n = 10⁶), C7 (why η = 1), C9–C11 (the block test): (rerun pending).

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
| conditional tests, correlation 0 | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| conditional tests, correlation 0.7 | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat partial r joint, correlation 0.7 | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R joint, correlation 0.7 | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R-squared joint, correlation 0.7 | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |

**Calibration at 3**, averaged over everything except the method and k:

| | k = 2 | k = 5 | k = 10 | k = 20 |
|---|---|---|---|---|
| True null found: flat partial r conditional | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| True null found: flat R conditional | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| True null found: flat R-squared conditional | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| True null found: JZS | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| Small effect (0.1) called null: flat partial r conditional | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| Small effect (0.1) called null: flat R conditional | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| Small effect (0.1) called null: flat R-squared conditional | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| Small effect (0.1) called null: JZS | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| Effect 0.2 found: flat partial r / flat R / flat R-squared / JZS | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |

False effects for a true null, largest over methods: (rerun pending).

## 07: estimation (24 designs, 100 datasets each)

Excess prediction error for a new observation, relative to OLS (below 1 =
better), mean over correlation and pattern; coverage of the 95% intervals;
mean bias of real effects:

| | prediction error, k = 3: n = 30 / 60 / 200 | k = 8: n = 30 / 60 / 200 | coverage, real effects, k = 3 / 8 | lowest cell, real effects | mean bias, real effects |
|---|---|---|---|---|---|
| flat partial r | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R-squared | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| JZS default | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| Zellner–Siow | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| OLS | 1 | 1 | (rerun pending) | (rerun pending) | (rerun pending) |

Divergent transitions over all Stan fits: (rerun pending).

## 08: k growing with n

Median log BF₀₁ for the true null model against the full model (200 datasets
per n):

| | n = 200 | n = 2000 |
|---|---|---|
| k = n/4: flat partial r / flat R / flat R-squared / JZS / Zellner–Siow | (rerun pending) | (rerun pending) |
| k = √n: flat partial r / flat R / flat R-squared / JZS / Zellner–Siow | (rerun pending) | (rerun pending) |

When the small model is true, the evidence for it at n = 2000, k = n/4:
(rerun pending).

## 09: variable selection (12 candidates, 4 real: 0.3, 0.2, 0.15, 0.1; 100 datasets per cell)

At the 0.5 cutoff (n = 500, uncorrelated, beta-binomial model prior):

| | real predictors missed (of 4) | null predictors included (of 8) | inclusion of the 0.1 effect | true model found (uniform model prior) |
|---|---|---|---|---|
| flat partial r | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R-squared | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| JZS default | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| Zellner–Siow (1) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| BIC | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |

**Whole trade-off curves.** Sliding the cutoff on the inclusion probabilities
traces a curve of missed against false inclusions per method (figure
`09_trade_off_curves.png`); the AUC for "real vs null predictor" summarises it.
Difference from flat R, beta-binomial model prior, 95% bootstrap over datasets:

| | n | JZS | flat partial r | flat R-squared | Zellner–Siow | BIC |
|---|---|---|---|---|---|---|
| uncorrelated | 50 / 100 / 200 / 500 | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |
| correlation 0.5 | 50 / 100 / 200 / 500 | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) | (rerun pending) |

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
| flat partial r | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R | (rerun pending) | (rerun pending) | (rerun pending) |
| flat R-squared | (rerun pending) | (rerun pending) | (rerun pending) |

**Coherence (P3–P5).** The conditional test puts flat priors on the other
coefficients, so it is not a ratio of marginal likelihoods under the joint
prior. The Savage–Dickey ratio is: the joint model's prior density of x1's
partial correlation at 0 over its posterior density at 0 (`fit_regression()`
draws, logspline). Conditional and joint tests minus Savage–Dickey, on 20
datasets (k = 3, 8; n = 50, 200):

| | conditional: mean absolute / largest | joint: mean absolute / largest |
|---|---|---|
| flat partial r | (rerun pending) | (rerun pending) |
| flat R | (rerun pending) | (rerun pending) |
| flat R-squared | (rerun pending) | (rerun pending) |

Figures: `results/figures/10_calibration.png`, `10_coherence.png`.

## Not yet done

The list of open work, for the whole project, is in
[`docs/0_overview.md`](../docs/0_overview.md) ("Open").

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
