# Validation

Each script in `validation/` answers one question, for the three priors side
by side: flat partial r (`flat_partial_r`, $R^2 \sim$ Beta(*k*/2, 1)), flat R
(`flat_r`, Beta(½, 1)) and flat R-squared (`flat_r_squared`, Beta(1, 1)).
Study 00 comes first because it is the result that started the comparison.
Run the scripts from the project folder:

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

The scripts load the package from the project folder with
`pkgload::load_all()`.

`Rscript validation/run_all.R` runs every study, 00 to 11, each in its own R
process. It writes `validation/results/run_summary.md` (the commit, the R and
CmdStan versions, the number of cores, study 00's table and every PASS/FAIL
check) and a log per study in `cache/logs/`. `quarto render
validation/report.qmd` then builds `validation/report.html`, which presents
every figure and check of the run with a short explanation of how to read
each. `run_all.R quick` runs study 00 and the correctness checks 01, 02 and
05; `run_all.R 03 04 09` runs a selection; `VALIDATION_CORES=32` sets the
number of cores.

### Running it on another machine

The environment is pinned in two files: `renv.lock` records the R packages,
and `validation/Dockerfile` builds R, those packages and CmdStan into one
image. renv is not activated in this folder; the lockfile is the record, and
the image installs from it.

| Where | How |
|---|---|
| Any machine with Docker (a cloud VM, a SURF Research Cloud workspace) | `git clone` the repository, then `docker build -f validation/Dockerfile -t flatcor-validation .` (once, 20–40 minutes) and `docker run --rm -v "$PWD":/project -e VALIDATION_CORES=32 flatcor-validation`. Results land in `validation/results/` of the clone |
| A machine with R but no Docker | `Rscript -e 'install.packages("renv"); renv::restore(library = .libPaths()[1], prompt = FALSE)'` (into your normal library), then `Rscript -e 'cmdstanr::install_cmdstan(version = "2.40.0")'`, then `Rscript validation/run_all.R` |
| A cluster without Docker (Snellius) | the same as the previous row, or the image through Apptainer: `apptainer build flatcor.sif docker-daemon://flatcor-validation:latest` |

To update the lockfile after installing or upgrading packages, run
`Rscript -e 'renv::snapshot(type = "implicit", exclude = c("flatcor", "rstan"), prompt = FALSE)'`.
rstan is left out because it is only the fallback when cmdstanr is missing,
and leaving it out keeps the image small.

The run times below are for a run without caches (the first full run, on 4
November 2026). With the caches in `cache/`, a rerun takes about 20 minutes.

| Script | Question | How | Run time |
|---|---|---|---|
| `00_initial_simulation.R` | The first result: what happens to partial r, R and R² as k grows, in valid correlation matrices? | 200,000 matrices per k = 1–10; three versions (partial r flat, R flat, R² flat); F1–F3: each version flat where it should be; writes `validation/results/00_first_result.md` | about a minute |
| `01_prior_and_posterior.R` | Does the Stan model (`inst/stan/flatcor.stan`) get the prior and the posterior right, for each prior? | P1–P4 (flat R), P5–P7 (flat partial r), P8–P10 (flat R-squared): sample the prior alone and compare with its distribution. S: simulation-based calibration, k = 8, n = 60, 400 datasets, all three. Contrast: flat R sampled in coefficient coordinates | about 8 minutes |
| `02_bayes_factor.R` | Is the whole-model Bayes factor right, for each prior? | One loop over the three priors: A1 against `integrate()` by a route that is not the package's; A2 invariance; A3 against a brute-force Bayes factor; A4 a second derivation; B1–B4 predictive matching, information and model-selection consistency; B5 intrinsic consistency | about a minute |
| `03_simulations.R` | Testing one coefficient across many datasets | 4 correlated predictors, x1's slope varies, x3 null; all one-coefficient methods (below); U1–U2 | 10–15 minutes |
| `04_matched_comparison.R` | Is a difference in 03 discrimination or calibration? | AUC (threshold-free), each method minus JZS with a paired bootstrap, matched thresholds, and why (M1–M4) | about 3 minutes |
| `05_conditional_test.R` | Are the conditional tests right, and why η = 1? | C1–C5 for every prior; C7 why η = 1; C8 against Kucharský et al. (2026); C9–C11 the block test | under a minute |
| `06_broad_simulation.R` | Does the picture of 03–04 hold beyond one design? | k = 2–20, n = 50–1000, correlation 0–0.7, sparse and dense other coefficients; all one-coefficient methods; every outcome (says null, undecided, says effect) at k = 2 and 20 | 20–30 minutes |
| `07_estimation.R` | How good are the estimates? | bias, coefficient error and prediction error relative to `lm()`, coverage of 95% intervals; the three priors, JZS, Zellner–Siow, `lm()`; resumes if interrupted | about 1.5 hours |
| `08_growing_k.R` | Consistency when k grows with n (k = n/4, k = √n) | null true and small model true; all whole-model methods (G1–G2) | about a minute |
| `09_variable_selection.R` | Choosing among 4,096 subsets of 12 candidates (4 real) | all whole-model methods and BIC; uniform and beta-binomial model priors; trade-off curves and AUC; V-checks | 3–5 minutes |
| `10_calibration.R` | Is the conditional test calibrated on its own prior, and how far is it from one prior for everything? | P1–P2 calibration per prior, with data drawn from that prior and from JZS; P3–P5 against the Savage–Dickey ratio under the joint model; P6–P7 what a zero partial correlation means under each joint prior | about 70 minutes |
| `11_factor_terms.R` | Testing a factor (a term of q columns) | the block conditional test of the three priors against their joint tests, JZS (lmBF) and Zellner–Siow; q = 2–5, n = 50–1000, what the other predictors explain 0–0.9; F1–F3 | about an hour |
| `comparison_methods.R` | helper: the common method set (names, Bayes factors), the figure style (`method_scales()`), the LKJ sampler, the g-mixture check | – | – |
| `settings.R` | helper: parallel cores, caching, `record_check()`, `save_figure()` | – | – |

## Which criterion from the literature is checked where

| Criterion | Source | Study and checks |
|---|---|---|
| Basic: proper prior on the tested coefficients, no arbitrary scale | Bayarri, Berger, Forte & García-Donato (2012) | by construction (`docs/2_how_it_calculates.md`, section 11); 01, P-checks (the prior integrates to its stated distribution) |
| Model selection consistency, fixed k | Bayarri et al. (2012) | 02, B3, B3b (the rate, with and without an effect in the smaller model), B4 |
| Model selection consistency, k growing with n | Bayarri et al. (2012) | 08, G1–G2 |
| Information consistency | Bayarri et al. (2012) | 02, B2, B2b (the boundary n = k + 3) |
| Intrinsic consistency | Bayarri et al. (2012) | 02, B5 (the prior with the sample covariance tends to the prior with the population covariance) |
| Predictive matching (null and dimensional) | Bayarri et al. (2012) | 02, B1 (null: n = k + 1 gives BF = 1). Dimensional matching (two models of the same size at their minimal sample give BF = 1) holds by construction: the Bayes factor depends only on R², n and k, and at n = k + 1 both fit exactly |
| Measurement and group invariance | Bayarri et al. (2012) | 02, A2; by construction |
| Correctness of the computation | the g-prior closed form: Liang, Paulo, Molina, Clyde & Berger (2008) | 02, A1, A3, A4; 05, C1–C5, C9–C11 |
| Simulation-based calibration of the posterior | Talts, Betancourt, Simpson, Vehtari & Gelman (2018) | 01, S |
| Frequentist operating characteristics: nulls found, small effects called null, effects found, false effects, discrimination (AUC) | Schönbrodt & Wagenmakers (2018), Bayes factor design analysis | 03, 04, 06; variable selection 09; estimation 07 |
| Calibration of the Bayes factor (share of nulls = 1/(1 + BF₁₀) under the prior) | Schad, Nicenboim, Bürkner, Betancourt & Vasishth (2023) | 10, P1–P2 |
| Coherence with one prior for everything (Savage–Dickey) | Wagenmakers, Lodewyckx, Kuriyal & Grasman (2010); the conditional Lindley paradox it avoids: Som, Hans & MacEachern (2016) | 10, P3–P7 |

## The methods compared (from study 03 on)

Every study from 03 on uses one set of methods, defined once in
`comparison_methods.R`. In the figures each prior keeps its colour from study
00, with its conditional test drawn solid and its joint test dashed; the two
reference methods are grey.

| Method | What it says | Whole models | One coefficient |
|---|---|---|---|
| **flat partial r** | R² ~ Beta(k/2, 1): every correlation matrix equally likely; `stan_lm`'s prior with η = 1 | ✓ | joint and conditional |
| **flat R** | R² ~ Beta(½, 1) for every k | ✓ | joint and conditional |
| **flat R-squared** | R² ~ Beta(1, 1) for every k | ✓ | joint and conditional |
| JZS default (reference) | Cauchy, scale 0.354 (BayesFactor, JASP) | ✓ | joint |
| Zellner–Siow (reference) | Cauchy, scale 1; R² ~ Beta(k/2, ½) | ✓ | joint |

- The **joint test** (= the model with the coefficient against the model
  without it) compares two whole models.
- The **conditional test** (= a test of the coefficient's partial
  correlation alone, with the prior that the joint prior implies for one
  coefficient) looks only at that coefficient. Under flat partial r that
  implied prior is a flat partial correlation for every *k*.
- Study 07 (estimation) fits the three priors with `fit_regression()` (one
  Stan model with the shape *a* as data), JZS, Zellner–Siow and `lm()`. The
  conditional tests have no estimation version.

Each check prints PASS or FAIL together with its criterion; studies 00, 01,
02 and 05 stop with an error if any check fails. Simulation results, Stan
fits and logs are cached in `cache/` at the project root, which git does not
track. `validation/results/` holds only tracked files: the run summary, the
first result, the session information and the figures.

**Simulation-based calibration** (study 01, check S) works as follows: draw
the true parameters from the prior, simulate a dataset, fit it, and record
the rank of the true value among the posterior draws. If the sampler is
correct, every rank is equally likely, so the counts per bin should be flat
(Talts et al., 2018). In this run every bin test passes (χ² *p* > .01); the
closest call is flat partial r's σ, with *p* = .020.

The Stan model samples *R* and a direction, not the coefficients, because the
prior has a sharp peak at zero when written as a density on the coefficients
(`docs/2_how_it_calculates.md`, section 3). A sampler that moves through the
coefficients directly can miss that peak without any warning. Checks P1–P10
of study 01 show that the model reproduces each prior (KS ≤ .012). In a
separate diagnostic run of the same calibration setting (*k* = 8, *n* = 60,
400 datasets per prior), the package's model had no divergent transitions in
any of its 1,200 fits, whereas the contrast model, sampled in coefficient
coordinates, had every transition divergent in 165 of its 400 fits.

## Last results

The run of 4 November 2026 at commit 21d42dc (started 12:24) covered all
studies (00–11) on 7 cores, and **137 of 138 checks passed**
(`validation/results/run_summary.md`; logs in `cache/logs/`). Studies 03, 06,
07, 09, 10 (parts 1 and 2) and 11 reused the simulations cached by the first
full run; because each design cell has its own seed, these are the same
simulations. Studies 00, 01, 02, 04, 05, 08 and 10 (part 3) were computed
afresh. The one failure, check P4 of study 10 for flat R-squared, records a
property of that prior rather than a bug; it is explained after the table
and in section 10.

| Check | Criterion | Flat partial r | Flat R | Flat R-squared |
|---|---|---|---|---|
| 00 F1–F3 | KS distance < 0.01 | pass (F1, 0.0026) | pass (F2, 0.0024) | pass (F3, 0.0038) |
| 01 P-checks | KS distance < 0.03 | pass (P5–P7, 0.003–0.012) | pass (P1–P4, 0.004–0.009) | pass (P8–P10, 0.006–0.009) |
| 01 S | rank bins, χ² p > 0.01 (R², b1, σ) | pass (p = 0.64, 0.21, 0.020) | pass (0.38, 0.63, 0.42) | pass (0.51, 0.82, 0.22) |
| 02 A1 | log BF within 0.001 of `integrate()`, 162 settings | pass, largest error 2.6e-07 | pass, 4.1e-05 | pass, 8.4e-05 |
| 02 A2 | identical up to 1e-8 | pass, 2.4e-14 | pass, 1.9e-14 | pass, 2.2e-14 |
| 02 A3 | within 0.05 of brute force (400,000 draws) | pass: 0.817 vs 0.797 (k = 3), −0.692 vs −0.689 (k = 1), 0.284 vs 0.308 (k = 5) | pass: 1.644 vs 1.629, −0.692 vs −0.689, 2.507 vs 2.524 | pass: 1.325 vs 1.308, −1.730 vs −1.727, 2.155 vs 2.175 |
| 02 A4 | within 0.002 of a second derivation | pass, 5.1e-05 | pass, 3.5e-04 | pass, 2.4e-05 |
| 02 B1 | BF within 0.001 of 1 at n = k + 1 | pass (k = 1–6) | pass (k = 1–6) | pass (k = 2–6) |
| 02 B2, B2b | increasing, > 50 at R² = 0.9999; boundary n = k + 3 | pass (13.5 29.7 45.8 61.9) | pass (12.6 28.6 44.7 60.9) | pass (13.2 29.3 45.4 61.5) |
| 02 B3 | slope of median log BF₀₁ in 0.35–0.65 | pass (0.49) | pass (0.53) | pass (0.51) |
| 02 B3b | slopes within 20% of q/2 (smaller model x1), within 30% of a (intercept only) | pass: 0.50 1.47 2.95; 0.51 1.44 2.90 | pass: 0.56 1.58 2.97; 0.51 0.48 0.50 | pass: 0.55 1.55 2.95; 1.01 0.96 0.99 |
| 02 B4 | median log BF₁₀ increasing, > 20 at n = 5000 | pass (−0.0, 8.5, 53.3) | pass (0.9, 9.3, 54.1) | pass (1.1, 9.5, 54.3) |
| 02 B5 | KS distance at n = 30, 300, 10⁵: < 0.01 at 10⁵, and at n = 30 more than twice that | pass (0.034, 0.005, 0.005) | pass (0.026, 0.005, 0.005) | pass (0.032, 0.005, 0.005) |
| 03 U1–U2 | k = 1 agreement within 1e-3; g-prior Monte Carlo within 0.01 | pass (U1 4.9e-06; U2 0.0012) | pass (U1, shared with flat partial r) | pass (U2 0.0015; U1 not applicable: at k = 1 a different prior) |
| 05 C1–C5 | see the script | pass (C1 1.8e-06 for all; C4a KS ≤ 0.0025; C4b 1.8e-05; C5 3.5e-06) | pass (C2 KS 0.0024; C5 6.7e-05) | pass (C2 KS 0.0021; C5 5.0e-05) |
| 05 C7–C11 | see the script; C10 (1.0e-04) and C11 (1.6e-14) cover all three | pass (C7, C8; C9 KS ≤ 0.017) | pass (C9 KS 0.0033) | pass (C9 KS 0.0038) |
| 08 G1–G2 | median log BF for the truth rises with n (slope > 0.1) | pass (G1 slopes 98.9, 20.0; G2 108.2, 19.7) | pass (G1 0.25, 0.40; G2 42.7, 12.9) | pass (G1 0.50, 0.80; G2 42.6, 12.8) |
| 09 V1–V5 | tables and draws (shared by the methods) | pass (V3, V4: KS ≤ 0.013) | pass (V2: 7.7e-06) | pass (V1 tables: 2.8e-06) |
| 10 P1 | calibrated on its own prior: χ² p > 0.01, slope near 1 | pass (p 0.40–0.79; slopes 0.98–1.03) | pass (0.52–0.98; 0.99–1.02) | pass (0.37–0.94; 0.98–1.02) |
| 10 P2 | miscalibration shows with JZS data | pass (χ² 503–556) | pass (233–1198) | pass (50–265) |
| 10 P3–P5 | Savage–Dickey estimable and agreeing; conditional within a factor 2 | pass (P4 largest 0.092, factor 1.10) | not estimable: infinite prior density at 0 (see P7) | P3, P5 pass; **P4 fail** (largest 0.874, factor 2.40) |
| 10 P6–P7 | the joint prior given \|ρ₁\| < ε, ε = 0.1, 0.02, 0.005: density near 0 settles (P6) or grows and the others' R² goes to 0 (P7) | pass (P6: density 0.498 0.502 0.485 at k = 3, 0.504 0.497 0.491 at k = 8; others' R² 0.37, 0.73) | pass (P7: density 1.51 2.30 2.99 at k = 3, 2.23 3.83 5.20 at k = 8; others' median R² down to 0.001 and 0.003) | pass (P6: density 0.741 0.767 0.781 at k = 3, 1.342 1.537 1.558 at k = 8; others' R² 0.17–0.22) |
| 11 F1–F3 | block conditional test unchanged by what the others explain; same AUC for the three; JZS (lmBF) precise | pass (F1 1.4e-12 over all three; F2 identical AUC in 72 cells) | pass | pass (F3, shared: largest error 0.032) |

**Why check P4 of study 10 fails for flat R-squared.** On the 19 datasets
where the Savage–Dickey ratio can be estimated, flat R-squared's conditional
test differs from it by 0.283 on the log scale on average (largest 0.874, a
factor of 2.40). The difference is small at *k* = 3 (mean −0.04 for both
*n* = 50 and *n* = 200) and grows with *k* (mean −0.48 at *k* = 8, *n* = 50;
−0.34 at *k* = 8, *n* = 200). For flat partial r it is at most 0.092 (a factor
of 1.10). The Savage–Dickey ratio equals the conditional test only when the
prior of the other parameters, given a zero partial correlation, is the
joint prior's conditional. Under flat partial r the partial correlations are
independent (the LKJ vine), so this holds. Under flat R-squared they are
linked through *R*², so the joint posterior of *x*₁'s partial correlation also
draws on the other coefficients, and more so the more of them there are. The
check stays as written, because it records that flat R-squared's conditional
test does not correspond to one prior for the whole model.

## Summary of what the studies show about the three priors

The comparison criterion by criterion is in
[`docs/3_comparison.md`](3_comparison.md).

- **Correctness and the Bayarri et al. (2012) criteria** (studies 01, 02,
  05). All three priors pass every check: the sampler, simulation-based
  calibration, the Bayes factors against independent routes, predictive
  matching, invariance, information consistency from *n* = *k* + 3,
  model-selection consistency at fixed *k*, and intrinsic consistency.
  Against the intercept-only model, the evidence for a true null grows like
  *a* log *n*: ½ log *n* for flat R, log *n* for flat R-squared and
  *k*/2 log *n* for flat partial r.
- **Number of predictors growing with the sample size** (study 08). All three
  are consistent, but for a true null with *k* = *n*/4 the evidence grows
  slowly under flat R (median log BF₀₁ of 1.9 at *n* = 2,000) and flat
  R-squared (3.7), and quickly under flat partial r (408) and Zellner–Siow
  (411); JZS lies in between (149).
- **Testing one coefficient** (studies 03, 04, 06). The three conditional
  tests rank datasets identically (the same AUC) and differ only in where
  the Bayes factor lands. Flat partial r confirms true nulls most often, and
  its rate does not depend on *k*. Flat R and flat R-squared confirm fewer
  true nulls as *k* grows (at *k* = 20: .88, .06 and .26) and call a small
  effect null less often. Their joint tests depend on how much the other
  predictors explain (study 04, M4: −5.2 and −5.1 log BF per unit of the
  others' *R*², against −0.5 for flat partial r), and with correlated
  predictors they separate nulls from small effects less well than JZS.
- **Calibration and coherence** (study 10). Each prior is calibrated on data
  drawn from itself. Only flat partial r's conditional test matches the
  Savage–Dickey ratio under its joint model; flat R-squared's moves away from
  it as *k* grows. Under flat R, a zero partial correlation for *x*₁ drags
  every other coefficient to zero (check P7), so no joint model corresponds
  to its test of one coefficient.
- **Variable selection and estimation** (studies 07, 09). Flat R and flat
  R-squared include more predictors (fewer real ones missed, more nulls
  included) and shrink the estimates more. This gives them a lower
  prediction error at small *n* (*k* = 8, *n* = 30: 0.51 and 0.53 times that of
  ordinary least squares, against 0.74 for flat partial r), but more bias in
  real effects and lower interval coverage in sparse designs (lowest cell .76
  and .82, against .92). Flat partial r behaves much like Zellner–Siow.

| Question | Flat partial r | Flat R | Flat R-squared | Study |
|---|---|---|---|---|
| Separating nulls from small effects, one coefficient (AUC minus JZS, conditional, correlation 0.7, mean over k) | +0.038 | +0.038 | +0.038 | 04, 06 |
| Confirming a true null, one coefficient (BF₀₁ > 3), conditional, k = 2 / 20 | 0.86 / 0.88 | 0.55 / 0.06 | 0.86 / 0.26 | 06 |
| Calling a small effect (0.1) null, conditional, k = 2 / 20 | 0.52 / 0.60 | 0.19 / 0.005 | 0.52 / 0.04 | 06 |
| Evidence for the null when k grows with n (median log BF₀₁, k = n/4, n = 2000) | 407.6 | 1.9 | 3.7 | 08 |
| Variable selection (n = 500, uncorrelated): real missed of 4 / nulls included of 8 | 0.84 / 0.13 | 0.45 / 0.91 | 0.44 / 0.92 | 09 |
| Estimation: prediction error against OLS, k = 8, n = 30 / 200 | 0.74 / 0.95 | 0.51 / 0.94 | 0.53 / 0.92 | 07 |
| Estimation: 95% interval coverage, real effects, k = 3 / 8 (lowest cell) | 0.957 / 0.940 (0.92) | 0.947 / 0.928 (0.76) | 0.957 / 0.942 (0.82) | 07 |

## 00: what valid correlation matrices force

What happens to one partial correlation and to *R*² as predictors are added,
if the correlation matrix of (*y*, *x*₁ … *x*ₖ) must stay valid? The study
draws 200,000 valid matrices per *k* and keeps a different quantity flat in
each of three versions. Keeping each partial correlation flat makes *R*² grow
with *k*; keeping *R* or *R*² flat makes each partial correlation shrink.

| k | P(\|partial r\| > 0.3) | P(R² > 0.5) | P(\|partial r\| > 0.3) | P(R² > 0.5) | P(\|partial r\| > 0.3) | P(R² > 0.5) |
|---|---|---|---|---|---|---|
| | partial r flat | partial r flat | R flat | R flat | R² flat | R² flat |
| 1 | 0.70 | 0.29 | 0.70 | 0.29 | 0.91 | 0.50 |
| 2 | 0.70 | 0.50 | 0.50 | 0.29 | 0.70 | 0.50 |
| 3 | 0.70 | 0.65 | 0.41 | 0.29 | 0.60 | 0.50 |
| 4 | 0.70 | 0.75 | 0.36 | 0.29 | 0.54 | 0.50 |
| 5 | 0.70 | 0.82 | 0.32 | 0.29 | 0.50 | 0.50 |
| 6 | 0.70 | 0.88 | 0.29 | 0.29 | 0.46 | 0.50 |
| 7 | 0.70 | 0.91 | 0.27 | 0.29 | 0.43 | 0.50 |
| 8 | 0.70 | 0.94 | 0.25 | 0.29 | 0.40 | 0.50 |
| 9 | 0.70 | 0.96 | 0.23 | 0.29 | 0.38 | 0.50 |
| 10 | 0.70 | 0.97 | 0.22 | 0.29 | 0.36 | 0.50 |

- At *k* = 1, "partial r flat" and "R flat" are the same prior (Jeffreys');
  "R² flat" is not flat on *r* even then.
- Figures: `validation/results/figures/00_partial_r_R_R2_by_k.png` (the three
  distributions overlaid per *k*) and `00_three_versions_against_k.png`.
- Because the three versions cannot be reconciled, every later study compares
  three priors.

## 03: how often is the test of one coefficient right? (200 datasets per cell)

How often does each method reach the correct conclusion about one
coefficient? The design has four predictors, uncorrelated or with
correlation 0.5; *x*₁'s slope varies, *x*₃ is null, and *n* runs from 20 to
500. The table gives the share of datasets in which the null predictor is
confirmed, P(BF₀₁ > 3), and in which a small real effect (*x*₁ = 0.1) is
called null, P(BF₀₁ > 3):

| Method | null found, n = 200 unc. / corr. | null found, n = 500 unc. / corr. | 0.1 called null, n = 200 unc. / corr. | 0.1 called null, n = 500 unc. / corr. |
|---|---|---|---|---|
| flat partial r conditional | 0.91 / 0.90 | 0.94 / 0.94 | 0.62 / 0.64 | 0.34 / 0.52 |
| flat R conditional | 0.41 / 0.40 | 0.70 / 0.72 | 0.16 / 0.20 | 0.12 / 0.22 |
| flat R-squared conditional | 0.80 / 0.79 | 0.89 / 0.88 | 0.46 / 0.50 | 0.23 / 0.38 |
| flat partial r joint | 0.82 / 0.83 | 0.91 / 0.90 | 0.47 / 0.51 | 0.25 / 0.40 |
| flat R joint | 0.40 / 0.54 | 0.73 / 0.78 | 0.08 / 0.15 | 0.08 / 0.22 |
| flat R-squared joint | 0.42 / 0.54 | 0.73 / 0.78 | 0.08 / 0.17 | 0.09 / 0.22 |
| JZS default | 0.61 / 0.65 | 0.80 / 0.82 | 0.22 / 0.26 | 0.13 / 0.27 |
| Zellner–Siow | 0.85 / 0.85 | 0.92 / 0.91 | 0.52 / 0.55 | 0.28 / 0.45 |

## 04: are the differences in discrimination or in where the test decides?

Do the methods in study 03 differ in how well they separate nulls from small
effects, or only in where they put the threshold? The **AUC** (= the share of
pairs in which a null dataset gets a larger BF₀₁ than a dataset with a small
effect of 0.1) answers this without any threshold. The table gives each
method's AUC minus that of JZS, with a paired bootstrap 95% CI:

| | n = 100 | n = 200 | n = 500 |
|---|---|---|---|
| conditional tests (all three), uncorrelated | +0.003 (−0.018, 0.023) | +0.003 (−0.010, 0.015) | −0.000 (−0.004, 0.003) |
| conditional tests (all three), correlated (0.5) | +0.047 (0.027, 0.067) | +0.027 (0.014, 0.040) | +0.028 (0.018, 0.039) |
| flat partial r joint, correlated | +0.029 (0.015, 0.042) | +0.018 (0.009, 0.028) | +0.022 (0.014, 0.031) |
| flat R joint, correlated | −0.039 (−0.058, −0.020) | −0.034 (−0.049, −0.019) | −0.028 (−0.036, −0.021) |
| flat R-squared joint, correlated | −0.032 (−0.047, −0.016) | −0.029 (−0.042, −0.017) | −0.027 (−0.035, −0.020) |

- The three conditional tests are expected to have identical AUCs. All three
  rank datasets by the coefficient's partial correlation, so they can differ
  only in where their Bayes factors fall relative to 3.
- Check M4 explains why the joint tests differ. It measures how much the log
  BF for *x*₁ changes per unit of *R*² explained by the other predictors, at a
  fixed *t*-statistic. The change is about zero for the three conditional
  tests (0.002, −0.011 and −0.001), −0.51 for the flat partial r joint test,
  −5.19 for the flat R joint test and −5.10 for the flat R-squared joint
  test; it is −1.93 for JZS and −0.52 for Zellner–Siow.

## 05: are the conditional tests computed correctly?

Are the conditional tests implemented correctly, for every prior and for
blocks of columns? All checks pass:

- C1–C5 for every prior: the test equals the joint test at *k* = 1; the
  implied prior for one coefficient agrees with correlated joint draws; the
  share mixture with *a* = *k*/2 gives a flat partial correlation; and the
  grids agree with `integrate()` up to *n* = 10⁶. C1 has a largest difference
  of 1.8e-06; C2 a KS distance of .0024 (flat R) and .0021 (flat R-squared);
  C4a KS ≤ .0025; C4b 1.8e-05; C5 a largest difference of 3.5e-06 (flat
  partial r), 6.7e-05 (flat R) and 5.0e-05 (flat R-squared).
- C7 shows why η = 1.
- C9–C11 check the block test: C9 KS ≤ .017; C10 1.0e-04; C11 1.6e-14.

**C8: the flat partial r conditional test differs from the Kucharský et al.
(2026) test in one assumption only.** (Numbers from the run of 4 November
2026.) Their exact Bayes factor for a partial correlation with a flat prior
(their eq. 2, with α = 1) is

$$\mathrm{BF}_{10} = \frac{B(\tfrac12, 1 + h)}{B(\tfrac12, 1)}\,
{}_2F_1\!\left(h, h;\ \tfrac32 + h;\ r^2\right),\qquad h = \frac{n - k}{2},$$

where $k$ is our number of predictors (their $k - 1$ controlling variables).
Our implementation of it reproduces their worked example (BF₁₀ = 45.88). Our
own test gives a lower value in 47 of the 48 cases on the grid *r* = 0, 0.1,
0.3, 0.5; *n* = 20, 50, 200, 1000; *k* = 1, 3, 8 (differences in log BF):

| | n = 20 | 50 | 200 | 1000 |
|---|---|---|---|---|
| k = 1, r = 0 | −0.054 | −0.023 | −0.006 | −0.001 |
| k = 3, r = 0 | −0.105 | −0.044 | −0.011 | −0.002 |
| k = 8, r = 0 | −0.262 | −0.098 | −0.024 | −0.005 |
| k = 8, r = 0.5 | −0.128 | −0.029 | −0.005 | −0.001 |

The difference is largest at *r* = 0, shrinks as 1/*n* and grows with *k*. It
has one source. Their model is multivariate normal, so the predictors are
random and the tested predictor's residual SD (given the others) is an
unknown parameter with a 1/σ prior. Ours conditions on *X*, so the effect
size *f* is expressed in units of the observed residual SD. Adding their
uncertainty to ours (*f* in sample units = *f* in population units ×
√(*W*/*n*), with *W* ~ χ²(*n* − *k*)) makes the two agree to within 5e-7 in
all 48 cases.

## 06: does the picture hold across designs? (216 cells)

Do the results of studies 03 and 04 hold beyond one design? The study varies
*k* (2–20), *n* (50–1,000), the predictor correlation (0–0.7) and the
pattern of the other coefficients (sparse or dense), with 200 datasets per
cell.

**Discrimination.** AUC for null against a slope of 0.1, minus that of JZS,
averaged over *n* and the pattern of the other coefficients:

| k | 2 | 5 | 10 | 20 |
|---|---|---|---|---|
| conditional tests, correlation 0 | +0.001 | +0.004 | +0.002 | +0.007 |
| conditional tests, correlation 0.7 | +0.035 | +0.045 | +0.041 | +0.031 |
| flat partial r joint, correlation 0.7 | +0.022 | +0.027 | +0.020 | +0.012 |
| flat R joint, correlation 0.7 | −0.064 | −0.035 | −0.022 | −0.013 |
| flat R-squared joint, correlation 0.7 | −0.059 | −0.033 | −0.019 | −0.012 |

**Calibration at a Bayes factor of 3**, averaged over everything except the
method and *k*:

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

The rate of false effects for a true null (BF₁₀ > 3), averaged in the same
way, is at most .047 over all methods (flat R-squared joint, *k* = 2), and at
most .018 for the conditional tests. In a single design cell the largest is
.105 (flat R-squared joint, *k* = 2, *n* = 200, correlation 0.3, dense).

**What each test says** at a Bayes factor of 3, as a share of datasets, for
the conditional tests and JZS, averaged over *n*, correlation and pattern.
Correct answers are in bold.

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

The three conditional tests rank the datasets identically (AUC for null
against 0.1 of .728 at *k* = 2 and .701 at *k* = 20 for all three; JZS .713 and
.684), so none of them is more accurate; they differ only in where they
decide. At *k* = 20, flat R almost never calls a small effect null (.005)
because it almost never calls anything null: 92% of true nulls stay
undecided, and it finds effects no more often than the other priors (.21 for
a slope of 0.1). Flat partial r reaches a decision on most datasets. It
confirms 88% of true nulls, but it also calls 60% of 0.1 effects and 31% of
0.2 effects null, so with a small effect and little data it leans towards the
null. Neither behaviour is a gain in accuracy: flat R stays silent, and flat
partial r leans.

## 07: how good are the estimates? (24 designs, 100 datasets each)

How close are the estimates to the truth, and how often do the intervals
contain it? The table gives the excess prediction error for a new
observation relative to ordinary least squares (below 1 is better), averaged
over correlation and pattern; the coverage of the 95% intervals; and the mean
bias of real effects:

| | prediction error, k = 3: n = 30 / 60 / 200 | k = 8: n = 30 / 60 / 200 | coverage, real effects, k = 3 / 8 | lowest cell, real effects | mean bias, real effects |
|---|---|---|---|---|---|
| flat partial r | 0.83 / 0.89 / 0.97 | 0.74 / 0.84 / 0.95 | 0.957 / 0.940 | 0.916 | −0.016 (−8%) |
| flat R | 0.72 / 0.85 / 1.02 | 0.51 / 0.73 / 0.94 | 0.947 / 0.928 | 0.760 | −0.054 (−28%) |
| flat R-squared | 0.74 / 0.84 / 0.97 | 0.53 / 0.72 / 0.92 | 0.957 / 0.942 | 0.815 | −0.040 (−21%) |
| JZS default | 0.70 / 0.80 / 0.95 | 0.56 / 0.72 / 0.91 | 0.959 / 0.951 | 0.870 | −0.037 (−19%) |
| Zellner–Siow | 0.86 / 0.91 / 0.97 | 0.77 / 0.85 / 0.95 | 0.952 / 0.939 | 0.910 | −0.014 (−7%) |
| OLS | 1 | 1 | 0.953 / 0.945 | 0.920 | −0.000 (0%) |

The lowest-coverage cell is *k* = 8, *n* = 60, uncorrelated, sparse for
flat R, flat R-squared and JZS; *k* = 8, *n* = 30, uncorrelated, dense for
flat partial r; *k* = 8, *n* = 30, uncorrelated, sparse for Zellner–Siow; and
*k* = 8, *n* = 60, correlation 0.5, sparse for OLS. The mean bias of real
effects is the mean error, with the error relative to the true value in
brackets.

None of the Stan fits had a divergent transition (0 of 2,400 fits for each
of the three priors).

## 08: is each method consistent when the number of predictors grows with *n*?

Does the evidence for the true model keep growing when *k* grows with *n*
(*k* = *n*/4 and *k* = √*n*)? The table gives the median log BF₀₁ for the
true null model against the full model (200 datasets per *n*):

| | n = 200 | n = 2000 |
|---|---|---|
| k = n/4: flat partial r / flat R / flat R-squared / JZS / Zellner–Siow | 39.5 / 1.4 / 2.6 / 15.3 / 41.5 | 407.6 / 1.9 / 3.7 / 149.3 / 410.8 |
| k = √n: flat partial r / flat R / flat R-squared / JZS / Zellner–Siow | 17.7 / 1.6 / 3.2 / 8.5 / 19.2 | 83.0 / 2.6 / 5.1 / 44.6 / 85.1 |

When a small model is true, the median log BF for it against the full model
at *n* = 2,000 and *k* = *n*/4 is 443.8 under flat partial r, 174.1 under
flat R, 174.1 under flat R-squared, 244.9 under JZS and 446.1 under
Zellner–Siow. Every method's evidence keeps rising, so checks G1 and G2 pass
for all three priors. Under flat R and flat R-squared, however, the evidence
for a true null grows slowly: their prior often puts the full model's *R*²
near zero, so an empty model with 500 predictors is hardly penalised.

## 09: how well does each method select variables? (12 candidates, 4 real: 0.3, 0.2, 0.15, 0.1; 100 datasets per cell)

How well does each method pick out the four real predictors among twelve
candidates? At an inclusion-probability cutoff of 0.5 (*n* = 500,
uncorrelated predictors, beta-binomial model prior):

| | real predictors missed (of 4) | null predictors included (of 8) | inclusion of the 0.1 effect | true model found (uniform model prior) |
|---|---|---|---|---|
| flat partial r | 0.84 | 0.13 | 0.45 | 0.37 |
| flat R | 0.45 | 0.91 | 0.64 | 0.36 |
| flat R-squared | 0.44 | 0.92 | 0.65 | 0.35 |
| JZS default | 0.58 | 0.42 | 0.58 | 0.37 |
| Zellner–Siow (1) | 0.89 | 0.12 | 0.43 | 0.35 |
| BIC | 1.01 | 0.05 | 0.34 | 0.27 |

**Trade-off curves.** Sliding the cutoff on the inclusion probabilities
traces, for each method, a curve of real predictors missed against null
predictors included (figure `09_trade_off_curves.png`). The AUC for
separating real from null predictors summarises each curve. The table gives
the difference from flat R under the beta-binomial model prior, with a 95%
bootstrap over datasets:

| | n | JZS | flat partial r | flat R-squared | Zellner–Siow | BIC |
|---|---|---|---|---|---|---|
| uncorrelated | 50 / 100 / 200 / 500 | −0.025* / −0.029* / −0.006 / +0.002 | −0.019* / −0.021* / +0.006 / +0.003 | −0.017* / −0.018* / −0.004 / +0.000* | −0.021* / −0.030* / +0.001 / +0.002 | +0.005 / +0.012 / +0.028* / +0.009* |
| correlation 0.5 | 50 / 100 / 200 / 500 | +0.001 / +0.006 / +0.007* / +0.003* | +0.011 / +0.020* / +0.022* / +0.011* | +0.000 / −0.000 / −0.000 / +0.000 | +0.011 / +0.022* / +0.023* / +0.011* | +0.030* / +0.039* / +0.035* / +0.019* |

\* the 95% interval excludes 0.

- Against the intercept-only model, the evidence for the null grows like
  *a* log *n* (½ log *n* for flat R, log *n* for flat R-squared, *k*/2 log *n*
  for flat partial r; study 02, check B3b). How inclusive each prior is
  therefore follows from its shape *a*; the prior over models determines the
  rest.

## 10: is the conditional test calibrated, and is it coherent with one prior for the whole model?

Two questions: is each conditional test calibrated on data from its own
prior, and does it equal the Bayes factor of its joint model?

**Calibration (checks P1–P2).** A Bayes factor is calibrated if, among
datasets that give a particular BF₁₀, the share of true nulls is
1/(1 + BF₁₀) (with prior odds of 1:1), *when the data come from the prior
that the Bayes factor assumes*. Datasets are drawn from each prior's joint
version, with *x*₁'s coefficient set to 0 in half of them, and from the JZS
prior, where miscalibration should show:

| | χ² p over bins (own prior) | logistic slope (own prior) | JZS data: χ² |
|---|---|---|---|
| flat partial r | 0.40–0.79 | 0.98–1.03 | 503–556 (slopes 0.89–0.95) |
| flat R | 0.52–0.98 | 0.99–1.02 | 233–1198 (slopes 1.16–1.41) |
| flat R-squared | 0.37–0.94 | 0.98–1.02 | 50–265 (slopes 0.95–1.12) |

The P1 ranges cover *k* = 3, 8 and *n* = 50, 200 (χ² over 8–11 bins); all
four cells pass for every prior.

**Coherence (checks P3–P5).** The conditional test puts flat priors on the
other coefficients, so it is not a ratio of marginal likelihoods under the
joint prior. The Savage–Dickey ratio is such a ratio: the joint model's prior
density of *x*₁'s partial correlation at 0, divided by its posterior density
at 0 (from `fit_regression()` draws, using logspline). The table gives the
conditional and joint tests minus the Savage–Dickey ratio, on 20 datasets
(*k* = 3, 8; *n* = 50, 200):

| | conditional: mean absolute / largest | joint: mean absolute / largest |
|---|---|---|
| flat partial r | 0.028 / 0.092 (16 of 20 estimable) | 0.226 / 0.422 |
| flat R | not estimable: infinite prior density at 0 (P7) | not estimable |
| flat R-squared | 0.283 / 0.874 (19 of 20 estimable) | 1.007 / 2.365 |

These numbers cover the datasets where the Savage–Dickey ratio can be
estimated (chain SE < 0.1); on the others the two agree to within 3 Monte
Carlo SE (check P5). For flat R-squared the conditional test lies below the
Savage–Dickey ratio, and further below at larger *k* (mean −0.04 at *k* = 3;
−0.48 and −0.34 at *k* = 8 for *n* = 50 and 200), so check P4 (agreement
within a factor of 2) fails. The reason is that its partial correlations are
linked through *R*², so the joint posterior of *x*₁'s partial correlation
also uses the other coefficients (see "Why check P4 of study 10 fails for
flat R-squared" under Last results).

**What a zero partial correlation means (checks P6–P7).** This part uses
prior draws only: 400,000 per prior and *k*, with independent predictors. The
joint prior is conditioned on |ρ₁| < ε for ε = 0.1, 0.02 and 0.005. The table
gives the density of ρ₁ near 0, P(|ρ₁| < ε) / (2ε), and the median *R*² of
the other predictors under that condition:

| | k | density near 0 (ε = 0.1 / 0.02 / 0.005) | others' median R² |
|---|---|---|---|
| flat partial r | 3 | 0.498 / 0.502 / 0.485 | 0.372 / 0.372 / 0.374 |
| flat partial r | 8 | 0.504 / 0.497 / 0.491 | 0.730 / 0.733 / 0.728 |
| flat R-squared | 3 | 0.741 / 0.767 / 0.781 | 0.185 / 0.169 / 0.180 |
| flat R-squared | 8 | 1.342 / 1.537 / 1.558 | 0.221 / 0.174 / 0.162 |
| flat R | 3 | 1.506 / 2.301 / 2.990 | 0.024 / 0.005 / 0.001 |
| flat R | 8 | 2.229 / 3.827 / 5.204 | 0.054 / 0.012 / 0.003 |

Under flat partial r (density ½, a flat partial correlation) and flat
R-squared, the density settles and the other predictors keep their *R*², so
"*x*₁'s partial correlation is 0" defines a proper restricted model and the
Savage–Dickey ratio applies (check P6). Under flat R the density keeps
growing as ε shrinks, because it is infinite at 0, and the other predictors'
*R*² goes to 0 (check P7). Conditioning on a zero partial correlation for
*x*₁ then drags every other coefficient to zero, so the restricted model
collapses to the intercept-only model. Under flat R, therefore, no joint
model corresponds to testing one coefficient, and its conditional test
stands on its own.

Figures: `validation/results/figures/10_calibration.png`, `10_coherence.png`.

## 11: does the block test of a factor stay unaffected by the other predictors? (q = 2–5 columns, 200 datasets per cell)

When a whole factor of *q* columns is tested, does the evidence for it depend
on how much the other predictors explain? Checks F1–F3 pass. The block
conditional tests give the same Bayes factor however much the other
predictors explain (largest difference 1.4e-12 over 14,400 datasets); the
three rank datasets identically (the same AUC in all 72 cells); and JZS, as
computed by lmBF, is precise enough (largest proportional error 0.032).

The table gives the median log BF₁₀ for a null factor at *n* = 200, as the
other predictors explain *R*² = 0, 0.5 or 0.9 (2 other predictors; 10 others
in brackets):

| | R² of the others 0 | 0.5 | 0.9 |
|---|---|---|---|
| flat partial r conditional | −5.85 (−5.87) | −5.85 (−5.87) | −5.85 (−5.87) |
| flat R conditional | −1.58 (−1.11) | −1.58 (−1.11) | −1.58 (−1.11) |
| flat R-squared conditional | −2.88 (−1.99) | −2.88 (−1.99) | −2.88 (−1.99) |
| flat partial r joint | −4.90 (−3.31) | −6.01 (−4.34) | −8.65 (−7.00) |
| flat R joint | 0.10 (0.00) | −5.70 (−3.40) | −9.43 (−7.08) |
| flat R-squared joint | 0.20 (0.02) | −5.70 (−3.39) | −9.44 (−7.08) |
| JZS (lmBF) | −2.66 (−2.60) | −2.66 (−2.63) | −2.67 (−2.63) |
| Zellner–Siow | −5.29 (−3.45) | −6.44 (−4.49) | −9.07 (−7.14) |

The conditional tests do not change along a row, whereas the joint tests do
(the conditional Lindley paradox; Som et al., 2016). The change is largest
for flat R and flat R-squared, which give no evidence either way for a null
factor when the other predictors explain nothing. JZS, as lmBF computes it,
gives the factor its own *g*, so it barely changes either. True null factors
are confirmed (BF₀₁ > 3), averaged over *n* and the other predictors, at the
following rates from *q* = 2 to *q* = 5: flat partial r conditional .96 to
.99, flat R conditional .51 to .63, and flat R-squared conditional .79 to .87.

Figures: `validation/results/figures/11_auc.png`, `11_calibration.png`.

## Not yet done

The list of open work for the whole project is in
[`docs/0_overview.md`](0_overview.md), section "Open".

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
