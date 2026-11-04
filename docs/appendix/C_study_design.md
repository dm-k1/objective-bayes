# Appendix C. The design of the validation studies (ADEMP)

**Draft for review.** One description per study of the redesigned validation
(see [`REDESIGN.md`](../working/REDESIGN.md) for the structure and the decisions), in the
ADEMP form of Morris, White and Crowther (2019) as adapted for psychology by
Siepe et al. (2024): **A**ims, **D**ata-generating mechanisms, **E**stimands,
**M**ethods and **P**erformance measures. Every number of simulated datasets
(*n*<sub>sim</sub>) is derived from a target Monte Carlo standard error (MCSE), fixed
here before anything is run. The package under test is flatcor 1.0.0
([dm-k1/flatcor](https://github.com/dm-k1/flatcor), pinned in `renv.lock`).

Notation: *n* is the sample size of one simulated dataset, *k* the number of
predictors, *n*<sub>sim</sub> the number of simulated datasets (or prior draws) per
design cell.

---

## How *n*<sub>sim</sub> is chosen

For a performance measure that is a proportion *p* (a rate of decisions,
interval coverage), the Monte Carlo standard error is

$$\mathrm{MCSE} = \sqrt{\frac{p(1-p)}{n_\text{sim}}}, \qquad\text{so}\qquad
n_\text{sim} = \frac{p(1-p)}{\mathrm{MCSE}^2},$$

taken at the worst case *p* = .5 where *p* is not known in advance (Morris et
al., 2019, Table 6). For a mean (a bias, a mean log difference), MCSE =
SD / √*n*<sub>sim</sub>, with the SD taken from the run of 4 November 2026 where one
exists.

| Target | MCSE | *n*<sub>sim</sub> (worst case) |
|---|---|---|
| A proportion reported per design cell | .02 | 625 |
| A proportion reported pooled over cells | .01 | 2,500 in total (met by pooling cells at 625 each) |
| Interval coverage near .95, per cell | .01 | 475 |
| A mean log Bayes factor difference | 0.05 | (SD / 0.05)² |
| A Kolmogorov–Smirnov check at the .01 level | critical distance ≤ .01 | 1.63² / .01² ≈ 26,600 (effective draws) |

Checks then compare an estimate with its target **relative to its MCSE**
(for example "within 2 MCSE of zero") instead of against a fixed number.

All results are reported with their MCSE, computed by `rsimsum` where it
applies.

---

## Shared data-generating mechanisms

Each dataset has an ID (data set, design cell, replication) and its own
L'Ecuyer-CMRG random stream derived from that ID, so the same ID always gives
the same data, on any machine and in any order. Datasets are regenerated from
their seed, not stored.

| Set | How the data arise | Grid | Used by |
|---|---|---|---|
| **A** | Fixed true effects. Predictors are multivariate normal with unit variances and equal correlation ρ; noise SD 1; slopes in units of the noise SD. | *k* ∈ {2, 5, 10, 20}; *n* ∈ {50, 200, 1,000}; ρ ∈ {0, .3, .7}; other slopes sparse (one at 0.3, the rest 0) or dense (all 0.1); slope of *x*₁ ∈ {0, 0.1, 0.2, 0.3}. That is 72 designs, each at 4 slopes. | 7, 9 (a subset), 8 (extended with factors) |
| **B** | Drawn from each prior: the coefficients and the noise SD from the prior (proper stand-ins for the intercept and noise), then the data. | *k* ∈ {3, 8}; *n* ∈ {50, 200} (and *n* = 60 at *k* = 8 for calibration of the posterior); one cell per prior. | 3, 4, 6 |
| **C** | Variable selection: 12 candidates, of which 4 are real (0.3, 0.2, 0.15, 0.1). | *n* ∈ {50, 100, 200, 500}; ρ ∈ {0, .5}. | 10 |
| **D** | *k* growing with *n*: *k* = *n*/4 and *k* = √*n*. Under the null only *R*² is drawn (it is Beta(*k*/2, (*n* − *k* − 1)/2)); with a small true model (3 slopes of 0.3) full datasets are drawn. | *n* ∈ {50, 100, 200, 500, 1,000, 2,000}. | 5 |

---

## The motivating simulation (outside the validation)

`R/docs/motivating_simulation.R`. Not a study: it shows why there are three
priors. Kept as it is (200,000 matrices per *k*), with its checks F1–F3 now
judged against the KS critical value for its *n*<sub>sim</sub>.

---

## 1. Is the prior sampled correctly?

- **Aim.** The Stan model reproduces each prior.
- **Data.** None; the model samples the prior only (`prior_only = 1`).
- **Estimands.** The distribution of *R*² (and, for flat partial r, of each
  partial correlation) under each prior, at *k* ∈ {1, 3, 8}, with uncorrelated
  and correlated predictors.
- **Methods.** flatcor's Stan model for each of the three priors.
- **Performance.** The KS distance between the draws and the exact
  distribution.
- ***n*<sub>sim</sub>.** 26,600 effective draws per setting (KS critical distance at the
  .01 level ≤ .01); the number of iterations is set from the effective sample
  size, not the raw draw count.

## 2. Are the Bayes factors computed correctly?

- **Aim.** Every Bayes factor flatcor reports equals the integral it stands
  for.
- **Data.** A fixed grid of summary statistics (*R*², *n*, *k*; the partial
  correlation and the number of columns tested), including extremes (*R*²
  near 1, *n* up to 10⁶), plus a few simulated datasets for the
  data-based routes.
- **Estimands.** The log Bayes factor of each route.
- **Methods.**
  - Closed forms and grids in flatcor, against three independent routes:
    adaptive integration (`integrate()`); a high-precision ground truth with
    `Rmpfr` (200-bit arithmetic, replacing the floating-point integral as the
    reference); and, for the joint model, `bridgesampling` on the Stan fit.
  - Invariance under rescaling and recoding.
- **Performance.** The absolute log difference from the ground truth.
- ***n*<sub>sim</sub>.** Deterministic except the brute-force and bridge-sampling
  routes: those run until their own MCSE on the log scale is below 0.01 (a
  stated precision for an estimator, not a stopping rule on the result).

## 3. Is the posterior calibrated?

- **Aim.** The Stan posterior is calibrated for each prior (simulation-based
  calibration; Talts et al., 2018), with the `SBC` package.
- **Data.** Set B at *k* = 8, *n* = 60, and *k* = 3, *n* = 50.
- **Estimands.** The rank of the true *R*², of one coefficient and of the noise
  SD among the posterior draws.
- **Methods.** flatcor's Stan model for each prior; the coefficient-coordinates
  model as a contrast, expected to fail.
- **Performance.** `SBC`'s ECDF-difference test and rank histograms.
- ***n*<sub>sim</sub>.** 1,000 fits per prior and setting, thinned to 100 posterior draws
  each. That is larger than the 400 of the current run, to give the ECDF test
  power against miscalibration of about 5% in the rank distribution.
  *To check against the guidance in the `SBC` documentation before fixing.*

## 4. Is the Bayes factor calibrated?

- **Aim.** On data from each prior, the share of true nulls among datasets with
  a given BF₁₀ equals 1 / (1 + BF₁₀) (Schad et al., 2023); on data from another
  prior (JZS), it does not.
- **Data.** Set B: half the datasets from the null, half from the prior's
  alternative; generating priors: the three, and JZS.
- **Estimands.** The observed share of nulls per bin of log BF₁₀.
- **Methods.** The conditional test under each prior.
- **Performance.** Observed against predicted share per bin, the χ² statistic
  over bins and the logistic calibration slope, with MCSE.
- ***n*<sub>sim</sub>.** 20,000 datasets per cell, as now: with bins holding at least 625
  datasets each, every bin's share has MCSE ≤ .02. Cheap (closed forms).

## 5. Does each prior meet the default-prior criteria?

- **Aim.** The criteria of Bayarri et al. (2012): predictive matching,
  invariance, information consistency, intrinsic consistency, and
  model-selection consistency for fixed and for growing *k*.
- **Data.** Mostly none (the criteria are statements about the Bayes factor as
  a function of *R*², *n* and *k*); simulated *R*² under the null for the rates;
  set D for growing *k*.
- **Estimands.**
  - Fixed *k*: the slope of the median log BF₀₁ against log *n* (expected *a* for
    the intercept-only model, *q*/2 between models).
  - Growing *k*: the slope against log *n* of the median log Bayes factor for
    the truth.
- **Methods.** The three priors; JZS and Zellner–Siow as references.
- **Performance.** The slopes with bootstrap MCSE.
- ***n*<sub>sim</sub>.**
  - Fixed *k*: 2,000 *R*² draws per (*n*, *q*) setting, which gives a slope MCSE
    of about 0.02. The current 200 gave slopes that moved by up to 0.1 between
    runs.
  - Growing *k* under the null: 2,000 *R*² draws per setting (cheap).
  - Growing *k* with a small true model: 400 datasets per setting. Bounded by
    the cost of the *k* = 500 regressions at *n* = 2,000.

## 6. Is the test of one coefficient coherent with one model?

- **Aim.** For each prior, does the conditional test equal the Bayes factor of
  the full model against its restriction by conditioning on the tested partial
  correlation being zero (the Savage–Dickey ratio), and is that restriction a
  proper model?
- **Data.** Set B, with *x*₁'s *t*-value set to 0, 1, 1.5, 2 and 2.5; for the
  restriction, prior draws only.
- **Estimands.**
  - The mean log difference between the conditional test and Savage–Dickey,
    per prior and *k*.
  - The density of the partial correlation near zero, and the other
    predictors' *R*² given it.
- **Methods.** flatcor's conditional test; Savage–Dickey from the Stan
  posterior (logspline, with a kernel estimate as a cross-check). Flat R is
  covered by the restriction check only: its prior density at zero is infinite.
- **Performance.**
  - The mean log difference with its MCSE. The check passes when it is within
    2 MCSE of zero; this replaces the current P4, which used the largest
    difference over about 10 datasets.
  - The restriction: the distribution of the others' *R*² against its
    predicted form (for flat partial r, Beta((*k*−1)/2, 3/2)), by KS distance.
- ***n*<sub>sim</sub>.**
  - Coherence: the SD of the difference was about 0.30 in the last run (flat
    R-squared, *k* = 8). For MCSE ≤ 0.05 that needs (0.30/0.05)² = 36
    estimable datasets per prior and *k*. Planned: 50 per prior, *k* and *n*
    (200 Stan fits per prior), since some fits are not estimable.
  - Restriction: 400,000 prior draws per prior and *k*, as now (cheap).

## 7. How does the test of one coefficient behave? (worked in full)

- **Aims.**
  1. How often each test gives each answer about one coefficient: supports the
     null (BF₀₁ > 3), supports an effect (BF₁₀ > 3), or undecided.
  2. How that depends on *k*.
  3. Whether differences between priors are differences in accuracy (the
     ranking of datasets) or only in where they decide.
- **Data.** Set A, all 72 designs, at slopes of *x*₁ 0, 0.1, 0.2 and 0.3. This
  replaces current studies 03, 04 and 06. Slope 0.3 is kept to show where all
  tests agree.
- **Estimands, per method and design cell:**
  - at slope 0: P(supports the null) and P(supports an effect, a false
    effect);
  - at slope > 0: P(supports the null) and P(supports an effect);
  - the AUC for separating slope 0 from slope 0.1 (and from 0.2), and its
    difference from JZS;
  - pooled over cells at each *k*: the same rates.
- **Methods.** The conditional test under each prior (3); the joint test under
  each prior (3); JZS (BayesFactor, scale 0.354) and Zellner–Siow (scale 1). All
  eight are applied to every dataset, so every comparison is paired.
- **Performance.**
  - Rates with MCSE √(*p*(1−*p*)/*n*<sub>sim</sub>).
  - The AUC difference from JZS with a paired bootstrap MCSE over datasets.
  - The check of the main claim: at every *k*, the three conditional tests'
    AUCs agree within 2 MCSE, and the difference in the null-confirmation rate
    between *k* = 2 and *k* = 20 is within 2 MCSE of zero for flat partial r.
- ***n*<sub>sim</sub>.** 625 datasets per design cell and slope, so every per-cell rate
  has MCSE ≤ .02. Pooled over the 18 cells at each *k*, that gives 11,250
  datasets and MCSE ≤ .005. The AUC per cell then has MCSE of about .018 at
  worst, and its paired difference from JZS less. In total 72 × 4 × 625 =
  180,000 datasets, each analysed by 8 methods. The cost is dominated by
  BayesFactor's JZS and Zellner–Siow (about 360,000 calls; a few core-hours on
  a cluster).

## 8. How does it behave for factors?

- **Aim.** The block (several-column) conditional test, which applies to
  factors, against the joint test and JZS. Its evidence should not depend on
  what the other predictors explain.
- **Data.** Set A extended: a factor of *q* ∈ {2, 3, 4, 5} columns with effect
  (Cohen's *f*²) 0, .02 or .15; 2 or 10 other predictors explaining *R*² of 0,
  .5 or .9; *n* ∈ {50, 200, 1,000}.
- **Estimands.** The same three answers as in study 7, and the AUC null vs
  small.
- **Methods.** The block conditional test (3 priors); the joint test (3); JZS
  (BayesFactor `lmBF`); Zellner–Siow.
- **Performance.** As in study 7, plus the slope of the median log BF against
  the others' *R*² (expected zero for the conditional tests).
- ***n*<sub>sim</sub>.** 625 datasets per cell (rates MCSE ≤ .02), against 200 now.

## 9. How good are the estimates?

- **Aim.** Bias, coverage and prediction error of the Stan posterior under each
  prior, against ordinary least squares and the BayesFactor posteriors (JZS,
  Zellner–Siow).
- **Data.** A subset of set A, chosen so that estimation and testing are judged
  on **the same datasets**: *k* ∈ {2, 10}, *n* ∈ {50, 200}, ρ ∈ {0, .7},
  sparse and dense, slope of *x*₁ 0.2. That is 16 cells. The current study used
  its own grid (*k* ∈ {3, 8}, *n* ∈ {30, 60, 200}).
- **Estimands, per method and cell:**
  - the bias of *x*₁'s estimate, as a share of its true value;
  - the coverage of *x*₁'s 95% interval;
  - the excess prediction error relative to OLS.
- **Methods.** flatcor's Stan model (3 priors); BayesFactor posterior draws
  (JZS, Zellner–Siow); OLS.
- **Performance.** Means and coverage with MCSE. Coverage is judged on one
  coefficient per dataset (*x*₁), so the replications are independent; the
  current study pooled correlated coefficients within datasets.
- ***n*<sub>sim</sub>.** 475 datasets per cell (coverage near .95 with MCSE ≤ .01), so
  16 × 475 × 3 ≈ 23,000 Stan fits. That is a cluster job: each fit takes
  seconds.

## 10. How does it select variables?

- **Aim.** Variable selection over all 4,096 subsets of 12 candidates, with a
  uniform and a beta-binomial prior over models.
- **Data.** Set C.
- **Estimands, per method and cell:**
  - the mean number of real predictors missed (of 4) and of null predictors
    included (of 8), at an inclusion probability of .5;
  - the mean inclusion probability of each true slope;
  - the AUC of inclusion probabilities, real against null predictors.
- **Methods.** The three priors (whole-model Bayes factors); JZS, Zellner–Siow;
  BIC; and, from `BAS`, hyper-*g*/*n* (and the robust and intrinsic priors if
  `bas.lm` offers them for linear models; to check).
- **Performance.** Means with MCSE; the AUC difference from flat R with a
  paired bootstrap MCSE.
- ***n*<sub>sim</sub>.** 625 datasets per cell, against 100 now: per-cell shares then have
  MCSE ≤ .02, and mean counts have MCSE ≤ about .05.

---

## Summary

| Study | *n*<sub>sim</sub> per cell (now) | *n*<sub>sim</sub> per cell (planned) | Main cost |
|---|---|---|---|
| 1 | 4,000–40,000 draws | 26,600 effective draws | Stan, prior only |
| 2 | deterministic | deterministic; brute force and bridge sampling to MCSE 0.01 | `Rmpfr` |
| 3 | 400 fits | 1,000 fits | Stan |
| 4 | 20,000 | 20,000 | none |
| 5 | 200 | 2,000 *R*² draws; 400 datasets (growing *k*, true model) | the *k* = 500 regressions |
| 6 | 5 per *t*-value (20 per prior) | 50 per prior, *k* and *n* | Stan |
| 7 | 200 | 625 | BayesFactor |
| 8 | 200 | 625 | BayesFactor |
| 9 | 100 | 475 | Stan |
| 10 | 100 | 625 | enumeration of 4,096 models |

## Open choices for review

1. The targets: per-cell proportions at MCSE .02, coverage at .01, coherence at
   0.05. Stricter targets multiply the cost by the square of the change (half
   the MCSE means four times the datasets).
2. Whether study 9 should move onto set A's grid (as proposed, to share
   datasets with study 7) or keep its own smaller *n* (30, 60), where shrinkage
   matters most.
3. Whether slope 0.3 in set A is worth its cost (25% of study 7's datasets) to
   show where the tests agree.
4. The SBC settings (study 3) against the `SBC` package's guidance.

## References

- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012).
  Criteria for Bayesian model choice with application to variable selection.
  *Annals of Statistics*, *40*(3), 1550–1577.
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
  studies to evaluate statistical methods. *Statistics in Medicine*, *38*(11),
  2074–2102. https://doi.org/10.1002/sim.8086
- Schad, D. J., Nicenboim, B., Bürkner, P.-C., Betancourt, M., & Vasishth, S.
  (2023). Workflow techniques for the robust use of Bayes factors.
  *Psychological Methods*, *28*(6), 1404–1426.
- Siepe, B. S., Bartoš, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W., &
  Pawel, S. (2024). Simulation studies for methodological research in
  psychology: A standardized template for planning, preregistration, and
  reporting. *Psychological Methods*. https://doi.org/10.1037/met0000695
- Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018).
  Validating Bayesian inference algorithms with simulation-based calibration.
  arXiv:1804.06788.
