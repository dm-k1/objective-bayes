# flatcor: objective Bayes for linear regression with flat correlations

> Define what's possible, not what you expect.

**What it shows.** For testing one regression coefficient with a default
Bayes factor, it matters what the default keeps flat as predictors are added.
Keep each partial correlation flat (**flat partial r**) and the test means the
same at every number of predictors. Keep the total flat (**flat R**, **flat
R²**) and the test drifts with the number of predictors, without becoming
more accurate: at 20 predictors all three rank datasets identically, but
flat R leaves 92% of true nulls undecided, where flat partial r confirms 88%.
Only flat partial r's test is also the Bayes factor of one joint model.

![What each test says about one coefficient as predictors are added](validation/results/figures/06_outcomes.png)

**The evidence in one page:** [`docs/findings.md`](docs/findings.md). All
three are built and validated side by side against the criteria from the
literature (131 of 132 checks pass; [report](validation/report.html),
[details](docs/5_validation.md)); the comparison and the recommendation of
flat partial r, with its costs, are in [`docs/4_comparison.md`](docs/4_comparison.md).

## The three priors

**Where it starts** ([`validation/00_initial_simulation.R`](validation/00_initial_simulation.R),
the first result): draw valid correlation matrices of $(y, x_1 \ldots x_k)$.
Keep each partial correlation flat and $R^2$ grows with $k$; keep $R$ or
$R^2$ flat and each partial correlation shrinks toward zero with $k$. The
three options that leave:

| | What is flat | Implied $R^2$ | One coefficient's prior |
|---|---|---|---|
| **flat partial r** (`prior = "flat_partial_r"`) | every valid correlation matrix of $(y, x_1 \ldots x_k)$ | Beta($k/2$, 1), grows with $k$ | a flat partial correlation, for every $k$ |
| **flat R** (`prior = "flat_r"`) | the multiple correlation $R$ | Beta(½, 1), for every $k$ | its share of the flat total; narrows as $k$ grows |
| **flat R²** (`prior = "flat_r_squared"`) | $R^2$ | Beta(1, 1), for every $k$ | its share of the flat total; narrows as $k$ grows |

All three are $R^2$ ~ Beta($a$, 1) with $a$ = $k/2$, ½ or 1. There is **no
default**: every function asks for `prior =`.

All three divide the total over the predictors with no preference. As a
density on the coefficients:

```
P(beta) = P(total) x P(this combination | total) x unit conversion
          R^2 ~ Beta(a, 1)     every split equally likely   noise units -> data units
```

- All three are mixtures of Zellner g-priors (flat R² from two predictors
  on), the family of the JZS / JASP default; intercept and noise get the same
  neutral priors as there (flat, and $1/\sigma$).
- Flat partial r is rstanarm `stan_lm`'s prior with η = 1, the one shape for
  which every coefficient's partial correlation is flat.

**What you get:** estimation with Stan (with `loo()`, `pp_check()`, bayesplot),
Bayes factors between any two linear models from a formula that needs only
$R^2$, $n$ and $k$, and a test per coefficient.

**Then [`docs/0_overview.md`](docs/0_overview.md)**: the map of the
project, what each claim rests on and how far it is shown, and what to watch
out for when using it.

## Quick start

Install once, from the project folder (or from GitHub once it is public):

```r
devtools::install()                    # or: remotes::install_github("dm-k1/objective-bayes")
```

Then:

```r
library(flatcor)

## every function asks for the prior: "flat_partial_r", "flat_r" or "flat_r_squared"
fit <- fit_regression(y ~ x1 + x2 + x3, my_data, prior = "flat_partial_r")   # estimation (Stan)
fit                                                                         # estimates and diagnostics
loo(fit); pp_check(fit)

bayes_factor_per_term(y ~ x1 + x2 + x3, my_data, prior = "flat_partial_r")   # BF per term (conditional test)
bayes_factor_models(fit, lm(y ~ x2 + x3, my_data), prior = "flat_partial_r") # any two models
bayes_factor_models(lm(y ~ x1 + x2, my_data), lm(y ~ x3, my_data),           # no Stan needed
                    prior = "flat_r")

## the same with another prior
bayes_factor_per_term(y ~ x1 + x2 + x3, my_data, prior = "flat_r_squared")
```

Requirements: R ≥ 4.3; MASS, posterior, loo and bayesplot are installed with
the package. For estimation also cmdstanr with CmdStan (or rstan). Optional:
BayesFactor (JZS comparison), ggplot2, dplyr, purrr, tidyr, patchwork,
logspline (example, figures and validation). Run the tests with `devtools::test()`. The exact package
versions are in `renv.lock`, and `validation/Dockerfile` builds the whole environment;
running the validation anywhere: `docs/5_validation.md`.

## What's where

| Path | What it is |
|---|---|
| `R/fit.R` | `fit_regression()`, with `print()`, `loo()` and `pp_check()` |
| `R/priors.R` | The three priors: their names and the shape $a$ of $R^2$ ~ Beta($a$, 1) |
| `R/bayes_factor.R` | `bayes_factor_models()`, `bayes_factor_per_term()` and `log_bayes_factor_vs_intercept_only()`, the formula behind them |
| `R/conditional.R` | `log_bayes_factor_conditional()`: the test for one coefficient or one term |
| `R/prior_draws.R` | `draw_coefficients_from_prior()`: the prior as a recipe |
| `R/support_curve.R` | Optional reporting next to a null Bayes factor: which effect sizes the data can't tell apart from zero |
| `inst/stan/flatcor.stan` | The Stan model: samples $R$ (with $R^2$ ~ Beta($a$, 1)) and a direction, computes the coefficients |
| `tests/` | Fast checks (`devtools::test()`); the full checks are in `validation/` |
| `DESCRIPTION`, `NAMESPACE`, `man/` | Package files (`man/` and `NAMESPACE` are generated by roxygen2 from the comments in `R/`) |
| `docs/example.R` | Use it on a dataset, start to finish |
| `docs/findings.md` | **What the project shows**, in one page: the claim, the three results behind it, the costs, what is new |
| `docs/0_overview.md` | **The map**: the fork from study 00, the three priors, the questions the package answers, the status of each claim, what to watch out for, what is open |
| `docs/1_from_behaviour.md` | Where the prior comes from: from what effects can do to the one sentence, with simulated figures |
| `docs/2_how_it_calculates.md` | How it calculates: the density factor by factor, and the Bayes factor |
| `docs/3_consistency.md` | Proofs and sketches for the consistency criteria |
| `docs/4_comparison.md` | The three priors compared, criterion by criterion, with the JZS default and Zellner–Siow as references |
| `docs/5_validation.md` | The validation report: how to run the studies, the criteria from the literature and where each is checked, every result, references |
| `docs/relation_to_others.md` | Which priors drift with the number of predictors; flat partial r as `stan_lm`'s prior; why this is not R2-D2 |
| `docs/figures/make_figures.R` | Makes every figure in the docs (each a simulation against a formula) |
| `validation/` | Study 00 (the first result: the fork), then for each prior: does the Stan model get prior and posterior right, is the Bayes factor right, how does it behave; `run_all.R` runs everything; see `docs/5_validation.md` |
| `validation/report.qmd`, `report.html` | The validation report: every figure and check, with what each score means (open `report.html` in a browser; rebuild with `quarto render validation/report.qmd`) |
| `renv.lock`, `validation/Dockerfile` | The pinned environment: R packages, and R + packages + CmdStan as a container |

## Scope

Supported: linear regression with a normal outcome, an intercept and fixed
effects, including continuous predictors, factors (tested as one block),
interactions, polynomials and spline columns.

Not (yet) supported, and `fit_regression()` stops for some of these:
group-level (random) effects (the shared $R^2$ total would need the group-level
variance, cf. R2D2M2), other families, distributional models, models without
an intercept.

What to watch out for when using it (shrinkage at small $n$, the conditional
test with many predictors, variable selection, Stan messages):
[`docs/0_overview.md`](docs/0_overview.md).

## History

- **v1 (March 2026)**, separate project *Computational Bayesian prior
  selection*: an automated prior predictive calibration in brms, with one normal
  scale shared by all predictors so that together they explain at most the
  outcome's variance (scale ≈ 1/√k).
- **v2 (October 2026, this project)**: the exact version, with a closed-form
  Bayes factor and Stan estimation in the prior's own coordinates; first
  flat R, then uniform correlations (now flat partial r) as the default.
- **v3 (November 2026)**: study 00 showed that not everything can stay flat,
  so the three options (flat partial r, flat R, flat R²) are built and
  validated as equals, with no default; the results recommend flat partial r
  (`docs/4_comparison.md`, section 8). The superseded scripts (by-hand
  tutorials, the brms version of the prior, earlier approaches) were in
  `archive/` and are in the git history up to commit a2886a7.

## Author

Diemithry M. Kloppenburg, University of Amsterdam.

## License

MIT (`LICENSE.md`).
