# Objective Bayes for linear regression: the possible-range prior

> Define what's possible, not what you expect.

Bayes factors and estimation for linear regression with a prior that only
encodes the **possible range** of effects, not an expectation about their size:

- the explained share of variance (R²) gets a Beta(½, 1) prior, the same for
  any number of predictors. With one predictor this is exactly a flat prior on
  the correlation from −1 to 1 (Jeffreys' default);
- every way of dividing that share over the predictors (every direction) is
  equally plausible, relative to how the predictors correlate;
- flat prior on the intercept, Jeffreys prior on the noise level (the same
  conventions as the JZS / JASP default).

Estimation runs through **brms**, so all the usual brms output works. Bayes
factors between any two linear models (nested or not) come from an exact
closed form that depends on the data only through R², n and k. It is checked
against exact integration over R² from 0.01 to 0.99, n from 10 to 10,000 and
k from 1 to 30: the largest error is 4 × 10⁻⁵ on the log scale
(`validation/01_checks.R`).

## Quick start

Set the working directory to the project folder (all scripts use paths
relative to it), then:

```r
source("R/possible_range.R")   # light: brms is only loaded when you fit
library(brms)                  # for summary(), plot(), ... on fitted models

fit  <- brm_possible(y ~ x1 + x2 + x3, my_data)    # estimation (brms)
summary(fit)

possible_bayes_factor_per_term(y ~ x1 + x2 + x3, my_data)   # BF per term
possible_bayes_factor(fit, brm_possible(y ~ x2 + x3, my_data))  # any two models
possible_bayes_factor(lm(y ~ x1 + x2, my_data), lm(y ~ x3, my_data))  # no Stan needed
```

Requirements: R ≥ 4.3, brms with the cmdstanr backend (CmdStan installed).
Optional: BayesFactor (JZS comparison), dplyr/purrr/tidyr/ggplot2 (tutorial
and validation scripts).

## Folder structure

| Path | What it is |
|---|---|
| `R/possible_range.R` | The method: `brm_possible()` (estimation), `possible_bayes_factor()` and `possible_bayes_factor_per_term()` (closed-form Bayes factors), `support_for_each_effect_size()` and `support_verdicts()` (reporting) |
| `vignettes/possible_range_vignette.R` | How to use it on a dataset, start to finish |
| `tutorial/tutorial_step_by_step.R` | Why the method works, step by step, without a math background |
| `tutorial/posterior_by_hand.R` | The posterior and Bayes factor computed by hand on a grid (one predictor), with a plot per step |
| `docs/derivation.md` | The math: the prior density on the coefficients and the closed-form Bayes factor, step by step |
| `docs/relation_to_r2d2.md` | Why this is not the R2-D2 prior, with a one-predictor comparison table |
| `validation/` | Checks of the closed form and simulation studies; see `validation/README.md` |
| `archive/` | Superseded or abandoned code, kept for history |

## Scope

Supported: linear regression with a normal outcome, an intercept and fixed
effects. That includes continuous predictors, factors (tested as one block),
interactions, polynomials and spline columns.

Not (yet) supported, and `brm_possible()` stops for some of these:

- group-level (random) effects: the R² budget would need to include the
  group-level variance (cf. the R2D2M2 prior of Aguilar & Bürkner, 2023);
- other families (logistic, Poisson, …) and distributional models;
- models without an intercept.

Points of attention:

- With several predictors at small n, the prior pulls estimates toward zero
  (about 20% for the vignette's x1 at n = 60). This follows from the fixed R²
  budget: each predictor's prior share is peaked at zero.
- The brms parameterization has a sharp (integrable) peak at all-zero
  coefficients. Occasional divergent transitions can result; use
  `control = list(adapt_delta = 0.99)`. If they persist, see
  `archive/stan_unit_vector_model.R` for a parameterization without that peak.
- `possible_bayes_factor_per_term()` removes one term at a time. Testing a
  main effect while its interaction stays in the model follows the usual
  marginality caveats.
- Do not use `brms::hypothesis()` for Bayes factors with this prior: brms
  cannot sample from a custom prior, so Savage–Dickey does not work. Use
  `possible_bayes_factor()`.

## Related work

- Jeffreys' uniform prior on the correlation; Ly, Verhagen & Wagenmakers on
  Jeffreys' default tests.
- Kucharský, Wagenmakers, van den Bergh & Ly (2026), *Test*: exact Bayes factor
  for one partial correlation with a (stretched) beta prior; α = 1 is the
  per-coefficient flat prior compared in `validation/02`.
- rstanarm `stan_lm`: R² prior with a uniform direction, but R² ~ Beta(K/2, η),
  which drifts with the number of predictors.
- R2-D2 (Zhang et al., 2022) and R2D2M2 (Aguilar & Bürkner, 2023): fixed prior
  on R², split by a Dirichlet; compared in `validation/03`.
- Liang et al. (2008), hyper-g; Maruyama & George (2011), generalized g-prior.

## History

- **v1 (March 2026)**, separate project *Computational Bayesian prior
  selection*: an automated prior predictive calibration in brms. One normal
  scale shared by all predictors, matched so that the predictors together can
  explain at most the outcome's variance; the scale shrinks with the number of
  predictors (≈ 1/√K).
- **v2 (October 2026, this project)**: the exact version. A fixed Beta(½, 1)
  prior on R² with a uniform direction, which reduces to Jeffreys' flat
  correlation prior with one predictor, a closed-form Bayes factor, and brms
  estimation.

## Author

Diemithry M. Kloppenburg, University of Amsterdam.

## License

Not chosen yet.
