# Objective Bayes for linear regression: the flat-R prior

> Define what's possible, not what you expect.

**In one sentence:** the total strength of all predictors together (the
multiple correlation $R$) is flat from 0 to 1, for any number of predictors,
and every way of dividing it over the predictors is equally plausible.

As a density on the coefficients:

```
P(beta) = P(total) x P(this combination | total) x unit conversion
          flat on R   every split equally likely   noise units -> data units
```

- With one predictor this is exactly Jeffreys' flat prior on the correlation.
- With more predictors the prior's claim about explained variance stays the
  same ($R^2 \sim$ Beta(½, 1)), instead of rising with every predictor as
  per-coefficient priors do. Each predictor's share then peaks at zero: a
  built-in correction for testing many predictors.
- Intercept and noise get the usual neutral priors (flat, and $1/\sigma$), as in
  the JZS / JASP default.
- It is a mixture of Zellner g-priors, the family of the JZS default: the member
  whose implied $R^2$ stays the same for every number of predictors
  ([`docs/2_how_it_calculates.md`](docs/2_how_it_calculates.md), section 8).

It does carry expectations, stated rather than hidden: flat on $R$ means
$E[R^2] = \tfrac13$, and the per-predictor shrinkage follows from the shared
budget ([`docs/1_from_behaviour.md`](docs/1_from_behaviour.md), last section).

**What you get:** estimation with Stan, sampled in the prior's own coordinates
($R$ and a direction), with the usual tools (`loo()`, `pp_check()`, bayesplot);
and Bayes factors between any two linear models, nested or not, from an exact
closed form that needs only $R^2$, $n$ and $k$.

## Quick start

Set the working directory to the project folder, then:

```r
source("R/flat_r.R")

fit <- fit_flat_r_model(y ~ x1 + x2 + x3, my_data)   # estimation (Stan)
fit                                                  # estimates and diagnostics
loo(fit); pp_check(fit)

bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, my_data)              # BF per term
bayes_factor_flat_r(fit, lm(y ~ x2 + x3, my_data))                   # any two models
bayes_factor_flat_r(lm(y ~ x1 + x2, my_data), lm(y ~ x3, my_data))   # no Stan needed
```

Requirements: R ≥ 4.3 and MASS. For estimation: cmdstanr with CmdStan (or
rstan), posterior, loo, bayesplot. Optional: BayesFactor (JZS comparison),
ggplot2, dplyr, purrr, tidyr (vignette, figures and validation).

## What's where

| Path | What it is |
|---|---|
| `R/flat_r.R` | The method: `fit_flat_r_model()` (estimation), `bayes_factor_flat_r()` and `bayes_factor_flat_r_per_term()` (Bayes factors), `draw_coefficients_from_flat_r_prior()` (the prior as a recipe) |
| `R/flat_r.stan` | The Stan model: samples $R$ and a direction, computes the coefficients |
| `R/support_curve.R` | Optional reporting next to a null Bayes factor: which effect sizes the data can't tell apart from zero |
| `vignette.R` | Use it on a dataset, start to finish |
| `docs/1_from_behaviour.md` | Where the prior comes from: from what effects can do to the one sentence, with simulated figures |
| `docs/2_how_it_calculates.md` | How it calculates: the density factor by factor, and the Bayes factor |
| `docs/relation_to_others.md` | Which priors drift with the number of predictors; why this is not R2-D2 |
| `docs/figures/make_figures.R` | Makes every figure in the docs (each a simulation against a formula) |
| `validation/` | Does the Stan model get prior and posterior right, is the Bayes factor right, how does it behave; see `validation/README.md` |
| `archive/` | Superseded code (by-hand tutorials, the brms version of the prior, earlier approaches), kept for history |

## Scope

Supported: linear regression with a normal outcome, an intercept and fixed
effects, including continuous predictors, factors (tested as one block),
interactions, polynomials and spline columns.

Not (yet) supported, and `fit_flat_r_model()` stops for some of these:
group-level (random) effects (the $R^2$ budget would need the group-level
variance, cf. R2D2M2), other families, distributional models, models without
an intercept.

Points of attention:

- With several predictors at small $n$ the prior pulls estimates toward zero;
  the vignette prints the size of the pull next to the ordinary estimates.
- At the start of warmup Stan may print a few "Informational Message ... proposal
  is about to be rejected" lines (the sampler tried R at its boundary). A few are
  harmless; `fit` reports divergent transitions, which is what to check.
- Why not brms: brms samples the coefficients directly, and for $k \ge 2$ the
  prior has a sharp peak at all-zero coefficients that such a sampler misses,
  without warnings. The Stan model here samples $R$ and a direction instead
  (`docs/2_how_it_calculates.md`, section 6).
- `bayes_factor_flat_r_per_term()` removes one term at a time; testing a main
  effect while its interaction stays in the model follows the usual
  marginality caveats.

## History

- **v1 (March 2026)**, separate project *Computational Bayesian prior
  selection*: an automated prior predictive calibration in brms, with one normal
  scale shared by all predictors so that together they explain at most the
  outcome's variance (scale ≈ 1/√k).
- **v2 (October 2026, this project)**: the exact version, with a closed-form
  Bayes factor and Stan estimation in the prior's own coordinates.

## Author

Diemithry M. Kloppenburg, University of Amsterdam.

## License

Not chosen yet.
