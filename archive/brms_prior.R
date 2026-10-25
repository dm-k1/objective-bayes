## ARCHIVED: the flat-R prior as a custom brms prior (a Stan density on the
## coefficients b). Correct as a density, but samplers moving in b miss the
## peak at zero for larger k (validation/01, contrast). Replaced by
## R/flat_r.stan, which samples R and a direction. Stands alone; brms is
## called as brms::...

## The prior as a density on the coefficients, in one sentence:
##   P(b) = P(total) x P(this combination | total) x unit conversion
##   1. P(total)                        flat on R, per unit of f
##                                      (f = length of b in noise units, Cohen's f)
##   2. P(this combination | total)     every split equally likely:
##                                      1 / size of all combinations with this total
##   3. unit conversion                 from noise units to the units of b
## Where each factor comes from: docs/1_from_behaviour.md and docs/2_how_it_calculates.md.
## Checked in brms: archive/check_brms_prior.R.

flat_r_stan_code <- "
  // P(b) = P(total) x P(this combination | total) x unit conversion
  real flat_r_lpdf(vector b, real sigma, matrix predictor_covariance) {
    int k = rows(b);

    // The total: length of b in noise units (Cohen's f)
    real f_squared = quad_form(predictor_covariance, b) / square(sigma);
    real f = sqrt(f_squared);

    // 1. P(total): flat on R, per unit of f
    real log_p_total = -1.5 * log1p(f_squared);

    // 2. P(this combination | total): 1 / size of all combinations with this total
    //    = 1 / surface of a sphere with radius f in k dimensions
    real log_size_of_all_combinations = log(2) + 0.5 * k * log(pi()) - lgamma(0.5 * k)
                                        + (k - 1) * log(f);
    real log_p_combination_given_total = -log_size_of_all_combinations;

    // 3. Unit conversion: from noise units to the units of b
    real log_conversion = 0.5 * log_determinant(predictor_covariance) - k * log(sigma);

    // log of (1 x 2 x 3)
    return log_p_total + log_p_combination_given_total + log_conversion;
  }
"

## Builds the brms prior for a formula: each model uses its own predictors'
## covariance, so the prior adapts to correlated predictors automatically.
build_flat_r_prior <- function(formula, data) {
  predictors <- model.matrix(formula, data)[, -1, drop = FALSE]
  neutral_priors <- brms::set_prior("", class = "Intercept") +   # flat intercept
                    brms::set_prior("", class = "sigma")         # Jeffreys term added below
  if (ncol(predictors) == 0)
    return(list(prior = neutral_priors + brms::set_prior("target += -log(sigma);", check = FALSE),
                stanvars = NULL))
  centered_predictors  <- scale(predictors, scale = FALSE)
  predictor_covariance <- crossprod(centered_predictors) / nrow(centered_predictors)
  list(
    prior = neutral_priors + brms::set_prior(
      "target += flat_r_lpdf(b | sigma, predictor_covariance) - log(sigma);",
      check = FALSE),
    stanvars = brms::stanvar(scode = flat_r_stan_code, block = "functions") +
               brms::stanvar(predictor_covariance, name = "predictor_covariance"))
}

## Fit exactly like brm(): summary(), plot(), pp_check(), bayes_R2(),
## loo() all work as usual. Extra arguments go to brm(). The family may
## be given, but only gaussian with the identity link is supported.
fit_flat_r_model <- function(formula, data, family = gaussian(), ...) {
  if (any(grepl("|", deparse(formula), fixed = TRUE)))
    stop("fit_flat_r_model() supports fixed effects only (no group-level terms like (1 | g)).")
  if (is.character(family)) family <- get(family, mode = "function")
  if (is.function(family))  family <- family()
  if (!identical(family$family, "gaussian") || !identical(family$link, "identity"))
    stop("fit_flat_r_model() supports the gaussian family with the identity link only.")
  prior_parts <- build_flat_r_prior(formula, data)
  brms::brm(formula, data = data, family = gaussian(),
            prior = prior_parts$prior, stanvars = prior_parts$stanvars,
            backend = getOption("brms.backend", "cmdstanr"), ...)
}
