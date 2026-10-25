## =====================================================================
##  Validation 01: does the brms model carry the prior the docs describe?
##
##  The Stan code brms generates for fit_flat_r_model() is sampled with its
##  likelihood switched off, so the draws come from flat_r_lpdf itself. They
##  are compared with what the prior should give. This checks the Stan code
##  directly, including that it compiles.
##
##    P1. one predictor:       the implied correlation is flat from -1 to 1
##    P2. 3 correlated:        the implied R-squared is Beta(1/2, 1)
##    P3. 8 predictors:        still Beta(1/2, 1): no drift with k
##    P4. 3 correlated:        each predictor's partial correlation matches the
##                             drawing recipe draw_coefficients_from_flat_r_prior()
##
##  Criterion: Kolmogorov-Smirnov distance below 0.03 (the largest gap
##  between the two cumulative distributions). With 20,000 MCMC draws a
##  correct prior lands well below that; a wrong factor lands far above it.
##  Divergent transitions are reported (the prior has a sharp peak at zero
##  for k >= 2, docs/2_how_it_calculates.md section 3).
##
##  Prior-only sampling needs proper priors on the intercept and the noise
##  (the default ones are improper). The checked quantities, the correlation
##  and R-squared, do not depend on the noise level, so any proper choice
##  gives the same answer.
##
##  Run from the project folder: Rscript validation/01_prior_in_brms.R
##  A few minutes (three models compile).
## =====================================================================
source("R/flat_r.R")
source("validation/settings.R")   # parallel cores, cache folders, record_check()
library(brms)

## ---- Prior-only fit ------------------------------------------------------
## brms refuses sample_prior = "only" when a prior is a custom target += line.
## So: take the exact Stan code and data brms generates, switch its own
## prior_only flag on (likelihood off), and sample that.
fit_prior_only <- function(formula, data) {
  prior_parts <- build_flat_r_prior(formula, data)
  priors <- set_prior("normal(0, 10)", class = "Intercept") +
            set_prior("lognormal(0, 1)", class = "sigma") +
            set_prior("target += flat_r_lpdf(b | sigma, predictor_covariance);", check = FALSE)
  stan_code <- make_stancode(formula, data = data, family = gaussian(), prior = priors,
                             stanvars = prior_parts$stanvars)
  stan_data <- make_standata(formula, data = data, family = gaussian(), prior = priors,
                             stanvars = prior_parts$stanvars)
  stan_data$prior_only <- 1L
  stan_data <- unclass(stan_data)
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_code))
    fit   <- model$sample(data = stan_data, chains = 4, iter_warmup = 5000,
                          iter_sampling = 5000, seed = 1, adapt_delta = 0.99, refresh = 0)
    list(draws = posterior::as_draws_matrix(fit$draws(c("b", "sigma"))),
         divergent = sum(fit$diagnostic_summary()$num_divergent))
  } else {
    fit <- rstan::stan(model_code = stan_code, data = stan_data, chains = 4, iter = 10000,
                       warmup = 5000, seed = 1, refresh = 0,
                       control = list(adapt_delta = 0.99))
    list(draws = posterior::as_draws_matrix(fit),
         divergent = rstan::get_num_divergent(fit))
  }
}

## Coefficients (one column per predictor) and sigma from a prior-only fit
coefficients_of <- function(fit) fit$draws[, grep("^b\\[", colnames(fit$draws)), drop = FALSE]
sigma_of        <- function(fit) as.numeric(fit$draws[, "sigma"])

## Implied f^2 (Cohen's f squared) per draw
implied_f_squared <- function(fit, predictor_covariance) {
  coefficients <- coefficients_of(fit)
  rowSums((coefficients %*% predictor_covariance) * coefficients) / sigma_of(fit)^2
}

report_divergences <- function(fit, label)
  message(sprintf("%s: %d divergent transitions out of 20000", label, fit$divergent))

## (MCMC repeats a value when a proposal is rejected; ks.test warns about those
## ties, which don't matter for the distance itself)
ks_distance <- function(draws, reference) suppressWarnings(unname(ks.test(draws, reference)$statistic))

## ---- Data (only the predictors matter: no likelihood is used) -----------
set.seed(42)
sample_size <- 60
three_correlated <- matrix(rnorm(sample_size * 3), sample_size) %*%
                    chol(matrix(c(1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1), 3, 3))
eight_uncorrelated <- matrix(rnorm(sample_size * 8), sample_size)
prior_data <- data.frame(y = rnorm(sample_size),
                         x  = rnorm(sample_size),
                         setNames(as.data.frame(three_correlated), paste0("c", 1:3)),
                         setNames(as.data.frame(eight_uncorrelated), paste0("u", 1:8)))

covariance_of <- function(columns) {
  centered <- scale(as.matrix(prior_data[, columns]), scale = FALSE)
  crossprod(centered) / nrow(centered)
}


## =====================================================================
## P1. One predictor: flat correlation
## =====================================================================
fit_one <- fit_prior_only(y ~ x, prior_data)
report_divergences(fit_one, "P1")
signed_f <- as.numeric(coefficients_of(fit_one)) * sqrt(covariance_of("x")[1, 1]) /
            sigma_of(fit_one)
implied_correlation <- signed_f / sqrt(1 + signed_f^2)
distance <- ks_distance(implied_correlation, function(r) punif(r, min = -1, max = 1))
record_check("P1. k = 1: correlation flat from -1 to 1", distance < 0.03,
             sprintf("KS distance %.3f", distance))


## =====================================================================
## P2. Three correlated predictors: R-squared ~ Beta(1/2, 1)
## =====================================================================
correlated_columns <- paste0("c", 1:3)
fit_three <- fit_prior_only(y ~ c1 + c2 + c3, prior_data)
report_divergences(fit_three, "P2")
f_squared_three <- implied_f_squared(fit_three, covariance_of(correlated_columns))
r_squared_three <- f_squared_three / (1 + f_squared_three)
distance <- ks_distance(r_squared_three, function(r2) pbeta(r2, shape1 = 0.5, shape2 = 1))
record_check("P2. k = 3 correlated: R-squared ~ Beta(1/2, 1)", distance < 0.03,
             sprintf("KS distance %.3f", distance))


## =====================================================================
## P3. Eight predictors: still Beta(1/2, 1)
## =====================================================================
uncorrelated_columns <- paste0("u", 1:8)
fit_eight <- fit_prior_only(reformulate(uncorrelated_columns, response = "y"), prior_data)
report_divergences(fit_eight, "P3")
f_squared_eight <- implied_f_squared(fit_eight, covariance_of(uncorrelated_columns))
r_squared_eight <- f_squared_eight / (1 + f_squared_eight)
distance <- ks_distance(r_squared_eight, function(r2) pbeta(r2, shape1 = 0.5, shape2 = 1))
record_check("P3. k = 8: R-squared ~ Beta(1/2, 1) (no drift)", distance < 0.03,
             sprintf("KS distance %.3f", distance))


## =====================================================================
## P4. Per-predictor shares: brms vs the drawing recipe
## =====================================================================
## Partial correlation of predictor j: its coefficient times the spread of x_j
## left after the other predictors, in noise units, turned into a correlation.
covariance_three <- covariance_of(correlated_columns)
spread_left_after_others <- sqrt(1 / diag(solve(covariance_three)))
partial_correlations <- function(coefficients_in_noise_units) {
  scaled <- sweep(coefficients_in_noise_units, 2, spread_left_after_others, "*")
  scaled / sqrt(1 + scaled^2)
}
brms_partials   <- partial_correlations(coefficients_of(fit_three) / sigma_of(fit_three))
set.seed(7)
recipe_partials <- partial_correlations(draw_coefficients_from_flat_r_prior(100000, covariance_three))
for (j in 1:3) {
  distance <- ks_distance(brms_partials[, j], recipe_partials[, j])
  record_check(sprintf("P4. partial correlation of c%d: brms vs recipe", j), distance < 0.03,
               sprintf("KS distance %.3f", distance))
}

report_checks()
