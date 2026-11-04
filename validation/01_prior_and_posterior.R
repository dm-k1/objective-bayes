## =====================================================================
##  Validation 01: does the Stan model get the prior and the posterior right?
##
##  inst/stan/flatcor.stan (used by fit_regression()) samples R, with
##  R-squared ~ Beta(a, 1), and a direction (uniform), and computes the
##  coefficients from them. This script checks that it gets the prior and the
##  posterior right, for each of the three priors (a = 1/2: flat_r; a = k/2:
##  flat_partial_r; a = 1: flat_r_squared). Simulation-based calibration
##  follows Talts et al. (2018).
##
##  Contrast (reported, not a check): the flat_r prior written as a density on
##  the coefficients and sampled directly, as a custom brms prior would do.
##  For larger k that sampler misses the peak of the prior at zero, without
##  divergent transitions (with brms at k = 8: KS distance 0.074).
##
##    P1-P4. flat_r: implied correlation flat (k = 1); R-squared ~
##        Beta(1/2, 1) (k = 3 correlated, k = 8); per-predictor shares match
##        the recipe (k = 3)
##    P5-P7. flat_partial_r: R-squared ~ Beta(k/2, 1) (k = 3 correlated,
##        k = 8); every coefficient's partial correlation flat on (-1, 1)
##        (k = 3 correlated)
##    P8-P10. flat_r_squared: R-squared flat (k = 3 correlated, k = 8); with
##        one predictor |r| has density 2|r| (k = 1)
##        Criterion: KS distance < 0.03
##    S. posterior: simulation-based calibration, k = 8, n = 60, all three priors. Draw the
##        parameters from the prior, simulate data, fit, and record where the
##        true value falls among the posterior draws (its rank). If the
##        sampler is right, every rank is equally likely.
##        Criterion: chi-squared test of the 10 rank bins, p > 0.01
##
##  Both use proper stand-in priors on the intercept and the noise (needed to
##  sample a prior or simulate from it); R, f and the direction don't depend
##  on them.
##
##  Run from the project folder: Rscript validation/01_prior_and_posterior.R
##  About 8 minutes on 4 cores.
## =====================================================================
pkgload::load_all(".", quiet = TRUE)   # the flatcor package, from this folder
source("validation/settings.R")   # parallel cores, record_check()

## ---- Sampling helper: cmdstanr if installed, otherwise rstan ------------
compile_stan <- function(stan_code) {
  if (requireNamespace("cmdstanr", quietly = TRUE))
    cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_code))
  else
    rstan::stan_model(model_code = stan_code)
}
sample_stan <- function(model, stan_data, chains, iter_warmup, iter_sampling, thin, seed) {
  if (inherits(model, "CmdStanModel")) {
    fit <- model$sample(data = stan_data, chains = chains, iter_warmup = iter_warmup,
                        iter_sampling = iter_sampling, thin = thin, seed = seed,
                        refresh = 0, show_messages = FALSE)
    unclass(posterior::as_draws_matrix(fit$draws(c("b", "sigma"))))   # plain matrix, as rstan gives
  } else {
    fit <- rstan::sampling(model, data = stan_data, chains = chains,
                           iter = iter_warmup + iter_sampling, warmup = iter_warmup,
                           thin = thin, seed = seed, refresh = 0)
    as.matrix(fit, pars = c("b", "sigma"))
  }
}
coefficients_of <- function(draws) draws[, grep("^b\\[", colnames(draws)), drop = FALSE]

## ---- The two models ------------------------------------------------------
## direction:    inst/stan/flatcor.stan (R and a direction as parameters)
## coefficients: the same prior as a density on b (docs/2_how_it_calculates.md,
##               section 5), sampled directly
direction_model <- compile_stan(paste(readLines("inst/stan/flatcor.stan"), collapse = "\n"))
flat_r_density_code <- "
  // P(b) = P(total) x P(this combination | total) x unit conversion
  real flat_r_lpdf(vector b, real sigma, matrix predictor_covariance) {
    int k = rows(b);
    real f_squared = quad_form(predictor_covariance, b) / square(sigma);
    real f = sqrt(f_squared);
    real log_p_total = -1.5 * log1p(f_squared);
    real log_p_combination_given_total = -(log(2) + 0.5 * k * log(pi()) - lgamma(0.5 * k)
                                           + (k - 1) * log(f));
    real log_conversion = 0.5 * log_determinant(predictor_covariance) - k * log(sigma);
    return log_p_total + log_p_combination_given_total + log_conversion;
  }
"
coefficients_model <- compile_stan(paste0("
functions {", flat_r_density_code, "}
data {
  int<lower=1> N; int<lower=1> K; matrix[N, K] X; vector[N] y;
  int<lower=0, upper=1> prior_only; int<lower=0, upper=1> proper_stand_ins;
}
transformed data {
  matrix[N, K] Xc;
  for (j in 1:K) Xc[, j] = X[, j] - mean(X[, j]);
  matrix[K, K] predictor_covariance = crossprod(Xc) / N;
}
parameters { real alpha; real<lower=0> sigma; vector[K] b; }
model {
  alpha ~ normal(0, 10);
  sigma ~ lognormal(0, 1);
  target += flat_r_lpdf(b | sigma, predictor_covariance);
  if (!prior_only) y ~ normal(alpha + Xc * b, sigma);
}
"))

sample_size <- 60
covariance_of <- function(predictors) {
  centered <- scale(predictors, scale = FALSE)
  crossprod(centered) / nrow(centered)
}
implied_r_squared <- function(draws, predictor_covariance) {
  coefficients <- coefficients_of(draws)
  f_squared <- rowSums((coefficients %*% predictor_covariance) * coefficients) /
               draws[, "sigma"]^2
  f_squared / (1 + f_squared)
}
ks_distance <- function(draws, reference) suppressWarnings(unname(ks.test(draws, reference)$statistic))


## =====================================================================
## P1-P4. The prior
## =====================================================================
set.seed(42)
one_predictor    <- matrix(rnorm(sample_size), sample_size)
three_correlated <- matrix(rnorm(sample_size * 3), sample_size) %*%
                    chol(matrix(c(1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1), 3, 3))
eight_predictors <- matrix(rnorm(sample_size * 8), sample_size)

sample_prior <- function(predictors, prior = "flat_r")
  sample_stan(direction_model,
              list(N = sample_size, K = ncol(predictors), X = predictors,
                   y = rnorm(sample_size), prior_only = 1L, proper_stand_ins = 1L,
                   r_squared_shape = prior_shape(prior, ncol(predictors))),
              chains = 4, iter_warmup = 1000, iter_sampling = 5000, thin = 1, seed = 1)

draws_one <- sample_prior(one_predictor)
signed_f  <- coefficients_of(draws_one)[, 1] * sqrt(covariance_of(one_predictor)[1, 1]) /
             draws_one[, "sigma"]
distance  <- ks_distance(signed_f / sqrt(1 + signed_f^2), function(r) punif(r, min = -1, max = 1))
record_check("P1. k = 1: correlation flat from -1 to 1", distance < 0.03,
             sprintf("KS distance %.3f", distance))

draws_three <- sample_prior(three_correlated)
distance <- ks_distance(implied_r_squared(draws_three, covariance_of(three_correlated)),
                        function(r2) pbeta(r2, shape1 = 0.5, shape2 = 1))
record_check("P2. k = 3 correlated: R-squared ~ Beta(1/2, 1)", distance < 0.03,
             sprintf("KS distance %.3f", distance))

draws_eight <- sample_prior(eight_predictors)
distance <- ks_distance(implied_r_squared(draws_eight, covariance_of(eight_predictors)),
                        function(r2) pbeta(r2, shape1 = 0.5, shape2 = 1))
record_check("P3. k = 8: R-squared ~ Beta(1/2, 1)", distance < 0.03,
             sprintf("KS distance %.3f", distance))

covariance_three <- covariance_of(three_correlated)
spread_left_after_others <- sqrt(1 / diag(solve(covariance_three)))
partial_correlations <- function(coefficients_in_noise_units) {
  scaled <- sweep(coefficients_in_noise_units, 2, spread_left_after_others, "*")
  scaled / sqrt(1 + scaled^2)
}
stan_partials   <- partial_correlations(coefficients_of(draws_three) / draws_three[, "sigma"])
set.seed(7)
recipe_partials <- partial_correlations(draw_coefficients_from_prior(100000, covariance_three, prior = "flat_r"))
distances <- sapply(1:3, function(j) ks_distance(stan_partials[, j], recipe_partials[, j]))
record_check("P4. k = 3: per-predictor shares match the recipe", all(distances < 0.03),
             sprintf("KS distances %s", paste(sprintf("%.3f", distances), collapse = " ")))

## P5-P7. flat_partial_r: R-squared ~ Beta(k/2, 1), and each coefficient's
## partial correlation flat (Jeffreys' flat correlation for every
## coefficient, whatever k is; validation/05, C4 and C7)
partial_three <- sample_prior(three_correlated, prior = "flat_partial_r")
distance <- ks_distance(implied_r_squared(partial_three, covariance_of(three_correlated)),
                        function(r2) pbeta(r2, shape1 = 3 / 2, shape2 = 1))
record_check("P5. flat_partial_r, k = 3 correlated: R-squared ~ Beta(3/2, 1)", distance < 0.03,
             sprintf("KS distance %.3f", distance))

partial_eight <- sample_prior(eight_predictors, prior = "flat_partial_r")
distance <- ks_distance(implied_r_squared(partial_eight, covariance_of(eight_predictors)),
                        function(r2) pbeta(r2, shape1 = 4, shape2 = 1))
record_check("P6. flat_partial_r, k = 8: R-squared ~ Beta(4, 1)", distance < 0.03,
             sprintf("KS distance %.3f", distance))

partial_r_draws <- partial_correlations(coefficients_of(partial_three) / partial_three[, "sigma"])
distances <- sapply(1:3, function(j)
  ks_distance(partial_r_draws[, j], function(r) punif(r, min = -1, max = 1)))
record_check("P7. flat_partial_r, k = 3 correlated: each partial correlation flat",
             all(distances < 0.03),
             sprintf("KS distances %s", paste(sprintf("%.3f", distances), collapse = " ")))

## P8-P10. flat_r_squared: R-squared flat for every k
squared_three <- sample_prior(three_correlated, prior = "flat_r_squared")
distance <- ks_distance(implied_r_squared(squared_three, covariance_of(three_correlated)), punif)
record_check("P8. flat_r_squared, k = 3 correlated: R-squared flat", distance < 0.03,
             sprintf("KS distance %.3f", distance))

squared_eight <- sample_prior(eight_predictors, prior = "flat_r_squared")
distance <- ks_distance(implied_r_squared(squared_eight, covariance_of(eight_predictors)), punif)
record_check("P9. flat_r_squared, k = 8: R-squared flat", distance < 0.03,
             sprintf("KS distance %.3f", distance))

squared_one <- sample_prior(one_predictor, prior = "flat_r_squared")
signed_f    <- coefficients_of(squared_one)[, 1] * sqrt(covariance_of(one_predictor)[1, 1]) /
               squared_one[, "sigma"]
distance <- ks_distance(signed_f^2 / (1 + signed_f^2), punif)
record_check("P10. flat_r_squared, k = 1: r^2 flat (|r| has density 2|r|)", distance < 0.03,
             sprintf("KS distance %.3f", distance))


## =====================================================================
## S. The posterior: simulation-based calibration, k = 8
## =====================================================================
number_of_predictors <- 8
replications <- 400
set.seed(3)
predictors  <- matrix(rnorm(sample_size * number_of_predictors), sample_size)
centered    <- scale(predictors, scale = FALSE)
covariance  <- covariance_of(predictors)

## Four fits per replication: the package's model with each of the three
## priors (data simulated from that prior), and flat_r in coefficient
## coordinates (the contrast). flat_r and the contrast share their datasets
## (seed offset 0); the other two priors get their own seeds. (Equal seeds for
## two fits of the same compiled model would also give them the same cmdstanr
## output file, which the first fit's clean-up deletes.)
calibration_fits <- list(
  flat_r         = list(model = direction_model,    prior = "flat_r",         seed_offset = 0),
  coefficients   = list(model = coefficients_model, prior = "flat_r",         seed_offset = 0),
  flat_partial_r = list(model = direction_model,    prior = "flat_partial_r", seed_offset = 100000),
  flat_r_squared = list(model = direction_model,    prior = "flat_r_squared", seed_offset = 200000))

one_replication <- function(replication) {
  sapply(calibration_fits, function(setup) {
    set.seed(1000 + replication + setup$seed_offset)
    true_alpha <- rnorm(1, mean = 0, sd = 10)
    true_sigma <- rlnorm(1, meanlog = 0, sdlog = 1)
    true_b     <- true_sigma * drop(draw_coefficients_from_prior(1, covariance, prior = setup$prior))
    outcome    <- true_alpha + drop(centered %*% true_b) + rnorm(sample_size, mean = 0, sd = true_sigma)
    true_f_squared <- sum((true_b %*% covariance) * true_b) / true_sigma^2
    true_r_squared <- true_f_squared / (1 + true_f_squared)
    stan_data <- list(N = sample_size, K = number_of_predictors, X = predictors, y = outcome,
                      prior_only = 0L, proper_stand_ins = 1L,
                      r_squared_shape = prior_shape(setup$prior, number_of_predictors))
    ## 100 draws per fit (1000 thinned by 10), so ranks run from 0 to 100
    draws <- suppressWarnings(sample_stan(setup$model, stan_data, chains = 1, iter_warmup = 1000,
                                          iter_sampling = 1000, thin = 10,
                                          seed = replication + setup$seed_offset))
    c(r_squared = sum(implied_r_squared(draws, covariance) < true_r_squared),
      b1        = sum(coefficients_of(draws)[, 1] < true_b[1]),
      sigma     = sum(draws[, "sigma"] < true_sigma))
  })
}
message("Simulation-based calibration: ", replications, " datasets x 4 fits ...")
ranks <- parallel_map(seq_len(replications), one_replication)

for (model_name in names(calibration_fits)) for (quantity in c("r_squared", "b1", "sigma")) {
  quantity_ranks <- sapply(ranks, function(r) r[quantity, model_name])
  bin_counts <- table(cut(quantity_ranks, breaks = seq(0, 100, by = 10), include.lowest = TRUE))
  p_value <- chisq.test(bin_counts)$p.value
  label <- if (model_name == "coefficients") sprintf("S. flat_r in coefficient coordinates, %s", quantity)
           else sprintf("S. %s, %s", model_name, quantity)
  details <- sprintf("p = %.3f; bins %s", p_value, paste(bin_counts, collapse = " "))
  if (model_name != "coefficients") {
    record_check(label, p_value > 0.01, details)
  } else {
    message(label, " (contrast, not a check): ", details)
  }
}

report_checks()
