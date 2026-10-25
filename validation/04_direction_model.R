## =====================================================================
##  Validation 04: sampling in the prior's own coordinates
##
##  Check 01 found that sampling the coefficients directly (what brms does)
##  under-samples the peak at zero for larger k. R/flat_r_direction.stan
##  samples R (flat) and a direction (uniform) instead, and computes the
##  coefficients from them. This script checks that it gets the prior and
##  the posterior right, with the coefficient-coordinates model as contrast.
##
##    D1. prior: implied correlation flat (k = 1); R-squared ~ Beta(1/2, 1)
##        (k = 3 correlated, k = 8); per-predictor shares match the recipe (k = 3)
##        Criterion: KS distance < 0.03
##    D2. posterior: simulation-based calibration, k = 8, n = 60. Draw the
##        parameters from the prior, simulate data, fit, and record where the
##        true value falls among the posterior draws (its rank). If the
##        sampler is right, every rank is equally likely.
##        Criterion: chi-squared test of the 10 rank bins, p > 0.01
##
##  Both use proper stand-in priors on the intercept and the noise (needed to
##  sample a prior or simulate from it); R, f and the direction don't depend
##  on them.
##
##  Run from the project folder: Rscript validation/04_direction_model.R
##  About 5 minutes on 4 cores.
## =====================================================================
source("R/flat_r.R")
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
    posterior::as_draws_matrix(fit$draws(c("b", "sigma")))
  } else {
    fit <- rstan::sampling(model, data = stan_data, chains = chains,
                           iter = iter_warmup + iter_sampling, warmup = iter_warmup,
                           thin = thin, seed = seed, refresh = 0)
    as.matrix(fit, pars = c("b", "sigma"))
  }
}
coefficients_of <- function(draws) draws[, grep("^b\\[", colnames(draws)), drop = FALSE]

## ---- The two models ------------------------------------------------------
## direction:    R/flat_r_direction.stan (R and a direction as parameters)
## coefficients: the same prior as a density on b, as brms samples it
direction_model <- compile_stan(paste(readLines("R/flat_r_direction.stan"), collapse = "\n"))
coefficients_model <- compile_stan(paste0("
functions {", flat_r_stan_code, "}
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
## D1. The prior
## =====================================================================
set.seed(42)
one_predictor    <- matrix(rnorm(sample_size), sample_size)
three_correlated <- matrix(rnorm(sample_size * 3), sample_size) %*%
                    chol(matrix(c(1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1), 3, 3))
eight_predictors <- matrix(rnorm(sample_size * 8), sample_size)

sample_prior <- function(predictors)
  sample_stan(direction_model,
              list(N = sample_size, K = ncol(predictors), X = predictors,
                   y = rnorm(sample_size), prior_only = 1L, proper_stand_ins = 1L),
              chains = 4, iter_warmup = 1000, iter_sampling = 5000, thin = 1, seed = 1)

draws_one <- sample_prior(one_predictor)
signed_f  <- coefficients_of(draws_one)[, 1] * sqrt(covariance_of(one_predictor)[1, 1]) /
             draws_one[, "sigma"]
distance  <- ks_distance(signed_f / sqrt(1 + signed_f^2), function(r) punif(r, min = -1, max = 1))
record_check("D1a. k = 1: correlation flat from -1 to 1", distance < 0.03,
             sprintf("KS distance %.3f", distance))

draws_three <- sample_prior(three_correlated)
distance <- ks_distance(implied_r_squared(draws_three, covariance_of(three_correlated)),
                        function(r2) pbeta(r2, shape1 = 0.5, shape2 = 1))
record_check("D1b. k = 3 correlated: R-squared ~ Beta(1/2, 1)", distance < 0.03,
             sprintf("KS distance %.3f", distance))

draws_eight <- sample_prior(eight_predictors)
distance <- ks_distance(implied_r_squared(draws_eight, covariance_of(eight_predictors)),
                        function(r2) pbeta(r2, shape1 = 0.5, shape2 = 1))
record_check("D1c. k = 8: R-squared ~ Beta(1/2, 1)", distance < 0.03,
             sprintf("KS distance %.3f (coefficient coordinates in brms: 0.074)", distance))

covariance_three <- covariance_of(three_correlated)
spread_left_after_others <- sqrt(1 / diag(solve(covariance_three)))
partial_correlations <- function(coefficients_in_noise_units) {
  scaled <- sweep(coefficients_in_noise_units, 2, spread_left_after_others, "*")
  scaled / sqrt(1 + scaled^2)
}
stan_partials   <- partial_correlations(coefficients_of(draws_three) / draws_three[, "sigma"])
set.seed(7)
recipe_partials <- partial_correlations(draw_coefficients_from_flat_r_prior(100000, covariance_three))
distances <- sapply(1:3, function(j) ks_distance(stan_partials[, j], recipe_partials[, j]))
record_check("D1d. k = 3: per-predictor shares match the recipe", all(distances < 0.03),
             sprintf("KS distances %s", paste(sprintf("%.3f", distances), collapse = " ")))


## =====================================================================
## D2. The posterior: simulation-based calibration, k = 8
## =====================================================================
number_of_predictors <- 8
replications <- 400
set.seed(3)
predictors  <- matrix(rnorm(sample_size * number_of_predictors), sample_size)
centered    <- scale(predictors, scale = FALSE)
covariance  <- covariance_of(predictors)

one_replication <- function(replication) {
  set.seed(1000 + replication)
  true_alpha <- rnorm(1, mean = 0, sd = 10)
  true_sigma <- rlnorm(1, meanlog = 0, sdlog = 1)
  true_b     <- true_sigma * drop(draw_coefficients_from_flat_r_prior(1, covariance))
  outcome    <- true_alpha + drop(centered %*% true_b) + rnorm(sample_size, mean = 0, sd = true_sigma)
  true_f_squared <- sum((true_b %*% covariance) * true_b) / true_sigma^2
  true_r_squared <- true_f_squared / (1 + true_f_squared)
  stan_data <- list(N = sample_size, K = number_of_predictors, X = predictors, y = outcome,
                    prior_only = 0L, proper_stand_ins = 1L)
  ## 100 draws per fit (1000 thinned by 10), so ranks run from 0 to 100
  sapply(list(direction = direction_model, coefficients = coefficients_model), function(model) {
    draws <- suppressWarnings(sample_stan(model, stan_data, chains = 1, iter_warmup = 1000,
                                          iter_sampling = 1000, thin = 10, seed = replication))
    c(r_squared = sum(implied_r_squared(draws, covariance) < true_r_squared),
      b1        = sum(coefficients_of(draws)[, 1] < true_b[1]),
      sigma     = sum(draws[, "sigma"] < true_sigma))
  })
}
message("Simulation-based calibration: ", replications, " datasets x 2 models ...")
ranks <- parallel::mclapply(seq_len(replications), one_replication, mc.cores = number_of_cores)

for (model_name in c("direction", "coefficients")) for (quantity in c("r_squared", "b1", "sigma")) {
  quantity_ranks <- sapply(ranks, function(r) r[quantity, model_name])
  bin_counts <- table(cut(quantity_ranks, breaks = seq(0, 100, by = 10), include.lowest = TRUE))
  p_value <- chisq.test(bin_counts)$p.value
  label <- sprintf("D2. %s coordinates, %s", model_name, quantity)
  details <- sprintf("p = %.3f; bins %s", p_value, paste(bin_counts, collapse = " "))
  if (model_name == "direction") {
    record_check(label, p_value > 0.01, details)
  } else {
    message(label, " (contrast, not a check): ", details)
  }
}

report_checks()
