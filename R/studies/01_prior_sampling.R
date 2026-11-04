## =====================================================================
##  Study 1: is the prior sampled correctly?
##  (docs/appendix/C_study_design.md, study 1; replaces checks P1-P10 of the old study 01)
##
##  Aims. flatcor's Stan model (inst/stan/flatcor.stan in the installed
##  package, used by fit_regression()) samples R, with R-squared ~ Beta(a, 1),
##  and a uniform direction, and computes the coefficients from them. Does it
##  reproduce each of the three priors?
##  Data. None: the model samples the prior only (prior_only = 1), with proper
##  stand-in priors on the intercept and the noise. The predictors are fixed
##  designs (n = 60): one predictor, three correlated predictors
##  (correlation .4) and eight uncorrelated ones, each drawn once from its own
##  stream (R/data/data_sets.R).
##  Estimands. The distribution of R-squared, of the correlation (k = 1) and of
##  each coefficient's partial correlation, against its exact form.
##  Methods. The Stan model for flat_r (a = 1/2), flat_partial_r (a = k/2) and
##  flat_r_squared (a = 1).
##  Performance. The Kolmogorov-Smirnov distance between the draws and the
##  exact distribution, judged against the KS critical value at the .01 level
##  for the effective number of draws, 1.63 / sqrt(ESS) (two-sample: for the
##  ESS and the size of the reference sample). The ESS is posterior::ess_bulk()
##  of the checked quantity, so autocorrelated draws count for what they are
##  worth.
##  n_sim. Enough iterations that every checked quantity reaches an ESS of
##  n_sim_of("ks_draws") = 26,569 (critical distance <= .01;
##  R/shared/precision.R): the chains are lengthened until it does.
##
##  Checks (old number in brackets):
##    1.1  flat_r, k = 1: the correlation is flat on (-1, 1)              (P1)
##    1.2  flat_r, k = 3 correlated: R-squared ~ Beta(1/2, 1)             (P2)
##    1.3  flat_r, k = 8: R-squared ~ Beta(1/2, 1)                        (P3)
##    1.4  flat_r, k = 3 correlated: each partial correlation matches the
##         recipe draw_coefficients_from_prior() (two-sample)            (P4)
##    1.5  flat_partial_r, k = 3 correlated: R-squared ~ Beta(3/2, 1)     (P5)
##    1.6  flat_partial_r, k = 8: R-squared ~ Beta(4, 1)                  (P6)
##    1.7  flat_partial_r, k = 3 correlated: each partial correlation
##         flat on (-1, 1)                                                (P7)
##    1.8  flat_r_squared, k = 3 correlated: R-squared flat               (P8)
##    1.9  flat_r_squared, k = 8: R-squared flat                          (P9)
##    1.10 flat_r_squared, k = 1: r^2 flat (|r| has density 2|r|)         (P10)
##  The old threshold was a fixed KS distance of 0.03; each check now passes
##  when the distance is below the critical value for its ESS.
##
##  Run from the project folder:  Rscript R/studies/01_prior_sampling.R
##  Smoke test (ESS target 500, about a minute):
##    VALIDATION_SMOKE=1 Rscript R/studies/01_prior_sampling.R
##  Raw results are cached in cache/study_01/ (delete them to rerun); the
##  tables and checks are saved to results/studies/study_01.rds.
## =====================================================================
library(posterior)
source("R/shared/settings.R")             # flatcor, cores, record_check()
source("R/shared/precision.R")            # n_sim_of(), smoke_mode
source("R/data/data_sets.R")            # with_stream_of()

study_folder <- file.path(cache_folder, if (smoke_mode) "study_01_smoke" else "study_01")
dir.create(study_folder, recursive = TRUE, showWarnings = FALSE)

## the ESS every checked quantity must reach (smoke mode: 500)
target_ess <- if (smoke_mode) 500L else precision_targets$ks_draws$n_sim
chains <- 4


## ---- Sampling helper: cmdstanr if installed, otherwise rstan -------------------
## Returns the draws of b and sigma as a posterior::draws_array (iterations x
## chains x variables), so the ESS of any derived quantity can be computed.
stan_model_code <- paste(readLines(system.file("stan", "flatcor.stan", package = "flatcor")),
                         collapse = "\n")
use_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
prior_model <- if (use_cmdstanr) {
  cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_model_code))
} else {
  rstan::stan_model(model_code = stan_model_code)
}
sample_prior_draws <- function(stan_data, iter_sampling, seed) {
  if (use_cmdstanr) {
    fit <- prior_model$sample(data = stan_data, chains = chains, parallel_chains = min(chains, number_of_cores),
                              iter_warmup = 1000, iter_sampling = iter_sampling, seed = seed,
                              refresh = 0, show_messages = FALSE)
    fit$draws(c("b", "sigma"))
  } else {
    fit <- rstan::sampling(prior_model, data = stan_data, chains = chains,
                           iter = 1000 + iter_sampling, warmup = 1000, seed = seed, refresh = 0)
    as_draws_array(as.array(fit, pars = c("b", "sigma")))
  }
}


## ---- The fixed designs ----------------------------------------------------------
sample_size <- 60
covariance_of <- function(predictors) {
  centered <- scale(predictors, scale = FALSE)
  crossprod(centered) / nrow(centered)
}
designs <- list(
  one = with_stream_of("study01-design-one", matrix(rnorm(sample_size), sample_size)),
  three_correlated = with_stream_of("study01-design-three-correlated",
    matrix(rnorm(sample_size * 3), sample_size) %*%
      chol(matrix(c(1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1), 3, 3))),
  eight = with_stream_of("study01-design-eight", matrix(rnorm(sample_size * 8), sample_size)))


## ---- Quantities computed from the draws (each an iterations x chains matrix) ----
coefficient_draws <- function(draws, j) as.matrix(extract_variable_array(draws, sprintf("b[%d]", j))[, , 1])
sigma_draws <- function(draws) as.matrix(extract_variable_array(draws, "sigma")[, , 1])

## explained share: f^2 = b' Sigma b / sigma^2, R^2 = f^2 / (1 + f^2)
implied_r_squared <- function(draws, predictor_covariance) {
  k <- ncol(predictor_covariance)
  coefficients <- sapply(seq_len(k), function(j) as.vector(coefficient_draws(draws, j)))
  f_squared <- rowSums((coefficients %*% predictor_covariance) * coefficients) /
               as.vector(sigma_draws(draws))^2
  matrix(f_squared / (1 + f_squared), ncol = chains)
}

## each coefficient's partial correlation: the coefficient in noise units times
## the predictor's spread left after the others, rho = f_j / sqrt(1 + f_j^2)
partial_correlation <- function(draws, predictor_covariance, j) {
  spread_left <- sqrt(1 / diag(solve(predictor_covariance)))[j]
  f_j <- coefficient_draws(draws, j) / sigma_draws(draws) * spread_left
  f_j / sqrt(1 + f_j^2)
}

## the signed correlation with one predictor
correlation_one <- function(draws, predictor_covariance) {
  f <- coefficient_draws(draws, 1) * sqrt(predictor_covariance[1, 1]) / sigma_draws(draws)
  f / sqrt(1 + f^2)
}


## ---- Sampling until every checked quantity reaches the target ESS --------------
## Each setting (design and prior) is sampled with a fixed seed. If the
## smallest ESS of its checked quantities falls short of the target, the
## chains are lengthened in proportion and sampled again (same seed); the rule
## concerns precision only, not the outcome of any check.
settings <- list(
  list(name = "flat_r-one",                design = "one",              prior = "flat_r",         seed = 101),
  list(name = "flat_r-three_correlated",   design = "three_correlated", prior = "flat_r",         seed = 102),
  list(name = "flat_r-eight",              design = "eight",            prior = "flat_r",         seed = 103),
  list(name = "flat_partial_r-three_correlated", design = "three_correlated", prior = "flat_partial_r", seed = 105),
  list(name = "flat_partial_r-eight",      design = "eight",            prior = "flat_partial_r", seed = 106),
  list(name = "flat_r_squared-three_correlated", design = "three_correlated", prior = "flat_r_squared", seed = 108),
  list(name = "flat_r_squared-eight",      design = "eight",            prior = "flat_r_squared", seed = 109),
  list(name = "flat_r_squared-one",        design = "one",              prior = "flat_r_squared", seed = 110))

## the quantities whose ESS must reach the target, per setting
checked_quantities <- function(draws, setting) {
  covariance <- covariance_of(designs[[setting$design]])
  k <- ncol(covariance)
  ## k = 1: flat_r checks r itself, flat_r_squared checks r^2 (its r is
  ## bimodal near -1 and 1, so the sign mixes slowly while r^2 does not)
  if (k == 1) return(list(if (setting$prior == "flat_r_squared") correlation_one(draws, covariance)^2
                          else correlation_one(draws, covariance)))
  c(list(implied_r_squared(draws, covariance)),
    if (setting$design == "three_correlated")
      lapply(seq_len(k), function(j) partial_correlation(draws, covariance, j)))
}

sample_setting <- function(setting) {
  cache_file <- file.path(study_folder, paste0(setting$name, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  predictors <- designs[[setting$design]]
  stan_data <- list(N = sample_size, K = ncol(predictors), X = predictors,
                    y = with_stream_of(paste0("study01-outcome-", setting$design), rnorm(sample_size)),
                    prior_only = 1L, proper_stand_ins = 1L,
                    r_squared_shape = prior_shape(setting$prior, ncol(predictors)))
  iter_sampling <- ceiling(1.2 * target_ess / chains)
  for (attempt in 1:4) {
    draws <- sample_prior_draws(stan_data, iter_sampling, setting$seed)
    smallest_ess <- min(vapply(checked_quantities(draws, setting), ess_bulk, numeric(1)))
    if (smallest_ess >= target_ess) break
    iter_sampling <- ceiling(iter_sampling * 1.2 * target_ess / smallest_ess)
  }
  result <- list(draws = draws, iter_sampling = iter_sampling, smallest_ess = smallest_ess)
  saveRDS(result, cache_file)
  result
}

samples <- setNames(lapply(settings, sample_setting), vapply(settings, `[[`, "", "name"))
cat("\nIterations per chain and smallest ESS per setting (target", target_ess, "):\n")
print(data.frame(setting = names(samples),
                 iterations_per_chain = vapply(samples, `[[`, numeric(1), "iter_sampling"),
                              smallest_ess = round(vapply(samples, `[[`, numeric(1), "smallest_ess"))),
      row.names = FALSE)


## ---- The KS check against the critical value for the ESS ------------------------
## one-sample: critical distance 1.63 / sqrt(ESS); two-sample against a
## reference sample of size m: 1.63 * sqrt((ESS + m) / (ESS * m))
ks_against <- function(quantity, reference, reference_size = NULL) {
  distance <- suppressWarnings(unname(ks.test(as.vector(quantity), reference)$statistic))
  ess <- ess_bulk(quantity)
  critical <- if (is.null(reference_size)) 1.63 / sqrt(ess) else
    1.63 * sqrt((ess + reference_size) / (ess * reference_size))
  c(distance = distance, ess = ess, critical = critical)
}
record_ks <- function(label, results) {
  results <- if (is.list(results)) results else list(results)
  distances <- vapply(results, `[[`, numeric(1), "distance")
  criticals <- vapply(results, `[[`, numeric(1), "critical")
  record_check(label, all(distances < criticals),
               sprintf("KS %s; critical %s (ESS %s)",
                       paste(sprintf("%.4f", distances), collapse = " "),
                       paste(sprintf("%.4f", criticals), collapse = " "),
                       paste(round(vapply(results, `[[`, numeric(1), "ess")), collapse = " ")))
}
covariance <- lapply(designs, covariance_of)

## 1.1 (old P1) flat_r, k = 1: the correlation is flat on (-1, 1)
record_ks("1.1 flat_r, k = 1: the correlation is flat on (-1, 1)",
          ks_against(correlation_one(samples[["flat_r-one"]]$draws, covariance$one),
                     function(r) punif(r, min = -1, max = 1)))

## 1.2, 1.3 (old P2, P3) flat_r: R-squared ~ Beta(1/2, 1)
record_ks("1.2 flat_r, k = 3 correlated: R-squared ~ Beta(1/2, 1)",
          ks_against(implied_r_squared(samples[["flat_r-three_correlated"]]$draws, covariance$three_correlated),
                     function(r2) pbeta(r2, 0.5, 1)))
record_ks("1.3 flat_r, k = 8: R-squared ~ Beta(1/2, 1)",
          ks_against(implied_r_squared(samples[["flat_r-eight"]]$draws, covariance$eight),
                     function(r2) pbeta(r2, 0.5, 1)))

## 1.4 (old P4) flat_r, k = 3 correlated: each partial correlation matches the
## prior as a recipe (draw_coefficients_from_prior(), in noise units)
reference_size <- if (smoke_mode) 20000 else 200000
recipe <- with_stream_of("study01-recipe-flat_r-three",
  draw_coefficients_from_prior(reference_size, covariance$three_correlated, prior = "flat_r"))
spread_left <- sqrt(1 / diag(solve(covariance$three_correlated)))
recipe_partials <- sweep(recipe, 2, spread_left, "*")
recipe_partials <- recipe_partials / sqrt(1 + recipe_partials^2)
record_ks("1.4 flat_r, k = 3 correlated: each partial correlation matches the recipe",
          lapply(1:3, function(j)
            ks_against(partial_correlation(samples[["flat_r-three_correlated"]]$draws,
                                           covariance$three_correlated, j),
                       recipe_partials[, j], reference_size = reference_size)))

## 1.5, 1.6 (old P5, P6) flat_partial_r: R-squared ~ Beta(k/2, 1)
record_ks("1.5 flat_partial_r, k = 3 correlated: R-squared ~ Beta(3/2, 1)",
          ks_against(implied_r_squared(samples[["flat_partial_r-three_correlated"]]$draws,
                                       covariance$three_correlated),
                     function(r2) pbeta(r2, 3 / 2, 1)))
record_ks("1.6 flat_partial_r, k = 8: R-squared ~ Beta(4, 1)",
          ks_against(implied_r_squared(samples[["flat_partial_r-eight"]]$draws, covariance$eight),
                     function(r2) pbeta(r2, 4, 1)))

## 1.7 (old P7) flat_partial_r, k = 3 correlated: each partial correlation flat
record_ks("1.7 flat_partial_r, k = 3 correlated: each partial correlation flat on (-1, 1)",
          lapply(1:3, function(j)
            ks_against(partial_correlation(samples[["flat_partial_r-three_correlated"]]$draws,
                                           covariance$three_correlated, j),
                       function(r) punif(r, min = -1, max = 1))))

## 1.8, 1.9 (old P8, P9) flat_r_squared: R-squared flat
record_ks("1.8 flat_r_squared, k = 3 correlated: R-squared flat",
          ks_against(implied_r_squared(samples[["flat_r_squared-three_correlated"]]$draws,
                                       covariance$three_correlated), punif))
record_ks("1.9 flat_r_squared, k = 8: R-squared flat",
          ks_against(implied_r_squared(samples[["flat_r_squared-eight"]]$draws, covariance$eight), punif))

## 1.10 (old P10) flat_r_squared, k = 1: r^2 flat (|r| has density 2|r|)
record_ks("1.10 flat_r_squared, k = 1: r^2 flat (|r| has density 2|r|)",
          ks_against(correlation_one(samples[["flat_r_squared-one"]]$draws, covariance$one)^2, punif))

save_study_results("01", list(iterations = vapply(samples, `[[`, numeric(1), "iter_sampling"),
                              smallest_ess = vapply(samples, `[[`, numeric(1), "smallest_ess")))
report_checks()
