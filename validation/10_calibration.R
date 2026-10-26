## =====================================================================
##  Validation 10: is the conditional test calibrated, and how far is it
##  from one prior for everything?
##
##  The conditional test (R/conditional.R, uniform correlations) puts a flat
##  prior on the tested coefficient's partial correlation and flat priors on
##  the other coefficients, so it is not a ratio of marginal likelihoods under
##  the joint prior. Two checks of what that costs:
##
##  Part 1, calibration. Datasets drawn from the joint uniform-correlations
##  prior (draw_coefficients_from_prior()), half of them with x1's coefficient
##  set to exactly 0. If the test is the Bayes factor of that generating
##  process, then among datasets with Bayes factor BF10 the share of nulls is
##  1 / (1 + BF10) (prior odds 1:1). Datasets are binned by log BF10 and the
##  share of nulls per bin is compared with the mean of 1 / (1 + BF10) in it,
##  with the Monte Carlo error of a sum of Bernoulli draws. Also a logistic
##  regression of "null" on log BF01: intercept 0 and slope 1 when calibrated.
##  The same with data from the JZS prior (Cauchy, scale sqrt(2)/4) shows that
##  calibration is relative to the prior that generated the data.
##    k = 3, 8;  n = 50, 200;  predictor correlation 0.3;  20,000 datasets
##    per cell and prior (10,000 nulls).
##
##  Part 2, coherence. On 20 datasets (k = 3, 8; n = 50, 200; five datasets
##  each, x1's t-value about 0, 1, 1.5, 2, 2.5; the other coefficients drawn
##  from the prior), the conditional Bayes factor against the Savage-Dickey
##  ratio for x1's partial correlation under the joint model fitted with
##  fit_regression(): prior density at 0 (1/2, flat on (-1, 1); validation
##  05, C4) over posterior density at 0 (logspline, kernel as a cross-check).
##  The Savage-Dickey ratio is the Bayes factor of the joint model against the
##  joint model restricted to a zero partial correlation: "one prior for
##  everything". Also shown: the joint test (bayes_factor_models(), the model
##  with against without x1).
##
##  Run from the project folder: Rscript validation/10_calibration.R
##  Part 1 a few minutes on 7 cores (cached); part 2 about 10-20 minutes
##  (each fit cached in validation/results/coherence_v1/). Needs cmdstanr
##  (or rstan) and logspline.
## =====================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
source("validation/settings.R")   # parallel_map(), caching, record_check()
pkgload::load_all(".", quiet = TRUE)


## ---- Shared: a fixed design and x1's partial correlation -----------------------
## Predictors are drawn once per dataset and then held fixed (the model
## conditions on them). The prior is set up in their sample covariance (/n),
## as the Stan model does, so x1's partial f is f1 = b1 * sd(x1 | others) / sigma.
draw_predictors <- function(n, k, correlation = 0.3) {
  covariance <- matrix(correlation, k, k); diag(covariance) <- 1
  predictors <- matrix(rnorm(n * k), n) %*% chol(covariance)
  scale(predictors, center = TRUE, scale = FALSE)
}
sample_covariance <- function(predictors) crossprod(predictors) / nrow(predictors)
sd_given_others   <- function(covariance, j = 1) sqrt(1 / solve(covariance)[j, j])

## x1's sample partial correlation, from its t-value
sample_partial_r <- function(predictors, outcome) {
  fit <- lm.fit(cbind(1, predictors), outcome)
  df  <- length(outcome) - ncol(predictors) - 1
  sigma_hat <- sqrt(sum(fit$residuals^2) / df)
  r_inverse <- chol2inv(qr.R(fit$qr))
  t_value   <- fit$coefficients[2] / (sigma_hat * sqrt(r_inverse[2, 2]))
  unname(t_value / sqrt(t_value^2 + df))
}


## =====================================================================
##  Part 1: calibration
## =====================================================================
jzs_scale <- sqrt(2) / 4

## Coefficients from the generating prior, for this design (noise SD 1)
draw_coefficients <- function(generating_prior, covariance) {
  if (generating_prior == "uniform correlations")
    return(drop(draw_coefficients_from_prior(1, covariance, prior = "uniform_correlations")))
  ## JZS (BayesFactor's g-prior): b | g ~ N(0, g Sigma^-1), g ~ Inverse-Gamma(1/2, s^2/2),
  ## so each coefficient's partial f is Cauchy(0, s)
  g <- 1 / rgamma(1, shape = 1 / 2, rate = jzs_scale^2 / 2)
  MASS::mvrnorm(1, rep(0, ncol(covariance)), g * solve(covariance))
}

calibration_design <- expand_grid(generating_prior     = c("uniform correlations", "JZS"),
                                  number_of_predictors = c(3, 8),
                                  sample_size          = c(50, 200),
                                  chunk                = 1:10)
datasets_per_chunk <- 2000   # 20,000 per cell, half nulls

calibration_runs <- run_simulation_grid(
  calibration_design,
  results_file = file.path(results_folder, "calibration.rds"),
  settings = list(datasets_per_chunk = datasets_per_chunk, jzs_scale = jzs_scale, version = 1),
  base_seed = 1010,
  simulate_cell = function(cell) {
    n <- cell$sample_size; k <- cell$number_of_predictors
    is_null <- rep(c(TRUE, FALSE), length.out = datasets_per_chunk)
    rows <- lapply(is_null, function(null) {
      predictors   <- draw_predictors(n, k)
      covariance   <- sample_covariance(predictors)
      coefficients <- draw_coefficients(cell$generating_prior, covariance)
      if (null) coefficients[1] <- 0
      outcome   <- drop(predictors %*% coefficients) + rnorm(n)
      partial_r <- sample_partial_r(predictors, outcome)
      true_f    <- coefficients[1] * sd_given_others(covariance)
      c(true_partial_correlation = true_f / sqrt(1 + true_f^2), partial_r = partial_r,
        log_BF10 = log_bayes_factor_uniform_correlations_conditional(partial_r, n, k))
    })
    cbind(as.data.frame(cell), is_null = is_null, as.data.frame(do.call(rbind, rows)))
  })

## ---- Per bin: observed share of nulls against the mean of 1 / (1 + BF10) --------
calibration_runs <- mutate(calibration_runs,
                           predicted_null = plogis(-log_BF10),   # 1 / (1 + BF10)
                           cell = sprintf("k = %d, n = %d", number_of_predictors, sample_size))
bin_edges <- c(-Inf, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 6, Inf)

calibration_bins <- calibration_runs |>
  mutate(bin = cut(log_BF10, bin_edges)) |>
  group_by(generating_prior, cell, bin) |>
  summarise(datasets = n(), observed = mean(is_null), predicted = mean(predicted_null),
            mean_log_BF10 = mean(log_BF10),
            ## Monte Carlo SD of the observed share if calibrated: Bernoulli(p_i) draws
            mc_sd = sqrt(sum(predicted_null * (1 - predicted_null))) / n(),
            .groups = "drop") |>
  mutate(z = (observed - predicted) / mc_sd)

## Chi-square over the bins with at least 5 expected nulls and 5 expected effects
calibration_tests <- calibration_bins |>
  filter(datasets * predicted >= 5, datasets * (1 - predicted) >= 5) |>
  group_by(generating_prior, cell) |>
  summarise(bins = n(), chi_square = sum(z^2), largest_abs_z = max(abs(z)), .groups = "drop") |>
  mutate(p_value = pchisq(chi_square, bins, lower.tail = FALSE))

## Logistic regression of "null" on log BF01: intercept 0, slope 1 when calibrated
logistic_fits <- calibration_runs |>
  group_by(generating_prior, cell) |>
  group_modify(function(rows, key) {
    ## (datasets with huge Bayes factors give fitted probabilities of 0 or 1: harmless)
    fit <- suppressWarnings(glm(is_null ~ I(-log_BF10), family = binomial, data = rows))
    estimates <- summary(fit)$coefficients
    data.frame(intercept = estimates[1, 1], intercept_se = estimates[1, 2],
               slope = estimates[2, 1], slope_se = estimates[2, 2])
  }) |>
  ungroup()

cat("\nCalibration per bin (share of nulls: observed, predicted = mean 1/(1 + BF10)):\n")
print(as.data.frame(calibration_bins |> select(generating_prior, cell, bin, datasets,
                                               observed, predicted, z)), digits = 3)
cat("\nChi-square over the bins, and the logistic calibration line:\n")
print(as.data.frame(left_join(calibration_tests, logistic_fits, by = c("generating_prior", "cell"))),
      digits = 3)

for (cell_name in unique(calibration_tests$cell)) {
  test <- filter(calibration_tests, generating_prior == "uniform correlations", cell == cell_name)
  line <- filter(logistic_fits, generating_prior == "uniform correlations", cell == cell_name)
  record_check(sprintf("P1. calibrated on its own prior, %s", cell_name),
               test$p_value > 0.01 && abs(line$slope - 1) < 3 * line$slope_se,
               sprintf("chi2 = %.1f on %d bins (p = %.2f); slope %.3f (SE %.3f)",
                       test$chi_square, test$bins, test$p_value, line$slope, line$slope_se))
}
for (cell_name in unique(calibration_tests$cell)) {
  test <- filter(calibration_tests, generating_prior == "JZS", cell == cell_name)
  line <- filter(logistic_fits, generating_prior == "JZS", cell == cell_name)
  record_check(sprintf("P2. data from JZS: miscalibration shows, %s", cell_name),
               test$p_value < 0.001,
               sprintf("chi2 = %.0f on %d bins; slope %.2f (SE %.2f)",
                       test$chi_square, test$bins, line$slope, line$slope_se))
}

## What the effects called null are, on the test's own prior (study 06: slope 0.1)
called_null <- calibration_runs |>
  filter(generating_prior == "uniform correlations", log_BF10 < -log(3)) |>
  group_by(cell) |>
  summarise(called_null = n(), share_true_nulls = mean(is_null),
            effects_with_abs_rho_below_0.1 = mean(abs(true_partial_correlation[!is_null]) < 0.1),
            effects_with_abs_rho_below_0.2 = mean(abs(true_partial_correlation[!is_null]) < 0.2),
            .groups = "drop")
cat("\nDatasets with BF01 > 3, data from the uniform-correlations prior:\n")
print(as.data.frame(called_null), digits = 3)

save_figure(
  ggplot(filter(calibration_bins, datasets >= 20),
         aes(predicted, observed, colour = generating_prior)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = predicted - 2 * mc_sd, ymax = predicted + 2 * mc_sd),
                  width = 0, colour = "grey70") +
    geom_point(size = 2) + geom_line() +
    facet_wrap(~ cell, nrow = 1) +
    scale_colour_manual(values = c(`uniform correlations` = "#1b9e77", JZS = "grey40")) +
    coord_equal() + theme(legend.position = "bottom") +
    labs(x = "predicted share of nulls, mean of 1 / (1 + BF10)", y = "observed share of nulls",
         colour = "data drawn from",
         title = "Calibration of the conditional test (bins of log BF10, grey: 2 Monte Carlo SD)"),
  "10_calibration", width = 12, height = 4.5)


## =====================================================================
##  Part 2: coherence with the joint model (Savage-Dickey)
## =====================================================================
coherence_folder <- file.path(results_folder, "coherence_v1")
dir.create(coherence_folder, recursive = TRUE, showWarnings = FALSE)
coherence_design <- expand_grid(number_of_predictors = c(3, 8), sample_size = c(50, 200),
                                target_t = c(0, 1, 1.5, 2, 2.5)) |>
  mutate(seed = 5000 + row_number())

coherence_one_dataset <- function(k, n, target_t, seed) {
  cache_file <- file.path(coherence_folder, sprintf("k%d_n%d_t%.1f.rds", k, n, target_t))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  set.seed(seed)
  predictors   <- draw_predictors(n, k)
  covariance   <- sample_covariance(predictors)
  coefficients <- drop(draw_coefficients_from_prior(1, covariance, prior = "uniform_correlations"))
  coefficients[1] <- target_t / sqrt(n) / sd_given_others(covariance)   # x1's f = t / sqrt(n)
  colnames(predictors) <- paste0("x", seq_len(k))
  data <- data.frame(predictors, y = drop(predictors %*% coefficients) + rnorm(n))
  formula <- reformulate(colnames(predictors), response = "y")

  ## the conditional test and the joint test
  partial_r   <- sample_partial_r(predictors, data$y)
  conditional <- log_bayes_factor_uniform_correlations_conditional(partial_r, n, k)
  joint <- log_bayes_factor_vs_intercept_only(sample_r_squared(predictors, data$y), n, k,
                                              "uniform_correlations") -
    log_bayes_factor_vs_intercept_only(sample_r_squared(predictors[, -1, drop = FALSE], data$y),
                                       n, k - 1, "uniform_correlations")

  ## Savage-Dickey: the joint model's posterior for x1's partial correlation
  fit <- suppressWarnings(fit_regression(formula, data, chains = 4, iter_warmup = 1000,
                                         iter_sampling = 5000, seed = seed))
  draws <- as.data.frame(fit$draws)
  partial_f   <- draws$x1 * sd_given_others(sample_covariance(scale(predictors, scale = FALSE))) /
    draws$sigma
  partial_rho <- partial_f / sqrt(1 + partial_f^2)
  log_density_at_0 <- function(values)
    log(logspline::dlogspline(0, logspline::logspline(values, lbound = -1, ubound = 1)))
  by_chain <- vapply(split(partial_rho, fit$chain_ids), log_density_at_0, numeric(1))
  kernel   <- density(partial_rho, bw = "SJ", n = 2^12, from = -1, to = 1)
  result <- data.frame(
    number_of_predictors = k, sample_size = n, target_t = target_t,
    partial_r = partial_r, conditional = conditional, joint = joint,
    savage_dickey        = log(1 / 2) - log_density_at_0(partial_rho),
    savage_dickey_se     = sd(by_chain) / sqrt(length(by_chain)),
    savage_dickey_kernel = log(1 / 2) - log(approx(kernel$x, kernel$y, 0)$y),
    divergent = fit$divergent,
    max_rhat  = max(posterior::summarise_draws(fit$draws, rhat = posterior::rhat)$rhat, na.rm = TRUE))
  saveRDS(result, cache_file)
  result
}

message("Part 2: ", nrow(coherence_design), " fits (cached in ", coherence_folder, ") ...")
coherence <- bind_rows(lapply(seq_len(nrow(coherence_design)), function(row)
  with(coherence_design[row, ], coherence_one_dataset(number_of_predictors, sample_size,
                                                      target_t, seed)))) |>
  mutate(conditional_minus_sd = conditional - savage_dickey,
         joint_minus_sd       = joint - savage_dickey)

cat("\nLog BF10 for x1: conditional test, Savage-Dickey under the joint model, joint test\n")
print(as.data.frame(coherence |> select(number_of_predictors, sample_size, target_t, partial_r,
                                        conditional, savage_dickey, savage_dickey_se,
                                        savage_dickey_kernel, joint, conditional_minus_sd,
                                        joint_minus_sd, divergent, max_rhat)), digits = 3)
cat(sprintf("\nconditional - Savage-Dickey: mean %.3f, mean absolute %.3f, largest absolute %.3f\n",
            mean(coherence$conditional_minus_sd), mean(abs(coherence$conditional_minus_sd)),
            max(abs(coherence$conditional_minus_sd))))
cat(sprintf("joint - Savage-Dickey:       mean %.3f, mean absolute %.3f, largest absolute %.3f\n",
            mean(coherence$joint_minus_sd), mean(abs(coherence$joint_minus_sd)),
            max(abs(coherence$joint_minus_sd))))
print(as.data.frame(coherence |> group_by(number_of_predictors, sample_size) |>
  summarise(mean_conditional_minus_sd = mean(conditional_minus_sd),
            largest_abs = max(abs(conditional_minus_sd)), .groups = "drop")), digits = 3)

## The posterior density at 0 can only be estimated from draws where there are
## draws near 0. Where x1's evidence is strong (log BF10 above about 4) the
## estimate is poor (Monte Carlo SE from the four chains above 0.1); those
## datasets are compared within their Monte Carlo error instead.
estimable <- coherence$savage_dickey_se < 0.1
cat(sprintf("\nSavage-Dickey estimable (chain SE < 0.1): %d of %d datasets\n",
            sum(estimable), nrow(coherence)))
cat(sprintf("conditional - Savage-Dickey, estimable: mean %.3f, mean absolute %.3f, largest absolute %.3f\n",
            mean(coherence$conditional_minus_sd[estimable]),
            mean(abs(coherence$conditional_minus_sd[estimable])),
            max(abs(coherence$conditional_minus_sd[estimable]))))
cat(sprintf("joint - Savage-Dickey, estimable:       mean %.3f, mean absolute %.3f, largest absolute %.3f\n",
            mean(coherence$joint_minus_sd[estimable]), mean(abs(coherence$joint_minus_sd[estimable])),
            max(abs(coherence$joint_minus_sd[estimable]))))

record_check("P3. Savage-Dickey: logspline and kernel agree (estimable ones)",
             max(abs(coherence$savage_dickey - coherence$savage_dickey_kernel)[estimable]) < 0.25,
             sprintf("%d of %d estimable; largest difference %.3f",
                     sum(estimable), nrow(coherence),
                     max(abs(coherence$savage_dickey - coherence$savage_dickey_kernel)[estimable])))
record_check("P4. conditional within a factor 2 of Savage-Dickey (estimable)",
             max(abs(coherence$conditional_minus_sd[estimable])) < log(2),
             sprintf("largest |difference| %.3f on the log scale (factor %.2f); mean |difference| %.3f",
                     max(abs(coherence$conditional_minus_sd[estimable])),
                     exp(max(abs(coherence$conditional_minus_sd[estimable]))),
                     mean(abs(coherence$conditional_minus_sd[estimable]))))
record_check("P5. the others: conditional within 3 Monte Carlo SE of Savage-Dickey",
             all(abs(coherence$conditional_minus_sd[!estimable]) <
                   3 * coherence$savage_dickey_se[!estimable]),
             paste(sprintf("%.2f vs %.2f (SE %.2f)", coherence$conditional[!estimable],
                           coherence$savage_dickey[!estimable], coherence$savage_dickey_se[!estimable]),
                   collapse = "; "))

save_figure(
  ggplot(coherence |> pivot_longer(c(conditional, joint), names_to = "test", values_to = "log_BF10"),
         aes(savage_dickey, log_BF10, colour = test, shape = factor(sample_size))) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = savage_dickey - 2 * savage_dickey_se,
                       xmax = savage_dickey + 2 * savage_dickey_se), height = 0, alpha = 0.5) +
    geom_point(size = 2.5) +
    facet_wrap(~ number_of_predictors, labeller = labeller(number_of_predictors = function(v) paste("k =", v))) +
    scale_colour_manual(values = c(conditional = "#1b9e77", joint = "#d95f02")) +
    labs(x = "Savage-Dickey log BF10 under the joint model", y = "log BF10",
         colour = NULL, shape = "n",
         title = "The conditional and the joint test against one prior for everything"),
  "10_coherence", width = 10, height = 4.5)

report_checks()
