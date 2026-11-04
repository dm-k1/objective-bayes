## =====================================================================
##  Validation 10: is the conditional test calibrated, and how far is it
##  from one prior for everything?
##
##  The conditional test (R/conditional.R) puts the prior that the joint prior
##  implies for one coefficient on the tested coefficient's partial
##  correlation, and flat priors on the other coefficients, so it is not a
##  ratio of marginal likelihoods under the joint prior. Two checks of what
##  that costs, for the three contenders (flat partial r, flat R, flat
##  R-squared; validation/comparison_methods.R):
##
##  Part 1, calibration. Datasets drawn from each contender's joint prior
##  (draw_coefficients_from_prior()), half of them with x1's coefficient set to
##  exactly 0, and tested with that contender's conditional test. If the test
##  is the Bayes factor of that generating process, then among datasets with
##  Bayes factor BF10 the share of nulls is 1 / (1 + BF10) (prior odds 1:1).
##  Datasets are binned by log BF10 and the share of nulls per bin is compared
##  with the mean of 1 / (1 + BF10) in it, with the Monte Carlo error of a sum
##  of Bernoulli draws. Also a logistic regression of "null" on log BF01:
##  intercept 0 and slope 1 when calibrated. Data from the JZS prior (Cauchy,
##  scale sqrt(2)/4), tested with each of the three, show that calibration is
##  relative to the prior that generated the data.
##    k = 3, 8;  n = 50, 200;  predictor correlation 0.3;  20,000 datasets
##    per cell and generating prior (10,000 nulls).
##
##  Part 2, coherence. On 20 datasets per prior (k = 3, 8; n = 50, 200; five
##  datasets each, x1's t-value about 0, 1, 1.5, 2, 2.5; the other
##  coefficients drawn from the prior), the conditional Bayes factor against
##  the Savage-Dickey ratio for x1's partial correlation under the joint model
##  fitted with fit_regression(): the prior density of the partial correlation
##  at 0 over its posterior density at 0 (logspline, kernel as a
##  cross-check). The Savage-Dickey ratio is the Bayes factor of the joint
##  model against the joint model restricted to a zero partial correlation:
##  "one prior for everything". Also shown: the joint test (the model with
##  against without x1). Prior density at 0:
##    flat partial r   1/2 (flat on (-1, 1); validation 05, C4)
##    flat R-squared   half the density of |f| at 0, from the share mixture
##                     (finite: pi/4 at k = 3)
##    flat R           not included: its implied prior for one partial
##                     correlation has an infinite density at 0 (for k >= 2 the
##                     share mixture E[s^(-1/2)] diverges), so the Savage-Dickey
##                     ratio is not defined.
##
##  Run from the project folder: Rscript validation/10_calibration.R
##  Part 1 about 20 minutes on 7 cores (cached); part 2 about 20-40 minutes
##  (each fit cached in cache/coherence_three_priors/). Needs
##  cmdstanr (or rstan) and logspline.
## =====================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
source("validation/settings.R")   # parallel_map(), caching, record_check()
pkgload::load_all(".", quiet = TRUE)
source("validation/comparison_methods.R")   # contenders, method_colours


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
  if (generating_prior %in% names(contenders))
    return(drop(draw_coefficients_from_prior(1, covariance, prior = contenders[[generating_prior]])))
  ## JZS (BayesFactor's g-prior): b | g ~ N(0, g Sigma^-1), g ~ Inverse-Gamma(1/2, s^2/2),
  ## so each coefficient's partial f is Cauchy(0, s)
  g <- 1 / rgamma(1, shape = 1 / 2, rate = jzs_scale^2 / 2)
  MASS::mvrnorm(1, rep(0, ncol(covariance)), g * solve(covariance))
}

calibration_design <- expand_grid(generating_prior     = c(names(contenders), "JZS"),
                                  number_of_predictors = c(3, 8),
                                  sample_size          = c(50, 200),
                                  chunk                = 1:10)
datasets_per_chunk <- 2000   # 20,000 per cell, half nulls

## Data from a contender are tested with that contender's conditional test;
## data from JZS with each of the three.
calibration_runs <- run_simulation_grid(
  calibration_design,
  results_file = file.path(cache_folder, "calibration.rds"),
  settings = list(datasets_per_chunk = datasets_per_chunk, jzs_scale = jzs_scale,
                  contenders = contenders, version = 2),
  base_seed = 1010,
  simulate_cell = function(cell) {
    n <- cell$sample_size; k <- cell$number_of_predictors
    tests <- if (cell$generating_prior == "JZS") names(contenders) else cell$generating_prior
    is_null <- rep(c(TRUE, FALSE), length.out = datasets_per_chunk)
    rows <- lapply(is_null, function(null) {
      predictors   <- draw_predictors(n, k)
      covariance   <- sample_covariance(predictors)
      coefficients <- draw_coefficients(cell$generating_prior, covariance)
      if (null) coefficients[1] <- 0
      outcome   <- drop(predictors %*% coefficients) + rnorm(n)
      partial_r <- sample_partial_r(predictors, outcome)
      true_f    <- coefficients[1] * sd_given_others(covariance)
      data.frame(is_null = null, test = tests,
                 true_partial_correlation = true_f / sqrt(1 + true_f^2), partial_r = partial_r,
                 log_BF10 = vapply(tests, function(test)
                   log_bayes_factor_conditional(partial_r, n, k, contenders[[test]]), numeric(1)))
    })
    cbind(as.data.frame(cell), do.call(rbind, rows), row.names = NULL)
  })

## ---- Per bin: observed share of nulls against the mean of 1 / (1 + BF10) --------
calibration_runs <- mutate(calibration_runs,
                           predicted_null = plogis(-log_BF10),   # 1 / (1 + BF10)
                           cell = sprintf("k = %d, n = %d", number_of_predictors, sample_size))
bin_edges <- c(-Inf, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 6, Inf)

calibration_bins <- calibration_runs |>
  mutate(bin = cut(log_BF10, bin_edges)) |>
  group_by(generating_prior, test, cell, bin) |>
  summarise(datasets = n(), observed = mean(is_null), predicted = mean(predicted_null),
            mean_log_BF10 = mean(log_BF10),
            ## Monte Carlo SD of the observed share if calibrated: Bernoulli(p_i) draws
            mc_sd = sqrt(sum(predicted_null * (1 - predicted_null))) / n(),
            .groups = "drop") |>
  mutate(z = (observed - predicted) / mc_sd)

## Chi-square over the bins with at least 5 expected nulls and 5 expected effects
calibration_tests <- calibration_bins |>
  filter(datasets * predicted >= 5, datasets * (1 - predicted) >= 5) |>
  group_by(generating_prior, test, cell) |>
  summarise(bins = n(), chi_square = sum(z^2), largest_abs_z = max(abs(z)), .groups = "drop") |>
  mutate(p_value = pchisq(chi_square, bins, lower.tail = FALSE))

## Logistic regression of "null" on log BF01: intercept 0, slope 1 when calibrated
logistic_fits <- calibration_runs |>
  group_by(generating_prior, test, cell) |>
  group_modify(function(rows, key) {
    ## (datasets with huge Bayes factors give fitted probabilities of 0 or 1: harmless)
    fit <- suppressWarnings(glm(is_null ~ I(-log_BF10), family = binomial, data = rows))
    estimates <- summary(fit)$coefficients
    data.frame(intercept = estimates[1, 1], intercept_se = estimates[1, 2],
               slope = estimates[2, 1], slope_se = estimates[2, 2])
  }) |>
  ungroup()

cat("\nCalibration per bin (share of nulls: observed, predicted = mean 1/(1 + BF10)):\n")
print(as.data.frame(calibration_bins |> select(generating_prior, test, cell, bin, datasets,
                                               observed, predicted, z)), digits = 3)
cat("\nChi-square over the bins, and the logistic calibration line:\n")
print(as.data.frame(left_join(calibration_tests, logistic_fits,
                              by = c("generating_prior", "test", "cell"))), digits = 3)

for (contender in names(contenders)) for (cell_name in unique(calibration_tests$cell)) {
  test <- filter(calibration_tests, generating_prior == contender, test == contender, cell == cell_name)
  line <- filter(logistic_fits, generating_prior == contender, test == contender, cell == cell_name)
  record_check(sprintf("P1. %s: calibrated on its own prior, %s", contender, cell_name),
               test$p_value > 0.01 && abs(line$slope - 1) < 3 * line$slope_se,
               sprintf("chi2 = %.1f on %d bins (p = %.2f); slope %.3f (SE %.3f)",
                       test$chi_square, test$bins, test$p_value, line$slope, line$slope_se))
}
for (contender in names(contenders)) for (cell_name in unique(calibration_tests$cell)) {
  test <- filter(calibration_tests, generating_prior == "JZS", test == contender, cell == cell_name)
  line <- filter(logistic_fits, generating_prior == "JZS", test == contender, cell == cell_name)
  record_check(sprintf("P2. %s, data from JZS: miscalibration shows, %s", contender, cell_name),
               test$p_value < 0.001,
               sprintf("chi2 = %.0f on %d bins; slope %.2f (SE %.2f)",
                       test$chi_square, test$bins, line$slope, line$slope_se))
}

## What the effects called null are, on the test's own prior (study 06: slope 0.1)
called_null <- calibration_runs |>
  filter(generating_prior == test, log_BF10 < -log(3)) |>
  group_by(test, cell) |>
  summarise(called_null = n(), share_true_nulls = mean(is_null),
            effects_with_abs_rho_below_0.1 = mean(abs(true_partial_correlation[!is_null]) < 0.1),
            effects_with_abs_rho_below_0.2 = mean(abs(true_partial_correlation[!is_null]) < 0.2),
            .groups = "drop")
cat("\nDatasets with BF01 > 3, each contender on data from its own prior:\n")
print(as.data.frame(called_null), digits = 3)

save_figure(
  ggplot(filter(calibration_bins, datasets >= 20) |>
           mutate(data_from = ifelse(generating_prior == "JZS", "JZS", "its own prior")),
         aes(predicted, observed, colour = data_from)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = predicted - 2 * mc_sd, ymax = predicted + 2 * mc_sd),
                  width = 0, colour = "grey70") +
    geom_point(size = 2) + geom_line() +
    facet_grid(test ~ cell) +
    scale_colour_manual(values = c(`its own prior` = "#0072B2", JZS = "grey40")) +
    coord_equal() + theme(legend.position = "bottom") +
    labs(x = "predicted share of nulls, mean of 1 / (1 + BF10)", y = "observed share of nulls",
         colour = "data drawn from",
         title = "Calibration of the conditional test (bins of log BF10, grey: 2 Monte Carlo SD)",
         subtitle = "rows: the contender whose conditional test is used"),
  "10_calibration", width = 12, height = 10)


## =====================================================================
##  Part 2: coherence with the joint model (Savage-Dickey)
## =====================================================================
coherence_folder <- file.path(cache_folder, "coherence_three_priors")
dir.create(coherence_folder, recursive = TRUE, showWarnings = FALSE)
coherence_priors <- c("flat partial r", "flat R-squared")   # flat R: see the header
coherence_design <- expand_grid(prior = coherence_priors, number_of_predictors = c(3, 8),
                                sample_size = c(50, 200), target_t = c(0, 1, 1.5, 2, 2.5)) |>
  mutate(seed = 5000 + row_number())

## Prior density of x1's partial correlation at 0 (both signs): flat partial r
## 1/2; otherwise half the density of |f| at 0 (d rho / d f = 1 there), the
## slope of the share mixture's CDF at a tiny f
log_prior_density_at_0 <- function(prior, k) {
  if (prior == "flat partial r") return(log(1 / 2))
  tiny <- 1e-6
  log(prior_cdf_of_f(tiny, k, points_share = 4000, shape = prior_shape(contenders[[prior]], k)) /
        tiny / 2)
}

coherence_one_dataset <- function(prior, k, n, target_t, seed) {
  cache_file <- file.path(coherence_folder, sprintf("%s_k%d_n%d_t%.1f.rds", contenders[[prior]], k, n,
                                                    target_t))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  set.seed(seed)
  predictors   <- draw_predictors(n, k)
  covariance   <- sample_covariance(predictors)
  coefficients <- drop(draw_coefficients_from_prior(1, covariance, prior = contenders[[prior]]))
  coefficients[1] <- target_t / sqrt(n) / sd_given_others(covariance)   # x1's f = t / sqrt(n)
  colnames(predictors) <- paste0("x", seq_len(k))
  data <- data.frame(predictors, y = drop(predictors %*% coefficients) + rnorm(n))
  formula <- reformulate(colnames(predictors), response = "y")

  ## the conditional test and the joint test
  partial_r   <- sample_partial_r(predictors, data$y)
  conditional <- log_bayes_factor_conditional(partial_r, n, k, contenders[[prior]])
  joint <- log_bayes_factor_vs_intercept_only(sample_r_squared(predictors, data$y), n, k,
                                              contenders[[prior]]) -
    log_bayes_factor_vs_intercept_only(sample_r_squared(predictors[, -1, drop = FALSE], data$y),
                                       n, k - 1, contenders[[prior]])

  ## Savage-Dickey: the joint model's posterior for x1's partial correlation
  fit <- suppressWarnings(fit_regression(formula, data, prior = contenders[[prior]], chains = 4,
                                         iter_warmup = 1000, iter_sampling = 5000, seed = seed))
  draws <- as.data.frame(fit$draws)
  partial_f   <- draws$x1 * sd_given_others(sample_covariance(scale(predictors, scale = FALSE))) /
    draws$sigma
  partial_rho <- partial_f / sqrt(1 + partial_f^2)
  log_density_at_0 <- function(values)
    log(logspline::dlogspline(0, logspline::logspline(values, lbound = -1, ubound = 1)))
  by_chain <- vapply(split(partial_rho, fit$chain_ids), log_density_at_0, numeric(1))
  kernel   <- density(partial_rho, bw = "SJ", n = 2^12, from = -1, to = 1)
  prior_at_0 <- log_prior_density_at_0(prior, k)
  result <- data.frame(
    prior = prior, number_of_predictors = k, sample_size = n, target_t = target_t,
    partial_r = partial_r, conditional = conditional, joint = joint,
    savage_dickey        = prior_at_0 - log_density_at_0(partial_rho),
    savage_dickey_se     = sd(by_chain) / sqrt(length(by_chain)),
    savage_dickey_kernel = prior_at_0 - log(approx(kernel$x, kernel$y, 0)$y),
    divergent = fit$divergent,
    max_rhat  = max(posterior::summarise_draws(fit$draws, rhat = posterior::rhat)$rhat, na.rm = TRUE))
  saveRDS(result, cache_file)
  result
}

message("Part 2: ", nrow(coherence_design), " fits (cached in ", coherence_folder, ") ...")
coherence <- bind_rows(lapply(seq_len(nrow(coherence_design)), function(row)
  with(coherence_design[row, ], coherence_one_dataset(prior, number_of_predictors, sample_size,
                                                      target_t, seed)))) |>
  mutate(conditional_minus_sd = conditional - savage_dickey,
         joint_minus_sd       = joint - savage_dickey)

cat("\nLog BF10 for x1: conditional test, Savage-Dickey under the joint model, joint test\n")
print(as.data.frame(coherence |> select(prior, number_of_predictors, sample_size, target_t, partial_r,
                                        conditional, savage_dickey, savage_dickey_se,
                                        savage_dickey_kernel, joint, conditional_minus_sd,
                                        joint_minus_sd, divergent, max_rhat)), digits = 3)
cat(sprintf("\nconditional - Savage-Dickey: mean %.3f, mean absolute %.3f, largest absolute %.3f\n",
            mean(coherence$conditional_minus_sd), mean(abs(coherence$conditional_minus_sd)),
            max(abs(coherence$conditional_minus_sd))))
cat(sprintf("joint - Savage-Dickey:       mean %.3f, mean absolute %.3f, largest absolute %.3f\n",
            mean(coherence$joint_minus_sd), mean(abs(coherence$joint_minus_sd)),
            max(abs(coherence$joint_minus_sd))))
print(as.data.frame(coherence |> group_by(prior, number_of_predictors, sample_size) |>
  summarise(mean_conditional_minus_sd = mean(conditional_minus_sd),
            largest_abs = max(abs(conditional_minus_sd)), .groups = "drop")), digits = 3)

## The posterior density at 0 can only be estimated from draws where there are
## draws near 0. Where x1's evidence is strong (log BF10 above about 4) the
## estimate is poor (Monte Carlo SE from the four chains above 0.1); those
## datasets are compared within their Monte Carlo error instead.
for (prior in coherence_priors) {
  of_prior  <- filter(coherence, prior == !!prior)
  estimable <- of_prior$savage_dickey_se < 0.1
  cat(sprintf("\n%s: Savage-Dickey estimable (chain SE < 0.1): %d of %d datasets\n",
              prior, sum(estimable), nrow(of_prior)))
  cat(sprintf("conditional - Savage-Dickey, estimable: mean %.3f, mean absolute %.3f, largest absolute %.3f\n",
              mean(of_prior$conditional_minus_sd[estimable]),
              mean(abs(of_prior$conditional_minus_sd[estimable])),
              max(abs(of_prior$conditional_minus_sd[estimable]))))
  cat(sprintf("joint - Savage-Dickey, estimable:       mean %.3f, mean absolute %.3f, largest absolute %.3f\n",
              mean(of_prior$joint_minus_sd[estimable]), mean(abs(of_prior$joint_minus_sd[estimable])),
              max(abs(of_prior$joint_minus_sd[estimable]))))

  record_check(sprintf("P3. %s: Savage-Dickey, logspline and kernel agree (estimable ones)", prior),
               max(abs(of_prior$savage_dickey - of_prior$savage_dickey_kernel)[estimable]) < 0.25,
               sprintf("%d of %d estimable; largest difference %.3f",
                       sum(estimable), nrow(of_prior),
                       max(abs(of_prior$savage_dickey - of_prior$savage_dickey_kernel)[estimable])))
  record_check(sprintf("P4. %s: conditional within a factor 2 of Savage-Dickey (estimable)", prior),
               max(abs(of_prior$conditional_minus_sd[estimable])) < log(2),
               sprintf("largest |difference| %.3f on the log scale (factor %.2f); mean |difference| %.3f",
                       max(abs(of_prior$conditional_minus_sd[estimable])),
                       exp(max(abs(of_prior$conditional_minus_sd[estimable]))),
                       mean(abs(of_prior$conditional_minus_sd[estimable]))))
  record_check(sprintf("P5. %s: the others, conditional within 3 Monte Carlo SE of Savage-Dickey", prior),
               all(abs(of_prior$conditional_minus_sd[!estimable]) <
                     3 * of_prior$savage_dickey_se[!estimable]),
               paste(sprintf("%.2f vs %.2f (SE %.2f)", of_prior$conditional[!estimable],
                             of_prior$savage_dickey[!estimable], of_prior$savage_dickey_se[!estimable]),
                     collapse = "; "))
}

save_figure(
  ggplot(coherence |> pivot_longer(c(conditional, joint), names_to = "test", values_to = "log_BF10"),
         aes(savage_dickey, log_BF10, colour = test, shape = factor(sample_size))) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = savage_dickey - 2 * savage_dickey_se,
                       xmax = savage_dickey + 2 * savage_dickey_se), height = 0, alpha = 0.5) +
    geom_point(size = 2.5) +
    facet_grid(prior ~ number_of_predictors,
               labeller = labeller(number_of_predictors = function(v) paste("k =", v))) +
    scale_colour_manual(values = c(conditional = "#0072B2", joint = "#56B4E9")) +
    labs(x = "Savage-Dickey log BF10 under the joint model", y = "log BF10",
         colour = NULL, shape = "n",
         title = "The conditional and the joint test against one prior for everything"),
  "10_coherence", width = 10, height = 8)

report_checks()
