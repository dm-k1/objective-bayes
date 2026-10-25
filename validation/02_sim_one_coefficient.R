## =====================================================================
##  Validation 02: one coefficient at a time
##  JZS default vs a flat prior on the partial correlation (the
##  one-predictor case of the flat-R idea), two predictors
##  (x and a covariate z).
##
##   A. simulation   true slope x sample size x predictor correlation
##   B. summaries    misleading / correct evidence, coverage, bias
##   C. confusion and calibration in different "worlds" of true effects
##   D. precision matching: which null threshold makes the flat partial r null
##      calls as reliable as JZS's?
##
##  Run from the project folder. Section A runs in parallel and is cached:
##  it only reruns when its settings change (or rerun_simulation = TRUE).
##  Figures are saved to validation/results/figures/.
## =====================================================================
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

source("R/flat_r.R")
source("validation/comparison_methods.R")   # JZS and flat partial correlation
source("validation/settings.R")             # parallel cores, cached simulation runner, save_figure()

results_file     <- file.path(results_folder, "sim_one_coefficient.rds")
rerun_simulation <- FALSE   # TRUE forces a rerun even if cached results match


## #####################################################################
##  A. SIMULATION
## #####################################################################

## ---- Settings ----------------------------------------------------------
true_slopes            <- c(0, 0.1, 0.2, 0.3, 0.5)   # slope of x (population SDs all 1)
sample_sizes           <- c(20, 50, 100, 200, 500)
predictor_correlations <- c(0, 0.5)                  # correlation between x and covariate z
replications_per_cell  <- 200
covariate_slope        <- 0.5                        # slope of z
jzs_scale              <- sqrt(2) / 4                # BayesFactor's "medium" for regression

slope_grid <- seq(-2, 2, length.out = 4001)          # for the JZS posterior

## True partial correlation of x and y given z
true_partial_correlation <- function(true_slope, predictor_correlation) {
  noise_variance <- 1 - (true_slope^2 + covariate_slope^2 +
                         2 * true_slope * covariate_slope * predictor_correlation)
  true_slope * sqrt(1 - predictor_correlation^2) /
    sqrt(true_slope^2 * (1 - predictor_correlation^2) + noise_variance)
}

## JZS posterior for the slope of x, approximated by a Cauchy prior on the
## slope (the Bayes factor itself is exact, from BayesFactor)
posterior_jzs_cauchy_approximation <- function(estimate, standard_error, degrees_of_freedom) {
  posterior  <- dt((estimate - slope_grid) / standard_error, degrees_of_freedom) / standard_error *
                dcauchy(slope_grid, 0, jzs_scale)
  cumulative <- cumsum(posterior) / sum(posterior)
  list(slope_mean  = sum(slope_grid * posterior) / sum(posterior),
       slope_lower = slope_grid[which.max(cumulative >= 0.025)],
       slope_upper = slope_grid[which.max(cumulative >= 0.975)])
}

## ---- One dataset, all methods ------------------------------------------
simulate_one_dataset <- function(true_slope, sample_size, predictor_correlation) {
  z <- rnorm(sample_size)
  x <- predictor_correlation * z + sqrt(1 - predictor_correlation^2) * rnorm(sample_size)
  noise_sd <- sqrt(1 - (true_slope^2 + covariate_slope^2 +
                        2 * true_slope * covariate_slope * predictor_correlation))   # Var(y) = 1
  y <- true_slope * x + covariate_slope * z + rnorm(sample_size, sd = noise_sd)

  fit             <- lm(y ~ x + z)
  coefficient_row <- summary(fit)$coefficients["x", ]
  estimate        <- coefficient_row[["Estimate"]]
  standard_error  <- coefficient_row[["Std. Error"]]
  degrees_of_freedom <- fit$df.residual
  r_squared_full     <- summary(fit)$r.squared
  r_squared_without  <- summary(lm(y ~ z))$r.squared
  t_critical         <- qt(0.975, degrees_of_freedom)

  ## 1. JZS default
  jzs_posterior <- posterior_jzs_cauchy_approximation(estimate, standard_error, degrees_of_freedom)
  jzs <- tibble(method          = "JZS default",
                BF10            = bayes_factor_jzs_one_coefficient(r_squared_full, r_squared_without,
                                                                   sample_size, 2, rscale = jzs_scale),
                posterior_mean  = jzs_posterior$slope_mean,
                interval_lower  = jzs_posterior$slope_lower,
                interval_upper  = jzs_posterior$slope_upper)

  ## 2. Flat prior on the partial correlation, from -1 to 1
  predictors     <- cbind(x, z)
  flat_posterior <- posterior_flat_partial_correlation(fit, 1, predictors)
  flat <- tibble(method            = "flat partial correlation",
                 BF10              = bayes_factor_flat_partial_correlation(fit, 1, predictors),
                 posterior_mean    = flat_posterior$slope_mean,
                 interval_lower    = flat_posterior$slope_lower,
                 interval_upper    = flat_posterior$slope_upper,
                 correlation_mean  = flat_posterior$correlation_mean,
                 correlation_lower = flat_posterior$correlation_lower,
                 correlation_upper = flat_posterior$correlation_upper)

  ## 3. Frequentist reference
  frequentist <- tibble(method         = "frequentist",
                        BF10           = NA_real_,
                        posterior_mean = estimate,
                        interval_lower = estimate - t_critical * standard_error,
                        interval_upper = estimate + t_critical * standard_error)

  bind_rows(jzs, flat, frequentist) |>
    mutate(p_value = coefficient_row[["Pr(>|t|)"]])
}

## ---- Run the grid (parallel; cached until the settings change) ----------
simulation_design <- expand_grid(true_slope            = true_slopes,
                                 sample_size           = sample_sizes,
                                 predictor_correlation = predictor_correlations)

simulation_runs <- run_simulation_grid(
  simulation_design, results_file = results_file, rerun = rerun_simulation,
  settings = list(replications_per_cell = replications_per_cell,
                  covariate_slope = covariate_slope, jzs_scale = jzs_scale),
  simulate_cell = function(cell) {
    map_dfr(seq_len(replications_per_cell), function(replication)
      mutate(simulate_one_dataset(cell$true_slope, cell$sample_size, cell$predictor_correlation),
             replication = replication)) |>
      mutate(true_slope = cell$true_slope, sample_size = cell$sample_size,
             predictor_correlation = cell$predictor_correlation, .before = 1)
  }) |>
  mutate(true_partial_correlation = true_partial_correlation(true_slope, predictor_correlation))


## #####################################################################
##  B. SUMMARIES
## #####################################################################

cell_summaries <- simulation_runs |>
  group_by(true_slope, sample_size, predictor_correlation, method) |>
  summarise(
    bf10_above_3  = mean(BF10 > 3),
    bf10_above_10 = mean(BF10 > 10),
    bf01_above_3  = mean(BF10 < 1/3),
    bf01_above_10 = mean(BF10 < 1/10),
    p_below_05    = mean(p_value < .05),
    coverage      = mean(interval_lower <= true_slope & true_slope <= interval_upper),
    correlation_coverage = mean(correlation_lower <= true_partial_correlation &   # flat method only
                                true_partial_correlation <= correlation_upper),
    interval_width = mean(interval_upper - interval_lower),
    bias           = mean(posterior_mean - true_slope),
    rmse           = sqrt(mean((posterior_mean - true_slope)^2)),
    .groups = "drop") |>
  mutate(
    misleading_3  = ifelse(true_slope == 0, bf10_above_3,  bf01_above_3),
    misleading_10 = ifelse(true_slope == 0, bf10_above_10, bf01_above_10),
    correct_3     = ifelse(true_slope == 0, bf01_above_3,  bf10_above_3),
    correct_10    = ifelse(true_slope == 0, bf01_above_10, bf10_above_10))

## Large-sample check: at n = 500 the flat-prior interval should match the
## frequentist one
cell_summaries |>
  filter(sample_size == 500, method != "JZS default") |>
  select(predictor_correlation, true_slope, method, coverage, correlation_coverage, interval_width) |>
  arrange(predictor_correlation, true_slope, method) |>
  print(n = Inf)

## Testing
cell_summaries |>
  filter(method != "frequentist") |>
  select(predictor_correlation, true_slope, sample_size, method,
         misleading_3, misleading_10, correct_3, correct_10) |>
  arrange(predictor_correlation, true_slope, sample_size, method) |>
  print(n = Inf)

## Estimation
cell_summaries |>
  select(predictor_correlation, true_slope, sample_size, method,
         coverage, correlation_coverage, interval_width, bias, rmse, p_below_05) |>
  arrange(predictor_correlation, true_slope, sample_size, method) |>
  print(n = Inf)

## ---- Figures -------------------------------------------------------------
method_colours <- c(`JZS default` = "grey40", `flat partial correlation` = "#1b9e77",
                    frequentist = "black")
bayesian_summaries <- filter(cell_summaries, method != "frequentist")
facets_by_cell     <- facet_grid(predictor_correlation ~ true_slope, labeller = label_both)
sample_size_axis   <- scale_x_log10(breaks = sample_sizes)

save_figure(
  ggplot(bayesian_summaries, aes(sample_size, misleading_3, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() + sample_size_axis + facets_by_cell +
    scale_colour_manual(values = method_colours) +
    labs(y = "P(BF > 3 for the wrong hypothesis)", title = "Misleading evidence"),
  "02_misleading_evidence")

save_figure(
  ggplot(bayesian_summaries, aes(sample_size, correct_3, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    geom_line(data = filter(cell_summaries, method == "frequentist"), aes(y = p_below_05),
              colour = "black", linetype = "dashed") +
    sample_size_axis + facets_by_cell + scale_colour_manual(values = method_colours) +
    labs(y = "P(BF > 3 for the true hypothesis)", title = "Correct evidence",
         subtitle = "dashed = P(p < .05)"),
  "02_correct_evidence")

save_figure(
  ggplot(cell_summaries, aes(sample_size, coverage, colour = method)) +
    geom_hline(yintercept = 0.95, linetype = "dotted") +
    geom_line(linewidth = 0.9) + geom_point() + sample_size_axis + facets_by_cell +
    scale_colour_manual(values = method_colours) +
    labs(y = "coverage of the true slope", title = "Coverage of 95% intervals"),
  "02_coverage")

save_figure(
  ggplot(cell_summaries, aes(sample_size, bias, colour = method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 0.9) + geom_point() + sample_size_axis + facets_by_cell +
    scale_colour_manual(values = method_colours) +
    labs(y = "posterior mean - true slope", title = "Bias"),
  "02_bias")


## #####################################################################
##  C. CONFUSION AND CALIBRATION
## #####################################################################

## How often each true slope occurs in the world you study (rows sum to 1)
effect_worlds <- tribble(
  ~world,                 ~`0`, ~`0.1`, ~`0.2`, ~`0.3`, ~`0.5`,
  "balanced",              0.5,  0.125,  0.125,  0.125,  0.125,
  "small effects common",  0.5,  0.25,   0.15,   0.07,   0.03,
  "large effects common",  0.5,  0.03,   0.07,   0.15,   0.25,
  "few true effects",      0.8,  0.05,   0.05,   0.05,   0.05) |>
  pivot_longer(-world, names_to = "true_slope", values_to = "world_weight") |>
  mutate(true_slope = as.numeric(true_slope))

bayesian_runs <- simulation_runs |>
  filter(method != "frequentist") |>
  mutate(truth   = ifelse(true_slope == 0, "null", "effect"),
         verdict = case_when(BF10 > 3   ~ "effect",
                             BF10 < 1/3 ~ "null",
                             TRUE       ~ "undecided")) |>
  group_by(true_slope, sample_size, predictor_correlation, method) |>
  mutate(run_weight = 1 / n()) |>          # each cell sums to 1
  ungroup()

weighted_runs <- bayesian_runs |>
  inner_join(effect_worlds, by = "true_slope", relationship = "many-to-many") |>
  mutate(weight = world_weight * run_weight)

## ---- 1. Confusion table per world, method, n ------------------------------
confusion <- weighted_runs |>
  group_by(world, predictor_correlation, sample_size, method) |>
  summarise(
    says_null        = sum(weight[verdict == "null"]),
    says_effect      = sum(weight[verdict == "effect"]),
    null_precision   = sum(weight[verdict == "null"   & truth == "null"])   / says_null,
    effect_precision = sum(weight[verdict == "effect" & truth == "effect"]) / says_effect,
    decisive         = says_null + says_effect,
    accuracy_decided = sum(weight[verdict == truth]) / decisive,
    wrong_null       = sum(weight[verdict == "null"   & truth == "effect"]),   # missed effects
    wrong_effect     = sum(weight[verdict == "effect" & truth == "null"]),     # false alarms
    .groups = "drop") |>
  mutate(null_lean = log(wrong_null / pmax(wrong_effect, 1e-4)))   # > 0: errors lean to null

confusion |>
  select(world, predictor_correlation, sample_size, method, decisive, null_precision,
         effect_precision, accuracy_decided, wrong_null, wrong_effect) |>
  arrange(world, predictor_correlation, sample_size, method) |>
  print(n = Inf)

## ---- 2. Calibration: when a method says BF = k, are the odds k:1? --------
## With 50% nulls the prior odds are 1, so observed odds of effect vs
## null should equal the stated BF if the method is calibrated for that world.
calibration <- weighted_runs |>
  filter(world != "few true effects") |>
  mutate(log_bf10 = log(BF10),
         bin      = cut(log_bf10, breaks = c(-Inf, seq(-4, 4, by = 1), Inf))) |>
  group_by(world, method, bin) |>
  summarise(stated   = weighted.mean(pmin(pmax(log_bf10, -5), 5), weight),
            observed = log(sum(weight[truth == "effect"]) / sum(weight[truth == "null"])),
            mass     = sum(weight), .groups = "drop") |>
  filter(is.finite(observed), mass > 0.002)

save_figure(
  ggplot(calibration, aes(stated, observed, colour = method)) +
    geom_abline(linetype = "dashed") +
    geom_line(linewidth = 0.9) + geom_point(aes(size = mass)) +
    scale_colour_manual(values = method_colours) +
    facet_wrap(~ world) +
    labs(x = "stated log BF10", y = "observed log odds (effect : null)",
         size = "share of runs",
         title = "Calibration: does the Bayes factor mean what it says?",
         subtitle = paste("on the dashed line = calibrated; above = understates evidence",
                          "for effects; below = overstates it")),
  "02_calibration")


## #####################################################################
##  D. PRECISION MATCHING
## #####################################################################

null_thresholds <- 10^seq(0, 2, by = 0.01)   # candidate BF01 thresholds, 1 to 100

## Null precision and number of calls for one method at a BF01 threshold
## (effect calls stay at BF10 > 3)
null_precision_at_threshold <- function(runs, threshold) {
  null_call   <- runs$BF10 < 1 / threshold
  effect_call <- runs$BF10 > 3
  tibble(threshold      = threshold,
         says_null      = sum(runs$weight[null_call]),
         null_precision = sum(runs$weight[null_call & runs$truth == "null"]) /
                          sum(runs$weight[null_call]),
         decisive       = sum(runs$weight[null_call | effect_call]))
}

precision_matched <- weighted_runs |>
  group_by(world, predictor_correlation, sample_size) |>
  group_modify(function(cell_runs, key) {
    jzs_runs  <- filter(cell_runs, method == "JZS default")
    flat_runs <- filter(cell_runs, method == "flat partial correlation")
    jzs_reference <- null_precision_at_threshold(jzs_runs, 3)
    no_match <- tibble(threshold = NA_real_, says_null = NA_real_,
                       null_precision = NA_real_, decisive = NA_real_)
    flat_match <- if (is.na(jzs_reference$null_precision)) no_match else {
      flat_curve   <- map_dfr(null_thresholds, ~ null_precision_at_threshold(flat_runs, .x))
      at_least_jzs <- filter(flat_curve, !is.na(null_precision),
                             null_precision >= jzs_reference$null_precision)
      if (nrow(at_least_jzs)) slice_min(at_least_jzs, threshold, n = 1, with_ties = FALSE)
      else no_match
    }
    tibble(jzs_null_precision  = jzs_reference$null_precision,
           jzs_says_null       = jzs_reference$says_null,
           jzs_decisive        = jzs_reference$decisive,
           threshold_needed    = flat_match$threshold,
           flat_null_precision = flat_match$null_precision,
           flat_says_null      = flat_match$says_null,
           flat_decisive       = flat_match$decisive)
  }) |>
  ungroup() |>
  mutate(extra_null_calls = flat_says_null - jzs_says_null,
         extra_decisive   = flat_decisive  - jzs_decisive)

precision_matched |>
  mutate(across(where(is.double), ~ round(.x, 3))) |>
  select(world, predictor_correlation, sample_size, threshold_needed,
         jzs_null_precision, flat_null_precision, jzs_says_null, flat_says_null,
         extra_null_calls, extra_decisive) |>
  arrange(world, predictor_correlation, sample_size) |>
  print(n = Inf)

save_figure(
  ggplot(precision_matched, aes(sample_size, extra_null_calls, colour = factor(predictor_correlation))) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 0.9) + geom_point() +
    sample_size_axis + facet_wrap(~ world) +
    labs(x = "n", y = "extra null calls (share of all studies)",
         colour = "predictor correlation",
         title = "Flat-prior null calls, made as reliable as JZS's",
         subtitle = "above 0: same reliability, more conclusions; below 0: fewer"),
  "02_precision_matched_extra_null_calls")

save_figure(
  ggplot(precision_matched, aes(sample_size, threshold_needed, colour = factor(predictor_correlation))) +
    geom_hline(yintercept = 3, linetype = "dotted") +
    geom_line(linewidth = 0.9) + geom_point() +
    sample_size_axis + scale_y_log10() + facet_wrap(~ world) +
    labs(x = "n", y = "BF01 threshold needed", colour = "predictor correlation",
         title = "Null threshold that matches JZS's null reliability",
         subtitle = "dotted = the usual 3"),
  "02_precision_matched_threshold")
