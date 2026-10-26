## =====================================================================
##  Validation 06: does the picture hold beyond one design?
##
##  Study 03/04 used one design (4 predictors, one effect pattern). Here:
##    k (predictors)          2, 5, 10, 20
##    n                       50, 200, 1000
##    predictor correlation   0, 0.3, 0.7 (all pairs equal)
##    the other coefficients  sparse: one other = 0.3, rest 0
##                            dense:  all others 0.1
##    tested coefficient x1   0 (null), 0.1, 0.2
##  200 datasets per cell, noise SD 1 (effects in noise units).
##
##  Methods: the common set of validation/comparison_methods.R, all testing x1:
##    JZS default, Zellner-Siow, flat-R joint, uniform correlations joint
##    (model with vs without x1), flat-R conditional and uniform correlations
##    conditional (x1's partial correlation only; the latter is a flat prior on
##    the partial correlation).
##
##  Per cell: discrimination (AUC null vs 0.1 and vs 0.2, threshold-free) and
##  calibration at 3 (nulls found, false effects, 0.1 called null, 0.2 found).
##
##  Run from the project folder: Rscript validation/06_broad_simulation.R
##  About 20-30 minutes on 7 cores; cached.
## =====================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
source("validation/settings.R")             # parallel_map(), caching, save_figure()
pkgload::load_all(".", quiet = TRUE)
source("validation/comparison_methods.R")   # the common method set

## ---- Design ----------------------------------------------------------------
design <- expand_grid(
  number_of_predictors  = c(2, 5, 10, 20),
  sample_size           = c(50, 200, 1000),
  predictor_correlation = c(0, 0.3, 0.7),
  other_coefficients    = c("sparse", "dense"),
  slope_of_x1           = c(0, 0.1, 0.2))
replications_per_cell <- 200

slopes_for <- function(k, slope_of_x1, pattern) {
  others <- if (pattern == "sparse") c(0.3, rep(0, k - 2)) else rep(0.1, k - 1)
  c(slope_of_x1, others)
}

one_dataset <- function(k, n, correlation, pattern, slope_of_x1) {
  covariance <- matrix(correlation, k, k); diag(covariance) <- 1
  predictors <- matrix(rnorm(n * k), n) %*% chol(covariance)
  outcome    <- drop(predictors %*% slopes_for(k, slope_of_x1, pattern)) + rnorm(n)
  full_fit   <- lm(outcome ~ predictors)
  log_bayes_factors_one_coefficient(full_fit, 1, outcome, predictors)
}

runs <- run_simulation_grid(
  design,
  results_file = file.path(results_folder, "sim_broad.rds"),
  settings = list(replications_per_cell = replications_per_cell,
                  methods = one_coefficient_methods, version = 2),
  simulate_cell = function(cell)
    do.call(rbind, lapply(seq_len(replications_per_cell), function(replication)
      cbind(as.data.frame(cell), replication = replication,
            one_dataset(cell$number_of_predictors, cell$sample_size, cell$predictor_correlation,
                        cell$other_coefficients, cell$slope_of_x1)))))


runs <- mutate(runs, method = factor(method, levels = one_coefficient_methods))

## ---- Summaries ---------------------------------------------------------------
area_under_curve <- function(higher, lower)
  mean(outer(higher, lower, ">")) + 0.5 * mean(outer(higher, lower, "=="))

cell_summaries <- runs |>
  group_by(number_of_predictors, sample_size, predictor_correlation, other_coefficients, method) |>
  summarise(
    auc_null_vs_0.1   = area_under_curve(-log_BF10[slope_of_x1 == 0], -log_BF10[slope_of_x1 == 0.1]),
    auc_null_vs_0.2   = area_under_curve(-log_BF10[slope_of_x1 == 0], -log_BF10[slope_of_x1 == 0.2]),
    null_found        = mean(log_BF10[slope_of_x1 == 0]   < -log(3)),
    false_effect      = mean(log_BF10[slope_of_x1 == 0]   >  log(3)),
    small_called_null = mean(log_BF10[slope_of_x1 == 0.1] < -log(3)),
    effect_0.2_found  = mean(log_BF10[slope_of_x1 == 0.2] >  log(3)),
    .groups = "drop")
saveRDS(cell_summaries, file.path(results_folder, "sim_broad_summaries.rds"))

## Each method against JZS, averaged over the cells with the same k and
## correlation (positive = better than JZS)
versus_jzs <- cell_summaries |>
  group_by(number_of_predictors, sample_size, predictor_correlation, other_coefficients) |>
  mutate(auc_minus_jzs = auc_null_vs_0.1 - auc_null_vs_0.1[method == "JZS default"]) |>
  ungroup() |>
  filter(method != "JZS default") |>
  group_by(method, number_of_predictors, predictor_correlation) |>
  summarise(mean_auc_minus_jzs = mean(auc_minus_jzs),
            worst = min(auc_minus_jzs), best = max(auc_minus_jzs), .groups = "drop")
cat("\nAUC (null vs 0.1) minus JZS, averaged over n and the pattern of the other coefficients:\n")
print(as.data.frame(versus_jzs), digits = 3)

cat("\nCalibration at 3, averaged over everything except method and k:\n")
print(as.data.frame(cell_summaries |>
  group_by(method, number_of_predictors) |>
  summarise(across(c(null_found, false_effect, small_called_null, effect_0.2_found), mean),
            .groups = "drop")), digits = 3)


## ---- Figures -------------------------------------------------------------------
method_colours <- method_colours[one_coefficient_methods]
by_k_and_correlation <- facet_grid(
  predictor_correlation ~ number_of_predictors,
  labeller = labeller(number_of_predictors  = function(v) paste("k =", v),
                      predictor_correlation = function(v) paste("correlation =", v)))

save_figure(
  ggplot(filter(cell_summaries, other_coefficients == "sparse"),
         aes(sample_size, auc_null_vs_0.1, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    scale_x_log10(breaks = c(50, 200, 1000)) + scale_colour_manual(values = method_colours) +
    by_k_and_correlation + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "AUC: null vs slope 0.1", colour = NULL,
         title = "Discrimination across designs (other coefficients sparse)"),
  "06_auc_sparse", width = 12, height = 8)

save_figure(
  ggplot(filter(cell_summaries, other_coefficients == "dense"),
         aes(sample_size, auc_null_vs_0.1, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    scale_x_log10(breaks = c(50, 200, 1000)) + scale_colour_manual(values = method_colours) +
    by_k_and_correlation + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "AUC: null vs slope 0.1", colour = NULL,
         title = "Discrimination across designs (other coefficients dense, all 0.1)"),
  "06_auc_dense", width = 12, height = 8)

save_figure(
  ggplot(cell_summaries |>
           pivot_longer(c(null_found, small_called_null), names_to = "rate", values_to = "value") |>
           mutate(rate = recode(rate, null_found = "true null found (BF01 > 3)",
                                small_called_null = "slope 0.1 called null (BF01 > 3)")) |>
           group_by(method, number_of_predictors, sample_size, rate) |>
           summarise(value = mean(value), .groups = "drop"),
         aes(sample_size, value, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    scale_x_log10(breaks = c(50, 200, 1000)) + scale_colour_manual(values = method_colours) +
    facet_grid(rate ~ number_of_predictors,
               labeller = labeller(number_of_predictors = function(v) paste("k =", v))) +
    theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "proportion of datasets", colour = NULL,
         title = "Calibration at 3, averaged over correlation and pattern"),
  "06_calibration", width = 12, height = 7)
