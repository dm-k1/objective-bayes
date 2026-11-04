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
##  Methods: the common set of validation/comparison_methods.R, all testing x1.
##  The three contenders, flat partial r, flat R and flat R-squared, each two
##  ways: joint (model with vs without x1) and conditional (x1's partial
##  correlation only; for flat partial r a flat prior on it). The references,
##  JZS default and Zellner-Siow, model with vs without x1.
##
##  Per cell: discrimination (AUC null vs 0.1 and vs 0.2, threshold-free) and
##  calibration at 3 (nulls found, false effects, 0.1 called null, 0.2 found).
##  And every outcome at 3 (says null, undecided, says effect) for the
##  conditional tests and JZS at k = 2 and 20: with equal AUC, a test that
##  calls fewer effects null may simply decide less often.
##
##  Run from the project folder: Rscript validation/06_broad_simulation.R
##  About 20-30 minutes on 7 cores (eight methods); cached.
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
  results_file = file.path(cache_folder, "sim_broad.rds"),
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
saveRDS(cell_summaries, file.path(cache_folder, "sim_broad_summaries.rds"))

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

## Every outcome, not only the correct one: at 3, a test says null, says
## effect, or stays undecided. Shown for the conditional tests and JZS at
## k = 2 and 20. With the same AUC, a test that "calls fewer effects null"
## may simply say less: the difference is where it decides, not how well.
cat("\nWhat each test says at 3 (share of datasets), k = 2 and 20:\n")
print(as.data.frame(runs |>
  filter(number_of_predictors %in% c(2, 20), grepl("conditional|JZS", method)) |>
  mutate(outcome = case_when(log_BF10 < -log(3) ~ "says null", log_BF10 > log(3) ~ "says effect",
                             TRUE ~ "undecided")) |>
  count(method, number_of_predictors, slope_of_x1, outcome) |>
  group_by(method, number_of_predictors, slope_of_x1) |> mutate(share = n / sum(n)) |> ungroup() |>
  select(-n) |>
  pivot_wider(names_from = outcome, values_from = share, values_fill = 0) |>
  select(method, number_of_predictors, slope_of_x1, `says null`, undecided, `says effect`)),
  digits = 3)


## ---- Figures -------------------------------------------------------------------
by_k_and_correlation <- facet_grid(
  predictor_correlation ~ number_of_predictors,
  labeller = labeller(number_of_predictors  = function(v) paste("k =", v),
                      predictor_correlation = function(v) paste("correlation =", v)))

## AUC minus JZS per cell: the AUC curves themselves lie on top of each other
auc_minus_jzs_per_cell <- cell_summaries |>
  group_by(number_of_predictors, sample_size, predictor_correlation, other_coefficients) |>
  mutate(auc_minus_jzs = auc_null_vs_0.1 - auc_null_vs_0.1[method == "JZS default"]) |>
  ungroup() |> filter(method != "JZS default")
for (pattern in c("sparse", "dense"))
  save_figure(
    ggplot(filter(auc_minus_jzs_per_cell, other_coefficients == pattern),
           aes(sample_size, auc_minus_jzs, colour = method, linetype = method)) +
      geom_hline(yintercept = 0, colour = "grey40") +
      geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
      scale_x_log10(breaks = c(50, 200, 1000)) +
      method_scales(setdiff(one_coefficient_methods, "JZS default")) +
      by_k_and_correlation +
      labs(x = "sample size (log scale)", y = "AUC (null vs slope 0.1) minus JZS",
           title = sprintf("Discrimination against JZS (other coefficients %s)",
                           if (pattern == "sparse") "sparse" else "dense, all 0.1"),
           subtitle = "above 0: separates nulls from small effects better than JZS; the three conditional tests coincide"),
    paste0("06_auc_", pattern), width = 12, height = 8)

## What each test says, as one bar per test: the share that says null, stays
## undecided or says effect (at 3), by k and the true slope of x1, averaged
## over n, correlation and pattern. Equal AUC, different bars: the tests
## differ in where they decide, not in how well they separate.
outcomes <- runs |>
  filter(grepl("conditional|JZS", method)) |>
  mutate(test = factor(sub(" conditional", "", method),
                       levels = c("flat partial r", "flat R", "flat R-squared", "JZS default")),
         outcome = factor(case_when(log_BF10 < -log(3) ~ "says null", log_BF10 > log(3) ~ "says effect",
                                    TRUE ~ "undecided"),
                          levels = c("says effect", "undecided", "says null")),
         truth = factor(paste("true slope", slope_of_x1),
                        levels = paste("true slope", c(0, 0.1, 0.2)))) |>
  count(test, number_of_predictors, truth, outcome) |>
  group_by(test, number_of_predictors, truth) |> mutate(share = n / sum(n)) |> ungroup()
save_figure(
  ggplot(outcomes, aes(test, share, fill = outcome)) +
    geom_col(width = 0.75) +
    geom_text(data = filter(outcomes, share >= 0.08),
              aes(label = sprintf("%.2f", share)), position = position_stack(vjust = 0.5),
              size = 2.7, colour = "white") +
    scale_fill_manual(values = c(`says null` = "#542788", undecided = "grey70", `says effect` = "#B35806"),
                      breaks = c("says null", "undecided", "says effect"), name = "at BF 3:") +
    facet_grid(truth ~ number_of_predictors,
               labeller = labeller(number_of_predictors = function(v) paste("k =", v))) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") +
    labs(x = NULL, y = "share of datasets",
         title = "What each conditional test (and JZS) says about x1, as k grows",
         subtitle = "averaged over n, correlation and pattern; correct: says null in the top row, says effect below"),
  "06_outcomes", width = 11, height = 8)
