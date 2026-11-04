## =====================================================================
##  Validation 11: testing a factor (a term with several columns)
##
##  bayes_factor_per_term() tests a term of q columns (a factor with q + 1
##  levels) with the block conditional test: only the block's partial R-squared
##  given all the other columns, with the prior the joint prior implies for it
##  (flat partial r: partial R-squared ~ Beta(q/2, 1); flat R and flat
##  R-squared: share Beta(q/2, (k - q)/2) of a total with R-squared ~ Beta(1/2, 1)
##  or Beta(1, 1)). Study 05 (C9-C11) checks that it is computed right; this
##  study compares it, for the three contenders, with their joint tests (the
##  model with against without the factor), and with two references: JZS as
##  BayesFactor computes it for a factor (lmBF: a separate g for the factor,
##  scale 1/2, and one for the continuous covariates, scale sqrt(2)/4) and
##  Zellner-Siow (one g for all columns, scale 1; joint).
##
##  Design:
##    q (columns of the factor)       2, 3, 4, 5 (3 to 6 equal-sized levels)
##    n                               50, 200, 1000
##    other predictors (continuous)   2, 10, uncorrelated, all with equal slopes
##    what the others explain         R-squared 0, 0.5, 0.9 (of the outcome
##                                    when the factor is null)
##    factor effect, Cohen's f^2      0 (null), 0.02 (small), 0.15 (medium):
##                                    the variance of the level means in noise
##                                    units; equally spaced level means
##  200 datasets per cell. The three values of what the others explain use the
##  same predictors, factor and noise, so they differ only in the others'
##  slopes: the conditional tests must give exactly the same answer for all
##  three (they see the factor only after the others are projected out), and
##  any difference in the joint tests and JZS is the conditional Lindley
##  paradox (the evidence for a term depending on how much the others explain).
##
##  Per cell: discrimination (AUC null vs small, threshold-free), nulls
##  confirmed (BF01 > 3), small effects called null (BF01 > 3), false effects
##  and medium effects found (BF10 > 3).
##    F1. the conditional tests: identical for the three values of what the
##        others explain (each dataset), for each of the three contenders.
##    F2. the three conditional tests have the same AUC in every cell (all
##        rank datasets by the block's partial R-squared).
##    F3. JZS (lmBF) is computed precisely enough: its reported proportional
##        error is below 0.05 in every dataset (recomputed with 10^5 and then
##        10^6 draws where it is above 0.02).
##
##  Run from the project folder: Rscript validation/11_factor_terms.R
##  About an hour on 7 cores (lmBF dominates); cached.
## =====================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
source("validation/settings.R")             # parallel_map(), caching, record_check()
library(flatcor)   # the package under test, at the version in renv.lock
source("validation/comparison_methods.R")   # log_bayes_factor_model(), method_scales()

## ---- Design ----------------------------------------------------------------
design <- expand_grid(
  columns_of_factor = 2:5,
  sample_size       = c(50, 200, 1000),
  other_predictors  = c(2, 10),
  factor_effect     = c("null", "small", "medium"))
replications_per_cell <- 200
explained_by_others   <- c(0, 0.5, 0.9)
cohen_f_squared       <- c(null = 0, small = 0.02, medium = 0.15)

## The three contenders, conditional (block) and joint, and the two references;
## colours and line types from method_scales()
factor_methods <- c(paste(names(contenders), "conditional (block)"),
                    paste(names(contenders), "joint"), "JZS (lmBF)", "Zellner-Siow")

## Level means: equally spaced, centred, with mean square f^2 (equal-sized
## levels), so the factor's variance is f^2 times the noise variance
level_means <- function(q, f_squared) {
  spacing <- seq_len(q + 1) - (q + 2) / 2
  spacing * sqrt(f_squared / mean(spacing^2))
}

## All methods for one dataset, for each value of what the others explain
one_dataset <- function(q, n, other_predictors, f_squared) {
  group  <- factor(sample(rep(seq_len(q + 1), length.out = n)))
  others <- matrix(rnorm(n * other_predictors), n,
                   dimnames = list(NULL, paste0("x", seq_len(other_predictors))))
  noise  <- rnorm(n)
  k <- q + other_predictors
  design_full    <- model.matrix(~ group + others)
  design_without <- cbind(1, others)
  formula_full    <- as.formula(paste("y ~ group +", paste(colnames(others), collapse = " + ")))
  formula_without <- as.formula(paste("y ~", paste(colnames(others), collapse = " + ")))

  do.call(rbind, lapply(explained_by_others, function(r_squared_others) {
    slopes  <- rep(sqrt(r_squared_others / (1 - r_squared_others) / other_predictors), other_predictors)
    outcome <- drop(others %*% slopes) + level_means(q, f_squared)[group] + noise
    total_ss         <- sum((outcome - mean(outcome))^2)
    residual_full    <- sum(qr.resid(qr(design_full), outcome)^2)
    residual_without <- sum(qr.resid(qr(design_without), outcome)^2)
    r_squared_full    <- 1 - residual_full / total_ss
    r_squared_without <- 1 - residual_without / total_ss
    block_partial_r   <- sqrt(1 - residual_full / residual_without)
    joint <- function(model_method)
      log_bayes_factor_model(model_method, r_squared_full, n, k) -
      log_bayes_factor_model(model_method, r_squared_without, n, k - q)
    data_frame  <- data.frame(y = outcome, group = group, others)
    jzs_full    <- BayesFactor::lmBF(formula_full,    data_frame, progress = FALSE)
    ## two g's (factor, covariates): Monte Carlo; more draws where it is imprecise
    for (iterations in c(1e5, 1e6))
      if (BayesFactor::extractBF(jzs_full)$error > 0.02)
        jzs_full <- BayesFactor::recompute(jzs_full, iterations = iterations, progress = FALSE)
    jzs_without <- BayesFactor::lmBF(formula_without, data_frame, progress = FALSE)
    data.frame(
      explained_by_others = r_squared_others,
      block_partial_r_squared = block_partial_r^2,
      method = factor_methods,
      log_BF10 = c(
        vapply(contenders, function(prior)
          log_bayes_factor_conditional(block_partial_r, n, k, prior, number_of_columns = q), numeric(1)),
        vapply(names(contenders), joint, numeric(1)),
        BayesFactor::extractBF(jzs_full, logbf = TRUE)$bf -
          BayesFactor::extractBF(jzs_without, logbf = TRUE)$bf,
        joint("Zellner-Siow")),
      jzs_error = ifelse(factor_methods == "JZS (lmBF)", BayesFactor::extractBF(jzs_full)$error, NA))
  }))
}

runs <- run_simulation_grid(
  design,
  results_file = file.path(cache_folder, "sim_factor_terms.rds"),
  settings = list(replications_per_cell = replications_per_cell, explained_by_others = explained_by_others,
                  cohen_f_squared = cohen_f_squared, methods = factor_methods, version = 2),
  simulate_cell = function(cell)
    do.call(rbind, lapply(seq_len(replications_per_cell), function(replication)
      cbind(as.data.frame(cell), replication = replication,
            one_dataset(cell$columns_of_factor, cell$sample_size, cell$other_predictors,
                        cohen_f_squared[[cell$factor_effect]])))))

runs <- mutate(runs, method = factor(method, levels = factor_methods))


## ---- F1-F3: checks -----------------------------------------------------------
conditional_spread <- runs |>
  filter(grepl("conditional", method)) |>
  group_by(columns_of_factor, sample_size, other_predictors, factor_effect, replication, method) |>
  summarise(spread = diff(range(log_BF10)), .groups = "drop")
record_check("F1 conditional: same for every R2 of the others",
             max(conditional_spread$spread) < 1e-8,
             sprintf("largest difference in log BF10 %.1e (%d datasets, three contenders)",
                     max(conditional_spread$spread), nrow(conditional_spread) / length(contenders)))

area_under_curve <- function(higher, lower)
  mean(outer(higher, lower, ">")) + 0.5 * mean(outer(higher, lower, "=="))

cell_summaries <- runs |>
  group_by(columns_of_factor, sample_size, other_predictors, explained_by_others, method) |>
  summarise(
    auc_null_vs_small  = area_under_curve(-log_BF10[factor_effect == "null"], -log_BF10[factor_effect == "small"]),
    auc_null_vs_medium = area_under_curve(-log_BF10[factor_effect == "null"], -log_BF10[factor_effect == "medium"]),
    null_found         = mean(log_BF10[factor_effect == "null"]   < -log(3)),
    false_effect       = mean(log_BF10[factor_effect == "null"]   >  log(3)),
    small_called_null  = mean(log_BF10[factor_effect == "small"]  < -log(3)),
    small_found        = mean(log_BF10[factor_effect == "small"]  >  log(3)),
    medium_found       = mean(log_BF10[factor_effect == "medium"] >  log(3)),
    .groups = "drop")
saveRDS(cell_summaries, file.path(cache_folder, "sim_factor_terms_summaries.rds"))

conditional_auc_difference <- cell_summaries |>
  filter(grepl("conditional", method)) |>
  group_by(columns_of_factor, sample_size, other_predictors, explained_by_others) |>
  summarise(difference = diff(range(auc_null_vs_small)), .groups = "drop")
record_check("F2 the three conditional tests: same AUC",
             max(conditional_auc_difference$difference) < 1e-12,
             sprintf("largest difference %.1e over %d cells", max(conditional_auc_difference$difference),
                     nrow(conditional_auc_difference)))

jzs_errors <- runs$jzs_error[runs$method == "JZS (lmBF)"]
record_check("F3 JZS (lmBF): proportional error < 0.05",
             max(jzs_errors, na.rm = TRUE) < 0.05,
             sprintf("median %.4f, 99%% %.4f, largest %.3f", median(jzs_errors, na.rm = TRUE),
                     quantile(jzs_errors, 0.99, na.rm = TRUE), max(jzs_errors, na.rm = TRUE)))


## ---- Summaries ---------------------------------------------------------------
print_table <- function(title, table) {
  cat("\n", title, "\n", sep = "")
  print(as.data.frame(table), digits = 3, row.names = FALSE)
}

print_table("AUC null vs small (f^2 = 0.02), averaged over q and the number of others; columns n_R2 of the others:",
  cell_summaries |>
    group_by(method, sample_size, explained_by_others) |>
    summarise(auc = mean(auc_null_vs_small), .groups = "drop") |>
    pivot_wider(names_from = c(sample_size, explained_by_others), values_from = auc,
                names_glue = "n{sample_size}_R2o{explained_by_others}"))

print_table("AUC null vs small minus JZS (lmBF), over the cells of each n: mean, worst, best:",
  cell_summaries |>
    group_by(columns_of_factor, sample_size, other_predictors, explained_by_others) |>
    mutate(auc_minus_jzs = auc_null_vs_small - auc_null_vs_small[method == "JZS (lmBF)"]) |>
    ungroup() |> filter(method != "JZS (lmBF)") |>
    group_by(method, sample_size) |>
    summarise(mean = mean(auc_minus_jzs), worst = min(auc_minus_jzs), best = max(auc_minus_jzs),
              .groups = "drop"))

for (rate in c("null_found", "small_called_null", "false_effect", "small_found", "medium_found")) {
  print_table(sprintf("%s, averaged over q, by n, the number of others and their R2:", rate),
    cell_summaries |>
      group_by(method, sample_size, other_predictors, explained_by_others) |>
      summarise(value = mean(.data[[rate]]), .groups = "drop") |>
      pivot_wider(names_from = explained_by_others, values_from = value, names_prefix = "R2o="))
}

print_table("null_found and small_called_null by q (averaged over n, the others and their R2):",
  cell_summaries |>
    group_by(method, columns_of_factor) |>
    summarise(null_found = mean(null_found), small_called_null = mean(small_called_null),
              .groups = "drop") |>
    pivot_wider(names_from = columns_of_factor, values_from = c(null_found, small_called_null),
                names_glue = "{.value}_q{columns_of_factor}"))

## The paradox in one number per method: the median log BF10 for a null factor
## as the others go from explaining 0 to 0.9 (n = 200)
print_table("Median log BF10 for a null factor (n = 200), by what the others explain:",
  runs |>
    filter(factor_effect == "null", sample_size == 200) |>
    group_by(method, other_predictors, explained_by_others) |>
    summarise(median_log_BF10 = median(log_BF10), .groups = "drop") |>
    pivot_wider(names_from = explained_by_others, values_from = median_log_BF10, names_prefix = "R2o="))


## ---- Figures -------------------------------------------------------------------
by_others <- facet_grid(other_predictors ~ explained_by_others,
  labeller = labeller(other_predictors    = function(v) paste(v, "other predictors"),
                      explained_by_others = function(v) paste("others explain R-squared =", v)))
averaged_over_q <- cell_summaries |>
  group_by(method, sample_size, other_predictors, explained_by_others) |>
  summarise(across(c(auc_null_vs_small, null_found, small_called_null), mean), .groups = "drop")

save_figure(
  ggplot(averaged_over_q |>
           group_by(sample_size, other_predictors, explained_by_others) |>
           mutate(auc_minus_jzs = auc_null_vs_small - auc_null_vs_small[method == "JZS (lmBF)"]) |>
           ungroup() |> filter(method != "JZS (lmBF)"),
         aes(sample_size, auc_minus_jzs, colour = method, linetype = method)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200, 1000)) + method_scales(setdiff(factor_methods, "JZS (lmBF)")) +
    by_others +
    labs(x = "sample size (log scale)", y = "AUC (null vs small factor) minus JZS", colour = NULL,
         title = "Testing a factor: separating nulls from small effects, against JZS (averaged over q = 2-5)",
         subtitle = "above 0: separates better than JZS"),
  "11_auc", width = 12, height = 7)

save_figure(
  ggplot(averaged_over_q |>
           pivot_longer(c(null_found, small_called_null), names_to = "rate", values_to = "value") |>
           mutate(rate = recode(rate, null_found = "true null confirmed (BF01 > 3)",
                                small_called_null = "small effect called null (BF01 > 3)")),
         aes(sample_size, value, colour = method, linetype = method)) +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200, 1000)) + method_scales(factor_methods) +
    facet_grid(rate + other_predictors ~ explained_by_others,
      labeller = labeller(other_predictors    = function(v) paste(v, "other predictors"),
                          explained_by_others = function(v) paste("others explain R-squared =", v))) +
    labs(x = "sample size (log scale)", y = "proportion of datasets", colour = NULL,
         title = "Testing a factor: calibration at 3 (averaged over q = 2-5)",
         subtitle = "The conditional tests do not change along a row; the joint tests and JZS do (the conditional Lindley paradox)"),
  "11_calibration", width = 12, height = 11)

report_checks()
