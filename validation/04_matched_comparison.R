## =====================================================================
##  Validation 04: compared at equal error rates
##
##  Study 03 counts how often BF01 > 3. That mixes two things: how well a
##  method separates nulls from small effects (discrimination), and where its
##  Bayes factors happen to sit relative to 3 (calibration). This script
##  separates them, for every method of the common set: the three contenders
##  (flat partial r, flat R, flat R-squared), joint and conditional, and the
##  references JZS default and Zellner-Siow.
##
##    M1. discrimination: the AUC, P(BF01 is larger for a null dataset than for
##        a dataset with a small effect), for x1 = 0 against x1 = 0.1 and 0.2,
##        from the datasets of study 03. Threshold-free: equal AUC means any
##        difference in 03 comes from calibration.
##    M2. each method minus JZS in AUC, with a paired bootstrap (all methods
##        were computed on the same datasets).
##    M3. matched threshold: each method's BF01 threshold at which a small
##        effect (x1 = 0.1) is called null in 10% of datasets, and how many
##        true nulls it then finds.
##    M4. why: how strongly each method's evidence for x1 depends on how much
##        the OTHER predictors explain, at a fixed t-statistic for x1.
##
##  Needs the results of study 03 (validation/results/sim_one_coefficient.rds).
##  Run from the project folder: Rscript validation/04_matched_comparison.R
##  About 3 minutes.
## =====================================================================
library(dplyr)
library(ggplot2)
source("validation/settings.R")   # results folder, save_figure()
pkgload::load_all(".", quiet = TRUE)
source("validation/comparison_methods.R")

simulation_runs <- readRDS(file.path(results_folder, "sim_one_coefficient.rds"))

## P(a value from `higher` exceeds a value from `lower`), ties count half
area_under_curve <- function(higher, lower)
  mean(outer(higher, lower, ">")) + 0.5 * mean(outer(higher, lower, "=="))


## ---- M1 and M3: discrimination and matched thresholds -------------------
x1_runs <- simulation_runs |>
  filter(coefficient == "x1", slope_of_x1 %in% c(0, 0.1, 0.2)) |>
  mutate(log_BF01 = -log_BF10, method = factor(method, levels = one_coefficient_methods))

discrimination <- x1_runs |>
  group_by(predictor_correlation, sample_size, method) |>
  summarise(
    auc_null_vs_0.1 = area_under_curve(log_BF01[slope_of_x1 == 0], log_BF01[slope_of_x1 == 0.1]),
    auc_null_vs_0.2 = area_under_curve(log_BF01[slope_of_x1 == 0], log_BF01[slope_of_x1 == 0.2]),
    matched_threshold_BF01 = exp(quantile(log_BF01[slope_of_x1 == 0.1], probs = 0.90)),
    nulls_found_matched = mean(log_BF01[slope_of_x1 == 0] > quantile(log_BF01[slope_of_x1 == 0.1], probs = 0.90)),
    .groups = "drop")
print(as.data.frame(discrimination), digits = 3)

save_figure(
  ggplot(discrimination, aes(sample_size, auc_null_vs_0.1, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    scale_x_log10(breaks = c(20, 100, 500)) +
    scale_colour_manual(values = method_colours[one_coefficient_methods]) +
    facet_wrap(~ predictor_correlation,
               labeller = labeller(predictor_correlation = function(v) paste("predictor correlation =", v))) +
    labs(x = "sample size (log scale)", y = "AUC: null vs slope 0.1", colour = NULL,
         title = "Separating nulls from small effects, at any threshold",
         subtitle = "0.5 = no better than chance; equal curves = equal discrimination") +
    theme(legend.position = "bottom"),
  "04_discrimination", width = 10, height = 5)


## ---- M2: each method minus JZS, paired bootstrap --------------------------
set.seed(1)
paired_differences <- do.call(rbind, lapply(c(0, 0.5), function(correlation)
  do.call(rbind, lapply(c(100, 200, 500), function(n) {
    cell <- filter(x1_runs, predictor_correlation == correlation, sample_size == n)
    log_bf01_of <- function(method_name, slope)
      cell$log_BF01[cell$method == method_name & cell$slope_of_x1 == slope][
        order(cell$replication[cell$method == method_name & cell$slope_of_x1 == slope])]
    replications <- length(log_bf01_of("JZS default", 0))
    bootstrap_rows <- replicate(1000, list(sample(replications, replace = TRUE),
                                           sample(replications, replace = TRUE)), simplify = FALSE)
    auc_of <- function(method_name, null_rows, effect_rows)
      area_under_curve(log_bf01_of(method_name, 0)[null_rows], log_bf01_of(method_name, 0.1)[effect_rows])
    do.call(rbind, lapply(setdiff(one_coefficient_methods, "JZS default"), function(method_name) {
      difference <- function(rows) auc_of(method_name, rows[[1]], rows[[2]]) -
                                   auc_of("JZS default", rows[[1]], rows[[2]])
      bootstrap <- vapply(bootstrap_rows, difference, numeric(1))
      data.frame(predictor_correlation = correlation, sample_size = n, method = method_name,
                 auc_minus_jzs = difference(list(seq_len(replications), seq_len(replications))),
                 lower_95 = quantile(bootstrap, 0.025), upper_95 = quantile(bootstrap, 0.975),
                 row.names = NULL)
    }))
  }))))
cat("\nAUC (null vs 0.1) minus JZS, paired bootstrap 95% interval:\n")
print(paired_differences, digits = 3)


## ---- M4: why: dependence on what the other predictors explain ------------
## Same design as study 03 (4 predictors, correlation 0.5, x2 = 0.3, n = 200),
## x1 null. Regress each method's log BF for x1 on the t-statistic of x1
## (squared) and on the R-squared of the model without x1.
set.seed(2)
sample_size <- 200
number_of_predictors <- 4
covariance <- matrix(0.5, number_of_predictors, number_of_predictors); diag(covariance) <- 1
one_dataset <- function() {
  slopes <- c(0, 0.3, 0, 0)
  predictors <- matrix(rnorm(sample_size * number_of_predictors), sample_size) %*% chol(covariance)
  outcome <- drop(predictors %*% slopes) +
             rnorm(sample_size, sd = sqrt(1 - drop(t(slopes) %*% covariance %*% slopes)))
  full_fit <- lm(outcome ~ predictors)
  log_bfs <- log_bayes_factors_one_coefficient(full_fit, 1, outcome, predictors)
  data.frame(t_squared = summary(full_fit)$coefficients[2, 3]^2,
             r_squared_reduced = summary(lm(outcome ~ predictors[, -1]))$r.squared,
             method = log_bfs$method, log_BF10 = log_bfs$log_BF10)
}
datasets <- do.call(rbind, lapply(1:1500, function(i) one_dataset()))
dependence <- datasets |>
  group_by(method = factor(method, levels = one_coefficient_methods)) |>
  summarise(per_unit_t_squared = coef(lm(log_BF10 ~ t_squared + r_squared_reduced))[2],
            per_unit_r_squared_others = coef(lm(log_BF10 ~ t_squared + r_squared_reduced))[3],
            rank_correlation_with_t_squared = cor(log_BF10, t_squared, method = "spearman"),
            .groups = "drop")
cat("\nM4. Dependence of log BF for x1 on t^2 and on the R-squared of the other predictors:\n")
print(as.data.frame(dependence), digits = 3)
