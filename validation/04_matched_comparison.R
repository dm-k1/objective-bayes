## =====================================================================
##  Validation 04: compared at equal error rates
##
##  Study 03 counts how often BF01 > 3. That mixes two things: how well a
##  method separates nulls from small effects (discrimination), and where its
##  Bayes factors happen to sit relative to 3 (calibration). This script
##  separates them.
##
##    M1. discrimination: the AUC, P(BF01 is larger for a null dataset than for
##        a dataset with a small effect), for x1 = 0 against x1 = 0.1 and 0.2,
##        from the datasets of study 03. Threshold-free: equal AUC means any
##        difference in 03 comes from calibration.
##    M2. flat-R minus JZS in AUC, with a paired bootstrap (both methods were
##        computed on the same datasets).
##    M3. matched threshold: each method's BF01 threshold at which a small
##        effect (x1 = 0.1) is called null in 10% of datasets, and how many
##        true nulls it then finds.
##    M4. why: how strongly each method's evidence for x1 depends on how much
##        the OTHER predictors explain, at a fixed t-statistic for x1.
##
##  Needs the results of study 03 (validation/results/sim_joint_prior.rds).
##  Run from the project folder: Rscript validation/04_matched_comparison.R
## =====================================================================
library(dplyr)
library(ggplot2)
source("validation/settings.R")   # results folder, save_figure()
pkgload::load_all(".", quiet = TRUE)

simulation_runs <- readRDS(file.path(results_folder, "sim_joint_prior.rds"))

## P(a value from `higher` exceeds a value from `lower`), ties count half
area_under_curve <- function(higher, lower)
  mean(outer(higher, lower, ">")) + 0.5 * mean(outer(higher, lower, "=="))


## ---- M1 and M3: discrimination and matched thresholds -------------------
x1_runs <- simulation_runs |>
  filter(coefficient == "x1", slope_of_x1 %in% c(0, 0.1, 0.2)) |>
  mutate(BF01 = 1 / BF10)

discrimination <- x1_runs |>
  group_by(predictor_correlation, sample_size, method) |>
  summarise(
    auc_null_vs_0.1 = area_under_curve(BF01[slope_of_x1 == 0], BF01[slope_of_x1 == 0.1]),
    auc_null_vs_0.2 = area_under_curve(BF01[slope_of_x1 == 0], BF01[slope_of_x1 == 0.2]),
    matched_threshold  = quantile(BF01[slope_of_x1 == 0.1], probs = 0.90),
    nulls_found_matched = mean(BF01[slope_of_x1 == 0] > quantile(BF01[slope_of_x1 == 0.1], probs = 0.90)),
    .groups = "drop")
print(as.data.frame(discrimination), digits = 3)

save_figure(
  ggplot(discrimination, aes(sample_size, auc_null_vs_0.1, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    scale_x_log10(breaks = c(20, 100, 500)) +
    scale_colour_manual(values = c(`JZS default` = "grey40", `joint flat-R` = "#d95f02",
                                   `flat partial correlation` = "#1b9e77",
                                   `R2-D2 matched` = "#7570b3")) +
    facet_wrap(~ predictor_correlation,
               labeller = labeller(predictor_correlation = function(v) paste("predictor correlation =", v))) +
    labs(x = "sample size (log scale)", y = "AUC: null vs slope 0.1", colour = NULL,
         title = "Separating nulls from small effects, at any threshold",
         subtitle = "0.5 = no better than chance; equal curves = equal discrimination") +
    theme(legend.position = "bottom"),
  "04_discrimination", width = 10, height = 5)


## ---- M2: flat-R minus JZS, paired bootstrap -----------------------------
set.seed(1)
paired_differences <- do.call(rbind, lapply(c(0, 0.5), function(correlation)
  do.call(rbind, lapply(c(100, 200, 500), function(n) {
    cell <- filter(x1_runs, predictor_correlation == correlation, sample_size == n)
    bf01_of <- function(method_name, slope)
      cell$BF01[cell$method == method_name & cell$slope_of_x1 == slope]
    difference <- function(null_rows, effect_rows)
      area_under_curve(bf01_of("joint flat-R", 0)[null_rows], bf01_of("joint flat-R", 0.1)[effect_rows]) -
      area_under_curve(bf01_of("JZS default", 0)[null_rows],  bf01_of("JZS default", 0.1)[effect_rows])
    replications <- length(bf01_of("JZS default", 0))
    bootstrap <- replicate(1000, difference(sample(replications, replace = TRUE),
                                            sample(replications, replace = TRUE)))
    data.frame(predictor_correlation = correlation, sample_size = n,
               auc_flat_r_minus_jzs = difference(seq_len(replications), seq_len(replications)),
               lower_95 = quantile(bootstrap, 0.025), upper_95 = quantile(bootstrap, 0.975),
               row.names = NULL)
  }))))
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
  full_fit <- summary(lm(outcome ~ predictors))
  r_squared_full    <- full_fit$r.squared
  r_squared_reduced <- summary(lm(outcome ~ predictors[, -1]))$r.squared
  c(t_squared  = full_fit$coefficients[2, 3]^2,
    r_squared_reduced = r_squared_reduced,
    log_bf_flat_r = log_bayes_factor_flat_r_vs_intercept_only(r_squared_full, sample_size, 4) -
                    log_bayes_factor_flat_r_vs_intercept_only(r_squared_reduced, sample_size, 3),
    log_bf_jzs =
      BayesFactor::linearReg.R2stat(sample_size, 4, r_squared_full, rscale = sqrt(2) / 4, simple = FALSE)$bf -
      BayesFactor::linearReg.R2stat(sample_size, 3, r_squared_reduced, rscale = sqrt(2) / 4, simple = FALSE)$bf)
}
datasets <- as.data.frame(t(replicate(1500, one_dataset())))
dependence <- rbind(
  `joint flat-R` = coef(lm(log_bf_flat_r ~ t_squared + r_squared_reduced, data = datasets)),
  `JZS default`  = coef(lm(log_bf_jzs ~ t_squared + r_squared_reduced, data = datasets)))
print(round(dependence, 2))
cat(sprintf("Rank correlation of log BF with t^2: flat-R %.2f, JZS %.2f\n",
            cor(datasets$log_bf_flat_r, datasets$t_squared, method = "spearman"),
            cor(datasets$log_bf_jzs, datasets$t_squared, method = "spearman")))
