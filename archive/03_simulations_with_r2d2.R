## Archived (October 2026): study 03 as it was before the common method set.
## Compared JZS, joint flat-R, an approximate flat partial correlation and R2-D2;
## its results are the 03 section of validation/README.md before that change.
## Needs bayes_factor_jzs_one_coefficient() and the R2-D2 functions from
## validation/comparison_methods.R.

## =====================================================================
##  Validation 03: how does the flat-R prior behave across many datasets?
##  4 correlated predictors; x1's slope varies, x2 = 0.3, x3 and x4 are null.
##  Tested: x1 (effect varies) and x3 (always null). 200 replications per cell.
##
##  Methods (all Bayes factors for one coefficient, model with vs without it):
##   - JZS default (BayesFactor / JASP)
##   - joint flat-R (this project)
##   - flat partial correlation, one coefficient at a time
##   - R2-D2 with the same R-squared prior, Dirichlet split
##
##  Uses the closed-form Bayes factor (checked in validation/02_bayes_factor.R).
##
##  Run from the project folder. The simulation runs in parallel and is
##  cached: it only reruns when its settings change (or rerun_simulation = TRUE).
##  Figures are saved to validation/results/figures/.
## =====================================================================
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

pkgload::load_all(".", quiet = TRUE)   # the flatr package, from this folder
source("validation/comparison_methods.R")
source("validation/settings.R")   # parallel cores, cached simulation runner, record_check()

results_file     <- file.path(results_folder, "sim_joint_prior.rds")
rerun_simulation <- FALSE   # TRUE forces a rerun even if cached results match

## ---- Settings ------------------------------------------------------------
number_of_predictors   <- 4
slopes_of_x1           <- c(0, 0.1, 0.2, 0.3, 0.5)
slopes_of_x2_to_x4     <- c(0.3, 0, 0)
tested_predictors      <- c(1, 3)
sample_sizes           <- c(20, 50, 100, 200, 500)
predictor_correlations <- c(0, 0.5)               # equicorrelation between all predictors
replications_per_cell  <- 200
jzs_scale              <- sqrt(2) / 4

## R2-D2 prior draws, shared by all datasets (full and reduced models)
set.seed(99)
r2d2_draws_full    <- draw_from_r2d2_prior(number_of_predictors)
r2d2_draws_reduced <- draw_from_r2d2_prior(number_of_predictors - 1)

## ---- Monte Carlo check: two independent draw sets, same dataset ---------
## The R2-D2 Bayes factor is a Monte Carlo average, so two independent sets
## of prior draws should give nearly the same answer. Criterion: within 0.1
## on the log scale.
set.seed(5)
check_predictors <- matrix(rnorm(100 * 4), 100)
check_outcome    <- drop(check_predictors %*% c(0.2, 0.3, 0, 0)) + rnorm(100)
centered_check_predictors <- scale(check_predictors, scale = FALSE)
centered_check_outcome    <- check_outcome - mean(check_outcome)
log_bf_draw_set_a <- log_bayes_factor_r2d2_vs_intercept_only(
  centered_check_predictors, centered_check_outcome, draw_from_r2d2_prior(4))
log_bf_draw_set_b <- log_bayes_factor_r2d2_vs_intercept_only(
  centered_check_predictors, centered_check_outcome, draw_from_r2d2_prior(4))
record_check("R2-D2 Monte Carlo error (two draw sets)",
             abs(log_bf_draw_set_a - log_bf_draw_set_b) < 0.1,
             sprintf("log BF %.3f vs %.3f", log_bf_draw_set_a, log_bf_draw_set_b))

## ---- One dataset ---------------------------------------------------------
simulate_one_dataset <- function(slope_of_x1, sample_size, predictor_correlation) {
  predictor_covariance <- matrix(predictor_correlation, number_of_predictors, number_of_predictors)
  diag(predictor_covariance) <- 1
  true_slopes <- c(slope_of_x1, slopes_of_x2_to_x4)
  predictors  <- matrix(rnorm(sample_size * number_of_predictors), sample_size) %*%
                 chol(predictor_covariance)
  noise_sd    <- sqrt(1 - drop(t(true_slopes) %*% predictor_covariance %*% true_slopes))
  outcome     <- drop(predictors %*% true_slopes) + rnorm(sample_size, sd = noise_sd)
  centered_predictors <- scale(predictors, scale = FALSE)
  centered_outcome    <- outcome - mean(outcome)

  full_fit       <- lm(outcome ~ predictors)
  r_squared_full <- summary(full_fit)$r.squared
  log_bf_flat_r_full <- log_bayes_factor_flat_r_vs_intercept_only(
    r_squared_full, sample_size, number_of_predictors)
  log_bf_r2d2_full <- log_bayes_factor_r2d2_vs_intercept_only(
    centered_predictors, centered_outcome, r2d2_draws_full)

  map_dfr(tested_predictors, function(tested) {
    r_squared_without <- summary(lm(outcome ~ predictors[, -tested]))$r.squared
    log_bf_flat_r_without <- log_bayes_factor_flat_r_vs_intercept_only(
      r_squared_without, sample_size, number_of_predictors - 1)
    log_bf_r2d2_without <- log_bayes_factor_r2d2_vs_intercept_only(
      centered_predictors[, -tested, drop = FALSE], centered_outcome, r2d2_draws_reduced)
    tibble(coefficient = paste0("x", tested),
           true_slope  = true_slopes[tested],
           method      = c("JZS default", "joint flat-R",
                           "flat partial correlation", "R2-D2 matched"),
           BF10        = c(bayes_factor_jzs_one_coefficient(r_squared_full, r_squared_without,
                                                            sample_size, number_of_predictors,
                                                            rscale = jzs_scale),
                           exp(log_bf_flat_r_full - log_bf_flat_r_without),
                           bayes_factor_flat_partial_correlation(full_fit, tested, predictors),
                           exp(log_bf_r2d2_full - log_bf_r2d2_without)))
  })
}

## ---- Run (parallel; cached until the settings change) --------------------
simulation_runs <- run_simulation_grid(
  expand_grid(slope_of_x1           = slopes_of_x1,
              sample_size           = sample_sizes,
              predictor_correlation = predictor_correlations),
  results_file = results_file, rerun = rerun_simulation,
  settings = list(replications_per_cell = replications_per_cell,
                  number_of_predictors = number_of_predictors,
                  slopes_of_x2_to_x4 = slopes_of_x2_to_x4,
                  tested_predictors = tested_predictors, jzs_scale = jzs_scale),
  simulate_cell = function(cell) {
    map_dfr(seq_len(replications_per_cell), function(replication)
      mutate(simulate_one_dataset(cell$slope_of_x1, cell$sample_size, cell$predictor_correlation),
             replication = replication)) |>
      mutate(slope_of_x1 = cell$slope_of_x1, sample_size = cell$sample_size,
             predictor_correlation = cell$predictor_correlation, .before = 1)
  })

## ---- Summaries -------------------------------------------------------------
cell_summaries <- simulation_runs |>
  mutate(correct           = ifelse(true_slope == 0, BF10 < 1/3, BF10 > 3),
         misleading        = ifelse(true_slope == 0, BF10 > 3,   BF10 < 1/3),
         misleading_strong = ifelse(true_slope == 0, BF10 > 10,  BF10 < 1/10)) |>
  group_by(coefficient, predictor_correlation, slope_of_x1, sample_size, method) |>
  summarise(correct_3     = mean(correct),
            misleading_3  = mean(misleading),
            misleading_10 = mean(misleading_strong),
            runs          = n(), .groups = "drop")

## x1: slope varies
cell_summaries |>
  filter(coefficient == "x1") |>
  select(predictor_correlation, slope_of_x1, sample_size, method,
         correct_3, misleading_3, misleading_10) |>
  arrange(predictor_correlation, slope_of_x1, sample_size, method) |>
  print(n = Inf)

## x3: always null (pooled over the slopes of x1)
simulation_runs |>
  filter(coefficient == "x3") |>
  group_by(predictor_correlation, sample_size, method) |>
  summarise(null_found          = mean(BF10 < 1/3),
            false_effect        = mean(BF10 > 3),
            strong_false_effect = mean(BF10 > 10), .groups = "drop") |>
  arrange(predictor_correlation, sample_size, method) |>
  print(n = Inf)

## ---- Figures ---------------------------------------------------------------
method_colours <- c(`JZS default` = "grey40", `joint flat-R` = "#d95f02",
                    `flat partial correlation` = "#1b9e77", `R2-D2 matched` = "#7570b3")
x1_summaries   <- filter(cell_summaries, coefficient == "x1")
facets_by_cell <- facet_grid(
  predictor_correlation ~ slope_of_x1,
  labeller = labeller(slope_of_x1           = function(v) paste("slope of x1 =", v),
                      predictor_correlation = function(v) paste("predictor correlation =", v)))
sample_size_axis <- scale_x_log10(breaks = c(20, 100, 500))
bottom_legend    <- theme(legend.position = "bottom")

save_figure(
  ggplot(x1_summaries, aes(sample_size, correct_3, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    sample_size_axis + scale_colour_manual(values = method_colours) +
    facets_by_cell + bottom_legend +
    labs(x = "sample size (log scale)", y = "P(BF > 3 for the true hypothesis)", colour = NULL,
         title = "x1: correct evidence",
         subtitle = "slope 0: evidence for no effect; slope > 0: evidence for an effect"),
  "03_correct_evidence_x1", width = 12, height = 6)

save_figure(
  ggplot(x1_summaries, aes(sample_size, misleading_3, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    sample_size_axis + scale_colour_manual(values = method_colours) +
    facets_by_cell + bottom_legend +
    labs(x = "sample size (log scale)", y = "P(BF > 3 for the wrong hypothesis)", colour = NULL,
         title = "x1: misleading evidence",
         subtitle = "slope 0: a false effect; slope > 0: a real effect called null"),
  "03_misleading_evidence_x1", width = 12, height = 6)

## The trade-off in one picture: how often each method finds the true null
## (x3) against how often it calls a small real effect (x1 = 0.1) null
trade_off <- inner_join(
  simulation_runs |> filter(coefficient == "x3") |>
    group_by(predictor_correlation, sample_size, method) |>
    summarise(null_found = mean(BF10 < 1/3), .groups = "drop"),
  x1_summaries |> filter(slope_of_x1 == 0.1) |>
    select(predictor_correlation, sample_size, method, small_effect_called_null = misleading_3),
  by = c("predictor_correlation", "sample_size", "method"))

save_figure(
  ggplot(trade_off, aes(small_effect_called_null, null_found, colour = method)) +
    geom_path(linewidth = 0.8) +
    geom_point(aes(size = sample_size)) +
    scale_size_continuous(breaks = c(20, 100, 500), range = c(1, 4)) +
    scale_colour_manual(values = method_colours) +
    facet_wrap(~ predictor_correlation,
               labeller = labeller(predictor_correlation = function(v) paste("predictor correlation =", v))) +
    bottom_legend +
    labs(x = "P(small real effect, slope 0.1, called null: BF01 > 3)",
         y = "P(true null found: BF01 > 3)", colour = NULL, size = "n",
         title = "Finding nulls vs calling small effects null",
         subtitle = "top left is best; each line runs from n = 20 to n = 500"),
  "03_trade_off", width = 11, height = 5.5)
## the docs show this one too
file.copy(file.path(figures_folder, "03_trade_off.png"),
          file.path("docs", "figures", "7_trade_off.png"), overwrite = TRUE)

report_checks()
