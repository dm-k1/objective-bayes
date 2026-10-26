## =====================================================================
##  Validation 03: testing one coefficient, across many datasets
##
##  4 correlated predictors; x1's slope varies, x2 = 0.3, x3 and x4 are null.
##  Tested: x1 (effect varies) and x3 (always null). 200 replications per cell.
##
##  Methods: the common set of validation/comparison_methods.R, all testing one
##  coefficient (effect against none):
##    JZS default, Zellner-Siow                   model with vs without it
##    flat-R joint, uniform correlations joint    model with vs without it
##    flat-R conditional,                         the coefficient's partial
##    uniform correlations conditional            correlation only
##  (The earlier version, with R2-D2 and an approximate flat partial
##  correlation, is archive/03_simulations_with_r2d2.R.)
##
##  Checks first: U1 the four flat-on-r methods agree with one predictor;
##  U2 the uniform-correlations Bayes factor against a Monte Carlo average.
##
##  Run from the project folder: Rscript validation/03_simulations.R
##  About 5-10 minutes; cached (reruns only when its settings change).
##  Figures are saved to validation/results/figures/.
## =====================================================================
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

pkgload::load_all(".", quiet = TRUE)   # the flatcor package, from this folder
source("validation/settings.R")        # parallel cores, cached simulation runner, record_check()
source("validation/comparison_methods.R")

results_file     <- file.path(results_folder, "sim_one_coefficient.rds")
rerun_simulation <- FALSE   # TRUE forces a rerun even if cached results match

## ---- Settings ------------------------------------------------------------
number_of_predictors   <- 4
slopes_of_x1           <- c(0, 0.1, 0.2, 0.3, 0.5)
slopes_of_x2_to_x4     <- c(0.3, 0, 0)
tested_predictors      <- c(1, 3)
sample_sizes           <- c(20, 50, 100, 200, 500)
predictor_correlations <- c(0, 0.5)               # equicorrelation between all predictors
replications_per_cell  <- 200


## ---- Checks --------------------------------------------------------------
## U1. With one predictor, flat-R and uniform correlations are the same prior
##     (flat on r), and the conditional tests equal the joint ones.
u1_differences <- sapply(list(c(0.2, 60), c(0.05, 100), c(0.5, 30)), function(s) {
  values <- c(log_bayes_factor_model("flat-R", s[1], s[2], 1),
              log_bayes_factor_model("uniform correlations", s[1], s[2], 1),
              log_bayes_factor_flat_r_conditional(sqrt(s[1]), s[2], 1),
              log_bayes_factor_uniform_correlations_conditional(sqrt(s[1]), s[2], 1))
  diff(range(values))
})
record_check("U1. k = 1: flat-R and uniform correlations, joint and conditional, agree",
             max(u1_differences) < 1e-3, sprintf("largest spread %.1e", max(u1_differences)))

## U2. Uniform correlations is the g-prior with g/n ~ Inverse-Gamma(1, 1/2):
##     the grid against a Monte Carlo average over 2,000,000 draws of g.
set.seed(1)
u2_differences <- sapply(list(c(0.15, 100, 4), c(0.3, 50, 8), c(0.02, 500, 3)), function(s) {
  g <- s[2] / (2 * rexp(2e6))
  log_integrand <- (s[2] - 1 - s[3]) / 2 * log1p(g) - (s[2] - 1) / 2 * log1p(g * (1 - s[1]))
  monte_carlo <- max(log_integrand) + log(mean(exp(log_integrand - max(log_integrand))))
  abs(log_bayes_factor_model("uniform correlations", s[1], s[2], s[3]) - monte_carlo)
})
record_check("U2. uniform correlations: grid vs Monte Carlo", max(u2_differences) < 0.01,
             sprintf("largest difference %.4f", max(u2_differences)))


## ---- One dataset ---------------------------------------------------------
simulate_one_dataset <- function(slope_of_x1, sample_size, predictor_correlation) {
  predictor_covariance <- matrix(predictor_correlation, number_of_predictors, number_of_predictors)
  diag(predictor_covariance) <- 1
  true_slopes <- c(slope_of_x1, slopes_of_x2_to_x4)
  predictors  <- matrix(rnorm(sample_size * number_of_predictors), sample_size) %*%
                 chol(predictor_covariance)
  noise_sd    <- sqrt(1 - drop(t(true_slopes) %*% predictor_covariance %*% true_slopes))
  outcome     <- drop(predictors %*% true_slopes) + rnorm(sample_size, sd = noise_sd)
  full_fit    <- lm(outcome ~ predictors)
  map_dfr(tested_predictors, function(tested)
    tibble(coefficient = paste0("x", tested), true_slope = true_slopes[tested],
           log_bayes_factors_one_coefficient(full_fit, tested, outcome, predictors)))
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
                  tested_predictors = tested_predictors, methods = one_coefficient_methods,
                  version = 2),
  simulate_cell = function(cell) {
    map_dfr(seq_len(replications_per_cell), function(replication)
      mutate(simulate_one_dataset(cell$slope_of_x1, cell$sample_size, cell$predictor_correlation),
             replication = replication)) |>
      mutate(slope_of_x1 = cell$slope_of_x1, sample_size = cell$sample_size,
             predictor_correlation = cell$predictor_correlation, .before = 1)
  })
simulation_runs <- mutate(simulation_runs, method = factor(method, levels = one_coefficient_methods))

## ---- Summaries -------------------------------------------------------------
cell_summaries <- simulation_runs |>
  mutate(correct           = ifelse(true_slope == 0, log_BF10 < -log(3),  log_BF10 > log(3)),
         misleading        = ifelse(true_slope == 0, log_BF10 >  log(3),  log_BF10 < -log(3)),
         misleading_strong = ifelse(true_slope == 0, log_BF10 >  log(10), log_BF10 < -log(10))) |>
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
  summarise(null_found          = mean(log_BF10 < -log(3)),
            false_effect        = mean(log_BF10 >  log(3)),
            strong_false_effect = mean(log_BF10 >  log(10)), .groups = "drop") |>
  arrange(predictor_correlation, sample_size, method) |>
  print(n = Inf)

## ---- Figures ---------------------------------------------------------------
x1_summaries   <- filter(cell_summaries, coefficient == "x1")
facets_by_cell <- facet_grid(
  predictor_correlation ~ slope_of_x1,
  labeller = labeller(slope_of_x1           = function(v) paste("slope of x1 =", v),
                      predictor_correlation = function(v) paste("predictor correlation =", v)))
sample_size_axis <- scale_x_log10(breaks = c(20, 100, 500))
bottom_legend    <- theme(legend.position = "bottom")
method_scale     <- scale_colour_manual(values = method_colours[one_coefficient_methods])

save_figure(
  ggplot(x1_summaries, aes(sample_size, correct_3, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    sample_size_axis + method_scale + facets_by_cell + bottom_legend +
    labs(x = "sample size (log scale)", y = "P(BF > 3 for the true hypothesis)", colour = NULL,
         title = "x1: correct evidence",
         subtitle = "slope 0: evidence for no effect; slope > 0: evidence for an effect"),
  "03_correct_evidence_x1", width = 12, height = 6)

save_figure(
  ggplot(x1_summaries, aes(sample_size, misleading_3, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    sample_size_axis + method_scale + facets_by_cell + bottom_legend +
    labs(x = "sample size (log scale)", y = "P(BF > 3 for the wrong hypothesis)", colour = NULL,
         title = "x1: misleading evidence",
         subtitle = "slope 0: a false effect; slope > 0: a real effect called null"),
  "03_misleading_evidence_x1", width = 12, height = 6)

## The trade-off in one picture: how often each method finds the true null
## (x3) against how often it calls a small real effect (x1 = 0.1) null
trade_off <- inner_join(
  simulation_runs |> filter(coefficient == "x3") |>
    group_by(predictor_correlation, sample_size, method) |>
    summarise(null_found = mean(log_BF10 < -log(3)), .groups = "drop"),
  x1_summaries |> filter(slope_of_x1 == 0.1) |>
    select(predictor_correlation, sample_size, method, small_effect_called_null = misleading_3),
  by = c("predictor_correlation", "sample_size", "method"))

save_figure(
  ggplot(trade_off, aes(small_effect_called_null, null_found, colour = method)) +
    geom_path(linewidth = 0.8) +
    geom_point(aes(size = sample_size)) +
    scale_size_continuous(breaks = c(20, 100, 500), range = c(1, 4)) +
    method_scale +
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
