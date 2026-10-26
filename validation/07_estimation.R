## =====================================================================
##  Validation 07: how good are the estimates?
##
##  The prior shrinks estimates toward zero when several predictors share the
##  budget. Does that help or hurt? Fixed true coefficients, many datasets,
##  and for every coefficient: bias, error (RMSE), and whether the 95%
##  interval covers the truth; per dataset: the excess prediction error for a
##  new observation, (b_hat - b)' Sigma (b_hat - b). Methods (the estimation
##  versions of the common set in validation/comparison_methods.R):
##    - flat-R                fit_flat_r_model(), posterior mean and 95% interval
##    - uniform correlations  rstanarm stan_lm() with R-squared ~ Beta(k/2, 1)
##                            (R2(location = k/(k+2), what = "mean"), i.e. eta = 1)
##    - JZS default           BayesFactor posterior (Cauchy, scale 0.354)
##    - Zellner-Siow          BayesFactor posterior (Cauchy, scale 1)
##    - ordinary              lm(): estimate and 95% confidence interval
##  (The conditional tests have no estimation version: they are tests.)
##
##  Design (noise SD 1, effects in noise units):
##    k (predictors)    3, 8
##    n                 30, 60, 200
##    correlation       0, 0.5 (all pairs equal)
##    true coefficients sparse: 0.4, 0.2, and the rest 0
##                      dense:  all 0.15
##  100 datasets per cell: 2,400 flat-R fits and 2,400 stan_lm fits.
##  Needs cmdstanr (or rstan), rstanarm and BayesFactor.
##
##  Run from the project folder, in a terminal, and leave it running:
##    Rscript validation/07_estimation.R
##  Roughly 2.5-4 hours. Each design cell is saved when it finishes
##  (validation/results/estimation_v2/), so an interrupted run resumes where
##  it stopped. (estimation/, from the version without stan_lm, is not used.)
## =====================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
source("validation/settings.R")   # results folder, save_figure()
pkgload::load_all(".", quiet = TRUE)

cell_folder <- file.path(results_folder, "estimation_v2")
dir.create(cell_folder, recursive = TRUE, showWarnings = FALSE)

design <- expand_grid(number_of_predictors  = c(3, 8),
                      sample_size           = c(30, 60, 200),
                      predictor_correlation = c(0, 0.5),
                      true_coefficients     = c("sparse", "dense"))
replications_per_cell <- 100

true_slopes_for <- function(k, pattern)
  if (pattern == "sparse") c(0.4, 0.2, rep(0, k - 2)) else rep(0.15, k)


## ---- One dataset: estimates and intervals from every method -----------------
estimate_one_dataset <- function(k, n, correlation, pattern, seed) {
  set.seed(seed)
  covariance <- matrix(correlation, k, k); diag(covariance) <- 1
  true_slopes <- true_slopes_for(k, pattern)
  predictors  <- matrix(rnorm(n * k), n) %*% chol(covariance)
  colnames(predictors) <- paste0("x", seq_len(k))
  data <- data.frame(predictors, y = drop(predictors %*% true_slopes) + rnorm(n))
  formula <- reformulate(colnames(predictors), response = "y")

  ## flat-R
  flat_r_fit <- suppressWarnings(fit_flat_r_model(formula, data, chains = 4, iter_warmup = 1000,
                                                  iter_sampling = 1000, seed = seed))
  flat_r_draws <- as.data.frame(flat_r_fit$draws)[, colnames(predictors), drop = FALSE]

  ## uniform correlations: rstanarm stan_lm, R-squared ~ Beta(k/2, eta) with
  ## mean (k/2) / (k/2 + eta); mean k / (k + 2) gives eta = 1
  stan_lm_fit <- suppressWarnings(rstanarm::stan_lm(
    formula, data = data, prior = rstanarm::R2(location = k / (k + 2), what = "mean"),
    chains = 4, iter = 2000, seed = seed, refresh = 0))
  stan_lm_draws <- as.data.frame(as.matrix(stan_lm_fit, pars = colnames(predictors)))
  stan_lm_divergent <- sum(sapply(rstan::get_sampler_params(stan_lm_fit$stanfit, inc_warmup = FALSE),
                                  function(chain) sum(chain[, "divergent__"])))

  ## JZS default and Zellner-Siow (posterior samples from BayesFactor)
  bayesfactor_draws <- function(scale)
    as.data.frame(BayesFactor::lmBF(formula, data, rscaleCont = scale, posterior = TRUE,
                                    iterations = 4000, progress = FALSE))[, colnames(predictors), drop = FALSE]
  jzs_draws <- bayesfactor_draws(sqrt(2) / 4)
  zs_draws  <- bayesfactor_draws(1)

  ## ordinary least squares
  ols_fit <- lm(formula, data)
  ols_interval <- confint(ols_fit)[colnames(predictors), ]

  summarise_draws <- function(draws, method) data.frame(
    method = method, coefficient = colnames(predictors), true_value = true_slopes,
    estimate = colMeans(draws),
    lower = apply(draws, 2, quantile, probs = 0.025),
    upper = apply(draws, 2, quantile, probs = 0.975))
  rbind(summarise_draws(flat_r_draws, "flat-R"),
        summarise_draws(stan_lm_draws, "uniform correlations"),
        summarise_draws(jzs_draws, "JZS default"),
        summarise_draws(zs_draws, "Zellner-Siow"),
        data.frame(method = "ordinary (lm)", coefficient = colnames(predictors),
                   true_value = true_slopes, estimate = coef(ols_fit)[colnames(predictors)],
                   lower = ols_interval[, 1], upper = ols_interval[, 2])) |>
    group_by(method) |>
    ## excess error predicting a new observation: (b_hat - b)' Sigma (b_hat - b)
    mutate(excess_prediction_error = drop(t(estimate - true_value) %*% covariance %*%
                                            (estimate - true_value))) |>
    ungroup() |>
    mutate(divergent_transitions = ifelse(method == "flat-R", flat_r_fit$divergent,
                                          ifelse(method == "uniform correlations", stan_lm_divergent, NA)))
}


## ---- Run, one design cell at a time (each saved when done) ------------------
started <- Sys.time()
for (row in seq_len(nrow(design))) {
  cell <- design[row, ]
  cell_file <- file.path(cell_folder, sprintf("cell_%02d.rds", row))
  if (file.exists(cell_file)) next
  message(sprintf("[%s] cell %d of %d: k = %d, n = %d, correlation = %.1f, %s",
                  format(Sys.time(), "%H:%M"), row, nrow(design), cell$number_of_predictors,
                  cell$sample_size, cell$predictor_correlation, cell$true_coefficients))
  cell_results <- do.call(rbind, lapply(seq_len(replications_per_cell), function(replication)
    cbind(as.data.frame(cell), replication = replication,
          estimate_one_dataset(cell$number_of_predictors, cell$sample_size,
                               cell$predictor_correlation, cell$true_coefficients,
                               seed = 10000 * row + replication))))
  saveRDS(cell_results, cell_file)
}
message("Done in ", format(round(Sys.time() - started, 1)))

results <- do.call(rbind, lapply(sort(list.files(cell_folder, pattern = "^cell_.*\\.rds$",
                                                 full.names = TRUE)), readRDS))


## ---- Summaries -----------------------------------------------------------------
results <- results |>
  mutate(error = estimate - true_value,
         covered = lower <= true_value & true_value <= upper,
         interval_width = upper - lower,
         kind = ifelse(true_value == 0, "true zero", "true effect"))

for (stan_method in c("flat-R", "uniform correlations")) {
  divergences_per_fit <- results |>
    filter(method == stan_method) |>
    distinct(number_of_predictors, sample_size, predictor_correlation, true_coefficients,
             replication, divergent_transitions)
  cat(sprintf("\n%s fits with divergent transitions: %d of %d (total %d transitions)\n", stan_method,
              sum(divergences_per_fit$divergent_transitions > 0), nrow(divergences_per_fit),
              sum(divergences_per_fit$divergent_transitions)))
}

per_cell <- results |>
  group_by(number_of_predictors, sample_size, predictor_correlation, true_coefficients,
           method, kind) |>
  summarise(bias     = mean(error),
            rmse     = sqrt(mean(error^2)),
            coverage = mean(covered),
            width    = mean(interval_width), .groups = "drop")
saveRDS(per_cell, file.path(results_folder, "estimation_summaries.rds"))

## Total squared error over all coefficients, relative to ordinary least squares
## (below 1: the method estimates the whole coefficient vector better)
relative_error <- results |>
  group_by(number_of_predictors, sample_size, predictor_correlation, true_coefficients,
           method, replication) |>
  summarise(squared_error = sum(error^2), .groups = "drop") |>
  group_by(number_of_predictors, sample_size, predictor_correlation, true_coefficients, method) |>
  summarise(mean_squared_error = mean(squared_error), .groups = "drop") |>
  group_by(number_of_predictors, sample_size, predictor_correlation, true_coefficients) |>
  mutate(relative_to_ordinary = mean_squared_error /
           mean_squared_error[method == "ordinary (lm)"]) |>
  ungroup()

cat("\nTotal squared error relative to ordinary least squares (below 1 = better):\n")
print(as.data.frame(relative_error |>
  select(-mean_squared_error) |>
  pivot_wider(names_from = method, values_from = relative_to_ordinary)), digits = 3)

## Excess prediction error for a new observation, relative to ordinary least squares
relative_prediction_error <- results |>
  distinct(number_of_predictors, sample_size, predictor_correlation, true_coefficients,
           method, replication, excess_prediction_error) |>
  group_by(number_of_predictors, sample_size, predictor_correlation, true_coefficients, method) |>
  summarise(mean_excess = mean(excess_prediction_error), .groups = "drop") |>
  group_by(number_of_predictors, sample_size, predictor_correlation, true_coefficients) |>
  mutate(relative_to_ordinary = mean_excess / mean_excess[method == "ordinary (lm)"]) |>
  ungroup()

cat("\nExcess prediction error relative to ordinary least squares (below 1 = better):\n")
print(as.data.frame(relative_prediction_error |>
  select(-mean_excess) |>
  pivot_wider(names_from = method, values_from = relative_to_ordinary)), digits = 3)

cat("\nCoverage of the 95% intervals (should be about 0.95 for the ordinary interval):\n")
print(as.data.frame(per_cell |>
  select(number_of_predictors, sample_size, predictor_correlation, true_coefficients,
         method, kind, coverage) |>
  pivot_wider(names_from = method, values_from = coverage)), digits = 3)


## ---- Figures ---------------------------------------------------------------------
source("validation/comparison_methods.R")   # method_colours
method_colours <- method_colours[c("flat-R", "uniform correlations", "JZS default", "Zellner-Siow",
                                   "ordinary (lm)")]
cells_facet <- facet_grid(
  number_of_predictors + predictor_correlation ~ true_coefficients,
  labeller = labeller(number_of_predictors  = function(v) paste("k =", v),
                      predictor_correlation = function(v) paste("cor =", v)))

save_figure(
  ggplot(relative_error, aes(sample_size, relative_to_ordinary, colour = method)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_line(linewidth = 0.9) + geom_point() +
    scale_x_log10(breaks = c(30, 60, 200)) + scale_colour_manual(values = method_colours) +
    cells_facet + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "total squared error / ordinary least squares",
         colour = NULL, title = "Does the shrinkage pay off?",
         subtitle = "below 1: closer to the true coefficients than ordinary regression"),
  "07_relative_error", width = 9, height = 10)

save_figure(
  ggplot(relative_prediction_error, aes(sample_size, relative_to_ordinary, colour = method)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_line(linewidth = 0.9) + geom_point() +
    scale_x_log10(breaks = c(30, 60, 200)) + scale_colour_manual(values = method_colours) +
    cells_facet + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "excess prediction error / ordinary least squares",
         colour = NULL, title = "Does the shrinkage pay off for prediction?",
         subtitle = "below 1: predicts new observations better than ordinary regression"),
  "07_prediction_error", width = 9, height = 10)

save_figure(
  ggplot(per_cell, aes(sample_size, coverage, colour = method, linetype = kind)) +
    geom_hline(yintercept = 0.95, linetype = "dotted") +
    geom_line(linewidth = 0.8) + geom_point() +
    scale_x_log10(breaks = c(30, 60, 200)) + scale_colour_manual(values = method_colours) +
    cells_facet + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "coverage of the 95% interval", colour = NULL,
         linetype = NULL, title = "Do the intervals contain the true coefficient?"),
  "07_coverage", width = 9, height = 10)

save_figure(
  ggplot(per_cell, aes(sample_size, bias, colour = method, linetype = kind)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 0.8) + geom_point() +
    scale_x_log10(breaks = c(30, 60, 200)) + scale_colour_manual(values = method_colours) +
    cells_facet + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "mean of (estimate - truth)", colour = NULL,
         linetype = NULL, title = "Bias: how far estimates are pulled toward zero"),
  "07_bias", width = 9, height = 10)
