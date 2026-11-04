## =====================================================================
##  Validation 07: how good are the estimates?
##
##  The priors shrink estimates toward zero when several predictors share the
##  total. Does that help or hurt? Fixed true coefficients, many datasets,
##  and for every coefficient: bias, error (RMSE), and whether the 95%
##  interval covers the truth; per dataset: the excess prediction error for a
##  new observation, (b_hat - b)' Sigma (b_hat - b). Methods (the estimation
##  versions of the common set in validation/comparison_methods.R):
##    the three contenders, fit_regression() with the same Stan model,
##    posterior mean and 95% interval:
##    - flat partial r        prior = "flat_partial_r": R-squared ~ Beta(k/2, 1)
##                            (the prior of rstanarm stan_lm with eta = 1)
##    - flat R                prior = "flat_r": R-squared ~ Beta(1/2, 1)
##    - flat R-squared        prior = "flat_r_squared": R-squared ~ Beta(1, 1)
##    two references:
##    - JZS default           BayesFactor posterior (Cauchy, scale 0.354)
##    - Zellner-Siow          BayesFactor posterior (Cauchy, scale 1)
##    and ordinary lm(): estimate and 95% confidence interval
##  (The conditional tests have no estimation version: they are tests.)
##
##  Design (noise SD 1, effects in noise units):
##    k (predictors)    3, 8
##    n                 30, 60, 200
##    correlation       0, 0.5 (all pairs equal)
##    true coefficients sparse: 0.4, 0.2, and the rest 0
##                      dense:  all 0.15
##  100 datasets per cell: 2,400 Stan fits per contender.
##  Needs cmdstanr (or rstan) and BayesFactor.
##
##  Run from the project folder, in a terminal, and leave it running:
##    Rscript validation/07_estimation.R
##  Roughly 1.5 hours on 8 cores (the last two-prior run took 56 minutes).
##  Each design cell is saved when it finishes
##  (cache/estimation_three_priors/), so an interrupted run
##  resumes where it stopped; a saved cell is reused only if it holds exactly
##  the current method set.
## =====================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
source("validation/settings.R")   # results folder, save_figure()
pkgload::load_all(".", quiet = TRUE)
source("validation/comparison_methods.R")   # contenders, method_scales()

cell_folder <- file.path(cache_folder, "estimation_three_priors")
estimation_methods <- c(names(contenders), "JZS default", "Zellner-Siow", "ordinary (lm)")
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

  ## the three contenders: the same Stan model, R-squared ~ Beta(a, 1)
  stan_fits <- lapply(contenders, function(prior)
    suppressWarnings(fit_regression(formula, data, prior = prior, chains = 4, iter_warmup = 1000,
                                    iter_sampling = 1000, seed = seed)))

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
  stan_summaries <- do.call(rbind, lapply(names(stan_fits), function(method)
    cbind(summarise_draws(as.data.frame(stan_fits[[method]]$draws)[, colnames(predictors), drop = FALSE],
                          method),
          divergent_transitions = stan_fits[[method]]$divergent)))
  rbind(stan_summaries[setdiff(names(stan_summaries), "divergent_transitions")],
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
    mutate(divergent_transitions = stan_summaries$divergent_transitions[
      match(method, stan_summaries$method)])
}


## ---- Run, one design cell at a time (each saved when done) ------------------
started <- Sys.time()
for (row in seq_len(nrow(design))) {
  cell <- design[row, ]
  cell_file <- file.path(cell_folder, sprintf("cell_%02d.rds", row))
  if (file.exists(cell_file) && setequal(unique(readRDS(cell_file)$method), estimation_methods)) next
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

results <- do.call(rbind, lapply(file.path(cell_folder, sprintf("cell_%02d.rds", seq_len(nrow(design)))),
                                 readRDS))


## ---- Summaries -----------------------------------------------------------------
results <- results |>
  mutate(error = estimate - true_value,
         covered = lower <= true_value & true_value <= upper,
         interval_width = upper - lower,
         kind = ifelse(true_value == 0, "true zero", "true effect"))

for (stan_method in names(contenders)) {
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
saveRDS(per_cell, file.path(cache_folder, "estimation_summaries.rds"))

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
cells_facet <- facet_grid(
  number_of_predictors + predictor_correlation ~ true_coefficients,
  labeller = labeller(number_of_predictors  = function(v) paste("k =", v),
                      predictor_correlation = function(v) paste("cor =", v)))

save_figure(
  ggplot(relative_error, aes(sample_size, relative_to_ordinary, colour = method)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = c(30, 60, 200)) + method_scales(estimation_methods, linetype = FALSE) +
    cells_facet + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "total squared error / ordinary least squares",
         colour = NULL, title = "Does the shrinkage pay off?",
         subtitle = "below 1: closer to the true coefficients than ordinary regression"),
  "07_relative_error", width = 9, height = 10)

save_figure(
  ggplot(relative_prediction_error, aes(sample_size, relative_to_ordinary, colour = method)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = c(30, 60, 200)) + method_scales(estimation_methods, linetype = FALSE) +
    cells_facet + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "excess prediction error / ordinary least squares",
         colour = NULL, title = "Does the shrinkage pay off for prediction?",
         subtitle = "below 1: predicts new observations better than ordinary regression"),
  "07_prediction_error", width = 9, height = 10)

save_figure(
  ggplot(per_cell, aes(sample_size, coverage, colour = method, linetype = kind)) +
    geom_hline(yintercept = 0.95, linetype = "dotted") +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = c(30, 60, 200)) + method_scales(estimation_methods, linetype = FALSE) +
    cells_facet + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "coverage of the 95% interval", colour = NULL,
         linetype = NULL, title = "Do the intervals contain the true coefficient?"),
  "07_coverage", width = 9, height = 10)

save_figure(
  ggplot(per_cell, aes(sample_size, bias, colour = method, linetype = kind)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = c(30, 60, 200)) + method_scales(estimation_methods, linetype = FALSE) +
    cells_facet + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "mean of (estimate - truth)", colour = NULL,
         linetype = NULL, title = "Bias: how far estimates are pulled toward zero"),
  "07_bias", width = 9, height = 10)
