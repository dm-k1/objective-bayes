## =====================================================================
##  Study 9: how good are the estimates?
##  (validation/ADEMP.md, study 9; replaces old study 07_estimation)
##
##  Aims. Bias, interval coverage and prediction error of the Stan posterior
##  under each prior, against ordinary least squares and the BayesFactor
##  posteriors (JZS, Zellner-Siow), on the same datasets as study 7.
##  Data. A subset of set A (validation/data_sets.R, generate_A, so the
##  datasets are exactly study 7's for these cells): k in {2, 10}, n in
##  {50, 200}, predictor correlation in {0, .7}, the other slopes sparse or
##  dense, slope of x1 0.2. That is 16 cells.
##  Estimands, per method and cell:
##    - the bias of x1's estimate (posterior mean), as a share of its true
##      value 0.2;
##    - the coverage of x1's 95% interval (one coefficient per dataset, so the
##      replications are independent);
##    - the excess prediction error (b_hat - b)' Sigma (b_hat - b), Sigma the
##      population predictor covariance, relative to OLS.
##  Methods. flatcor's fit_regression() with each prior (Stan; the seed of
##  every fit follows from its dataset's ID); BayesFactor posterior draws for
##  JZS (scale 0.354) and Zellner-Siow (scale 1); OLS (lm, confint).
##  Performance. Relative bias and coverage with their Monte Carlo SE from
##  rsimsum (the posterior SD, or OLS's standard error, as the model SE);
##  the prediction-error ratio with a paired bootstrap SE over datasets.
##  n_sim. 475 datasets per cell (validation/precision.R, coverage_per_cell):
##  coverage near .95 then has MCSE <= .01.
##
##  Checks, set relative to the Monte Carlo error:
##    9.1 Flat partial r's 95% intervals for x1 cover near nominally: in every
##        cell the coverage is above .93 or within 2 MCSE of .95.
##    9.2 Flat partial r predicts at least as well as OLS: in every cell its
##        prediction-error ratio is below 1 + 2 MCSE.
##    9.3 The Stan model samples without divergent transitions: none over all
##        fits of the three priors.
##  The old study had no PASS/FAIL checks; 9.1 and 9.2 test what its docs
##  reported (flat partial r's coverage of at least .916, better than OLS in
##  every cell), 9.3 its count of divergent transitions.
##
##  Run from the project folder:  Rscript validation/09_estimation.R
##  Smoke test (2 cells, 10 datasets, short chains):
##    VALIDATION_SMOKE=1 Rscript validation/09_estimation.R
##  Results are cached in cache/study_09/ (delete them to rerun). Each fit of
##  fit_regression() runs its 4 chains in parallel, so SimDesign uses one
##  worker per 4 cores.
## =====================================================================
library(SimDesign)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rsimsum)
source("validation/settings.R")             # flatcor, cores, record_check(), save_figure()
source("validation/comparison_methods.R")   # contenders, method_scales()
source("validation/precision.R")            # n_sim_of()
source("validation/data_sets.R")            # generate_A(), cell_id_A(), slopes_A(), with_stream_of()

study_folder <- file.path(cache_folder, if (smoke_mode) "study_09_smoke" else "study_09")
dir.create(study_folder, recursive = TRUE, showWarnings = FALSE)

design <- expand.grid(k = c(2, 10), n = c(50, 200), rho = c(0, 0.7),
                      others = c("sparse", "dense"), slope_x1 = 0.2, stringsAsFactors = FALSE)
if (smoke_mode) design <- subset(design, k == 2 & n == 50 & rho == 0)
n_sim <- n_sim_of("coverage_per_cell")
true_slope_x1 <- 0.2

## Stan settings: 4 chains of 1,000 draws; short chains in smoke mode
stan_chains   <- if (smoke_mode) 2 else 4
stan_warmup   <- if (smoke_mode) 200 else 1000
stan_sampling <- if (smoke_mode) 200 else 1000
bayesfactor_iterations <- if (smoke_mode) 1000 else 4000
workers <- max(1, number_of_cores %/% stan_chains)

estimation_methods <- c(names(contenders), "JZS default", "Zellner-Siow", "ordinary (lm)")
method_key <- setNames(make.names(estimation_methods), estimation_methods)


## ---- The simulation ------------------------------------------------------------
simulate_data <- function(condition, fixed_objects) generate_A(condition, condition$REPLICATION)

analyse_dataset <- function(condition, dat, fixed_objects) {
  dataset_id <- paste(cell_id_A(condition), condition$REPLICATION, sep = "-r")
  k <- condition$k
  covariance  <- matrix(condition$rho, k, k); diag(covariance) <- 1
  true_slopes <- slopes_A(k, condition$others, condition$slope_x1)
  data <- data.frame(dat$predictors, y = dat$outcome)
  formula <- reformulate(colnames(dat$predictors), response = "y")
  ## every fit's seed follows from its dataset's ID, so results do not depend
  ## on the number of cores or the order of the fits
  seed_of <- function(label) abs(digest::digest2int(paste(dataset_id, label))) %% .Machine$integer.max

  summary_of <- function(draws) {      # posterior mean and 95% interval of x1, and the excess prediction error
    estimates <- colMeans(draws)
    error <- estimates - true_slopes
    c(estimate = estimates[["x1"]], sd = sd(draws[, "x1"]),
      lower = quantile(draws[, "x1"], 0.025, names = FALSE),
      upper = quantile(draws[, "x1"], 0.975, names = FALSE),
      prediction_error = drop(t(error) %*% covariance %*% error))
  }
  per_method <- list()
  divergent  <- 0
  for (method in names(contenders)) {
    fit <- suppressWarnings(fit_regression(formula, data, prior = contenders[[method]],
                                           chains = stan_chains, iter_warmup = stan_warmup,
                                           iter_sampling = stan_sampling, seed = seed_of(method)))
    divergent <- divergent + fit$divergent
    per_method[[method]] <- summary_of(as.matrix(as.data.frame(fit$draws)[, colnames(dat$predictors),
                                                                         drop = FALSE]))
  }
  bayesfactor_draws <- function(scale, label)
    with_stream_of(paste(dataset_id, label), as.matrix(as.data.frame(
      BayesFactor::lmBF(formula, data, rscaleCont = scale, posterior = TRUE,
                        iterations = bayesfactor_iterations, progress = FALSE))[, colnames(dat$predictors),
                                                                               drop = FALSE]))
  per_method[["JZS default"]]  <- summary_of(bayesfactor_draws(sqrt(2) / 4, "JZS"))
  per_method[["Zellner-Siow"]] <- summary_of(bayesfactor_draws(1, "ZS"))
  ols_fit <- lm(formula, data)
  ols_error <- coef(ols_fit)[colnames(dat$predictors)] - true_slopes
  per_method[["ordinary (lm)"]] <- c(estimate = coef(ols_fit)[["x1"]],
                                     sd = summary(ols_fit)$coefficients["x1", "Std. Error"],
                                     lower = confint(ols_fit)["x1", 1], upper = confint(ols_fit)["x1", 2],
                                     prediction_error = drop(t(ols_error) %*% covariance %*% ols_error))
  values <- unlist(lapply(estimation_methods, function(method)
    setNames(per_method[[method]], paste(names(per_method[[method]]), method_key[[method]], sep = "."))))
  c(values, divergent = divergent)
}

## per cell: the divergent transitions of the Stan fits (the estimands are
## computed from the stored results below)
summarise_cell <- function(condition, results, fixed_objects) c(divergent = sum(results[, "divergent"]))

results_file <- file.path(study_folder, "simulation.rds")
if (file.exists(results_file)) {
  simulation <- readRDS(results_file)
  message("Using cached results: ", results_file)
} else {
  simulation <- runSimulation(design, replications = n_sim, generate = simulate_data,
                              analyse = analyse_dataset, summarise = summarise_cell,
                              parallel = workers > 1, ncores = workers,
                              packages = c("flatcor", "BayesFactor"), store_results = TRUE,
                              ## the resume file of an interrupted run, in this study's cache
                              save_details = list(out_rootdir = study_folder,
                                                  tmpfilename = "resume_simulation.qs2"),
                              control = list(include_replication_index = TRUE), progress = FALSE)
  saveRDS(simulation, results_file)
}

## one row per dataset and method: x1's estimate and interval, and the excess
## prediction error
per_dataset <- SimExtract(simulation, what = "results")
estimates <- per_dataset |>
  select(-divergent) |>
  pivot_longer(matches("^(estimate|sd|lower|upper|prediction_error)[.]"),
               names_to = c(".value", "key"), names_pattern = "^([a-z_]+)[.](.*)$") |>
  mutate(method = names(method_key)[match(key, method_key)],
         cell = sprintf("k = %d, n = %d, rho = %.1f, %s", k, n, rho, others))


## ---- Bias and coverage of x1, with their Monte Carlo SE (rsimsum) -------------
x1_summary <- simsum(as.data.frame(rename(estimates, ci_lower = lower, ci_upper = upper)),
                     estvarname = "estimate", se = "sd", true = true_slope_x1, methodvar = "method",
                     by = c("k", "n", "rho", "others"), ci.limits = c("ci_lower", "ci_upper"))
## rbias: the bias as a share of the true value 0.2 (rsimsum's relative bias)
bias_and_coverage <- tidy(summary(x1_summary, stats = c("rbias", "cover"))) |>
  mutate(stat = recode(stat, rbias = "bias"),
         across(c(k, n, rho), ~ as.numeric(as.character(.x))))   # rsimsum returns the by-variables as factors

cat("\nBias of x1 (share of its true value 0.2) and coverage of its 95% interval, with MCSE:\n")
print(as.data.frame(bias_and_coverage |>
  mutate(value = sprintf("%.3f (%.3f)", est, mcse)) |>
  select(stat, k, n, rho, others, method, value) |>
  pivot_wider(names_from = method, values_from = value)), row.names = FALSE)


## ---- Excess prediction error relative to OLS, with a paired bootstrap SE --------
bootstrap_replicates <- if (smoke_mode) 100 else 1000
RNGkind("L'Ecuyer-CMRG"); set.seed(9)   # the bootstrap's own, fixed stream
prediction_ratio <- estimates |>
  group_by(k, n, rho, others) |>
  group_modify(function(cell, keys) {
    error_of <- function(method_name)                       # in replication order
      cell$prediction_error[cell$method == method_name][order(cell$REPLICATION[cell$method == method_name])]
    replications <- sum(cell$method == "ordinary (lm)")
    resamples <- replicate(bootstrap_replicates, sample(replications, replace = TRUE), simplify = FALSE)
    ols <- error_of("ordinary (lm)")
    bind_rows(lapply(setdiff(estimation_methods, "ordinary (lm)"), function(method_name) {
      errors <- error_of(method_name)
      tibble(method = method_name, ratio = mean(errors) / mean(ols),
             mcse = sd(vapply(resamples, function(rows) mean(errors[rows]) / mean(ols[rows]), numeric(1))))
    }))
  }) |>
  ungroup()

cat("\nExcess prediction error relative to OLS (below 1: better than OLS), with MCSE:\n")
print(as.data.frame(prediction_ratio |>
  mutate(value = sprintf("%.3f (%.3f)", ratio, mcse)) |> select(-ratio, -mcse) |>
  pivot_wider(names_from = method, values_from = value)), row.names = FALSE)


## ---- Checks --------------------------------------------------------------------
coverage_partial_r <- filter(bias_and_coverage, stat == "cover", method == "flat partial r")
record_check("9.1 flat partial r's 95% intervals for x1 cover near nominally",
             all(coverage_partial_r$est > 0.93 | abs(coverage_partial_r$est - 0.95) <= 2 * coverage_partial_r$mcse),
             sprintf("lowest coverage %.3f (MCSE %.3f) over %d cells",
                     min(coverage_partial_r$est), coverage_partial_r$mcse[which.min(coverage_partial_r$est)],
                     nrow(coverage_partial_r)))

ratio_partial_r <- filter(prediction_ratio, method == "flat partial r")
record_check("9.2 flat partial r predicts at least as well as OLS",
             all(ratio_partial_r$ratio < 1 + 2 * ratio_partial_r$mcse),
             sprintf("highest ratio to OLS %.3f (MCSE %.3f) over %d cells",
                     max(ratio_partial_r$ratio), ratio_partial_r$mcse[which.max(ratio_partial_r$ratio)],
                     nrow(ratio_partial_r)))

record_check("9.3 the Stan model samples without divergent transitions",
             sum(per_dataset$divergent) == 0,
             sprintf("%d divergent transitions over %d fits (%d datasets x 3 priors)",
                     sum(per_dataset$divergent), 3 * nrow(per_dataset), nrow(per_dataset)))


## ---- Figures -------------------------------------------------------------------
cells_facet <- facet_grid(others + rho ~ k,
  labeller = labeller(k = function(v) paste("k =", v), rho = function(v) paste("correlation =", v)))
coverage_plot <- filter(bias_and_coverage, stat == "cover")
save_figure(
  ggplot(coverage_plot, aes(n, est, colour = method)) +
    geom_hline(yintercept = 0.95, linetype = "dotted") +
    geom_pointrange(aes(ymin = est - 2 * mcse, ymax = est + 2 * mcse),
                    position = position_dodge(width = 0.15), size = 0.25, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200)) + method_scales(estimation_methods, linetype = FALSE) +
    cells_facet +
    labs(x = "sample size (log scale)", y = "coverage of x1's 95% interval",
         title = "Do the intervals contain the true slope of x1 (0.2)?",
         subtitle = sprintf("%d datasets per cell; bars: 2 Monte Carlo SE", n_sim)),
  "study_09_coverage", width = 9, height = 10)

bias_plot <- filter(bias_and_coverage, stat == "bias")
save_figure(
  ggplot(bias_plot, aes(n, est, colour = method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_pointrange(aes(ymin = est - 2 * mcse, ymax = est + 2 * mcse),
                    position = position_dodge(width = 0.15), size = 0.25, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200)) + method_scales(estimation_methods, linetype = FALSE) +
    cells_facet +
    labs(x = "sample size (log scale)", y = "bias of x1, as a share of its true value",
         title = "Bias: how far estimates of x1 are pulled towards zero",
         subtitle = sprintf("%d datasets per cell; bars: 2 Monte Carlo SE", n_sim)),
  "study_09_bias", width = 9, height = 10)

save_figure(
  ggplot(prediction_ratio, aes(n, ratio, colour = method)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_pointrange(aes(ymin = ratio - 2 * mcse, ymax = ratio + 2 * mcse),
                    position = position_dodge(width = 0.15), size = 0.25, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200)) +
    method_scales(setdiff(estimation_methods, "ordinary (lm)"), linetype = FALSE) +
    cells_facet +
    labs(x = "sample size (log scale)", y = "excess prediction error / OLS",
         title = "Does the shrinkage pay off for prediction?",
         subtitle = "below 1: predicts new observations better than ordinary regression; bars: 2 Monte Carlo SE"),
  "study_09_prediction_error", width = 9, height = 10)

saveRDS(list(bias_and_coverage = bias_and_coverage, prediction_ratio = prediction_ratio,
             divergent = sum(per_dataset$divergent)),
        file.path(study_folder, "summaries.rds"))
report_checks()
