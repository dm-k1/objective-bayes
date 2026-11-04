## =====================================================================
##  Study 4: is the Bayes factor calibrated?
##  (validation/ADEMP.md, study 4; replaces part 1 of the old study 10)
##
##  Aims. On data drawn from each prior, the share of true nulls among
##  datasets with a given BF10 equals the posterior probability of the null,
##  1 / (1 + BF10), with equal prior odds (Schad et al., 2023); on data drawn
##  from another prior (JZS), it does not.
##  Data. Set B (validation/data_set_B.R): k = 3 and 8, n = 50 and 200;
##  generating priors the three and JZS; half of each cell's datasets have x1
##  null (odd replications), half draw it from the prior.
##  Estimands. The observed share of nulls per bin of log BF10, against the
##  mean predicted share; the logistic calibration line of "null" on log BF01
##  (intercept 0, slope 1 when calibrated).
##  Methods. The conditional test for x1 under each of the three priors: on
##  data from a prior, that prior's test; on JZS data, all three.
##  Performance. Per bin, the observed share with its Monte Carlo SD under
##  calibration; a chi-square over the bins (each bin's z = (observed -
##  predicted) / Monte Carlo SD under calibration); the logistic slope with its
##  SE. (The SBC package's miscalibration test of Dimitriadis et al. needs the
##  package monotone, which is not in renv.lock; its Brier-score test is
##  one-sided and misses underconfidence, so neither is used.)
##  n_sim. 20,000 datasets per cell (validation/precision.R, bf_calibration):
##  bins holding at least 625 datasets have shares with MCSE <= .02.
##
##  Checks:
##    4.1 (old P1) Each prior's conditional test is calibrated on data from its
##        own prior: the chi-square over the bins does not reject (p > .01) and
##        the logistic slope is within 3 SE of 1, in every cell.
##    4.2 (old P2) On JZS data, each test's miscalibration shows (p < .001).
##
##  Run from the project folder:  Rscript validation/R/04_bayes_factor_calibration.R
##  Smoke test: VALIDATION_SMOKE=1 Rscript validation/R/04_bayes_factor_calibration.R
##  Raw results are cached in cache/study_04/ (delete them to rerun); the
##  tables and checks are saved to validation/results/study_04.rds.
## =====================================================================
library(SimDesign)
library(dplyr)
library(tidyr)
library(ggplot2)
source("validation/R/settings.R")             # flatcor, cores, record_check(), save_figure()
source("validation/R/comparison_methods.R")   # contenders
source("validation/R/precision.R")            # n_sim_of(), mcse_of_proportion()
source("validation/R/data_sets.R")            # with_stream_of()
source("validation/R/data_set_B.R")           # generate_B(), sample_partial_r()

study_folder <- file.path(cache_folder, if (smoke_mode) "study_04_smoke" else "study_04")
dir.create(study_folder, recursive = TRUE, showWarnings = FALSE)

design <- createDesign(generating_prior = c(unname(contenders), "JZS"), k = c(3, 8), n = c(50, 200))
if (smoke_mode) design <- subset(design, k == 3 & n == 50)
n_sim <- if (smoke_mode) 200L else n_sim_of("bf_calibration")   # calibration needs a few hundred to bin


## ---- The simulation ------------------------------------------------------------
simulate_data <- function(condition, fixed_objects) {
  r <- condition$REPLICATION
  generate_B(condition$generating_prior, condition$k, condition$n,
             x1 = if (r %% 2 == 1) "null" else "prior", purpose = "calibration", replication = r)
}

## log BF10 for x1 under each prior's conditional test, and whether x1 was null
analyse_dataset <- function(condition, dat, fixed_objects) {
  partial_r <- sample_partial_r(dat$predictors, dat$outcome)
  log_bfs <- vapply(unname(contenders), function(prior)
    log_bayes_factor_conditional(partial_r, condition$n, condition$k, prior = prior), numeric(1))
  c(is_null = as.numeric(dat$coefficients[1] == 0),
    true_partial_correlation = dat$partial_correlation_x1, setNames(log_bfs, unname(contenders)))
}

summarise_cell <- function(condition, results, fixed_objects) c(share_null = mean(results$is_null))

results_file <- file.path(study_folder, "simulation.rds")
if (file.exists(results_file)) {
  simulation <- readRDS(results_file)
  message("Using cached results: ", results_file)
} else {
  simulation <- runSimulation(design, replications = n_sim, generate = simulate_data,
                              analyse = analyse_dataset, summarise = summarise_cell,
                              parallel = number_of_cores > 1, ncores = number_of_cores,
                              packages = c("flatcor"), store_results = TRUE,
                              control = list(include_replication_index = TRUE), progress = FALSE,
                              save_details = list(out_rootdir = study_folder, tmpfilename = "resume_simulation.qs2"))
  saveRDS(simulation, results_file)
}

## one row per dataset and test: on data from a prior only that prior's test,
## on JZS data all three
calibration_runs <- SimExtract(simulation, what = "results") |>
  pivot_longer(all_of(unname(contenders)), names_to = "test", values_to = "log_BF10") |>
  filter(generating_prior == "JZS" | generating_prior == test) |>
  mutate(is_null = is_null == 1,
         predicted_null = plogis(-log_BF10),          # 1 / (1 + BF10)
         cell = sprintf("k = %d, n = %d", k, n))


## ---- Per bin: observed share of nulls against the mean of 1 / (1 + BF10) --------
bin_edges <- c(-Inf, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 6, Inf)
calibration_bins <- calibration_runs |>
  mutate(bin = cut(log_BF10, bin_edges)) |>
  group_by(generating_prior, test, cell, bin) |>
  summarise(datasets = n(), observed = mean(is_null), predicted = mean(predicted_null),
            ## Monte Carlo SD of the observed share if calibrated (Bernoulli(p_i) draws)
            mc_sd = sqrt(sum(predicted_null * (1 - predicted_null))) / n(),
            mcse_observed = mcse_of_proportion(mean(is_null), n()),
            .groups = "drop") |>
  mutate(z = (observed - predicted) / mc_sd)

## Per cell: the chi-square over the bins with at least 5 expected nulls and 5
## expected effects, and the logistic calibration line
chi_square_tests <- calibration_bins |>
  filter(datasets * predicted >= 5, datasets * (1 - predicted) >= 5) |>
  group_by(generating_prior, test, cell) |>
  summarise(bins = n(), chi_square = sum(z^2), largest_abs_z = max(abs(z)), .groups = "drop") |>
  mutate(p_value = pchisq(chi_square, bins, lower.tail = FALSE))
logistic_lines <- calibration_runs |>
  group_by(generating_prior, test, cell) |>
  group_modify(function(rows, key) {
    ## (datasets with huge Bayes factors give fitted probabilities of 0 or 1: harmless)
    fit <- suppressWarnings(glm(is_null ~ I(-log_BF10), family = binomial, data = rows))
    estimates <- summary(fit)$coefficients
    tibble(datasets = nrow(rows), intercept = estimates[1, 1], intercept_se = estimates[1, 2],
           slope = estimates[2, 1], slope_se = estimates[2, 2])
  }) |>
  ungroup()
calibration_tests <- left_join(chi_square_tests, logistic_lines, by = c("generating_prior", "test", "cell"))

cat("\nCalibration per bin (share of nulls: observed with its MCSE, predicted = mean 1/(1 + BF10)):\n")
print(as.data.frame(calibration_bins |> select(generating_prior, test, cell, bin, datasets,
                                               observed, mcse_observed, predicted)), digits = 3)
cat("\nChi-square over the bins and the logistic calibration line per cell:\n")
print(as.data.frame(calibration_tests), digits = 3)


## ---- Checks --------------------------------------------------------------------
for (prior in unname(contenders)) for (cell_name in unique(calibration_tests$cell)) {
  row <- filter(calibration_tests, generating_prior == prior, test == prior, cell == cell_name)
  record_check(sprintf("4.1 %s: calibrated on its own prior, %s", prior, cell_name),
               row$p_value > 0.01 && abs(row$slope - 1) < 3 * row$slope_se,
               sprintf("chi2 = %.1f on %d bins (p = %.3f; %d datasets); slope %.3f (SE %.3f)",
                       row$chi_square, row$bins, row$p_value, row$datasets, row$slope, row$slope_se))
}
for (prior in unname(contenders)) for (cell_name in unique(calibration_tests$cell)) {
  row <- filter(calibration_tests, generating_prior == "JZS", test == prior, cell == cell_name)
  record_check(sprintf("4.2 %s, data from JZS: miscalibration shows, %s", prior, cell_name),
               row$p_value < 0.001,
               sprintf("chi2 = %.0f on %d bins (p = %.2g); slope %.2f (SE %.2f)",
                       row$chi_square, row$bins, row$p_value, row$slope, row$slope_se))
}

## What the effects called null are, on the test's own prior
called_null <- calibration_runs |>
  filter(generating_prior == test, log_BF10 < -log(3)) |>
  group_by(test, cell) |>
  summarise(called_null = n(), share_true_nulls = mean(is_null),
            mcse = mcse_of_proportion(mean(is_null), n()),
            effects_with_abs_rho_below_0.1 = mean(abs(true_partial_correlation[!is_null]) < 0.1),
            .groups = "drop")
cat("\nDatasets with BF01 > 3, each prior's test on data from its own prior:\n")
print(as.data.frame(called_null), digits = 3)


## ---- Figure --------------------------------------------------------------------
save_figure(
  ggplot(filter(calibration_bins, datasets >= 20) |>
           mutate(data_from = ifelse(generating_prior == "JZS", "JZS", "its own prior")),
         aes(predicted, observed, colour = data_from)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = predicted - 2 * mc_sd, ymax = predicted + 2 * mc_sd),
                  width = 0, colour = "grey70") +
    geom_point(size = 2) + geom_line() +
    facet_grid(test ~ cell) +
    scale_colour_manual(values = c(`its own prior` = "#0072B2", JZS = "grey40")) +
    coord_equal() + theme(legend.position = "bottom") +
    labs(x = "predicted share of nulls, mean of 1 / (1 + BF10)", y = "observed share of nulls",
         colour = "data drawn from",
         title = "Calibration of the conditional test (bins of log BF10, grey: 2 Monte Carlo SD)",
         subtitle = sprintf("%d datasets per cell; rows: the prior whose conditional test is used", n_sim)),
  "study_04_calibration", width = 12, height = 10)

save_study_results("04", list(calibration_bins = calibration_bins, calibration_tests = calibration_tests,
                              called_null = called_null))
report_checks()
