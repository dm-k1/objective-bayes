## =====================================================================
##  Validation 09: variable selection with many candidate predictors
##
##  Where the prior's claim about model size should matter: comparing models
##  of DIFFERENT sizes. Twelve candidate predictors, four of which matter.
##  Every subset (2^12 = 4,096 models) gets a Bayes factor against the
##  intercept-only model; with a prior over models that gives posterior model
##  probabilities and inclusion probabilities per predictor.
##
##  Design (noise SD 1, effects in noise units):
##    true coefficients   0.3, 0.2, 0.15, 0.1, and 8 zeros
##    correlation         0, 0.5 (all pairs equal)
##    n                   50, 100, 200, 500
##  100 datasets per cell.
##
##  Bayes factors per model (each depends only on R-squared, n and k):
##    flat-R, JZS default (scale 0.354), Zellner-Siow (scale 1), and BIC
##    (exp(-BIC/2) relative to the null) as a non-Bayesian-prior reference.
##  Priors over models:
##    uniform      every subset equally likely (favours models with ~6 predictors)
##    beta-binomial every model SIZE equally likely (Scott & Berger, 2010)
##
##  Per dataset: is the highest-probability model the true one; the median
##  probability model (predictors with inclusion probability > 0.5); number of
##  false inclusions and missed predictors; inclusion probabilities.
##
##  For speed, each method's log Bayes factor is tabulated once per (n, k) on a
##  grid of R-squared values and interpolated (check V1 below: the error is
##  reported against direct computation).
##
##  Run from the project folder: Rscript validation/09_variable_selection.R
##  About 10 minutes; cached.
## =====================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
source("validation/settings.R")   # parallel_map(), caching, record_check(), save_figure()
pkgload::load_all(".", quiet = TRUE)

number_of_candidates <- 12
true_slopes <- c(0.3, 0.2, 0.15, 0.1, rep(0, 8))
true_model  <- true_slopes != 0
design <- expand.grid(sample_size = c(50, 100, 200, 500), predictor_correlation = c(0, 0.5))
replications_per_cell <- 100
methods <- c("flat-R", "JZS default", "Zellner-Siow (1)", "BIC")

## all 4,096 subsets as a logical matrix, one row per model
all_models <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), number_of_candidates)))
model_sizes <- rowSums(all_models)


## ---- Log Bayes factor against the null, tabulated in R-squared ---------------
log_bf_direct <- function(method, r_squared, n, k) {
  if (k == 0) return(0)
  switch(method,
    `flat-R`           = log_bayes_factor_flat_r_vs_intercept_only(r_squared, n, k),
    ## simple = FALSE returns the log Bayes factor (simple = TRUE overflows to Inf)
    `JZS default`      = BayesFactor::linearReg.R2stat(n, k, r_squared, rscale = sqrt(2) / 4, simple = FALSE)$bf,
    `Zellner-Siow (1)` = BayesFactor::linearReg.R2stat(n, k, r_squared, rscale = 1, simple = FALSE)$bf,
    BIC                = -(n * log(1 - r_squared) + k * log(n)) / 2)
}
## grid in log-odds of R-squared, where log BF is smooth
log_odds_grid <- seq(-14, 9, length.out = 300)
build_tables <- function(n) lapply(setNames(methods, methods), function(method)
  lapply(0:number_of_candidates, function(k) {
    if (k == 0 || method == "BIC") return(NULL)
    r2 <- plogis(log_odds_grid)
    splinefun(log_odds_grid, sapply(r2, log_bf_direct, method = method, n = n, k = k))
  }))
log_bf_table <- function(tables, method, r_squared, n, k) {
  if (k == 0) return(0)
  if (method == "BIC") return(log_bf_direct("BIC", r_squared, n, k))
  tables[[method]][[k + 1]](pmin(pmax(qlogis(r_squared), min(log_odds_grid)), max(log_odds_grid)))
}

message("Tabulating Bayes factors for each n ...")
tables_per_n <- parallel_map(design$sample_size[!duplicated(design$sample_size)], build_tables)
names(tables_per_n) <- unique(design$sample_size)

## V1. interpolation error against direct computation, random R-squared values
set.seed(1)
interpolation_errors <- sapply(1:200, function(i) {
  sizes <- unique(design$sample_size)
  n <- sizes[sample.int(length(sizes), 1)]; k <- sample.int(number_of_candidates, 1)
  method <- sample(c("flat-R", "JZS default", "Zellner-Siow (1)"), 1); r2 <- runif(1, 0.001, 0.9)
  abs(log_bf_table(tables_per_n[[as.character(n)]], method, r2, n, k) - log_bf_direct(method, r2, n, k))
})
record_check("V1. tabulated log BF vs direct (200 random cases)", max(interpolation_errors) < 1e-3,
             sprintf("largest error %.1e", max(interpolation_errors)))


## ---- One dataset ------------------------------------------------------------
log_model_prior <- list(
  uniform        = rep(0, nrow(all_models)),
  `beta-binomial` = -lchoose(number_of_candidates, model_sizes))   # each size gets 1/13

analyse_dataset <- function(n, correlation, seed) {
  set.seed(seed)
  covariance <- matrix(correlation, number_of_candidates, number_of_candidates); diag(covariance) <- 1
  predictors <- matrix(rnorm(n * number_of_candidates), n) %*% chol(covariance)
  outcome <- drop(predictors %*% true_slopes) + rnorm(n)
  r_squared <- apply(all_models, 1, function(in_model)
    if (!any(in_model)) 0 else sample_r_squared(predictors[, in_model, drop = FALSE], outcome))
  tables <- tables_per_n[[as.character(n)]]
  do.call(rbind, lapply(methods, function(method) {
    log_bf <- mapply(function(r2, k) log_bf_table(tables, method, r2, n, k), r_squared, model_sizes)
    do.call(rbind, lapply(names(log_model_prior), function(model_prior) {
      log_posterior <- log_bf + log_model_prior[[model_prior]]
      posterior <- exp(log_posterior - max(log_posterior)); posterior <- posterior / sum(posterior)
      inclusion <- colSums(all_models * posterior)
      best_model <- all_models[which.max(posterior), ]
      median_model <- inclusion > 0.5
      data.frame(method = method, model_prior = model_prior,
                 best_is_true          = all(best_model == true_model),
                 median_model_is_true  = all(median_model == true_model),
                 false_inclusions      = sum(median_model & !true_model),
                 missed                = sum(!median_model & true_model),
                 size_of_median_model  = sum(median_model),
                 posterior_of_true     = posterior[which(apply(all_models, 1, function(m) all(m == true_model)))],
                 inclusion_true_mean   = mean(inclusion[true_model]),
                 inclusion_null_mean   = mean(inclusion[!true_model]),
                 inclusion_x4          = inclusion[4])   # the smallest real effect, 0.1
    }))
  }))
}

runs <- run_simulation_grid(
  design,
  results_file = file.path(results_folder, "sim_variable_selection.rds"),
  settings = list(replications_per_cell = replications_per_cell, true_slopes = true_slopes, version = 1),
  simulate_cell = function(cell)
    do.call(rbind, lapply(seq_len(replications_per_cell), function(replication)
      cbind(as.data.frame(cell), replication = replication,
            analyse_dataset(cell$sample_size, cell$predictor_correlation,
                            seed = 1000 * cell$sample_size + 10 * replication +
                                   round(10 * cell$predictor_correlation))))))


## ---- Summaries -----------------------------------------------------------------
summaries <- runs |>
  group_by(predictor_correlation, sample_size, model_prior, method) |>
  summarise(best_is_true         = mean(best_is_true),
            median_model_is_true = mean(median_model_is_true),
            false_inclusions     = mean(false_inclusions),
            missed               = mean(missed),
            size_of_median_model = mean(size_of_median_model),
            inclusion_true       = mean(inclusion_true_mean),
            inclusion_null       = mean(inclusion_null_mean),
            inclusion_x4         = mean(inclusion_x4), .groups = "drop")
print(as.data.frame(summaries), digits = 3)
saveRDS(summaries, file.path(results_folder, "variable_selection_summaries.rds"))

method_colours <- c(`flat-R` = "#d95f02", `JZS default` = "grey40",
                    `Zellner-Siow (1)` = "#7570b3", BIC = "#1b9e77")
long <- summaries |>
  pivot_longer(c(median_model_is_true, false_inclusions, missed),
               names_to = "measure", values_to = "value") |>
  mutate(measure = recode(measure,
    median_model_is_true = "median probability model = true model (share)",
    false_inclusions     = "false inclusions (of 8 null predictors)",
    missed               = "missed predictors (of 4 real ones)"))

save_figure(
  ggplot(long, aes(sample_size, value, colour = method, linetype = factor(predictor_correlation))) +
    geom_line(linewidth = 0.8) + geom_point() +
    scale_x_log10(breaks = c(50, 100, 200, 500)) + scale_colour_manual(values = method_colours) +
    facet_grid(measure ~ model_prior, scales = "free_y") +
    theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = NULL, colour = NULL, linetype = "correlation",
         title = "Variable selection: 12 candidates, 4 real (0.3, 0.2, 0.15, 0.1)",
         subtitle = "columns: prior over models"),
  "09_variable_selection", width = 10, height = 9)

save_figure(
  ggplot(summaries |> pivot_longer(c(inclusion_true, inclusion_null), names_to = "kind", values_to = "value") |>
           mutate(kind = recode(kind, inclusion_true = "real predictors", inclusion_null = "null predictors")),
         aes(sample_size, value, colour = method, linetype = kind)) +
    geom_line(linewidth = 0.8) + geom_point() +
    scale_x_log10(breaks = c(50, 100, 200, 500)) + scale_colour_manual(values = method_colours) +
    facet_grid(predictor_correlation ~ model_prior,
               labeller = labeller(predictor_correlation = function(v) paste("correlation =", v))) +
    theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "mean inclusion probability", colour = NULL, linetype = NULL,
         title = "Inclusion probabilities: real predictors should go up, null ones down"),
  "09_inclusion", width = 10, height = 7)

report_checks()
