## =====================================================================
##  Study 10: how does each prior select variables?
##  (docs/appendix/C_study_design.md, study 10; replaces old study 09)
##
##  Aims. Variable selection over all 2^12 = 4,096 subsets of twelve
##  candidate predictors, of which four matter: how many real predictors each
##  method misses and how many null ones it includes, and whether methods
##  differ in how well they rank real above null predictors or only in where
##  they cut.
##  Data. Set C (R/data/data_set_C.R): slopes 0.3, 0.2, 0.15, 0.1 and
##  eight zeros; n in {50, 100, 200, 500}; predictor correlation 0 or .5.
##  Estimands, per method, prior over models and cell: the mean number of
##  real predictors missed (of 4) and of null predictors included (of 8) by
##  the median probability model (inclusion probability > .5); the share of
##  datasets where the median probability model and the highest-probability
##  model are the true one; the mean inclusion probability of each real
##  predictor; the AUC of inclusion probabilities, real against null.
##  Methods. Whole-model Bayes factors against the intercept-only model
##  (R/shared/comparison_methods.R): the three priors (flat partial r, flat
##  R, flat R-squared), JZS (scale 0.354) and Zellner-Siow (scale 1), BIC
##  (exp(-BIC/2) relative to the null), and the hyper-g/n prior of Liang et al.
##  (2008) from BAS (bas.lm, alpha = 3). bas.lm offers no robust or intrinsic
##  prior for linear models (its priors: AIC, BIC, g-prior, JZS, ZS-null,
##  hyper-g, hyper-g-laplace, hyper-g-n, EB-local, EB-global), so those are
##  not included. Priors over models: uniform (every subset equally likely)
##  and beta-binomial(1, 1) (every model size equally likely; Scott & Berger,
##  2010). All methods see the same datasets.
##  Performance. Means with MCSE SD / sqrt(n_sim); shares with
##  mcse_of_proportion(); the AUC difference from flat R with a paired
##  bootstrap SE over datasets.
##  n_sim. 625 datasets per cell (precision.R, rate_per_cell): shares have
##  MCSE <= .02, and mean counts MCSE of about .05 or less.
##
##  Checks (the old V-checks, renumbered; thresholds for KS checks are the
##  critical distance at the .01 level for the number of draws, 1.63 / sqrt(n)):
##    10.1 (old V1) the tabulated log BF equals direct computation (200 random
##         cases, largest error < 1e-3).
##    10.2 (old V2) the Beta(a, 1) g-mixture with a = 1/2 equals the flat R
##         closed form (< 1e-3).
##    10.3 (old V3) the g-mixture with a = k/2 implies R-squared ~ Beta(k/2, 1).
##    10.4 (old V4) uniformly drawn correlation matrices imply R-squared ~
##         Beta(k/2, 1).
##    10.5 (old V5) W ~ Inverse-Gamma(b, 1/2) implies R-squared ~ Beta(1/2, b).
##
##  Run from the project folder:  Rscript R/studies/10_variable_selection.R
##  Smoke test (n = 50, correlation 0, 10 datasets):
##    VALIDATION_SMOKE=1 Rscript R/studies/10_variable_selection.R
##  Raw results are cached in cache/study_10/ (delete them to rerun); the
##  tables and checks are saved to results/study_10.rds.
## =====================================================================
library(SimDesign)
library(BAS)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/shared/settings.R")             # flatcor, cores, record_check(), save_figure()
source("R/shared/comparison_methods.R")   # log_bayes_factor_model(), draw_lkj_correlation_matrix()
source("R/shared/precision.R")            # n_sim_of(), mcse_of_proportion(), smoke_mode
source("R/data/data_sets.R")            # with_stream_of()
source("R/data/data_set_C.R")           # design_C(), generate_C(), true_slopes_C

study_folder <- file.path(cache_folder, if (smoke_mode) "study_10_smoke" else "study_10")
dir.create(study_folder, recursive = TRUE, showWarnings = FALSE)

design <- design_C()
if (smoke_mode) design <- subset(design, n == 50 & rho == 0)
n_sim <- n_sim_of("rate_per_cell")

number_of_candidates <- length(true_slopes_C)
true_model <- true_slopes_C != 0
tabulated_methods <- model_methods                     # Bayes factors from R-squared, n and k
methods <- c(tabulated_methods, "BIC", "hyper-g/n (BAS)")
model_priors <- c("uniform", "beta-binomial")

## all 4,096 subsets as a logical matrix, one row per model
all_models <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), number_of_candidates)))
model_sizes <- rowSums(all_models)
log_model_prior <- list(
  uniform         = rep(0, nrow(all_models)),
  `beta-binomial` = -lchoose(number_of_candidates, model_sizes))   # each size gets 1/13


## ---- Log Bayes factor against the null, tabulated in R-squared ---------------
## Each method's log BF depends only on R-squared, n and k; it is tabulated once
## per (n, method, k) on a grid in the log-odds of R-squared and interpolated.
log_odds_grid <- seq(-14, 9, length.out = 300)
table_tasks <- expand.grid(n = unique(design$n), method = tabulated_methods,
                           k = seq_len(number_of_candidates), stringsAsFactors = FALSE)
message("Tabulating Bayes factors for each n, method and k (", nrow(table_tasks), " tables) ...")
table_values <- parallel_map(seq_len(nrow(table_tasks)), function(row)
  suppressWarnings(sapply(plogis(log_odds_grid), log_bayes_factor_model,
                          method = table_tasks$method[row], n = table_tasks$n[row],
                          k = table_tasks$k[row])))
tables_per_n <- lapply(setNames(nm = as.character(unique(design$n))), function(n)
  lapply(setNames(nm = tabulated_methods), function(method)
    lapply(seq_len(number_of_candidates), function(k)
      splinefun(log_odds_grid, table_values[[which(table_tasks$n == as.numeric(n) &
                                                     table_tasks$method == method &
                                                     table_tasks$k == k)]]))))
log_bf_tabulated <- function(tables, method, r_squared, k)
  ifelse(k == 0, 0,
         vapply(seq_along(r_squared), function(i) if (k[i] == 0) 0 else
           tables[[method]][[k[i]]](pmin(pmax(qlogis(r_squared[i]), min(log_odds_grid)),
                                         max(log_odds_grid))), numeric(1)))


## ---- Checks of the ingredients (old V1-V5) --------------------------------------
ks_critical <- function(draws) 1.63 / sqrt(draws)
RNGkind("L'Ecuyer-CMRG")

## 10.1 interpolation error against direct computation, random R-squared values
set.seed(1)
interpolation_errors <- sapply(1:200, function(i) {
  n <- unique(design$n)[sample.int(length(unique(design$n)), 1)]
  k <- sample.int(number_of_candidates, 1)
  method <- sample(tabulated_methods, 1); r2 <- runif(1, 0.001, 0.9)
  abs(log_bf_tabulated(tables_per_n[[as.character(n)]], method, r2, k) -
      suppressWarnings(log_bayes_factor_model(method, r2, n, k)))
})
record_check("10.1 tabulated log BF equals direct computation (200 random cases)",
             max(interpolation_errors) < 1e-3, sprintf("largest error %.1e", max(interpolation_errors)))

## 10.2 the Beta(a, 1) g-mixture with a = 1/2 is the flat R closed form
v2_errors <- sapply(list(c(0.19, 60, 3), c(0.05, 100, 5), c(0.3, 200, 10), c(0.4, 50, 1)), function(s)
  abs(log_bayes_factor_beta_r_squared(s[1], s[2], s[3], 0.5) -
      log_bayes_factor_flat_r_vs_intercept_only(s[1], s[2], s[3])))
record_check("10.2 Beta(a, 1) g-mixture with a = 1/2 equals flat R", max(v2_errors) < 1e-3,
             sprintf("largest difference %.1e", max(v2_errors)))

## 10.3 the mixture with a = k/2 (g/n ~ Inverse-Gamma(1, 1/2)) implies R-squared ~ Beta(k/2, 1)
set.seed(2)
draws_v3 <- 2e5
v3_distances <- sapply(c(2, 5, 12), function(k) {
  f_squared <- rchisq(draws_v3, k) / (2 * rexp(draws_v3))
  suppressWarnings(unname(ks.test(f_squared / (1 + f_squared), "pbeta", k / 2, 1)$statistic))
})
record_check("10.3 g-mixture with a = k/2 implies R-squared ~ Beta(k/2, 1) (k = 2, 5, 12)",
             max(v3_distances) < ks_critical(draws_v3),
             sprintf("KS distances %s; critical %.4f", paste(sprintf("%.4f", v3_distances), collapse = " "),
                     ks_critical(draws_v3)))

## 10.4 uniformly drawn correlation matrices (LKJ, eta = 1) imply R-squared ~ Beta(k/2, 1)
set.seed(3)
draws_v4 <- 4000
v4_distances <- sapply(c(1, 2, 5), function(k) {
  r_squared <- replicate(draws_v4, {
    correlation <- draw_lkj_correlation_matrix(k + 1, eta = 1)   # variable 1 is the outcome
    rho <- correlation[1, -1]
    drop(rho %*% solve(correlation[-1, -1, drop = FALSE], rho))
  })
  suppressWarnings(unname(ks.test(r_squared, "pbeta", k / 2, 1)$statistic))
})
record_check("10.4 uniform correlation matrices imply R-squared ~ Beta(k/2, 1) (k = 1, 2, 5)",
             max(v4_distances) < ks_critical(draws_v4),
             sprintf("KS distances %s; critical %.4f", paste(sprintf("%.4f", v4_distances), collapse = " "),
                     ks_critical(draws_v4)))

## 10.5 with W ~ Inverse-Gamma(b, 1/2) the mixture implies R-squared ~ Beta(1/2, b) for every k
set.seed(4)
draws_v5 <- 2e5
v5_settings <- expand.grid(k = c(1, 5, 20), b = c(0.5, 2))
v5_distances <- mapply(function(k, b) {
  share <- if (k == 1) 1 else rbeta(draws_v5, 0.5, (k - 1) / 2)
  f_squared <- share * (0.5 / rgamma(draws_v5, b)) * rchisq(draws_v5, k)
  suppressWarnings(unname(ks.test(f_squared / (1 + f_squared), "pbeta", 0.5, b)$statistic))
}, v5_settings$k, v5_settings$b)
record_check("10.5 W ~ Inverse-Gamma(b, 1/2) implies R-squared ~ Beta(1/2, b) (b = 1/2, 2; k = 1, 5, 20)",
             max(v5_distances) < ks_critical(draws_v5),
             sprintf("KS distances %s; critical %.4f", paste(sprintf("%.4f", v5_distances), collapse = " "),
                     ks_critical(draws_v5)))


## ---- The simulation ------------------------------------------------------------
## R-squared of every subset, from the sample correlations
subset_r_squared <- function(predictors, outcome) {
  correlations <- cor(cbind(outcome, predictors))
  with_outcome <- correlations[1, -1]; among <- correlations[-1, -1]
  apply(all_models, 1, function(in_model) if (!any(in_model)) 0 else
    drop(with_outcome[in_model] %*% solve(among[in_model, in_model, drop = FALSE], with_outcome[in_model])))
}

## selection summaries from posterior model probabilities (log scale)
selection_from <- function(log_posterior) {
  posterior <- exp(log_posterior - max(log_posterior)); posterior <- posterior / sum(posterior)
  inclusion <- colSums(all_models * posterior)
  c(selection_from_inclusion(inclusion),
    best_is_true = as.numeric(all(all_models[which.max(posterior), ] == true_model)))
}
selection_from_inclusion <- function(inclusion) {
  median_model <- inclusion > 0.5
  c(median_is_true   = as.numeric(all(median_model == true_model)),
    missed           = sum(!median_model & true_model),
    false_inclusions = sum(median_model & !true_model),
    setNames(inclusion, paste0("inclusion_x", seq_along(inclusion))))
}

simulate_data <- function(condition, fixed_objects) generate_C(condition, condition$REPLICATION)

analyse_dataset <- function(condition, dat, fixed_objects) {
  n <- condition$n
  r_squared <- subset_r_squared(dat$predictors, dat$outcome)
  tables <- fixed_objects$tables_per_n[[as.character(n)]]
  log_bf <- c(lapply(setNames(nm = tabulated_methods), function(method)
                log_bf_tabulated(tables, method, r_squared, model_sizes)),
              list(BIC = -0.5 * (n * log(1 - pmin(r_squared, 1 - 1e-12)) + model_sizes * log(n))))
  ours <- unlist(lapply(names(log_bf), function(method)
    unlist(lapply(model_priors, function(model_prior) {
      values <- selection_from(log_bf[[method]] + log_model_prior[[model_prior]])
      setNames(values, paste(method, model_prior, names(values), sep = "|"))
    }))))
  ## hyper-g/n from BAS: the same 4,096 models, enumerated
  frame <- data.frame(y = dat$outcome, dat$predictors)
  bas_values <- unlist(lapply(model_priors, function(model_prior) {
    fit <- bas.lm(y ~ ., data = frame, prior = "hyper-g-n", alpha = 3, method = "BAS",
                  n.models = 2^number_of_candidates,
                  modelprior = if (model_prior == "uniform") uniform() else beta.binomial(1, 1))
    inclusion <- fit$probne0[-1]                                  # drop the intercept
    best <- logical(number_of_candidates)
    best[setdiff(fit$which[[which.max(fit$postprobs)]], 0)] <- TRUE
    values <- c(selection_from_inclusion(inclusion), best_is_true = as.numeric(all(best == true_model)))
    setNames(values, paste("hyper-g/n (BAS)", model_prior, names(values), sep = "|"))
  }))
  c(ours, bas_values)
}

results_file <- file.path(study_folder, "simulation.rds")
if (file.exists(results_file)) {
  simulation <- readRDS(results_file)
  message("Using cached results: ", results_file)
} else {
  simulation <- runSimulation(design, replications = n_sim, generate = simulate_data,
                              analyse = analyse_dataset,
                              fixed_objects = list(tables_per_n = tables_per_n),
                              parallel = number_of_cores > 1, ncores = number_of_cores,
                              packages = c("flatcor", "BayesFactor", "BAS"), store_results = TRUE,
                              control = list(include_replication_index = TRUE), progress = FALSE,
                              ## the resume file of an interrupted run, kept per study so studies
                              ## never pick up each other's state
                              save_details = list(out_rootdir = study_folder,
                                                  tmpfilename = "resume_simulation.qs2"))
  saveRDS(simulation, results_file)
}

## one row per dataset, method and prior over models
## (without a summarise function, runSimulation() returns the per-dataset rows
## themselves: one row per dataset, all methods side by side)
raw_results <- if (is(simulation, "SimDesign")) SimExtract(simulation, what = "results") else simulation
per_dataset <- raw_results |>
  group_by(n, rho) |> mutate(dataset = row_number()) |> ungroup() |>
  pivot_longer(contains("|"), names_to = "key", values_to = "value") |>
  separate_wider_delim(key, delim = "|", names = c("method", "model_prior", "measure")) |>
  pivot_wider(names_from = measure, values_from = value)


## ---- Summaries with their Monte Carlo SE -----------------------------------------
summaries <- per_dataset |>
  group_by(n, rho, model_prior, method) |>
  summarise(n_sim = n(),
            missed_mcse = sd(missed) / sqrt(n_sim), missed = mean(missed),
            false_inclusions_mcse = sd(false_inclusions) / sqrt(n_sim),
            false_inclusions = mean(false_inclusions),
            median_is_true = mean(median_is_true), best_is_true = mean(best_is_true),
            across(paste0("inclusion_x", 1:4), list(mean = mean, mcse = \(v) sd(v) / sqrt(length(v)))),
            .groups = "drop") |>
  mutate(median_is_true_mcse = mcse_of_proportion(median_is_true, n_sim),
         best_is_true_mcse = mcse_of_proportion(best_is_true, n_sim))

cat("\nMedian probability model: real predictors missed (of 4) and nulls included (of 8), with MCSE:\n")
print(as.data.frame(summaries |>
  transmute(n, rho, model_prior, method,
            missed = sprintf("%.2f (%.2f)", missed, missed_mcse),
            false_inclusions = sprintf("%.2f (%.2f)", false_inclusions, false_inclusions_mcse),
            median_is_true = sprintf("%.2f (%.2f)", median_is_true, median_is_true_mcse))), row.names = FALSE)


## ---- Ranking: AUC of inclusion probabilities, real against null -----------------
## pooled over the predictors of all datasets in a cell; the difference from
## flat R with a paired bootstrap over datasets (the same datasets for every method)
auc_real_vs_null <- function(inclusion_matrix) {   # rows: datasets, columns: 12 predictors
  real <- as.vector(inclusion_matrix[, true_model]); null <- as.vector(inclusion_matrix[, !true_model])
  ranks <- rank(c(real, null))
  (sum(ranks[seq_along(real)]) - length(real) * (length(real) + 1) / 2) / (length(real) * length(null))
}
bootstrap_replicates <- if (smoke_mode) 100 else 1000
inclusion_columns <- paste0("inclusion_x", seq_len(number_of_candidates))
RNGkind("L'Ecuyer-CMRG"); set.seed(10)   # the bootstrap's own, fixed stream
auc_table <- per_dataset |>
  group_by(n, rho, model_prior) |>
  group_modify(function(cell, keys) {
    matrices <- lapply(split(cell, cell$method), function(rows)
      as.matrix(rows[order(rows$dataset), inclusion_columns]))
    reference <- matrices[["flat R"]]
    resamples <- replicate(bootstrap_replicates, sample.int(nrow(reference), replace = TRUE),
                           simplify = FALSE)
    bind_rows(lapply(names(matrices), function(method) {
      m <- matrices[[method]]
      differences <- vapply(resamples, function(rows)
        auc_real_vs_null(m[rows, , drop = FALSE]) - auc_real_vs_null(reference[rows, , drop = FALSE]),
        numeric(1))
      tibble(method = method, auc = auc_real_vs_null(m),
             difference_from_flat_r = auc_real_vs_null(m) - auc_real_vs_null(reference),
             difference_mcse = sd(differences))
    }))
  }) |>
  ungroup()
cat("\nAUC of inclusion probabilities, real vs null predictors (difference from flat R, bootstrap MCSE):\n")
print(as.data.frame(auc_table |> mutate(across(c(auc, difference_from_flat_r, difference_mcse),
                                               \(v) round(v, 3)))), row.names = FALSE)


## ---- Figures ---------------------------------------------------------------------
long <- summaries |>
  select(n, rho, model_prior, method, median_is_true, false_inclusions, missed,
         median_is_true_mcse, false_inclusions_mcse, missed_mcse) |>
  pivot_longer(c(median_is_true, false_inclusions, missed), names_to = "measure", values_to = "value") |>
  mutate(mcse = case_when(measure == "median_is_true" ~ median_is_true_mcse,
                          measure == "false_inclusions" ~ false_inclusions_mcse,
                          TRUE ~ missed_mcse),
         measure = recode(measure,
           median_is_true   = "true model chosen (share)",
           false_inclusions = "null predictors included (of 8)",
           missed           = "real predictors missed (of 4)"))
save_figure(
  ggplot(long, aes(n, value, colour = method, linetype = factor(rho))) +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    geom_errorbar(aes(ymin = value - 2 * mcse, ymax = value + 2 * mcse), width = 0.03, alpha = 0.5) +
    scale_x_log10(breaks = c(50, 100, 200, 500)) + method_scales(methods, linetype = FALSE) +
    facet_grid(measure ~ model_prior, scales = "free_y") +
    theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = NULL, colour = NULL, linetype = "correlation",
         title = "Variable selection: 12 candidates, 4 real (0.3, 0.2, 0.15, 0.1)",
         subtitle = sprintf("median probability model; %d datasets per cell; bars: 2 MCSE; columns: prior over models",
                            n_sim)),
  "study_10_variable_selection", width = 10, height = 9)

inclusion_long <- summaries |>
  select(n, rho, model_prior, method, matches("^inclusion_x[1-4]_(mean|mcse)$")) |>
  pivot_longer(matches("^inclusion_x"), names_to = c("predictor", ".value"),
               names_pattern = "inclusion_(x[1-4])_(mean|mcse)") |>
  mutate(slope = factor(paste("true slope", true_slopes_C[as.integer(sub("x", "", predictor))]),
                        levels = paste("true slope", true_slopes_C[1:4])))
save_figure(
  ggplot(filter(inclusion_long, model_prior == "beta-binomial"),
         aes(n, mean, colour = method, linetype = factor(rho))) +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 100, 200, 500)) + method_scales(methods, linetype = FALSE) +
    facet_wrap(~ slope, nrow = 1) + theme(legend.position = "bottom") +
    labs(x = "sample size (log scale)", y = "mean inclusion probability", colour = NULL,
         linetype = "correlation",
         title = "Inclusion probability of each real predictor (beta-binomial prior over models)"),
  "study_10_inclusion", width = 12, height = 5)

save_study_results("10", list(summaries = summaries, auc = auc_table))
report_checks()
