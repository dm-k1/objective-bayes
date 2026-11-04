## =====================================================================
##  Study 8: how does the test behave for factors?
##  (docs/appendix/C_study_design.md, study 8; replaces study 11)
##
##  Aims. The block conditional test (a term of q columns, a factor with
##  q + 1 levels) against the joint test and JZS: how often each supports the
##  null, supports an effect or stays undecided, how well each separates a
##  null factor from a small one, and whether its evidence depends on what the
##  other predictors explain (it should not: the conditional Lindley paradox).
##  Data. Data set A extended with a factor (R/data/data_set_A_factors.R):
##  q = 2 to 5, n = 50, 200, 1,000, 2 or 10 other predictors, factor effect
##  Cohen's f^2 = 0, .02, .15. Each dataset is analysed three times, with the
##  others explaining R-squared 0, .5 and .9; the predictors, factor and noise
##  are the same, so only the others' slopes differ.
##  Estimands, per method, cell and value of what the others explain: the
##  share of each answer at BF 3; the AUC for a null factor against a small
##  one and its difference from JZS; the slope of the median log BF10 for a
##  null factor against what the others explain.
##  Methods. The block conditional test (3 priors), the joint test (3), JZS as
##  BayesFactor computes it for a factor (lmBF: a g for the factor, scale 1/2,
##  and one for the covariates, scale sqrt(2)/4) and Zellner-Siow (one g, scale
##  1, joint). All on the same datasets.
##  Performance. Shares with their Monte Carlo SE; AUC differences from JZS
##  and the slopes with paired bootstrap SEs.
##  n_sim. 625 datasets per cell (R/shared/precision.R, rate_per_cell):
##  every per-cell share has MCSE <= .02.
##
##  Checks:
##    8.1 (was F1) The conditional tests give the same log BF for the three
##        values of what the others explain, in every dataset (an exact
##        property: they see the factor only after the others are projected
##        out, so the threshold is numerical, 1e-8).
##    8.2 (was F2) The three conditional tests have the same AUC in every
##        cell (exact: all rank datasets by the block's partial R-squared).
##    8.3 (was F3) JZS (lmBF) is computed precisely enough: its reported
##        proportional error is below .05 in every dataset (recomputed with
##        10^5 and then 10^6 draws where it is above .02).
##    8.4 (new) For a null factor, the slope of the median log BF10 against
##        what the others explain is within 2 MCSE of zero for each
##        conditional test.
##
##  Run from the project folder:  Rscript R/studies/08_factors.R
##  Smoke test (one corner, 10 datasets per cell):
##    VALIDATION_SMOKE=1 Rscript R/studies/08_factors.R
##  Raw results are cached in cache/study_08/ (delete them to rerun); the
##  tables and checks are saved to results/study_08.rds.
## =====================================================================
library(SimDesign)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/shared/settings.R")             # flatcor, cores, record_check(), save_figure()
source("R/shared/comparison_methods.R")   # contenders, log_bayes_factor_model(), method_scales()
source("R/shared/precision.R")            # n_sim_of(), mcse_of_proportion()
source("R/data/data_sets.R")            # with_stream_of()
source("R/data/data_set_A_factors.R")   # design_A_factors(), generate_A_factors()

study_folder <- file.path(cache_folder, if (smoke_mode) "study_08_smoke" else "study_08")
dir.create(study_folder, recursive = TRUE, showWarnings = FALSE)

design <- design_A_factors()
if (smoke_mode) design <- subset(design, q == 2 & n == 50 & other_predictors == 2)
n_sim <- n_sim_of("rate_per_cell")

factor_methods <- c(paste(names(contenders), "conditional (block)"),
                    paste(names(contenders), "joint"), "JZS (lmBF)", "Zellner-Siow")
method_key <- setNames(make.names(factor_methods), factor_methods)


## ---- The simulation ------------------------------------------------------------
simulate_data <- function(condition, fixed_objects) generate_A_factors(condition, condition$REPLICATION)

## All methods for one dataset, for each value of what the others explain.
## Returns log BF10 per method and value (names "<method>.R2o<value>"), the
## JZS proportional error ("jzs_error.R2o<value>") and the replication number
## (SimDesign does not keep it with the stored results, and reserves names in
## capitals, hence "replication").
analyse_dataset <- function(condition, dat, fixed_objects) {
  q <- condition$q; n <- condition$n; p <- condition$other_predictors; k <- q + p
  others <- dat$others; group <- dat$group
  design_full    <- model.matrix(~ group + others)
  design_without <- cbind(1, others)
  formula_full    <- as.formula(paste("y ~ group +", paste(colnames(others), collapse = " + ")))
  formula_without <- as.formula(paste("y ~", paste(colnames(others), collapse = " + ")))
  dataset_id <- paste(cell_id_A_factors(condition), condition$REPLICATION, sep = "-r")
  c(replication = condition$REPLICATION, unlist(lapply(explained_by_others_values, function(r_squared_others) {
    outcome <- outcome_A_factors(dat, condition, r_squared_others)
    total_ss         <- sum((outcome - mean(outcome))^2)
    residual_full    <- sum(qr.resid(qr(design_full), outcome)^2)
    residual_without <- sum(qr.resid(qr(design_without), outcome)^2)
    r_squared_full    <- 1 - residual_full / total_ss
    r_squared_without <- 1 - residual_without / total_ss
    block_partial_r   <- sqrt(1 - residual_full / residual_without)
    joint <- function(model_method)
      log_bayes_factor_model(model_method, r_squared_full, n, k) -
      log_bayes_factor_model(model_method, r_squared_without, n, k - q)
    data_frame <- data.frame(y = outcome, group = group, others)
    ## lmBF integrates over two g's by Monte Carlo: give it the dataset's own
    ## stream, so its answer does not depend on the order or number of cores
    jzs <- with_stream_of(paste(dataset_id, "jzs", r_squared_others, sep = "-"), {
      jzs_full <- BayesFactor::lmBF(formula_full, data_frame, progress = FALSE)
      for (iterations in c(1e5, 1e6))
        if (BayesFactor::extractBF(jzs_full)$error > 0.02)
          jzs_full <- BayesFactor::recompute(jzs_full, iterations = iterations, progress = FALSE)
      jzs_without <- BayesFactor::lmBF(formula_without, data_frame, progress = FALSE)
      c(log_bf = BayesFactor::extractBF(jzs_full, logbf = TRUE)$bf -
          BayesFactor::extractBF(jzs_without, logbf = TRUE)$bf,
        error = BayesFactor::extractBF(jzs_full)$error)
    })
    values <- c(
      vapply(contenders, function(prior)
        log_bayes_factor_conditional(block_partial_r, n, k, prior, number_of_columns = q), numeric(1)),
      vapply(names(contenders), joint, numeric(1)),
      jzs[["log_bf"]],
      joint("Zellner-Siow"),
      jzs[["error"]])
    setNames(values, paste0(c(unname(method_key), "jzs_error"), ".R2o", r_squared_others))
  })))
}

## the share of each answer at BF 3, per method and value of what the others explain
summarise_cell <- function(condition, results, fixed_objects) {
  log_bfs <- results[, !startsWith(names(results), "jzs_error") & names(results) != "replication",
                     drop = FALSE]
  c(setNames(colMeans(log_bfs < -log(3)), paste0("null.", names(log_bfs))),
    setNames(colMeans(log_bfs > log(3)), paste0("effect.", names(log_bfs))))
}

results_file <- file.path(study_folder, "simulation.rds")
if (file.exists(results_file)) {
  simulation <- readRDS(results_file)
  message("Using cached results: ", results_file)
} else {
  simulation <- runSimulation(design, replications = n_sim, generate = simulate_data,
                              analyse = analyse_dataset, summarise = summarise_cell,
                              parallel = number_of_cores > 1, ncores = number_of_cores,
                              packages = c("flatcor", "BayesFactor"), store_results = TRUE,
                              control = list(include_replication_index = TRUE), progress = FALSE,
                              save_details = list(out_rootdir = study_folder,
                                                  tmpfilename = "resume_simulation.qs2"))
  saveRDS(simulation, results_file)
}

## every dataset's log BF10, one row per dataset, method and value of what the
## others explain
raw <- SimExtract(simulation, what = "results")
jzs_errors <- raw |> select(starts_with("jzs_error")) |> unlist(use.names = FALSE)
log_bf <- raw |>
  select(-starts_with("jzs_error"), -any_of("REPLICATION")) |>   # SimDesign's own copy, when present
  rename(REPLICATION = replication) |>
  pivot_longer(contains(".R2o"), names_to = c("key", "explained_by_others"), names_sep = "\\.R2o",
               values_to = "log_BF10") |>
  mutate(explained_by_others = as.numeric(explained_by_others),
         method = factor(names(method_key)[match(key, method_key)], levels = factor_methods),
         answer = case_when(log_BF10 < -log(3) ~ "says null", log_BF10 > log(3) ~ "says effect",
                            TRUE ~ "undecided"))


## ---- Shares of each answer, with their Monte Carlo SE ---------------------------
shares_per_cell <- log_bf |>
  count(q, n, other_predictors, factor_effect, explained_by_others, method, answer, name = "datasets") |>
  complete(nesting(q, n, other_predictors, factor_effect, explained_by_others, method), answer,
           fill = list(datasets = 0)) |>
  group_by(q, n, other_predictors, factor_effect, explained_by_others, method) |>
  mutate(n_sim = sum(datasets), share = datasets / n_sim, mcse = mcse_of_proportion(share, n_sim)) |>
  ungroup()

## averaged over q (independent cells), by n, the number of others and their R-squared
shares_over_q <- shares_per_cell |>
  group_by(n, other_predictors, factor_effect, explained_by_others, method, answer) |>
  summarise(cells = n(), share = mean(share), mcse = sqrt(sum(mcse^2)) / cells, .groups = "drop")

cat("\nShare of each answer at BF 3, averaged over q (with MCSE), null and small factors:\n")
print(as.data.frame(shares_over_q |> filter(factor_effect %in% c("null", "small")) |>
  mutate(share = sprintf("%.3f (%.3f)", share, mcse)) |> select(-mcse, -cells) |>
  pivot_wider(names_from = answer, values_from = share)), row.names = FALSE)


## ---- Discrimination: AUC null vs small, and its difference from JZS -------------
## P(a value from higher exceeds one from lower), ties counting half (Mann-Whitney)
area_under_curve <- function(higher, lower) {
  ranks <- rank(c(higher, lower))
  (sum(ranks[seq_along(higher)]) - length(higher) * (length(higher) + 1) / 2) /
    (length(higher) * length(lower))
}
bootstrap_replicates <- if (smoke_mode) 100 else 1000

RNGkind("L'Ecuyer-CMRG"); set.seed(8)   # the bootstrap's own, fixed stream
auc <- log_bf |>
  filter(factor_effect %in% c("null", "small")) |>
  group_by(q, n, other_predictors, explained_by_others) |>
  group_modify(function(cell, keys) {
    evidence_for_null <- function(method_name, effect)          # log BF01, in replication order
      -cell$log_BF10[cell$method == method_name & cell$factor_effect == effect][
        order(cell$REPLICATION[cell$method == method_name & cell$factor_effect == effect])]
    replications <- sum(cell$method == "JZS (lmBF)" & cell$factor_effect == "null")
    resamples <- replicate(bootstrap_replicates,
                           list(sample(replications, replace = TRUE), sample(replications, replace = TRUE)),
                           simplify = FALSE)
    auc_of <- function(method_name, rows_null = seq_len(replications), rows_small = seq_len(replications))
      area_under_curve(evidence_for_null(method_name, "null")[rows_null],
                       evidence_for_null(method_name, "small")[rows_small])
    bind_rows(lapply(factor_methods, function(method_name) {
      versus_jzs <- vapply(resamples, function(rows)
        auc_of(method_name, rows[[1]], rows[[2]]) - auc_of("JZS (lmBF)", rows[[1]], rows[[2]]), numeric(1))
      tibble(method = method_name, auc = auc_of(method_name),
             auc_minus_jzs = auc_of(method_name) - auc_of("JZS (lmBF)"),
             mcse_minus_jzs = sd(versus_jzs))
    }))
  }) |>
  ungroup() |>
  mutate(method = factor(method, levels = factor_methods))

cat("\nAUC (null vs small factor) minus JZS, mean over q, by n and what the others explain:\n")
print(as.data.frame(auc |> group_by(method, n, explained_by_others) |>
  summarise(auc_minus_jzs = mean(auc_minus_jzs), mcse = sqrt(sum(mcse_minus_jzs^2)) / n(),
            .groups = "drop")), digits = 3, row.names = FALSE)


## ---- The paradox in one number: slope of the median log BF10 for a null factor --
## against what the others explain, with a bootstrap SE over datasets
RNGkind("L'Ecuyer-CMRG"); set.seed(84)
slope_of_medians <- function(values, explained) {
  medians <- tapply(values, explained, median)
  unname(coef(lm(medians ~ as.numeric(names(medians))))[2])
}
paradox <- log_bf |>
  filter(factor_effect == "null") |>
  group_by(q, n, other_predictors, method) |>
  group_modify(function(cell, keys) {
    replications <- sort(unique(cell$REPLICATION))
    resampled <- vapply(seq_len(bootstrap_replicates), function(b) {
      chosen <- sample(replications, replace = TRUE)
      rows <- unlist(lapply(chosen, function(r) which(cell$REPLICATION == r)))
      slope_of_medians(cell$log_BF10[rows], cell$explained_by_others[rows])
    }, numeric(1))
    tibble(slope = slope_of_medians(cell$log_BF10, cell$explained_by_others), mcse = sd(resampled))
  }) |>
  ungroup()

cat("\nSlope of the median log BF10 for a null factor against the others' R-squared,\n",
    "mean over cells (with MCSE):\n", sep = "")
print(as.data.frame(paradox |> group_by(method, n) |>
  summarise(slope = mean(slope), mcse = sqrt(sum(mcse^2)) / n(), .groups = "drop")),
  digits = 3, row.names = FALSE)


## ---- Checks --------------------------------------------------------------------
conditional_spread <- log_bf |>
  filter(grepl("conditional", method)) |>
  group_by(q, n, other_predictors, factor_effect, REPLICATION, method) |>
  summarise(spread = diff(range(log_BF10)), .groups = "drop")
record_check("8.1 conditional tests: the same for every R-squared of the others",
             max(conditional_spread$spread) < 1e-8,
             sprintf("largest difference in log BF10 %.1e (%d datasets x three contenders)",
                     max(conditional_spread$spread), nrow(conditional_spread) / length(contenders)))

conditional_auc_spread <- auc |>
  filter(grepl("conditional", method)) |>
  group_by(q, n, other_predictors, explained_by_others) |>
  summarise(spread = diff(range(auc)), .groups = "drop")
record_check("8.2 the three conditional tests: the same AUC in every cell",
             max(conditional_auc_spread$spread) < 1e-12,
             sprintf("largest difference %.1e over %d cells", max(conditional_auc_spread$spread),
                     nrow(conditional_auc_spread)))

record_check("8.3 JZS (lmBF): proportional error < .05",
             max(jzs_errors, na.rm = TRUE) < 0.05,
             sprintf("median %.4f, 99%% %.4f, largest %.3f", median(jzs_errors, na.rm = TRUE),
                     quantile(jzs_errors, 0.99, na.rm = TRUE), max(jzs_errors, na.rm = TRUE)))

conditional_paradox <- filter(paradox, grepl("conditional", method))
record_check("8.4 conditional tests: no slope against the others' R-squared (null factor)",
             all(abs(conditional_paradox$slope) <= 2 * conditional_paradox$mcse + 1e-8),
             sprintf("largest |slope| %.1e; slopes of the joint tests and JZS: %s",
                     max(abs(conditional_paradox$slope)),
                     paste(sprintf("%s %.2f", levels(paradox$method)[4:8],
                                   tapply(paradox$slope, paradox$method, mean)[4:8]), collapse = "; ")))


## ---- Figures -------------------------------------------------------------------
by_others <- facet_grid(other_predictors ~ explained_by_others,
  labeller = labeller(other_predictors    = function(v) paste(v, "other predictors"),
                      explained_by_others = function(v) paste("others explain R-squared =", v)))

save_figure(
  ggplot(auc |> filter(method != "JZS (lmBF)") |>
           group_by(method, n, other_predictors, explained_by_others) |>
           summarise(auc_minus_jzs = mean(auc_minus_jzs), .groups = "drop"),
         aes(n, auc_minus_jzs, colour = method, linetype = method)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200, 1000)) + method_scales(setdiff(factor_methods, "JZS (lmBF)")) +
    by_others +
    labs(x = "sample size (log scale)", y = "AUC (null vs small factor) minus JZS", colour = NULL,
         title = "Testing a factor: separating nulls from small effects, against JZS (mean over q = 2-5)",
         subtitle = "above 0: separates better than JZS"),
  "study_08_auc", width = 12, height = 7)

save_figure(
  ggplot(shares_over_q |>
           filter((factor_effect == "null" & answer == "says null") |
                    (factor_effect == "small" & answer == "says null")) |>
           mutate(rate = ifelse(factor_effect == "null", "true null confirmed (BF01 > 3)",
                                "small effect called null (BF01 > 3)")),
         aes(n, share, colour = method, linetype = method)) +
    geom_line(linewidth = 0.8, alpha = 0.8) + geom_point(size = 1.3, alpha = 0.8) +
    scale_x_log10(breaks = c(50, 200, 1000)) + method_scales(factor_methods) +
    facet_grid(rate + other_predictors ~ explained_by_others,
      labeller = labeller(other_predictors    = function(v) paste(v, "other predictors"),
                          explained_by_others = function(v) paste("others explain R-squared =", v))) +
    labs(x = "sample size (log scale)", y = "share of datasets", colour = NULL,
         title = "Testing a factor: calibration at 3 (mean over q = 2-5)",
         subtitle = "The conditional tests do not change along a row; the joint tests and JZS do (the conditional Lindley paradox)"),
  "study_08_calibration", width = 12, height = 11)

save_study_results("08", list(shares_per_cell = shares_per_cell, shares_over_q = shares_over_q, auc = auc,
                              paradox = paradox))
report_checks()
