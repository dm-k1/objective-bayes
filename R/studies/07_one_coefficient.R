## =====================================================================
##  Study 7: how does the test of one coefficient behave?
##  (docs/ADEMP.md, study 7; replaces studies 03, 04 and 06)
##
##  Aims. How often each test supports the null (BF01 > 3), supports an
##  effect (BF10 > 3) or stays undecided about one coefficient (x1); how that
##  depends on the number of predictors k; and whether differences between
##  priors are differences in accuracy (how datasets are ranked) or only in
##  where the test decides.
##  Data. Set A (R/data/data_sets.R): 72 designs (k, n, predictor
##  correlation, the other slopes) at slopes of x1 of 0, 0.1, 0.2 and 0.3.
##  Estimands, per method and design cell: the share of each answer; the AUC
##  for separating slope 0 from slope 0.1 and its difference from JZS.
##  Methods. The conditional and the joint test under each of the three
##  priors, JZS (scale 0.354) and Zellner-Siow (scale 1), all on the same
##  datasets.
##  Performance. Shares with their Monte Carlo SE sqrt(p (1 - p) / n_sim);
##  AUC differences from JZS with a paired bootstrap SE.
##  n_sim. 625 datasets per cell (R/shared/precision.R): every per-cell share
##  has MCSE <= .02, and shares pooled over the 18 cells at each k <= .005.
##
##  Checks, set relative to the Monte Carlo error:
##    7.1 At every k, the three conditional tests separate nulls from small
##        effects equally well: their AUCs agree within 2 MCSE in every cell.
##    7.2 Flat partial r confirms true nulls at the same rate at k = 2 and
##        k = 20: the difference lies within the equivalence margin of +-.05,
##        judged by its 95% Monte Carlo interval.
##    7.3 Flat R confirms fewer true nulls at k = 20 than at k = 2: the drop
##        exceeds .05, judged by its 95% Monte Carlo interval.
##
##  Run from the project folder:  Rscript R/studies/07_one_coefficient.R
##  Smoke test (a corner of the design, 10 datasets per cell, about a minute):
##    VALIDATION_SMOKE=1 Rscript R/studies/07_one_coefficient.R
##  Raw results are cached in cache/study_07/ (delete them to rerun); the
##  tables and checks are saved to results/study_07.rds.
## =====================================================================
library(SimDesign)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/shared/settings.R")             # flatcor, cores, record_check(), save_figure()
source("R/shared/comparison_methods.R")   # the methods, log_bayes_factors_one_coefficient()
source("R/shared/precision.R")            # n_sim_of(), mcse_of_proportion()
source("R/data/data_sets.R")            # design_A(), generate_A()

study_folder <- file.path(cache_folder, if (smoke_mode) "study_07_smoke" else "study_07")
dir.create(study_folder, recursive = TRUE, showWarnings = FALSE)

design <- design_A()
if (smoke_mode) design <- subset(design, k %in% c(2, 20) & n == 50 & rho == 0 & others == "sparse")
n_sim <- n_sim_of("rate_per_cell")

## method names without spaces, for SimDesign's result columns
method_key <- setNames(make.names(one_coefficient_methods), one_coefficient_methods)


## ---- The simulation ------------------------------------------------------------
simulate_data <- function(condition, fixed_objects) generate_A(condition, condition$REPLICATION)

analyse_dataset <- function(condition, dat, fixed_objects) {
  full_fit <- lm(dat$outcome ~ dat$predictors)
  log_bfs  <- log_bayes_factors_one_coefficient(full_fit, 1, dat$outcome, dat$predictors)
  setNames(log_bfs$log_BF10, method_key[log_bfs$method])   # log BF10 for x1, one per method
}

## the share of each answer at BF 3, per method
summarise_cell <- function(condition, results, fixed_objects) {
  c(setNames(colMeans(results < -log(3)), paste0("null.", names(results))),
    setNames(colMeans(results > log(3)), paste0("effect.", names(results))))
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
                              ## its own resume file, so concurrent studies never collide
                              save_details = list(out_rootdir = study_folder,
                                                  tmpfilename = "resume_simulation.qs2"))
  saveRDS(simulation, results_file)
}

## every dataset's log BF10, one row per dataset and method
log_bf <- SimExtract(simulation, what = "results") |>
  pivot_longer(all_of(unname(method_key)), names_to = "key", values_to = "log_BF10") |>
  mutate(method = names(method_key)[match(key, method_key)],
         answer = case_when(log_BF10 < -log(3) ~ "says null", log_BF10 > log(3) ~ "says effect",
                            TRUE ~ "undecided"))


## ---- Shares of each answer, with their Monte Carlo SE ---------------------------
shares_per_cell <- log_bf |>
  count(k, n, rho, others, slope_x1, method, answer, name = "datasets") |>
  complete(nesting(k, n, rho, others, slope_x1, method), answer, fill = list(datasets = 0)) |>
  group_by(k, n, rho, others, slope_x1, method) |>
  mutate(n_sim = sum(datasets), share = datasets / n_sim, mcse = mcse_of_proportion(share, n_sim)) |>
  ungroup()

## pooled over the cells at each k: the mean of the per-cell shares, whose
## MCSE follows from the per-cell MCSEs (independent cells)
shares_per_k <- shares_per_cell |>
  group_by(k, slope_x1, method, answer) |>
  summarise(cells = n(), share = mean(share), mcse = sqrt(sum(mcse^2)) / cells, .groups = "drop")

cat("\nShare of each answer at BF 3, pooled over the cells at each k (with MCSE):\n")
print(as.data.frame(shares_per_k |> filter(grepl("conditional|JZS", method)) |>
  mutate(share = sprintf("%.3f (%.3f)", share, mcse)) |> select(-mcse, -cells) |>
  pivot_wider(names_from = answer, values_from = share)), row.names = FALSE)


## ---- Discrimination: AUC null vs slope 0.1, and its difference from JZS ---------
## P(a value from higher exceeds one from lower), ties counting half: the
## Mann-Whitney form, from ranks, so the bootstrap stays fast at 625 per group
area_under_curve <- function(higher, lower) {
  ranks <- rank(c(higher, lower))
  (sum(ranks[seq_along(higher)]) - length(higher) * (length(higher) + 1) / 2) /
    (length(higher) * length(lower))
}
bootstrap_replicates <- if (smoke_mode) 100 else 1000

RNGkind("L'Ecuyer-CMRG"); set.seed(7)   # the bootstrap's own, fixed stream
auc <- log_bf |>
  filter(slope_x1 %in% c(0, 0.1)) |>
  group_by(k, n, rho, others) |>
  group_modify(function(cell, keys) {
    evidence_for_null <- function(method_name, slope)        # log BF01, in replication order
      -cell$log_BF10[cell$method == method_name & cell$slope_x1 == slope][
        order(cell$REPLICATION[cell$method == method_name & cell$slope_x1 == slope])]
    replications <- sum(cell$method == "JZS default" & cell$slope_x1 == 0)
    resamples <- replicate(bootstrap_replicates,
                           list(sample(replications, replace = TRUE), sample(replications, replace = TRUE)),
                           simplify = FALSE)
    auc_of <- function(method_name, rows_null = seq_len(replications), rows_effect = seq_len(replications))
      area_under_curve(evidence_for_null(method_name, 0)[rows_null],
                       evidence_for_null(method_name, 0.1)[rows_effect])
    bind_rows(lapply(one_coefficient_methods, function(method_name) {
      versus_jzs <- vapply(resamples, function(rows)
        auc_of(method_name, rows[[1]], rows[[2]]) - auc_of("JZS default", rows[[1]], rows[[2]]), numeric(1))
      tibble(method = method_name, auc = auc_of(method_name),
             auc_minus_jzs = auc_of(method_name) - auc_of("JZS default"),
             mcse_minus_jzs = sd(versus_jzs),
             mcse_auc = sd(vapply(resamples, function(rows) auc_of(method_name, rows[[1]], rows[[2]]),
                                  numeric(1))))
    }))
  }) |>
  ungroup()

cat("\nAUC (null vs slope 0.1) minus JZS, mean over cells at each k:\n")
print(as.data.frame(auc |> group_by(k, method) |>
  summarise(auc = mean(auc), auc_minus_jzs = mean(auc_minus_jzs), .groups = "drop")), digits = 3)


## ---- Checks --------------------------------------------------------------------
conditional_methods <- paste(names(contenders), "conditional")

## 7.1 the three conditional tests rank datasets alike: their AUCs agree within
## 2 MCSE in every cell
auc_spread <- auc |>
  filter(method %in% conditional_methods) |>
  group_by(k, n, rho, others) |>
  summarise(spread = max(auc) - min(auc), mcse = max(mcse_auc), .groups = "drop")
record_check("7.1 the three conditional tests separate nulls from small effects equally well",
             all(auc_spread$spread <= 2 * auc_spread$mcse),
             sprintf("largest AUC spread %.4f; largest 2 MCSE %.4f (%d cells)",
                     max(auc_spread$spread), 2 * max(auc_spread$mcse), nrow(auc_spread)))

null_confirmed <- function(method_name, k_value)
  filter(shares_per_k, method == method_name, k == k_value, slope_x1 == 0, answer == "says null")
difference_interval <- function(method_name) {
  at_2 <- null_confirmed(method_name, 2); at_20 <- null_confirmed(method_name, 20)
  difference <- at_20$share - at_2$share
  mcse <- sqrt(at_2$mcse^2 + at_20$mcse^2)
  c(difference = difference, lower = difference - 1.96 * mcse, upper = difference + 1.96 * mcse)
}

## 7.2 flat partial r: the same rate at k = 2 and k = 20 (equivalence margin .05)
stable <- difference_interval("flat partial r conditional")
record_check("7.2 flat partial r confirms true nulls at the same rate at k = 2 and k = 20",
             stable[["lower"]] > -0.05 && stable[["upper"]] < 0.05,
             sprintf("difference %.3f, 95%% Monte Carlo interval [%.3f, %.3f]; margin +-.05",
                     stable[["difference"]], stable[["lower"]], stable[["upper"]]))

## 7.3 flat R: fewer at k = 20 than at k = 2, by more than .05
drift <- difference_interval("flat R conditional")
record_check("7.3 flat R confirms fewer true nulls at k = 20 than at k = 2",
             drift[["upper"]] < -0.05,
             sprintf("difference %.3f, 95%% Monte Carlo interval [%.3f, %.3f]",
                     drift[["difference"]], drift[["lower"]], drift[["upper"]]))


## ---- Figure: what each test says, by k and true slope ---------------------------
outcome_shares <- shares_per_k |>
  filter(grepl("conditional|JZS", method), slope_x1 %in% c(0, 0.1, 0.2)) |>
  mutate(test = factor(sub(" conditional", "", method),
                       levels = c("flat partial r", "flat R", "flat R-squared", "JZS default")),
         answer = factor(answer, levels = c("says effect", "undecided", "says null")),
         truth = factor(paste("true slope", slope_x1), levels = paste("true slope", c(0, 0.1, 0.2))))
save_figure(
  ggplot(outcome_shares, aes(test, share, fill = answer)) +
    geom_col(width = 0.75) +
    geom_text(data = filter(outcome_shares, share >= 0.08), aes(label = sprintf("%.2f", share)),
              position = position_stack(vjust = 0.5), size = 2.7, colour = "white") +
    scale_fill_manual(values = c(`says null` = "#542788", undecided = "grey70", `says effect` = "#B35806"),
                      breaks = c("says null", "undecided", "says effect"), name = "at BF 3:") +
    facet_grid(truth ~ k, labeller = labeller(k = function(v) paste("k =", v))) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") +
    labs(x = NULL, y = "share of datasets",
         title = "What each conditional test (and JZS) says about x1, as k grows",
         subtitle = sprintf("%d datasets per cell, pooled over n, correlation and pattern; correct: says null in the top row",
                            n_sim)),
  "study_07_outcomes", width = 11, height = 8)

save_study_results("07", list(shares_per_cell = shares_per_cell, shares_per_k = shares_per_k, auc = auc))
report_checks()
