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
##  Bayes factors per model (each depends only on R-squared, n and k): the
##  whole-model set of validation/comparison_methods.R,
##    flat-R                 R-squared ~ Beta(1/2, 1) for every k
##    uniform correlations   R-squared ~ Beta(k/2, 1): every valid correlation
##                           matrix of (y, x1 ... xk) equally likely (LKJ, eta = 1;
##                           check V4); rstanarm stan_lm's prior with eta = 1
##    JZS default (scale 0.354) and Zellner-Siow (scale 1)
##  and BIC (exp(-BIC/2) relative to the null) as a non-Bayesian reference.
##  (An earlier experiment, R-squared ~ Beta(k/(k+1), 1), barely differed from
##  flat-R and was dropped; validation/README.md, 09.)
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
##  About 3-5 minutes; cached.
## =====================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
source("validation/settings.R")   # parallel_map(), caching, record_check(), save_figure()
pkgload::load_all(".", quiet = TRUE)
source("validation/comparison_methods.R")   # the common method set

number_of_candidates <- 12
true_slopes <- c(0.3, 0.2, 0.15, 0.1, rep(0, 8))
true_model  <- true_slopes != 0
design <- expand.grid(sample_size = c(50, 100, 200, 500), predictor_correlation = c(0, 0.5))
replications_per_cell <- 100
methods <- c(model_methods, "BIC")


## all 4,096 subsets as a logical matrix, one row per model
all_models <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), number_of_candidates)))
model_sizes <- rowSums(all_models)


## ---- Log Bayes factor against the null, tabulated in R-squared ---------------
log_bf_direct <- function(method, r_squared, n, k) log_bayes_factor_model(method, r_squared, n, k)
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
  method <- sample(setdiff(methods, "BIC"), 1); r2 <- runif(1, 0.001, 0.9)
  abs(log_bf_table(tables_per_n[[as.character(n)]], method, r2, n, k) - log_bf_direct(method, r2, n, k))
})
record_check("V1. tabulated log BF vs direct (200 random cases)", max(interpolation_errors) < 1e-3,
             sprintf("largest error %.1e", max(interpolation_errors)))

## V2. the Beta(a, 1) Bayes factor (comparison_methods.R) with a = 1/2 is the
##     flat-R closed form: the same code then gives uniform correlations at a = k/2
v2_errors <- sapply(list(c(0.19, 60, 3), c(0.05, 100, 5), c(0.3, 200, 10), c(0.4, 50, 1)), function(s)
  abs(log_bayes_factor_beta_r_squared(s[1], s[2], s[3], 0.5) -
      log_bayes_factor_flat_r_vs_intercept_only(s[1], s[2], s[3])))
record_check("V2. Beta(a, 1) version with a = 1/2 equals flat-R", max(v2_errors) < 1e-3,
             sprintf("largest difference %.1e", max(v2_errors)))

## V3. the mixture with a = k/2 (B = 1, g/n ~ Inverse-Gamma(1, 1/2)) implies
##     R-squared ~ Beta(k/2, 1)
set.seed(2)
v3_distances <- sapply(c(2, 5, 12), function(k) {
  f_squared <- rchisq(2e5, k) / (2 * rexp(2e5))
  suppressWarnings(unname(ks.test(f_squared / (1 + f_squared), "pbeta", k / 2, 1)$statistic))
})
record_check("V3. g-mixture with a = k/2 implies R-squared ~ Beta(k/2, 1) (k = 2, 5, 12)", max(v3_distances) < 0.01,
             sprintf("KS distances %s", paste(sprintf("%.4f", v3_distances), collapse = " ")))

## V4. "uniform correlations": draw whole correlation matrices of (y, x1 ... xk) uniformly
## (LKJ with eta = 1, C-vine sampler of Lewandowski, Kurowicka & Joe, 2009) and compute
## the population R-squared of y on the x's; it should be Beta(k/2, 1)
## (draw_lkj_correlation_matrix() is in comparison_methods.R)
set.seed(3)
v4_distances <- sapply(c(1, 2, 5), function(k) {
  r_squared <- replicate(4000, {
    correlation <- draw_lkj_correlation_matrix(k + 1, eta = 1)   # variable 1 is the outcome
    rho <- correlation[1, -1]
    drop(rho %*% solve(correlation[-1, -1, drop = FALSE], rho))
  })
  suppressWarnings(unname(ks.test(r_squared, "pbeta", k / 2, 1)$statistic))
})
record_check("V4. uniform correlation matrices imply R-squared ~ Beta(k/2, 1) (k = 1, 2, 5)", max(v4_distances) < 0.03,
             sprintf("KS distances %s", paste(sprintf("%.4f", v4_distances), collapse = " ")))

## V5. constant in k does not single out flat-R: with the same B but
##     W ~ Inverse-Gamma(b, 1/2), the mixture implies R-squared ~ Beta(1/2, b) for
##     every k (docs/2_how_it_calculates.md, section 8). b = 1 is flat-R; only
##     b = 1 is Jeffreys' flat correlation at k = 1.
set.seed(4)
v5_settings <- expand.grid(k = c(1, 5, 20), b = c(0.5, 2))
v5_distances <- mapply(function(k, b) {
  share <- if (k == 1) 1 else rbeta(2e5, 0.5, (k - 1) / 2)
  f_squared <- share * (0.5 / rgamma(2e5, b)) * rchisq(2e5, k)
  suppressWarnings(unname(ks.test(f_squared / (1 + f_squared), "pbeta", 0.5, b)$statistic))
}, v5_settings$k, v5_settings$b)
record_check("V5. W ~ Inverse-Gamma(b, 1/2) implies R-squared ~ Beta(1/2, b) (b = 1/2, 2; k = 1, 5, 20)",
             max(v5_distances) < 0.01,
             sprintf("KS distances %s", paste(sprintf("%.4f", v5_distances), collapse = " ")))


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
                 inclusion_x4          = inclusion[4],   # the smallest real effect, 0.1
                 inclusion             = I(list(unname(inclusion))))   # all 12, for the curves below
    }))
  }))
}

runs <- run_simulation_grid(
  design,
  results_file = file.path(results_folder, "sim_variable_selection.rds"),
  settings = list(replications_per_cell = replications_per_cell, true_slopes = true_slopes,
                  methods = methods, version = 5),
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

method_colours <- method_colours[methods]
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

## ---- Equal or worse? Whole trade-off curves, not one point ----------------------
## The 0.5 cutoff gives each method one point (misses vs false inclusions). Sliding
## the cutoff traces a curve. If two methods' curves coincide, the prior only chose
## where on the curve you sit; if one lies below, that method separates real from
## null predictors better. Summary: AUC of "real vs null" over all predictors in a
## cell (the cutoff is shared across datasets, as in practice), with a bootstrap
## over datasets for the difference from flat-R (the same datasets for every method).
auc_real_vs_null <- function(inclusion_matrix) {   # rows: datasets, columns: 12 predictors
  real <- as.vector(inclusion_matrix[, true_model]); null <- as.vector(inclusion_matrix[, !true_model])
  ranks <- rank(c(real, null))
  (sum(ranks[seq_along(real)]) - length(real) * (length(real) + 1) / 2) / (length(real) * length(null))
}
inclusion_matrix_of <- function(rows) do.call(rbind, rows$inclusion[order(rows$replication)])

set.seed(4)
curve_cells <- split(runs, list(runs$predictor_correlation, runs$sample_size, runs$model_prior), drop = TRUE)
auc_table <- do.call(rbind, lapply(curve_cells, function(cell) {
  per_method <- lapply(split(cell, cell$method), inclusion_matrix_of)
  reference <- per_method[["flat-R"]]
  boot_index <- replicate(1000, sample.int(nrow(reference), replace = TRUE))
  do.call(rbind, lapply(names(per_method), function(method) {
    m <- per_method[[method]]
    differences <- apply(boot_index, 2, function(i) auc_real_vs_null(m[i, ]) - auc_real_vs_null(reference[i, ]))
    data.frame(predictor_correlation = cell$predictor_correlation[1], sample_size = cell$sample_size[1],
               model_prior = cell$model_prior[1], method = method, auc = auc_real_vs_null(m),
               difference_from_flat_r = auc_real_vs_null(m) - auc_real_vs_null(reference),
               lower = quantile(differences, 0.025), upper = quantile(differences, 0.975))
  }))
}))
rownames(auc_table) <- NULL
cat("\nAUC of inclusion probabilities, real vs null predictors (difference from flat-R, 95% bootstrap):\n")
print(auc_table |> mutate(across(c(auc, difference_from_flat_r, lower, upper), \(v) round(v, 3))) |>
        arrange(model_prior, predictor_correlation, sample_size, method), row.names = FALSE)
saveRDS(auc_table, file.path(results_folder, "variable_selection_auc.rds"))

## the curves themselves: mean misses (of 4) against mean false inclusions (of 8) per cutoff
cutoffs <- c(seq(0.01, 0.99, by = 0.01))
curves <- do.call(rbind, lapply(curve_cells, function(cell)
  do.call(rbind, lapply(split(cell, cell$method), function(rows) {
    m <- inclusion_matrix_of(rows)
    data.frame(predictor_correlation = rows$predictor_correlation[1], sample_size = rows$sample_size[1],
               model_prior = rows$model_prior[1], method = rows$method[1], cutoff = cutoffs,
               false_inclusions = sapply(cutoffs, function(c) mean(rowSums(m[, !true_model] > c))),
               missed           = sapply(cutoffs, function(c) mean(rowSums(m[, true_model] <= c))))
  }))))
at_half <- curves |> filter(abs(cutoff - 0.5) < 1e-9)
save_figure(
  ggplot(curves |> filter(model_prior == "beta-binomial"),
         aes(false_inclusions, missed, colour = method)) +
    geom_path(linewidth = 0.7) +
    geom_point(data = at_half |> filter(model_prior == "beta-binomial"), size = 2) +
    scale_colour_manual(values = method_colours) +
    facet_grid(predictor_correlation ~ sample_size,
               labeller = labeller(predictor_correlation = function(v) paste("correlation =", v),
                                   sample_size = function(v) paste("n =", v))) +
    theme(legend.position = "bottom") +
    labs(x = "false inclusions (of 8 null predictors)", y = "missed (of 4 real predictors)", colour = NULL,
         title = "Trade-off curves: every cutoff on the inclusion probability (beta-binomial model prior)",
         subtitle = "dots: the usual 0.5 cutoff; a curve lower and further left is better"),
  "09_trade_off_curves", width = 12, height = 7)

report_checks()
