## =====================================================================
##  Validation 08: consistency when the number of predictors grows with n
##
##  Model selection consistency (criterion 2 of Bayarri et al., 2012) is
##  usually shown for a fixed number of predictors. With many predictors,
##  k grows with n, and some default priors then fail: the evidence for the
##  true model stops growing, or goes the wrong way.
##
##  Two growth rates: k = n / 4 (proportional) and k = sqrt(n).
##    G1. The null is true (no predictor matters). Compare the k-predictor
##        model with the null: log BF01 should grow without limit.
##        R-squared of a null model is Beta(k/2, (n-k-1)/2), so no data
##        needs to be simulated.
##    G2. A small model is true (3 predictors, 0.3 each, noise SD 1); the
##        other k - 3 predictors do nothing. Compare the true model with the
##        full k-predictor model: log BF(true vs full) should grow without limit.
##  Methods: the whole-model set of validation/comparison_methods.R: the three
##  contenders, flat partial r (R-squared ~ Beta(k/2, 1)), flat R
##  (Beta(1/2, 1)) and flat R-squared (Beta(1, 1)), each with a check; and the
##  two references, the JZS default (scale 0.354) and Zellner-Siow (scale 1),
##  reported for comparison.
##
##  Run from the project folder: Rscript validation/08_growing_k.R
##  About 10-20 minutes (the n = 2000, k = 500 regressions in G2 dominate).
## =====================================================================
library(dplyr)
library(ggplot2)
source("validation/settings.R")   # parallel_map(), record_check(), save_figure()
pkgload::load_all(".", quiet = TRUE)
source("validation/comparison_methods.R")   # the common method set

sample_sizes <- c(50, 100, 200, 500, 1000, 2000)
growth_rates <- list(`k = n / 4` = function(n) floor(n / 4), `k = sqrt(n)` = function(n) floor(sqrt(n)))
replications <- 200

log_bf_vs_null <- log_bayes_factor_model
methods <- model_methods


## ---- G1: the null is true ---------------------------------------------------
set.seed(1)
null_cells <- expand.grid(n = sample_sizes, growth = names(growth_rates), stringsAsFactors = FALSE)
null_results <- do.call(rbind, parallel_map(seq_len(nrow(null_cells)), function(row) {
  n <- null_cells$n[row]; k <- growth_rates[[null_cells$growth[row]]](n)
  set.seed(100 + row)
  r_squared <- rbeta(replications, k / 2, (n - k - 1) / 2)
  do.call(rbind, lapply(methods, function(method)
    data.frame(growth = null_cells$growth[row], n = n, k = k, method = method,
               median_log_bf01 = -median(sapply(r_squared, log_bf_vs_null, method = method,
                                                n = n, k = k)))))
}))


## ---- G2: a small model is true, extra predictors do nothing -----------------
true_slopes <- c(0.3, 0.3, 0.3)
alternative_results <- do.call(rbind, parallel_map(seq_len(nrow(null_cells)), function(row) {
  n <- null_cells$n[row]; k <- growth_rates[[null_cells$growth[row]]](n)
  set.seed(200 + row)
  log_bf_true_vs_full <- replicate(replications, {
    predictors <- matrix(rnorm(n * k), n)
    outcome <- drop(predictors[, 1:3] %*% true_slopes) + rnorm(n)
    r_squared_true <- sample_r_squared(predictors[, 1:3], outcome)
    r_squared_full <- sample_r_squared(predictors, outcome)
    sapply(methods, function(method)
      log_bf_vs_null(method, r_squared_true, n, 3) - log_bf_vs_null(method, r_squared_full, n, k))
  })
  data.frame(growth = null_cells$growth[row], n = n, k = k, method = methods,
             median_log_bf_true_vs_full = apply(log_bf_true_vs_full, 1, median), row.names = NULL)
}))


## ---- Results -------------------------------------------------------------------
print(null_results, digits = 3)
print(alternative_results, digits = 3)

## A method passes when the median evidence for the truth keeps rising with n:
## slope against log(n) above 0.1, and higher at the largest n than at the
## smallest. (Not a fixed height: with k growing with n the evidence can rise
## slowly and still be consistent; the slope is reported.)
growth_slope <- function(values) coef(lm(values ~ log(sample_sizes)))[[2]]
check_growth <- function(values) growth_slope(values) > 0.1 && tail(values, 1) > values[1]
for (growth in names(growth_rates)) for (method in methods) {
  g1 <- null_results$median_log_bf01[null_results$growth == growth & null_results$method == method]
  g2 <- alternative_results$median_log_bf_true_vs_full[alternative_results$growth == growth &
                                                         alternative_results$method == method]
  if (method %in% names(contenders)) {
    record_check(sprintf("G1. %s, %s, null true: evidence for the null grows", method, growth),
                 check_growth(g1),
                 sprintf("median log BF01 at n = 50 ... 2000: %s; slope %.2f",
                         paste(sprintf("%.1f", g1), collapse = " "), growth_slope(g1)))
    record_check(sprintf("G2. %s, %s, small model true: evidence for it grows", method, growth),
                 check_growth(g2),
                 sprintf("median log BF at n = 50 ... 2000: %s; slope %.2f",
                         paste(sprintf("%.1f", g2), collapse = " "), growth_slope(g2)))
  } else {
    message(sprintf("%s, %s (comparison): G1 %s; G2 %s", method, growth,
                    paste(sprintf("%.1f", g1), collapse = " "), paste(sprintf("%.1f", g2), collapse = " ")))
  }
}

save_figure(
  ggplot(rbind(transmute(null_results, growth, n, method, value = median_log_bf01,
                         question = "G1: null true, log BF(null vs k-predictor model)"),
               transmute(alternative_results, growth, n, method, value = median_log_bf_true_vs_full,
                         question = "G2: 3-predictor model true, log BF(true vs full)")),
         aes(n, value, colour = method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 0.9) + geom_point() +
    scale_x_log10(breaks = sample_sizes) +
    scale_colour_manual(values = method_colours[model_methods]) +
    facet_grid(question ~ growth, scales = "free_y") + theme(legend.position = "bottom") +
    labs(x = "n (log scale)", y = "median log Bayes factor for the truth", colour = NULL,
         title = "Consistency when the number of predictors grows with n",
         subtitle = "consistent: keeps rising; below 0: favours the wrong model"),
  "08_growing_k", width = 11, height = 7)

report_checks()
