## =====================================================================
##  Validation 05 (experiment): a per-coefficient test with the flat-R prior
##
##  Study 04 found that the joint flat-R test for ONE coefficient depends on
##  how much the other predictors explain, which costs discrimination when
##  predictors correlate. This script tries the alternative:
##
##    test one coefficient's partial correlation only, with the prior that the
##    joint flat-R prior implies for it.
##
##  In one sentence: same expectation per coefficient as the joint prior (the
##  multiplicity correction stays), but the evidence depends only on that
##  coefficient's own sample partial correlation, n and k.
##
##  The implied prior for one coefficient (derivation below the functions):
##    f_j^2 = share * R^2 / (1 - R^2),  R flat on (0, 1),
##    share ~ Beta(1/2, (k - 1)/2),     sign + or - with equal probability,
##  and the partial correlation is f_j / sqrt(1 + f_j^2). It does not depend on
##  how the predictors correlate (check C2). With k = 1 it is flat.
##
##  Checks first (C1-C3), then the comparison with JZS and the joint test on
##  the design of study 03 (x1 tested; slope 0, 0.1, 0.2, 0.3).
##
##  Run from the project folder: Rscript validation/05_conditional_test.R
##  About 10 minutes on 4 cores; the simulation is cached.
## =====================================================================
library(dplyr)
library(ggplot2)
source("validation/settings.R")   # parallel_map(), caching, record_check(), save_figure()
pkgload::load_all(".", quiet = TRUE)
source("validation/comparison_methods.R")


## ---- The test ------------------------------------------------------------
## Model: y = alpha + X b + noise. The other coefficients and alpha get flat
## priors, the noise 1/sigma (the same under both hypotheses, so they cancel).
## Integrating them out leaves, for the tested coefficient in noise units f,
##
##   p(y | f) / p(y | 0) = exp(-n f^2 / 2) * I(sqrt(n) r f) / I(0),
##   I(b) = integral_0^inf v^m exp(-v^2/2 + b v) dv,   m = n - k - 1,
##
## with r the sample partial correlation of the tested predictor. The Bayes
## factor is this ratio averaged over the prior for f.

## log I(b), for a vector of b: a grid around the peak of the integrand
log_peak_integral <- function(b, m, half_width = 12, points = 241) {
  peak    <- (b + sqrt(b^2 + 4 * m)) / 2
  width   <- 1 / sqrt(1 + m / peak^2)
  offsets <- seq(-half_width, half_width, length.out = points)
  v <- matrix(peak, length(b), points) + outer(width, offsets)
  v[v <= 0] <- NA
  log_integrand <- m * log(v) - v^2 / 2 + b * v
  largest <- apply(log_integrand, 1, max, na.rm = TRUE)
  largest + log(rowSums(exp(log_integrand - largest), na.rm = TRUE)) +
    log(width * (offsets[2] - offsets[1]))
}

## Prior points for |f|: R and the share on evenly spaced quantiles
## (deterministic, equal weights)
## Accuracy: about 1e-4 on the log scale against 4x as many points.
log_bayes_factor_flat_r_conditional <- function(partial_r, sample_size, number_of_predictors,
                                                points_R = 100, points_share = 50) {
  m <- sample_size - number_of_predictors - 1
  R <- (seq_len(points_R) - 0.5) / points_R
  share <- if (number_of_predictors == 1) 1 else
    qbeta((seq_len(points_share) - 0.5) / points_share, 0.5, (number_of_predictors - 1) / 2)
  f <- sqrt(as.vector(outer(R^2 / (1 - R^2), share)))
  b <- sqrt(sample_size) * partial_r * f
  log_null <- log_peak_integral(0, m)
  log_ratio <- c(-sample_size * f^2 / 2 + log_peak_integral( b, m) - log_null,
                 -sample_size * f^2 / 2 + log_peak_integral(-b, m) - log_null)
  largest <- max(log_ratio)
  largest + log(mean(exp(log_ratio - largest)))
}

## Sample partial correlation of predictor j, from its t-statistic
partial_correlation_of <- function(lm_fit, coefficient_row) {
  t_value <- summary(lm_fit)$coefficients[coefficient_row, "t value"]
  t_value / sqrt(t_value^2 + lm_fit$df.residual)
}


## ---- C1-C3: checks ---------------------------------------------------------
## C1. With one predictor the test is the joint flat-R test (both are flat on r).
c1_differences <- sapply(list(c(0.3, 30), c(0.1, 100), c(0.5, 50), c(0.05, 500)), function(s)
  abs(log_bayes_factor_flat_r_conditional(s[1], s[2], 1, points_R = 4000) -
      log_bayes_factor_flat_r_vs_intercept_only(s[1]^2, s[2], 1)))
record_check("C1. k = 1: equals the joint flat-R Bayes factor", max(c1_differences) < 1e-3,
             sprintf("largest difference %.1e", max(c1_differences)))

## C2. The implied prior for one coefficient's partial correlation does not
##     depend on how the predictors correlate: draws from the joint prior with
##     correlated predictors against the share formula.
set.seed(1)
correlated <- matrix(0.6, 4, 4); diag(correlated) <- 1
joint_draws <- draw_coefficients_from_flat_r_prior(200000, correlated)
spread_left <- sqrt(1 / diag(solve(correlated)))[2]
from_joint <- joint_draws[, 2] * spread_left / sqrt(1 + (joint_draws[, 2] * spread_left)^2)
f_formula  <- sqrt(rbeta(200000, 0.5, 1.5) * (function(R) R^2 / (1 - R^2))(runif(200000)))
from_formula <- sample(c(-1, 1), 200000, replace = TRUE) * f_formula / sqrt(1 + f_formula^2)
c2_distance <- suppressWarnings(unname(ks.test(from_joint, from_formula)$statistic))
record_check("C2. implied prior: correlated joint draws = share formula (k = 4)", c2_distance < 0.01,
             sprintf("KS distance %.4f", c2_distance))

## C3. With a flat prior on the partial correlation (share = 1), the test should
##     come close to the approximate flat partial r method of comparison_methods.R
##     (which plugs in the estimated noise). Criterion: within 0.05 on the log scale.
log_bf_flat_partial_exact <- function(partial_r, n, k, points = 4000) {
  m <- n - k - 1; R <- (seq_len(points) - 0.5) / points; f <- R / sqrt(1 - R^2)
  b <- sqrt(n) * partial_r * f; log_null <- log_peak_integral(0, m)
  log_ratio <- c(-n * f^2 / 2 + log_peak_integral(b, m) - log_null,
                 -n * f^2 / 2 + log_peak_integral(-b, m) - log_null)
  largest <- max(log_ratio); largest + log(mean(exp(log_ratio - largest)))
}
set.seed(3)
c3_differences <- sapply(c(60, 100, 200, 400), function(n) {
  X <- matrix(rnorm(n * 3), n) %*% chol(0.5 + 0.5 * diag(3))
  y <- drop(X %*% c(0.15, 0.3, 0)) + rnorm(n)
  fit <- lm(y ~ X)
  abs(log_bf_flat_partial_exact(partial_correlation_of(fit, 2), n, 3) -
      log(bayes_factor_flat_partial_correlation(fit, 1, X)))
})
record_check("C3. flat prior version vs the approximate flat partial r", max(c3_differences) < 0.05,
             sprintf("largest difference %.3f", max(c3_differences)))


## ---- The comparison: design of study 03 ------------------------------------
number_of_predictors   <- 4
slopes_of_x1           <- c(0, 0.1, 0.2, 0.3)
slopes_of_x2_to_x4     <- c(0.3, 0, 0)
sample_sizes           <- c(50, 100, 200, 500)
predictor_correlations <- c(0, 0.5)
replications_per_cell  <- 200

one_dataset <- function(slope_of_x1, sample_size, predictor_correlation) {
  covariance <- matrix(predictor_correlation, number_of_predictors, number_of_predictors)
  diag(covariance) <- 1
  slopes     <- c(slope_of_x1, slopes_of_x2_to_x4)
  predictors <- matrix(rnorm(sample_size * number_of_predictors), sample_size) %*% chol(covariance)
  outcome    <- drop(predictors %*% slopes) +
                rnorm(sample_size, sd = sqrt(1 - drop(t(slopes) %*% covariance %*% slopes)))
  full_fit <- lm(outcome ~ predictors)
  r_squared_full    <- summary(full_fit)$r.squared
  r_squared_reduced <- summary(lm(outcome ~ predictors[, -1]))$r.squared
  data.frame(
    method = c("JZS default", "joint flat-R", "conditional flat-R"),
    log_BF10 = c(
      log(bayes_factor_jzs_one_coefficient(r_squared_full, r_squared_reduced, sample_size,
                                           number_of_predictors, rscale = sqrt(2) / 4)),
      log_bayes_factor_flat_r_vs_intercept_only(r_squared_full, sample_size, number_of_predictors) -
        log_bayes_factor_flat_r_vs_intercept_only(r_squared_reduced, sample_size, number_of_predictors - 1),
      log_bayes_factor_flat_r_conditional(partial_correlation_of(full_fit, 2), sample_size,
                                          number_of_predictors)))
}

runs <- run_simulation_grid(
  tidyr::expand_grid(slope_of_x1 = slopes_of_x1, sample_size = sample_sizes,
                     predictor_correlation = predictor_correlations),
  results_file = file.path(results_folder, "sim_conditional_test.rds"),
  settings = list(replications_per_cell = replications_per_cell,
                  slopes_of_x2_to_x4 = slopes_of_x2_to_x4, version = 1),
  simulate_cell = function(cell)
    do.call(rbind, lapply(seq_len(replications_per_cell), function(replication)
      cbind(slope_of_x1 = cell$slope_of_x1, sample_size = cell$sample_size,
            predictor_correlation = cell$predictor_correlation, replication = replication,
            one_dataset(cell$slope_of_x1, cell$sample_size, cell$predictor_correlation)))))


## ---- Discrimination (threshold-free) and calibration (at 3) -----------------
area_under_curve <- function(higher, lower)
  mean(outer(higher, lower, ">")) + 0.5 * mean(outer(higher, lower, "=="))

summary_table <- runs |>
  group_by(predictor_correlation, sample_size, method) |>
  summarise(
    auc_null_vs_0.1     = area_under_curve(-log_BF10[slope_of_x1 == 0], -log_BF10[slope_of_x1 == 0.1]),
    null_found          = mean(log_BF10[slope_of_x1 == 0] < -log(3)),
    false_effect        = mean(log_BF10[slope_of_x1 == 0] >  log(3)),
    small_called_null   = mean(log_BF10[slope_of_x1 == 0.1] < -log(3)),
    effect_0.3_found    = mean(log_BF10[slope_of_x1 == 0.3] >  log(3)),
    .groups = "drop")
print(as.data.frame(summary_table), digits = 3)

## Paired bootstrap: conditional minus JZS, and conditional minus joint, in AUC
set.seed(4)
paired <- do.call(rbind, lapply(split(runs, list(runs$predictor_correlation, runs$sample_size)), function(cell) {
  get <- function(m, s) -cell$log_BF10[cell$method == m & cell$slope_of_x1 == s]
  difference <- function(i0, i1, other)
    area_under_curve(get("conditional flat-R", 0)[i0], get("conditional flat-R", 0.1)[i1]) -
    area_under_curve(get(other, 0)[i0], get(other, 0.1)[i1])
  reps <- length(get("JZS default", 0))
  do.call(rbind, lapply(c("JZS default", "joint flat-R"), function(other) {
    boot <- replicate(500, difference(sample(reps, replace = TRUE), sample(reps, replace = TRUE), other))
    data.frame(predictor_correlation = cell$predictor_correlation[1], sample_size = cell$sample_size[1],
               compared_with = other, auc_difference = difference(seq_len(reps), seq_len(reps), other),
               lower_95 = quantile(boot, 0.025), upper_95 = quantile(boot, 0.975), row.names = NULL)
  }))
}))
print(paired[order(paired$compared_with, paired$predictor_correlation, paired$sample_size), ], digits = 3)

save_figure(
  ggplot(summary_table, aes(sample_size, auc_null_vs_0.1, colour = method)) +
    geom_line(linewidth = 0.9) + geom_point() +
    scale_x_log10(breaks = sample_sizes) +
    scale_colour_manual(values = c(`JZS default` = "grey40", `joint flat-R` = "#d95f02",
                                   `conditional flat-R` = "#e7298a")) +
    facet_wrap(~ predictor_correlation,
               labeller = labeller(predictor_correlation = function(v) paste("predictor correlation =", v))) +
    labs(x = "sample size (log scale)", y = "AUC: null vs slope 0.1", colour = NULL,
         title = "Per-coefficient test: joint vs conditional flat-R",
         subtitle = "higher = separates nulls from small effects better, at any threshold") +
    theme(legend.position = "bottom"),
  "05_conditional_discrimination", width = 10, height = 5)

report_checks()
