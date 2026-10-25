## =====================================================================
##  Validation 02: is the Bayes factor right?
##  One script, one pass/fail summary. Every check has an explicit criterion;
##  the script stops with an error if anything failed.
##
##  PART A. Correctness
##    A1. the grid computes the formula: against adaptive integration
##        (integrate()) over a sweep of R-squared, n and k, extremes included
##    A2. invariance: rescaling or recoding predictors leaves the BF unchanged
##    A3. the formula is the right one: closed form vs brms bridge sampling
##        with the Stan prior (continuous predictor; factor tested as a block)
##
##  PART B. Properties of a good default (Bayarri et al., 2012)
##    B1. predictive matching (n = k + 1 gives BF = 1)
##    B2. information consistency (BF grows without limit as R2 -> 1)
##    B3. consistency under the null (log BF01 grows like 1/2 log n), simulated
##
##  Run from the project folder: Rscript validation/02_bayes_factor.R
##  A few minutes; the brms fits (A3) are cached after the first run.
## =====================================================================
source("R/flat_r.R")
source("validation/settings.R")             # parallel cores, cache folders, record_check()

run_bridge_sampling_check <- TRUE    # set FALSE to skip A3 (needs CmdStan)

## brms fits from A3 are cached in validation/results/fits/ and reused on the
## next run. They are refitted automatically when the model code, the prior or
## the data change (file_refit = "on_change"). Delete the files to force a refit.
cache_folder <- fits_folder


## =====================================================================
## Ground truth: the same Bayes factor by adaptive integration
## =====================================================================
## R/flat_r.R averages over the scale g and over directions with
## fixed grids (fast). Here both averages are computed with integrate():
## adaptive quadrature that adds points where the curve changes quickly and
## reports its own error estimate. Slow, but accurate, so it serves as the
## judge in A1. Everything is on the log scale, scaled by the largest value,
## so that astronomically large Bayes factors do not overflow.

## ---- Average over directions, for one value of g ----------------------------
## Directions enter through the angle phi between a direction and the data's
## direction; only (0, pi/2) is needed (phi and pi - phi give the same alignment c^2).
## The weight of each angle is sin(phi)^(k-2), normalized exactly.
log_direction_average_by_integration <- function(scaled_share, sample_size, number_of_predictors) {
  half_n_minus_one <- (sample_size - 1) / 2
  log_value_at_zero <- -half_n_minus_one * log1p(-scaled_share)   # largest possible value
  if (number_of_predictors == 1) return(log_value_at_zero)          # only c^2 = 1

  log_total_weight <- 0.5 * log(pi) - log(2) +
    lgamma((number_of_predictors - 1) / 2) - lgamma(number_of_predictors / 2)
  scaled_integrand <- function(angle)
    exp(-half_n_minus_one * log1p(-scaled_share * cos(angle)^2) - log_value_at_zero +
        (number_of_predictors - 2) * log(sin(angle)) - log_total_weight)

  ## When the evidence is strong the integrand is a narrow peak near angle 0,
  ## of width about sqrt((1 - z) / (2 m z)). Splitting the interval at a few
  ## multiples of that width makes sure integrate() sees the peak.
  peak_width <- if (scaled_share > 0)
    sqrt((1 - scaled_share) / (2 * half_n_minus_one * scaled_share)) else Inf
  inner_breaks <- c(3, 10, 40) * peak_width
  break_points <- c(0, inner_breaks[inner_breaks < pi / 2], pi / 2)
  pieces <- vapply(seq_len(length(break_points) - 1), function(i)
    integrate(scaled_integrand, break_points[i], break_points[i + 1],
              rel.tol = 1e-10, subdivisions = 1000L)$value, numeric(1))
  log_value_at_zero + log(sum(pieces))
}

## ---- Full Bayes factor vs intercept only -----------------------------------------
log_bayes_factor_flat_r_by_integration <- function(explained_share, sample_size,
                                                           number_of_predictors) {
  if (number_of_predictors == 0) return(0)

  ## log of the integrand over log(g) (includes the factor g, because the step is in log g)
  log_integrand <- function(log_g_values) vapply(log_g_values, function(log_g) {
    g <- exp(log_g)
    log(sample_size / 2) - 2 * log_g - sample_size / (2 * g) + log_g -   # IG(1, n/2) density, times g
      0.5 * log1p(g) +
      log_direction_average_by_integration(explained_share * g / (1 + g), sample_size,
                                           number_of_predictors)
  }, numeric(1))

  ## find the peak, then integrate exp(log_integrand - peak) around it
  peak_location <- optimize(log_integrand, c(-20, 40), maximum = TRUE)$maximum
  peak_height   <- log_integrand(peak_location)
  scaled <- function(log_g) exp(log_integrand(log_g) - peak_height)
  break_points <- c(-40, peak_location - 3, peak_location + 3, 60)
  pieces <- vapply(1:3, function(i)
    integrate(scaled, break_points[i], break_points[i + 1],
              rel.tol = 1e-10, subdivisions = 1000L)$value, numeric(1))
  peak_height + log(sum(pieces))
}


## ---- Example data (same as the vignette) ---------------------
set.seed(42)
sample_size <- 60
predictor_correlations <- matrix(0.4, 3, 3); diag(predictor_correlations) <- 1
predictors <- matrix(rnorm(sample_size * 3), sample_size) %*% chol(predictor_correlations)
colnames(predictors) <- c("x1", "x2", "x3")
outcome <- drop(predictors %*% c(0.3, 0.2, 0)) + rnorm(sample_size)
example_data <- data.frame(predictors, y = outcome)


## #####################################################################
##  PART A. CORRECTNESS
## #####################################################################

## =====================================================================
## A1. Against the ground truth, over a sweep of inputs
## =====================================================================
## Criterion: the code's log Bayes factor within 0.001 of the ground truth for
## every setting. The worst settings are listed, so you can see where the
## grid is not accurate enough.

sweep_settings <- expand.grid(explained_share      = c(0.01, 0.1, 0.3, 0.6, 0.9, 0.99),
                              sample_size          = c(10, 30, 100, 1000, 10000),
                              number_of_predictors = c(1, 2, 3, 5, 10, 30))
sweep_settings <- subset(sweep_settings, sample_size > number_of_predictors + 1)

message("Ground truth for ", nrow(sweep_settings), " settings (adaptive integration) ...")
sweep_values <- parallel::mclapply(seq_len(nrow(sweep_settings)), function(row) {
  s <- sweep_settings[row, ]
  c(code         = log_bayes_factor_flat_r_vs_intercept_only(s$explained_share, s$sample_size,
                                                                     s$number_of_predictors),
    ground_truth = log_bayes_factor_flat_r_by_integration(s$explained_share, s$sample_size,
                                                                  s$number_of_predictors))
}, mc.cores = number_of_cores)
sweep_settings <- cbind(sweep_settings, do.call(rbind, sweep_values))
sweep_settings$error <- abs(sweep_settings$code - sweep_settings$ground_truth)

ordinary <- abs(sweep_settings$ground_truth) < 50     # |log BF| < 50: the range that matters in practice
record_check(sprintf("A1a. ground truth, ordinary evidence (%d settings)", sum(ordinary)),
             max(sweep_settings$error[ordinary]) < 1e-3,
             sprintf("largest error %.1e", max(sweep_settings$error[ordinary])))
record_check(sprintf("A1b. ground truth, all settings (%d)", nrow(sweep_settings)),
             max(sweep_settings$error) < 1e-3,
             sprintf("largest error %.1e", max(sweep_settings$error)))

cat("\nSettings with the largest errors:\n")
print(head(sweep_settings[order(-sweep_settings$error), ], 6), row.names = FALSE, digits = 6)
cat("\n")


## =====================================================================
## A2. Invariance
## =====================================================================
## The BF should not depend on units or on how the same predictors are
## written. Criterion: identical up to 1e-8 (relative).

bf_original <- bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, example_data)$BF10

rescaled_data <- transform(example_data, x1 = 100 * x1 + 7, x2 = -0.01 * x2, x3 = x3 - 50)
bf_rescaled <- bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, rescaled_data)$BF10
record_check("A2a. rescaling and shifting predictors", max(abs(bf_rescaled / bf_original - 1)) < 1e-8,
             sprintf("largest relative difference %.1e", max(abs(bf_rescaled / bf_original - 1))))

outcome_rescaled <- transform(example_data, y = 1000 * y - 3)
bf_outcome <- bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, outcome_rescaled)$BF10
record_check("A2b. rescaling the outcome", max(abs(bf_outcome / bf_original - 1)) < 1e-8,
             sprintf("largest relative difference %.1e", max(abs(bf_outcome / bf_original - 1))))

recoded_data <- transform(example_data, sum_12 = x1 + x2, difference_12 = x1 - x2)
bf_full_original <- bayes_factor_flat_r(lm(y ~ x1 + x2 + x3, example_data),
                                                lm(y ~ 1, example_data))
bf_full_recoded  <- bayes_factor_flat_r(lm(y ~ sum_12 + difference_12 + x3, recoded_data),
                                                lm(y ~ 1, recoded_data))
record_check("A2c. recoding predictors (same model, other columns)",
             abs(bf_full_recoded / bf_full_original - 1) < 1e-8,
             sprintf("original %.4f, recoded %.4f", bf_full_original, bf_full_recoded))


## =====================================================================
## A3. Closed form vs brms bridge sampling
## =====================================================================
## The brms prior (a Stan density on the coefficients) and the closed form
## must describe the same prior. This is the check of the formula itself
## (A1 only checks that the grid computes the formula accurately). Bridge sampling is noisy: criterion is
## agreement within 10% (0.1 on the log scale).
## Case A: a continuous predictor (x1). Case B: a 3-level factor tested
## as a block, so dummy-coded columns are covered too.

if (run_bridge_sampling_check) {
  library(brms)   # save_pars(), bayes_factor()
  set.seed(5)
  example_data$group   <- factor(sample(c("a", "b", "c"), sample_size, replace = TRUE))
  example_data$y_group <- example_data$y + c(a = 0, b = 0.5, c = -0.3)[as.character(example_data$group)]

  bridge_cases <- list(
    list(label = "A3a. bridge sampling, continuous x1", cache_name = "checkA3a",
         full = y ~ x1 + x2 + x3, reduced = y ~ x2 + x3),
    list(label = "A3b. bridge sampling, factor as a block", cache_name = "checkA3b",
         full = y_group ~ x1 + group, reduced = y_group ~ x1))

  ## same sampler settings for every fit; cached per case and model
  fit_for_bridge_sampling <- function(formula, cache_file)
    fit_flat_r_model(formula, example_data, iter = 10000, seed = 2, refresh = 0,
                             save_pars = save_pars(all = TRUE),
                             control = list(adapt_delta = 0.99),
                             file = file.path(cache_folder, cache_file), file_refit = "on_change")

  for (case in bridge_cases) {
    full_fit    <- fit_for_bridge_sampling(case$full,    paste0(case$cache_name, "_full"))
    reduced_fit <- fit_for_bridge_sampling(case$reduced, paste0(case$cache_name, "_reduced"))
    bridge <- bayes_factor(full_fit, reduced_fit, silent = TRUE)$bf
    closed <- bayes_factor_flat_r(lm(case$full, example_data), lm(case$reduced, example_data))
    record_check(case$label, abs(log(bridge) - log(closed)) < 0.1,
                 sprintf("bridge %.3f, closed %.3f", bridge, closed))
  }
}





## #####################################################################
##  PART B. PROPERTIES OF A GOOD DEFAULT
## #####################################################################

## =====================================================================
## B. Properties of a good default
## =====================================================================

## B1. Predictive matching: with one more observation than predictors the
##    data can't tell the models apart, so BF = 1. Criterion: within 0.001.
matching <- sapply(1:6, function(k)
  exp(log_bayes_factor_flat_r_vs_intercept_only(1 - 1e-12, k + 1, k)))
record_check("B1. predictive matching (k = 1-6)", all(abs(matching - 1) < 1e-3),
             paste(sprintf("%.4f", matching), collapse = " "))

## B2. Information consistency: as the fit becomes perfect, evidence for an
##    effect grows without limit. Criterion: strictly increasing, and above
##    log BF = 50 at R2 = 0.9999.
information <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(r2)
  log_bayes_factor_flat_r_vs_intercept_only(r2, 20, 3))
record_check("B2. information consistency", all(diff(information) > 0) && information[4] > 50,
             paste(sprintf("%.1f", information), collapse = " "))

## B3. Consistency under the null, by simulation: x3 truly does nothing.
##    The median log BF01 for x3 should grow with n at a rate of about
##    1/2 log n. Criterion: increasing in n, and slope against log n
##    between 0.35 and 0.65.
set.seed(4)
sample_sizes <- c(50, 200, 1000, 5000)
median_log_bf01 <- sapply(sample_sizes, function(n) {
  median(replicate(200, {
    x <- matrix(rnorm(n * 3), n)
    y <- drop(x %*% c(0.3, 0.2, 0)) + rnorm(n)
    -(log_bayes_factor_flat_r_vs_intercept_only(sample_r_squared(x, y), n, 3) -
      log_bayes_factor_flat_r_vs_intercept_only(sample_r_squared(x[, 1:2], y), n, 2))
  }))
})
slope_against_log_n <- coef(lm(median_log_bf01 ~ log(sample_sizes)))[[2]]
record_check("B3. consistency under the null (simulated)",
             all(diff(median_log_bf01) > 0) && slope_against_log_n > 0.35 && slope_against_log_n < 0.65,
             sprintf("median log BF01 %s; slope %.2f",
                     paste(sprintf("%.2f", median_log_bf01), collapse = " "), slope_against_log_n))


## =====================================================================
## Summary
## =====================================================================
report_checks()
