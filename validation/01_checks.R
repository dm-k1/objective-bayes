## =====================================================================
##  Validation 01: is everything correct?
##  One script, one pass/fail summary. Every check has an explicit criterion;
##  the script stops with an error if anything failed.
##
##  PART A. Calculator: is the closed-form Bayes factor computed correctly?
##    The Bayes factor depends only on R-squared, n and k, so these inputs are
##    swept directly (extremes included) instead of simulating datasets.
##    A1. against the ground truth (adaptive integration with integrate())
##        over a sweep of R-squared, n and k; lists where the code's grid is
##        not accurate enough
##    A2. against brute force: draw from the prior and average (an
##        independent check of the formula itself), k = 1, 3, 5
##    A3. invariance: rescaling or recoding predictors leaves the BF unchanged
##    A4. reproduces the tutorial's Bayes factors (catches accidental changes)
##    A5. the brms prior and the closed form describe the same prior (bridge
##        sampling; optional, needs CmdStan): continuous predictor, and a
##        factor tested as a block
##
##  PART B. Default properties: does the method behave as a default should?
##    B1. one predictor: agrees with the per-coefficient flat-r Bayes factor
##    B2a. predictive matching (n = k + 1 gives BF = 1)
##    B2b. information consistency (BF grows without limit as R2 -> 1)
##    B2c. consistency under the null (log BF01 grows like 1/2 log n),
##         by simulating null data
##
##  Run from the project folder: Rscript validation/01_checks.R
##  A few minutes; the brms part (A5) is cached after the first run.
## =====================================================================
source("R/possible_range.R")
source("validation/comparison_methods.R")   # flat-r Bayes factor for B1
source("validation/settings.R")             # parallel cores, cache folders, check()

run_bridge_sampling_check <- TRUE    # set FALSE to skip A5 (needs CmdStan)

## brms fits from A5 are cached in validation/results/fits/ and reused on the
## next run. They are refitted automatically when the model code, the prior or
## the data change (file_refit = "on_change"). Delete the files to force a refit.
cache_folder <- fits_folder


## =====================================================================
## Ground truth: the same Bayes factor by adaptive integration
## =====================================================================
## R/possible_range.R averages over the scale g and over directions with
## fixed grids (fast). Here both averages are computed with integrate():
## adaptive quadrature that adds points where the curve changes quickly and
## reports its own error estimate. Slow, but accurate, so it serves as the
## judge in A1. Everything is on the log scale, scaled by the largest value,
## so that astronomically large Bayes factors do not overflow.

## ---- Average over directions, for one value of g ----------------------------
## Directions enter through the angle phi between a direction and the data's
## direction; only (0, pi/2) is needed (phi and pi - phi give the same t^2).
## The weight of each angle is sin(phi)^(k-2), normalized exactly.
log_direction_average_exact <- function(scaled_share, sample_size, number_of_predictors) {
  half_n_minus_one <- (sample_size - 1) / 2
  log_value_at_zero <- -half_n_minus_one * log1p(-scaled_share)   # largest possible value
  if (number_of_predictors == 1) return(log_value_at_zero)          # only t^2 = 1

  log_total_weight <- 0.5 * log(pi) - log(2) +
    lgamma((number_of_predictors - 1) / 2) - lgamma(number_of_predictors / 2)
  scaled_integrand <- function(angle)
    exp(-half_n_minus_one * log1p(-scaled_share * cos(angle)^2) - log_value_at_zero +
        (number_of_predictors - 2) * log(sin(angle)) - log_total_weight)

  ## When the evidence is strong the integrand is a narrow peak near angle 0,
  ## of width about sqrt((1 - z) / (2 m z)). Splitting the interval at a few
  ## multiples of that width makes sure integrate() sees the peak.
  peak_width <- if (scaled_share > 0) sqrt((1 - scaled_share) / (2 * half_n_minus_one * scaled_share)) else Inf
  inner_breaks <- c(3, 10, 40) * peak_width
  break_points <- c(0, inner_breaks[inner_breaks < pi / 2], pi / 2)
  pieces <- vapply(seq_len(length(break_points) - 1), function(i)
    integrate(scaled_integrand, break_points[i], break_points[i + 1],
              rel.tol = 1e-10, subdivisions = 1000L)$value, numeric(1))
  log_value_at_zero + log(sum(pieces))
}

## ---- Full Bayes factor vs intercept only -----------------------------------------
log_bayes_factor_ground_truth <- function(explained_share, sample_size, number_of_predictors) {
  if (number_of_predictors == 0) return(0)

  ## log of the integrand over u = log(g) (includes the Jacobian g)
  log_integrand <- function(log_g) vapply(log_g, function(u) {
    g <- exp(u)
    log(sample_size / 2) - 2 * u - sample_size / (2 * g) + u -      # IG(1, n/2) density, times g
      0.5 * log1p(g) +
      log_direction_average_exact(explained_share * g / (1 + g), sample_size, number_of_predictors)
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


## ---- Example data (same as the tutorial and vignette) ---------------------
set.seed(42)
sample_size <- 60
predictor_correlations <- matrix(0.4, 3, 3); diag(predictor_correlations) <- 1
predictors <- matrix(rnorm(sample_size * 3), sample_size) %*% chol(predictor_correlations)
colnames(predictors) <- c("x1", "x2", "x3")
outcome <- drop(predictors %*% c(0.3, 0.2, 0)) + rnorm(sample_size)
example_data <- data.frame(predictors, y = outcome)




## #####################################################################
##  PART A. CALCULATOR
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
  c(code         = log_bayes_factor_vs_intercept_only(s$explained_share, s$sample_size,
                                                      s$number_of_predictors),
    ground_truth = log_bayes_factor_ground_truth(s$explained_share, s$sample_size,
                                                 s$number_of_predictors))
}, mc.cores = number_of_cores)
sweep_settings <- cbind(sweep_settings, do.call(rbind, sweep_values))
sweep_settings$error <- abs(sweep_settings$code - sweep_settings$ground_truth)

ordinary <- abs(sweep_settings$ground_truth) < 50     # |log BF| < 50: the range that matters in practice
check(sprintf("A1a. ground truth, ordinary evidence (%d settings)", sum(ordinary)),
      max(sweep_settings$error[ordinary]) < 1e-3,
      sprintf("largest error %.1e", max(sweep_settings$error[ordinary])))
check(sprintf("A1b. ground truth, all settings (%d)", nrow(sweep_settings)),
      max(sweep_settings$error) < 1e-3,
      sprintf("largest error %.1e", max(sweep_settings$error)))

cat("\nSettings with the largest errors:\n")
print(head(sweep_settings[order(-sweep_settings$error), ], 6), row.names = FALSE, digits = 6)
cat("\n")


## =====================================================================
## A2. Closed form vs brute force
## =====================================================================
## Brute force uses no formula: draw coefficient sets from the prior, see
## how well each predicts the outcome (noise level averaged on a grid), and
## average. Criterion: agreement within 0.05 on the log scale.

draw_effects_from_prior <- function(number_of_draws, predictor_covariance) {
  number_of_predictors <- ncol(predictor_covariance)
  explained_share  <- rbeta(number_of_draws, 0.5, 1)
  overall_strength <- sqrt(explained_share / (1 - explained_share))
  random_numbers   <- matrix(rnorm(number_of_draws * number_of_predictors),
                             number_of_draws, number_of_predictors)
  directions       <- random_numbers / sqrt(rowSums(random_numbers^2))
  covariance_root  <- chol(predictor_covariance)
  t(backsolve(covariance_root, t(directions * overall_strength)))
}

log_bayes_factor_brute_force <- function(predictors, outcome,
                                         number_of_draws = 40000, draws_per_batch = 2000) {
  sample_size            <- length(outcome)
  centered_predictors    <- scale(predictors, scale = FALSE)
  centered_outcome       <- outcome - mean(outcome)
  outcome_sum_of_squares <- sum(centered_outcome^2)
  predictor_covariance   <- crossprod(centered_predictors) / sample_size

  ## the noise level is unknown: average over it on a grid (as 1 / noise)
  inverse_noise_values <- seq(1e-4, 6 * sqrt(sample_size / outcome_sum_of_squares), length.out = 800)
  inverse_noise_step   <- inverse_noise_values[2] - inverse_noise_values[1]
  log_inverse_noise    <- log(inverse_noise_values)

  log_evidence_null <- log_sum_exp((sample_size - 2) * log_inverse_noise -
                                   inverse_noise_values^2 * outcome_sum_of_squares / 2) +
                       log(inverse_noise_step)

  log_evidence_per_draw <- numeric(0)
  for (batch in seq_len(number_of_draws / draws_per_batch)) {
    effects          <- draw_effects_from_prior(draws_per_batch, predictor_covariance)
    predicted_values <- centered_predictors %*% t(effects)
    overlap_with_outcome     <- drop(crossprod(centered_outcome, predicted_values))
    predicted_sum_of_squares <- colSums(predicted_values^2)
    log_integrand <- (sample_size - 2) *
        matrix(log_inverse_noise, draws_per_batch, length(inverse_noise_values), byrow = TRUE) -
      (outer(rep(outcome_sum_of_squares, draws_per_batch), inverse_noise_values^2) -
       2 * outer(overlap_with_outcome, inverse_noise_values) + predicted_sum_of_squares) / 2
    log_evidence_per_draw <- c(log_evidence_per_draw,
                               log_sum_exp_per_row(log_integrand) + log(inverse_noise_step))
  }
  log_sum_exp(log_evidence_per_draw) - log(number_of_draws) - log_evidence_null
}

brute_force_cases <- list(
  list(label = "k = 3, tutorial data", predictors = predictors, outcome = outcome),
  local({ set.seed(11); x <- matrix(rnorm(80), 80, 1)
          list(label = "k = 1, n = 80", predictors = x, outcome = drop(0.3 * x) + rnorm(80)) }),
  local({ set.seed(12); x <- matrix(rnorm(100 * 5), 100) %*% chol(0.3 + 0.7 * diag(5))
          list(label = "k = 5, n = 100, correlated", predictors = x,
               outcome = drop(x %*% c(0.2, 0, 0.1, 0, 0)) + rnorm(100)) }))

set.seed(1)
for (case in brute_force_cases) {
  closed_form <- log_bayes_factor_vs_intercept_only(
    explained_share_of(case$predictors, case$outcome), length(case$outcome), ncol(case$predictors))
  brute_force <- log_bayes_factor_brute_force(case$predictors, case$outcome)
  check(paste("A2. brute force,", case$label),
        abs(closed_form - brute_force) < 0.05,
        sprintf("closed %.3f, brute %.3f", closed_form, brute_force))
}


## =====================================================================
## A3. Invariance
## =====================================================================
## The BF should not depend on units or on how the same predictors are
## written. Criterion: identical up to 1e-8 (relative).

bf_original <- possible_bayes_factor_per_term(y ~ x1 + x2 + x3, example_data)$BF10

rescaled_data <- transform(example_data, x1 = 100 * x1 + 7, x2 = -0.01 * x2, x3 = x3 - 50)
bf_rescaled <- possible_bayes_factor_per_term(y ~ x1 + x2 + x3, rescaled_data)$BF10
check("A3a. rescaling and shifting predictors", max(abs(bf_rescaled / bf_original - 1)) < 1e-8,
      sprintf("largest relative difference %.1e", max(abs(bf_rescaled / bf_original - 1))))

outcome_rescaled <- transform(example_data, y = 1000 * y - 3)
bf_outcome <- possible_bayes_factor_per_term(y ~ x1 + x2 + x3, outcome_rescaled)$BF10
check("A3b. rescaling the outcome", max(abs(bf_outcome / bf_original - 1)) < 1e-8,
      sprintf("largest relative difference %.1e", max(abs(bf_outcome / bf_original - 1))))

recoded_data <- transform(example_data, sum_12 = x1 + x2, difference_12 = x1 - x2)
bf_full_original <- possible_bayes_factor(lm(y ~ x1 + x2 + x3, example_data), lm(y ~ 1, example_data))
bf_full_recoded  <- possible_bayes_factor(lm(y ~ sum_12 + difference_12 + x3, recoded_data),
                                          lm(y ~ 1, recoded_data))
check("A3c. recoding predictors (same model, other columns)",
      abs(bf_full_recoded / bf_full_original - 1) < 1e-8,
      sprintf("original %.4f, recoded %.4f", bf_full_original, bf_full_recoded))


## =====================================================================
## A4. Reproduces the tutorial
## =====================================================================
## Criterion: within 0.1% of the recorded values 18.696, 0.3827, 1.5872
## (computed with the packed direction grid).
recorded <- c(x1 = 18.696, x2 = 0.3827, x3 = 1.5872)
check("A4. tutorial Bayes factors", all(abs(bf_original / recorded - 1) < 1e-3),
      paste(sprintf("%.4f", bf_original), collapse = " "))


## =====================================================================
## A5. Closed form vs brms bridge sampling
## =====================================================================
## The brms prior (a Stan density on the coefficients) and the closed form
## must describe the same prior. Bridge sampling is noisy: criterion is
## agreement within 10% (0.1 on the log scale).
## Case A: a continuous predictor (x1). Case B: a 3-level factor tested
## as a block, so dummy-coded columns are covered too.

if (run_bridge_sampling_check) {
  library(brms)   # save_pars(), bayes_factor()
  set.seed(5)
  example_data$group   <- factor(sample(c("a", "b", "c"), sample_size, replace = TRUE))
  example_data$y_group <- example_data$y + c(a = 0, b = 0.5, c = -0.3)[as.character(example_data$group)]

  bridge_cases <- list(
    list(label = "A5a. bridge sampling, continuous x1", cache_name = "checkA5A",
         full = y ~ x1 + x2 + x3, reduced = y ~ x2 + x3),
    list(label = "A5b. bridge sampling, factor as a block", cache_name = "checkA5B",
         full = y_group ~ x1 + group, reduced = y_group ~ x1))

  for (case in bridge_cases) {
    full_fit <- brm_possible(case$full, example_data, iter = 10000, seed = 2, refresh = 0,
                             save_pars = save_pars(all = TRUE),
                             control = list(adapt_delta = 0.99),
                             file = file.path(cache_folder, paste0(case$cache_name, "_full")),
                             file_refit = "on_change")
    reduced_fit <- brm_possible(case$reduced, example_data, iter = 10000, seed = 2, refresh = 0,
                                save_pars = save_pars(all = TRUE),
                                control = list(adapt_delta = 0.99),
                                file = file.path(cache_folder, paste0(case$cache_name, "_reduced")),
                                file_refit = "on_change")
    bridge <- bayes_factor(full_fit, reduced_fit, silent = TRUE)$bf
    closed <- possible_bayes_factor(lm(case$full, example_data), lm(case$reduced, example_data))
    check(case$label, abs(log(bridge) - log(closed)) < 0.1,
          sprintf("bridge %.3f, closed %.3f", bridge, closed))
  }
}



## #####################################################################
##  PART B. DEFAULT PROPERTIES
## #####################################################################

## =====================================================================
## B1. One predictor: agrees with the per-coefficient flat-r Bayes factor
## =====================================================================
## The closed form with k = 1 is flat on the correlation. The per-coefficient
## version (validation/comparison_methods.R) is an independent calculation of
## the same idea, but it plugs in the estimated noise level, so it is an
## approximation. Criterion: agreement within 3% (relative).

set.seed(3)
for (true_slope in c(0, 0.2, 0.4)) {
  x <- rnorm(100); y <- true_slope * x + rnorm(100)
  one_predictor_fit <- lm(y ~ x)
  closed_form <- exp(log_bayes_factor_vs_intercept_only(summary(one_predictor_fit)$r.squared, 100, 1))
  flat_r      <- bayes_factor_flat_partial_r(one_predictor_fit, 1, matrix(x))
  check(sprintf("B1. k = 1 vs flat partial r, slope %.1f", true_slope),
        abs(closed_form / flat_r - 1) < 0.03,
        sprintf("closed %.3f, flat r %.3f", closed_form, flat_r))
}


## =====================================================================
## B2. Properties of a good default
## =====================================================================

## a) Predictive matching: with one more observation than predictors the
##    data can't tell the models apart, so BF = 1. Criterion: within 0.001.
matching <- sapply(1:6, function(k)
  exp(log_bayes_factor_vs_intercept_only(1 - 1e-12, k + 1, k)))
check("B2a. predictive matching (k = 1-6)", all(abs(matching - 1) < 1e-3),
      paste(sprintf("%.4f", matching), collapse = " "))

## b) Information consistency: as the fit becomes perfect, evidence for an
##    effect grows without limit. Criterion: strictly increasing, and above
##    log BF = 50 at R2 = 0.9999.
information <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(r2)
  log_bayes_factor_vs_intercept_only(r2, 20, 3))
check("B2b. information consistency", all(diff(information) > 0) && information[4] > 50,
      paste(sprintf("%.1f", information), collapse = " "))

## c) Consistency under the null, by simulation: x3 truly does nothing.
##    The median log BF01 for x3 should grow with n at a rate of about
##    1/2 log n. Criterion: increasing in n, and slope against log n
##    between 0.35 and 0.65.
set.seed(4)
sample_sizes <- c(50, 200, 1000, 5000)
median_log_bf01 <- sapply(sample_sizes, function(n) {
  median(replicate(200, {
    x <- matrix(rnorm(n * 3), n)
    y <- drop(x %*% c(0.3, 0.2, 0)) + rnorm(n)
    -(log_bayes_factor_vs_intercept_only(explained_share_of(x, y), n, 3) -
      log_bayes_factor_vs_intercept_only(explained_share_of(x[, 1:2], y), n, 2))
  }))
})
slope_against_log_n <- coef(lm(median_log_bf01 ~ log(sample_sizes)))[[2]]
check("B2c. consistency under the null (simulated)",
      all(diff(median_log_bf01) > 0) && slope_against_log_n > 0.35 && slope_against_log_n < 0.65,
      sprintf("median log BF01 %s; slope %.2f",
              paste(sprintf("%.2f", median_log_bf01), collapse = " "), slope_against_log_n))


## =====================================================================
## Summary
## =====================================================================
report_checks()
