## =====================================================================
##  Validation 02: is the Bayes factor right?
##  One script, one pass/fail summary. Every check has an explicit criterion;
##  the script stops with an error if anything failed.
##
##  PART A. Correctness
##    A1. the grid computes the formula: against adaptive integration
##        (integrate()) over a sweep of R-squared, n and k, extremes included
##    A2. invariance: rescaling or recoding predictors leaves the BF unchanged
##    A3. the formula is the right one: closed form vs brute force (draw
##        coefficients from the prior, average how well they predict; no formula)
##    A4. a second derivation: the prior as a mixture of g-priors, whose Bayes
##        factor is the standard g-prior one averaged over the hyperprior
##
##  PART B. Properties of a good default (Bayarri et al., 2012)
##    B1. predictive matching (n = k + 1 gives BF = 1)
##    B2. information consistency (BF grows without limit as R2 -> 1)
##    B2b. ... exactly from n = k + 3 on (the analytic boundary)
##    B3. consistency under the null (log BF01 grows like 1/2 log n), simulated
##    B3b. ... the 1/2 log n holds when the smaller model is the null model;
##         with a real effect in the smaller model it is q/2 log n (q extra nulls)
##    B4. consistency under the alternative (log BF10 grows with n), simulated
##
##  Run from the project folder: Rscript validation/02_bayes_factor.R
##  About a minute. No Stan needed.
## =====================================================================
pkgload::load_all(".", quiet = TRUE)   # the flatr package, from this folder
source("validation/settings.R")             # parallel cores, cache folders, record_check()


## =====================================================================
## Ground truth: the same Bayes factor by adaptive integration
## =====================================================================
## R/bayes_factor.R averages over the scale g and over directions with
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
sweep_values <- parallel_map(seq_len(nrow(sweep_settings)), function(row) {
  s <- sweep_settings[row, ]
  c(code         = log_bayes_factor_flat_r_vs_intercept_only(s$explained_share, s$sample_size,
                                                                     s$number_of_predictors),
    ground_truth = log_bayes_factor_flat_r_by_integration(s$explained_share, s$sample_size,
                                                                  s$number_of_predictors))
})
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
## A3. Closed form vs brute force
## =====================================================================
## A1 checks that the grid computes the formula. This checks the formula
## itself, without using it: draw coefficient sets from the prior (the recipe
## of docs/2 section 5), see how well each predicts the outcome (the noise
## level averaged on a grid), and average. Monte Carlo, so noisy: with
## 400,000 draws its SD is about 0.02 on the log scale (k = 5).
## Criterion: agreement within 0.05 on the log scale.

## The same log Bayes factor as log_bayes_factor_flat_r_vs_intercept_only(),
## without the formula: draw coefficient sets from the prior, see how well
## each predicts the outcome (the noise level averaged on a grid), and
## average. Slow and noisy (Monte Carlo); agrees to about 0.05 on the log scale.
log_bayes_factor_flat_r_by_sampling <- function(predictors, outcome,
                                                        number_of_draws = 40000,
                                                        draws_per_batch = 2000) {
  sample_size            <- length(outcome)
  centered_predictors    <- scale(predictors, scale = FALSE)
  centered_outcome       <- outcome - mean(outcome)
  outcome_sum_of_squares <- sum(centered_outcome^2)
  predictor_covariance   <- crossprod(centered_predictors) / sample_size

  ## the noise level is unknown: average over it on a grid (written as
  ## 1 / noise, which makes the arithmetic simpler)
  inverse_noise_values <- seq(1e-4, 6 * sqrt(sample_size / outcome_sum_of_squares), length.out = 800)
  inverse_noise_step   <- inverse_noise_values[2] - inverse_noise_values[1]
  log_inverse_noise    <- log(inverse_noise_values)

  ## how well "no effect" predicts the data
  log_evidence_null <- log_sum_exp((sample_size - 2) * log_inverse_noise -
                                   inverse_noise_values^2 * outcome_sum_of_squares / 2) +
                       log(inverse_noise_step)

  ## how well each drawn coefficient set predicts the data
  log_evidence_per_draw <- numeric(0)
  for (batch in seq_len(number_of_draws / draws_per_batch)) {
    coefficients     <- draw_coefficients_from_flat_r_prior(draws_per_batch, predictor_covariance)
    predicted_values <- centered_predictors %*% t(coefficients)
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
  list(label = "k = 3, correlated (vignette data)", predictors = predictors, outcome = outcome),
  local({ set.seed(11); x <- matrix(rnorm(80), 80, 1)
          list(label = "k = 1, n = 80", predictors = x, outcome = drop(0.3 * x) + rnorm(80)) }),
  local({ set.seed(12); x <- matrix(rnorm(100 * 5), 100) %*% chol(0.3 + 0.7 * diag(5))
          list(label = "k = 5, n = 100, correlated", predictors = x,
               outcome = drop(x %*% c(0.2, 0, 0.1, 0, 0)) + rnorm(100)) }))

set.seed(1)
for (case in brute_force_cases) {
  closed_form <- log_bayes_factor_flat_r_vs_intercept_only(
    sample_r_squared(case$predictors, case$outcome), length(case$outcome), ncol(case$predictors))
  brute_force <- log_bayes_factor_flat_r_by_sampling(case$predictors, case$outcome,
                                                      number_of_draws = 400000)
  record_check(paste("A3. brute force,", case$label),
               abs(closed_form - brute_force) < 0.05,
               sprintf("closed %.3f, brute force %.3f", closed_form, brute_force))
}


## =====================================================================
## A4. The same prior as a mixture of g-priors
## =====================================================================
## The flat-R prior is a mixture of Zellner g-priors with
## g / n = B * W, B ~ Beta(1/2, (k-1)/2), W ~ Inverse-Gamma(1, 1/2)
## (docs/2_how_it_calculates.md, section 8). Its Bayes factor is then the
## standard g-prior Bayes factor (Liang et al., 2008) averaged over that
## hyperprior: a different route to the same number.
## Criterion: within 0.002 on the log scale (nested numerical integration).

## Density of t = g / n: p(t) = E[B exp(-B / (2t))] / (2 t^2)
hyperprior_density_g_over_n <- function(t, k) {
  if (k == 1) return(exp(-1 / (2 * t)) / (2 * t^2))
  weighted_beta <- function(b) b * dbeta(b, 0.5, (k - 1) / 2) * exp(-b / (2 * t))
  split <- min(1, 60 * t)   # the integrand sits near b = 0 when t is small
  pieces <- integrate(weighted_beta, 0, split, rel.tol = 1e-10)$value +
    (if (split < 1) integrate(weighted_beta, split, 1, rel.tol = 1e-10)$value else 0)
  pieces / (2 * t^2)
}
log_bayes_factor_by_g_mixture <- function(r_squared, sample_size, k) {
  integrand <- function(log_t) sapply(log_t, function(lt) {
    t <- exp(lt); g <- sample_size * t
    exp((sample_size - 1 - k) / 2 * log1p(g) - (sample_size - 1) / 2 * log1p(g * (1 - r_squared))) *
      hyperprior_density_g_over_n(t, k) * t          # times t: the step is in log t
  })
  log(integrate(integrand, -20, 25, rel.tol = 1e-8, subdivisions = 2000)$value)
}
g_mixture_cases <- list(c(0.19, 60, 3), c(0.05, 100, 5), c(0.30, 200, 10), c(0.10, 30, 2), c(0.40, 50, 1))
g_mixture_errors <- sapply(g_mixture_cases, function(case)
  abs(log_bayes_factor_by_g_mixture(case[1], case[2], case[3]) -
      log_bayes_factor_flat_r_vs_intercept_only(case[1], case[2], case[3])))
record_check("A4. closed form vs mixture of g-priors (5 settings)", max(g_mixture_errors) < 0.002,
             sprintf("largest difference %.1e", max(g_mixture_errors)))


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

##    The exact boundary (docs/2_how_it_calculates.md, section 11): information
##    consistent if and only if n >= k + 3. At n = k + 2 the evidence stays
##    bounded as R2 -> 1; at n = k + 3 it keeps growing.
boundary <- sapply(1:3, function(k) {
  at <- function(n) sapply(c(1 - 1e-6, 1 - 1e-8), log_bayes_factor_flat_r_vs_intercept_only,
                           sample_size = n, number_of_predictors = k)
  c(bounded_at_k_plus_2 = diff(at(k + 2)) < 0.01, growing_at_k_plus_3 = diff(at(k + 3)) > 0.2)
})
record_check("B2b. information consistency starts exactly at n = k + 3 (k = 1-3)", all(boundary),
             "bounded at n = k + 2, growing at n = k + 3")

## B3. Consistency under the null, by simulation: x3 truly does nothing.
##    The median log BF01 for x3 should grow with n at a rate of about
##    1/2 log n (one extra null predictor, q = 1; B3b separates the cases
##    where q matters). Criterion: increasing in n, and slope against log n
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


## B3b. The rate depends on which smaller model is true (docs/3_consistency.md).
##     Against the intercept-only model the evidence for the null grows like
##     1/2 log n whatever the number q of extra null predictors: flat on R puts
##     prior mass of order n^(-1/2) near beta = 0, for every k. Between two
##     non-null models (y = 0.5 x1 + e; x1 against x1 plus q nulls) it grows
##     like q/2 log n, the usual penalty for q extra dimensions.
##     Criterion: slope of the median log BF01 against log n within 20% of
##     q/2 for every q (smaller model x1), and between 0.35 and 0.65 for every q
##     (smaller model intercept only).
set.seed(6)
extra_null_predictors <- c(1, 3, 6)
nested_design <- expand.grid(extra_nulls = extra_null_predictors,
                             smaller_model = c("x1", "intercept only"),
                             stringsAsFactors = FALSE)
nested_medians <- parallel_map(seq_len(nrow(nested_design)), function(row) {
  q <- nested_design$extra_nulls[row]
  has_effect <- nested_design$smaller_model[row] == "x1"
  k_small <- as.integer(has_effect)
  sapply(sample_sizes, function(n) median(replicate(200, {
    x <- matrix(rnorm(n * (k_small + q)), n)
    y <- (if (has_effect) 0.5 * x[, 1] else 0) + rnorm(n)
    log_bf_small <- if (has_effect)
      log_bayes_factor_flat_r_vs_intercept_only(sample_r_squared(x[, 1, drop = FALSE], y), n, 1) else 0
    log_bf_small - log_bayes_factor_flat_r_vs_intercept_only(sample_r_squared(x, y), n, k_small + q)
  })))
})
nested_design$slope <- sapply(nested_medians, function(m) coef(lm(m ~ log(sample_sizes)))[[2]])
from_x1   <- nested_design$smaller_model == "x1"
from_null <- !from_x1
record_check("B3b. nested, smaller model x1: slope q/2 (q = 1, 3, 6)",
             all(abs(nested_design$slope[from_x1] / (nested_design$extra_nulls[from_x1] / 2) - 1) < 0.2),
             sprintf("slopes %s", paste(sprintf("%.2f", nested_design$slope[from_x1]), collapse = " ")))
record_check("B3b. nested, smaller model intercept only: slope 1/2 (q = 1, 3, 6)",
             all(nested_design$slope[from_null] > 0.35 & nested_design$slope[from_null] < 0.65),
             sprintf("slopes %s", paste(sprintf("%.2f", nested_design$slope[from_null]), collapse = " ")))


## B4. Model selection consistency when there IS an effect: x1 = 0.15 is real,
##     so the median log BF10 for x1 should grow without limit with n.
##     Criterion: increasing in n, and above 20 at n = 5000.
set.seed(5)
median_log_bf10 <- sapply(c(200, 1000, 5000), function(n) {
  median(replicate(100, {
    x <- matrix(rnorm(n * 3), n)
    y <- drop(x %*% c(0.15, 0.2, 0)) + rnorm(n)
    log_bayes_factor_flat_r_vs_intercept_only(sample_r_squared(x, y), n, 3) -
      log_bayes_factor_flat_r_vs_intercept_only(sample_r_squared(x[, 2:3], y), n, 2)
  }))
})
record_check("B4. consistency under the alternative (simulated)",
             all(diff(median_log_bf10) > 0) && median_log_bf10[3] > 20,
             sprintf("median log BF10 at n = 200, 1000, 5000: %s",
                     paste(sprintf("%.1f", median_log_bf10), collapse = " ")))


## =====================================================================
## Summary
## =====================================================================
report_checks()
