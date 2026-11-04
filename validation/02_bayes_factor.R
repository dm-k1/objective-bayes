## =====================================================================
##  Validation 02: is the Bayes factor right?
##  One script, one pass/fail summary, the same checks for each of the three
##  priors (flat_r, flat_partial_r, flat_r_squared). Every check has an
##  explicit criterion; the script stops with an error if anything failed.
##
##  PART A. Correctness
##    A1. the grid computes the formula: against adaptive integration
##        (integrate()) over a sweep of R-squared, n and k, extremes included,
##        by a route that is not the package's
##    A2. invariance: rescaling or recoding predictors leaves the BF unchanged
##    A3. the formula is the right one: the package vs brute force (draw
##        coefficients from the prior, average how well they predict; no formula)
##    A4. a second derivation: flat_r as a mixture of g-priors; for every prior
##        the block route (the k columns tested against none,
##        R/conditional.R), which is how flat_r_squared is computed, against
##        the closed forms
##
##  PART B. Properties of a good default (Bayarri, Berger, Forte &
##  Garcia-Donato, 2012)
##    B1. predictive matching (n = k + 1 gives BF = 1)
##    B2. information consistency (BF grows without limit as R2 -> 1)
##    B2b. ... exactly from n = k + 3 on (the analytic boundary)
##    B3. consistency under the null (log BF01 grows with n), simulated
##    B3b. the rate: against the intercept-only model it is a log n, with a the
##         shape of R-squared ~ Beta(a, 1) (1/2, k/2 or 1: the prior's mass
##         near zero, P(R^2 < 1/n) ~ n^(-a)); with a real effect in the
##         smaller model it is q/2 log n (q extra nulls), the usual penalty
##    B4. consistency under the alternative (log BF10 grows with n), simulated
##    B5. intrinsic consistency: the prior uses the sample covariance of the
##        predictors, so it changes with n; it must converge to a proper prior
##        (the one with the population covariance) as n grows
##
##  Run from the project folder: Rscript validation/02_bayes_factor.R
##  A few minutes. No Stan needed.
## =====================================================================
pkgload::load_all(".", quiet = TRUE)   # the flatcor package, from this folder
source("validation/settings.R")             # parallel cores, cache folders, record_check()
source("validation/comparison_methods.R")   # log_bayes_factor_beta_r_squared() (A4)

priors <- c("flat_r", "flat_partial_r", "flat_r_squared")
log_bf <- function(prior, r_squared, n, k) log_bayes_factor_vs_intercept_only(r_squared, n, k, prior)


## =====================================================================
## Ground truth: the same Bayes factor by adaptive integration
## =====================================================================
## R/bayes_factor.R uses fixed grids (fast). Here each prior's Bayes factor is
## computed with integrate(): adaptive quadrature that adds points where the
## curve changes quickly and reports its own error estimate. Slow, but
## accurate, so it serves as the judge in A1. Everything is on the log scale,
## scaled by the largest value, so that astronomically large Bayes factors do
## not overflow. A different route per prior, none of them the package's:
##   flat_r          over g and over directions (the one-predictor prior as a
##                   mixture of normals, averaged over every direction)
##   flat_partial_r  over g: a g-prior mixture with g ~ Inverse-Gamma(1, n/2)
##   flat_r_squared  k >= 3: a g-prior mixture, g/n = B W, B ~ Beta(1, k/2 - 1),
##                   W ~ Inverse-Gamma(1, 1/2) (docs/2_how_it_calculates.md,
##                   section 8); k = 2: the same prior as flat_partial_r
##                   (a = k/2 = 1); k = 1: over the effect f directly

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



## ---- flat_partial_r: one integral over g ---------------------------------------
log_bf_g_inverse_gamma_by_integration <- function(r_squared, n, k) {
  log_integrand <- function(log_g) log(n / 2) - log_g - n / (2 * exp(log_g)) +
    (n - 1 - k) / 2 * log1p(exp(log_g)) - (n - 1) / 2 * log1p(exp(log_g) * (1 - r_squared))
  peak <- optimize(log_integrand, c(-20, 40), maximum = TRUE)$maximum
  height <- log_integrand(peak)
  break_points <- c(-40, peak - 3, peak + 3, 60)
  height + log(sum(vapply(1:3, function(i)
    integrate(function(x) exp(log_integrand(x) - height), break_points[i], break_points[i + 1],
              rel.tol = 1e-10, subdivisions = 1000L)$value, numeric(1))))
}

## ---- R-squared ~ Beta(a, 1) as a g-prior mixture, a < k/2 -----------------------
## g/n = t = B W, B ~ Beta(a, k/2 - a), W ~ Inverse-Gamma(1, 1/2), so
## p(t) = E[B exp(-B / (2t))] / (2 t^2). a = 1/2 is flat_r (A4), a = 1 is
## flat_r_squared for k >= 3.
log_hyperprior_density_g_over_n <- function(t, k, a) {
  weighted_beta <- function(b) b * dbeta(b, a, k / 2 - a) * exp(-b / (2 * t))
  split <- min(1, 60 * t)   # the integrand sits near b = 0 when t is small
  pieces <- integrate(weighted_beta, 0, split, rel.tol = 1e-10)$value +
    (if (split < 1) integrate(weighted_beta, split, 1, rel.tol = 1e-10)$value else 0)
  log(pieces) - log(2) - 2 * log(t)
}
log_bf_g_mixture_by_integration <- function(r_squared, n, k, a) {
  log_integrand <- function(log_t) vapply(log_t, function(lt) {
    g <- n * exp(lt)
    value <- (n - 1 - k) / 2 * log1p(g) - (n - 1) / 2 * log1p(g * (1 - r_squared)) +
      log_hyperprior_density_g_over_n(exp(lt), k, a) + lt      # + lt: the step is in log t
    if (is.finite(value)) value else -Inf   # far out in the tails the density underflows
  }, numeric(1))
  peak <- optimize(log_integrand, c(-20, 25), maximum = TRUE)$maximum
  height <- log_integrand(peak)
  break_points <- sort(c(-20, max(-19, peak - 3), min(24, peak + 3), 25))
  height + log(sum(vapply(1:3, function(i)
    integrate(function(x) exp(log_integrand(x) - height), break_points[i], break_points[i + 1],
              rel.tol = 1e-8, subdivisions = 2000L)$value, numeric(1))))
}

## ---- flat_r_squared with one predictor: over the effect f ---------------------
## |f| has density 2 f / (1 + f^2)^2 (r^2 = f^2 / (1 + f^2) flat). With the
## intercept flat and the noise 1/sigma integrated out,
##   p(y | f) / p(y | 0) = exp(-n f^2/2) (I(b) + I(-b)) / (2 I(0)),
##   I(b) = integral_0^inf v^m exp(-v^2/2 + b v) dv,  b = sqrt(n) r f,  m = n - 2,
## each I computed by integrate() around its peak.
## (At n = 10,000 a tolerance of 1e-10 can hit round-off; then 1e-8.)
log_I_by_integration <- function(b, m) {
  peak     <- (b + sqrt(b^2 + 4 * m)) / 2
  log_peak <- m * log(peak) - peak^2 / 2 + b * peak
  scaled   <- function(v) exp(m * log(v) - v^2 / 2 + b * v - log_peak)
  value <- tryCatch(integrate(scaled, max(0, peak - 40), peak + 40, rel.tol = 1e-10)$value,
                    error = function(e) integrate(scaled, max(0, peak - 40), peak + 40,
                                                  rel.tol = 1e-8)$value)
  log_peak + log(value)
}
log_bf_one_predictor_flat_r_squared_by_integration <- function(r_squared, n) {
  m <- n - 2; r <- sqrt(r_squared)
  log_I0 <- log_I_by_integration(0, m)
  log_integrand <- function(f) vapply(f, function(f1) {
    b <- sqrt(n) * r * f1
    terms <- c(log_I_by_integration(b, m), log_I_by_integration(-b, m))
    -n * f1^2 / 2 + max(terms) + log(sum(exp(terms - max(terms)))) - log(2) - log_I0 +
      log(2 * f1) - 2 * log1p(f1^2)
  }, numeric(1))
  peak <- optimize(log_integrand, c(1e-8, 50), maximum = TRUE)$maximum
  height <- log_integrand(peak)
  width <- 1 / sqrt(n)
  break_points <- sort(unique(c(1e-12, max(1e-12, peak - 30 * width), peak, peak + 30 * width, 200)))
  height + log(sum(vapply(seq_len(length(break_points) - 1), function(i)
    integrate(function(x) exp(log_integrand(x) - height), break_points[i], break_points[i + 1],
              rel.tol = 1e-8, subdivisions = 2000L)$value, numeric(1))))
}

## ---- Each prior's ground truth --------------------------------------------------
log_bf_by_integration <- function(prior, r_squared, n, k) {
  if (k == 0) return(0)
  switch(prior,
    flat_r         = log_bayes_factor_flat_r_by_integration(r_squared, n, k),
    flat_partial_r = log_bf_g_inverse_gamma_by_integration(r_squared, n, k),
    flat_r_squared = if (k == 1) log_bf_one_predictor_flat_r_squared_by_integration(r_squared, n) else
                     if (k == 2) log_bf_g_inverse_gamma_by_integration(r_squared, n, k) else
                     log_bf_g_mixture_by_integration(r_squared, n, k, a = 1))
}

## ---- Example data (same as the vignette) ---------------------
set.seed(42)
sample_size <- 60
predictor_correlations <- matrix(0.4, 3, 3); diag(predictor_correlations) <- 1
predictors <- matrix(rnorm(sample_size * 3), sample_size) %*% chol(predictor_correlations)
colnames(predictors) <- c("x1", "x2", "x3")
outcome <- drop(predictors %*% c(0.3, 0.2, 0)) + rnorm(sample_size)
example_data <- data.frame(predictors, y = outcome)


## The same log Bayes factor as log_bayes_factor_vs_intercept_only(), without
## the formula: draw coefficient sets from the prior, see how well each
## predicts the outcome (the noise level averaged on a grid), and average.
## Slow and noisy (Monte Carlo); agrees to about 0.05 on the log scale.
log_bayes_factor_by_sampling <- function(predictors, outcome, prior, number_of_draws = 40000,
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
    ## flat_partial_r and flat_r_squared put little mass at small R, where the
    ## data usually are, so plain Monte Carlo would rest on a few rare draws.
    ## Instead: draw from the flat_r recipe (R flat) and weight each draw by
    ## the ratio of the two densities of R, 2a R^(2a - 1) (importance sampling;
    ## the weights are bounded by 2a). For flat_r every weight is 1.
    coefficients     <- draw_coefficients_from_prior(draws_per_batch, predictor_covariance, prior = "flat_r")
    f_squared_drawn  <- rowSums((coefficients %*% predictor_covariance) * coefficients)
    a <- prior_shape(prior, ncol(predictors))
    log_weight <- log(2 * a) + (2 * a - 1) / 2 * log(f_squared_drawn / (1 + f_squared_drawn))
    predicted_values <- centered_predictors %*% t(coefficients)
    overlap_with_outcome     <- drop(crossprod(centered_outcome, predicted_values))
    predicted_sum_of_squares <- colSums(predicted_values^2)
    log_integrand <- (sample_size - 2) *
        matrix(log_inverse_noise, draws_per_batch, length(inverse_noise_values), byrow = TRUE) -
      (outer(rep(outcome_sum_of_squares, draws_per_batch), inverse_noise_values^2) -
       2 * outer(overlap_with_outcome, inverse_noise_values) + predicted_sum_of_squares) / 2
    log_evidence_per_draw <- c(log_evidence_per_draw,
                               log_sum_exp_per_row(log_integrand) + log(inverse_noise_step) + log_weight)
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

## The brute-force cases of A3
brute_force_cases <- list(
  list(label = "k = 3, correlated (vignette data)", predictors = predictors, outcome = outcome),
  local({ set.seed(11); x <- matrix(rnorm(80), 80, 1)
          list(label = "k = 1, n = 80", predictors = x, outcome = drop(0.3 * x) + rnorm(80)) }),
  local({ set.seed(12); x <- matrix(rnorm(100 * 5), 100) %*% chol(0.3 + 0.7 * diag(5))
          list(label = "k = 5, n = 100, correlated", predictors = x,
               outcome = drop(x %*% c(0.2, 0, 0.1, 0, 0)) + rnorm(100)) }))

## Settings of A1, B3b
sweep_settings <- subset(expand.grid(explained_share      = c(0.01, 0.1, 0.3, 0.6, 0.9, 0.99),
                                     sample_size          = c(10, 30, 100, 1000, 10000),
                                     number_of_predictors = c(1, 2, 3, 5, 10, 30)),
                         sample_size > number_of_predictors + 1)
sample_sizes <- c(50, 200, 1000, 5000)
nested_design <- expand.grid(extra_nulls = c(1, 3, 6), smaller_model = c("x1", "intercept only"),
                             stringsAsFactors = FALSE)
from_x1 <- nested_design$smaller_model == "x1"


for (prior in priors) {
  cat("\n## ========================= ", prior, " =========================\n\n", sep = "")

  ## ###################################################################
  ##  PART A. CORRECTNESS
  ## ###################################################################

  ## A1. Against the ground truth, over a sweep of inputs. Criterion: within
  ##     0.001 for every setting. The worst settings are listed, so you can see
  ##     where the grid is not accurate enough.
  message(prior, ": ground truth for ", nrow(sweep_settings), " settings (adaptive integration) ...")
  sweep_values <- parallel_map(seq_len(nrow(sweep_settings)), function(row) {
    s <- sweep_settings[row, ]
    c(code         = log_bf(prior, s$explained_share, s$sample_size, s$number_of_predictors),
      ground_truth = log_bf_by_integration(prior, s$explained_share, s$sample_size,
                                           s$number_of_predictors))
  })
  sweep <- cbind(sweep_settings, do.call(rbind, sweep_values))
  sweep$error <- abs(sweep$code - sweep$ground_truth)
  ordinary <- abs(sweep$ground_truth) < 50     # |log BF| < 50: the range that matters in practice
  record_check(sprintf("A1a. %s: ground truth, ordinary evidence (%d settings)", prior, sum(ordinary)),
               max(sweep$error[ordinary]) < 1e-3,
               sprintf("largest error %.1e", max(sweep$error[ordinary])))
  record_check(sprintf("A1b. %s: ground truth, all settings (%d)", prior, nrow(sweep)),
               max(sweep$error) < 1e-3, sprintf("largest error %.1e", max(sweep$error)))
  cat("\nSettings with the largest errors:\n")
  print(head(sweep[order(-sweep$error), ], 6), row.names = FALSE, digits = 6)
  cat("\n")

  ## A2. Invariance: the BF should not depend on units or on how the same
  ##     predictors are written. Criterion: identical up to 1e-8 (relative).
  bf_original <- bayes_factor_per_term(y ~ x1 + x2 + x3, example_data, prior = prior)$BF10
  rescaled_data <- transform(example_data, x1 = 100 * x1 + 7, x2 = -0.01 * x2, x3 = x3 - 50,
                             y = 1000 * y - 3)
  bf_rescaled <- bayes_factor_per_term(y ~ x1 + x2 + x3, rescaled_data, prior = prior)$BF10
  recoded_data <- transform(example_data, sum_12 = x1 + x2, difference_12 = x1 - x2)
  bf_full_original <- bayes_factor_models(lm(y ~ x1 + x2 + x3, example_data),
                                          lm(y ~ 1, example_data), prior = prior)
  bf_full_recoded  <- bayes_factor_models(lm(y ~ sum_12 + difference_12 + x3, recoded_data),
                                          lm(y ~ 1, recoded_data), prior = prior)
  largest <- max(abs(bf_rescaled / bf_original - 1), abs(bf_full_recoded / bf_full_original - 1))
  record_check(sprintf("A2. %s: rescaling, shifting and recoding", prior), largest < 1e-8,
               sprintf("largest relative difference %.1e", largest))

  ## A3. Package vs brute force: draw coefficient sets from the prior (the
  ##     recipe of docs/2 section 5), see how well each predicts, average.
  ##     Monte Carlo, so noisy: with 400,000 draws its SD is about 0.02 on the
  ##     log scale (k = 5). Criterion: within 0.05 on the log scale.
  set.seed(1)
  for (case in brute_force_cases) {
    package_value <- log_bf(prior, sample_r_squared(case$predictors, case$outcome),
                            length(case$outcome), ncol(case$predictors))
    brute_force <- log_bayes_factor_by_sampling(case$predictors, case$outcome, prior,
                                                number_of_draws = 400000)
    record_check(sprintf("A3. %s: brute force, %s", prior, case$label),
                 abs(package_value - brute_force) < 0.05,
                 sprintf("package %.3f, brute force %.3f", package_value, brute_force))
  }

  ## A4. A second derivation. flat_r: the prior as a mixture of Zellner
  ##     g-priors, g/n = B W, B ~ Beta(1/2, (k-1)/2), W ~ Inverse-Gamma(1, 1/2)
  ##     (docs/2_how_it_calculates.md, section 8), by nested integrate().
  ##     flat_r and flat_partial_r: the block route, which computes
  ##     flat_r_squared, against their closed forms. flat_r_squared: the
  ##     the block route against its closed form (a g-mixture with the
  ##     hyperprior of t = g/n computed once per k), and for k >= 3 the
  ##     g-mixture by nested integrate().
  ##     Criterion: within 0.002 on the log scale.
  ##     (The nested integral for flat_r runs on the first five cases, those of
  ##     the original check; at very small R-squared its inner integral over a
  ##     Beta(1/2, .) density is unstable. The block route covers all seven,
  ##     extreme evidence included.)
  g_mixture_cases <- list(c(0.19, 60, 3), c(0.05, 100, 5), c(0.30, 200, 10), c(0.10, 30, 2),
                          c(0.40, 50, 1), c(0.9, 30, 4), c(0.001, 500, 3), c(0.99, 10000, 30))
  second_route <- function(case) {
    r_squared <- case[1]; n <- case[2]; k <- case[3]
    original_case <- any(vapply(g_mixture_cases[1:5], identical, logical(1), case))
    switch(prior,
      flat_r = c(g_mixture = if (!original_case) NULL else
                             if (k == 1) log_bf_g_inverse_gamma_by_integration(r_squared, n, 1) else
                             log_bf_g_mixture_by_integration(r_squared, n, k, a = 1 / 2),
                 block = log_bayes_factor_shape_vs_intercept_only(r_squared, n, k, 1 / 2)),
      flat_partial_r = c(block = log_bayes_factor_shape_vs_intercept_only(r_squared, n, k, k / 2)),
      flat_r_squared = c(block = log_bayes_factor_shape_vs_intercept_only(r_squared, n, k, 1),
                         if (k >= 3 && r_squared < 0.99)
                           c(g_mixture = log_bf_g_mixture_by_integration(r_squared, n, k, a = 1))))
  }
  a4_differences <- unlist(lapply(g_mixture_cases, function(case)
    abs(second_route(case) - log_bf(prior, case[1], case[2], case[3]))))
  record_check(sprintf("A4. %s: second derivation (%d comparisons)", prior, length(a4_differences)),
               max(a4_differences) < 0.002, sprintf("largest difference %.1e", max(a4_differences)))


  ## ###################################################################
  ##  PART B. PROPERTIES OF A GOOD DEFAULT
  ## ###################################################################

  ## B1. Predictive matching: with one more observation than predictors the
  ##    data can't tell the models apart, so BF = 1. Criterion: within 0.001.
  ##    (flat_r_squared with k = 1 goes through the block route, which returns
  ##    exactly 1 at n = k + 1 by this property, so k = 1 is not a test there.)
  matching <- sapply(if (prior == "flat_r_squared") 2:6 else 1:6, function(k)
    exp(log_bf(prior, 1 - 1e-12, k + 1, k)))
  record_check(sprintf("B1. %s: predictive matching (k = %s)", prior,
                       if (prior == "flat_r_squared") "2-6" else "1-6"),
               all(abs(matching - 1) < 1e-3), paste(sprintf("%.4f", matching), collapse = " "))

  ## B2. Information consistency: as the fit becomes perfect, evidence for an
  ##    effect grows without limit. Criterion: strictly increasing, and above
  ##    log BF = 50 at R2 = 0.9999.
  information <- sapply(c(0.9, 0.99, 0.999, 0.9999), function(r2) log_bf(prior, r2, 20, 3))
  record_check(sprintf("B2. %s: information consistency", prior),
               all(diff(information) > 0) && information[4] > 50,
               paste(sprintf("%.1f", information), collapse = " "))

  ##    The exact boundary (docs/2_how_it_calculates.md, section 11): information
  ##    consistent if and only if n >= k + 3. At n = k + 2 the evidence stays
  ##    bounded as R2 -> 1; at n = k + 3 it keeps growing. All three priors
  ##    have the same tail, P(f > t) ~ a / t^2.
  boundary <- sapply(1:3, function(k) {
    at <- function(n) sapply(c(1 - 1e-6, 1 - 1e-8), function(r2) log_bf(prior, r2, n, k))
    c(bounded_at_k_plus_2 = diff(at(k + 2)) < 0.01, growing_at_k_plus_3 = diff(at(k + 3)) > 0.2)
  })
  record_check(sprintf("B2b. %s: information consistency starts exactly at n = k + 3 (k = 1-3)", prior),
               all(boundary), "bounded at n = k + 2, growing at n = k + 3")

  ## B3. Consistency under the null, by simulation: x3 truly does nothing.
  ##    The median log BF01 for x3 should grow with n (one extra null
  ##    predictor next to two real ones, so about 1/2 log n). Criterion:
  ##    increasing in n, and slope against log n between 0.35 and 0.65.
  set.seed(4)
  median_log_bf01 <- sapply(sample_sizes, function(n) median(replicate(200, {
    x <- matrix(rnorm(n * 3), n)
    y <- drop(x %*% c(0.3, 0.2, 0)) + rnorm(n)
    -(log_bf(prior, sample_r_squared(x, y), n, 3) - log_bf(prior, sample_r_squared(x[, 1:2], y), n, 2))
  })))
  slope_against_log_n <- coef(lm(median_log_bf01 ~ log(sample_sizes)))[[2]]
  record_check(sprintf("B3. %s: consistency under the null (simulated)", prior),
               all(diff(median_log_bf01) > 0) && slope_against_log_n > 0.35 && slope_against_log_n < 0.65,
               sprintf("median log BF01 %s; slope %.2f",
                       paste(sprintf("%.2f", median_log_bf01), collapse = " "), slope_against_log_n))

  ## B3b. The rate depends on which smaller model is true (docs/3_consistency.md).
  ##     Smaller model x1 (y = 0.5 x1 + e; x1 against x1 plus q nulls): q/2 log n.
  ##     Smaller model intercept only: a log n, the prior's mass near zero
  ##     (flat_r 1/2 and flat_r_squared 1 for every q; flat_partial_r q/2).
  ##     Criterion: slope within 20% of q/2 (smaller model x1) and within 30%
  ##     of a (intercept only), for every q. (The 30% band is that of the
  ##     original flat_r check, 0.35-0.65 around 1/2.)
  set.seed(6)
  nested_medians <- parallel_map(seq_len(nrow(nested_design)), function(row) {
    q <- nested_design$extra_nulls[row]
    has_effect <- nested_design$smaller_model[row] == "x1"
    k_small <- as.integer(has_effect)
    sapply(sample_sizes, function(n) median(replicate(200, {
      x <- matrix(rnorm(n * (k_small + q)), n)
      y <- (if (has_effect) 0.5 * x[, 1] else 0) + rnorm(n)
      log_bf_small <- if (has_effect) log_bf(prior, sample_r_squared(x[, 1, drop = FALSE], y), n, 1) else 0
      log_bf_small - log_bf(prior, sample_r_squared(x, y), n, k_small + q)
    })))
  })
  slopes <- sapply(nested_medians, function(m) coef(lm(m ~ log(sample_sizes)))[[2]])
  expected_slope <- ifelse(from_x1, nested_design$extra_nulls / 2,
                           sapply(nested_design$extra_nulls, function(q) prior_shape(prior, q)))
  relative_error <- abs(slopes / expected_slope - 1)
  record_check(sprintf("B3b. %s: smaller model x1, slope q/2 (q = 1, 3, 6)", prior),
               all(relative_error[from_x1] < 0.2),
               sprintf("slopes %s (expected %s)", paste(sprintf("%.2f", slopes[from_x1]), collapse = " "),
                       paste(sprintf("%.2f", expected_slope[from_x1]), collapse = " ")))
  record_check(sprintf("B3b. %s: smaller model intercept only, slope a (q = 1, 3, 6)", prior),
               all(relative_error[!from_x1] < 0.3),
               sprintf("slopes %s (expected %s)", paste(sprintf("%.2f", slopes[!from_x1]), collapse = " "),
                       paste(sprintf("%.2f", expected_slope[!from_x1]), collapse = " ")))

  ## B4. Model selection consistency when there IS an effect: x1 = 0.15 is
  ##     real, so the median log BF10 for x1 should grow without limit with n.
  ##     Criterion: increasing in n, and above 20 at n = 5000.
  set.seed(5)
  median_log_bf10 <- sapply(c(200, 1000, 5000), function(n) median(replicate(100, {
    x <- matrix(rnorm(n * 3), n)
    y <- drop(x %*% c(0.15, 0.2, 0)) + rnorm(n)
    log_bf(prior, sample_r_squared(x, y), n, 3) - log_bf(prior, sample_r_squared(x[, 2:3], y), n, 2)
  })))
  record_check(sprintf("B4. %s: consistency under the alternative (simulated)", prior),
               all(diff(median_log_bf10) > 0) && median_log_bf10[3] > 20,
               sprintf("median log BF10 at n = 200, 1000, 5000: %s",
                       paste(sprintf("%.1f", median_log_bf10), collapse = " ")))

  ## B5. Intrinsic consistency (Bayarri et al., 2012, criterion 5): the prior
  ##     on the coefficients is defined through the sample covariance Sigma_n of
  ##     the predictors, so it depends on n. As n grows it must tend to a proper
  ##     prior, the one with the population covariance. Predictors drawn from
  ##     three correlated normals; the implied prior of the first coefficient
  ##     (noise units) with Sigma_n against the one with the population
  ##     covariance. Criterion: below 0.01 at n = 100,000 (100,000 draws each:
  ##     sampling noise about 0.005), and at n = 30 at least twice that.
  ##     (Set after the first run, which asked for a distance shrinking at every
  ##     step: by n = 300 it is already at the noise floor, so n = 300 and
  ##     100,000 cannot be ordered.)
  set.seed(8)
  population <- matrix(0.5, 3, 3); diag(population) <- 1
  limit_draws <- draw_coefficients_from_prior(100000, population, prior = prior)[, 1]
  intrinsic <- sapply(c(30, 300, 100000), function(n) {
    x <- matrix(rnorm(n * 3), n) %*% chol(population)
    sample_covariance <- crossprod(scale(x, scale = FALSE)) / n
    suppressWarnings(unname(ks.test(draw_coefficients_from_prior(100000, sample_covariance,
                                                                 prior = prior)[, 1],
                                    limit_draws)$statistic))
  })
  record_check(sprintf("B5. %s: intrinsic consistency (prior with Sigma_n -> prior with Sigma)", prior),
               intrinsic[3] < 0.01 && intrinsic[1] > 2 * intrinsic[3],
               sprintf("KS distance at n = 30, 300, 100000: %s",
                       paste(sprintf("%.4f", intrinsic), collapse = " ")))
}


## =====================================================================
## Summary
## =====================================================================
report_checks()
