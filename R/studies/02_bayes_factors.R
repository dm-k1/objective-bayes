## =====================================================================
##  Study 2: are the Bayes factors computed correctly?
##  (docs/appendix/C_study_design.md, study 2; replaces part A of the old study 02 and
##  the old study 05)
##
##  Aims. Every Bayes factor flatcor reports equals the integral it stands
##  for: the whole-model Bayes factor (R-squared ~ Beta(a, 1) with a uniform
##  direction; flat_r a = 1/2, flat_partial_r a = k/2, flat_r_squared a = 1)
##  and the conditional test of one coefficient or a block of q columns.
##  Data. Fixed grids of summary statistics (R-squared, n, k; the partial
##  correlation and the number of columns tested), extremes included, and a
##  few simulated datasets for the routes that need data.
##  Estimands. The log Bayes factor of each route.
##  Methods. flatcor's closed forms and grids, against independent routes that
##  share no code with the package:
##    adaptive integration (stats::integrate()), as in the old studies;
##    a high-precision ground truth with Rmpfr (200-bit arithmetic): a
##      trapezoid rule on log g, which converges exponentially for these smooth
##      integrands, checked by halving the step;
##    brute force (draw coefficients from the prior, average how well they
##      predict; no formula);
##    bridge sampling (bridgesampling) on draws from the package's Stan model.
##  Performance. The absolute log difference from the reference. For the
##  Monte Carlo routes (brute force, bridge sampling) the difference is judged
##  against their own Monte Carlo error, and each runs until that error is
##  below 0.01 on the log scale (a stated precision for an estimator, not a
##  stopping rule on the result).
##  n_sim. Deterministic, except: the Monte Carlo routes (to MCSE 0.01); and
##  the distribution checks (2.8, 2.10, 2.13, 2.15), which use
##  precision_targets$ks_draws draws (KS critical distance at the .01 level
##  of .01), judged against the KS critical value for their number of draws,
##  with the .01 level split over the checks of a family (Bonferroni).
##
##  Checks (old check in brackets):
##    Whole model, for each prior
##    2.1 the grid against integrate(), ordinary evidence (|log BF| < 50) and
##        all settings: within 1e-3 [A1a, A1b]
##    2.2 the closed forms against a 200-bit ground truth (flat_partial_r for
##        every k; flat_r_squared at k = 2, where it is the same prior), on a
##        grid with extremes (R-squared to .99999, n to 1e6, k to 30), as 2.1
##        is split: ordinary evidence (|log BF| < 50) within 1e-6 [2.2a]; all
##        settings within a relative error of 1e-6 [2.2b]. integrate()'s own
##        error against the same truth is reported (new). The split was made
##        after the first full run (5 November 2026), in which one absolute
##        tolerance of 1e-6 failed only where the evidence is astronomical
##        (log BF up to 5.8 million, off by up to 3.1: a relative error of
##        5e-7), while every setting with |log BF| < 1,000 was within 1e-6.
##    2.3 invariance under rescaling, shifting and recoding: within 1e-8 [A2]
##    2.4 brute force, run to MCSE 0.01: within 3 MCSE [A3, which used a fixed
##        400,000 draws and a fixed 0.05]
##    2.5 a second derivation (g-prior mixtures; the block route): within
##        0.002 [A4]
##    2.6 bridge sampling on the package's Stan fit, run to an estimated error
##        of 0.01: within 3 of its errors (new). bridgesampling's error allows
##        for autocorrelation (the spectral density of the posterior draws;
##        Fruhwirth-Schnatter, 2004), but not for a biased sampler, so the
##        fits use adapt_delta = 0.99. In the first full run (flatcor 1.0.0),
##        flat_r_squared at k = 1 had 13 divergent transitions and missed by
##        0.0047 against 3 x 0.0012: R had the density |R| on (-1, 1), which
##        vanishes at 0. flatcor 1.0.1 samples |R| and marginalises its sign
##    One coefficient and blocks
##    2.7 k = 1: the conditional test equals the joint test [C1]
##    2.8 the implied prior for one coefficient does not depend on how the
##        predictors correlate (KS) [C2]
##    2.9 the exact flat partial correlation against its approximation [C3]
##    2.10 flat_partial_r: one partial correlation is flat, any k (KS) [C4a]
##    2.11 the share mixture with a = k/2 is that flat partial correlation [C4b]
##    2.12 the conditional grid against integrate() [C5]
##    2.13 LKJ(eta): one partial correlation is Beta(eta, eta); flat only at
##         eta = 1 (KS) [C7]
##    2.14 against Kucharský et al. (2026), with their series in 200-bit
##         arithmetic as the reference [C8]
##    2.15 the implied prior of a block of q columns (KS) [C9]
##    2.16 the block test against integrate() [C10]
##    2.17 the block likelihood ratio from the data themselves [C11]
##
##  Run from the project folder:  Rscript R/studies/02_bayes_factors.R
##  Smoke test (small grids, few draws, a few minutes):
##    VALIDATION_SMOKE=1 Rscript R/studies/02_bayes_factors.R
##  Needs cmdstanr with CmdStan (check 2.6), Rmpfr and bridgesampling.
##  Raw results are cached in cache/study_02/ (delete them to rerun); the
##  tables and checks are saved to results/studies/study_02.rds.
## =====================================================================
library(dplyr)
source("R/shared/settings.R")             # flatcor, cores, record_check(), parallel_map()
source("R/shared/comparison_methods.R")   # draw_lkj_correlation_matrix(), partial_correlation_of()
source("R/shared/precision.R")            # precision_targets, smoke_mode

study_folder <- file.path(cache_folder, if (smoke_mode) "study_02_smoke" else "study_02")
dir.create(study_folder, recursive = TRUE, showWarnings = FALSE)

priors <- c("flat_r", "flat_partial_r", "flat_r_squared")
log_bf <- function(prior, r_squared, n, k) log_bayes_factor_vs_intercept_only(r_squared, n, k, prior)

## draws for the distribution checks (a smoke run uses 2,000: enough to run)
ks_draws <- if (smoke_mode) 2000L else precision_targets$ks_draws$n_sim

## Kolmogorov-Smirnov critical distance at level alpha (asymptotic), one-sample
## with n draws or two-sample with n and m draws:
## sqrt(-log(alpha / 2) / 2) * sqrt((n + m) / (n m))
ks_critical <- function(alpha, n, m = Inf) sqrt(-log(alpha / 2) / 2) * sqrt(1 / n + 1 / m)
ks_distance <- function(...) suppressWarnings(unname(ks.test(...)$statistic))

## a cached computation: the expression is evaluated once and its value saved
cached <- function(name, expression) {
  file <- file.path(study_folder, paste0(name, ".rds"))
  if (file.exists(file)) return(readRDS(file))
  value <- expression
  saveRDS(value, file)
  value
}


## =====================================================================
## Ground truths by adaptive integration (from the old study 02)
## =====================================================================
## A different route per prior, none of them the package's:
##   flat_r          over g and over directions (the one-predictor prior as a
##                   mixture of normals, averaged over every direction)
##   flat_partial_r  over g: a g-prior mixture with g ~ Inverse-Gamma(1, n/2)
##   flat_r_squared  k >= 3: a g-prior mixture, g/n = B W, B ~ Beta(1, k/2 - 1),
##                   W ~ Inverse-Gamma(1, 1/2) (docs/appendix/A_derivations.md,
##                   section 8); k = 2: the same prior as flat_partial_r
##                   (a = k/2 = 1); k = 1: over the effect f directly
## Everything on the log scale, scaled by the largest value.

## average over directions for one value of g: the angle phi between a
## direction and the data's direction, weight sin(phi)^(k - 2)
log_direction_average_by_integration <- function(scaled_share, sample_size, number_of_predictors) {
  half_n_minus_one <- (sample_size - 1) / 2
  log_value_at_zero <- -half_n_minus_one * log1p(-scaled_share)
  if (number_of_predictors == 1) return(log_value_at_zero)
  log_total_weight <- 0.5 * log(pi) - log(2) +
    lgamma((number_of_predictors - 1) / 2) - lgamma(number_of_predictors / 2)
  scaled_integrand <- function(angle)
    exp(-half_n_minus_one * log1p(-scaled_share * cos(angle)^2) - log_value_at_zero +
        (number_of_predictors - 2) * log(sin(angle)) - log_total_weight)
  ## strong evidence makes a narrow peak near angle 0: split the interval there
  peak_width <- if (scaled_share > 0)
    sqrt((1 - scaled_share) / (2 * half_n_minus_one * scaled_share)) else Inf
  inner_breaks <- c(3, 10, 40) * peak_width
  break_points <- c(0, inner_breaks[inner_breaks < pi / 2], pi / 2)
  pieces <- vapply(seq_len(length(break_points) - 1), function(i)
    integrate(scaled_integrand, break_points[i], break_points[i + 1],
              rel.tol = 1e-10, subdivisions = 1000L)$value, numeric(1))
  log_value_at_zero + log(sum(pieces))
}

log_bayes_factor_flat_r_by_integration <- function(explained_share, sample_size, number_of_predictors) {
  if (number_of_predictors == 0) return(0)
  ## the integrand over log(g): IG(1, n/2) density times g, the g-prior factor,
  ## and the average over directions
  log_integrand <- function(log_g_values) vapply(log_g_values, function(log_g) {
    g <- exp(log_g)
    log(sample_size / 2) - 2 * log_g - sample_size / (2 * g) + log_g - 0.5 * log1p(g) +
      log_direction_average_by_integration(explained_share * g / (1 + g), sample_size,
                                           number_of_predictors)
  }, numeric(1))
  peak_location <- optimize(log_integrand, c(-20, 40), maximum = TRUE)$maximum
  peak_height   <- log_integrand(peak_location)
  scaled <- function(log_g) exp(log_integrand(log_g) - peak_height)
  break_points <- c(-40, peak_location - 3, peak_location + 3, 60)
  pieces <- vapply(1:3, function(i)
    integrate(scaled, break_points[i], break_points[i + 1], rel.tol = 1e-10, subdivisions = 1000L)$value,
    numeric(1))
  peak_height + log(sum(pieces))
}

## flat_partial_r: one integral over g, with g ~ Inverse-Gamma(1, n/2)
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

## R-squared ~ Beta(a, 1) as a g-prior mixture, a < k/2: g/n = t = B W,
## B ~ Beta(a, k/2 - a), W ~ Inverse-Gamma(1, 1/2), p(t) = E[B exp(-B/(2t))] / (2 t^2)
log_hyperprior_density_g_over_n <- function(t, k, a) {
  weighted_beta <- function(b) b * dbeta(b, a, k / 2 - a) * exp(-b / (2 * t))
  split <- min(1, 60 * t)
  pieces <- integrate(weighted_beta, 0, split, rel.tol = 1e-10)$value +
    (if (split < 1) integrate(weighted_beta, split, 1, rel.tol = 1e-10)$value else 0)
  log(pieces) - log(2) - 2 * log(t)
}
log_bf_g_mixture_by_integration <- function(r_squared, n, k, a) {
  log_integrand <- function(log_t) vapply(log_t, function(lt) {
    g <- n * exp(lt)
    value <- (n - 1 - k) / 2 * log1p(g) - (n - 1) / 2 * log1p(g * (1 - r_squared)) +
      log_hyperprior_density_g_over_n(exp(lt), k, a) + lt
    if (is.finite(value)) value else -Inf
  }, numeric(1))
  peak <- optimize(log_integrand, c(-20, 25), maximum = TRUE)$maximum
  height <- log_integrand(peak)
  break_points <- sort(c(-20, max(-19, peak - 3), min(24, peak + 3), 25))
  height + log(sum(vapply(1:3, function(i)
    integrate(function(x) exp(log_integrand(x) - height), break_points[i], break_points[i + 1],
              rel.tol = 1e-8, subdivisions = 2000L)$value, numeric(1))))
}

## log of the integral of v^m exp(-v^2/2 + b v) over v > 0, around its peak
log_I_by_integration <- function(b, m) {
  peak     <- (b + sqrt(b^2 + 4 * m)) / 2
  log_peak <- m * log(peak) - peak^2 / 2 + b * peak
  scaled   <- function(v) exp(m * log(v) - v^2 / 2 + b * v - log_peak)
  value <- tryCatch(integrate(scaled, max(0, peak - 40), peak + 40, rel.tol = 1e-10)$value,
                    error = function(e) integrate(scaled, max(0, peak - 40), peak + 40,
                                                  rel.tol = 1e-8)$value)
  log_peak + log(value)
}

## flat_r_squared with one predictor, over the effect f: |f| has density
## 2 f / (1 + f^2)^2, and p(y | f) / p(y | 0) = exp(-n f^2/2) (I(b) + I(-b)) / (2 I(0))
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

log_bf_by_integration <- function(prior, r_squared, n, k) {
  if (k == 0) return(0)
  switch(prior,
    flat_r         = log_bayes_factor_flat_r_by_integration(r_squared, n, k),
    flat_partial_r = log_bf_g_inverse_gamma_by_integration(r_squared, n, k),
    flat_r_squared = if (k == 1) log_bf_one_predictor_flat_r_squared_by_integration(r_squared, n) else
                     if (k == 2) log_bf_g_inverse_gamma_by_integration(r_squared, n, k) else
                     log_bf_g_mixture_by_integration(r_squared, n, k, a = 1))
}


## =====================================================================
## A 200-bit ground truth with Rmpfr (new)
## =====================================================================
## The flat_partial_r Bayes factor as one integral over u = log g (the same
## integrand as log_bf_g_inverse_gamma_by_integration(), but evaluated in
## 200-bit arithmetic):
##   BF = integral of (n/2) exp(-u - n e^-u / 2) (1 + e^u)^((n-1-k)/2)
##                    (1 + e^u (1 - R^2))^(-(n-1)/2) du.
## The integrand is smooth and decays fast on both sides, so the trapezoid
## rule on a wide interval converges exponentially in the step; the value is
## accepted when halving the step changes it by less than 1e-12 on the log
## scale. The peak height (in double precision) only scales the sum.
log_bf_flat_partial_r_high_precision <- function(r_squared, n, k, bits = 200, step = 0.01) {
  log_integrand_double <- function(u) log(n / 2) - u - n / 2 * exp(-u) +
    (n - 1 - k) / 2 * log1p(exp(u)) - (n - 1) / 2 * log1p(exp(u) * (1 - r_squared))
  peak <- optimize(log_integrand_double, c(-30, 60), maximum = TRUE)$maximum
  height <- log_integrand_double(peak)
  one <- Rmpfr::mpfr(1, bits)
  n_m <- n * one; k_m <- k * one; share_m <- Rmpfr::mpfr(r_squared, bits)
  trapezoid <- function(h) {
    u <- Rmpfr::mpfr(seq(peak - 60, peak + 60, by = h), bits)
    log_values <- log(n_m / 2) - u - n_m / 2 * exp(-u) +
      (n_m - 1 - k_m) / 2 * log1p(exp(u)) - (n_m - 1) / 2 * log1p(exp(u) * (1 - share_m)) - height
    values <- exp(log_values)
    total <- sum(values) - (values[1] + values[length(values)]) / 2
    height + as.numeric(log(total * h))
  }
  coarse <- trapezoid(step); fine <- trapezoid(step / 2)
  c(value = fine, step_change = abs(fine - coarse))
}


## =====================================================================
## Brute force, run until its Monte Carlo error is below 0.01 (from A3)
## =====================================================================
## The log Bayes factor without the formula: draw coefficient sets from the
## prior, see how well each predicts the outcome (the noise level averaged on
## a grid of 1 / noise), and average. Importance sampling from the flat_r
## recipe with weights 2a R^(2a - 1) for the other priors (bounded by 2a).
## Returns the per-draw log evidence, so the Monte Carlo error can be
## estimated: for log(mean(w)), MCSE = sd(w) / (mean(w) sqrt(N)).
log_evidence_per_draw_by_sampling <- function(predictors, outcome, prior, number_of_draws,
                                              draws_per_batch = 2000) {
  sample_size            <- length(outcome)
  centered_predictors    <- scale(predictors, scale = FALSE)
  centered_outcome       <- outcome - mean(outcome)
  outcome_sum_of_squares <- sum(centered_outcome^2)
  predictor_covariance   <- crossprod(centered_predictors) / sample_size
  inverse_noise_values <- seq(1e-4, 6 * sqrt(sample_size / outcome_sum_of_squares), length.out = 800)
  inverse_noise_step   <- inverse_noise_values[2] - inverse_noise_values[1]
  log_inverse_noise    <- log(inverse_noise_values)
  log_evidence_null <- log_sum_exp((sample_size - 2) * log_inverse_noise -
                                     inverse_noise_values^2 * outcome_sum_of_squares / 2) +
                       log(inverse_noise_step)
  a <- prior_shape(prior, ncol(predictors))
  unlist(lapply(seq_len(ceiling(number_of_draws / draws_per_batch)), function(batch) {
    coefficients    <- draw_coefficients_from_prior(draws_per_batch, predictor_covariance, prior = "flat_r")
    f_squared_drawn <- rowSums((coefficients %*% predictor_covariance) * coefficients)
    log_weight <- log(2 * a) + (2 * a - 1) / 2 * log(f_squared_drawn / (1 + f_squared_drawn))
    predicted_values <- centered_predictors %*% t(coefficients)
    overlap_with_outcome     <- drop(crossprod(centered_outcome, predicted_values))
    predicted_sum_of_squares <- colSums(predicted_values^2)
    log_integrand <- (sample_size - 2) *
        matrix(log_inverse_noise, draws_per_batch, length(inverse_noise_values), byrow = TRUE) -
      (outer(rep(outcome_sum_of_squares, draws_per_batch), inverse_noise_values^2) -
       2 * outer(overlap_with_outcome, inverse_noise_values) + predicted_sum_of_squares) / 2
    log_sum_exp_per_row(log_integrand) + log(inverse_noise_step) + log_weight - log_evidence_null
  }))
}
estimate_from_log_evidence <- function(log_evidence) {
  largest <- max(log_evidence)
  weights <- exp(log_evidence - largest)
  c(log_bf = largest + log(mean(weights)), mcse = sd(weights) / (mean(weights) * sqrt(length(weights))))
}
## draws in rounds of 100,000 (in parallel) until the MCSE is below 0.01, at most
## 6.4 million; a smoke run takes one round of 20,000
log_bayes_factor_by_sampling_to_precision <- function(predictors, outcome, prior, target = 0.01) {
  round_size <- if (smoke_mode) 20000 else 100000
  maximum <- if (smoke_mode) round_size else 6.4e6
  log_evidence <- numeric(0)
  repeat {
    log_evidence <- c(log_evidence, unlist(parallel_map(seq_len(number_of_cores), function(core)
      log_evidence_per_draw_by_sampling(predictors, outcome, prior,
                                        ceiling(round_size / number_of_cores)))))
    estimate <- estimate_from_log_evidence(log_evidence)
    if (estimate[["mcse"]] < target || length(log_evidence) >= maximum) break
  }
  c(estimate, draws = length(log_evidence))
}


## =====================================================================
## Bridge sampling on the package's Stan fit (new)
## =====================================================================
## The joint-model Bayes factor against the intercept-only model as a ratio
## of marginal likelihoods. The full model's posterior draws come from the
## package's Stan model (system.file("stan", "flatcor.stan", package =
## "flatcor"), fitted with cmdstanr). bridgesampling's matrix method then
## needs the log posterior with every constant, which is written here in R
## (the Stan program drops constants):
##   likelihood    normal, with all constants
##   alpha, sigma  flat and 1/sigma (the package's improper priors; the same in
##                 both models, so their arbitrary scale cancels in the ratio)
##   R_magnitude   the size of R, density 2a R^(2a - 1) on (0, 1); with k = 1 the
##                 sign of R is + or - with probability 1/2, and the Stan model
##                 (flatcor 1.0.1) marginalises it: the likelihood is averaged
##                 over both signs, here as in the Stan program
##   z             k > 1: the auxiliary vector of the direction, density
##                 exp(-((|z| - 1)/0.1)^2 / 2) / C_k on R^k, with C_k = surface
##                 area of the unit sphere times the integral over the radius
##   coefficients  b = sigma f U^-1 (z/|z|), f = R / sqrt(1 - R^2), U the upper
##                 Cholesky factor of the predictors' sample covariance (as in
##                 the Stan model)
## The intercept-only model's marginal likelihood is in closed form:
##   n^(-1/2) (2 pi)^(-(n-1)/2) Gamma((n-1)/2) / 2 * (SST/2)^(-(n-1)/2).
## Its error: bridgesampling's coefficient of variation of the marginal
## likelihood, which is the MCSE of its log. Draws are doubled until it is
## below 0.01 (at most 4 x 32,000).
log_normalizing_constant_of_direction <- function(k) {
  radial <- integrate(function(r) exp(-((r - 1) / 0.1)^2 / 2 + (k - 1) * log(r)), 0, Inf,
                      rel.tol = 1e-12)$value
  log(2) + k / 2 * log(pi) - lgamma(k / 2) + log(radial)
}
log_marginal_likelihood_intercept_only <- function(outcome) {
  n <- length(outcome); sst <- sum((outcome - mean(outcome))^2)
  -0.5 * log(n) - (n - 1) / 2 * log(2 * pi) + lgamma((n - 1) / 2) - log(2) - (n - 1) / 2 * log(sst / 2)
}
bridge_sampling_log_bf <- function(predictors, outcome, prior, stan_model, seed, target = 0.01) {
  n <- length(outcome); k <- ncol(predictors); a <- prior_shape(prior, k)
  centered <- scale(predictors, scale = FALSE)
  upper_root <- chol(crossprod(centered) / n)
  log_c <- if (k > 1) log_normalizing_constant_of_direction(k) else 0
  parameter_names <- c("alpha", "sigma", "R_magnitude", if (k > 1) paste0("z[", seq_len(k), "]"))
  log_posterior <- function(pars, data) {
    sigma <- pars[["sigma"]]; R <- pars[["R_magnitude"]]
    direction <- if (k > 1) { z <- pars[paste0("z[", seq_len(k), "]")]; z / sqrt(sum(z^2)) } else 1
    b <- sigma * R / sqrt(1 - R^2) * backsolve(upper_root, direction)
    log_prior <- -log(sigma) + log(2 * a) + (2 * a - 1) * log(R) +
      (if (k > 1) -((sqrt(sum(z^2)) - 1) / 0.1)^2 / 2 - log_c else 0)
    log_likelihood <- function(sign) sum(dnorm(outcome, pars[["alpha"]] + sign * drop(centered %*% b), sigma, log = TRUE))
    if (k > 1) return(log_likelihood(1) + log_prior)
    log_sum_exp(c(log_likelihood(1), log_likelihood(-1))) - log(2) + log_prior
  }
  lower <- setNames(c(-Inf, 0, 0, rep(-Inf, if (k > 1) k else 0)), parameter_names)
  upper <- setNames(c(Inf, Inf, 1, rep(Inf, if (k > 1) k else 0)), parameter_names)
  draws_per_chain <- if (smoke_mode) 1000 else 4000
  repeat {
    fit <- stan_model$sample(data = list(N = n, K = k, X = predictors, y = outcome, prior_only = 0L,
                                         proper_stand_ins = 0L, r_squared_shape = a),
                             chains = 4, parallel_chains = min(4, number_of_cores),
                             iter_warmup = 1000, iter_sampling = draws_per_chain,
                             adapt_delta = 0.99, seed = seed, refresh = 0, show_messages = FALSE)
    samples <- posterior::as_draws_matrix(fit$draws(parameter_names))
    samples <- unclass(samples)[, parameter_names, drop = FALSE]
    bridge <- suppressWarnings(bridgesampling::bridge_sampler(samples = samples,
                log_posterior = log_posterior, data = NULL, lb = lower, ub = upper, silent = TRUE))
    error <- bridgesampling::error_measures(bridge)$cv
    if (error < target || smoke_mode || draws_per_chain >= 32000) break
    draws_per_chain <- 2 * draws_per_chain
  }
  c(log_bf = bridge$logml - log_marginal_likelihood_intercept_only(outcome), error = error,
    draws = 4 * draws_per_chain, divergent = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent))
}


## =====================================================================
## Data: the example (as docs) and the brute-force cases (from A3)
## =====================================================================
set.seed(42)
example_n <- 60
example_correlations <- matrix(0.4, 3, 3); diag(example_correlations) <- 1
example_predictors <- matrix(rnorm(example_n * 3), example_n) %*% chol(example_correlations)
colnames(example_predictors) <- c("x1", "x2", "x3")
example_outcome <- drop(example_predictors %*% c(0.3, 0.2, 0)) + rnorm(example_n)
example_data <- data.frame(example_predictors, y = example_outcome)

data_cases <- list(
  list(label = "k = 3, correlated (example data)", predictors = example_predictors,
       outcome = example_outcome),
  local({ set.seed(11); x <- matrix(rnorm(80), 80, 1)
          list(label = "k = 1, n = 80", predictors = x, outcome = drop(0.3 * x) + rnorm(80)) }),
  local({ set.seed(12); x <- matrix(rnorm(100 * 5), 100) %*% chol(0.3 + 0.7 * diag(5))
          list(label = "k = 5, n = 100, correlated", predictors = x,
               outcome = drop(x %*% c(0.2, 0, 0.1, 0, 0)) + rnorm(100)) }))
if (smoke_mode) data_cases <- data_cases[1]


## =====================================================================
## WHOLE MODEL
## =====================================================================
sweep_settings <- subset(expand.grid(explained_share      = c(0.01, 0.1, 0.3, 0.6, 0.9, 0.99),
                                     sample_size          = c(10, 30, 100, 1000, 10000),
                                     number_of_predictors = c(1, 2, 3, 5, 10, 30)),
                         sample_size > number_of_predictors + 1)
if (smoke_mode) sweep_settings <- subset(sweep_settings, sample_size %in% c(30, 1000) &
                                           number_of_predictors %in% c(1, 3))

## 2.2's grid: extremes included
high_precision_settings <- subset(expand.grid(explained_share = c(0.001, 0.1, 0.5, 0.9, 0.999, 0.99999),
                                              sample_size = c(20, 200, 1e4, 1e6),
                                              number_of_predictors = c(2, 3, 10, 30)),
                                  sample_size > number_of_predictors + 2)
if (smoke_mode) high_precision_settings <- subset(high_precision_settings,
                                                  sample_size == 200 & number_of_predictors == 3)

g_mixture_cases <- list(c(0.19, 60, 3), c(0.05, 100, 5), c(0.30, 200, 10), c(0.10, 30, 2),
                        c(0.40, 50, 1), c(0.9, 30, 4), c(0.001, 500, 3), c(0.99, 10000, 30))
if (smoke_mode) g_mixture_cases <- g_mixture_cases[c(1, 5)]

stan_model <- local({
  ## compiled from a copy of the package's Stan file in the cache, so the
  ## executable is not written into the installed package
  stan_file <- file.path(study_folder, "flatcor.stan")
  file.copy(system.file("stan", "flatcor.stan", package = "flatcor"), stan_file, overwrite = TRUE)
  cmdstanr::cmdstan_model(stan_file, quiet = TRUE)
})

for (prior in priors) {
  cat("\n## ========================= ", prior, " =========================\n\n", sep = "")

  ## 2.1 [A1] the grid against adaptive integration
  sweep <- cached(paste0("sweep_", prior), {
    message(prior, ": integrate() for ", nrow(sweep_settings), " settings ...")
    values <- parallel_map(seq_len(nrow(sweep_settings)), function(row) {
      s <- sweep_settings[row, ]
      c(code         = log_bf(prior, s$explained_share, s$sample_size, s$number_of_predictors),
        ground_truth = log_bf_by_integration(prior, s$explained_share, s$sample_size,
                                             s$number_of_predictors))
    })
    cbind(sweep_settings, do.call(rbind, values))
  })
  sweep$error <- abs(sweep$code - sweep$ground_truth)
  ordinary <- abs(sweep$ground_truth) < 50
  record_check(sprintf("2.1a %s: grid = integrate(), ordinary evidence (%d settings)", prior, sum(ordinary)),
               max(sweep$error[ordinary]) < 1e-3, sprintf("largest error %.1e", max(sweep$error[ordinary])))
  record_check(sprintf("2.1b %s: grid = integrate(), all settings (%d)", prior, nrow(sweep)),
               max(sweep$error) < 1e-3, sprintf("largest error %.1e", max(sweep$error)))

  ## 2.2 (new) the 200-bit ground truth, where it applies: flat_partial_r (every
  ## k) and flat_r_squared at k = 2 (the same prior, a = k/2 = 1)
  if (prior != "flat_r") {
    settings <- if (prior == "flat_partial_r") high_precision_settings else
      transform(unique(high_precision_settings[, c("explained_share", "sample_size")]),
                number_of_predictors = 2)
    high <- cached(paste0("high_precision_", prior), {
      message(prior, ": 200-bit ground truth for ", nrow(settings), " settings ...")
      values <- parallel_map(seq_len(nrow(settings)), function(row) {
        s <- settings[row, ]
        truth <- log_bf_flat_partial_r_high_precision(s$explained_share, s$sample_size,
                                                      s$number_of_predictors)
        c(truth, code = log_bf(prior, s$explained_share, s$sample_size, s$number_of_predictors),
          integrate = log_bf_g_inverse_gamma_by_integration(s$explained_share, s$sample_size,
                                                            s$number_of_predictors))
      })
      cbind(settings, do.call(rbind, values))
    })
    ordinary <- abs(high$value) < 50
    absolute_error <- abs(high$code - high$value)
    relative_error <- absolute_error / pmax(1, abs(high$value))
    label <- if (prior == "flat_r_squared") ", k = 2" else ""
    record_check(sprintf("2.2a %s: closed form = 200-bit ground truth, ordinary evidence (%d settings%s)",
                         prior, sum(ordinary), label),
                 max(absolute_error[ordinary]) < 1e-6 && max(high$step_change) < 1e-12,
                 sprintf("largest error %.1e; integrate() against the same truth: largest %.1e; step check %.0e",
                         max(absolute_error[ordinary]), max(abs(high$integrate - high$value)[ordinary]),
                         max(high$step_change)))
    record_check(sprintf("2.2b %s: closed form = 200-bit ground truth, all settings, relative (%d settings%s)",
                         prior, nrow(high), label),
                 max(relative_error) < 1e-6,
                 sprintf("largest relative error %.1e (absolute %.1e at log BF %.3g)", max(relative_error),
                         absolute_error[which.max(relative_error)], high$value[which.max(relative_error)]))
  }

  ## 2.3 [A2] invariance
  bf_original <- bayes_factor_per_term(y ~ x1 + x2 + x3, example_data, prior = prior)$BF10
  rescaled_data <- transform(example_data, x1 = 100 * x1 + 7, x2 = -0.01 * x2, x3 = x3 - 50,
                             y = 1000 * y - 3)
  bf_rescaled <- bayes_factor_per_term(y ~ x1 + x2 + x3, rescaled_data, prior = prior)$BF10
  recoded_data <- transform(example_data, sum_12 = x1 + x2, difference_12 = x1 - x2)
  bf_full_original <- bayes_factor_models(lm(y ~ x1 + x2 + x3, example_data), lm(y ~ 1, example_data),
                                          prior = prior)
  bf_full_recoded  <- bayes_factor_models(lm(y ~ sum_12 + difference_12 + x3, recoded_data),
                                          lm(y ~ 1, recoded_data), prior = prior)
  largest <- max(abs(bf_rescaled / bf_original - 1), abs(bf_full_recoded / bf_full_original - 1))
  record_check(sprintf("2.3 %s: rescaling, shifting and recoding", prior), largest < 1e-8,
               sprintf("largest relative difference %.1e", largest))

  ## 2.4 [A3] brute force, to MCSE 0.01
  brute <- cached(paste0("brute_force_", prior), {
    RNGkind("L'Ecuyer-CMRG"); set.seed(21)
    lapply(data_cases, function(case) {
      message(prior, ": brute force, ", case$label, " ...")
      log_bayes_factor_by_sampling_to_precision(case$predictors, case$outcome, prior)
    })
  })
  for (i in seq_along(data_cases)) {
    case <- data_cases[[i]]; estimate <- brute[[i]]
    package_value <- log_bf(prior, sample_r_squared(case$predictors, case$outcome),
                            length(case$outcome), ncol(case$predictors))
    record_check(sprintf("2.4 %s: brute force, %s", prior, case$label),
                 abs(package_value - estimate[["log_bf"]]) < 3 * estimate[["mcse"]] &&
                   (smoke_mode || estimate[["mcse"]] < 0.01),
                 sprintf("package %.4f, brute force %.4f (MCSE %.4f, %s draws)", package_value,
                         estimate[["log_bf"]], estimate[["mcse"]], format(estimate[["draws"]], big.mark = ",")))
  }

  ## 2.5 [A4] a second derivation: flat_r as a g-prior mixture by nested
  ## integrate() (on the first five cases; at very small R-squared the inner
  ## integral is unstable), and for every prior the block route against the
  ## closed form; flat_r_squared also the g-mixture for k >= 3
  second_route <- function(case) {
    r_squared <- case[1]; n <- case[2]; k <- case[3]
    original_case <- any(vapply(g_mixture_cases[seq_len(min(5, length(g_mixture_cases)))],
                                identical, logical(1), case))
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
  a5_differences <- unlist(lapply(g_mixture_cases, function(case)
    abs(second_route(case) - log_bf(prior, case[1], case[2], case[3]))))
  record_check(sprintf("2.5 %s: second derivation (%d comparisons)", prior, length(a5_differences)),
               max(a5_differences) < 0.002, sprintf("largest difference %.1e", max(a5_differences)))

  ## 2.6 (new) bridge sampling on the package's Stan fit, to an error of 0.01
  bridged <- cached(paste0("bridge_sampling_", prior), lapply(seq_along(data_cases), function(i) {
    message(prior, ": bridge sampling, ", data_cases[[i]]$label, " ...")
    bridge_sampling_log_bf(data_cases[[i]]$predictors, data_cases[[i]]$outcome, prior, stan_model,
                           seed = 600 + i)
  }))
  for (i in seq_along(data_cases)) {
    case <- data_cases[[i]]; estimate <- bridged[[i]]
    package_value <- log_bf(prior, sample_r_squared(case$predictors, case$outcome),
                            length(case$outcome), ncol(case$predictors))
    record_check(sprintf("2.6 %s: bridge sampling on the Stan fit, %s", prior, case$label),
                 abs(package_value - estimate[["log_bf"]]) < 3 * estimate[["error"]] &&
                   (smoke_mode || estimate[["error"]] < 0.01),
                 sprintf("package %.4f, bridge sampling %.4f (error %.4f, %s draws, %d divergent)",
                         package_value, estimate[["log_bf"]], estimate[["error"]],
                         format(estimate[["draws"]], big.mark = ","), as.integer(estimate[["divergent"]])))
  }
}


## =====================================================================
## ONE COEFFICIENT AND BLOCKS (from the old study 05)
## =====================================================================
cat("\n## ========================= one coefficient and blocks =========================\n\n")

## 2.7 [C1] with one predictor the conditional test is the joint test
c1_differences <- unlist(lapply(priors, function(prior)
  sapply(list(c(0.3, 30), c(0.1, 100), c(0.5, 50), c(0.05, 500)), function(s)
    abs(log_bayes_factor_conditional(s[1], s[2], 1, prior) -
        log_bayes_factor_vs_intercept_only(s[1]^2, s[2], 1, prior)))))
record_check("2.7 k = 1: conditional = joint Bayes factor (every prior)", max(c1_differences) < 1e-3,
             sprintf("largest difference %.1e", max(c1_differences)))

## 2.8 [C2] the implied prior for one coefficient's partial correlation does
## not depend on how the predictors correlate: correlated joint draws against
## the share formula (two-sample KS, two tests at a family level of .01)
RNGkind("L'Ecuyer-CMRG"); set.seed(1)
c2_draws <- if (smoke_mode) 20000 else 200000
correlated <- matrix(0.6, 4, 4); diag(correlated) <- 1
spread_left <- sqrt(1 / diag(solve(correlated)))[2]
c2_distances <- sapply(c(flat_r = 1 / 2, flat_r_squared = 1), function(a) {
  joint_draws <- draw_coefficients_from_prior(c2_draws, correlated,
                                              prior = if (a == 1) "flat_r_squared" else "flat_r")
  from_joint <- joint_draws[, 2] * spread_left / sqrt(1 + (joint_draws[, 2] * spread_left)^2)
  r_squared  <- rbeta(c2_draws, a, 1)
  f_formula  <- sqrt(rbeta(c2_draws, 0.5, 1.5) * r_squared / (1 - r_squared))
  from_formula <- sample(c(-1, 1), c2_draws, replace = TRUE) * f_formula / sqrt(1 + f_formula^2)
  ks_distance(from_joint, from_formula)
})
c2_critical <- ks_critical(0.01 / 2, c2_draws, c2_draws)
record_check("2.8 implied prior: correlated joint draws = share formula (k = 4; flat_r, flat_r_squared)",
             all(c2_distances < c2_critical),
             sprintf("KS distances %s; critical %.4f", paste(sprintf("%.4f", c2_distances), collapse = " "),
                     c2_critical))

## 2.9 [C3] the exact flat partial correlation (flat_partial_r conditional
## test) against the approximation that plugs in the estimated noise: close,
## not equal (within 0.05)
set.seed(3)
c3_differences <- sapply(c(60, 100, 200, 400), function(n) {
  X <- matrix(rnorm(n * 3), n) %*% chol(0.5 + 0.5 * diag(3))
  y <- drop(X %*% c(0.15, 0.3, 0)) + rnorm(n)
  fit <- lm(y ~ X)
  abs(log_bayes_factor_conditional(partial_correlation_of(fit, 2), n, 3, "flat_partial_r") -
      log(bayes_factor_flat_partial_correlation(fit, 1, X)))
})
record_check("2.9 exact vs approximate flat partial correlation", max(c3_differences) < 0.05,
             sprintf("largest difference %.3f", max(c3_differences)))

## 2.10 [C4a] flat_partial_r as the g-prior b ~ N(0, W Sigma^-1), W ~
## Inverse-Gamma(1, 1/2): one partial correlation is flat on (-1, 1) for any k
RNGkind("L'Ecuyer-CMRG"); set.seed(4)
c4_k <- c(2, 5, 12)
c4_distances <- sapply(c4_k, function(k) {
  covariance <- matrix(0.6, k, k); diag(covariance) <- 1
  precision  <- solve(covariance)
  W <- 1 / (2 * rexp(ks_draws))
  b <- sqrt(W) * MASS::mvrnorm(ks_draws, rep(0, k), precision)
  f <- b[, 1] / sqrt(precision[1, 1])
  ks_distance(f / sqrt(1 + f^2), "punif", -1, 1)
})
c4_critical <- ks_critical(0.01 / length(c4_k), ks_draws)
record_check("2.10 flat_partial_r: one partial correlation is flat (k = 2, 5, 12)",
             all(c4_distances < c4_critical),
             sprintf("KS distances %s; critical %.4f", paste(sprintf("%.4f", c4_distances), collapse = " "),
                     c4_critical))

## 2.11 [C4b] the share mixture with a = k/2 is the flat partial correlation,
## P(|f| < f) = f / sqrt(1 + f^2), for every k; for a block, partial R^2 ~
## Beta(q/2, 1). Within 1e-4.
f_values <- c(0.003, 0.03, 0.1, 0.3, 1, 3, 10)
c4b_differences <- c(
  unlist(lapply(c(2, 3, 5, 12, 30), function(k)
    abs(prior_cdf_of_f(f_values, k, 500, 1, shape = k / 2) - f_values / sqrt(1 + f_values^2)))),
  unlist(lapply(list(c(5, 2), c(8, 3)), function(kq)
    abs(prior_cdf_of_f(f_values, kq[1], 500, kq[2], shape = kq[1] / 2) -
          (f_values^2 / (1 + f_values^2))^(kq[2] / 2)))))
record_check("2.11 share mixture with a = k/2 = flat partial correlation (k = 2-30; blocks)",
             max(c4b_differences) < 1e-4, sprintf("largest difference %.1e", max(c4b_differences)))

## 2.12 [C5] the conditional test against integrate(): BF = integral of
## L(f) p(f) df, L(f) = exp(-n f^2/2) (I(b) + I(-b)) / (2 I(0)), b = sqrt(n) r f,
## p(f) the density of |f| = sqrt(s) R / sqrt(1 - R^2) mixed over the share
## s ~ Beta(1/2, (k-1)/2); P(|f| < f | s) = (f^2 / (s + f^2))^a. Within 2e-4.
prior_density_given_share <- function(f, s, a)
  a * (f^2 / (s + f^2))^(a - 1) * 2 * f * s / (s + f^2)^2
log_bayes_factor_conditional_by_integration <- function(partial_r, n, k, prior) {
  m <- n - k - 1
  likelihood_ratio <- function(f) {
    b <- sqrt(n) * partial_r * f
    log_ratio_by_sign <- -n * f^2 / 2 - log_I_by_integration(0, m) +
      c(log_I_by_integration(b, m), log_I_by_integration(-b, m))
    mean(exp(log_ratio_by_sign))
  }
  a <- prior_shape(prior, k)
  prior_density <- function(f) {
    if (prior == "flat_partial_r") return((1 + f^2)^(-3 / 2))
    if (k == 1) return(prior_density_given_share(f, 1, a))
    given_share <- function(s) stats::dbeta(s, 0.5, (k - 1) / 2) * prior_density_given_share(f, s, a)
    integrate(function(u) given_share(exp(u)) * exp(u), log(f^2) - 50, log(0.5),
              rel.tol = 1e-10)$value + integrate(given_share, 0.5, 1, rel.tol = 1e-10)$value
  }
  integrand <- Vectorize(function(f) likelihood_ratio(f) * prior_density(f))
  peak  <- sqrt(m / n) * abs(partial_r) / sqrt(1 - partial_r^2)
  width <- 1 / sqrt(n * (1 - partial_r^2))
  breaks <- unique(pmax(0, peak + c(-15, -5, -2, 0, 2, 5, 15) * width))
  log(sum(sapply(seq_len(length(breaks) - 1), function(i)
    integrate(integrand, breaks[i], breaks[i + 1], rel.tol = 1e-8)$value)))
}
c5_cases <- expand.grid(t_value = c(0, 1, 2, 4, 6), k = c(1, 2, 3, 5, 10, 30), n = c(50, 200, 1e4, 1e6),
                        prior = priors, stringsAsFactors = FALSE)
c5_cases <- subset(c5_cases, k > 1 | prior != "flat_partial_r")
if (smoke_mode) c5_cases <- subset(c5_cases, k %in% c(1, 3) & n %in% c(50, 1e4) & t_value %in% c(0, 2))
c5_cases <- cached("conditional_grid", {
  c5_cases$difference <- unlist(parallel_map(seq_len(nrow(c5_cases)), function(row) {
    case <- c5_cases[row, ]
    partial_r <- case$t_value / sqrt(case$t_value^2 + case$n - case$k - 1)
    abs(log_bayes_factor_conditional(partial_r, case$n, case$k, case$prior) -
        log_bayes_factor_conditional_by_integration(partial_r, case$n, case$k, case$prior))
  }))
  c5_cases
})
for (tested_prior in unique(c5_cases$prior)) {
  of_prior <- c5_cases[c5_cases$prior == tested_prior, ]
  worst <- of_prior[which.max(of_prior$difference), ]
  record_check(sprintf("2.12 %s: conditional grid = integrate() (n 50 to 1e6, k 1 to 30, t 0 to 6)",
                       tested_prior),
               max(of_prior$difference) < 2e-4,
               sprintf("largest difference %.1e (n = %g, k = %d, t = %g)", max(of_prior$difference),
                       worst$n, worst$k, worst$t_value))
}

## 2.13 [C7] why eta = 1: under LKJ(eta) the partial correlation of y with x1
## given all the others is Beta(eta, eta) on (-1, 1) for every k, so flat for
## every k only at eta = 1. Six tests against Beta(eta, eta) (family level
## .01), and at eta = 2 the distance to flat must exceed the critical value.
c7_design <- expand.grid(k = c(1, 3, 8), eta = c(1, 2))
c7_partials <- cached("lkj_partial_correlations", {
  RNGkind("L'Ecuyer-CMRG"); set.seed(8)
  parallel_map(seq_len(nrow(c7_design)), function(row) replicate(ks_draws, {
    precision <- solve(draw_lkj_correlation_matrix(c7_design$k[row] + 1, c7_design$eta[row]))
    -precision[1, 2] / sqrt(precision[1, 1] * precision[2, 2])
  }))
})
c7_design$to_beta <- mapply(function(partial, eta) ks_distance((partial + 1) / 2, "pbeta", eta, eta),
                            c7_partials, c7_design$eta)
c7_design$to_flat <- sapply(c7_partials, function(partial) ks_distance(partial, "punif", -1, 1))
c7_critical <- ks_critical(0.01 / nrow(c7_design), ks_draws)
record_check("2.13 LKJ(eta): one partial correlation ~ Beta(eta, eta) for k = 1, 3, 8; flat only at eta = 1",
             all(c7_design$to_beta < c7_critical) && all(c7_design$to_flat[c7_design$eta == 2] > c7_critical),
             sprintf("KS to Beta(eta, eta): %s; to flat at eta = 2: %s; critical %.4f",
                     paste(sprintf("%.3f", c7_design$to_beta), collapse = " "),
                     paste(sprintf("%.3f", c7_design$to_flat[c7_design$eta == 2]), collapse = " "),
                     c7_critical))

## 2.14 [C8] against the exact partial-correlation Bayes factor of Kucharský,
## Wagenmakers, van den Bergh & Ly (2026, TEST, eq. 2) with a flat prior
## (alpha = 1; gamma = delta = 0): with k - 1 controls and h = (n - k)/2,
##   BF = B(1/2, alpha + h) / B(1/2, alpha) 2F1(h, h; alpha + h + 1/2; r^2).
## Their predictors are random (the tested predictor's residual SD is a
## parameter); ours conditions on X. Their hypergeometric series is now
## summed in 200-bit arithmetic as the reference (it was double precision).
## (a) their formula reproduces their example (r = .51, n = 40, no controls:
##     BF10 = 45.88) within 0.01;
## (b) ours differs from theirs: the largest difference is reported;
## (c) ours with their prior on f (the flat partial correlation mixed over the
##     predictor's residual variance) equals theirs within 1e-5.
log_bayes_factor_kucharsky_high_precision <- function(partial_r, sample_size, number_of_controls,
                                                      alpha = 1, terms = 20000, bits = 200) {
  half_df <- (sample_size - number_of_controls - 1) / 2
  if (partial_r == 0) return(lbeta(0.5, alpha + half_df) - lbeta(0.5, alpha))
  j <- Rmpfr::mpfr(0:(terms - 1), bits)
  h <- Rmpfr::mpfr(half_df, bits); c_term <- Rmpfr::mpfr(alpha + half_df + 0.5, bits)
  log_terms <- 2 * (lgamma(h + j) - lgamma(h)) - lgamma(c_term + j) + lgamma(c_term) -
    lgamma(j + 1) + j * log(Rmpfr::mpfr(partial_r, bits)^2)
  largest <- max(log_terms)
  lbeta(0.5, alpha + half_df) - lbeta(0.5, alpha) + as.numeric(largest + log(sum(exp(log_terms - largest))))
}
prior_cdf_of_f_random_predictor <- function(f_values, sample_size, number_of_predictors) {
  df <- sample_size - number_of_predictors
  sapply(f_values, function(f) {
    if (f == 0) return(0)
    integrate(function(w) stats::dchisq(w, df) * f * sqrt(sample_size / w) / sqrt(1 + f^2 * sample_size / w),
              stats::qchisq(1e-14, df), stats::qchisq(1e-14, df, lower.tail = FALSE), rel.tol = 1e-12)$value
  })
}
c8_paper_example <- exp(log_bayes_factor_kucharsky_high_precision(0.51, 40, 0))
c8_cases <- expand.grid(partial_r = c(0, 0.1, 0.3, 0.5), n = c(20, 50, 200, 1000), k = c(1, 3, 8))
if (smoke_mode) c8_cases <- subset(c8_cases, n %in% c(50, 200) & k %in% c(1, 3))
c8_cases <- cached("kucharsky", {
  c8_results <- do.call(rbind, parallel_map(seq_len(nrow(c8_cases)), function(row) {
    partial_r <- c8_cases$partial_r[row]; n <- c8_cases$n[row]; k <- c8_cases$k[row]
    theirs <- log_bayes_factor_kucharsky_high_precision(partial_r, n, k - 1)
    c(ours_minus_theirs = log_bayes_factor_conditional(partial_r, n, k, "flat_partial_r") - theirs,
      random_predictor_minus_theirs = log_bayes_factor_from_prior_cdf_of_f(partial_r, n, k,
        function(f_values) prior_cdf_of_f_random_predictor(f_values, n, k), points_f = 2000) - theirs)
  }))
  cbind(c8_cases, c8_results)
})
worst <- c8_cases[which.max(abs(c8_cases$ours_minus_theirs)), ]
record_check(paste("2.14 Kucharský et al. (2026), flat prior: reproduces their BF 45.88; equals ours",
                   "with the predictor's residual SD unknown; differs from ours (fixed X)"),
             abs(c8_paper_example - 45.88) < 0.01 && max(abs(c8_cases$random_predictor_minus_theirs)) < 1e-5,
             sprintf(paste("their example %.3f; ours with random predictor - theirs: largest %.1e;",
                           "ours - theirs: largest %.3f (n = %g, k = %d, r = %g)"),
                     c8_paper_example, max(abs(c8_cases$random_predictor_minus_theirs)),
                     worst$ours_minus_theirs, worst$n, worst$k, worst$partial_r))
cat("2.14 ours - theirs, log BF (rows k and n, columns the partial correlation):\n")
c8_cases$k_and_n <- sprintf("k = %d, n = %4d", c8_cases$k, c8_cases$n)
print(round(xtabs(ours_minus_theirs ~ k_and_n + partial_r, c8_cases), 4))

## 2.15 [C9] the implied prior for a block of q columns:
## (a) flat_partial_r from LKJ(1) matrices: the block's partial R^2 is
##     Beta(q/2, 1) for every k;
## (b) flat_partial_r from g-prior draws with correlated predictors: the
##     block's f^2 / (1 + f^2) ~ Beta(q/2, 1);
## (c) flat_r and flat_r_squared: the block's share of the total is
##     Beta(q/2, (k - q)/2).
## Every test against the KS critical value at a family level of .01.
c9_lkj <- subset(expand.grid(k = c(2, 3, 5, 8), q = c(2, 3, 4)), q <= k)
c9_lkj$distance <- cached("block_lkj", {
  RNGkind("L'Ecuyer-CMRG"); set.seed(9)
  unlist(parallel_map(seq_len(nrow(c9_lkj)), function(row) {
    k <- c9_lkj$k[row]; q <- c9_lkj$q[row]
    others <- setdiff(seq_len(k + 1), 1 + seq_len(q))
    partial_r_squared <- replicate(ks_draws, {
      correlation <- draw_lkj_correlation_matrix(k + 1, 1)
      1 - solve(correlation[others, others])[1, 1] / solve(correlation)[1, 1]
    })
    ks_distance(partial_r_squared, "pbeta", q / 2, 1)
  }))
})
RNGkind("L'Ecuyer-CMRG"); set.seed(10)
c9_draws <- expand.grid(k = c(4, 6, 10), q = c(2, 3))
c9_draws <- cbind(c9_draws, do.call(rbind, lapply(seq_len(nrow(c9_draws)), function(row) {
  k <- c9_draws$k[row]; q <- c9_draws$q[row]; block <- seq_len(q)
  covariance <- matrix(0.6, k, k); diag(covariance) <- 1
  block_given_others <- covariance[block, block] - covariance[block, -block, drop = FALSE] %*%
    solve(covariance[-block, -block], covariance[-block, block, drop = FALSE])
  block_f_squared <- function(b) rowSums((b[, block, drop = FALSE] %*% block_given_others) *
                                           b[, block, drop = FALSE])
  partial <- draw_coefficients_from_prior(ks_draws, covariance, "flat_partial_r")
  share_distance <- function(draws)
    ks_distance(block_f_squared(draws) / rowSums((draws %*% covariance) * draws), "pbeta", q / 2, (k - q) / 2)
  c(flat_partial_r = ks_distance(block_f_squared(partial) / (1 + block_f_squared(partial)), "pbeta", q / 2, 1),
    flat_r_share = share_distance(draw_coefficients_from_prior(ks_draws, covariance, "flat_r")),
    flat_r_squared_share = share_distance(draw_coefficients_from_prior(ks_draws, covariance, "flat_r_squared")))
})))
c9_tests <- nrow(c9_lkj) + 3 * nrow(c9_draws)
c9_critical <- ks_critical(0.01 / c9_tests, ks_draws)
record_check(paste("2.15 block of q columns: flat_partial_r gives partial R^2 ~ Beta(q/2, 1)",
                   "(LKJ and g-prior draws); flat_r and flat_r_squared give share ~ Beta(q/2, (k-q)/2)"),
             all(c(c9_lkj$distance, c9_draws$flat_partial_r, c9_draws$flat_r_share,
                   c9_draws$flat_r_squared_share) < c9_critical),
             sprintf("largest KS distance: LKJ %.4f, g-prior draws %.4f, shares %.4f, %.4f; critical %.4f",
                     max(c9_lkj$distance), max(c9_draws$flat_partial_r), max(c9_draws$flat_r_share),
                     max(c9_draws$flat_r_squared_share), c9_critical))

## 2.16 [C10] the block test against integrate(), sharing no code with the
## package: BF = integral of L(f) p(f) df, L(f) = exp(-n f^2/2)
## E_theta[I(b cos theta)] / I(0), b = sqrt(n R2) f, theta the angle between
## the effect and the data's direction (density ~ sin^(q-2) theta),
## m = n - k + q - 2; flat_partial_r p(f) = q f^(q-1) (1 + f^2)^(-q/2-1), the
## others the share mixture with s ~ Beta(q/2, (k-q)/2). Within 2e-4.
log_bayes_factor_block_by_integration <- function(partial_r_squared, n, k, q, prior) {
  m <- n - k + q - 2
  log_I_0 <- log_I_by_integration(0, m)
  likelihood_ratio <- function(f) {
    b <- sqrt(n * partial_r_squared) * f
    if (b == 0) return(exp(-n * f^2 / 2))
    log_scale <- log_I_by_integration(b, m)
    direction_weight <- function(theta) sin(theta)^(q - 2) / beta(0.5, (q - 1) / 2)
    pieces <- unique(pmin(pi, c(0, 2, 8, 30) / sqrt(b)))
    pieces <- c(pieces[pieces < pi], pi)
    average <- sum(sapply(seq_len(length(pieces) - 1), function(i)
      integrate(Vectorize(function(theta)
        direction_weight(theta) * exp(log_I_by_integration(b * cos(theta), m) - log_scale)),
        pieces[i], pieces[i + 1], rel.tol = 1e-10)$value))
    exp(-n * f^2 / 2 + log_scale - log_I_0) * average
  }
  a <- prior_shape(prior, k)
  prior_density <- function(f) {
    if (prior == "flat_partial_r") return(q * f^(q - 1) * (1 + f^2)^(-q / 2 - 1))
    if (q == k) return(prior_density_given_share(f, 1, a))
    given_share <- function(s) stats::dbeta(s, q / 2, (k - q) / 2) * prior_density_given_share(f, s, a)
    integrate(function(u) given_share(exp(u)) * exp(u), log(f^2) - 50, log(0.5),
              rel.tol = 1e-10)$value + integrate(given_share, 0.5, 1, rel.tol = 1e-10)$value
  }
  integrand <- Vectorize(function(f) likelihood_ratio(f) * prior_density(f))
  partial_r <- sqrt(partial_r_squared)
  peak  <- sqrt((n - k - 1) / n) * partial_r / sqrt(1 - partial_r^2)
  width <- 1 / sqrt(n * (1 - partial_r_squared))
  breaks <- sort(unique(pmax(0, peak + c(-30, -15, -5, -2, 0, 2, 5, 15) * width)))
  log(sum(sapply(seq_len(length(breaks) - 1), function(i)
    integrate(integrand, breaks[i], breaks[i + 1], rel.tol = 1e-8)$value)))
}
c10_cases <- expand.grid(f_value = c(0, 1, 3, 6), q = c(2, 3, 5), k = c(5, 12), n = c(50, 1e4, 1e6),
                         prior = priors, stringsAsFactors = FALSE)
c10_cases <- subset(c10_cases, !(q == 5 & k == 5) | prior != "flat_partial_r")
if (smoke_mode) c10_cases <- subset(c10_cases, q == 2 & k == 5 & n == 50)
c10_cases <- cached("block_grid", {
  c10_cases$difference <- unlist(parallel_map(seq_len(nrow(c10_cases)), function(row) {
    case <- c10_cases[row, ]
    ## F statistic F_value for the block: partial R^2 = q F / (q F + n - k - 1)
    partial_r_squared <- case$q * case$f_value / (case$q * case$f_value + case$n - case$k - 1)
    ours <- log_bayes_factor_conditional(sqrt(partial_r_squared), case$n, case$k, case$prior,
                                         number_of_columns = case$q)
    abs(ours - log_bayes_factor_block_by_integration(partial_r_squared, case$n, case$k, case$q, case$prior))
  }))
  c10_cases
})
worst <- c10_cases[which.max(c10_cases$difference), ]
record_check("2.16 block conditional test (every prior) = integrate() (q 2 to 5, n 50 to 1e6, F 0 to 6)",
             max(c10_cases$difference) < 2e-4,
             sprintf("largest difference %.1e (%s, n = %g, k = %d, q = %d, F = %g)", max(c10_cases$difference),
                     worst$prior, worst$n, worst$k, worst$q, worst$f_value))

## 2.17 [C11] the block likelihood ratio from the data themselves: integrate
## over sigma (prior 1/sigma) the likelihood of the residualised outcome given
## the effect (whitened by the block's sample covariance given the others),
## and average over the effect's direction (q = 2). Within 1e-8.
set.seed(11)
c11_differences <- unlist(lapply(list(c(40, 4), c(80, 5), c(30, 2)), function(design) {
  n <- design[1]; k <- design[2]; q <- 2
  X <- matrix(rnorm(n * k), n) %*% chol(0.5 + 0.5 * diag(k))
  y <- drop(X %*% seq(0.1, 0.4, length.out = k)) + rnorm(n)
  others <- qr(cbind(1, X[, seq_len(k - q), drop = FALSE]))
  y_left <- qr.resid(others, y); Z <- qr.resid(others, X[, (k - q + 1):k, drop = FALSE])
  partial_r_squared <- 1 - sum(qr.resid(qr(cbind(1, X)), y)^2) / sum(y_left^2)
  eigen_Z <- eigen(crossprod(Z) / n)
  unwhiten <- eigen_Z$vectors %*% diag(1 / sqrt(eigen_Z$values)) %*% t(eigen_Z$vectors)
  degrees <- n - (k - q) - 1
  sigma_hat <- sqrt(sum(y_left^2) / degrees)
  log_marginal <- function(gamma) {
    log_integrand <- function(sigma) -(degrees + 1) * log(sigma) -
      sum((y_left - Z %*% (sigma * unwhiten %*% gamma))^2) / (2 * sigma^2)
    log_integrand(sigma_hat) + log(integrate(Vectorize(function(sigma)
      exp(log_integrand(sigma) - log_integrand(sigma_hat))), 0, Inf, rel.tol = 1e-12)$value)
  }
  log_null <- log_marginal(c(0, 0))
  sapply(c(0.05, 0.2, 0.5), function(f) {
    by_integration <- log(integrate(Vectorize(function(angle)
      exp(log_marginal(f * c(cos(angle), sin(angle))) - log_null)), 0, 2 * pi,
      rel.tol = 1e-10)$value / (2 * pi))
    abs(log_likelihood_ratio_block(f, partial_r_squared, n, k, q) - by_integration)
  })
}))
record_check("2.17 block likelihood ratio = integration over sigma and direction from the data (q = 2)",
             max(c11_differences) < 1e-8, sprintf("largest difference %.1e", max(c11_differences)))


save_study_results("02", list())
report_checks()
