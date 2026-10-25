## =====================================================================
##  Comparison methods used in the validation simulations
##  (the flat-R method itself lives in R/flat_r.R)
##
##   - JZS default (BayesFactor / JASP)
##   - per-coefficient flat prior on one partial correlation
##   - R2-D2 with the same R-squared prior, Dirichlet split
##
##  Names follow R/flat_r.R: bayes_factor_<method>[_<scope>].
## =====================================================================

## ---- JZS default: one coefficient, model with vs without it ----------
## rscale sqrt(2)/4 is BayesFactor's "medium" for regression.
bayes_factor_jzs_one_coefficient <- function(explained_share_full, explained_share_without,
                                             sample_size, number_of_predictors,
                                             rscale = sqrt(2) / 4) {
  with_coefficient <- BayesFactor::linearReg.R2stat(
    sample_size, number_of_predictors, explained_share_full,
    rscale = rscale, simple = TRUE)
  without_coefficient <- if (number_of_predictors == 1) 1 else
    BayesFactor::linearReg.R2stat(sample_size, number_of_predictors - 1,
                                  explained_share_without, rscale = rscale, simple = TRUE)
  with_coefficient / without_coefficient
}

## ---- Per-coefficient flat prior on the partial correlation -----------
## Flat prior on one partial correlation from -1 to 1 (Jeffreys' default
## applied to one coefficient at a time). The likelihood is evaluated on
## the slope via the partial-correlation <-> slope mapping, using the
## estimated noise level (so it is not exact; see Kucharsky et al. 2026
## for the exact version).
partial_correlation_grid <- seq(-1, 1, length.out = 8003)[-c(1, 8003)]
partial_correlation_step <- partial_correlation_grid[2] - partial_correlation_grid[1]

## Shared by the Bayes factor and the posterior below: the slope that
## belongs to each grid correlation, and the likelihood of each slope.
likelihood_over_partial_correlation_grid <- function(lm_fit, predictor_index, predictors) {
  coefficient_row <- summary(lm_fit)$coefficients[predictor_index + 1, ]
  estimate        <- coefficient_row[[1]]
  standard_error  <- coefficient_row[[2]]
  degrees_of_freedom <- lm_fit$df.residual
  spread_left_after_others <- if (ncol(predictors) > 1)
    sd(resid(lm(predictors[, predictor_index] ~ predictors[, -predictor_index]))) else
    sd(predictors[, predictor_index])
  conversion_to_slope <- spread_left_after_others / summary(lm_fit)$sigma
  slopes <- partial_correlation_grid /
            (conversion_to_slope * sqrt(1 - partial_correlation_grid^2))
  list(slopes     = slopes,
       likelihood = dt((estimate - slopes) / standard_error, degrees_of_freedom) / standard_error,
       likelihood_at_zero = dt(estimate / standard_error, degrees_of_freedom) / standard_error)
}

bayes_factor_flat_partial_correlation <- function(lm_fit, predictor_index, predictors) {
  grid <- likelihood_over_partial_correlation_grid(lm_fit, predictor_index, predictors)
  evidence_effect <- sum(grid$likelihood * 0.5) * partial_correlation_step   # flat: density 1/2
  evidence_effect / grid$likelihood_at_zero
}

## Posterior mean and 95% interval, for the slope and the partial correlation
## (used by archive/sim_one_coefficient_flat_partial_r.R)
posterior_flat_partial_correlation <- function(lm_fit, predictor_index, predictors) {
  grid <- likelihood_over_partial_correlation_grid(lm_fit, predictor_index, predictors)
  posterior_weights <- grid$likelihood / sum(grid$likelihood)   # flat prior: posterior ~ likelihood
  cumulative        <- cumsum(posterior_weights)
  lower_index <- which.max(cumulative >= 0.025)
  upper_index <- which.max(cumulative >= 0.975)
  list(slope_mean        = sum(grid$slopes * posterior_weights),
       slope_lower       = grid$slopes[lower_index],
       slope_upper       = grid$slopes[upper_index],
       correlation_mean  = sum(partial_correlation_grid * posterior_weights),
       correlation_lower = partial_correlation_grid[lower_index],
       correlation_upper = partial_correlation_grid[upper_index])
}

## ---- R2-D2 with the same R-squared prior ------------------------------
## coefficient_j | sigma ~ N(0, sigma^2 * omega * share_j / var(x_j)),
## omega = R2 / (1 - R2), R2 ~ Beta(1/2, 1), shares ~ Dirichlet(1/2).
## Bayes factor by Monte Carlo over shared prior draws (normal kernel is
## integrated exactly per draw).
draw_from_r2d2_prior <- function(number_of_predictors, number_of_draws = 1000,
                                 r2_shape_1 = 0.5, r2_shape_2 = 1, dirichlet_concentration = 0.5) {
  explained_share <- rbeta(number_of_draws, r2_shape_1, r2_shape_2)
  shares <- matrix(rgamma(number_of_draws * number_of_predictors, dirichlet_concentration),
                   number_of_draws, number_of_predictors)
  list(signal_to_noise = explained_share / (1 - explained_share),
       shares = shares / rowSums(shares))
}

## Log Bayes factor of k centered predictors vs intercept only
log_bayes_factor_r2d2_vs_intercept_only <- function(centered_predictors, centered_outcome,
                                                    prior_draws) {
  number_of_predictors <- ncol(centered_predictors)
  sample_size          <- length(centered_outcome)
  outcome_sum_of_squares <- sum(centered_outcome^2)
  cross_products   <- crossprod(centered_predictors)
  predictor_outcome <- drop(crossprod(centered_predictors, centered_outcome))
  predictor_variances <- colSums(centered_predictors^2) / sample_size

  log_evidence_per_draw <- vapply(seq_along(prior_draws$signal_to_noise), function(m) {
    prior_sd <- sqrt(prior_draws$signal_to_noise[m] * prior_draws$shares[m, ] / predictor_variances)
    root <- chol(diag(number_of_predictors) + (prior_sd %o% prior_sd) * cross_products)
    projected <- backsolve(root, prior_sd * predictor_outcome, transpose = TRUE)
    -sum(log(diag(root))) - (sample_size - 1) / 2 * log1p(-sum(projected^2) / outcome_sum_of_squares)
  }, numeric(1))
  log_sum_exp(log_evidence_per_draw) - log(length(log_evidence_per_draw))
}
