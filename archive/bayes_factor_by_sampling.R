## ARCHIVED: brute-force Bayes factor by sampling from the prior. Was an
## independent check of the closed form; that role is now played by brms bridge
## sampling (validation/02_bayes_factor.R, check A3). Needs R/flat_r.R sourced.

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
