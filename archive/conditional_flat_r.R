## =====================================================================
##  Helper: the per-coefficient (conditional) flat-R test
##  Used by validation/06_broad_simulation.R. Derivation and checks:
##  validation/05_conditional_test.R (the same functions, defined there too).
## =====================================================================

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

## The same test with a FLAT prior on the partial correlation (share = 1): the
## exact version of the per-coefficient flat partial r in comparison_methods.R
log_bayes_factor_flat_partial_exact <- function(partial_r, sample_size, number_of_predictors,
                                                points_R = 2000) {
  m <- sample_size - number_of_predictors - 1
  R <- (seq_len(points_R) - 0.5) / points_R
  f <- R / sqrt(1 - R^2)
  b <- sqrt(sample_size) * partial_r * f
  log_null <- log_peak_integral(0, m)
  log_ratio <- c(-sample_size * f^2 / 2 + log_peak_integral( b, m) - log_null,
                 -sample_size * f^2 / 2 + log_peak_integral(-b, m) - log_null)
  largest <- max(log_ratio)
  largest + log(mean(exp(log_ratio - largest)))
}

## Sample partial correlation of one coefficient, from its t-statistic
partial_correlation_of <- function(lm_fit, coefficient_row) {
  t_value <- summary(lm_fit)$coefficients[coefficient_row, "t value"]
  t_value / sqrt(t_value^2 + lm_fit$df.residual)
}
