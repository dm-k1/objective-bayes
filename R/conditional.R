## =====================================================================
##  Testing one coefficient: the conditional flat-R test
## =====================================================================
## In one sentence: test one coefficient's partial correlation, with the prior
## that the joint flat-R prior implies for one coefficient, so the evidence
## depends only on that coefficient's own contribution (its sample partial
## correlation), n and k. This avoids the conditional Lindley paradox of priors
## with one shared scale (Som, Hans & MacEachern, 2016), which affects the
## joint test as it does the JZS default. docs/2_how_it_calculates.md, section 9.
##
## Prior for one coefficient, in noise units f:
##   f^2 = share * R^2 / (1 - R^2),  R flat on (0, 1),  share ~ Beta(1/2, (k-1)/2),
## sign + or - with equal probability. It depends only on k, not on how the
## predictors correlate. The other coefficients and the intercept get flat
## priors and the noise 1/sigma (shared by both hypotheses). Integrating them
## out leaves
##   p(y | f) / p(y | 0) = exp(-n f^2 / 2) * I(sqrt(n) r f) / I(0),
##   I(b) = integral_0^inf v^m exp(-v^2/2 + b v) dv,   m = n - k - 1,
## with r the sample partial correlation. Checks: validation/05.

## log I(b) for a vector of b: a grid around the peak of the integrand
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

#' Log Bayes factor for one coefficient (conditional flat-R test)
#'
#' Tests one coefficient's partial correlation, with the prior that the joint
#' flat-R prior implies for one of \code{number_of_predictors} coefficients.
#' Depends only on the sample partial correlation, n and k. With one predictor
#' it equals the joint test.
#' @param partial_r Sample partial correlation of the tested predictor.
#' @param sample_size Number of observations, n.
#' @param number_of_predictors Number of predictor columns in the model, k.
#' @param points_R,points_share Grid sizes of the prior (accuracy about 1e-4 on
#'   the log scale with the defaults).
#' @return The log Bayes factor for an effect against no effect.
#' @examples
#' log_bayes_factor_flat_r_conditional(0.2, sample_size = 100, number_of_predictors = 3)
#' @export
log_bayes_factor_flat_r_conditional <- function(partial_r, sample_size, number_of_predictors,
                                                points_R = 100, points_share = 50) {
  if (sample_size < number_of_predictors + 2)
    stop("Need at least two more observations than predictors.")
  m <- sample_size - number_of_predictors - 1
  R <- (seq_len(points_R) - 0.5) / points_R
  share <- if (number_of_predictors == 1) 1 else
    stats::qbeta((seq_len(points_share) - 0.5) / points_share, 0.5, (number_of_predictors - 1) / 2)
  f <- sqrt(as.vector(outer(R^2 / (1 - R^2), share)))
  b <- sqrt(sample_size) * partial_r * f
  log_null <- log_peak_integral(0, m)
  log_ratio <- c(-sample_size * f^2 / 2 + log_peak_integral( b, m) - log_null,
                 -sample_size * f^2 / 2 + log_peak_integral(-b, m) - log_null)
  log_sum_exp(log_ratio) - log(length(log_ratio))
}
