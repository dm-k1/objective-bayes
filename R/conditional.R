## =====================================================================
##  Testing one coefficient: the conditional test, for both priors
## =====================================================================
## In one sentence: test one coefficient's partial correlation, with the prior
## that the joint flat-R prior implies for one coefficient, so the evidence
## depends only on that coefficient's own contribution (its sample partial
## correlation), n and k. This avoids the conditional Lindley paradox of priors
## with one shared scale (Som, Hans & MacEachern, 2016), which affects the
## joint test as it does the JZS default. docs/2_how_it_calculates.md, section 10.
##
## Prior for one coefficient, in noise units f:
##   f^2 = share * R^2 / (1 - R^2),  R flat on (0, 1),  share ~ Beta(1/2, (k-1)/2),
## sign + or - with equal probability. It depends only on k, not on how the
## predictors correlate. The other coefficients and the intercept get flat
## priors and the noise 1/sigma (shared by both hypotheses). Integrating them
## out leaves
##   p(y | f) / p(y | 0) = exp(-n f^2 / 2) * I(sqrt(n) r f) / I(0),
##   I(b) = integral_0^inf v^m exp(-v^2/2 + b v) dv,   m = n - k - 1,
## with r the sample partial correlation. The average over the prior is taken
## on a grid in f around the peak of this ratio (width about 1/sqrt(n)), each
## cell weighted by its exact prior probability. Checks: validation/05.

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

## Prior probability that |f| is below each of f_values (the prior's CDF).
## For one share s: |f| = sqrt(s) R / sqrt(1 - R^2) with R flat, so
## P(|f| < f | s) = f / sqrt(s + f^2) exactly. This is averaged over the share
## ~ Beta(1/2, (k-1)/2) on a grid even in logit(share): it reaches down to
## shares far below the smallest f^2, where the CDF changes fast, and up to
## shares near 1, where the Beta density is infinite when k = 2. Each grid cell
## carries its exact Beta probability (pbeta); the cells below and above the
## grid get share 0 and share 1. (At f = 0 the CDF is 0; f / sqrt(0 + 0) is not.)
prior_cdf_of_f <- function(f_values, number_of_predictors, points_share) {
  if (number_of_predictors == 1) return(f_values / sqrt(1 + f_values^2))
  smallest_share <- 1e-4 * min(f_values[f_values > 0])^2
  logit_edges <- seq(stats::qlogis(smallest_share), 25, length.out = points_share)
  share_probability <- diff(c(0, stats::pbeta(stats::plogis(logit_edges), 0.5,
                                              (number_of_predictors - 1) / 2), 1))
  share_values <- c(0, stats::plogis((logit_edges[-1] + logit_edges[-points_share]) / 2), 1)
  cdf_given_share <- outer(share_values, f_values, function(share, f)
    ifelse(f == 0, 0, f / sqrt(share + f^2)))
  colSums(share_probability * cdf_given_share)
}

## Log Bayes factor for one coefficient, for any prior on |f| given by its
## CDF (a function of a vector of f values). Used by both conditional tests.
log_bayes_factor_from_prior_cdf_of_f <- function(partial_r, sample_size, number_of_predictors,
                                                 prior_cdf, points_f = 600) {
  m <- sample_size - number_of_predictors - 1

  ## 1. A grid in f around the peak of the likelihood ratio, which is at
  ##    sqrt(m / n) r / sqrt(1 - r^2) (where its log has slope 0) with width
  ##    about 1 / sqrt(n (1 - r^2)). Past 15 widths the ratio is negligible.
  ##    The grid follows the peak, so it stays fine enough when n is large.
  peak_f  <- sqrt(m / sample_size) * abs(partial_r) / sqrt(1 - partial_r^2)
  width_f <- 1 / sqrt(sample_size * (1 - partial_r^2))
  f_edges <- seq(max(0, peak_f - 15 * width_f), peak_f + 15 * width_f,
                 length.out = points_f + 1)
  f_middles <- (f_edges[-1] + f_edges[-length(f_edges)]) / 2

  ## 2. The likelihood ratio in the middle of each cell, for both signs of f
  b <- sqrt(sample_size) * partial_r * f_middles
  log_null <- log_peak_integral(0, m)
  log_ratio_by_sign <- cbind(log_peak_integral( b, m), log_peak_integral(-b, m)) - log_null
  log_likelihood_ratio <- -sample_size * f_middles^2 / 2 +
    log_sum_exp_per_row(log_ratio_by_sign) - log(2)

  ## 3. Averaged over the prior: each cell weighs with its exact prior
  ##    probability (a difference of the prior's CDF), so where the prior
  ##    piles up (for flat-R, near f = 0 for small shares) needs no fine grid
  prior_probability <- diff(prior_cdf(f_edges))
  log_sum_exp(log_likelihood_ratio + log(prior_probability))
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
#' @param points_f,points_share Grid sizes: for f, around the peak of the
#'   likelihood, and for the share, inside the prior. With the defaults the
#'   error is below about 2e-4 on the log scale for n = 50 to 1e6, k = 1 to 30
#'   (against \code{stats::integrate()}, validation/05, check C5).
#' @return The log Bayes factor for an effect against no effect.
#' @examples
#' log_bayes_factor_flat_r_conditional(0.2, sample_size = 100, number_of_predictors = 3)
#' @export
log_bayes_factor_flat_r_conditional <- function(partial_r, sample_size, number_of_predictors,
                                                points_f = 600, points_share = 500) {
  if (sample_size < number_of_predictors + 2)
    stop("Need at least two more observations than predictors.")
  log_bayes_factor_from_prior_cdf_of_f(
    partial_r, sample_size, number_of_predictors,
    function(f_values) prior_cdf_of_f(f_values, number_of_predictors, points_share), points_f)
}

#' Log Bayes factor for one coefficient (conditional test, uniform correlations)
#'
#' Under uniform correlations (every correlation matrix of the outcome and the
#' predictors equally likely) one coefficient's partial correlation is flat on
#' (-1, 1) whatever the number of predictors: Jeffreys' flat correlation for
#' every coefficient (validation/05, checks C4 and C7). So this is the
#' conditional test with the prior P(|f| < f) = f / sqrt(1 + f^2), for any k.
#' The other coefficients get flat priors, as in the flat-R version.
#' @inheritParams log_bayes_factor_flat_r_conditional
#' @return The log Bayes factor for an effect against no effect.
#' @examples
#' log_bayes_factor_uniform_correlations_conditional(0.2, sample_size = 100, number_of_predictors = 3)
#' @export
log_bayes_factor_uniform_correlations_conditional <- function(partial_r, sample_size,
                                                              number_of_predictors, points_f = 600) {
  if (sample_size < number_of_predictors + 2)
    stop("Need at least two more observations than predictors.")
  log_bayes_factor_from_prior_cdf_of_f(partial_r, sample_size, number_of_predictors,
                                       function(f_values) f_values / sqrt(1 + f_values^2), points_f)
}
