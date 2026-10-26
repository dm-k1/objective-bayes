## =====================================================================
##  Testing one coefficient or one term: the conditional test, for both priors
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
##
## A term with q columns (a factor with q + 1 levels) is tested as a block, in
## the same way: f is the length of its effect in noise units, its direction
## within the block is uniform, and r^2 is the block's sample partial R^2.
##   flat-R                share ~ Beta(q/2, (k-q)/2)
##   uniform correlations  partial R^2 ~ Beta(q/2, 1), for every k
## The likelihood ratio becomes exp(-n f^2/2) E_c[I(sqrt(n) r f c)] / I(0) with
## m = n - k + q - 2 and c the cosine between the effect and the data's
## direction (for q = 1, c = +-1: the formula above). Checks: validation/05,
## C9 and C10.

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
## ~ Beta(q/2, (k-q)/2) (q columns tested; q = 1 for one coefficient) on a grid
## even in logit(share): it reaches down to shares far below the smallest f^2,
## where the CDF changes fast, and up to shares near 1, where the Beta density
## is infinite when k - q = 1. Each grid cell carries its exact Beta
## probability (pbeta); the cells below and above the grid get share 0 and
## share 1. (At f = 0 the CDF is 0; f / sqrt(0 + 0) is not.) For a block,
## |f| is the length of the block's effect in noise units.
prior_cdf_of_f <- function(f_values, number_of_predictors, points_share, number_of_columns = 1) {
  if (number_of_predictors == number_of_columns) return(f_values / sqrt(1 + f_values^2))
  smallest_share <- 1e-4 * min(f_values[f_values > 0])^2
  logit_edges <- seq(stats::qlogis(smallest_share), 25, length.out = points_share)
  share_probability <- diff(c(0, stats::pbeta(stats::plogis(logit_edges), number_of_columns / 2,
                                              (number_of_predictors - number_of_columns) / 2), 1))
  share_values <- c(0, stats::plogis((logit_edges[-1] + logit_edges[-points_share]) / 2), 1)
  cdf_given_share <- outer(share_values, f_values, function(share, f)
    ifelse(f == 0, 0, f / sqrt(share + f^2)))
  colSums(share_probability * cdf_given_share)
}

## Block of q columns: the effect is a vector gamma in noise units (whitened by
## the block's covariance given the other predictors), |gamma| = f, with a
## uniform direction. Integrating out sigma as for one column leaves
## exp(-n f^2/2) E_c[I(sqrt(n) r f c)] / I(0), with r^2 the block's sample
## partial R^2, c the cosine between gamma and the data's direction, and
## m = n - k + q - 2. The average over directions is the 0F1 below:
##   E_c[exp(t c)] = 0F1(; q/2; t^2/4) = Gamma(q/2) (t/2)^(1 - q/2) I_{q/2-1}(t)
## (Bessel I; cosh(t) when q = 1, the two signs). Its log, for t >= 0:
log_direction_average <- function(t, number_of_columns) {
  if (number_of_columns == 1) return(t + log1p(exp(-2 * t)) - log(2))
  order <- number_of_columns / 2 - 1
  small <- t < 1e-6
  t_large <- ifelse(small, 1, t)
  ifelse(small, t^2 / (2 * number_of_columns),
         lgamma(number_of_columns / 2) - order * log(t_large / 2) +
           log(besselI(t_large, order, expon.scaled = TRUE)) + t_large)
}

## log of integral_0^inf v^m exp(-v^2/2) E_c[exp(b v c)] dv for a vector of
## b >= 0: a grid around the peak of the integrand, as log_peak_integral().
## For large b v the direction average grows like exp(b v) v^(-(q-1)/2), so the
## peak is that of I(b) with m lowered by (q - 1)/2.
log_peak_integral_block <- function(b, m, number_of_columns, half_width = 12, points = 241) {
  m_peak  <- pmax(1, m - (number_of_columns - 1) / 2 * (b > 0))
  peak    <- (b + sqrt(b^2 + 4 * m_peak)) / 2
  width   <- 1 / sqrt(1 + m_peak / peak^2)
  offsets <- seq(-half_width, half_width, length.out = points)
  v <- matrix(peak, length(b), points) + outer(width, offsets)
  v[v <= 0] <- NA
  log_integrand <- m * log(v) - v^2 / 2 +
    matrix(log_direction_average(as.vector(b * v), number_of_columns), length(b), points)
  largest <- apply(log_integrand, 1, max, na.rm = TRUE)
  largest + log(rowSums(exp(log_integrand - largest), na.rm = TRUE)) +
    log(width * (offsets[2] - offsets[1]))
}

## Log likelihood ratio p(y | f) / p(y | 0) for a block of q columns, for a
## vector of f: depends only on the block's sample partial R^2, n, k and q.
## (In closed form exp(-n f^2/2) 1F1((m+1)/2; q/2; n r^2 f^2 / 2), which is hard
## to evaluate for large n; the integral over v is not.)
log_likelihood_ratio_block <- function(f_values, partial_r_squared, sample_size,
                                       number_of_predictors, number_of_columns) {
  m <- sample_size - number_of_predictors + number_of_columns - 2
  b <- sqrt(sample_size * partial_r_squared) * f_values
  -sample_size * f_values^2 / 2 + log_peak_integral_block(b, m, number_of_columns) -
    log_peak_integral_block(0, m, number_of_columns)
}

## Log Bayes factor for one coefficient or a block of q columns, for any prior
## on |f| given by its CDF (a function of a vector of f values). Used by both
## conditional tests. partial_r: the sample partial correlation, or for a block
## the square root of its sample partial R^2.
log_bayes_factor_from_prior_cdf_of_f <- function(partial_r, sample_size, number_of_predictors,
                                                 prior_cdf, points_f = 600, number_of_columns = 1) {
  m <- sample_size - number_of_predictors - 1

  ## 1. A grid in f around the peak of the likelihood ratio, which for one
  ##    column is at sqrt(m / n) r / sqrt(1 - r^2) (where its log has slope 0),
  ##    with width about 1 / sqrt(n (1 - r^2)). Past 15 widths the ratio is
  ##    negligible. The grid follows the peak, so it stays fine enough when n
  ##    is large. For a block the peak lies lower, by at most sqrt(q) widths
  ##    (the noncentral chi-square's mean exceeds its noncentrality by q).
  peak_f  <- sqrt(m / sample_size) * abs(partial_r) / sqrt(1 - partial_r^2)
  width_f <- 1 / sqrt(sample_size * (1 - partial_r^2))
  f_edges <- seq(max(0, peak_f - (15 + 2 * sqrt(number_of_columns - 1)) * width_f),
                 peak_f + 15 * width_f, length.out = points_f + 1)
  f_middles <- (f_edges[-1] + f_edges[-length(f_edges)]) / 2

  ## 2. The likelihood ratio in the middle of each cell (for one column: both
  ##    signs of f; for a block: every direction)
  log_likelihood_ratio <- log_likelihood_ratio_block(f_middles, partial_r^2, sample_size,
                                                     number_of_predictors, number_of_columns)

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
#'
#' A block of q columns (a factor with q + 1 levels, a multi-column
#' interaction) is tested the same way: its share of the flat total is
#' Beta(q/2, (k - q)/2), its direction within the block uniform, and the
#' evidence depends only on the block's sample partial R-squared, n, k and q.
#' @param partial_r Sample partial correlation of the tested predictor; for a
#'   block, the square root of the block's sample partial R-squared.
#' @param sample_size Number of observations, n.
#' @param number_of_predictors Number of predictor columns in the model, k.
#' @param points_f,points_share Grid sizes: for f, around the peak of the
#'   likelihood, and for the share, inside the prior. With the defaults the
#'   error is below about 2e-4 on the log scale for n = 50 to 1e6, k = 1 to 30
#'   (against \code{stats::integrate()}, validation/05, check C5).
#' @param number_of_columns Number of columns tested together, q (1 for one
#'   coefficient).
#' @return The log Bayes factor for an effect against no effect.
#' @examples
#' log_bayes_factor_flat_r_conditional(0.2, sample_size = 100, number_of_predictors = 3)
#' @export
log_bayes_factor_flat_r_conditional <- function(partial_r, sample_size, number_of_predictors,
                                                points_f = 600, points_share = 500,
                                                number_of_columns = 1) {
  if (sample_size < number_of_predictors + 2)
    stop("Need at least two more observations than predictors.")
  log_bayes_factor_from_prior_cdf_of_f(
    partial_r, sample_size, number_of_predictors,
    function(f_values) prior_cdf_of_f(f_values, number_of_predictors, points_share,
                                      number_of_columns),
    points_f, number_of_columns)
}

#' Log Bayes factor for one coefficient (conditional test, uniform correlations)
#'
#' Under uniform correlations (every correlation matrix of the outcome and the
#' predictors equally likely) one coefficient's partial correlation is flat on
#' (-1, 1) whatever the number of predictors: Jeffreys' flat correlation for
#' every coefficient (validation/05, checks C4 and C7). So this is the
#' conditional test with the prior P(|f| < f) = f / sqrt(1 + f^2), for any k.
#' The other coefficients get flat priors, as in the flat-R version.
#'
#' For a block of q columns the block's partial R-squared is Beta(q/2, 1),
#' again for any k (validation/05, check C9), so
#' P(|f| < f) = (f^2 / (1 + f^2))^(q/2), with a uniform direction within the block.
#' @inheritParams log_bayes_factor_flat_r_conditional
#' @return The log Bayes factor for an effect against no effect.
#' @examples
#' log_bayes_factor_uniform_correlations_conditional(0.2, sample_size = 100, number_of_predictors = 3)
#' @export
log_bayes_factor_uniform_correlations_conditional <- function(partial_r, sample_size,
                                                              number_of_predictors, points_f = 600,
                                                              number_of_columns = 1) {
  if (sample_size < number_of_predictors + 2)
    stop("Need at least two more observations than predictors.")
  log_bayes_factor_from_prior_cdf_of_f(
    partial_r, sample_size, number_of_predictors,
    function(f_values) (f_values^2 / (1 + f_values^2))^(number_of_columns / 2),
    points_f, number_of_columns)
}
