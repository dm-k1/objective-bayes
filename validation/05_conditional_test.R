## =====================================================================
##  Validation 05: checks of the per-coefficient (conditional) tests
##
##  The conditional test looks only at one coefficient's sample partial
##  correlation, with the prior that the joint prior implies for one
##  coefficient (R/conditional.R; docs/2_how_it_calculates.md, section 10):
##    flat-R                f_j^2 = share * R^2 / (1 - R^2), R flat,
##                          share ~ Beta(1/2, (k - 1)/2)
##    uniform correlations  f_j^2 = W chi^2_1, W ~ Inverse-Gamma(1, 1/2):
##                          a flat partial correlation, whatever k is
##  How the tests compare with the other methods is in studies 03, 04 and 06,
##  which use the common method set of validation/comparison_methods.R.
##
##    C1. k = 1: the conditional flat-R test equals the joint flat-R test.
##    C2. the implied flat-R prior for one coefficient does not depend on how
##        the predictors correlate: correlated joint draws against the formula.
##    C3. the exact flat-partial-correlation test against the approximate one
##        that plugs in the estimated noise (close, not equal).
##    C4. under uniform correlations the implied prior for one coefficient's
##        partial correlation is flat, for any k and correlated predictors.
##    C5. the conditional flat-R grid against stats::integrate() over
##        n = 50 to 1e6, k = 1 to 30 and t-values 0 to 6.
##    C6. the same for the uniform-correlations conditional test.
##    C7. why eta = 1: under LKJ(eta) the partial correlation of y with one
##        predictor, given all the others, is Beta(eta, eta) on (-1, 1) for any
##        k. It is flat (Jeffreys) for every k only at eta = 1.
##    C8. against the exact partial-correlation Bayes factor of Kucharský et al.
##        (2026) with a flat prior: not equal (theirs treats the predictors as
##        random), and equal once the predictor's residual SD is unknown.
##    C9. a block of q columns (a factor with q + 1 levels): the implied prior.
##        Uniform correlations: its partial R^2 is Beta(q/2, 1) for any k (LKJ
##        draws and g-prior draws); flat-R: its share is Beta(q/2, (k-q)/2).
##    C10. the block test (both priors) against stats::integrate().
##    C11. the block likelihood ratio against integration over sigma and the
##        effect's direction from simulated data.
##
##  Run from the project folder: Rscript validation/05_conditional_test.R
##  Under a minute.
## =====================================================================
source("validation/settings.R")   # record_check()
pkgload::load_all(".", quiet = TRUE)
source("validation/comparison_methods.R")


## C1. With one predictor the test is the joint flat-R test (both are flat on r).
c1_differences <- sapply(list(c(0.3, 30), c(0.1, 100), c(0.5, 50), c(0.05, 500)), function(s)
  abs(log_bayes_factor_flat_r_conditional(s[1], s[2], 1) -
      log_bayes_factor_flat_r_vs_intercept_only(s[1]^2, s[2], 1)))
record_check("C1. k = 1: equals the joint flat-R Bayes factor", max(c1_differences) < 1e-3,
             sprintf("largest difference %.1e", max(c1_differences)))

## C2. The implied prior for one coefficient's partial correlation does not
##     depend on how the predictors correlate.
set.seed(1)
correlated <- matrix(0.6, 4, 4); diag(correlated) <- 1
joint_draws <- draw_coefficients_from_prior(200000, correlated, prior = "flat_r")
spread_left <- sqrt(1 / diag(solve(correlated)))[2]
from_joint <- joint_draws[, 2] * spread_left / sqrt(1 + (joint_draws[, 2] * spread_left)^2)
f_formula  <- sqrt(rbeta(200000, 0.5, 1.5) * (function(R) R^2 / (1 - R^2))(runif(200000)))
from_formula <- sample(c(-1, 1), 200000, replace = TRUE) * f_formula / sqrt(1 + f_formula^2)
c2_distance <- suppressWarnings(unname(ks.test(from_joint, from_formula)$statistic))
record_check("C2. flat-R implied prior: correlated joint draws = share formula (k = 4)", c2_distance < 0.01,
             sprintf("KS distance %.4f", c2_distance))

## C3. Exact flat partial correlation (the uniform-correlations conditional test)
##     against the approximate version. Criterion: within 0.05 on the log scale.
set.seed(3)
c3_differences <- sapply(c(60, 100, 200, 400), function(n) {
  X <- matrix(rnorm(n * 3), n) %*% chol(0.5 + 0.5 * diag(3))
  y <- drop(X %*% c(0.15, 0.3, 0)) + rnorm(n)
  fit <- lm(y ~ X)
  abs(log_bayes_factor_uniform_correlations_conditional(partial_correlation_of(fit, 2), n, 3) -
      log(bayes_factor_flat_partial_correlation(fit, 1, X)))
})
record_check("C3. exact vs approximate flat partial correlation", max(c3_differences) < 0.05,
             sprintf("largest difference %.3f", max(c3_differences)))

## C4. Uniform correlations is the g-prior b ~ N(0, W Sigma^-1) (noise SD 1),
##     W ~ Inverse-Gamma(1, 1/2). One coefficient's partial correlation,
##     f_j / sqrt(1 + f_j^2) with f_j = b_j / sqrt((Sigma^-1)_jj), should be flat
##     on (-1, 1) for any k, with correlated predictors.
set.seed(4)
c4_distances <- sapply(c(2, 5, 12), function(k) {
  covariance <- matrix(0.6, k, k); diag(covariance) <- 1
  precision  <- solve(covariance)
  W <- 1 / (2 * rexp(100000))
  b <- sqrt(W) * MASS::mvrnorm(100000, rep(0, k), precision)
  f <- b[, 1] / sqrt(precision[1, 1])
  suppressWarnings(unname(ks.test(f / sqrt(1 + f^2), "punif", -1, 1)$statistic))
})
record_check("C4. uniform correlations: one partial correlation is flat (k = 2, 5, 12)",
             max(c4_distances) < 0.01,
             sprintf("KS distances %s", paste(sprintf("%.4f", c4_distances), collapse = " ")))

## C5. The conditional flat-R test against adaptive integration (stats::integrate),
##     which shares no code with the package: BF = integral of L(f) p(f) df, with
##       L(f) = exp(-n f^2/2) (I(b) + I(-b)) / (2 I(0)),  b = sqrt(n) r f,
##       p(f) = integral over the share s ~ Beta(1/2, (k-1)/2) of s (s + f^2)^(-3/2)
##     (the density of |f| = sqrt(s) R / sqrt(1 - R^2) given s, R flat).
##     The likelihood ratio peaks with width about 1/sqrt(n), so large n tests
##     whether the grid follows the peak. Criterion: within 2e-4 on the log scale.
log_I_by_integration <- function(b, m) {        # log of integral v^m exp(-v^2/2 + b v)
  peak     <- (b + sqrt(b^2 + 4 * m)) / 2
  log_peak <- m * log(peak) - peak^2 / 2 + b * peak
  log_peak + log(integrate(function(v) exp(m * log(v) - v^2 / 2 + b * v - log_peak),
                           max(0, peak - 40), peak + 40, rel.tol = 1e-10)$value)
}
##     The uniform-correlations prior for |f| is the flat-R one with k = 1
##     (share 1), whatever k is; only m = n - k - 1 changes with k.
log_bayes_factor_conditional_by_integration <- function(partial_r, n, k,
                                                        prior = c("flat-R", "uniform correlations")) {
  prior <- match.arg(prior)
  m <- n - k - 1
  likelihood_ratio <- function(f) {
    b <- sqrt(n) * partial_r * f
    log_ratio_by_sign <- -n * f^2 / 2 - log_I_by_integration(0, m) +
      c(log_I_by_integration(b, m), log_I_by_integration(-b, m))
    mean(exp(log_ratio_by_sign))
  }
  prior_density <- function(f) {
    if (k == 1 || prior == "uniform correlations") return((1 + f^2)^(-3 / 2))
    given_share <- function(s) stats::dbeta(s, 0.5, (k - 1) / 2) * s * (s + f^2)^(-3 / 2)
    ## shares below 1/2 on the log scale (the density changes near s = f^2),
    ## above 1/2 as they are (the Beta density is infinite at 1 when k = 2)
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
c5_cases <- expand.grid(t_value = c(0, 1, 2, 4, 6), k = c(1, 2, 3, 5, 10, 30),
                        n = c(50, 200, 1e4, 1e6))
c5_differences <- unlist(parallel_map(seq_len(nrow(c5_cases)), function(row) {
  n <- c5_cases$n[row]; k <- c5_cases$k[row]; t_value <- c5_cases$t_value[row]
  partial_r <- t_value / sqrt(t_value^2 + n - k - 1)
  abs(log_bayes_factor_flat_r_conditional(partial_r, n, k) -
      log_bayes_factor_conditional_by_integration(partial_r, n, k))
}))
worst <- c5_cases[which.max(c5_differences), ]
record_check("C5. conditional flat-R grid = stats::integrate (n 50 to 1e6, k 1 to 30, t 0 to 6)",
             max(c5_differences) < 2e-4,
             sprintf("largest difference %.1e (n = %g, k = %d, t = %g)", max(c5_differences),
                     worst$n, worst$k, worst$t_value))

## C6. The uniform-correlations conditional test (R/conditional.R)
##     against adaptive integration, on the same cases with k > 1.
c6_cases <- subset(c5_cases, k > 1)
c6_differences <- unlist(parallel_map(seq_len(nrow(c6_cases)), function(row) {
  n <- c6_cases$n[row]; k <- c6_cases$k[row]; t_value <- c6_cases$t_value[row]
  partial_r <- t_value / sqrt(t_value^2 + n - k - 1)
  abs(log_bayes_factor_uniform_correlations_conditional(partial_r, n, k) -
      log_bayes_factor_conditional_by_integration(partial_r, n, k, "uniform correlations"))
}))
worst <- c6_cases[which.max(c6_differences), ]
record_check("C6. uniform-correlations conditional grid = stats::integrate (same cases, k > 1)",
             max(c6_differences) < 2e-4,
             sprintf("largest difference %.1e (n = %g, k = %d, t = %g)", max(c6_differences),
                     worst$n, worst$k, worst$t_value))

## C7. Why eta = 1. Draw correlation matrices of (y, x1 ... xk) from LKJ(eta)
##     and compute the partial correlation of y with x1 given all the others
##     (from the inverse matrix). It should be Beta(eta, eta) on (-1, 1) for
##     every k, so flat for every k only when eta = 1.
##     Criterion: KS distance to Beta(eta, eta) below 0.03 (4,000 matrices), and
##     to the flat distribution above 0.05 when eta = 2.
set.seed(8)
c7_design <- expand.grid(k = c(1, 3, 8), eta = c(1, 2))
c7_partials <- parallel_map(seq_len(nrow(c7_design)), function(row) {
  replicate(4000, {
    precision <- solve(draw_lkj_correlation_matrix(c7_design$k[row] + 1, c7_design$eta[row]))
    -precision[1, 2] / sqrt(precision[1, 1] * precision[2, 2])
  })
})
c7_design$to_beta <- mapply(function(partial, eta)
  suppressWarnings(unname(ks.test((partial + 1) / 2, "pbeta", eta, eta)$statistic)),
  c7_partials, c7_design$eta)
c7_design$to_flat <- sapply(c7_partials, function(partial)
  suppressWarnings(unname(ks.test(partial, "punif", -1, 1)$statistic)))
record_check("C7. LKJ(eta): one partial correlation ~ Beta(eta, eta) for k = 1, 3, 8; flat only at eta = 1",
             all(c7_design$to_beta < 0.03) && all(c7_design$to_flat[c7_design$eta == 2] > 0.05),
             sprintf("KS to Beta(eta, eta): %s; to flat at eta = 2: %s",
                     paste(sprintf("%.3f", c7_design$to_beta), collapse = " "),
                     paste(sprintf("%.3f", c7_design$to_flat[c7_design$eta == 2]), collapse = " ")))

## C8. Against the exact partial-correlation Bayes factor of Kucharský,
##     Wagenmakers, van den Bergh & Ly (2026, TEST, eq. 2) with a flat prior
##     (stretched beta, alpha = 1; their gamma = delta = 0). With k - 1
##     controlling variables (our k counts the tested predictor too) and
##     h = (n - k) / 2 (their (n - controls - 1) / 2):
##       BF = B(1/2, alpha + h) / B(1/2, alpha) * 2F1(h, h; alpha + h + 1/2; r^2).
##     Their model is multivariate normal: the predictors are random, so the
##     tested predictor's residual SD (given the others) is an unknown
##     parameter. Ours conditions on X, so f uses the sample residual SD. The
##     only difference should be that: their f in sample units is the
##     population f times sqrt(W / n), W ~ chi^2(n - k) (the posterior of the
##     predictor's residual variance, with its 1/sigma prior). So:
##     (a) their formula reproduces their example (r = .51, n = 40, no
##         controls: BF10 = 45.88, section 4);
##     (b) ours differs from theirs: report the largest difference;
##     (c) ours with that prior on f (the flat partial correlation mixed over W)
##         equals theirs. Criteria: (a) within 0.01, (c) within 1e-5 on the log scale.
log_hypergeometric_2f1_series <- function(a, b, c, z, terms = 20000) {
  if (z == 0) return(0)
  j <- 0:(terms - 1)
  log_terms <- lgamma(a + j) - lgamma(a) + lgamma(b + j) - lgamma(b) -
    lgamma(c + j) + lgamma(c) - lgamma(j + 1) + j * log(z)
  largest <- max(log_terms)
  largest + log(sum(exp(log_terms - largest)))
}
log_bayes_factor_kucharsky_partial_correlation <- function(partial_r, sample_size,
                                                           number_of_controls, alpha = 1) {
  half_df <- (sample_size - number_of_controls - 1) / 2
  lbeta(0.5, alpha + half_df) - lbeta(0.5, alpha) +
    log_hypergeometric_2f1_series(half_df, half_df, alpha + half_df + 0.5, partial_r^2)
}
## P(|f| < f) in sample units when the partial correlation is flat in
## population units and the predictor's residual variance is unknown
prior_cdf_of_f_random_predictor <- function(f_values, sample_size, number_of_predictors) {
  df <- sample_size - number_of_predictors
  sapply(f_values, function(f) {
    if (f == 0) return(0)
    integrate(function(w) stats::dchisq(w, df) * f * sqrt(sample_size / w) /
                sqrt(1 + f^2 * sample_size / w),
              stats::qchisq(1e-14, df), stats::qchisq(1e-14, df, lower.tail = FALSE),
              rel.tol = 1e-12)$value
  })
}
c8_paper_example <- exp(log_bayes_factor_kucharsky_partial_correlation(0.51, 40, 0))
c8_cases <- expand.grid(partial_r = c(0, 0.1, 0.3, 0.5), n = c(20, 50, 200, 1000), k = c(1, 3, 8))
c8_results <- do.call(rbind, parallel_map(seq_len(nrow(c8_cases)), function(row) {
  partial_r <- c8_cases$partial_r[row]; n <- c8_cases$n[row]; k <- c8_cases$k[row]
  theirs <- log_bayes_factor_kucharsky_partial_correlation(partial_r, n, k - 1)
  c(ours_minus_theirs = log_bayes_factor_uniform_correlations_conditional(partial_r, n, k) - theirs,
    random_predictor_minus_theirs = log_bayes_factor_from_prior_cdf_of_f(partial_r, n, k,
      function(f_values) prior_cdf_of_f_random_predictor(f_values, n, k), points_f = 2000) - theirs)
}))
c8_cases <- cbind(c8_cases, c8_results)
worst <- c8_cases[which.max(abs(c8_cases$ours_minus_theirs)), ]
largest_at_n_1000 <- max(abs(c8_cases$ours_minus_theirs[c8_cases$n == 1000]))
record_check(paste("C8. Kucharský et al. (2026), flat prior: reproduces their BF 45.88; equals ours",
                   "with the predictor's residual SD unknown; differs from ours (fixed X)"),
             abs(c8_paper_example - 45.88) < 0.01 &&
               max(abs(c8_cases$random_predictor_minus_theirs)) < 1e-5,
             sprintf(paste("their example %.3f; ours with random predictor - theirs: largest %.1e;",
                           "ours - theirs: largest %.3f (n = %g, k = %d, r = %g), %.4f at n = 1000"),
                     c8_paper_example, max(abs(c8_cases$random_predictor_minus_theirs)),
                     worst$ours_minus_theirs, worst$n, worst$k, worst$partial_r, largest_at_n_1000))
cat("C8. ours - theirs, log BF (rows k and n, columns the partial correlation):\n")
c8_cases$k_and_n <- sprintf("k = %d, n = %4d", c8_cases$k, c8_cases$n)
print(round(xtabs(ours_minus_theirs ~ k_and_n + partial_r, c8_cases), 4))

## C9. The implied prior for a block of q columns (a term such as a factor with
##     q + 1 levels), given the other k - q predictors:
##     (a) uniform correlations, from LKJ(1) correlation matrices of
##         (y, x1 ... xk): the block's population partial R^2,
##         1 - (P_others)_yy / P_yy with P the inverse correlation matrix (and
##         P_others that of y and the others), is Beta(q/2, 1) for every k;
##     (b) uniform correlations, from g-prior coefficient draws with correlated
##         predictors: the block's f^2 = b' Sigma_(block | others) b gives
##         f^2 / (1 + f^2) ~ Beta(q/2, 1);
##     (c) flat-R, from the same draws: the block's share of the total f^2 is
##         Beta(q/2, (k - q)/2).
##     Criterion: KS distance below 0.03 for (a) (4,000 matrices, as C7) and below
##     0.01 for (b) and (c) (100,000 draws).
set.seed(9)
c9_lkj <- expand.grid(k = c(2, 3, 5, 8), q = c(2, 3, 4))
c9_lkj <- subset(c9_lkj, q <= k)
c9_lkj$distance <- unlist(parallel_map(seq_len(nrow(c9_lkj)), function(row) {
  k <- c9_lkj$k[row]; q <- c9_lkj$q[row]
  others <- setdiff(seq_len(k + 1), 1 + seq_len(q))      # y and the predictors not in the block
  partial_r_squared <- replicate(4000, {
    correlation <- draw_lkj_correlation_matrix(k + 1, 1)
    1 - solve(correlation[others, others])[1, 1] / solve(correlation)[1, 1]
  })
  suppressWarnings(unname(ks.test(partial_r_squared, "pbeta", q / 2, 1)$statistic))
}))
c9_draws <- expand.grid(k = c(4, 6, 10), q = c(2, 3))
c9_draws <- cbind(c9_draws, do.call(rbind, lapply(seq_len(nrow(c9_draws)), function(row) {
  k <- c9_draws$k[row]; q <- c9_draws$q[row]; block <- seq_len(q)
  covariance <- matrix(0.6, k, k); diag(covariance) <- 1
  block_given_others <- covariance[block, block] - covariance[block, -block, drop = FALSE] %*%
    solve(covariance[-block, -block], covariance[-block, block, drop = FALSE])
  block_f_squared <- function(b) rowSums((b[, block, drop = FALSE] %*% block_given_others) *
                                           b[, block, drop = FALSE])
  uniform <- draw_coefficients_from_prior(100000, covariance, "uniform_correlations")
  flat    <- draw_coefficients_from_prior(100000, covariance, "flat_r")
  c(uniform_correlations = suppressWarnings(unname(ks.test(
      block_f_squared(uniform) / (1 + block_f_squared(uniform)), "pbeta", q / 2, 1)$statistic)),
    flat_r_share = suppressWarnings(unname(ks.test(
      block_f_squared(flat) / rowSums((flat %*% covariance) * flat), "pbeta", q / 2, (k - q) / 2)$statistic)))
})))
record_check(paste("C9. block of q columns: uniform correlations gives partial R^2 ~ Beta(q/2, 1)",
                   "(LKJ and g-prior draws); flat-R gives share ~ Beta(q/2, (k-q)/2)"),
             all(c9_lkj$distance < 0.03) && all(c9_draws$uniform_correlations < 0.01) &&
               all(c9_draws$flat_r_share < 0.01),
             sprintf("largest KS distance: LKJ %.3f, g-prior draws %.4f, flat-R share %.4f",
                     max(c9_lkj$distance), max(c9_draws$uniform_correlations),
                     max(c9_draws$flat_r_share)))

## C10. The block test against adaptive integration, sharing no code with the
##      package (no Bessel function, no series): BF = integral of L(f) p(f) df,
##        L(f) = exp(-n f^2/2) E_theta[I(b cos theta)] / I(0),  b = sqrt(n R2) f,
##      theta the angle between the effect and the data's direction (density
##      proportional to sin^(q-2) theta on (0, pi)), m = n - k + q - 2, and
##        uniform correlations  p(f) = q f^(q-1) (1 + f^2)^(-q/2 - 1)
##        flat-R                p(f) = integral over s ~ Beta(q/2, (k-q)/2) of s (s + f^2)^(-3/2)
##      Criterion: within 2e-4 on the log scale.
##      (log_I_by_integration() is the one from C5.)
log_bayes_factor_block_by_integration <- function(partial_r_squared, n, k, q,
                                                  prior = c("flat-R", "uniform correlations")) {
  prior <- match.arg(prior)
  m <- n - k + q - 2
  log_I_0 <- log_I_by_integration(0, m)
  likelihood_ratio <- function(f) {
    b <- sqrt(n * partial_r_squared) * f
    if (b == 0) return(exp(-n * f^2 / 2))
    log_scale <- log_I_by_integration(b, m)
    direction_weight <- function(theta) sin(theta)^(q - 2) / beta(0.5, (q - 1) / 2)
    ## the integrand peaks at theta = 0 with width about 1/sqrt(b): split there
    pieces <- unique(pmin(pi, c(0, 2, 8, 30) / sqrt(b)))
    pieces <- c(pieces[pieces < pi], pi)
    average <- sum(sapply(seq_len(length(pieces) - 1), function(i)
      integrate(Vectorize(function(theta)
        direction_weight(theta) * exp(log_I_by_integration(b * cos(theta), m) - log_scale)),
        pieces[i], pieces[i + 1], rel.tol = 1e-10)$value))
    exp(-n * f^2 / 2 + log_scale - log_I_0) * average
  }
  prior_density <- function(f) {
    if (prior == "uniform correlations") return(q * f^(q - 1) * (1 + f^2)^(-q / 2 - 1))
    if (q == k) return((1 + f^2)^(-3 / 2))
    given_share <- function(s) stats::dbeta(s, q / 2, (k - q) / 2) * s * (s + f^2)^(-3 / 2)
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
c10_cases <- expand.grid(f_value = c(0, 1, 3, 6), q = c(2, 3, 5), k = c(5, 12),
                         n = c(50, 1e4, 1e6), prior = c("flat-R", "uniform correlations"),
                         stringsAsFactors = FALSE)
c10_cases <- subset(c10_cases, !(q == 5 & k == 5) | prior == "flat-R")
c10_differences <- unlist(parallel_map(seq_len(nrow(c10_cases)), function(row) {
  case <- c10_cases[row, ]
  ## F statistic F_value for the block: partial R^2 = q F / (q F + n - k - 1)
  partial_r_squared <- case$q * case$f_value / (case$q * case$f_value + case$n - case$k - 1)
  ours <- switch(case$prior,
    `flat-R` = log_bayes_factor_flat_r_conditional(sqrt(partial_r_squared), case$n, case$k,
                                                   number_of_columns = case$q),
    log_bayes_factor_uniform_correlations_conditional(sqrt(partial_r_squared), case$n, case$k,
                                                      number_of_columns = case$q))
  abs(ours - log_bayes_factor_block_by_integration(partial_r_squared, case$n, case$k, case$q,
                                                   case$prior))
}))
worst <- c10_cases[which.max(c10_differences), ]
record_check("C10. block conditional test (both priors) = stats::integrate (q 2 to 5, n 50 to 1e6, F 0 to 6)",
             max(c10_differences) < 2e-4,
             sprintf("largest difference %.1e (%s, n = %g, k = %d, q = %d, F = %g)",
                     max(c10_differences), worst$prior, worst$n, worst$k, worst$q, worst$f_value))

## C11. The block likelihood ratio from the data themselves: for simulated
##      data, integrate over sigma (prior 1/sigma) the likelihood of the
##      residualised outcome given the effect gamma (in noise units, whitened by
##      the block's sample covariance given the others), and average over the
##      direction of gamma (an angle, q = 2). This uses the data, not their
##      partial R^2, so it checks that the ratio depends on nothing else, and
##      the derivation (m = n - k + q - 2, the average over directions).
##      Criterion: within 1e-8 on the log scale.
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
  degrees <- n - (k - q) - 1                     # sigma^-degrees from the flat priors, times 1/sigma
  sigma_hat <- sqrt(sum(y_left^2) / degrees)
  log_marginal <- function(gamma) {              # log of the integral over sigma
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
record_check("C11. block likelihood ratio = integration over sigma and direction from the data (q = 2)",
             max(c11_differences) < 1e-8, sprintf("largest difference %.1e", max(c11_differences)))


report_checks()
