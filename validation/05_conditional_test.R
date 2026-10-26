## =====================================================================
##  Validation 05: checks of the per-coefficient (conditional) tests
##
##  The conditional test looks only at one coefficient's sample partial
##  correlation, with the prior that the joint prior implies for one
##  coefficient (R/conditional.R; docs/2_how_it_calculates.md, section 9):
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
joint_draws <- draw_coefficients_from_flat_r_prior(200000, correlated)
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

## C6. The uniform-correlations conditional test (validation/comparison_methods.R)
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

report_checks()
