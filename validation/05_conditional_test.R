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
##
##  Run from the project folder: Rscript validation/05_conditional_test.R
##  Under a minute.
## =====================================================================
source("validation/settings.R")   # record_check()
pkgload::load_all(".", quiet = TRUE)
source("validation/comparison_methods.R")


## C1. With one predictor the test is the joint flat-R test (both are flat on r).
c1_differences <- sapply(list(c(0.3, 30), c(0.1, 100), c(0.5, 50), c(0.05, 500)), function(s)
  abs(log_bayes_factor_flat_r_conditional(s[1], s[2], 1, points_R = 4000) -
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
  abs(log_bayes_factor_uniform_correlations_conditional(partial_correlation_of(fit, 2), n, 3,
                                                        points_R = 4000) -
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

report_checks()
