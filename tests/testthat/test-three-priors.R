## The three priors as one family, R-squared ~ Beta(a, 1) with a uniform
## direction (R/priors.R): fast versions of the family checks in validation/02
## and validation/05

test_that("each prior has its shape, and flat_partial_r is the default", {
  expect_equal(prior_shape("flat_partial_r", 6), 3)
  expect_equal(prior_shape("flat_r", 6), 1 / 2)
  expect_equal(prior_shape("flat_r_squared", 6), 1)
  expect_error(prior_shape(NULL, 3), "Choose a prior")
  expect_error(log_bayes_factor_vs_intercept_only(0.2, 50, 3, prior = "flat"), "Unknown prior")
  ## without a prior, the default is used and announced once per session
  default_prior_announced$done <- NULL
  expect_message(by_default <- log_bayes_factor_vs_intercept_only(0.2, 50, 3), "default prior")
  expect_equal(by_default, log_bayes_factor_vs_intercept_only(0.2, 50, 3, prior = "flat_partial_r"))
  expect_silent(log_bayes_factor_conditional(0.2, 50, 3))
  expect_equal(log_bayes_factor_conditional(0.2, 50, 3),
               log_bayes_factor_conditional(0.2, 50, 3, prior = "flat_partial_r"))
})

test_that("rows dropped for missing values are reported", {
  with_missing <- mtcars
  with_missing$hp[c(2, 5)] <- NA
  expect_warning(bayes_factor_per_term(mpg ~ wt + hp, with_missing, prior = "flat_partial_r"),
                 "2 of 32 rows")
  expect_warning(bayes_factor_models(lm(mpg ~ wt + hp, with_missing),
                                     lm(mpg ~ wt, with_missing[-c(2, 5), ]),
                                     prior = "flat_partial_r"), "2 rows with a missing value")
  expect_silent(bayes_factor_per_term(mpg ~ wt + hp, mtcars, prior = "flat_partial_r"))
})

test_that("the block route reproduces the closed forms of flat R and flat partial r", {
  for (case in list(c(0.1, 50, 1), c(0.3, 50, 3), c(0.05, 200, 5), c(0.5, 100, 8),
                    c(0.02, 1000, 10), c(0.9, 30, 4))) {
    r_squared <- case[1]; n <- case[2]; k <- case[3]
    expect_equal(log_bayes_factor_shape_vs_intercept_only(r_squared, n, k, shape = 1 / 2),
                 log_bayes_factor_flat_r_vs_intercept_only(r_squared, n, k), tolerance = 1e-4)
    expect_equal(log_bayes_factor_shape_vs_intercept_only(r_squared, n, k, shape = k / 2),
                 log_bayes_factor_flat_partial_r_vs_intercept_only(r_squared, n, k),
                 tolerance = 2e-4)
  }
})

test_that("flat R-squared: the whole-model Bayes factor agrees with adaptive integration", {
  ## the block likelihood ratio integrated against the prior density of |f|,
  ## 2 f / (1 + f^2)^2 (R^2 = f^2 / (1 + f^2) flat), by integrate()
  by_integration <- function(r_squared, n, k) {
    integrand <- function(f) exp(log_likelihood_ratio_block(f, r_squared, n, k, k)) *
      2 * f / (1 + f^2)^2
    peak <- sqrt(r_squared / (1 - r_squared)); width <- 1 / sqrt(n)
    breaks <- sort(unique(c(0, max(0, peak - 20 * width), peak, peak + 20 * width, Inf)))
    log(sum(sapply(seq_len(length(breaks) - 1), function(i)
      integrate(integrand, breaks[i], breaks[i + 1], rel.tol = 1e-10, subdivisions = 2000)$value)))
  }
  for (case in list(c(0.1, 50, 1), c(0.3, 50, 3), c(0.001, 500, 3), c(0.05, 200, 5)))
    expect_equal(log_bayes_factor_vs_intercept_only(case[1], case[2], case[3], prior = "flat_r_squared"),
                 by_integration(case[1], case[2], case[3]), tolerance = 1e-4)
})

test_that("flat R-squared: the closed form (k >= 3) agrees with the block route", {
  for (case in list(c(0.3, 50, 3), c(0.05, 200, 5), c(0.5, 100, 8), c(0.9, 30, 4)))
    expect_equal(log_bayes_factor_vs_intercept_only(case[1], case[2], case[3], prior = "flat_r_squared"),
                 log_bayes_factor_shape_vs_intercept_only(case[1], case[2], case[3], shape = 1),
                 tolerance = 2e-3)
})

test_that("flat R-squared: predictive matching and information consistency", {
  for (k in 1:5)
    expect_equal(log_bayes_factor_vs_intercept_only(1, k + 1, k, prior = "flat_r_squared"), 0)
  evidence <- sapply(c(0.9, 0.99, 0.999, 0.9999), log_bayes_factor_vs_intercept_only,
                     sample_size = 20, number_of_predictors = 3, prior = "flat_r_squared")
  expect_true(all(diff(evidence) > 0))
  expect_gt(evidence[4], 50)
})

test_that("conditional test: the share mixture with a = k/2 is a flat partial correlation", {
  f_values <- c(0.01, 0.1, 0.3, 1, 3)
  for (k in c(2, 3, 5, 10))
    expect_equal(prior_cdf_of_f(f_values, k, points_share = 500, shape = k / 2),
                 f_values / sqrt(1 + f_values^2), tolerance = 1e-4)
  ## the same for a block of q columns: partial R^2 ~ Beta(q/2, 1)
  expect_equal(prior_cdf_of_f(f_values, 6, points_share = 500, number_of_columns = 2, shape = 3),
               f_values^2 / (1 + f_values^2), tolerance = 1e-4)
})

test_that("conditional test: one predictor equals the joint test, for every prior", {
  for (prior in c("flat_partial_r", "flat_r", "flat_r_squared"))
    for (case in list(c(0.3, 30), c(0.1, 100)))
      expect_equal(log_bayes_factor_conditional(case[1], case[2], 1, prior = prior),
                   log_bayes_factor_vs_intercept_only(case[1]^2, case[2], 1, prior = prior),
                   tolerance = 1e-3)
})

test_that("conditional test: flat R and flat R-squared narrow with k, flat partial r does not", {
  ## at a fixed partial correlation, how far the log Bayes factor moves from
  ## k = 1 to k = 10
  spread <- function(prior)
    diff(range(sapply(c(1, 3, 10), function(k) log_bayes_factor_conditional(0.2, 200, k, prior))))
  expect_lt(spread("flat_partial_r"), 0.25)
  expect_gt(spread("flat_r"), spread("flat_partial_r"))
  expect_gt(spread("flat_r_squared"), spread("flat_partial_r"))
})

test_that("the block route holds at very strong evidence (large-argument Bessel function)", {
  for (case in list(c(0.99, 10000, 3), c(0.99, 10000, 30), c(0.999, 10000, 10)))
    expect_equal(log_bayes_factor_shape_vs_intercept_only(case[1], case[2], case[3], shape = 1 / 2),
                 log_bayes_factor_flat_r_vs_intercept_only(case[1], case[2], case[3]), tolerance = 1e-6)
  expect_true(is.finite(log_bayes_factor_vs_intercept_only(0.99, 10000, 30, prior = "flat_r_squared")))
})
