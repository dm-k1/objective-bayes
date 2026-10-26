## Fast versions of validation/02_bayes_factor.R

example_data <- local({
  set.seed(42)
  correlations <- matrix(0.4, 3, 3); diag(correlations) <- 1
  predictors <- matrix(rnorm(60 * 3), 60) %*% chol(correlations)
  data <- data.frame(x1 = predictors[, 1], x2 = predictors[, 2], x3 = predictors[, 3])
  data$y <- 0.3 * data$x1 + 0.2 * data$x2 + rnorm(60)
  data
})

test_that("the vignette's Bayes factors are unchanged", {
  per_term <- bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, example_data)
  expect_equal(per_term$BF10, c(18.696, 0.3827, 1.5872), tolerance = 1e-3)
})

test_that("rescaling or recoding predictors leaves the Bayes factor unchanged", {
  original <- bayes_factor_flat_r(lm(y ~ x1 + x2 + x3, example_data), lm(y ~ 1, example_data))
  rescaled <- transform(example_data, x1 = 100 * x1 + 7, y = 1000 * y - 3)
  recoded  <- transform(example_data, sum_12 = x1 + x2, difference_12 = x1 - x2)
  expect_equal(bayes_factor_flat_r(lm(y ~ x1 + x2 + x3, rescaled), lm(y ~ 1, rescaled)), original)
  expect_equal(bayes_factor_flat_r(lm(y ~ sum_12 + difference_12 + x3, recoded),
                                   lm(y ~ 1, recoded)), original)
})

test_that("the grid agrees with adaptive integration (one predictor)", {
  ## one predictor: the g-prior Bayes factor averaged over g ~ Inverse-Gamma(1, n/2)
  r_squared <- 0.09; n <- 40
  by_integration <- integrate(function(g)
    (1 + g)^((n - 2) / 2) * (1 + g * (1 - r_squared))^(-(n - 1) / 2) *
      (n / 2) * g^(-2) * exp(-n / (2 * g)), 0, Inf)$value
  expect_equal(log_bayes_factor_flat_r_vs_intercept_only(r_squared, n, 1), log(by_integration),
               tolerance = 1e-4)
})

test_that("predictive matching: with n = k + 1 the Bayes factor is 1", {
  for (k in 1:5)
    expect_equal(exp(log_bayes_factor_flat_r_vs_intercept_only(1 - 1e-12, k + 1, k)), 1,
                 tolerance = 1e-3)
})

test_that("information consistency: evidence grows without limit as R-squared goes to 1", {
  evidence <- sapply(c(0.9, 0.99, 0.999, 0.9999), log_bayes_factor_flat_r_vs_intercept_only,
                     sample_size = 20, number_of_predictors = 3)
  expect_true(all(diff(evidence) > 0))
  expect_gt(evidence[4], 50)
})

test_that("models must share outcome and rows", {
  expect_error(bayes_factor_flat_r(lm(y ~ x1, example_data), lm(y ~ x1, example_data[-1, ])),
               "same outcome")
})
