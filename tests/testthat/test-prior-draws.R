## The prior as a recipe gives R-squared ~ Beta(a, 1) (docs/2_how_it_calculates.md, section 5)
test_that("flat R: implied R-squared is Beta(1/2, 1) for k = 1, 3, 10, correlated or not", {
  set.seed(1)
  correlated <- matrix(0.5, 3, 3); diag(correlated) <- 1
  for (covariance in list(matrix(1), correlated, diag(10))) {
    draws <- draw_coefficients_from_prior(50000, covariance, prior = "flat_r")
    f_squared <- rowSums((draws %*% covariance) * draws)
    r_squared <- f_squared / (1 + f_squared)
    distance <- suppressWarnings(ks.test(r_squared, "pbeta", 0.5, 1)$statistic)
    expect_lt(distance, 0.01)
  }
})

test_that("with one predictor the correlation is flat from -1 to 1", {
  set.seed(2)
  effect <- draw_coefficients_from_prior(50000, matrix(1), prior = "flat_r")[, 1]
  correlation <- effect / sqrt(1 + effect^2)
  expect_lt(suppressWarnings(ks.test(correlation, "punif", -1, 1)$statistic), 0.01)
})

test_that("flat partial r: implied R-squared is Beta(k/2, 1), correlated or not", {
  set.seed(3)
  correlated <- matrix(0.5, 3, 3); diag(correlated) <- 1
  for (covariance in list(matrix(1), correlated, diag(10))) {
    draws <- draw_coefficients_from_prior(50000, covariance, prior = "flat_partial_r")
    f_squared <- rowSums((draws %*% covariance) * draws)
    distance <- suppressWarnings(ks.test(f_squared / (1 + f_squared), "pbeta",
                                         ncol(covariance) / 2, 1)$statistic)
    expect_lt(distance, 0.01)
  }
})

test_that("flat partial r: every coefficient's partial correlation is flat, for any k", {
  set.seed(4)
  correlated <- matrix(0.5, 5, 5); diag(correlated) <- 1
  draws <- draw_coefficients_from_prior(50000, correlated, prior = "flat_partial_r")
  partial_f <- draws[, 1] * sqrt(1 / solve(correlated)[1, 1])   # effect / spread left after the others
  expect_lt(suppressWarnings(ks.test(partial_f / sqrt(1 + partial_f^2), "punif", -1, 1)$statistic), 0.01)
})

test_that("flat R-squared: implied R-squared is flat for every k, correlated or not", {
  set.seed(5)
  correlated <- matrix(0.5, 3, 3); diag(correlated) <- 1
  for (covariance in list(matrix(1), correlated, diag(10))) {
    draws <- draw_coefficients_from_prior(50000, covariance, prior = "flat_r_squared")
    f_squared <- rowSums((draws %*% covariance) * draws)
    expect_lt(suppressWarnings(ks.test(f_squared / (1 + f_squared), "punif")$statistic), 0.01)
  }
})

test_that("the default is flat_partial_r, and an unknown prior is refused", {
  set.seed(1); by_default <- suppressMessages(draw_coefficients_from_prior(10, diag(2)))
  set.seed(1); expect_equal(by_default, draw_coefficients_from_prior(10, diag(2), prior = "flat_partial_r"))
  expect_error(draw_coefficients_from_prior(10, diag(2), prior = "uniform_correlations"),
               "Unknown prior")
})
