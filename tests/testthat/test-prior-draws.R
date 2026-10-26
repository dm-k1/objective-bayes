## The prior as a recipe gives R-squared ~ Beta(1/2, 1) for every k
## (docs/2_how_it_calculates.md, section 5)
test_that("implied R-squared is Beta(1/2, 1) for k = 1, 3, 10, correlated or not", {
  set.seed(1)
  correlated <- matrix(0.5, 3, 3); diag(correlated) <- 1
  for (covariance in list(matrix(1), correlated, diag(10))) {
    draws <- draw_coefficients_from_flat_r_prior(50000, covariance)
    f_squared <- rowSums((draws %*% covariance) * draws)
    r_squared <- f_squared / (1 + f_squared)
    distance <- suppressWarnings(ks.test(r_squared, "pbeta", 0.5, 1)$statistic)
    expect_lt(distance, 0.01)
  }
})

test_that("with one predictor the correlation is flat from -1 to 1", {
  set.seed(2)
  effect <- draw_coefficients_from_flat_r_prior(50000, matrix(1))[, 1]
  correlation <- effect / sqrt(1 + effect^2)
  expect_lt(suppressWarnings(ks.test(correlation, "punif", -1, 1)$statistic), 0.01)
})
