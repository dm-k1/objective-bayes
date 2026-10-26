test_that("fit_flat_r_model() runs and works with print, loo and pp_check", {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE) || requireNamespace("rstan", quietly = TRUE),
              "needs cmdstanr or rstan")
  set.seed(1)
  data <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  data$y <- 0.4 * data$x1 + rnorm(50)
  fit <- fit_flat_r_model(y ~ x1 + x2, data, chains = 2, iter_warmup = 500,
                          iter_sampling = 500, seed = 1)
  expect_s3_class(fit, "flat_r_fit")
  expect_setequal(posterior::variables(fit$draws), c("Intercept", "x1", "x2", "sigma", "R_squared"))
  expect_output(print(fit), "flat-R prior")
  expect_s3_class(loo(fit), "loo")
  expect_s3_class(pp_check(fit), "ggplot")
  expect_error(fit_flat_r_model(y ~ 0 + x1, data), "intercept")
  uniform <- fit_flat_r_model(y ~ x1 + x2, data, chains = 2, iter_warmup = 500,
                              iter_sampling = 500, seed = 1, prior = "uniform_correlations")
  expect_output(print(uniform), "uniform-correlations prior")
})
