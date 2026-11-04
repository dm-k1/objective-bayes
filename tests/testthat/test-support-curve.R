## The support curve is the likelihood ratio of the conditional test, so the
## conditional Bayes factor is the curve averaged over the prior.

support_data <- local({
  set.seed(7)
  data <- data.frame(x1 = rnorm(80), x2 = rnorm(80), group = factor(rep(c("a", "b", "c", "d"), 20)))
  data$y <- 0.3 * data$x1 + 0.1 * data$x2 + 0.3 * (data$group == "c") + rnorm(80)
  data
})
support_fit <- lm(y ~ x1 + x2 + group, support_data)

test_that("one column: the curve averaged over a flat partial correlation is the conditional BF", {
  ## flat partial r: the partial correlation is flat on (-1, 1)
  curve <- function(rho) support_per_partial_correlation(support_fit, "x2", rho)$support_vs_zero
  averaged <- integrate(curve, -1 + 1e-9, 1 - 1e-9, rel.tol = 1e-8)$value / 2
  expect_equal(averaged, bayes_factor_per_term(y ~ x1 + x2 + group, support_data, prior = "flat_partial_r")$BF10[2],
               tolerance = 1e-4)
})

test_that("block: the curve averaged over Beta(q/2, 1) is the block conditional BF", {
  curve <- function(r2) support_per_partial_r_squared(support_fit, "group", r2)$support_vs_zero *
    dbeta(r2, 3 / 2, 1)
  averaged <- integrate(curve, 0, 1 - 1e-9, rel.tol = 1e-8)$value
  expect_equal(averaged, bayes_factor_per_term(y ~ x1 + x2 + group, support_data, prior = "flat_partial_r")$BF10[3],
               tolerance = 1e-4)
})

test_that("a one-column term: the block curve is the signed curve averaged over both signs", {
  rho <- c(0.05, 0.2, 0.4)
  signed <- support_per_partial_correlation(support_fit, "x1", c(-rho, rho))$support_vs_zero
  block  <- support_per_partial_r_squared(support_fit, "x1", rho^2)$support_vs_zero
  expect_equal(block, (signed[1:3] + signed[4:6]) / 2, tolerance = 1e-10)
})

test_that("the curve is 1 at zero, signed, and accepts a term label or a design column", {
  support <- support_per_partial_correlation(support_fit, "x1", c(-0.3, 0, 0.3))
  expect_equal(support$support_vs_zero[2], 1)
  expect_gt(support$support_vs_zero[3], support$support_vs_zero[1])   # x1's effect is positive
  expect_equal(support_per_partial_r_squared(support_fit, "group"),
               support_per_partial_r_squared(support_fit, c("groupb", "groupc", "groupd")))
  expect_error(support_per_partial_correlation(support_fit, "group"), "has 3 columns")
  expect_error(support_per_partial_correlation(support_fit, "x9"), "not a term")
})

test_that("at large n it agrees with the plug-in approximation (known noise, t-density)", {
  set.seed(8)
  n <- 20000
  data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  data$y <- 0.02 * data$x1 + rnorm(n)
  fit <- lm(y ~ x1 + x2, data)
  rho <- seq(-0.01, 0.04, by = 0.005)
  slope_per_rho <- sd(qr.resid(qr(cbind(1, data$x2)), data$x1)) / summary(fit)$sigma
  coefficient <- summary(fit)$coefficients["x1", ]
  slopes <- rho / (slope_per_rho * sqrt(1 - rho^2))
  plug_in <- dt((coefficient[1] - slopes) / coefficient[2], fit$df.residual) /
             dt(coefficient[1] / coefficient[2], fit$df.residual)
  expect_equal(support_per_partial_correlation(fit, "x1", rho)$support_vs_zero, unname(plug_in),
               tolerance = 0.01)
})

test_that("support_regions reads either table", {
  regions <- support_regions(support_per_partial_r_squared(support_fit, "group"))
  expect_equal(regions$from[1], 0)
  expect_true(all(regions$verdict %in% c("explains the data better than zero",
                                         "explains the data worse than zero",
                                         "can't be told apart from zero")))
})
