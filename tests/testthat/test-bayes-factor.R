## Fast versions of validation/02_bayes_factor.R

example_data <- local({
  set.seed(42)
  correlations <- matrix(0.4, 3, 3); diag(correlations) <- 1
  predictors <- matrix(rnorm(60 * 3), 60) %*% chol(correlations)
  data <- data.frame(x1 = predictors[, 1], x2 = predictors[, 2], x3 = predictors[, 3])
  data$y <- 0.3 * data$x1 + 0.2 * data$x2 + rnorm(60)
  data
})

test_that("the vignette's Bayes factors are unchanged (joint and conditional)", {
  joint <- bayes_factor_per_term(y ~ x1 + x2 + x3, example_data, test = "joint", prior = "flat_r")
  expect_equal(joint$BF10, c(18.696, 0.3827, 1.5872), tolerance = 1e-3)
  conditional <- bayes_factor_per_term(y ~ x1 + x2 + x3, example_data, prior = "flat_r")
  expect_equal(conditional$BF10, c(17.489, 0.3911, 1.1833), tolerance = 1e-3)
  expect_equal(unique(conditional$test), "conditional")
})

test_that("conditional test: one predictor equals the joint test", {
  for (case in list(c(0.3, 30), c(0.1, 100), c(0.5, 50)))
    expect_equal(log_bayes_factor_flat_r_conditional(case[1], case[2], 1),
                 log_bayes_factor_flat_r_vs_intercept_only(case[1]^2, case[2], 1), tolerance = 1e-3)
})

test_that("conditional test: the grid agrees with adaptive integration, also for large n", {
  ## k = 3: share ~ Beta(1/2, 1) gives the prior of |f| the density
  ## asinh(1/f) - 1/sqrt(1 + f^2). The likelihood ratio by integrate() too:
  ## exp(-n f^2/2) (I(b) + I(-b)) / (2 I(0)), each I(b) scaled by its peak.
  log_I <- function(b, m) {
    peak     <- (b + sqrt(b^2 + 4 * m)) / 2
    log_peak <- m * log(peak) - peak^2 / 2 + b * peak
    log_peak + log(integrate(function(v) exp(m * log(v) - v^2 / 2 + b * v - log_peak),
                             max(0, peak - 40), peak + 40, rel.tol = 1e-10)$value)
  }
  by_integration <- function(partial_r, n) {
    m <- n - 4
    likelihood_ratio_times_prior <- Vectorize(function(f) {
      b <- sqrt(n) * partial_r * f
      mean(exp(-n * f^2 / 2 - log_I(0, m) + c(log_I(b, m), log_I(-b, m)))) *
        (asinh(1 / f) - 1 / sqrt(1 + f^2))
    })
    peak  <- sqrt(m / n) * abs(partial_r) / sqrt(1 - partial_r^2)
    width <- 1 / sqrt(n * (1 - partial_r^2))
    breaks <- unique(pmax(0, peak + c(-15, -5, -2, 0, 2, 5, 15) * width))
    log(sum(sapply(seq_len(length(breaks) - 1), function(i)
      integrate(likelihood_ratio_times_prior, breaks[i], breaks[i + 1], rel.tol = 1e-8)$value)))
  }
  for (case in list(c(3, 200), c(0, 1e4), c(2, 1e6))) {     # t-value, n
    partial_r <- case[1] / sqrt(case[1]^2 + case[2] - 4)
    expect_lt(abs(log_bayes_factor_flat_r_conditional(partial_r, case[2], 3) -
                  by_integration(partial_r, case[2])), 2e-4)
  }
})

test_that("the joint test depends on what the others explain; the conditional test does not", {
  ## same evidence about x1 (partial correlation 0.2, n = 100, k = 3); the
  ## other predictors explain 5% or 60% (docs/2_how_it_calculates.md, section 10)
  joint_bf <- function(r2_without) {
    r2_with <- 1 - (1 - r2_without) * (1 - 0.2^2)
    exp(log_bayes_factor_flat_r_vs_intercept_only(r2_with, 100, 3) -
        log_bayes_factor_flat_r_vs_intercept_only(r2_without, 100, 2))
  }
  expect_gt(joint_bf(0.05), 2)     # leans toward an effect
  expect_lt(joint_bf(0.60), 1)     # leans toward no effect, with the same x1 evidence
  expect_equal(exp(log_bayes_factor_flat_r_conditional(0.2, 100, 3)), 1.37, tolerance = 0.01)
})

test_that("multi-column terms fall back to the joint test", {
  data <- transform(example_data, group = factor(rep(c("a", "b", "c"), 20)))
  result <- bayes_factor_per_term(y ~ x1 + group, data)
  expect_equal(result$test, c("conditional", "joint"))
})

test_that("an aliased column stops with a clear message", {
  aliased <- transform(example_data, x4 = x1 + x2)
  expect_error(bayes_factor_per_term(y ~ x1 + x2 + x4, aliased), "aliased): x4")
  expect_error(bayes_factor_models(lm(y ~ x1 + x2 + x4, aliased), lm(y ~ 1, aliased)),
               "linearly dependent")
})

test_that("a term with missing values: the model without it uses the same rows", {
  with_missing <- transform(example_data, x3 = replace(x3, c(3, 8), NA))
  complete_rows <- with_missing[-c(3, 8), ]
  expect_equal(bayes_factor_per_term(y ~ x1 + x2 + x3, with_missing, test = "joint"),
               bayes_factor_per_term(y ~ x1 + x2 + x3, complete_rows, test = "joint"))
})

test_that("rescaling or recoding predictors leaves the Bayes factor unchanged", {
  original <- bayes_factor_models(lm(y ~ x1 + x2 + x3, example_data), lm(y ~ 1, example_data))
  rescaled <- transform(example_data, x1 = 100 * x1 + 7, y = 1000 * y - 3)
  recoded  <- transform(example_data, sum_12 = x1 + x2, difference_12 = x1 - x2)
  expect_equal(bayes_factor_models(lm(y ~ x1 + x2 + x3, rescaled), lm(y ~ 1, rescaled)), original)
  expect_equal(bayes_factor_models(lm(y ~ sum_12 + difference_12 + x3, recoded),
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
  expect_error(bayes_factor_models(lm(y ~ x1, example_data), lm(y ~ x1, example_data[-1, ])),
               "same outcome")
})

## Uniform correlations (R-squared ~ Beta(k/2, 1)): fast versions of validation/02, part C

test_that("uniform correlations: the grid agrees with adaptive integration", {
  r_squared <- 0.2; n <- 50; k <- 3
  by_integration <- integrate(function(g)
    (1 + g)^((n - 1 - k) / 2) * (1 + g * (1 - r_squared))^(-(n - 1) / 2) *
      (n / 2) * g^(-2) * exp(-n / (2 * g)), 0, Inf)$value
  expect_equal(log_bayes_factor_uniform_correlations_vs_intercept_only(r_squared, n, k),
               log(by_integration), tolerance = 1e-4)
})

test_that("with one predictor the two priors are the same", {
  expect_equal(log_bayes_factor_uniform_correlations_vs_intercept_only(0.09, 40, 1),
               log_bayes_factor_flat_r_vs_intercept_only(0.09, 40, 1), tolerance = 1e-4)
  expect_equal(log_bayes_factor_uniform_correlations_conditional(0.3, 40, 1),
               log_bayes_factor_flat_r_conditional(0.3, 40, 1), tolerance = 1e-4)
})

test_that("uniform correlations: predictive matching and invariance", {
  for (k in 1:5)
    expect_equal(exp(log_bayes_factor_uniform_correlations_vs_intercept_only(1 - 1e-12, k + 1, k)), 1,
                 tolerance = 1e-3)
  original <- bayes_factor_models(lm(y ~ x1 + x2 + x3, example_data), lm(y ~ 1, example_data),
                                  prior = "uniform_correlations")
  recoded  <- transform(example_data, sum_12 = x1 + x2, difference_12 = x1 - x2)
  expect_equal(bayes_factor_models(lm(y ~ sum_12 + difference_12 + x3, recoded), lm(y ~ 1, recoded),
                                   prior = "uniform_correlations"), original)
})

test_that("uniform correlations: the conditional test does not depend on k (flat partial r for every k)", {
  ## the prior is the same for every k; only the residual degrees of freedom
  ## n - k - 1 change, so at a fixed partial correlation the log Bayes factor
  ## moves little with k (1.57, 1.53, 1.39 for k = 1, 3, 10)
  by_k <- sapply(c(1, 3, 10), function(k) log_bayes_factor_uniform_correlations_conditional(0.2, 200, k))
  expect_lt(diff(range(by_k)), 0.25)
  ## flat-R's per-coefficient prior narrows with k, so its Bayes factor moves
  ## more (1.57, 1.94, 1.90)
  flat_r_by_k <- sapply(c(1, 3, 10), function(k) log_bayes_factor_flat_r_conditional(0.2, 200, k))
  expect_gt(diff(range(flat_r_by_k)), diff(range(by_k)))
})

test_that("per-term table works with both priors", {
  uniform <- bayes_factor_per_term(y ~ x1 + x2 + x3, example_data, prior = "uniform_correlations")
  expect_equal(unique(uniform$test), "conditional")
  expect_true(all(uniform$BF10 > 0))
  expect_false(isTRUE(all.equal(uniform$BF10, bayes_factor_per_term(y ~ x1 + x2 + x3, example_data, prior = "flat_r")$BF10)))
})
