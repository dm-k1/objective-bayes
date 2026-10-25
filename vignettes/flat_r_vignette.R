## =====================================================================
##  Vignette: the flat-R prior in practice
##  Fit with brms, then Bayes factors for any model comparison.
##  Run with the project folder as working directory.
## =====================================================================
source("R/flat_r.R")
library(brms)   # for summary(), bayes_R2(), pp_check(), plot() on the fitted models
options(mc.cores = max(1, parallel::detectCores() - 1))   # chains in parallel, one core left free

## ---- 1. Data ---------------------------------------------------------
## Three correlated predictors (x1 and x2 matter, x3 doesn't) and one
## unrelated extra predictor x4 for a non-nested comparison.
set.seed(42)
sample_size <- 60
predictor_correlations <- matrix(0.4, 3, 3); diag(predictor_correlations) <- 1
predictor_values <- matrix(rnorm(sample_size * 3), sample_size) %*% chol(predictor_correlations)
my_data <- data.frame(x1 = predictor_values[, 1],
                      x2 = predictor_values[, 2],
                      x3 = predictor_values[, 3])
my_data$y  <- 0.3 * my_data$x1 + 0.2 * my_data$x2 + rnorm(sample_size)
my_data$x4 <- rnorm(sample_size)

## ---- 2. Estimation: fit like any brms model -------------------------
full_model <- fit_flat_r_model(y ~ x1 + x2 + x3, my_data,
                                       chains = 4, iter = 2000, seed = 1, refresh = 0)
summary(full_model)              # estimates, intervals, Rhat, ESS
bayes_R2(full_model)             # explained share
pp_check(full_model)             # posterior predictive check
plot(full_model)                 # trace plots and posteriors

## For comparison: the ordinary confidence intervals
confint(lm(y ~ x1 + x2 + x3, my_data))

## ---- 3. Bayes factor per predictor ----------------------------------
bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, my_data)

## ---- 4. Comparing any two models ------------------------------------
without_x1 <- fit_flat_r_model(y ~ x2 + x3, my_data, chains = 4, seed = 1, refresh = 0)
other_set  <- fit_flat_r_model(y ~ x4,      my_data, chains = 4, seed = 1, refresh = 0)

bayes_factor_flat_r(full_model, without_x1)   # nested: evidence for x1
bayes_factor_flat_r(full_model, other_set)    # non-nested comparison

## Works on plain lm() fits too, without running Stan:
bayes_factor_flat_r(lm(y ~ x1 + x2 + x3, my_data), lm(y ~ 1, my_data))

## ---- 5. What a null Bayes factor rests on ---------------------------
## Which effect sizes can the data not tell apart from zero? (x2 here)
support_x2 <- support_per_partial_correlation(lm(y ~ x1 + x2 + x3, my_data), "x2")
plot(support_x2, type = "l", log = "y", lwd = 2,
     xlab = "true partial correlation", ylab = "support compared with zero (log scale)",
     main = "x2: which effect sizes can't be told apart from zero?")
abline(h = 1, lty = 3)
rect(-2, 1/3, 2, 3, col = rgb(0, 0, 0, 0.1), border = NA)   # within a factor 3 of zero
support_regions(support_x2)   # regions: better / can't tell / worse than zero

## ---- 6. Comparison with the default (JZS / JASP) --------------------
## Optional; needs the BayesFactor package
if (requireNamespace("BayesFactor", quietly = TRUE)) {
  default_bf_x1 <- BayesFactor::lmBF(y ~ x1 + x2 + x3, my_data) /
                   BayesFactor::lmBF(y ~ x2 + x3, my_data)
  print(default_bf_x1)
}

## The closed form is validated against brms bridge sampling in
## validation/01_checks.R (18.8 vs 18.7 for x1; 2.83 vs 2.82 for a factor), and
## against exact integration: accurate to 4e-05 on the log scale.
