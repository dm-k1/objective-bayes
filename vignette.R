## =====================================================================
##  Vignette: the flat-R prior on one dataset, start to finish
##  Estimation with brms, then Bayes factors for any model comparison.
##  Run with the project folder as working directory.
## =====================================================================
source("R/flat_r.R")
library(brms)      # summary(), bayes_R2(), pp_check(), plot() on the fitted models
library(ggplot2)
options(mc.cores = max(1, parallel::detectCores() - 1))   # chains in parallel, one core left free

## ---- 1. Data ---------------------------------------------------------
## Three correlated predictors (x1 and x2 matter, x3 doesn't) and one
## unrelated extra predictor x4 for a non-nested comparison.
set.seed(42)
sample_size <- 60
predictor_correlations <- matrix(0.4, nrow = 3, ncol = 3); diag(predictor_correlations) <- 1
predictor_values <- matrix(rnorm(sample_size * 3), sample_size) %*% chol(predictor_correlations)
my_data <- data.frame(x1 = predictor_values[, 1],
                      x2 = predictor_values[, 2],
                      x3 = predictor_values[, 3])
my_data$y  <- 0.3 * my_data$x1 + 0.2 * my_data$x2 + rnorm(sample_size)
my_data$x4 <- rnorm(sample_size)

## ---- 2. Estimation: fit like any brms model -------------------------
sampler_settings <- list(chains = 4, iter = 2000, warmup = 1000, seed = 1, refresh = 0)
fit_with <- function(formula)
  do.call(fit_flat_r_model, c(list(formula = formula, data = my_data), sampler_settings))

full_model <- fit_with(y ~ x1 + x2 + x3)
summary(full_model)              # estimates, intervals, Rhat, ESS
bayes_R2(full_model)             # explained share
pp_check(full_model)             # posterior predictive check
plot(full_model)                 # trace plots and posteriors

## Next to the ordinary estimates: at n = 60 with three predictors the prior
## pulls the estimates toward zero (docs/1_from_behaviour.md, step 4)
cbind(flat_r   = fixef(full_model)[c("x1", "x2", "x3"), "Estimate"],
      ordinary = coef(lm(y ~ x1 + x2 + x3, data = my_data))[c("x1", "x2", "x3")])

## ---- 3. Bayes factor per predictor ----------------------------------
## Closed form: no Stan needed, the data enter only through R-squared, n and k
bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, my_data)

## ---- 4. Comparing any two models ------------------------------------
without_x1 <- fit_with(y ~ x2 + x3)
other_set  <- fit_with(y ~ x4)

bayes_factor_flat_r(full_model, without_x1)   # nested: evidence for x1
bayes_factor_flat_r(full_model, other_set)    # non-nested comparison

## Works on plain lm() fits too:
bayes_factor_flat_r(lm(y ~ x1 + x2 + x3, data = my_data), lm(y ~ 1, data = my_data))

## ---- 5. What a null Bayes factor rests on ---------------------------
## Which effect sizes can the data not tell apart from zero? (x2 here)
support_x2 <- support_per_partial_correlation(lm(y ~ x1 + x2 + x3, data = my_data), "x2",
                                              partial_correlations = seq(-0.9, 0.9, by = 0.01))
support_regions(support_x2, band = 3)   # better / can't tell / worse than zero

ggplot(support_x2, aes(partial_correlation, support_vs_zero)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1 / 3, ymax = 3, alpha = 0.15) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  labs(x = "true partial correlation of x2", y = "support compared with zero (log scale)",
       title = "x2: which effect sizes can't be told apart from zero?",
       subtitle = "grey band: within a factor 3 of zero")

## ---- 6. Comparison with the default (JZS / JASP) --------------------
## Optional; needs the BayesFactor package
if (requireNamespace("BayesFactor", quietly = TRUE)) {
  default_bf_x1 <- BayesFactor::lmBF(y ~ x1 + x2 + x3, data = my_data) /
                   BayesFactor::lmBF(y ~ x2 + x3, data = my_data)
  print(default_bf_x1)
}
