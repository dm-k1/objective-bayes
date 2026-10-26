## =====================================================================
##  Vignette: the flat-R prior on one dataset, start to finish
##  Estimation (Stan), then Bayes factors for any model comparison.
## =====================================================================
library(flatr)   # install once, from the project folder: devtools::install()
library(ggplot2)

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

## ---- 2. Estimation ---------------------------------------------------
full_model <- fit_flat_r_model(y ~ x1 + x2 + x3, data = my_data, chains = 4,
                               iter_warmup = 1000, iter_sampling = 1000, seed = 1)
full_model                       # estimates, intervals, Rhat, ESS, divergences
loo(full_model)                  # leave-one-out cross-validation
pp_check(full_model)             # posterior predictive check
bayesplot::mcmc_trace(full_model$draws)   # trace plots (any bayesplot plot works on $draws)

## Next to the ordinary estimates: at n = 60 with three predictors the prior
## pulls the estimates toward zero (docs/1_from_behaviour.md, step 4)
cbind(flat_r   = colMeans(as.data.frame(full_model$draws)[, c("x1", "x2", "x3")]),
      ordinary = coef(lm(y ~ x1 + x2 + x3, data = my_data))[c("x1", "x2", "x3")])

## ---- 3. Bayes factor per predictor ----------------------------------
## No Stan needed. Each predictor is tested on its own partial correlation
## (the conditional test); test = "joint" compares the models with and without it
bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, my_data)
bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, my_data, test = "joint")
## The other prior (docs/0_overview.md): every correlation matrix equally
## likely, so each coefficient's partial correlation is flat, whatever k is
bayes_factor_flat_r_per_term(y ~ x1 + x2 + x3, my_data, prior = "uniform_correlations")

## ---- 4. Comparing any two models ------------------------------------
## Takes fit_flat_r_model() fits or plain lm() fits: only R-squared, n and k matter
bayes_factor_flat_r(full_model, lm(y ~ x2 + x3, data = my_data))   # nested: evidence for x1
bayes_factor_flat_r(full_model, lm(y ~ x4, data = my_data))        # non-nested comparison
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
