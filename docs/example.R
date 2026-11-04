## =====================================================================
##  Example: flatcor on one dataset, start to finish
##  Estimation (Stan), then Bayes factors for any model comparison, with the
##  three priors side by side: flat_partial_r, flat_r and flat_r_squared.
##  There is no default; every function asks for prior = (docs/0_overview.md).
## =====================================================================
library(flatcor)   # install once, from the project folder: devtools::install()
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
full_model <- fit_regression(y ~ x1 + x2 + x3, data = my_data, prior = "flat_partial_r",
                             chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1)
full_model                       # estimates, intervals, Rhat, ESS, divergences
loo(full_model)                  # leave-one-out cross-validation
pp_check(full_model)             # posterior predictive check
bayesplot::mcmc_trace(full_model$draws)   # trace plots (any bayesplot plot works on $draws)

## The three priors next to the ordinary estimates: at n = 60 with three
## predictors each prior pulls the estimates toward zero, by a different amount
## (docs/5_validation.md, study 07)
priors <- c("flat_partial_r", "flat_r", "flat_r_squared")
other_fits <- lapply(c(flat_r = "flat_r", flat_r_squared = "flat_r_squared"), function(prior)
  fit_regression(y ~ x1 + x2 + x3, data = my_data, prior = prior, chains = 4,
                 iter_warmup = 1000, iter_sampling = 1000, seed = 1))
posterior_means <- function(fit) colMeans(as.data.frame(fit$draws)[, c("x1", "x2", "x3")])
cbind(flat_partial_r = posterior_means(full_model),
      flat_r         = posterior_means(other_fits$flat_r),
      flat_r_squared = posterior_means(other_fits$flat_r_squared),
      ordinary       = coef(lm(y ~ x1 + x2 + x3, data = my_data))[c("x1", "x2", "x3")])

## ---- 3. Bayes factor per predictor ----------------------------------
## No Stan needed. Each predictor is tested on its own partial correlation
## (the conditional test); test = "joint" compares the models with and without it
bayes_factor_per_term(y ~ x1 + x2 + x3, my_data, prior = "flat_partial_r")
bayes_factor_per_term(y ~ x1 + x2 + x3, my_data, prior = "flat_partial_r", test = "joint")
## The three priors side by side. flat_partial_r keeps each coefficient's
## partial correlation flat for any k; flat_r and flat_r_squared keep the total
## flat, so each coefficient's prior narrows as predictors are added
per_term_bf10 <- sapply(priors, function(prior)
  bayes_factor_per_term(y ~ x1 + x2 + x3, my_data, prior = prior)$BF10)
rownames(per_term_bf10) <- c("x1", "x2", "x3")
per_term_bf10

## ---- 4. Comparing any two models ------------------------------------
## Takes fit_regression() fits or plain lm() fits: only R-squared, n and k matter
bayes_factor_models(full_model, lm(y ~ x2 + x3, data = my_data),            # nested: evidence for x1
                    prior = "flat_partial_r")
bayes_factor_models(full_model, lm(y ~ x4, data = my_data), prior = "flat_partial_r")   # non-nested
sapply(priors, function(prior)                                              # every prior
  bayes_factor_models(lm(y ~ x1 + x2 + x3, data = my_data), lm(y ~ 1, data = my_data), prior = prior))

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

## A term with several columns (a factor with more than two levels): one curve
## for the whole term, over its partial R-squared
my_data$group <- factor(rep(c("a", "b", "c"), length.out = nrow(my_data)))
support_group <- support_per_partial_r_squared(lm(y ~ x1 + x2 + x3 + group, data = my_data),
                                               "group")
support_regions(support_group, band = 3)

## ---- 6. Comparison with the default (JZS / JASP) --------------------
## Optional; needs the BayesFactor package
if (requireNamespace("BayesFactor", quietly = TRUE)) {
  default_bf_x1 <- BayesFactor::lmBF(y ~ x1 + x2 + x3, data = my_data) /
                   BayesFactor::lmBF(y ~ x2 + x3, data = my_data)
  print(default_bf_x1)
}
