## =================================================================
##  "Define what's possible, not what you expect"
##  A step-by-step walkthrough of the possible-range prior
##  for regression Bayes factors
## =================================================================
library(dplyr)
library(ggplot2)
library(BayesFactor)

## The method itself lives in one place: R/possible_range.R. The tutorial
## uses those functions, so it always shows the current version. (Run with
## the project folder as working directory.) The file also provides
## log_sum_exp(): adds up very large or very small numbers safely via their
## logarithms, since probabilities here can be like 1e-300.
source("R/possible_range.R")

## -----------------------------------------------------------------
## STEP 0. The problem
## -----------------------------------------------------------------
## A Bayes factor compares "no effect" with "an effect". But "an
## effect" has to be spelled out as a range of effect sizes with
## weights: the prior. Default tests (JZS / JASP) use a bell curve with
## a fixed width, which is an EXPECTATION about how big effects usually are.
##
## The idea here: don't state what you expect, only what is POSSIBLE,
## and treat every possible value as equally plausible.

## Example data: 3 correlated predictors; x1 and x2 matter, x3 doesn't
set.seed(42)
sample_size            <- 60
predictor_correlations <- matrix(0.4, 3, 3); diag(predictor_correlations) <- 1
predictors <- matrix(rnorm(sample_size * 3), sample_size) %*% chol(predictor_correlations)
colnames(predictors) <- c("x1", "x2", "x3")
outcome    <- drop(predictors %*% c(0.3, 0.2, 0)) + rnorm(sample_size)
regression_fit <- lm(outcome ~ predictors)
summary(regression_fit)

## -----------------------------------------------------------------
## STEP 1. One predictor: what is "possible"?
## -----------------------------------------------------------------
## With one predictor, the effect is fully described by the correlation
## between x and y, and a correlation can only lie between -1 and 1.
## "Every possible value equally plausible" = a flat prior on the
## correlation. This is Jeffreys' classic default for correlations.
##
## Translated to the effect-size scale (signal relative to noise),
## that flat prior becomes a heavy-tailed bell curve. It looks like the
## default Cauchy, but nobody chose its width: it follows from the fact
## that a correlation can't go beyond -1 or 1.

effect_size_values <- seq(-4, 4, length.out = 801)
prior_curves <- bind_rows(
  tibble(effect_size = effect_size_values,
         density     = 0.5 * (1 + effect_size_values^2)^(-1.5),
         prior       = "flat on the correlation (possible range)"),
  tibble(effect_size = effect_size_values,
         density     = dcauchy(effect_size_values, 0, sqrt(2) / 4),
         prior       = "JZS default (Cauchy, width 0.354)"))

ggplot(prior_curves, aes(effect_size, density, colour = prior)) +
  geom_line(linewidth = 1) +
  labs(x = "effect size (signal relative to noise)", y = "prior weight",
       title = "Step 1: a flat correlation prior, seen on the effect-size scale",
       subtitle = "the default puts more weight on small effects: that is its built-in expectation")

## -----------------------------------------------------------------
## STEP 2. Several predictors: the possible range is a shared limit
## -----------------------------------------------------------------
## With several predictors, no single coefficient has a fixed limit,
## but together they can't explain more than all of y's variance:
## the explained share (R-squared) can't exceed 1.
##
## ------------------------------------------------------------------
## Why the explained share (R-squared) gets its own prior
## ------------------------------------------------------------------
## Giving every coefficient its own independent prior sounds neutral,
## but each predictor then adds its own expected share of explained
## variance. The more predictors, the closer the implied R-squared
## gets to 1 before any data is seen. With 10 predictors the
## "every combination within the limit equally plausible" version
## expects R-squared = 0.83 (and gives only 3% to R-squared < 0.5).
## That's a hidden strong claim, and it grows with every predictor added.
##
## Prior frameworks WITH this drift:
##   - independent normal / Student-t priors per coefficient
##       (the usual brms/Stan setup, e.g. normal(0, 1) on standardized
##        betas; rstanarm's default autoscaled normal(0, 2.5))
##   - ridge regression with a fixed penalty (= independent normals)
##   - a uniform box per coefficient, e.g. uniform(-1, 1) on each
##       standardized beta (worse than the ellipsoid: mass sits in corners)
##   - uniform over all combinations within R-squared <= 1 (the ellipsoid)
##   - Zellner's g-prior with fixed g, e.g. g = n (unit information / BIC)
##   - JZS / Zellner-Siow, the multi-predictor default in BayesFactor
##       and JASP (milder, because of the heavy-tailed mixture over g,
##        but the implied R-squared still rises with the number of
##        predictors)
##   - spike-and-slab with a fixed slab width (drifts with the number
##       of included predictors)
##   - independent flat partial correlations per predictor, if combined
##       into one joint model (implied R-squared about 0.98 at 10
##       predictors); fine when each test is a separate one-coefficient
##       comparison
##
## Prior frameworks built to AVOID it (prior on R-squared or on total
## shrinkage, then split over predictors):
##   - R2-D2 (Zhang et al., 2022) and R2D2M2 (Aguilar & Buerkner, 2023)
##   - rstanarm::stan_lm with its R-squared prior
##   - regularized horseshoe (Piironen & Vehtari, 2017), with the global
##       scale tied to the expected number of relevant predictors
##   - hyper-g (Liang et al., 2008) and Maruyama & George (2011),
##       partially: prior on the shrinkage factor g/(1+g)
##   - the "possible range" prior used here: fixed Beta(1/2, 1) on
##       R-squared for any number of predictors, every direction
##       equally plausible
## ------------------------------------------------------------------
##
## The fix: give the explained share ITSELF the same prior, whatever the
## number of predictors, and let every WAY of dividing it over the
## predictors be equally plausible. The particular prior chosen for the
## explained share makes the one-predictor case exactly the flat
## correlation prior of Step 1.

explained_share_draws <- bind_rows(lapply(c(1, 4, 10), function(number_of_predictors) bind_rows(
  tibble(explained_share      = rbeta(20000, number_of_predictors / 2, 1),
         number_of_predictors = number_of_predictors,
         prior = "every combination equally plausible"),
  tibble(explained_share      = rbeta(20000, 0.5, 1),
         number_of_predictors = number_of_predictors,
         prior = "fixed prior on explained share + every direction equally plausible"))))

explained_share_draws |>
  group_by(prior, number_of_predictors) |>
  summarise(average_explained_share = mean(explained_share), .groups = "drop") |>
  print()

ggplot(explained_share_draws, aes(explained_share, colour = prior)) +
  geom_density(linewidth = 1, bounds = c(0, 1)) +
  facet_wrap(~ number_of_predictors, labeller = label_both) +
  labs(x = "explained share of variance (R-squared)",
       title = "Step 2: why the explained share gets its own prior",
       subtitle = "the 'every combination' version drifts toward 1 as predictors are added")

## ---- Why Beta(1/2, 1)? It is "flat on the correlation" on another ruler --
## With one predictor, R-squared is just the correlation squared. So a
## prior on the correlation IS a prior on R-squared: one statement on two
## scales. "Every correlation equally plausible" (flat from -1 to 1)
## becomes Beta(1/2, 1) on the R-squared scale.
##
## The curve looks tall near zero, but nothing is favoured there: squaring
## squeezes the small correlations together. Half of all correlations
## (|r| < 0.5) land in R-squared < 0.25, a quarter of the axis; the other
## half spread over the remaining three quarters. Like a square's side and
## its area: sides picked evenly from 0 to 1 give areas piled up near 0.
##
## The same event gets the same probability on both rulers:
## |r| < 0.5  <=>  R-squared < 0.25  (both 50%).
## This is also why "every value equally plausible" needs a chosen ruler:
## flat on r and flat on R-squared are different priors. The correlation
## is the natural "possible range" ruler: bounded, signed, and the same in
## every study.

correlation_draws <- runif(100000, -1, 1)   # every correlation equally plausible
r_squared_draws   <- correlation_draws^2    # the same draws, on the R-squared ruler

par(mfrow = c(1, 2))
hist(correlation_draws, breaks = 40, freq = FALSE, col = "grey85", border = "white",
     main = "Flat on the correlation", xlab = "correlation (one predictor)")
abline(v = c(-0.5, 0.5), lty = 2)
mtext(sprintf("share with |r| < 0.5: %.2f", mean(abs(correlation_draws) < 0.5)), cex = 0.8)

hist(r_squared_draws, breaks = 40, freq = FALSE, col = "grey85", border = "white",
     main = "...is Beta(1/2, 1) on R-squared", xlab = "R-squared", ylim = c(0, 5))
curve(dbeta(x, 0.5, 1), from = 0.002, to = 1, add = TRUE, lwd = 2)   # the formula
abline(v = 0.25, lty = 2)
mtext(sprintf("share with R-squared < 0.25: %.2f", mean(r_squared_draws < 0.25)), cex = 0.8)
par(mfrow = c(1, 1))

## -----------------------------------------------------------------
## STEP 3. Building the prior, and seeing what it implies
## -----------------------------------------------------------------
## To draw one set of coefficients (in units of the noise):
##   1. draw how much variance is explained (the explained share)
##   2. turn that into an overall effect strength
##   3. draw a random direction: how the effect is divided over predictors
##   4. reshape that direction to the predictors' correlations, so that
##      correlated predictors are handled automatically
## draw_coefficients_from_possible_range_prior() in R/possible_range.R
## (Part 4) does exactly these four steps; open the file to read them.

## The partial correlation a coefficient implies: its unique link with
## the outcome after removing what the other predictors explain
implied_partial_correlation <- function(effects, predictor_covariance, predictor_index) {
  spread_left_after_others <- sqrt(1 / solve(predictor_covariance)[predictor_index, predictor_index])
  scaled_effect <- effects[, predictor_index] * spread_left_after_others
  scaled_effect / sqrt(1 + scaled_effect^2)
}

## Check: with one predictor, the implied correlation is flat from -1 to 1
prior_draws_one_predictor <- draw_coefficients_from_possible_range_prior(50000, matrix(1))
hist(implied_partial_correlation(prior_draws_one_predictor, matrix(1), 1),
     breaks = 40, main = "Step 3a: one predictor gives a flat correlation prior",
     xlab = "correlation")

## With 4 correlated predictors, each one's share is narrower: if
## together they explain at most everything, one alone rarely explains
## a lot. Nobody assumed small effects; it follows from the shared
## limit. The price: each single test now behaves more like the JZS
## default (slower to find "no effect", less often fooled by small
## real effects) than the one-predictor flat prior does.
correlations_four_predictors <- matrix(0.4, 4, 4); diag(correlations_four_predictors) <- 1
prior_draws_four_predictors  <- draw_coefficients_from_possible_range_prior(
  50000, correlations_four_predictors)
partial_correlations_x1 <- implied_partial_correlation(
  prior_draws_four_predictors, correlations_four_predictors, 1)
hist(partial_correlations_x1, breaks = 60,
     main = "Step 3b: one predictor's share when there are 4", xlab = "partial correlation")
quantile(abs(partial_correlations_x1), c(0.25, 0.5, 0.75))

## ---- Step 3c: from flat to peaked as predictors are added ---------------
## The total (R-squared) keeps the same Beta(1/2, 1) prior for any number
## of predictors. What changes is how it is shared: a random direction
## splits the total over the predictors, so each one gets a random slice,
## on average 1/k of it. Like a cake whose size never changes: one person
## gets all of it; with four, each slice is usually smaller (and sometimes
## one person gets nearly everything).
##
## So one predictor's partial correlation goes from flat (k = 1) to more
## and more peaked at zero. Nobody stated "I expect small effects": it
## follows from "the total can't exceed 100%" plus "I don't know how it
## is split". (Uncorrelated predictors here, to show the pure effect of k.)

par(mfrow = c(1, 4))
for (number_of_predictors in c(1, 2, 4, 10)) {
  uncorrelated    <- diag(number_of_predictors)
  draws           <- draw_coefficients_from_possible_range_prior(50000, uncorrelated)
  first_predictor <- implied_partial_correlation(draws, uncorrelated, 1)
  hist(first_predictor, breaks = 50, xlim = c(-1, 1), col = "grey85", border = "white",
       main = paste(number_of_predictors, "predictor(s)"), xlab = "partial correlation of x1")
  mtext(sprintf("median |r| = %.2f", median(abs(first_predictor))), cex = 0.8)
}
par(mfrow = c(1, 1))

## -----------------------------------------------------------------
## STEP 4. The Bayes factor as a formula
## -----------------------------------------------------------------
## With the usual neutral priors on the intercept and the noise level
## (the same ones the JZS default uses), the Bayes factor of a model
## with predictors against the intercept-only model turns out to
## depend on just three things: the sample size, the number of
## predictors, and the observed explained share (R-squared).
##
## The formula averages over two things:
##  - an overall scale (like the JZS family does)
##  - every direction the effect could take, compared with the
##    direction the data actually point in: this is the
##    "every direction equally plausible" part

## The function that computes this is
## log_bayes_factor_possible_range_vs_intercept_only() in R/possible_range.R.
## Open that file to read the code: each line is commented, and
## docs/derivation.md (section 4) explains every step.
## Its inputs are exactly the three numbers above:
log_bayes_factor_possible_range_vs_intercept_only(explained_share = 0.19, sample_size = 60,
                                                  number_of_predictors = 3)

## -----------------------------------------------------------------
## STEP 5. Check the formula against brute force
## -----------------------------------------------------------------
## A Bayes factor means: how well did the prior predict the data,
## compared with "no effect"? Brute force does exactly that, with no
## formula: draw many coefficient sets from the prior (Step 3), see how
## well each predicts the outcome, and average. If this matches Step 4,
## the formula is right.

## log_bayes_factor_possible_range_by_sampling() in R/possible_range.R
## (Part 4) does this: the noise level is unknown, so it is averaged over
## on a grid, for "no effect" and for each drawn coefficient set.

explained_share_full_model <- summary(regression_fit)$r.squared
c(formula     = log_bayes_factor_possible_range_vs_intercept_only(explained_share_full_model,
                                                                  sample_size, 3),
  brute_force = log_bayes_factor_possible_range_by_sampling(predictors, outcome,
                                                            number_of_draws = 20000))
## These should agree to about 0.05 (on the log scale).

## -----------------------------------------------------------------
## STEP 6. Testing one predictor
## -----------------------------------------------------------------
## Evidence for predictor j = how well the model WITH it predicts,
## divided by how well the model WITHOUT it predicts, both using the
## same kind of prior. Same structure as the default test, so the two
## are directly comparable.

test_results <- lapply(1:3, function(predictor_index) {
  explained_share_without <- summary(lm(outcome ~ predictors[, -predictor_index]))$r.squared
  tibble(
    predictor = colnames(predictors)[predictor_index],
    estimate  = coef(regression_fit)[predictor_index + 1],
    p_value   = summary(regression_fit)$coefficients[predictor_index + 1, 4],
    evidence_for_effect_possible =
      exp(log_bayes_factor_possible_range_vs_intercept_only(explained_share_full_model,
                                                            sample_size, 3) -
          log_bayes_factor_possible_range_vs_intercept_only(explained_share_without,
                                                            sample_size, 2)),
    evidence_for_effect_default =
      linearReg.R2stat(sample_size, 3, explained_share_full_model, simple = TRUE) /
      linearReg.R2stat(sample_size, 2, explained_share_without,    simple = TRUE))
}) |>
  bind_rows() |>
  mutate(evidence_for_no_effect_possible = 1 / evidence_for_effect_possible)
print(test_results)

## -----------------------------------------------------------------
## STEP 7. What the evidence rests on: support for each effect size
## -----------------------------------------------------------------
## A Bayes factor that favours "no effect" mostly says LARGE effects are
## ruled out. To see which effect sizes the data can and can't tell
## apart from zero, compare each possible partial correlation directly
## with zero. Always report this next to a Bayes factor that favours
## "no effect".
##
## Note: the peak of each curve is the most a Bayes factor could ever
## be (you'd get it only by picking the best effect size AFTER seeing
## the data). The Bayes factor itself is this curve averaged over
## everything the prior says is possible, so it is always lower.


## support_per_partial_correlation() is in R/possible_range.R. It takes the
## fitted model and the name of one coefficient (here "predictorsx1", ...).
support_chart_data <- bind_rows(lapply(colnames(predictors), function(predictor_name)
  mutate(support_per_partial_correlation(regression_fit, paste0("predictors", predictor_name)),
         predictor = predictor_name)))

## Best-supported effect size per predictor
support_chart_data |>
  group_by(predictor) |>
  summarise(best_supported_value = partial_correlation[which.max(support_vs_zero)],
            best_support_vs_zero = max(support_vs_zero)) |>
  print()

## Split the axis into regions: where does each effect size explain the
## data better than zero, worse than zero, or about the same (within a
## factor 3)? Read each predictor from left to right.
## support_regions() (also in R/possible_range.R) does this for one predictor.
support_regions_per_predictor <- bind_rows(lapply(
  split(support_chart_data, support_chart_data$predictor),
  function(one_predictor) data.frame(predictor = one_predictor$predictor[1],
                                     support_regions(one_predictor))))
print(support_regions_per_predictor)

ggplot(support_chart_data, aes(partial_correlation, support_vs_zero)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1/3, ymax = 3, alpha = 0.15) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e-4, NA)) +
  facet_wrap(~ predictor) +
  labs(x = "true partial correlation", y = "support compared with zero (log scale)",
       title = "Step 7: which effect sizes do the data support?",
       subtitle = "grey band: data can't tell this value apart from zero (within a factor 3)")

## -----------------------------------------------------------------
## STEP 8. Estimation
## -----------------------------------------------------------------
## For estimates and intervals, use fit_possible_range_model() in R/possible_range.R.
## Note: with several predictors at a small sample, the prior pulls the
## estimates toward zero (in this example about 20% for x1: 0.31 vs 0.39),
## because each predictor's share of the fixed R-squared budget is peaked
## at zero (Step 3b). At larger samples the estimates approach the
## ordinary confidence interval below.
confint(regression_fit)

## -----------------------------------------------------------------
## STEP 9. Checks every good default prior should pass
## -----------------------------------------------------------------
## a) With the smallest possible sample (one more observation than
##    predictors), the data can't tell models apart, so the Bayes
##    factor should be exactly 1.
sapply(1:4, function(number_of_predictors)
  exp(log_bayes_factor_possible_range_vs_intercept_only(1 - 1e-12, number_of_predictors + 1,
                                                        number_of_predictors)))

## b) As the fit becomes perfect, evidence for an effect should grow
##    without limit.
sapply(c(0.9, 0.99, 0.999, 0.9999), function(explained_share)
  log_bayes_factor_possible_range_vs_intercept_only(explained_share, 20, 3))

## c) When a predictor truly does nothing, evidence for "no effect"
##    should keep growing as the sample grows.
sapply(c(50, 200, 1000, 5000), function(sample_size_check) {
  explained_share_without <- 0.2
  explained_share_with    <- explained_share_without + (1 - explained_share_without) / sample_size_check
  -(log_bayes_factor_possible_range_vs_intercept_only(explained_share_with,    sample_size_check, 3) -
    log_bayes_factor_possible_range_vs_intercept_only(explained_share_without, sample_size_check, 2))
})
