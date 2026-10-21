## =====================================================================
##  The posterior, step by step: one predictor, everything on a grid
##
##  With one predictor the possible values can be listed one by one, so
##  every step of the Bayesian calculation can be printed and plotted:
##  prior, likelihood, posterior, and the Bayes factor. MCMC (brms) does the
##  same thing by sampling when there are too many predictors for a grid.
##
##  Run with the project folder as working directory (step 9 uses the
##  method file to compare). Base R only.
## =====================================================================

## ---- Step 0. Data -------------------------------------------------------
set.seed(1)
sample_size <- 30
x <- rnorm(sample_size)
y <- 0.5 * x + rnorm(sample_size)

plot(x, y, pch = 16, main = "Step 0: the data", xlab = "x", ylab = "y")
abline(lm(y ~ x), lwd = 2)

## Center both variables. With a flat prior on the intercept, averaging over
## every possible intercept is the same as working with centered data and
## n - 1 instead of n (docs/derivation.md, section 4).
x_centered <- x - mean(x)
y_centered <- y - mean(y)
x_spread   <- sqrt(mean(x_centered^2))   # SD of x (n in the denominator, as in the prior)

## Three sums are all the data the calculation needs
sum_yy <- sum(y_centered^2)
sum_xy <- sum(x_centered * y_centered)
sum_xx <- sum(x_centered^2)


## ---- Step 1. Every possible correlation, with its prior ------------------
## The possible range is -1 to 1, every value equally plausible: a flat prior
## with density 1/2. We check correlations in steps of 0.01.
correlation_step <- 0.01
correlation_grid <- seq(-0.995, 0.995, by = correlation_step)
prior_weight_correlation <- 0.5 * correlation_step   # density x step = prior probability per point
sum(rep(prior_weight_correlation, length(correlation_grid)))   # adds up to (almost) 1

plot(correlation_grid, rep(0.5, length(correlation_grid)), type = "l", lwd = 2,
     ylim = c(0, 1), xlab = "correlation", ylab = "prior density",
     main = "Step 1: prior on the correlation\nevery possible value equally plausible")


## ---- Step 2. Every possible noise level, with its prior ------------------
## The noise prior is 1/sigma: every ORDER OF MAGNITUDE equally plausible.
## On a grid evenly spaced in log(sigma), every point gets the same weight.
log_noise_grid <- seq(log(sd(y)) - 3, log(sd(y)) + 1, length.out = 300)
noise_grid     <- exp(log_noise_grid)
log_noise_step <- log_noise_grid[2] - log_noise_grid[1]

par(mfrow = c(1, 2))
plot(noise_grid, 1 / noise_grid, type = "l", lwd = 2,
     xlab = "noise SD", ylab = "prior density (unnormalized)",
     main = "Step 2: prior on the noise\n1/sigma on the ordinary scale ...")
plot(log_noise_grid, rep(1, length(log_noise_grid)), type = "l", lwd = 2,
     ylim = c(0, 2), xlab = "log(noise SD)", ylab = "prior density (unnormalized)",
     main = "... is flat on the log scale")
par(mfrow = c(1, 1))


## ---- Step 3. All combinations ------------------------------------------
grid <- expand.grid(correlation = correlation_grid, noise_sd = noise_grid)
nrow(grid)   # 200 correlations x 300 noise levels = 60,000 combinations


## ---- Step 4. Turn each correlation into a slope -------------------------
## From "R-squared can't exceed 1" to a slope, in four small steps:
##  1. The outcome's variance = explained part + noise. The explained share
##     (R-squared) is explained / total, so it lies between 0 and 1.
##  2. With one predictor, R-squared is the correlation squared, so the
##     correlation r lies between -1 and 1. That is the possible range.
##  3. Explained part divided by noise = R2 / (1 - R2) = r^2 / (1 - r^2).
##     Its square root, w = r / sqrt(1 - r^2), is the effect in noise units:
##     how many noise-SDs y moves when x moves one SD. (This is the
##     standardized effect size of Step 4b.)
##  4. The slope is the same thing in the data's own units:
##     slope = w * noise_sd / SD(x).
## Example: r = 0.5 -> w = 0.5 / sqrt(0.75) = 0.58 -> with noise SD 1 and
## SD(x) 2, slope = 0.58 * 1 / 2 = 0.29.
grid$slope <- grid$noise_sd * grid$correlation / sqrt(1 - grid$correlation^2) / x_spread

## The same correlation means a bigger slope when there is more noise
plot(NA, xlim = c(-1, 1), ylim = c(-3, 3), xlab = "correlation", ylab = "slope",
     main = "Step 4: from correlation to slope\n(one curve per noise level)")
for (noise_level in c(0.5, 1, 2)) {
  lines(correlation_grid,
        noise_level * correlation_grid / sqrt(1 - correlation_grid^2) / x_spread, lwd = 2,
        lty = match(noise_level, c(0.5, 1, 2)))
}
abline(h = 0, v = 0, col = "grey")
legend("topleft", legend = paste("noise SD =", c(0.5, 1, 2)), lty = 1:3, lwd = 2, bty = "n")

## The same prior on the standardized effect-size scale, to compare with
## the Cauchy priors of the JZS default.
## The effect in noise units is w = r / sqrt(1 - r^2) (slope x SD(x) / noise SD).
## If r is flat from -1 to 1, w has density 0.5 * (1 + w^2)^(-3/2):
## a Student-t with 2 degrees of freedom whose width nobody chose (it follows
## from the bound of the correlation). The Cauchy is a Student-t with 1
## degree of freedom whose width IS chosen.
##
## Units: JASP's t-test widths are on Cohen's d. For two equal groups the
## effect in noise units is w = d / 2, so a width on d is halved here:
##   JASP t-test "medium"    0.707 on d  ->  0.354 on w  (= BayesFactor regression "medium")
##   JASP t-test "ultrawide" 1.414 on d  ->  0.707 on w
## On d, the flat correlation prior is a Student-t (2 df) with scale 1.41.
## Check by simulation: draw correlations from the flat prior, transform them,
## and compare the histogram with the formula.
correlation_draws <- runif(200000, -1, 1)
effect_size_draws <- correlation_draws / sqrt(1 - correlation_draws^2)

effect_size_values <- seq(-4, 4, length.out = 801)
## Histogram of the draws within -4 to 4, scaled to the full range (about 3%
## of the draws lie beyond |w| = 4, so heights are multiplied by the share kept)
shown <- abs(effect_size_draws) < 4
effect_size_histogram <- hist(effect_size_draws[shown], breaks = 160, plot = FALSE)
effect_size_histogram$density <- effect_size_histogram$density * mean(shown)
plot(effect_size_histogram, freq = FALSE, col = "grey90", border = "white",
     xlim = c(-4, 4), ylim = c(0, 1.35),
     xlab = "standardized effect size (effect in noise units)", ylab = "prior density",
     main = "Step 4b: the flat correlation prior on the effect-size scale\nvs the Cauchy defaults")
lines(effect_size_values, 0.5 * (1 + effect_size_values^2)^(-1.5), lwd = 3)
lines(effect_size_values, dcauchy(effect_size_values, 0, sqrt(2) / 4), lwd = 2, lty = 2)
lines(effect_size_values, dcauchy(effect_size_values, 0, sqrt(2) / 2), lwd = 2, lty = 3)
legend("top", bty = "n", lwd = c(3, 2, 2), lty = 1:3, cex = 0.75, seg.len = 3,
       legend = c("flat on r (formula; histogram = simulated)",
                  "Cauchy 0.354: JZS default (regression 'medium' = t-test 0.707 on d)",
                  "Cauchy 0.707: JASP t-test 'ultrawide' (1.41 on d)"))

## How much prior weight on small, medium and large effects?
## (Cohen's rough labels for the effect in noise units: 0.2, 0.5, 0.8)
prior_share_below <- function(cutoff) c(
  flat_r           = mean(abs(effect_size_draws) < cutoff),
  jzs_default      = 2 * pcauchy(cutoff, 0, sqrt(2) / 4) - 1,
  jasp_ultrawide   = 2 * pcauchy(cutoff, 0, sqrt(2) / 2) - 1)
round(rbind(`|effect| < 0.2` = prior_share_below(0.2),
            `|effect| < 0.5` = prior_share_below(0.5),
            `|effect| < 0.8` = prior_share_below(0.8),
            `|effect| < 2`   = prior_share_below(2)), 2)


## ---- Step 5. How well does each combination fit the data? ----------------
## Residual sum of squares for each slope, written with the three sums:
## sum((y - slope * x)^2) = sum_yy - 2 * slope * sum_xy + slope^2 * sum_xx
grid$residual_ss <- sum_yy - 2 * grid$slope * sum_xy + grid$slope^2 * sum_xx

## Log likelihood (normal errors, intercept averaged out). Constants that are
## the same for every combination are left out; they cancel later.
grid$log_likelihood <- -(sample_size - 1) * log(grid$noise_sd) -
                       grid$residual_ss / (2 * grid$noise_sd^2)

## One combination in detail: r = 0.40, noise SD closest to 1
example_row <- which.min(abs(grid$correlation - 0.40) + abs(grid$noise_sd - 1))
grid[example_row, ]
## the same residual sum of squares, computed the long way:
sum((y_centered - grid$slope[example_row] * x_centered)^2)

## Picture: the likelihood over all combinations (darker = fits better)
as_grid_matrix <- function(values) matrix(values, nrow = length(correlation_grid))
likelihood_surface <- as_grid_matrix(exp(grid$log_likelihood - max(grid$log_likelihood)))
image(correlation_grid, log_noise_grid, likelihood_surface,
      col = gray.colors(50, start = 1, end = 0),
      xlab = "correlation", ylab = "log(noise SD)",
      main = "Step 5: likelihood of every combination\n(darker = fits the data better)")
contour(correlation_grid, log_noise_grid, likelihood_surface, add = TRUE,
        levels = c(0.05, 0.25, 0.5, 0.75), drawlabels = FALSE)
points(grid$correlation[example_row], log(grid$noise_sd[example_row]), pch = 4, cex = 2, lwd = 2)


## ---- Step 6. Prior x likelihood -----------------------------------------
## Every grid point has the same prior weight (steps 1 and 2), so on the log
## scale the prior only adds a constant: the posterior has exactly the shape
## of the likelihood. (With a non-flat prior this is where the shape changes.)
grid$log_prior     <- log(prior_weight_correlation) + log(log_noise_step)
grid$log_posterior <- grid$log_likelihood + grid$log_prior   # not yet normalized


## ---- Step 7. Normalize: the posterior adds up to 1 ----------------------
## Subtract the largest value first so exp() doesn't overflow.
grid$posterior <- exp(grid$log_posterior - max(grid$log_posterior))
grid$posterior <- grid$posterior / sum(grid$posterior)
sum(grid$posterior)   # 1

image(correlation_grid, log_noise_grid, as_grid_matrix(grid$posterior),
      col = gray.colors(50, start = 1, end = 0),
      xlab = "correlation", ylab = "log(noise SD)",
      main = "Steps 6-7: posterior over all combinations\n(same shape as the likelihood: the prior is flat)")


## ---- Step 8. Read off what you want -------------------------------------
## Posterior for the correlation: add up over all noise levels
posterior_correlation <- tapply(grid$posterior, grid$correlation, sum)

plot(correlation_grid, posterior_correlation / correlation_step, type = "l", lwd = 2,
     xlab = "correlation", ylab = "density",
     main = "Step 8: posterior of the correlation\n(dashed: the flat prior it started from)")
lines(correlation_grid, rep(0.5, length(correlation_grid)), lty = 2)
abline(v = cor(x, y), col = "grey")   # the sample correlation

## Summaries for any quantity: weighted by the posterior
weighted_quantile <- function(values, weights, probabilities) {
  order_of_values <- order(values)
  cumulative <- cumsum(weights[order_of_values]) / sum(weights)
  sapply(probabilities, function(p) values[order_of_values][which(cumulative >= p)[1]])
}
c(posterior_mean_correlation = sum(grid$correlation * grid$posterior),
  lower_95 = weighted_quantile(grid$correlation, grid$posterior, 0.025),
  upper_95 = weighted_quantile(grid$correlation, grid$posterior, 0.975))

slope_summary <- c(posterior_median_slope = weighted_quantile(grid$slope, grid$posterior, 0.5),
                   lower_95 = weighted_quantile(grid$slope, grid$posterior, 0.025),
                   upper_95 = weighted_quantile(grid$slope, grid$posterior, 0.975))
slope_summary

## For comparison: the ordinary estimate and confidence interval
ordinary_fit <- lm(y ~ x)
coef(ordinary_fit)["x"]; confint(ordinary_fit)["x", ]

## Picture: posterior of the slope, with the ordinary interval
slope_breaks <- seq(-1, 2, by = 0.025)
inside <- grid$slope > min(slope_breaks) & grid$slope < max(slope_breaks)
slope_bins <- cut(grid$slope[inside], slope_breaks)
posterior_slope <- tapply(grid$posterior[inside], slope_bins, sum)
posterior_slope[is.na(posterior_slope)] <- 0
plot(head(slope_breaks, -1) + 0.0125, posterior_slope / 0.025, type = "l", lwd = 2,
     xlab = "slope", ylab = "density",
     main = "Step 8: posterior of the slope\n(grey: ordinary 95% confidence interval)")
abline(v = confint(ordinary_fit)["x", ], col = "grey", lty = 2)
abline(v = slope_summary[c("lower_95", "upper_95")], lty = 3)


## ---- Step 9. The Bayes factor from the same grid -------------------------
## Evidence = prior-weighted AVERAGE of the likelihood (not normalized).
## Compared with the model without an effect (slope = 0).
log_sum_exp <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

log_likelihood_null <- -(sample_size - 1) * log(noise_grid) - sum_yy / (2 * noise_grid^2)
log_evidence_null   <- log_sum_exp(log_likelihood_null + log(log_noise_step))

## For each correlation: how much better than "no effect" does it explain the
## data (averaged over the noise level)? This is the support curve.
log_support_per_correlation <- tapply(grid$log_likelihood + log(log_noise_step),
                                      grid$correlation, log_sum_exp) - log_evidence_null
support_per_correlation <- exp(log_support_per_correlation)

## The Bayes factor is this curve averaged over the flat prior
grid_bayes_factor <- sum(support_per_correlation * prior_weight_correlation)

plot(correlation_grid, support_per_correlation, type = "l", lwd = 2, log = "y",
     ylim = c(1e-4, max(support_per_correlation) * 2),
     xlab = "correlation", ylab = "support compared with no effect (log scale)",
     main = "Step 9: the Bayes factor is the support curve averaged over the prior")
abline(h = 1, lty = 3)
abline(h = grid_bayes_factor, lty = 2, lwd = 2)
legend("bottomright", bty = "n", lty = c(1, 2, 3), lwd = c(2, 2, 1),
       legend = c("support per correlation", sprintf("average = Bayes factor = %.2f", grid_bayes_factor),
                  "no difference from 'no effect'"))

## Same number from the closed form in the method file
source("R/possible_range.R")
closed_form_bayes_factor <- exp(log_bayes_factor_possible_range_vs_intercept_only(cor(x, y)^2,
                                                                                  sample_size, 1))
c(grid = grid_bayes_factor, closed_form = closed_form_bayes_factor)
