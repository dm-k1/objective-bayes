## =====================================================================
##  Figures for docs/1_from_behaviour.md and docs/2_how_it_calculates.md
##  Every figure is a simulation compared with a formula, so each claim in
##  the docs can be checked by running this file. The three priors keep the
##  colours of validation/00: flat partial r blue, flat R orange, flat
##  R-squared green.
##  Run from the project folder: Rscript docs/figures/make_figures.R
## =====================================================================
library(ggplot2)
pkgload::load_all(".", quiet = TRUE)   # the flatcor package, from this folder

figure_folder <- file.path("docs", "figures")
save_doc_figure <- function(figure, name, width = 8, height = 4.5) {
  ggsave(file.path(figure_folder, paste0(name, ".png")), figure,
         width = width, height = height, dpi = 150)
  invisible(figure)
}
theme_set(theme_minimal(base_size = 12))
prior_colours <- c(`flat partial r` = "#0072B2", `flat R` = "#D55E00", `flat R-squared` = "#009E73")
prior_names   <- c(`flat partial r` = "flat_partial_r", `flat R` = "flat_r", `flat R-squared` = "flat_r_squared")
set.seed(2026)
number_of_draws <- 200000


## ---- 1. One predictor: flat on r, seen on the effect-size scale ---------
## r flat from -1 to 1. The effect in noise units is f = r / sqrt(1 - r^2)
## (Cohen's f with a sign). Formula: 0.5 * (1 + f^2)^(-3/2). With one
## predictor flat partial r and flat R are both this prior; flat R-squared
## (r^2 flat) is the curve |f| / (1 + f^2)^2, which is zero at f = 0.
correlation_draws <- runif(number_of_draws, min = -1, max = 1)
signed_f_draws    <- correlation_draws / sqrt(1 - correlation_draws^2)
f_values          <- seq(-4, 4, by = 0.01)
shown_draws       <- signed_f_draws[abs(signed_f_draws) < 4]

curves_one_predictor <- rbind(
  data.frame(f = f_values, density = 0.5 * (1 + f_values^2)^(-3 / 2),
             prior = "flat on r (formula)"),
  data.frame(f = f_values, density = abs(f_values) / (1 + f_values^2)^2,
             prior = "flat on r^2 (formula)"),
  data.frame(f = f_values, density = dcauchy(f_values, location = 0, scale = sqrt(2) / 4),
             prior = "JZS default: Cauchy, scale 0.354"))

save_doc_figure(
  ggplot() +
    geom_histogram(aes(x = shown_draws, y = after_stat(density) * mean(abs(signed_f_draws) < 4)),
                   breaks = seq(-4, 4, by = 0.1), fill = "grey85", colour = "white") +
    geom_line(data = curves_one_predictor, aes(f, density, colour = prior), linewidth = 1) +
    scale_colour_manual(values = c("flat on r (formula)" = "black",
                                   "flat on r^2 (formula)" = "#009E73",
                                   "JZS default: Cauchy, scale 0.354" = "grey50")) +
    labs(x = "effect in noise units (f, signed)", y = "prior density", colour = NULL,
         title = "One predictor: a flat correlation is a heavy-tailed curve on the effect scale",
         subtitle = "grey bars: 200,000 simulated flat correlations, transformed; nobody chose its width") +
    theme(legend.position = "bottom"),
  "1_one_predictor")


## ---- 2. What each prior claims about the total as k grows ----------------
## Implied R-squared before any data, as predictors are added (uncorrelated
## predictors). Points: the median of draws from each prior's recipe
## (draw_coefficients_from_prior()); lines: the formula, the median of
## Beta(a, 1), 0.5^(1/a).
##  - flat partial r: a = k/2, the claim grows with k
##  - flat R: a = 1/2, and flat R-squared: a = 1, the same claim for every k
##  - JZS default (reference): Zellner-Siow, f^2 = scale^2 * chi2_k / chi2_1
numbers_of_predictors <- c(1, 2, 3, 5, 10, 20, 30)
drift_draws <- 50000
implied_r_squared <- do.call(rbind, lapply(numbers_of_predictors, function(k) {
  from_recipe <- sapply(prior_names, function(prior) {
    draws <- draw_coefficients_from_prior(drift_draws, diag(k), prior = prior)
    f_squared <- rowSums(draws^2)
    median(f_squared / (1 + f_squared))
  })
  jzs_f_squared <- (sqrt(2) / 4)^2 * rchisq(drift_draws, df = k) / rchisq(drift_draws, df = 1)
  data.frame(number_of_predictors = k,
             prior = c(names(prior_names), "JZS default (scale 0.354)"),
             simulated = c(from_recipe, median(jzs_f_squared / (1 + jzs_f_squared))),
             formula = c(0.5^(1 / sapply(prior_names, prior_shape, number_of_predictors = k)), NA))
}))
print(implied_r_squared)

save_doc_figure(
  ggplot(implied_r_squared, aes(number_of_predictors, colour = prior)) +
    geom_line(aes(y = formula), linewidth = 1, na.rm = TRUE) +
    geom_line(data = subset(implied_r_squared, is.na(formula)), aes(y = simulated),
              linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = simulated), size = 2) +
    scale_x_continuous(breaks = numbers_of_predictors) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = c(prior_colours, "JZS default (scale 0.354)" = "grey50")) +
    labs(x = "number of predictors (k)", y = "median implied R-squared, before any data",
         colour = NULL,
         title = "Flat partial r claims more explained variance with every predictor",
         subtitle = "flat R and flat R-squared keep the same claim for every k; points: draws, lines: Beta(a, 1)") +
    theme(legend.position = "bottom"),
  "2_drift")


## ---- 3. What follows: one predictor's partial correlation as k grows -----
## Nobody assumed small effects; under flat R and flat R-squared it follows
## from a fixed total, split evenly. Under flat partial r it stays flat.
partial_correlation_of_first <- function(coefficients) {
  ## uncorrelated predictors, noise 1: the partial correlation of x1 is
  ## b1 / sqrt(b1^2 + 1)
  coefficients[, 1] / sqrt(coefficients[, 1]^2 + 1)
}
share_draws <- do.call(rbind, lapply(c(1, 2, 4, 10), function(k)
  do.call(rbind, lapply(names(prior_names), function(prior)
    data.frame(number_of_predictors = paste("k =", k), prior = prior,
               partial_correlation = partial_correlation_of_first(
                 draw_coefficients_from_prior(50000, diag(k), prior = prior_names[[prior]])))))))
share_draws$number_of_predictors <- factor(share_draws$number_of_predictors,
                                           levels = paste("k =", c(1, 2, 4, 10)))
share_draws$prior <- factor(share_draws$prior, levels = names(prior_names))

save_doc_figure(
  ggplot(share_draws, aes(partial_correlation, fill = prior)) +
    geom_histogram(aes(y = after_stat(density)), breaks = seq(-1, 1, by = 0.05),
                   colour = "white") +
    facet_grid(prior ~ number_of_predictors) +
    scale_fill_manual(values = prior_colours, guide = "none") +
    labs(x = "partial correlation of x1", y = "prior density",
         title = "One predictor's partial correlation: flat for every k only under flat partial r",
         subtitle = "flat R and flat R-squared keep the total fixed, so a predictor that shares gets a peak at zero"),
  "3_share_per_predictor", width = 10, height = 6.5)


## ---- 4. The three factors along one line (flat R) -------------------------
## Factor 1: P(total) per unit of f, here for flat R (for R-squared ~ Beta(a, 1)
## it is 2a R^(2a - 1) (1 + f^2)^(-3/2)); factor 2: 1 / size of all combinations
## (2, 2 pi f, 4 pi f^2 for k = 1, 2, 3); together the density per point.
f_grid <- seq(0.02, 3, by = 0.01)
sphere_size <- function(f, k) 2 * pi^(k / 2) / gamma(k / 2) * f^(k - 1)
factor_curves <- do.call(rbind, lapply(1:3, function(k) rbind(
  data.frame(f = f_grid, k = paste("k =", k), value = (1 + f_grid^2)^(-3 / 2),
             factor = "1. P(total), per unit of f"),
  data.frame(f = f_grid, k = paste("k =", k), value = 1 / sphere_size(f_grid, k),
             factor = "2. P(combination | total) = 1 / size"),
  data.frame(f = f_grid, k = paste("k =", k),
             value = (1 + f_grid^2)^(-3 / 2) / sphere_size(f_grid, k),
             factor = "1 x 2: density per point"))))
factor_curves$factor <- factor(factor_curves$factor, levels = unique(factor_curves$factor))

save_doc_figure(
  ggplot(factor_curves, aes(f, value, colour = k, linetype = k)) +
    geom_line(linewidth = 1) +
    ## k in greys, so the colours stay reserved for the three priors
    scale_colour_manual(values = c(`k = 1` = "grey65", `k = 2` = "grey35", `k = 3` = "black")) +
    scale_linetype_manual(values = c(`k = 1` = "solid", `k = 2` = "dashed", `k = 3` = "solid")) +
    facet_wrap(~ factor, scales = "free_y") +
    coord_cartesian(ylim = c(0, 2.5)) +
    labs(x = "f (length of the coefficients in noise units)", y = NULL, colour = NULL, linetype = NULL,
         title = "The prior as a product: P(total) x P(combination | total), flat R",
         subtitle = "factor 1 is the same for every k; factor 2 makes the peak at zero for k >= 2") +
    theme(legend.position = "bottom"),
  "4_three_factors", width = 10, height = 4)


## ---- 5. Check: the density formula against counted draws (k = 2) --------
## Draw from the prior (the recipe), count how many land in a small square
## around a point, divide by (number of draws x square area), compare with
## the formula (flat R). Correlated predictors, so factor 3 (unit conversion) matters.
predictor_covariance <- matrix(c(1, 0.5, 0.5, 2), 2, 2)
coefficient_draws    <- draw_coefficients_from_prior(4e6, predictor_covariance, prior = "flat_r")
density_formula <- function(b) {
  f_squared <- drop(t(b) %*% predictor_covariance %*% b)
  (1 + f_squared)^(-3 / 2) / (2 * pi * sqrt(f_squared)) * sqrt(det(predictor_covariance))
}
square_half_width <- 0.025
check_points <- expand.grid(b1 = c(-0.6, -0.2, 0.1, 0.4, 0.8), b2 = c(-0.5, -0.15, 0.2, 0.6))
check_points$formula <- apply(check_points, 1, density_formula)
check_points$counted <- apply(check_points[, c("b1", "b2")], 1, function(point)
  mean(abs(coefficient_draws[, 1] - point[1]) < square_half_width &
       abs(coefficient_draws[, 2] - point[2]) < square_half_width) / (2 * square_half_width)^2)
print(check_points)

save_doc_figure(
  ggplot(check_points, aes(formula, counted)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(size = 2.5) +
    labs(x = "density from the formula (3 factors)", y = "density from counting draws",
         title = "The formula and the recipe describe the same prior",
         subtitle = "flat R, k = 2, correlated predictors, unequal spreads;\n20 points, 4 million draws"),
  "5_formula_vs_draws", width = 6, height = 5)


## ---- 6. The three priors among the mixtures of g-priors ------------------
## Under a g-prior, f^2 = (g / n) * chi2_k. R-squared ~ Beta(a, 1) is the
## mixture with g / n = B * W, B ~ Beta(a, k/2 - a), W ~ Inverse-Gamma(1, 1/2),
## for a <= k/2: flat R (a = 1/2), flat partial r (a = k/2, B = 1) and flat
## R-squared (a = 1, k >= 2; with one predictor it is not a g-mixture).
## Implied median R-squared before any data, as k grows (n = 100).
n_for_g <- 100
g_draws <- 400000
g_over_n_draws <- list(
  `Zellner-Siow (scale 1)` = function(k) 1 / rgamma(g_draws, shape = 0.5, rate = 0.5),
  `hyper-g/n (a = 3)`      = function(k) { s <- rbeta(g_draws, 1, 0.5); s / (1 - s) },
  `robust (Bayarri et al.)` = function(k) {
    lower <- (1 + n_for_g) / (k + 1)
    (lower / runif(g_draws)^2 - 1) / n_for_g },
  `intrinsic`              = function(k) (1 / (k + 1)) / rbeta(g_draws, 0.5, 0.5),
  `flat partial r`         = function(k) 1 / (2 * rexp(g_draws)),
  `flat R`                 = function(k)
    (if (k == 1) rep(1, g_draws) else rbeta(g_draws, 0.5, (k - 1) / 2)) / (2 * rexp(g_draws)),
  `flat R-squared`         = function(k) if (k == 1) NULL else
    (if (k == 2) rep(1, g_draws) else rbeta(g_draws, 1, k / 2 - 1)) / (2 * rexp(g_draws)))
g_family_r_squared <- do.call(rbind, lapply(names(g_over_n_draws), function(prior)
  do.call(rbind, lapply(c(1, 2, 3, 5, 10, 20, 30), function(k) {
    g_over_n <- g_over_n_draws[[prior]](k)
    if (is.null(g_over_n)) return(NULL)
    f_squared <- g_over_n * rchisq(g_draws, df = k)
    data.frame(prior = prior, number_of_predictors = k,
               median_r_squared = median(f_squared / (1 + f_squared)),
               ks_to_beta_half_one = suppressWarnings(unname(
                 ks.test(f_squared / (1 + f_squared), "pbeta", 0.5, 1)$statistic)))
  }))))
print(g_family_r_squared)

save_doc_figure(
  ggplot(g_family_r_squared, aes(number_of_predictors, median_r_squared, colour = prior)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_x_continuous(breaks = c(1, 2, 3, 5, 10, 20, 30)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = c(`Zellner-Siow (scale 1)` = "grey50",
                                   `hyper-g/n (a = 3)` = "#e7298a",
                                   `robust (Bayarri et al.)` = "#7570b3",
                                   `intrinsic` = "#A6761D", prior_colours)) +
    labs(x = "number of predictors (k)", y = "median implied R-squared, before any data",
         colour = NULL,
         title = "Mixtures of g-priors: only a prior on the total keeps the same claim for every k",
         subtitle = "n = 100; flat R and flat R-squared stay put, robust and intrinsic scale g with k but still drift") +
    theme(legend.position = "bottom"),
  "6_g_prior_family")
