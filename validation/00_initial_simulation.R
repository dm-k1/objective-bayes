## =====================================================================
##  Validation 00: an initial simulation of partial r, R and R-squared
##
##  The starting question. With one predictor, a flat prior on r (Jeffreys)
##  is the same as drawing a valid 2 x 2 correlation matrix of (y, x) at
##  random, and it gives R^2 ~ Beta(1/2, 1). Does that stay true with more
##  predictors?
##
##  Simulation: for k = 1 ... 10, draw valid (k + 1) x (k + 1) correlation
##  matrices of (y, x1 ... xk), every valid matrix equally likely. From each
##  matrix compute
##    - the partial correlation of y with x1, given the other predictors,
##    - R and R^2 of y on all k predictors.
##  Then compare R^2 with Beta(1/2, 1), the flat-on-R answer of k = 1.
##
##  What it shows: each partial r stays flat on (-1, 1) for every k, but R^2
##  moves away from Beta(1/2, 1) as k grows, towards high values. So once
##  there is more than one predictor, a flat prior on R no longer describes
##  "any valid correlation matrix equally likely".
##
##  Run from the project folder: Rscript validation/00_initial_simulation.R
##  Under a minute; no Stan needed.
## =====================================================================
library(ggplot2)
source("validation/settings.R")   # save_figure()

set.seed(2026)
number_of_draws       <- 20000
numbers_of_predictors <- 1:10


## ---- Random valid correlation matrices --------------------------------------
## The vine method of Lewandowski, Kurowicka and Joe (2009) with eta = 1:
## every valid d x d correlation matrix is equally likely. It builds the
## matrix from partial correlations, so every result is positive definite.
draw_correlation_matrix <- function(d) {
  correlations <- diag(d)
  partials     <- matrix(0, d, d)
  beta_shape   <- 1 + (d - 1) / 2
  for (level in seq_len(d - 1)) {
    beta_shape <- beta_shape - 1 / 2
    for (i in (level + 1):d) {
      partials[level, i] <- 2 * rbeta(1, beta_shape, beta_shape) - 1
      correlation <- partials[level, i]
      if (level > 1) for (l in (level - 1):1)
        correlation <- correlation * sqrt((1 - partials[l, i]^2) * (1 - partials[l, level]^2)) +
          partials[l, i] * partials[l, level]
      correlations[level, i] <- correlations[i, level] <- correlation
    }
  }
  correlations
}


## ---- Simulate ------------------------------------------------------------
## Variable 1 is y, variables 2 ... k + 1 are the predictors.
##   R^2       = 1 - 1 / (inverse correlation matrix)[y, y]
##   partial r = -inverse[y, x1] / sqrt(inverse[y, y] * inverse[x1, x1])
draws <- do.call(rbind, lapply(numbers_of_predictors, function(k) {
  one_matrix <- function(draw) {
    inverse <- solve(draw_correlation_matrix(k + 1))
    c(partial_r = -inverse[1, 2] / sqrt(inverse[1, 1] * inverse[2, 2]),
      r_squared = 1 - 1 / inverse[1, 1])
  }
  values <- do.call(rbind, parallel_map(seq_len(number_of_draws), one_matrix))
  data.frame(k = k, partial_r = values[, "partial_r"],
             R = sqrt(values[, "r_squared"]), r_squared = values[, "r_squared"])
}))


## ---- Table: simulated R^2 against Beta(1/2, 1) --------------------------------
summary_table <- do.call(rbind, lapply(split(draws, draws$k), function(d) {
  data.frame(
    k                     = d$k[1],
    median_R2             = median(d$r_squared),
    median_R2_beta_half   = qbeta(0.5, 1 / 2, 1),
    P_R2_above_0.5        = mean(d$r_squared > 0.5),
    P_0.5_beta_half       = 1 - pbeta(0.5, 1 / 2, 1),
    P_R2_above_0.8        = mean(d$r_squared > 0.8),
    P_0.8_beta_half       = 1 - pbeta(0.8, 1 / 2, 1),
    sd_partial_r          = sd(d$partial_r))   # flat on (-1, 1): 0.577
}))
cat("\nRandom valid correlation matrices (", number_of_draws, " per k)",
    " against Beta(1/2, 1)\n", sep = "")
print(format(summary_table, digits = 3), row.names = FALSE)


## ---- Figure 1: partial r, R and R^2 for k = 1 ... 10 --------------------------
## Black curves: what Beta(1/2, 1) on R^2 (flat on R) would give.
long_draws <- rbind(
  data.frame(k = draws$k, quantity = "partial r (y, x1 | rest)", value = draws$partial_r),
  data.frame(k = draws$k, quantity = "R",                        value = draws$R),
  data.frame(k = draws$k, quantity = "R-squared",                value = draws$r_squared))
quantity_levels <- c("partial r (y, x1 | rest)", "R", "R-squared")
long_draws$quantity <- factor(long_draws$quantity, levels = quantity_levels)

grid_values <- seq(0.005, 0.995, by = 0.005)
flat_r_curves <- do.call(rbind, lapply(numbers_of_predictors, function(k) rbind(
  data.frame(k = k, quantity = "partial r (y, x1 | rest)", value = c(-1, 1), density = 0.5),
  data.frame(k = k, quantity = "R",         value = grid_values, density = 1),
  data.frame(k = k, quantity = "R-squared", value = grid_values,
             density = dbeta(grid_values, 1 / 2, 1)))))
flat_r_curves$quantity <- factor(flat_r_curves$quantity, levels = quantity_levels)

save_figure(
  ggplot(long_draws, aes(value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = "grey70", colour = "white", boundary = 0) +
    geom_line(data = flat_r_curves, aes(value, density), linewidth = 0.7) +
    facet_grid(quantity ~ k, scales = "free_y", labeller = labeller(k = function(k) paste("k =", k))) +
    coord_cartesian(ylim = c(0, 4)) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    labs(x = NULL, y = "density",
         title = "Random valid correlation matrices: partial r stays flat, R-squared leaves Beta(1/2, 1)",
         subtitle = paste("grey: simulated; black: flat on R (R-squared ~ Beta(1/2, 1)),",
                          "which is exact for k = 1 only")) +
    theme_minimal(base_size = 11),
  "00_partial_r_R_R2_by_k", width = 16, height = 7)


## ---- Figure 2: chance of a high R^2 against k ----------------------------------
chance_of_high_r_squared <- rbind(
  data.frame(k = summary_table$k, threshold = "R-squared > 0.5", source = "simulated",
             probability = summary_table$P_R2_above_0.5),
  data.frame(k = summary_table$k, threshold = "R-squared > 0.8", source = "simulated",
             probability = summary_table$P_R2_above_0.8),
  data.frame(k = summary_table$k, threshold = "R-squared > 0.5", source = "Beta(1/2, 1)",
             probability = summary_table$P_0.5_beta_half),
  data.frame(k = summary_table$k, threshold = "R-squared > 0.8", source = "Beta(1/2, 1)",
             probability = summary_table$P_0.8_beta_half))

save_figure(
  ggplot(chance_of_high_r_squared, aes(k, probability, colour = threshold, linetype = source)) +
    geom_line(linewidth = 0.9) +
    geom_point(data = subset(chance_of_high_r_squared, source == "simulated")) +
    scale_x_continuous(breaks = numbers_of_predictors) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = c("R-squared > 0.5" = "black", "R-squared > 0.8" = "#d95f02")) +
    scale_linetype_manual(values = c("simulated" = "solid", "Beta(1/2, 1)" = "dashed")) +
    labs(x = "number of predictors (k)", y = "probability", colour = NULL, linetype = NULL,
         title = "With random correlation matrices, a high R-squared gets more likely as k grows",
         subtitle = "dashed: Beta(1/2, 1), the same for every k") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", legend.box = "vertical"),
  "00_chance_of_high_R2", width = 8, height = 5)
