## =====================================================================
##  Validation 00: the first result, the one that started the comparison
##
##  The starting question: what happens to partial r, R and R-squared as the
##  number of predictors k grows? With one predictor, a flat r (Jeffreys) is
##  all there is. With more, not everything can stay flat. This simulation
##  shows the fork, and studies 01-11 test the three options it leaves.
##
##  Three versions, k = 1 ... 10, all of them valid correlation matrices of
##  (y, x1 ... xk):
##    A. partial r flat: every valid correlation matrix equally likely,
##       which leaves each partial correlation flat on (-1, 1),
##    B. R flat:         the same matrices, with R rescaled to be flat on (0, 1),
##    C. R-squared flat: the same matrices, with R^2 rescaled to be flat on (0, 1).
##  These are the package's priors flat_partial_r, flat_r and flat_r_squared.
##
##  B and C keep the predictor correlations and the direction of y's
##  correlations with the predictors from A, and only change how large y's
##  correlations are, so that R takes the drawn value. Rescaling keeps the
##  matrix valid as long as R < 1, so the three versions differ only in what
##  is kept flat.
##
##  From each matrix: the partial correlation of y with x1 given the other
##  predictors, and R and R^2 of y on all k predictors.
##
##  What it shows: keep each partial r flat (A) and R-squared grows with k;
##  keep the total flat (B or C) and each partial r shrinks with k. With one
##  predictor A and B are the same.
##
##    F1-F3. each version is flat where it should be: A's partial r, B's R and
##        C's R^2 against the flat distribution, for every k (KS distance < 0.01)
##
##  Run from the project folder: Rscript validation/00_initial_simulation.R
##  About a minute on 8 cores; no Stan needed. Writes its table to
##  validation/results/00_first_result.md (included in run_summary.md).
## =====================================================================
library(ggplot2)
library(patchwork)
source("validation/settings.R")   # parallel_map(), save_figure(), record_check()
source("validation/comparison_methods.R")   # draw_lkj_correlation_matrix()

set.seed(2026)
number_of_draws       <- 200000
numbers_of_predictors <- 1:10
versions <- c("A. partial r flat", "B. R flat", "C. R-squared flat")
version_colours <- setNames(c("#0072B2", "#D55E00", "#009E73"), versions)


## ---- Random valid correlation matrices --------------------------------------
## draw_lkj_correlation_matrix() (validation/comparison_methods.R) with eta = 1:
## the vine method of Lewandowski, Kurowicka and Joe (2009), every valid
## correlation matrix equally likely, every result positive definite.

## Rescale y's correlations with the predictors so that R becomes new_R.
## R^2 = rho' solve(Sigma_xx) rho, so multiplying rho by c multiplies R by c.
with_total_R <- function(correlations, new_R) {
  current_R <- sqrt(1 - 1 / solve(correlations)[1, 1])
  correlations[1, -1] <- correlations[-1, 1] <- correlations[1, -1] * new_R / current_R
  correlations
}

## Partial r of y with x1 given the rest, and R^2, from one matrix.
## Variable 1 is y, variables 2 ... k + 1 are the predictors.
##   R^2       = 1 - 1 / (inverse correlation matrix)[y, y]
##   partial r = -inverse[y, x1] / sqrt(inverse[y, y] * inverse[x1, x1])
summarise_matrix <- function(correlations) {
  inverse <- solve(correlations)
  c(-inverse[1, 2] / sqrt(inverse[1, 1] * inverse[2, 2]), 1 - 1 / inverse[1, 1])
}


## ---- Simulate ------------------------------------------------------------
## One matrix per draw, used for all three versions.
draws <- do.call(rbind, lapply(numbers_of_predictors, function(k) {
  one_matrix <- function(draw) {
    correlations <- draw_lkj_correlation_matrix(k + 1)
    c(summarise_matrix(correlations),
      summarise_matrix(with_total_R(correlations, runif(1))),
      summarise_matrix(with_total_R(correlations, sqrt(runif(1)))))
  }
  values <- do.call(rbind, parallel_map(seq_len(number_of_draws), one_matrix))
  do.call(rbind, lapply(seq_along(versions), function(v) {
    data.frame(k = k, version = versions[v],
               partial_r = values[, 2 * v - 1], r_squared = values[, 2 * v],
               R = sqrt(values[, 2 * v]))
  }))
}))
draws$version <- factor(draws$version, levels = versions)


## ---- Table: partial r and R^2 as k grows, per version --------------------------
summary_table <- do.call(rbind, lapply(split(draws, list(draws$version, draws$k)), function(d) {
  data.frame(
    version              = d$version[1],
    k                    = d$k[1],
    sd_partial_r         = sd(d$partial_r),
    P_abs_partial_r_0.3  = mean(abs(d$partial_r) > 0.3),
    median_R2            = median(d$r_squared),
    P_R2_above_0.5       = mean(d$r_squared > 0.5),
    P_R2_above_0.8       = mean(d$r_squared > 0.8))
}))
summary_table <- summary_table[order(summary_table$version, summary_table$k), ]
cat("\nRandom valid correlation matrices (", number_of_draws, " per k)\n", sep = "")
for (v in versions) {
  cat("\n", v, "\n", sep = "")
  print(format(subset(summary_table, version == v, -version), digits = 3), row.names = FALSE)
}

## The first result, as a table for run_summary.md: per k, how much room one
## coefficient gets and how large the total is claimed to be
first_result <- reshape(summary_table[c("version", "k", "P_abs_partial_r_0.3", "P_R2_above_0.5")],
                        idvar = "k", timevar = "version", direction = "wide")
writeLines(c(
  "## 00: the first result",
  "",
  sprintf("%s valid correlation matrices of (y, x1 ... xk) per k. Keep each partial r flat",
          format(number_of_draws, big.mark = ",", scientific = FALSE)),
  "and R-squared grows with k; keep R or R-squared flat and each partial r shrinks with k.",
  "",
  paste0("| k | ", paste(rep(c("P(\\|partial r\\| > 0.3)", "P(R² > 0.5)"), 3), collapse = " | "), " |"),
  paste0("| | ", paste(rep(sub("^[ABC]\\. ", "", versions), each = 2), collapse = " | "), " |"),
  paste0(paste(rep("|---", 7), collapse = ""), "|"),
  apply(first_result, 1, function(row)
    paste0("| ", row[["k"]], " | ", paste(sprintf("%.2f", as.numeric(row[c(
      paste0(c("P_abs_partial_r_0.3.", "P_R2_above_0.5."), versions[1]),
      paste0(c("P_abs_partial_r_0.3.", "P_R2_above_0.5."), versions[2]),
      paste0(c("P_abs_partial_r_0.3.", "P_R2_above_0.5."), versions[3]))])), collapse = " | "), " |"))),
  file.path(results_folder, "00_first_result.md"))


## ---- Checks: each version is flat where it should be ----------------------------
flat_distance <- function(version, column, lower) sapply(numbers_of_predictors, function(k)
  suppressWarnings(ks.test(draws[[column]][draws$version == version & draws$k == k],
                           "punif", lower, 1)$statistic))
f1 <- flat_distance(versions[1], "partial_r", -1)
f2 <- flat_distance(versions[2], "R", 0)
f3 <- flat_distance(versions[3], "r_squared", 0)
record_check("F1. A: partial r flat on (-1, 1), every k", max(f1) < 0.01,
             sprintf("largest KS distance %.4f", max(f1)))
record_check("F2. B: R flat on (0, 1), every k", max(f2) < 0.01,
             sprintf("largest KS distance %.4f", max(f2)))
record_check("F3. C: R-squared flat on (0, 1), every k", max(f3) < 0.01,
             sprintf("largest KS distance %.4f", max(f3)))


## ---- Figure 1: partial r, R and R^2, the three versions overlaid ---------------
## One row per quantity, each on its own range; a few k to keep it readable
## (figure 2 has every k). Outlines instead of filled bars, so the three
## versions stay visible where they overlap.
shown_k <- c(1, 2, 3, 5, 10)

## Histogram densities per version and k, drawn as a line through the bin
## centres (geom_freqpoly would add a drop to zero at both ends).
histogram_row <- function(column, quantity_label, range) {
  breaks <- seq(range[1], range[2], length.out = 41)
  bins <- do.call(rbind, lapply(split(subset(draws, k %in% shown_k), ~ version + k, drop = TRUE),
                                function(d) {
    histogram <- hist(d[[column]], breaks = breaks, plot = FALSE)
    data.frame(version = d$version[1], k = d$k[1], value = histogram$mids, density = histogram$density)
  }))
  ggplot(bins, aes(value, density, colour = version)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ k, nrow = 1, labeller = labeller(k = function(k) paste("k =", k))) +
    coord_cartesian(xlim = range, ylim = c(0, 4)) +
    scale_colour_manual(values = version_colours) +
    labs(x = quantity_label, y = "density", colour = NULL) +
    theme_minimal(base_size = 11)
}

save_figure(
  (histogram_row("partial_r", "partial r of y with x1, given the other predictors", c(-1, 1)) /
     histogram_row("R",         "R", c(0, 1)) /
     histogram_row("r_squared", "R-squared", c(0, 1))) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Valid correlation matrices: not everything can stay flat as k grows",
      subtitle = paste(format(number_of_draws, big.mark = ",", scientific = FALSE), "matrices of (y, x1 ... xk) per k; the same matrices in every version,",
                       "only the size of y's correlations differs")) &
    theme(legend.position = "bottom"),
  "00_partial_r_R_R2_by_k", width = 13, height = 9)


## ---- Figure 2: the three versions against k ---------------------------------
trend <- rbind(
  data.frame(summary_table[c("version", "k")], measure = "P(|partial r| > 0.3)",
             value = summary_table$P_abs_partial_r_0.3),
  data.frame(summary_table[c("version", "k")], measure = "P(R-squared > 0.5)",
             value = summary_table$P_R2_above_0.5))

save_figure(
  ggplot(trend, aes(k, value, colour = version)) +
    geom_line(linewidth = 0.9) +
    geom_point() +
    facet_wrap(~ measure) +
    scale_x_continuous(breaks = numbers_of_predictors) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = version_colours) +
    labs(x = "number of predictors (k)", y = "share of correlation matrices", colour = NULL,
         title = "Keep partial r flat and R-squared grows; keep R or R-squared flat and partial r shrinks") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom"),
  "00_three_versions_against_k", width = 11, height = 5)

report_checks()
