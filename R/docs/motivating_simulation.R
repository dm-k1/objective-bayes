## =====================================================================
##  The motivating simulation: what happens as the number of predictors grows
##
##  The starting question: with one predictor, a flat correlation r
##  (Jeffreys) is the obvious default. What happens to the correlations when
##  predictors are added? This simulation answers that without choosing a
##  prior: it draws correlation matrices of (y, x1 ... xk) with every valid
##  matrix equally likely, and looks at what k does. Only then does it show
##  the three intuitive ways of deciding what should stay flat, which are the
##  three priors the validation (R/studies/) compares. It explains the
##  project and validates nothing, so it lives here with the code behind the
##  docs.
##
##  Part 1. What k does. Every valid correlation matrix of (y, x1 ... xk)
##    equally likely (LKJ with eta = 1), k = 1 ... 10, 15, 20, 30, 40, 50.
##    From each matrix: the correlation of y with x1 (marginal r), the partial
##    correlation of y with x1 given the other predictors (partial r), and R
##    and R^2 of y on all k predictors. At k = 1 these are one quantity
##    (partial r = r, R = |r|, R^2 = r^2); as k grows they come apart.
##
##  Part 2. Three choices of what stays flat. Each keeps one quantity flat at
##    every k, and the others move:
##    partial r flat: the matrices of part 1 as they are (prior flat_partial_r),
##    R flat:         the same matrices, rescaled so that R is flat (flat_r),
##    R-squared flat: the same matrices, rescaled so that R^2 is flat
##                    (flat_r_squared).
##    The rescaling keeps the predictor correlations and the direction of y's
##    correlations with the predictors, and only changes their size, so the
##    matrix stays valid (R < 1) and the three differ only in what is flat.
##
##  Checks (KS distance against the exact distribution at every k, judged at
##  the .01 level split over the values of k of a check, Bonferroni):
##    M0. the fast sampler used here (the correlation matrix of a Wishart draw)
##        gives the same matrices as the vine method of Lewandowski, Kurowicka
##        and Joe (2009), at k = 1, 3, 8
##    M1. the marginal r is Beta(d/2, d/2) on (-1, 1), d = k + 1 variables
##    M2. the partial r is flat on (-1, 1)
##    M3. R^2 is Beta(k/2, 1)
##    M4. R is flat in the "R flat" version, R^2 in the "R-squared flat" one
##
##  Run from the project folder: Rscript R/docs/motivating_simulation.R
##  About 8 minutes on 8 cores; no Stan needed. Writes its figures and its
##  table to docs/figures/ (motivation_*.png, motivation_table.md).
## =====================================================================
library(ggplot2)
library(patchwork)
source("R/shared/settings.R")             # parallel_map(), save_figure(), record_check()
source("R/shared/comparison_methods.R")   # draw_lkj_correlation_matrix()
figures_folder <- file.path("docs", "figures")   # save_figure() writes here, next to the docs' figures

set.seed(2026)
number_of_draws       <- 200000
numbers_of_predictors <- c(1:10, 15, 20, 30, 40, 50)
quantities <- c(marginal_r = "marginal r of y with x1", partial_r = "partial r of y with x1, given the others",
                R = "R", r_squared = "R-squared")
versions <- c("partial r flat", "R flat", "R-squared flat")
version_colours <- setNames(c("#0072B2", "#D55E00", "#009E73"), versions)   # Okabe-Ito


## ---- Random valid correlation matrices --------------------------------------
## Every valid correlation matrix of d variables equally likely: the
## correlation matrix of a Wishart draw with d + 1 degrees of freedom and
## identity scale has density proportional to det(C)^((nu - d - 1) / 2), which
## is flat for nu = d + 1. rWishart() is compiled code, so k = 50 stays cheap;
## check M0 compares it with the vine method (draw_lkj_correlation_matrix()).
draw_uniform_correlation_matrix <- function(dimension) {
  scatter <- stats::rWishart(1, dimension + 1, diag(dimension))[, , 1]
  stats::cov2cor(scatter)
}

## Rescale y's correlations with the predictors so that R becomes new_R.
## R^2 = rho' solve(Sigma_xx) rho, so multiplying rho by c multiplies R by c.
with_total_R <- function(correlations, new_R) {
  current_R <- sqrt(1 - 1 / solve(correlations)[1, 1])
  correlations[1, -1] <- correlations[-1, 1] <- correlations[1, -1] * new_R / current_R
  correlations
}

## Marginal r, partial r of y with x1 given the rest, and R^2, from one
## matrix. Variable 1 is y, variables 2 ... k + 1 are the predictors.
##   R^2       = 1 - 1 / (inverse correlation matrix)[y, y]
##   partial r = -inverse[y, x1] / sqrt(inverse[y, y] * inverse[x1, x1])
summarise_matrix <- function(correlations) {
  inverse <- solve(correlations)
  c(correlations[1, 2], -inverse[1, 2] / sqrt(inverse[1, 1] * inverse[2, 2]), 1 - 1 / inverse[1, 1])
}


## ---- Simulate ------------------------------------------------------------
## One matrix per draw, used for all three versions; in chunks of draws, so
## each parallel task is worth its overhead
chunks <- split(seq_len(number_of_draws), ceiling(seq_len(number_of_draws) / 5000))
draws <- do.call(rbind, lapply(numbers_of_predictors, function(k) {
  values <- do.call(rbind, parallel_map(chunks, function(chunk) t(vapply(chunk, function(draw) {
    correlations <- draw_uniform_correlation_matrix(k + 1)
    c(summarise_matrix(correlations),
      summarise_matrix(with_total_R(correlations, runif(1))),
      summarise_matrix(with_total_R(correlations, sqrt(runif(1)))))
  }, numeric(9)))))
  do.call(rbind, lapply(seq_along(versions), function(v) {
    data.frame(k = k, version = versions[v], marginal_r = values[, 3 * v - 2],
               partial_r = values[, 3 * v - 1], r_squared = values[, 3 * v],
               R = sqrt(values[, 3 * v]))
  }))
}))
draws$version <- factor(draws$version, levels = versions)
uniform <- subset(draws, version == "partial r flat")   # part 1: the matrices as drawn


## ---- Checks ------------------------------------------------------------------------
## The KS critical distance at level alpha (asymptotic): sqrt(-log(alpha / 2) / 2) / sqrt(n)
critical_for <- function(tests) sqrt(-log(0.01 / tests / 2) / 2) / sqrt(number_of_draws)
critical_distance <- critical_for(length(numbers_of_predictors))
ks_distance <- function(values, cdf, ...) suppressWarnings(ks.test(values, cdf, ...)$statistic)
per_k <- function(data, statistic) vapply(numbers_of_predictors, function(k)
  statistic(data[data$k == k, ], k), numeric(1))

## M0. the Wishart route against the vine method, on the quantities used here
## (a two-sample KS test, 20,000 vine matrices per k)
m0 <- sapply(c(1, 3, 8), function(k) {
  vine <- do.call(rbind, parallel_map(split(1:20000, rep(1:4, each = 5000)), function(chunk)
    t(vapply(chunk, function(draw) summarise_matrix(draw_lkj_correlation_matrix(k + 1)), numeric(3)))))
  wishart <- uniform[uniform$k == k, c("marginal_r", "partial_r", "r_squared")]
  min(vapply(1:3, function(j) suppressWarnings(ks.test(vine[, j], wishart[[j]])$p.value), numeric(1)))
})
record_check("M0. Wishart route = vine method (k = 1, 3, 8; marginal r, partial r, R^2)",
             min(m0) > 0.01 / 9, sprintf("smallest two-sample KS p = %.3f (Bonferroni over 9: .0011)", min(m0)))

m1 <- per_k(uniform, function(d, k) ks_distance((d$marginal_r + 1) / 2, "pbeta", (k + 1) / 2, (k + 1) / 2))
m2 <- per_k(uniform, function(d, k) ks_distance(d$partial_r, "punif", -1, 1))
m3 <- per_k(uniform, function(d, k) ks_distance(d$r_squared, "pbeta", k / 2, 1))
m4 <- c(per_k(subset(draws, version == "R flat"), function(d, k) ks_distance(d$R, "punif")),
        per_k(subset(draws, version == "R-squared flat"), function(d, k) ks_distance(d$r_squared, "punif")))
record_check("M1. marginal r ~ Beta(d/2, d/2) on (-1, 1), every k", max(m1) < critical_distance,
             sprintf("largest KS distance %.4f; critical %.4f", max(m1), critical_distance))
record_check("M2. partial r flat on (-1, 1), every k", max(m2) < critical_distance,
             sprintf("largest KS distance %.4f; critical %.4f", max(m2), critical_distance))
record_check("M3. R-squared ~ Beta(k/2, 1), every k", max(m3) < critical_distance,
             sprintf("largest KS distance %.4f; critical %.4f", max(m3), critical_distance))
record_check("M4. R flat (R flat), R-squared flat (R-squared flat), every k",
             max(m4) < critical_for(2 * length(numbers_of_predictors)),
             sprintf("largest KS distance %.4f; critical %.4f", max(m4), critical_for(2 * length(numbers_of_predictors))))


## ---- Table ---------------------------------------------------------------------
summary_table <- do.call(rbind, lapply(split(draws, list(draws$version, draws$k)), function(d) data.frame(
  version              = d$version[1],
  k                    = d$k[1],
  P_abs_marginal_r_0.3 = mean(abs(d$marginal_r) > 0.3),
  P_abs_partial_r_0.3  = mean(abs(d$partial_r) > 0.3),
  median_R2            = median(d$r_squared),
  P_R2_above_0.5       = mean(d$r_squared > 0.5))))
summary_table <- summary_table[order(summary_table$version, summary_table$k), ]
rownames(summary_table) <- NULL
cat("\nRandom valid correlation matrices (", number_of_draws, " per k)\n", sep = "")
for (v in versions) {
  cat("\n", v, "\n", sep = "")
  print(format(subset(summary_table, version == v, -version), digits = 3), row.names = FALSE)
}

share <- function(v, column) sprintf("%.2f", summary_table[summary_table$version == v, column])
table_k <- summary_table$k[summary_table$version == versions[1]]
writeLines(c(
  "## The motivating simulation",
  "",
  sprintf("%s correlation matrices of (y, x1 ... xk) per k.",
          format(number_of_draws, big.mark = ",", scientific = FALSE)),
  "",
  "**Part 1. What k does**, with every valid matrix equally likely: the marginal r",
  "narrows towards 0, the partial r stays flat, and R-squared grows towards 1.",
  "",
  "| k | P(\\|marginal r\\| > .3) | P(\\|partial r\\| > .3) | median R² | P(R² > .5) |",
  "|---|---|---|---|---|",
  sprintf("| %d | %s | %s | %s | %s |", table_k, share(versions[1], "P_abs_marginal_r_0.3"),
          share(versions[1], "P_abs_partial_r_0.3"), share(versions[1], "median_R2"),
          share(versions[1], "P_R2_above_0.5")),
  "",
  "**Part 2. Three choices of what stays flat.** Keep the partial r flat and R-squared",
  "grows with k; keep R or R-squared flat and the partial r narrows with k.",
  "",
  "| k | P(\\|partial r\\| > .3), partial r flat | R flat | R-squared flat | P(R² > .5), partial r flat | R flat | R-squared flat |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %d | %s | %s | %s | %s | %s | %s |", table_k,
          share(versions[1], "P_abs_partial_r_0.3"), share(versions[2], "P_abs_partial_r_0.3"),
          share(versions[3], "P_abs_partial_r_0.3"), share(versions[1], "P_R2_above_0.5"),
          share(versions[2], "P_R2_above_0.5"), share(versions[3], "P_R2_above_0.5"))),
  file.path(figures_folder, "motivation_table.md"))


## ---- Histogram densities, drawn as a line through the bin centres ----------------
densities <- function(data, column, range, group) {
  breaks <- seq(range[1], range[2], length.out = 41)
  do.call(rbind, lapply(split(data, data[group], drop = TRUE), function(d) {
    histogram <- hist(d[[column]], breaks = breaks, plot = FALSE)
    cbind(d[1, group, drop = FALSE], value = histogram$mids, density = histogram$density)
  }))
}
ranges <- list(marginal_r = c(-1, 1), partial_r = c(-1, 1), R = c(0, 1), r_squared = c(0, 1))


## ---- Figure 1: what k does to each quantity (part 1) --------------------------------
## One panel per quantity, one curve per k, light (k = 1) to dark (k = 50)
shown_k <- c(1, 2, 5, 10, 20, 50)
k_colours <- setNames(c("#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b"), shown_k)
part_1_panel <- function(column) {
  bins <- densities(subset(uniform, k %in% shown_k), column, ranges[[column]], "k")
  bins$k <- factor(bins$k, levels = shown_k)
  ggplot(bins, aes(value, density, colour = k)) +
    geom_line(linewidth = 0.8) +
    coord_cartesian(xlim = ranges[[column]], ylim = c(0, 6)) +
    scale_colour_manual(values = k_colours, name = "k") +
    labs(x = quantities[[column]], y = "density") +
    theme_minimal(base_size = 11)
}
save_figure(
  wrap_plots(lapply(names(quantities), part_1_panel), nrow = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Every valid correlation matrix equally likely: what adding predictors does",
      subtitle = paste0(format(number_of_draws, big.mark = ",", scientific = FALSE),
                        " matrices of (y, x1 ... xk) per k. At k = 1 the four are one quantity; ",
                        "as k grows the marginal r narrows, the partial r stays flat and R and R-squared grow")) &
    theme(legend.position = "bottom"),
  "motivation_what_k_does", width = 14, height = 4.5)


## ---- Figure 2: the three choices overlaid (part 2) ----------------------------------
## One row per quantity, each on its own range; a few k to keep it readable
version_panel <- function(column) {
  bins <- densities(subset(draws, k %in% c(1, 2, 5, 10, 50)), column, ranges[[column]], c("version", "k"))
  ggplot(bins, aes(value, density, colour = version)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ k, nrow = 1, labeller = labeller(k = function(k) paste("k =", k))) +
    coord_cartesian(xlim = ranges[[column]], ylim = c(0, 4)) +
    scale_colour_manual(values = version_colours) +
    labs(x = quantities[[column]], y = "density", colour = NULL) +
    theme_minimal(base_size = 11)
}
save_figure(
  (version_panel("partial_r") / version_panel("R") / version_panel("r_squared")) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Three choices of what stays flat: the others move as k grows",
      subtitle = "The same matrices in every version; only the size of y's correlations differs") &
    theme(legend.position = "bottom"),
  "motivation_partial_r_R_R2_by_k", width = 13, height = 9)


## ---- Figure 3: the three choices against k (part 2) ---------------------------------
trend <- rbind(
  data.frame(summary_table[c("version", "k")], measure = "P(|partial r| > .3)",
             value = summary_table$P_abs_partial_r_0.3),
  data.frame(summary_table[c("version", "k")], measure = "P(R-squared > .5)",
             value = summary_table$P_R2_above_0.5))
save_figure(
  ggplot(trend, aes(k, value, colour = version)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ measure) +
    scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = version_colours) +
    labs(x = "number of predictors (k, log scale)", y = "share of correlation matrices", colour = NULL,
         title = "Keep the partial r flat and R-squared grows; keep R or R-squared flat and the partial r narrows") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom"),
  "motivation_three_versions_against_k", width = 11, height = 5)

report_checks()
