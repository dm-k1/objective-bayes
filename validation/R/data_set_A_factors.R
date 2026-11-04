## =====================================================================
##  Data set A, extended with a factor (validation/ADEMP.md, study 8)
##
##  A factor with q + 1 equal-sized levels (q columns in the model), equally
##  spaced level means with Cohen's f^2 = 0 (null), .02 (small) or .15
##  (medium), plus 2 or 10 uncorrelated continuous predictors. What the
##  others explain is set afterwards from the same predictors, factor and
##  noise (validation/08_factors.R), so the values of it differ only in the
##  others' slopes. Every dataset has its own stream (validation/data_sets.R,
##  with_stream_of()), with an ID that starts with "A-factor".
## =====================================================================

design_A_factors <- function() {
  expand.grid(q = 2:5, n = c(50, 200, 1000), other_predictors = c(2, 10),
              factor_effect = c("null", "small", "medium"), stringsAsFactors = FALSE)
}

cohen_f_squared_of_factor <- c(null = 0, small = 0.02, medium = 0.15)

## What the other predictors explain of the outcome (when the factor is null)
explained_by_others_values <- c(0, 0.5, 0.9)

cell_id_A_factors <- function(condition)
  with(condition, sprintf("A-factor-q%d-n%d-others%d-%s", q, n, other_predictors, factor_effect))

## Level means: equally spaced and centred, with mean square f^2 (equal-sized
## levels), so the factor's variance is f^2 times the noise variance
factor_level_means <- function(q, f_squared) {
  spacing <- seq_len(q + 1) - (q + 2) / 2
  spacing * sqrt(f_squared / mean(spacing^2))
}

## One dataset of the extended set A: replication r of the cell in condition.
## Returns the factor, the other predictors and the noise; the outcome follows
## from them for each value of what the others explain (outcome_A_factors()).
generate_A_factors <- function(condition, replication) {
  n <- condition$n; p <- condition$other_predictors
  with_stream_of(paste(cell_id_A_factors(condition), replication, sep = "-r"), {
    group  <- factor(sample(rep(seq_len(condition$q + 1), length.out = n)))
    others <- matrix(rnorm(n * p), n, dimnames = list(NULL, paste0("x", seq_len(p))))
    noise  <- rnorm(n)
    list(group = group, others = others, noise = noise)
  })
}

## The outcome of a dataset when the others explain R-squared = r_squared_others
outcome_A_factors <- function(dataset, condition, r_squared_others) {
  p <- condition$other_predictors
  slopes <- rep(sqrt(r_squared_others / (1 - r_squared_others) / p), p)
  drop(dataset$others %*% slopes) +
    factor_level_means(condition$q, cohen_f_squared_of_factor[[condition$factor_effect]])[dataset$group] +
    dataset$noise
}
