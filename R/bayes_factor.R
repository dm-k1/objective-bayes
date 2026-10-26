## =====================================================================
##  Bayes factors (closed form), for the two priors
## =====================================================================
## Two ways to extend Jeffreys' flat correlation to k predictors
## (docs/0_overview.md):
##   flat-R                R flat for every k: R^2 ~ Beta(1/2, 1)
##   uniform correlations  every correlation matrix of (y, x1 ... xk) equally
##                         likely: R^2 ~ Beta(k/2, 1) (stan_lm with eta = 1)
## Both take a uniform direction for how the total is divided. Every model
## gets an exact Bayes factor against the intercept-only model; any two models
## then compare by dividing. The data enter only through the sample R-squared,
## n and k.

log_sum_exp <- function(log_values) {
  largest <- max(log_values)
  largest + log(sum(exp(log_values - largest)))
}
log_sum_exp_per_row <- function(log_matrix) {
  largest_per_row <- apply(log_matrix, 1, max)
  largest_per_row + log(rowSums(exp(log_matrix - largest_per_row)))
}

## Sample R-squared (always with intercept). Pure arithmetic, no fitting.
## A column that is a combination of the others (aliased) adds nothing to
## R-squared but would still count in k, so it stops instead.
sample_r_squared <- function(predictor_columns, outcome) {
  if (is.null(predictor_columns) || NCOL(predictor_columns) == 0) return(0)
  decomposition <- qr(cbind(1, predictor_columns))
  if (decomposition$rank < ncol(decomposition$qr))
    stop("The predictor columns are linearly dependent (one is a combination of ",
         "the others, or constant). Remove the redundant columns.")
  residuals_left <- qr.resid(decomposition, outcome)
  1 - sum(residuals_left^2) / sum((outcome - mean(outcome))^2)
}

#' Log Bayes factor of a model with k predictors against the intercept-only model
#'
#' Closed form: the one-predictor Bayes factor averaged over every direction
#' the effect could take. The data enter only through R-squared, n and k.
#' @param explained_share Sample R-squared of the model.
#' @param sample_size Number of observations, n.
#' @param number_of_predictors Number of predictor columns, k.
#' @param grid_points_g,grid_points_direction Grid sizes of the numerical integration.
#' @return The log Bayes factor.
#' @export
log_bayes_factor_flat_r_vs_intercept_only <- function(explained_share, sample_size,
                                                      number_of_predictors,
                                                      grid_points_g = 800,
                                                      grid_points_direction = 200) {
  if (number_of_predictors == 0) return(0)
  if (sample_size < number_of_predictors + 1)
    stop("Need at least one more observation than predictors.")

  ## 1. Average over g: the one-predictor prior written as a mixture of
  ##    normal priors with variance g (g ~ Inverse-Gamma(1, n/2)), on a log-g grid.
  g_grid   <- inverse_gamma_g_grid(sample_size, grid_points_g)
  g_values <- g_grid$g_values

  ## 2. Average over directions: only the alignment matters, the cosine of the
  ##    angle between a direction and the direction the data point in
  if (number_of_predictors == 1) {
    alignment_squared     <- 1
    log_direction_weights <- 0
  } else {
    ## Angle between a direction and the data's direction. Only (0, pi/2) is
    ## needed: angles phi and pi - phi give the same alignment. The grid packs
    ## its points near angle 0, where the integrand becomes a narrow peak when
    ## the evidence is strong (packing_power = 1 would be an even grid).
    grid_position <- (seq_len(grid_points_direction) - 0.5) / grid_points_direction
    packing_power <- 3
    angles        <- (pi / 2) * grid_position^packing_power
    angle_steps   <- (pi / 2) * packing_power * grid_position^(packing_power - 1) /
                     grid_points_direction
    ## exact normalizing constant: integral of sin^(k-2) over (0, pi/2)
    log_total_weight <- 0.5 * log(pi) - log(2) +
      lgamma((number_of_predictors - 1) / 2) - lgamma(number_of_predictors / 2)
    alignment_squared     <- cos(angles)^2
    log_direction_weights <- (number_of_predictors - 2) * log(sin(angles)) +
                             log(angle_steps) - log_total_weight
  }

  shrink_factor  <- g_values / (1 + g_values)
  log_fit_matrix <- -(sample_size - 1) / 2 *
    log1p(-pmin(explained_share, 1 - 1e-12) * outer(shrink_factor, alignment_squared)) +
    matrix(log_direction_weights, grid_points_g, length(alignment_squared), byrow = TRUE)

  ## 3. The one-predictor Bayes factor for each (g, alignment), averaged
  log_sum_exp(g_grid$log_weights - 0.5 * log1p(g_values) +
              log_sum_exp_per_row(log_fit_matrix))
}

## g ~ Inverse-Gamma(1, n/2) on a log-g grid: the scale of both priors
## (docs/2_how_it_calculates.md, section 8: g/n = share x W, W ~ IG(1, 1/2)).
## Each weight is the IG density times g times the step (the step is in log g).
inverse_gamma_g_grid <- function(sample_size, grid_points_g) {
  log_g_values <- seq(log(1e-6), log(1e10), length.out = grid_points_g)
  g_values     <- exp(log_g_values)
  list(g_values    = g_values,
       log_weights = log(sample_size / 2) - log_g_values - sample_size / (2 * g_values) +
                     log(log_g_values[2] - log_g_values[1]))
}

#' Log Bayes factor of a model with k predictors against the intercept-only model, uniform correlations
#'
#' The prior with every correlation matrix of the outcome and the k predictors
#' equally likely: R-squared ~ Beta(k/2, 1) with a uniform direction (rstanarm
#' \code{stan_lm} with eta = 1). As a mixture of g-priors its scale is the same
#' for every k, g ~ Inverse-Gamma(1, n/2), so the Bayes factor is the standard
#' g-prior Bayes factor averaged over g: one integral, on a log-g grid.
#' @inheritParams log_bayes_factor_flat_r_vs_intercept_only
#' @return The log Bayes factor.
#' @export
log_bayes_factor_uniform_correlations_vs_intercept_only <- function(explained_share, sample_size,
                                                                    number_of_predictors,
                                                                    grid_points_g = 800) {
  if (number_of_predictors == 0) return(0)
  if (sample_size < number_of_predictors + 1)
    stop("Need at least one more observation than predictors.")
  g_grid   <- inverse_gamma_g_grid(sample_size, grid_points_g)
  g_values <- g_grid$g_values
  ## the g-prior Bayes factor (Liang et al., 2008) for each g, averaged
  log_g_prior_bayes_factor <- (sample_size - 1 - number_of_predictors) / 2 * log1p(g_values) -
    (sample_size - 1) / 2 * log1p(g_values * (1 - pmin(explained_share, 1 - 1e-12)))
  log_sum_exp(g_grid$log_weights + log_g_prior_bayes_factor)
}

## Either prior, by name
log_bayes_factor_vs_intercept_only <- function(explained_share, sample_size, number_of_predictors,
                                               prior) {
  switch(prior,
    flat_r = log_bayes_factor_flat_r_vs_intercept_only(explained_share, sample_size,
                                                       number_of_predictors),
    uniform_correlations = log_bayes_factor_uniform_correlations_vs_intercept_only(
      explained_share, sample_size, number_of_predictors))
}

## Pull predictors and outcome out of a flatcor_fit or lm fit
extract_predictors_and_outcome <- function(model) {
  if (inherits(model, "flatcor_fit")) {
    predictors <- model$predictors
    outcome    <- model$outcome
  } else {
    predictors <- model.matrix(model)[, -1, drop = FALSE]
    outcome    <- as.numeric(model.response(model.frame(model)))
  }
  list(predictors = predictors, outcome = outcome)
}

#' Bayes factor between any two linear models
#'
#' Nested or not, as long as both use the same outcome and data rows.
#' @param model_a,model_b Fits from \code{fit_regression()} or \code{lm()}.
#' @param prior \code{"uniform_correlations"} (the default: every correlation
#'   matrix equally likely, R-squared ~ Beta(k/2, 1)) or \code{"flat_r"}
#'   (R-squared ~ Beta(1/2, 1) for every k). See \code{docs/0_overview.md}.
#' @return The Bayes factor of model A against model B.
#' @examples
#' bayes_factor_models(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ wt, mtcars))
#' bayes_factor_models(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ wt, mtcars),
#'                     prior = "flat_r")
#' @export
bayes_factor_models <- function(model_a, model_b, prior = c("uniform_correlations", "flat_r")) {
  prior <- match.arg(prior)
  parts_a <- extract_predictors_and_outcome(model_a)
  parts_b <- extract_predictors_and_outcome(model_b)
  if (!isTRUE(all.equal(parts_a$outcome, parts_b$outcome)))
    stop("Both models need the same outcome and the same data rows.")
  sample_size <- length(parts_a$outcome)
  log_bayes_factor_of <- function(parts)
    log_bayes_factor_vs_intercept_only(
      sample_r_squared(parts$predictors, parts$outcome), sample_size, ncol(parts$predictors), prior)
  exp(log_bayes_factor_of(parts_a) - log_bayes_factor_of(parts_b))
}

#' Bayes factor for every term of a model
#'
#' With against without each term. By default every term gets the conditional
#' test (\code{log_bayes_factor_flat_r_conditional()} or
#' \code{log_bayes_factor_uniform_correlations_conditional()}): its evidence
#' depends only on its own contribution, which avoids the conditional Lindley
#' paradox. A term with one column (a continuous predictor, a two-level factor)
#' is tested by its partial correlation; a term with several columns (a factor
#' with more levels, a multi-column interaction) as a block, by its partial
#' R-squared given all the other columns.
#' @param formula A model formula.
#' @param data A data frame.
#' @param test \code{"conditional"} (default) or \code{"joint"} (compare the
#'   models with and without each term, for every term).
#' @param prior \code{"uniform_correlations"} (the default) or \code{"flat_r"}; see
#'   \code{bayes_factor_models()}.
#' @return A data frame with BF10 (evidence for an effect), BF01, and the test used.
#' @examples
#' bayes_factor_per_term(mpg ~ wt + hp + qsec, mtcars)
#' bayes_factor_per_term(mpg ~ wt + hp + qsec, mtcars, prior = "flat_r")
#' @export
bayes_factor_per_term <- function(formula, data, test = c("conditional", "joint"),
                                  prior = c("uniform_correlations", "flat_r")) {
  test  <- match.arg(test)
  prior <- match.arg(prior)
  log_bayes_factor_conditional <- switch(prior,
    flat_r               = log_bayes_factor_flat_r_conditional,
    uniform_correlations = log_bayes_factor_uniform_correlations_conditional)
  full_model  <- lm(formula, data)
  ## an aliased column (e.g. x3 = x1 + x2) gets no coefficient and would still
  ## count in k
  aliased_columns <- names(which(is.na(coef(full_model))))
  if (length(aliased_columns) > 0)
    stop("These columns are linear combinations of the other columns (aliased): ",
         paste(aliased_columns, collapse = ", "), ". Remove them from the formula.")
  ## the reduced models use the full model's rows: lm() drops rows with a
  ## missing value in any term, so without the full model's rows a term's
  ## missing values would change the rows of the model without it
  rows_of_full_model <- data[row.names(model.frame(full_model)), , drop = FALSE]
  term_labels <- attr(terms(formula), "term.labels")
  design      <- model.matrix(full_model)
  columns_of_term <- split(colnames(design)[-1], attr(design, "assign")[-1])
  outcome     <- model.response(model.frame(full_model))
  sample_size <- nrow(design)
  number_of_predictors <- ncol(design) - 1
  results <- lapply(seq_along(term_labels), function(index) {
    columns <- columns_of_term[[as.character(index)]]
    if (test == "conditional") {
      ## the term's sample partial R-squared given all the other columns (one
      ## column: its squared partial correlation; the sign does not matter)
      residuals_without <- qr.resid(qr(design[, setdiff(colnames(design), columns), drop = FALSE]),
                                    outcome)
      partial_r <- sqrt(1 - sum(residuals(full_model)^2) / sum(residuals_without^2))
      list(bf = exp(log_bayes_factor_conditional(partial_r, sample_size, number_of_predictors,
                                                 number_of_columns = length(columns))),
           used = if (length(columns) == 1) "conditional" else "conditional (block)")
    } else {
      reduced_model <- lm(update(formula, paste(". ~ . -", term_labels[index])), rows_of_full_model)
      list(bf = bayes_factor_models(full_model, reduced_model, prior = prior), used = "joint")
    }
  })
  evidence_for_effect <- vapply(results, `[[`, numeric(1), "bf")
  data.frame(term = term_labels,
             BF10 = evidence_for_effect,          # evidence for an effect
             BF01 = 1 / evidence_for_effect,      # evidence for no effect
             test = vapply(results, `[[`, character(1), "used"),
             row.names = NULL)
}
