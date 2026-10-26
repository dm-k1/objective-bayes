## =====================================================================
##  Bayes factors (closed form)
## =====================================================================
## In one sentence: the Bayes factor with k predictors is the one-predictor
## Bayes factor, averaged over every direction the effect could take
## (docs/2_how_it_calculates.md, section 7). Every model gets an exact Bayes
## factor against the intercept-only model; any two models then compare by
## dividing. The data enter only through the sample R-squared, n and k.

log_sum_exp <- function(log_values) {
  largest <- max(log_values)
  largest + log(sum(exp(log_values - largest)))
}
log_sum_exp_per_row <- function(log_matrix) {
  largest_per_row <- apply(log_matrix, 1, max)
  largest_per_row + log(rowSums(exp(log_matrix - largest_per_row)))
}

## Sample R-squared (always with intercept). Pure arithmetic, no fitting.
sample_r_squared <- function(predictor_columns, outcome) {
  if (is.null(predictor_columns) || NCOL(predictor_columns) == 0) return(0)
  residuals_left <- qr.resid(qr(cbind(1, predictor_columns)), outcome)
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
  ##    The weight is the IG density times g (the step is in log g).
  log_g_values  <- seq(log(1e-6), log(1e10), length.out = grid_points_g)
  g_values      <- exp(log_g_values)
  log_g_step    <- log_g_values[2] - log_g_values[1]
  log_g_weights <- log(sample_size / 2) - log_g_values - sample_size / (2 * g_values)

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
  log_sum_exp(log_g_weights - 0.5 * log1p(g_values) +
              log_sum_exp_per_row(log_fit_matrix)) + log(log_g_step)
}

## Pull predictors and outcome out of a flat_r_fit or lm fit
extract_predictors_and_outcome <- function(model) {
  if (inherits(model, "flat_r_fit")) {
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
#' @param model_a,model_b Fits from \code{fit_flat_r_model()} or \code{lm()}.
#' @return The Bayes factor of model A against model B.
#' @examples
#' bayes_factor_flat_r(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ wt, mtcars))
#' @export
bayes_factor_flat_r <- function(model_a, model_b) {
  parts_a <- extract_predictors_and_outcome(model_a)
  parts_b <- extract_predictors_and_outcome(model_b)
  if (!isTRUE(all.equal(parts_a$outcome, parts_b$outcome)))
    stop("Both models need the same outcome and the same data rows.")
  sample_size <- length(parts_a$outcome)
  log_bayes_factor_of <- function(parts)
    log_bayes_factor_flat_r_vs_intercept_only(
      sample_r_squared(parts$predictors, parts$outcome), sample_size, ncol(parts$predictors))
  exp(log_bayes_factor_of(parts_a) - log_bayes_factor_of(parts_b))
}

#' Bayes factor for every term of a model
#'
#' With against without each term; factors and interactions are tested as a
#' block. Uses \code{lm()} only to build the design matrices.
#' @param formula A model formula.
#' @param data A data frame.
#' @return A data frame with BF10 (evidence for an effect) and BF01 per term.
#' @examples
#' bayes_factor_flat_r_per_term(mpg ~ wt + hp + qsec, mtcars)
#' @export
bayes_factor_flat_r_per_term <- function(formula, data) {
  full_model <- lm(formula, data)
  term_labels <- attr(terms(formula), "term.labels")
  evidence_for_effect <- sapply(term_labels, function(term) {
    reduced_model <- lm(update(formula, paste(". ~ . -", term)), data)
    bayes_factor_flat_r(full_model, reduced_model)
  })
  data.frame(term = term_labels,
             BF10 = evidence_for_effect,          # evidence for an effect
             BF01 = 1 / evidence_for_effect,      # evidence for no effect
             row.names = NULL)
}
