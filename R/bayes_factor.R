## =====================================================================
##  Bayes factors for whole models, for the three priors
## =====================================================================
## The three priors are R^2 ~ Beta(a, 1) with a uniform direction (R/priors.R):
##   flat_partial_r   a = k/2: every correlation matrix of (y, x1 ... xk)
##                    equally likely (stan_lm with eta = 1)
##   flat_r           a = 1/2: R flat for every k
##   flat_r_squared   a = 1:   R^2 flat for every k
## Every model gets an exact Bayes factor against the intercept-only model; any
## two models then compare by dividing. The data enter only through the sample
## R-squared, n and k.
##
## One route works for every a: the model's k columns tested as one block
## against none (R/conditional.R), with the prior P(|f| < f) = (f^2/(1+f^2))^a
## on the length f of the effect in noise units. flat_r and flat_partial_r also
## have a closed form through g-priors, and so does flat_r_squared for k >= 3
## (a g-prior mixture whose hyperprior is computed once per k); these are
## faster. The block route reproduces all of them (validation/02, check A4),
## and is the one used for flat_r_squared with one predictor.

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

## flat_r, closed form: the one-predictor Bayes factor averaged over every
## direction the effect could take.
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

## flat_partial_r, closed form: R-squared ~ Beta(k/2, 1) with a uniform
## direction (rstanarm stan_lm with eta = 1). As a mixture of g-priors its
## scale is the same for every k, g ~ Inverse-Gamma(1, n/2), so the Bayes factor
## is the standard g-prior Bayes factor averaged over g: one integral, on a
## log-g grid.
log_bayes_factor_flat_partial_r_vs_intercept_only <- function(explained_share, sample_size,
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

## flat_r_squared for k >= 3, closed form: R^2 ~ Beta(1, 1) with a uniform
## direction is the g-prior mixture g/n = t = B W, B ~ Beta(1, k/2 - 1),
## W ~ Inverse-Gamma(1, 1/2) (docs/2_how_it_calculates.md, section 8). The
## density of t, p(t) = E[B exp(-B / (2t))] / (2 t^2), depends only on k: it is
## computed once per k (integrate() over B at each point of a log-t grid) and
## kept, so each Bayes factor is one sum over the grid, as for flat_partial_r.
## (k = 2: B = 1, the flat_partial_r prior; k = 1: the block route.)
flat_r_squared_hyperprior <- new.env()
log_t_grid_flat_r_squared <- function(number_of_predictors, points_t = 1600) {
  key <- as.character(number_of_predictors)
  if (is.null(flat_r_squared_hyperprior[[key]])) {
    log_t <- seq(log(1e-9), log(1e9), length.out = points_t)
    other_shape <- number_of_predictors / 2 - 1
    log_density <- vapply(exp(log_t), function(t) {
      weighted_beta <- function(b) b * stats::dbeta(b, 1, other_shape) * exp(-b / (2 * t))
      split <- min(1, 60 * t)   # the integrand sits near b = 0 when t is small
      value <- stats::integrate(weighted_beta, 0, split, rel.tol = 1e-10)$value +
        (if (split < 1) stats::integrate(weighted_beta, split, 1, rel.tol = 1e-10)$value else 0)
      log(value) - log(2) - 2 * log(t)
    }, numeric(1))
    flat_r_squared_hyperprior[[key]] <- list(
      log_t = log_t, log_weights = log_density + log_t + log(log_t[2] - log_t[1]))
  }
  flat_r_squared_hyperprior[[key]]
}
log_bayes_factor_flat_r_squared_vs_intercept_only <- function(explained_share, sample_size,
                                                              number_of_predictors) {
  if (number_of_predictors == 0) return(0)
  if (number_of_predictors <= 2)
    return(if (number_of_predictors == 2)
      log_bayes_factor_flat_partial_r_vs_intercept_only(explained_share, sample_size, 2) else
      log_bayes_factor_shape_vs_intercept_only(explained_share, sample_size, 1, shape = 1))
  if (sample_size < number_of_predictors + 1)
    stop("Need at least one more observation than predictors.")
  grid <- log_t_grid_flat_r_squared(number_of_predictors)
  g_values <- sample_size * exp(grid$log_t)
  log_g_prior_bayes_factor <- (sample_size - 1 - number_of_predictors) / 2 * log1p(g_values) -
    (sample_size - 1) / 2 * log1p(g_values * (1 - pmin(explained_share, 1 - 1e-12)))
  keep <- is.finite(grid$log_weights)
  log_sum_exp(grid$log_weights[keep] + log_g_prior_bayes_factor[keep])
}

## Any a, through the block route: the k columns against none, with
## P(|f| < f) = (f^2 / (1 + f^2))^a. At n = k + 1 every model fits the data
## exactly (sample R^2 = 1) and the data say nothing about the effect: the
## Bayes factor is 1 for any proper prior on it (predictive matching, Bayarri
## et al., 2012), which the closed forms give too. The grid needs n >= k + 2.
log_bayes_factor_shape_vs_intercept_only <- function(explained_share, sample_size,
                                                     number_of_predictors, shape) {
  if (number_of_predictors == 0) return(0)
  if (sample_size < number_of_predictors + 1)
    stop("Need at least one more observation than predictors.")
  if (sample_size == number_of_predictors + 1) return(0)
  ## The grid: fine around the peak of the likelihood ratio, as in the
  ## conditional test. With few residual degrees of freedom (n - k < 10) the
  ## ratio is broad and, when R^2 is near 1, the peak grid alone misses its
  ## tail (validation/02, check B2b): then also 40 points per decade of f from
  ## 1e-8 to 1e8.
  r_squared <- pmin(explained_share, 1 - 1e-12)
  m       <- sample_size - 2
  peak_f  <- sqrt(m / sample_size) * sqrt(r_squared / (1 - r_squared))
  width_f <- 1 / sqrt(sample_size * (1 - r_squared))
  f_edges <- sort(unique(c(
    seq(max(0, peak_f - (15 + 2 * sqrt(number_of_predictors - 1)) * width_f),
        peak_f + 15 * width_f, length.out = 601),
    0, if (sample_size - number_of_predictors < 10) 10^seq(-8, 8, length.out = 641))))
  f_middles <- (f_edges[-1] + f_edges[-length(f_edges)]) / 2
  log_likelihood_ratio <- log_likelihood_ratio_block(f_middles, r_squared, sample_size,
                                                     number_of_predictors, number_of_predictors)
  prior_probability <- diff((f_edges^2 / (1 + f_edges^2))^shape)
  keep <- prior_probability > 0
  log_sum_exp(log_likelihood_ratio[keep] + log(prior_probability[keep]))
}

#' Log Bayes factor of a model with k predictors against the intercept-only model
#'
#' The data enter only through the sample R-squared, n and k. The three priors
#' are R-squared ~ Beta(a, 1) with a uniform direction: \code{"flat_partial_r"}
#' (a = k/2, every correlation matrix of the outcome and predictors equally
#' likely), \code{"flat_r"} (a = 1/2, R flat) and \code{"flat_r_squared"}
#' (a = 1, R-squared flat).
#' @param explained_share Sample R-squared of the model.
#' @param sample_size Number of observations, n.
#' @param number_of_predictors Number of predictor columns, k.
#' @param prior \code{"flat_partial_r"} (the default), \code{"flat_r"} or
#'   \code{"flat_r_squared"}.
#' @return The log Bayes factor.
#' @examples
#' log_bayes_factor_vs_intercept_only(0.2, 50, 3, prior = "flat_partial_r")
#' log_bayes_factor_vs_intercept_only(0.2, 50, 3, prior = "flat_r_squared")
#' @export
log_bayes_factor_vs_intercept_only <- function(explained_share, sample_size, number_of_predictors,
                                               prior = "flat_partial_r") {
  switch(check_prior(prior, defaulted = missing(prior)),
    flat_r         = log_bayes_factor_flat_r_vs_intercept_only(explained_share, sample_size,
                                                               number_of_predictors),
    flat_partial_r = log_bayes_factor_flat_partial_r_vs_intercept_only(
      explained_share, sample_size, number_of_predictors),
    flat_r_squared = log_bayes_factor_flat_r_squared_vs_intercept_only(
      explained_share, sample_size, number_of_predictors))
}

## lm() and model.frame() drop every row with a missing value in a variable of
## the formula. The sample size n enters every Bayes factor, so say so.
warn_dropped_rows <- function(data, rows_used) {
  rows_dropped <- nrow(data) - rows_used
  if (rows_dropped > 0)
    warning(sprintf(paste("%d of %d rows have a missing value in a variable of the formula",
                          "and were left out; the analysis uses n = %d."),
                    rows_dropped, nrow(data), rows_used), call. = FALSE)
}

## Pull predictors and outcome out of a flatcor_fit or lm fit
extract_predictors_and_outcome <- function(model) {
  if (inherits(model, "flatcor_fit")) {
    predictors <- model$predictors
    outcome    <- model$outcome
  } else {
    if (!is.null(model$na.action))
      warning(sprintf("%d rows with a missing value were left out of %s; it uses n = %d.",
                      length(model$na.action), deparse(formula(model)), nobs(model)),
              call. = FALSE)
    predictors <- model.matrix(model)[, -1, drop = FALSE]
    outcome    <- as.numeric(model.response(model.frame(model)))
  }
  list(predictors = predictors, outcome = outcome)
}

#' Bayes factor between any two linear models
#'
#' Nested or not, as long as both use the same outcome and data rows.
#' @param model_a,model_b Fits from \code{fit_regression()} or \code{lm()}.
#' @param prior \code{"flat_partial_r"} (every correlation matrix equally
#'   likely, R-squared ~ Beta(k/2, 1)), \code{"flat_r"} (R-squared ~
#'   Beta(1/2, 1) for every k) or \code{"flat_r_squared"} (R-squared ~
#'   Beta(1, 1) for every k). The default is \code{"flat_partial_r"}, the
#'   prior the validation recommends; the comparison is in
#'   \code{docs/3_comparison.md}.
#' @return The Bayes factor of model A against model B.
#' @examples
#' bayes_factor_models(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ wt, mtcars),
#'                     prior = "flat_partial_r")
#' bayes_factor_models(lm(mpg ~ wt + hp, mtcars), lm(mpg ~ wt, mtcars),
#'                     prior = "flat_r")
#' @export
bayes_factor_models <- function(model_a, model_b, prior = "flat_partial_r") {
  prior <- check_prior(prior, defaulted = missing(prior))
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
#' test (\code{log_bayes_factor_conditional()}): its evidence
#' depends only on its own contribution, which avoids the conditional Lindley
#' paradox. A term with one column (a continuous predictor, a two-level factor)
#' is tested by its partial correlation; a term with several columns (a factor
#' with more levels, a multi-column interaction) as a block, by its partial
#' R-squared given all the other columns.
#' @param formula A model formula.
#' @param data A data frame.
#' @param test \code{"conditional"} (default) or \code{"joint"} (compare the
#'   models with and without each term, for every term).
#' @param prior \code{"flat_partial_r"} (the default), \code{"flat_r"} or
#'   \code{"flat_r_squared"}. See \code{bayes_factor_models()}.
#' @return A data frame with BF10 (evidence for an effect), BF01, and the test used.
#' @examples
#' bayes_factor_per_term(mpg ~ wt + hp + qsec, mtcars, prior = "flat_partial_r")
#' bayes_factor_per_term(mpg ~ wt + hp + qsec, mtcars, prior = "flat_r")
#' @export
bayes_factor_per_term <- function(formula, data, prior = "flat_partial_r",
                                  test = c("conditional", "joint")) {
  prior <- check_prior(prior, defaulted = missing(prior))
  test  <- match.arg(test)
  full_model  <- lm(formula, data)
  warn_dropped_rows(data, nobs(full_model))
  full_model$na.action <- NULL   # warned once here, not again per term
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
      list(bf = exp(log_bayes_factor_conditional(partial_r, sample_size, number_of_predictors, prior,
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
