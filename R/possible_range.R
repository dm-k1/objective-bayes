## =====================================================================
##  The possible-range prior for linear regression
##  "Define what's possible, not what you expect"
##
##  Prior:
##    - explained share of variance (R-squared) ~ Beta(1/2, 1), the same
##      for any number of predictors (one predictor: flat on the
##      correlation from -1 to 1, Jeffreys' default)
##    - every way of dividing the effect over the predictors equally
##      plausible, relative to how the predictors correlate
##    - flat intercept, Jeffreys prior on the noise level
##
##  PART 1. Estimation      fit_possible_range_model()            (brms, standard output)
##  PART 2. Bayes factors   bayes_factor_possible_range()         (closed form, any two models)
##                          bayes_factor_possible_range_per_term()
##  PART 3. Reporting       support_per_partial_correlation(), support_regions()
##  PART 4. Prior sampling  draw_coefficients_from_possible_range_prior(),
##                          log_bayes_factor_possible_range_by_sampling() (slow check)
##
##  Naming: bayes_factor_<method>[_<scope>] for Bayes factors (the
##  comparison methods in validation/comparison_methods.R follow the same
##  pattern); verbs (fit_, build_, draw_) for functions that do something.
##
##  Supported: linear regression with a normal outcome, an intercept and
##  fixed effects (continuous predictors, factors, interactions,
##  polynomials, spline columns). Not yet: group-level (random) effects,
##  other families, distributional models. See README.md.
## =====================================================================
## brms is only needed for estimation (fit_possible_range_model). It is
## called as brms::... instead of library(brms), so sourcing this file for
## the Bayes factor functions is fast and does not need brms or Stan. For
## summary(), plot(), bayes_R2() etc. on a fitted model, load brms
## yourself: library(brms).


## =====================================================================
##  PART 1. ESTIMATION
## =====================================================================
## The prior written as a density on the coefficients. Whiten the
## coefficients (noise units, predictors decorrelated); call the length
## of the result s. Then density ~ (1 + s^2)^(-3/2) / s^(k-1):
## the first part is R-squared ~ Beta(1/2, 1), the second makes every
## direction equally plausible.

possible_range_stan_code <- "
  real possible_range_lpdf(vector b, real sigma, matrix covariance_root) {
    int k = rows(b);
    vector[k] whitened = covariance_root' * b / sigma;
    real length_squared = dot_self(whitened);
    return -1.5 * log1p(length_squared)
           - 0.5 * (k - 1) * log(length_squared)
           - (log(2) + 0.5 * k * log(pi()) - lgamma(0.5 * k))
           + sum(log(diagonal(covariance_root))) - k * log(sigma);
  }
"

## Builds the brms prior for a formula: each model uses its own predictors'
## covariance, so the prior adapts to correlated predictors automatically.
build_possible_range_prior <- function(formula, data) {
  predictors <- model.matrix(formula, data)[, -1, drop = FALSE]
  neutral_priors <- brms::set_prior("", class = "Intercept") +   # flat intercept
                    brms::set_prior("", class = "sigma")         # Jeffreys term added below
  if (ncol(predictors) == 0)
    return(list(prior = neutral_priors + brms::set_prior("target += -log(sigma);", check = FALSE),
                stanvars = NULL))
  centered_predictors <- scale(predictors, scale = FALSE)
  covariance_root <- t(chol(crossprod(centered_predictors) / nrow(centered_predictors)))
  list(
    prior = neutral_priors + brms::set_prior(
      "target += possible_range_lpdf(b | sigma, covariance_root) - log(sigma);",
      check = FALSE),
    stanvars = brms::stanvar(scode = possible_range_stan_code, block = "functions") +
               brms::stanvar(covariance_root, name = "covariance_root"))
}

## Fit exactly like brm(): summary(), plot(), pp_check(), bayes_R2(),
## loo() all work as usual. Extra arguments go to brm(). The family may
## be given, but only gaussian with the identity link is supported.
fit_possible_range_model <- function(formula, data, family = gaussian(), ...) {
  if (any(grepl("|", deparse(formula), fixed = TRUE)))
    stop("fit_possible_range_model() supports fixed effects only (no group-level terms like (1 | g)).")
  if (is.character(family)) family <- get(family, mode = "function")
  if (is.function(family))  family <- family()
  if (!identical(family$family, "gaussian") || !identical(family$link, "identity"))
    stop("fit_possible_range_model() supports the gaussian family with the identity link only.")
  prior_parts <- build_possible_range_prior(formula, data)
  brms::brm(formula, data = data, family = gaussian(),
            prior = prior_parts$prior, stanvars = prior_parts$stanvars,
            backend = getOption("brms.backend", "cmdstanr"), ...)
}


## =====================================================================
##  PART 2. BAYES FACTORS (closed form)
## =====================================================================
## Every model gets an exact Bayes factor against the intercept-only
## model; any two models then compare by dividing. Works for nested and
## non-nested models, as long as both use the same outcome and rows.
## The data enter only through the sample R-squared, n and k.

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

## Log Bayes factor: model with k predictors vs intercept only
log_bayes_factor_possible_range_vs_intercept_only <- function(explained_share, sample_size,
                                                              number_of_predictors,
                                                              grid_points_scale = 800,
                                                              grid_points_direction = 200) {
  if (number_of_predictors == 0) return(0)
  if (sample_size < number_of_predictors + 1)
    stop("Need at least one more observation than predictors.")

  ## average over an overall scale, weighted by the prior
  log_scale_values  <- seq(log(1e-6), log(1e10), length.out = grid_points_scale)
  scale_values      <- exp(log_scale_values)
  log_scale_step    <- log_scale_values[2] - log_scale_values[1]
  log_scale_weights <- log(sample_size / 2) - log_scale_values - sample_size / (2 * scale_values)

  ## average over every direction, compared with the direction the data point in
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

  shrink_factor  <- scale_values / (1 + scale_values)
  log_fit_matrix <- -(sample_size - 1) / 2 *
    log1p(-pmin(explained_share, 1 - 1e-12) * outer(shrink_factor, alignment_squared)) +
    matrix(log_direction_weights, grid_points_scale, length(alignment_squared), byrow = TRUE)

  log_sum_exp(log_scale_weights - 0.5 * log1p(scale_values) +
              log_sum_exp_per_row(log_fit_matrix)) + log(log_scale_step)
}

## Pull predictors and outcome out of a brms or lm fit
extract_predictors_and_outcome <- function(model) {
  if (inherits(model, "brmsfit")) {
    stan_data  <- brms::standata(model)
    predictors <- stan_data$X[, colnames(stan_data$X) != "Intercept", drop = FALSE]
    outcome    <- as.numeric(stan_data$Y)
  } else {
    predictors <- model.matrix(model)[, -1, drop = FALSE]
    outcome    <- as.numeric(model.response(model.frame(model)))
  }
  list(predictors = predictors, outcome = outcome)
}

## Bayes factor of model A against model B (fit_possible_range_model() or lm() fits)
bayes_factor_possible_range <- function(model_a, model_b) {
  parts_a <- extract_predictors_and_outcome(model_a)
  parts_b <- extract_predictors_and_outcome(model_b)
  if (!isTRUE(all.equal(parts_a$outcome, parts_b$outcome)))
    stop("Both models need the same outcome and the same data rows.")
  sample_size <- length(parts_a$outcome)
  log_bayes_factor_of <- function(parts)
    log_bayes_factor_possible_range_vs_intercept_only(
      sample_r_squared(parts$predictors, parts$outcome), sample_size, ncol(parts$predictors))
  exp(log_bayes_factor_of(parts_a) - log_bayes_factor_of(parts_b))
}

## Bayes factor for every term of a model: with vs without that term.
## Factors and interactions are tested as a block. Uses lm() only to
## build the reduced models' design matrices; nothing is estimated.
bayes_factor_possible_range_per_term <- function(formula, data) {
  full_model <- lm(formula, data)
  term_labels <- attr(terms(formula), "term.labels")
  evidence_for_effect <- sapply(term_labels, function(term) {
    reduced_model <- lm(update(formula, paste(". ~ . -", term)), data)
    bayes_factor_possible_range(full_model, reduced_model)
  })
  data.frame(term = term_labels,
             BF10 = evidence_for_effect,          # evidence for an effect
             BF01 = 1 / evidence_for_effect,      # evidence for no effect
             row.names = NULL)
}


## =====================================================================
##  PART 3. REPORTING
## =====================================================================
## For one coefficient: how much better each possible partial correlation
## explains the data than zero does. Report next to any Bayes factor that
## favours "no effect": it shows which effect sizes were ruled out and
## which weren't. The peak of the curve is the most any Bayes factor could
## be; the Bayes factor is this curve averaged over everything possible.
## Treats the noise level and predictor spread as known, so boundaries are
## approximate (about +/- 0.01-0.02 at n = 60).
##
## column_name is a column of the design matrix, not a term: for a factor
## that is one dummy column such as "groupb" (see colnames(model.matrix(lm_fit))).
support_per_partial_correlation <- function(lm_fit, column_name,
                                            partial_correlations = seq(-0.9, 0.9, by = 0.01)) {
  design_matrix <- model.matrix(lm_fit)
  if (!column_name %in% colnames(design_matrix))
    stop(sprintf("'%s' is not a design-matrix column. Available: %s", column_name,
                 paste(colnames(design_matrix)[-1], collapse = ", ")))
  coefficient_row <- summary(lm_fit)$coefficients[column_name, ]
  other_columns   <- design_matrix[, colnames(design_matrix) != column_name, drop = FALSE]
  spread_left_after_others <- sd(qr.resid(qr(other_columns), design_matrix[, column_name]))
  conversion_to_slope <- spread_left_after_others / summary(lm_fit)$sigma
  slopes <- partial_correlations / (conversion_to_slope * sqrt(1 - partial_correlations^2))
  data.frame(partial_correlation = partial_correlations,
             support_vs_zero = dt((coefficient_row[1] - slopes) / coefficient_row[2],
                                  lm_fit$df.residual) /
                               dt(coefficient_row[1] / coefficient_row[2], lm_fit$df.residual))
}

## Split the effect-size axis into regions: better than zero, worse than
## zero, or can't be told apart (within a factor `band`). Read left to right.
support_regions <- function(support_table, band = 3) {
  verdict <- ifelse(support_table$support_vs_zero > band,
                    "explains the data better than zero",
             ifelse(support_table$support_vs_zero < 1 / band,
                    "explains the data worse than zero",
                    "can't be told apart from zero"))
  region_number <- cumsum(c(TRUE, verdict[-1] != verdict[-length(verdict)]))
  do.call(rbind, lapply(split(seq_along(verdict), region_number), function(rows)
    data.frame(verdict = verdict[rows[1]],
               from    = min(support_table$partial_correlation[rows]),
               to      = max(support_table$partial_correlation[rows]))))
}


## =====================================================================
##  PART 4. PRIOR SAMPLING
## =====================================================================
## Used by the tutorial (to show what the prior implies) and by
## validation/01_checks.R (as an independent check of the closed form).

## One coefficient set per row, in units of the noise:
##   1. draw the explained share (R-squared) from Beta(1/2, 1)
##   2. turn it into an overall effect strength
##   3. draw a random direction: how the effect is divided over predictors
##   4. reshape that direction to the predictors' covariance, so that
##      correlated predictors are handled automatically
draw_coefficients_from_possible_range_prior <- function(number_of_draws, predictor_covariance) {
  number_of_predictors <- ncol(predictor_covariance)
  explained_share      <- rbeta(number_of_draws, 0.5, 1)
  overall_strength     <- sqrt(explained_share / (1 - explained_share))
  random_numbers <- matrix(rnorm(number_of_draws * number_of_predictors),
                           number_of_draws, number_of_predictors)
  directions     <- random_numbers / sqrt(rowSums(random_numbers^2))   # all directions equally likely
  covariance_root <- chol(predictor_covariance)
  t(backsolve(covariance_root, t(directions * overall_strength)))
}

## The same log Bayes factor as log_bayes_factor_possible_range_vs_intercept_only(),
## without the formula: draw coefficient sets from the prior, see how well
## each predicts the outcome (the noise level averaged on a grid), and
## average. Slow and noisy (Monte Carlo); agrees to about 0.05 on the log scale.
log_bayes_factor_possible_range_by_sampling <- function(predictors, outcome,
                                                        number_of_draws = 40000,
                                                        draws_per_batch = 2000) {
  sample_size            <- length(outcome)
  centered_predictors    <- scale(predictors, scale = FALSE)
  centered_outcome       <- outcome - mean(outcome)
  outcome_sum_of_squares <- sum(centered_outcome^2)
  predictor_covariance   <- crossprod(centered_predictors) / sample_size

  ## the noise level is unknown: average over it on a grid (written as
  ## 1 / noise, which makes the arithmetic simpler)
  inverse_noise_values <- seq(1e-4, 6 * sqrt(sample_size / outcome_sum_of_squares), length.out = 800)
  inverse_noise_step   <- inverse_noise_values[2] - inverse_noise_values[1]
  log_inverse_noise    <- log(inverse_noise_values)

  ## how well "no effect" predicts the data
  log_evidence_null <- log_sum_exp((sample_size - 2) * log_inverse_noise -
                                   inverse_noise_values^2 * outcome_sum_of_squares / 2) +
                       log(inverse_noise_step)

  ## how well each drawn coefficient set predicts the data
  log_evidence_per_draw <- numeric(0)
  for (batch in seq_len(number_of_draws / draws_per_batch)) {
    coefficients     <- draw_coefficients_from_possible_range_prior(draws_per_batch, predictor_covariance)
    predicted_values <- centered_predictors %*% t(coefficients)
    overlap_with_outcome     <- drop(crossprod(centered_outcome, predicted_values))
    predicted_sum_of_squares <- colSums(predicted_values^2)
    log_integrand <- (sample_size - 2) *
        matrix(log_inverse_noise, draws_per_batch, length(inverse_noise_values), byrow = TRUE) -
      (outer(rep(outcome_sum_of_squares, draws_per_batch), inverse_noise_values^2) -
       2 * outer(overlap_with_outcome, inverse_noise_values) + predicted_sum_of_squares) / 2
    log_evidence_per_draw <- c(log_evidence_per_draw,
                               log_sum_exp_per_row(log_integrand) + log(inverse_noise_step))
  }
  log_sum_exp(log_evidence_per_draw) - log(number_of_draws) - log_evidence_null
}
