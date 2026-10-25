## =====================================================================
##  The flat-R prior for linear regression
##  "Define what's possible, not what you expect"
##
##  P(b) = P(total) x P(this combination | total) x unit conversion
##         flat on R   every split equally likely   noise units -> data units
##
##  PART 1. Estimation      fit_flat_r_model()   (Stan; print, loo, pp_check)
##  PART 2. Bayes factors   bayes_factor_flat_r(), bayes_factor_flat_r_per_term()
##  PART 3. Prior draws     draw_coefficients_from_flat_r_prior()
##
##  Why it works: docs/1_from_behaviour.md. How it calculates:
##  docs/2_how_it_calculates.md. Supported: linear regression with a normal
##  outcome, an intercept and fixed effects (see README.md).
##
##  Needs: MASS. For estimation also cmdstanr (or rstan), posterior, loo,
##  bayesplot. The Bayes factors need no Stan.
## =====================================================================


## =====================================================================
##  PART 1. ESTIMATION
## =====================================================================
## The Stan model (R/flat_r.stan) samples the prior in its own coordinates:
## R (flat) and a direction (uniform), with the coefficients computed from
## them. Sampling the coefficients directly misses the peak of the prior at
## zero for larger k (validation/01_prior_and_posterior.R).

flat_r_stan_file <- file.path("R", "flat_r.stan")
compiled_flat_r_models <- new.env()   # compile once per session

compile_flat_r_model <- function(backend) {
  if (is.null(compiled_flat_r_models[[backend]]))
    compiled_flat_r_models[[backend]] <- if (backend == "cmdstanr")
      cmdstanr::cmdstan_model(flat_r_stan_file) else
      rstan::stan_model(flat_r_stan_file)
  compiled_flat_r_models[[backend]]
}

## Fit a linear regression with the flat-R prior. Returns a "flat_r_fit":
## print() gives estimates and diagnostics, loo() and pp_check() work as
## usual, fit$draws holds all posterior draws (posterior::draws_df).
fit_flat_r_model <- function(formula, data, chains = 4, iter_warmup = 1000,
                             iter_sampling = 1000, adapt_delta = 0.95, seed = NULL,
                             backend = if (requireNamespace("cmdstanr", quietly = TRUE))
                                         "cmdstanr" else "rstan") {
  model_frame <- model.frame(formula, data)
  design      <- model.matrix(formula, model_frame)
  if (!"(Intercept)" %in% colnames(design))
    stop("fit_flat_r_model() needs a model with an intercept.")
  predictors <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  if (ncol(predictors) == 0)
    stop("fit_flat_r_model() needs at least one predictor.")
  outcome <- as.numeric(model.response(model_frame))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  stan_data <- list(N = nrow(predictors), K = ncol(predictors), X = predictors, y = outcome,
                    prior_only = 0L, proper_stand_ins = 0L)
  model <- compile_flat_r_model(backend)
  if (backend == "cmdstanr") {
    fit <- model$sample(data = stan_data, chains = chains, parallel_chains = chains,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = adapt_delta, seed = seed, refresh = 0,
                        show_messages = FALSE)
    all_draws <- posterior::as_draws_df(fit$draws())
    divergent <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  } else {
    fit <- rstan::sampling(model, data = stan_data, chains = chains, cores = chains,
                           iter = iter_warmup + iter_sampling, warmup = iter_warmup,
                           control = list(adapt_delta = adapt_delta), seed = seed, refresh = 0)
    all_draws <- posterior::as_draws_df(as.array(fit))
    divergent <- rstan::get_num_divergent(fit)
  }
  if (divergent > 0) warning(divergent, " divergent transitions; check the fit.")

  ## Coefficients under their own names; the intercept back on the data's scale
  ## (the model's alpha is the mean of y, at the predictors' means)
  coefficients <- unclass(posterior::as_draws_matrix(posterior::subset_draws(all_draws, variable = "b")))
  colnames(coefficients) <- colnames(predictors)
  intercept <- all_draws$alpha - drop(coefficients %*% colMeans(predictors))
  draws <- posterior::as_draws_df(cbind(
    data.frame(Intercept = intercept, check.names = FALSE),
    as.data.frame(coefficients, check.names = FALSE),
    data.frame(sigma = all_draws$sigma, R_squared = all_draws$R_squared)))

  structure(list(draws        = draws,
                 log_lik      = unclass(posterior::as_draws_matrix(posterior::subset_draws(all_draws, variable = "log_lik"))),
                 y_rep        = unclass(posterior::as_draws_matrix(posterior::subset_draws(all_draws, variable = "y_rep"))),
                 chain_ids    = all_draws$.chain,
                 outcome      = outcome,
                 predictors   = predictors,
                 formula      = formula,
                 divergent    = divergent),
            class = "flat_r_fit")
}

print.flat_r_fit <- function(x, digits = 3, ...) {
  cat("Linear regression with the flat-R prior\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Observations:", length(x$outcome), "  Predictors:", ncol(x$predictors),
      "  Divergent transitions:", x$divergent, "\n\n")
  estimates <- posterior::summarise_draws(
    x$draws, mean, sd,
    ~ quantile(.x, probs = c(0.025, 0.975)),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
  print(as.data.frame(estimates), digits = digits, row.names = FALSE)
  cat("\nR_squared is the population explained share; its prior is Beta(1/2, 1).\n")
  invisible(x)
}

## loo() and pp_check() as for any Stan model
loo.flat_r_fit <- function(x, ...) {
  relative_efficiency <- loo::relative_eff(exp(x$log_lik), chain_id = x$chain_ids)
  loo::loo(x$log_lik, r_eff = relative_efficiency, ...)
}
pp_check.flat_r_fit <- function(object, ndraws = 50, ...) {
  shown <- sample(nrow(object$y_rep), ndraws)
  bayesplot::ppc_dens_overlay(object$outcome, object$y_rep[shown, , drop = FALSE], ...)
}
loo     <- function(x, ...) UseMethod("loo")
pp_check <- function(object, ...) UseMethod("pp_check")
loo.default      <- function(x, ...) loo::loo(x, ...)
pp_check.default <- function(object, ...) bayesplot::pp_check(object, ...)


## =====================================================================
##  PART 2. BAYES FACTORS (closed form)
## =====================================================================
## In one sentence: the Bayes factor with k predictors is the one-predictor
## Bayes factor, averaged over every direction the effect could take
## (docs/2_how_it_calculates.md, section 7).
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

## Bayes factor of model A against model B (fit_flat_r_model() or lm() fits)
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

## Bayes factor for every term of a model: with vs without that term.
## Factors and interactions are tested as a block. Uses lm() only to
## build the reduced models' design matrices; nothing is estimated.
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


## =====================================================================
##  PART 3. PRIOR DRAWS
## =====================================================================
## The prior as a recipe, the same sentence read forwards. Used by
## docs/figures/make_figures.R (to show what the prior implies) and by
## validation/ (to simulate from the prior and as the reference to check against).

## One coefficient set per row, in units of the noise (sigma = 1):
##   1. P(total):  draw the overall correlation R, flat; Cohen's f^2 = R^2 / (1 - R^2)
##   2. P(combination | total): draw a direction with no preference. Normal
##      draws with covariance solve(Sigma) have no preferred direction in terms
##      of explained variance, so correlated predictors are handled
##   3. scale each direction so its f^2 is exactly the drawn amount
## (Drawing R from -1 to 1 or from 0 to 1 gives the same R^2; the sign is
## carried by the direction.)
draw_coefficients_from_flat_r_prior <- function(number_of_draws, predictor_covariance) {
  number_of_predictors <- ncol(predictor_covariance)
  overall_correlation  <- runif(number_of_draws, min = -1, max = 1)
  f_squared            <- overall_correlation^2 / (1 - overall_correlation^2)
  directions <- MASS::mvrnorm(number_of_draws, mu = rep(0, number_of_predictors),
                              Sigma = solve(predictor_covariance))
  directions <- matrix(directions, number_of_draws, number_of_predictors)   # keep a matrix when k = 1
  f_squared_of_direction <- rowSums((directions %*% predictor_covariance) * directions)
  sqrt(f_squared / f_squared_of_direction) * directions
}
