## =====================================================================
##  Estimation: fit_flat_r_model()
## =====================================================================
## The Stan model samples the prior in its own coordinates: R (flat) and a
## direction (uniform), with the coefficients computed from them. Sampling
## the coefficients directly misses the prior's peak at zero for larger k
## (validation/01_prior_and_posterior.R).

compiled_flat_r_models <- new.env()   # compile once per session

## The Stan model ships in inst/stan. cmdstanr writes its executable to a
## user cache folder, because the installed package folder may not be writable.
compile_flat_r_model <- function(backend) {
  if (is.null(compiled_flat_r_models[[backend]])) {
    stan_file <- system.file("stan", "flat_r.stan", package = "flatr", mustWork = TRUE)
    compiled_flat_r_models[[backend]] <- if (backend == "cmdstanr") {
      executable_folder <- tools::R_user_dir("flatr", which = "cache")
      dir.create(executable_folder, recursive = TRUE, showWarnings = FALSE)
      cmdstanr::cmdstan_model(stan_file, dir = executable_folder)
    } else {
      rstan::stan_model(stan_file)
    }
  }
  compiled_flat_r_models[[backend]]
}

#' Fit a linear regression with the flat-R prior
#'
#' Samples the prior in its own coordinates, the overall correlation R (flat)
#' and a direction (uniform), with the coefficients computed from them
#' (\code{inst/stan/flat_r.stan}). Intercept flat, noise level \eqn{1/\sigma}.
#'
#' @param formula A model formula with an intercept and at least one predictor.
#' @param data A data frame.
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed Sampler settings.
#' @param backend \code{"cmdstanr"} (default when installed) or \code{"rstan"}.
#' @return A \code{flat_r_fit}: \code{print()} shows estimates and diagnostics,
#'   \code{loo()} and \code{pp_check()} work as usual, and \code{$draws}
#'   holds all posterior draws as a \code{posterior::draws_df}.
#' @examples
#' \dontrun{
#' fit <- fit_flat_r_model(mpg ~ wt + hp, mtcars)
#' fit
#' loo(fit)
#' }
#' @export
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

#' @export
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

#' Leave-one-out cross-validation for a flat-R fit
#'
#' Runs on one core by default: fast enough here, and forking the R session is
#' not allowed in Positron or RStudio.
#' @param x A \code{flat_r_fit}.
#' @param cores Number of cores.
#' @param ... Passed to \code{loo::loo()}.
#' @exportS3Method loo::loo
loo.flat_r_fit <- function(x, cores = 1, ...) {
  relative_efficiency <- loo::relative_eff(exp(x$log_lik), chain_id = x$chain_ids, cores = cores)
  loo::loo(x$log_lik, r_eff = relative_efficiency, cores = cores, ...)
}
#' Posterior predictive check for a flat-R fit
#' @param object A \code{flat_r_fit}.
#' @param ndraws Number of predicted datasets to draw.
#' @param ... Passed to \code{bayesplot::ppc_dens_overlay()}.
#' @exportS3Method bayesplot::pp_check
pp_check.flat_r_fit <- function(object, ndraws = 50, ...) {
  shown <- sample(nrow(object$y_rep), ndraws)
  bayesplot::ppc_dens_overlay(object$outcome, object$y_rep[shown, , drop = FALSE], ...)
}

#' @importFrom loo loo
#' @export
loo::loo

#' @importFrom bayesplot pp_check
#' @export
bayesplot::pp_check
