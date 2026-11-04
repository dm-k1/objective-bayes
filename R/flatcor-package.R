#' flatcor: objective Bayes for linear regression with flat correlations
#'
#' Define what's possible, not what you expect. With one predictor, Jeffreys'
#' flat correlation is the only choice; with several, not everything can stay
#' flat. Three priors, built and validated side by side, all R-squared ~
#' Beta(a, 1) with a uniform direction:
#' \itemize{
#'   \item \code{"flat_partial_r"}: every correlation matrix of the outcome and
#'     the predictors equally likely; each partial correlation flat on (-1, 1)
#'     whatever the number of predictors (a = k/2).
#'   \item \code{"flat_r"}: the multiple correlation R flat for every k (a = 1/2).
#'   \item \code{"flat_r_squared"}: R-squared flat for every k (a = 1).
#' }
#' There is no default: every function asks for \code{prior =}. All three divide
#' the total over the predictors with no preference:
#' P(total) x P(this combination | total) x unit conversion.
#'
#' Main functions: \code{fit_regression()} (estimation),
#' \code{bayes_factor_models()} and \code{bayes_factor_per_term()}
#' (Bayes factors), \code{log_bayes_factor_vs_intercept_only()} and
#' \code{log_bayes_factor_conditional()} (the formulas behind them),
#' \code{draw_coefficients_from_prior()} (the prior as a recipe).
#'
#' @importFrom stats model.frame model.matrix model.response lm update terms
#'   runif quantile dt sd
#' @keywords internal
"_PACKAGE"
