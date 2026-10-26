#' flatcor: objective Bayes for linear regression with flat correlations
#'
#' Define what's possible, not what you expect. The default prior makes every
#' correlation matrix of the outcome and the predictors equally likely
#' ("uniform correlations"): R-squared ~ Beta(k/2, 1), and each coefficient's
#' partial correlation is flat on (-1, 1), Jeffreys' prior, whatever the number
#' of predictors. The alternative, prior = "flat_r", is flat on the multiple
#' correlation R for every k. Both divide the total over the predictors with no
#' preference: P(total) x P(this combination | total) x unit conversion.
#'
#' Main functions: \code{fit_regression()} (estimation),
#' \code{bayes_factor_models()} and \code{bayes_factor_per_term()}
#' (Bayes factors), \code{draw_coefficients_from_prior()} (the prior as
#' a recipe).
#'
#' @importFrom stats model.frame model.matrix model.response lm update terms
#'   runif quantile dt sd
#' @keywords internal
"_PACKAGE"
