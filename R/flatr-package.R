#' flatr: objective Bayes for linear regression with the flat-R prior
#'
#' The total strength of all predictors together, the multiple correlation R,
#' is flat from 0 to 1 for any number of predictors, and every way of dividing
#' it over the predictors is equally plausible. As a density on the coefficients:
#' P(total) x P(this combination | total) x unit conversion.
#'
#' Main functions: \code{fit_flat_r_model()} (estimation),
#' \code{bayes_factor_flat_r()} and \code{bayes_factor_flat_r_per_term()}
#' (Bayes factors), \code{draw_coefficients_from_flat_r_prior()} (the prior as
#' a recipe).
#'
#' @importFrom stats model.frame model.matrix model.response lm update terms
#'   runif quantile dt sd
#' @keywords internal
"_PACKAGE"
