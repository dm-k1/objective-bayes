## =====================================================================
##  The prior as a recipe
## =====================================================================

#' Draw coefficient sets from the flat-R prior
#'
#' The prior as a recipe, in units of the noise (sigma = 1):
#' 1. P(total): draw the overall correlation R flat; Cohen's f^2 = R^2 / (1 - R^2).
#' 2. P(combination | total): draw a direction with no preference in explained
#'    variance (normal draws with covariance solve(Sigma)).
#' 3. Scale each direction so its f^2 is exactly the drawn amount.
#' @param number_of_draws Number of coefficient sets.
#' @param predictor_covariance Covariance matrix of the predictors.
#' @return A matrix with one coefficient set per row.
#' @examples
#' draws <- draw_coefficients_from_flat_r_prior(1000, diag(3))
#' r_squared <- rowSums(draws^2) / (1 + rowSums(draws^2))   # ~ Beta(1/2, 1)
#' @export
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
