## =====================================================================
##  The three priors: one family, R-squared ~ Beta(a, 1), uniform direction
## =====================================================================
## With one predictor, Jeffreys' flat correlation is the only choice. With k
## predictors, not everything can stay flat (validation/00): keeping each
## partial correlation flat makes R-squared grow with k, keeping the total flat
## makes each partial correlation shrink with k. The three options compared:
##
##   flat_partial_r   every valid correlation matrix of (y, x1 ... xk) equally
##                    likely; each partial correlation flat on (-1, 1) for any k
##                    R-squared ~ Beta(k/2, 1)
##   flat_r           the multiple correlation R flat on (0, 1) for any k
##                    R-squared ~ Beta(1/2, 1)
##   flat_r_squared   R-squared flat on (0, 1) for any k
##                    R-squared ~ Beta(1, 1)
##
## All three divide the total over the predictors with no preference (a
## uniform direction), so they differ only in the shape a of R-squared ~
## Beta(a, 1). On R that is the density 2a R^(2a - 1); with one predictor,
## flat_partial_r and flat_r are both flat on r.
## There is no default: every function asks for the prior by name.

prior_names <- c("flat_partial_r", "flat_r", "flat_r_squared")

## Stops unless prior is one of the three names.
check_prior <- function(prior) {
  if (missing(prior) || is.null(prior))
    stop("Choose a prior: prior = \"", paste(prior_names, collapse = "\", \""), "\".",
         call. = FALSE)
  if (length(prior) != 1 || !prior %in% prior_names)
    stop("Unknown prior \"", paste(prior, collapse = ", "), "\". Choose one of \"",
         paste(prior_names, collapse = "\", \""), "\".", call. = FALSE)
  prior
}

## The shape a of R-squared ~ Beta(a, 1) for k predictor columns.
prior_shape <- function(prior, number_of_predictors) {
  switch(check_prior(prior),
    flat_partial_r = number_of_predictors / 2,
    flat_r         = 1 / 2,
    flat_r_squared = 1)
}
