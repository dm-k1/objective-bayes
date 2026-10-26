## =====================================================================
##  Support curve: which effect sizes can the data tell apart from zero?
## =====================================================================
## For one coefficient: how much better each possible partial correlation
## explains the data than zero does. Report next to any Bayes factor that
## favours "no effect": it shows which effect sizes were ruled out and
## which weren't. The peak of the curve is the most any Bayes factor could
## be; the conditional Bayes factor is exactly this curve averaged over the
## prior (R/conditional.R).
##
## Exact: it is the likelihood ratio of the conditional test, with the other
## coefficients and the intercept integrated out over flat priors and the
## noise over 1/sigma, so no plug-in noise level. A partial correlation rho
## means an effect of f = rho / sqrt(1 - rho^2) in noise units (with the
## predictor's spread left after the others measured in the sample).
##   one column:  exp(-n f^2/2) I(sqrt(n) r f) / I(0),  m = n - k - 1,
##                signed: the curve tells positive from negative effects
##   a term with q columns (a factor with q + 1 levels): its partial R^2,
##                exp(-n f^2/2) E_c[I(sqrt(n) r f c)] / I(0),  m = n - k + q - 2,
##                averaged over the effect's direction within the term
## Both depend on the data only through the sample partial correlation (or
## partial R^2), n and k. docs/2_how_it_calculates.md, section 10.

## The design columns of a term, or one design column, of an lm() fit
columns_of_name <- function(lm_fit, name) {
  design_matrix <- model.matrix(lm_fit)
  if (all(name %in% colnames(design_matrix)[-1])) return(name)
  term_labels <- attr(terms(lm_fit), "term.labels")
  if (length(name) == 1 && name %in% term_labels)
    return(colnames(design_matrix)[attr(design_matrix, "assign") == match(name, term_labels)])
  stop(sprintf("'%s' is not a term or a design-matrix column. Terms: %s. Columns: %s",
               paste(name, collapse = ", "), paste(term_labels, collapse = ", "),
               paste(colnames(design_matrix)[-1], collapse = ", ")))
}

## The sample partial R^2 of some design columns given all the others, with
## n and k: all the conditional test needs from the data
partial_fit_of_columns <- function(lm_fit, columns) {
  design_matrix <- model.matrix(lm_fit)
  if (colnames(design_matrix)[1] != "(Intercept)")
    stop("The model needs an intercept.")
  if (anyNA(coef(lm_fit)))
    stop("Some columns are linear combinations of the others (aliased). Remove them.")
  outcome <- model.response(model.frame(lm_fit))
  residuals_without <- qr.resid(qr(design_matrix[, setdiff(colnames(design_matrix), columns),
                                                 drop = FALSE]), outcome)
  list(partial_r_squared = 1 - sum(residuals(lm_fit)^2) / sum(residuals_without^2),
       sample_size = nrow(design_matrix), number_of_predictors = ncol(design_matrix) - 1)
}

#' Support for each possible partial correlation, compared with zero
#'
#' Reporting next to a Bayes factor that favours "no effect": which effect
#' sizes the data can and cannot tell apart from zero. Exact: the likelihood
#' ratio of the conditional test, with the noise level and the other
#' coefficients integrated out (no plug-in values). The conditional Bayes
#' factor is this curve averaged over the prior. For a term with several
#' columns (a factor with more than two levels) use
#' \code{support_per_partial_r_squared()}.
#' @param lm_fit A fitted \code{lm()} with an intercept.
#' @param column_name A one-column term or a design-matrix column (see
#'   \code{colnames(model.matrix(lm_fit))}).
#' @param partial_correlations The partial correlations to evaluate.
#' @return A data frame with \code{partial_correlation} and \code{support_vs_zero}.
#' @examples
#' support <- support_per_partial_correlation(lm(mpg ~ wt + hp + qsec, mtcars), "qsec")
#' support_regions(support)
#' @export
support_per_partial_correlation <- function(lm_fit, column_name,
                                            partial_correlations = seq(-0.9, 0.9, by = 0.01)) {
  column <- columns_of_name(lm_fit, column_name)
  if (length(column) != 1)
    stop(sprintf("'%s' has %d columns: use support_per_partial_r_squared() for it.",
                 column_name, length(column)))
  fit <- partial_fit_of_columns(lm_fit, column)
  n <- fit$sample_size
  m <- n - fit$number_of_predictors - 1
  partial_r <- sign(coef(lm_fit)[[column]]) * sqrt(fit$partial_r_squared)
  f_values  <- partial_correlations / sqrt(1 - partial_correlations^2)
  data.frame(partial_correlation = partial_correlations,
             support_vs_zero = exp(-n * f_values^2 / 2 +
                                   log_peak_integral(sqrt(n) * partial_r * f_values, m) -
                                   log_peak_integral(0, m)))
}

#' Support for each possible partial R-squared of a term, compared with zero
#'
#' The block version of \code{support_per_partial_correlation()}, for a term
#' with any number of columns (a factor with more than two levels, an
#' interaction): which sizes of the term's effect, as its partial R-squared
#' given the other predictors, the data can and cannot tell apart from zero.
#' Exact, and averaged over the direction of the effect within the term, as in
#' the block conditional test, whose Bayes factor is this curve averaged over
#' the prior. For a one-column term it is the signed curve averaged over the
#' two signs.
#' @param lm_fit A fitted \code{lm()} with an intercept.
#' @param term A term label (such as \code{"group"}) or design-matrix columns.
#' @param partial_r_squared The partial R-squared values to evaluate.
#' @return A data frame with \code{partial_r_squared} and \code{support_vs_zero}.
#' @examples
#' support <- support_per_partial_r_squared(lm(mpg ~ wt + factor(cyl), mtcars), "factor(cyl)")
#' support_regions(support)
#' @export
support_per_partial_r_squared <- function(lm_fit, term,
                                          partial_r_squared = seq(0, 0.8, by = 0.005)) {
  columns <- columns_of_name(lm_fit, term)
  fit <- partial_fit_of_columns(lm_fit, columns)
  f_values <- sqrt(partial_r_squared / (1 - partial_r_squared))
  data.frame(partial_r_squared = partial_r_squared,
             support_vs_zero = exp(log_likelihood_ratio_block(
               f_values, fit$partial_r_squared, fit$sample_size, fit$number_of_predictors,
               length(columns))))
}

## Split the effect-size axis into regions: better than zero, worse than
## zero, or can't be told apart (within a factor `band`). Read left to right.
#' Split the effect-size axis into regions of support
#' @param support_table Output of \code{support_per_partial_correlation()} or
#'   \code{support_per_partial_r_squared()}.
#' @param band Within this factor of zero counts as "can't be told apart".
#' @return A data frame of regions, read from left to right, on the scale of
#'   the table's first column.
#' @export
support_regions <- function(support_table, band = 3) {
  effect_sizes <- support_table[[1]]
  verdict <- ifelse(support_table$support_vs_zero > band,
                    "explains the data better than zero",
             ifelse(support_table$support_vs_zero < 1 / band,
                    "explains the data worse than zero",
                    "can't be told apart from zero"))
  region_number <- cumsum(c(TRUE, verdict[-1] != verdict[-length(verdict)]))
  do.call(rbind, lapply(split(seq_along(verdict), region_number), function(rows)
    data.frame(verdict = verdict[rows[1]],
               from    = min(effect_sizes[rows]),
               to      = max(effect_sizes[rows]))))
}
