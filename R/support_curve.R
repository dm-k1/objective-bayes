## Reporting next to a Bayes factor that favours "no effect". Independent of
## the flat-R prior: uses only the fitted lm() and its t-distribution.
## =====================================================================
##  Support curve: which effect sizes can the data tell apart from zero?
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
