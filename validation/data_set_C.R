## =====================================================================
##  Data set C: variable selection (validation/ADEMP.md)
##
##  Twelve candidate predictors, of which four matter (slopes 0.3, 0.2, 0.15
##  and 0.1); the other eight are null. Predictors are multivariate normal with
##  unit variances and equal correlation rho; noise SD 1, so slopes are in
##  units of the noise SD. Each dataset has its own random stream from its ID
##  (with_stream_of() in validation/data_sets.R).
## =====================================================================

true_slopes_C <- c(0.3, 0.2, 0.15, 0.1, rep(0, 8))

design_C <- function() expand.grid(n = c(50, 100, 200, 500), rho = c(0, 0.5))

cell_id_C <- function(condition) with(condition, sprintf("C-n%d-rho%.1f", n, rho))

## One dataset of set C: replication r of the cell in condition
generate_C <- function(condition, replication) {
  k <- length(true_slopes_C); n <- condition$n
  covariance <- matrix(condition$rho, k, k); diag(covariance) <- 1
  with_stream_of(paste(cell_id_C(condition), replication, sep = "-r"), {
    predictors <- matrix(rnorm(n * k), n) %*% chol(covariance)
    colnames(predictors) <- paste0("x", seq_len(k))
    list(predictors = predictors, outcome = drop(predictors %*% true_slopes_C) + rnorm(n))
  })
}
