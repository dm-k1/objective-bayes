## =====================================================================
##  Data set D: the number of predictors k grows with the sample size n
##  (docs/ADEMP.md; used by study 5)
##
##  Two growth rates, k = n / 4 and k = sqrt(n), at n = 50 ... 2,000.
##    Under the null (no predictor matters) only R-squared is needed, and it
##    is exactly Beta(k / 2, (n - k - 1) / 2), so it is drawn directly.
##    With a small true model (3 slopes of 0.3, noise SD 1, the other k - 3
##    predictors null) full datasets are drawn: standard normal predictors.
##  Every draw has its own random stream from its ID (R/data/data_sets.R).
## =====================================================================

design_D <- function() {
  expand.grid(n = c(50, 100, 200, 500, 1000, 2000), growth = c("k = n / 4", "k = sqrt(n)"),
              stringsAsFactors = FALSE)
}

k_of_D <- function(n, growth) if (growth == "k = n / 4") floor(n / 4) else floor(sqrt(n))

cell_id_D <- function(condition)
  sprintf("D-n%d-%s", condition$n, if (condition$growth == "k = n / 4") "quarter" else "root")

## R-squared of the k-predictor model when the null is true (replication r)
generate_D_null <- function(condition, replication) {
  n <- condition$n; k <- k_of_D(n, condition$growth)
  with_stream_of(paste0(cell_id_D(condition), "-null-r", replication),
                 rbeta(1, k / 2, (n - k - 1) / 2))
}

## One dataset with the small true model (replication r)
generate_D_true <- function(condition, replication) {
  n <- condition$n; k <- k_of_D(n, condition$growth)
  with_stream_of(paste0(cell_id_D(condition), "-true-r", replication), {
    predictors <- matrix(rnorm(n * k), n)
    outcome <- drop(predictors[, 1:3] %*% c(0.3, 0.3, 0.3)) + rnorm(n)
    list(predictors = predictors, outcome = outcome)
  })
}
