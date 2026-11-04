## =====================================================================
##  The shared data sets of the validation (validation/ADEMP.md)
##
##  Every dataset has an ID: the data set, its design cell and its
##  replication. Its random numbers come from its own L'Ecuyer-CMRG stream,
##  seeded from that ID, so the same ID always gives the same data, in any
##  study, on any machine, in any order, on any number of cores. Datasets are
##  regenerated from their ID, not stored.
## =====================================================================

## Evaluate code with the random stream of one dataset, given its ID (a
## character string), and restore the caller's random state afterwards, so
## generating a dataset never changes what the rest of a study draws. The seed
## is a hash of the ID, so neighbouring IDs do not get neighbouring seeds.
with_stream_of <- function(dataset_id, code) {
  caller_kind <- RNGkind()
  caller_seed <- get0(".Random.seed", envir = globalenv())
  on.exit({
    do.call(RNGkind, as.list(caller_kind))
    if (is.null(caller_seed)) rm(".Random.seed", envir = globalenv())
    else assign(".Random.seed", caller_seed, envir = globalenv())
  })
  RNGkind("L'Ecuyer-CMRG")
  set.seed(digest::digest2int(dataset_id))
  code
}

## ---- Data set A: fixed true effects -----------------------------------------
## Predictors multivariate normal with unit variances and equal correlation rho;
## noise SD 1, so slopes are in units of the noise SD. Of the other slopes,
## "sparse" has one at 0.3 and the rest 0, "dense" has all at 0.1.
design_A <- function() {
  expand.grid(k = c(2, 5, 10, 20), n = c(50, 200, 1000), rho = c(0, 0.3, 0.7),
              others = c("sparse", "dense"), slope_x1 = c(0, 0.1, 0.2, 0.3),
              stringsAsFactors = FALSE)
}

cell_id_A <- function(condition)
  with(condition, sprintf("A-k%d-n%d-rho%.1f-%s-slope%.1f", k, n, rho, others, slope_x1))

slopes_A <- function(k, others, slope_x1) {
  other_slopes <- if (others == "sparse") c(0.3, rep(0, k - 2)) else rep(0.1, k - 1)
  c(slope_x1, other_slopes)
}

## One dataset of set A: replication r of the cell in condition
generate_A <- function(condition, replication) {
  k <- condition$k; n <- condition$n
  covariance <- matrix(condition$rho, k, k); diag(covariance) <- 1
  with_stream_of(paste(cell_id_A(condition), replication, sep = "-r"), {
    predictors <- matrix(rnorm(n * k), n) %*% chol(covariance)
    colnames(predictors) <- paste0("x", seq_len(k))
    outcome <- drop(predictors %*% slopes_A(k, condition$others, condition$slope_x1)) + rnorm(n)
    list(predictors = predictors, outcome = outcome)
  })
}
