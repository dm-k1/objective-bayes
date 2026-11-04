## =====================================================================
##  Data set B: data drawn from each prior (docs/ADEMP.md)
##
##  Used by study 3 (calibration of the posterior), study 4 (calibration of
##  the Bayes factor) and study 6 (coherence with the joint model). Needs
##  R/data/data_sets.R (with_stream_of()) and flatcor.
##
##  One dataset: predictors with equal correlation .3 (centred, then held
##  fixed: the model conditions on them); the intercept from N(0, 10) and the
##  noise SD from lognormal(0, 1), the proper stand-ins of the Stan model;
##  the coefficients from the generating prior, set up in the predictors'
##  sample covariance (/ n), as the Stan model does. The generating prior is
##  one of the three ("flat_partial_r", "flat_r", "flat_r_squared") or "JZS"
##  (BayesFactor's Cauchy prior, scale 0.354).
##
##  x1 can be set instead of drawn:
##    x1 = "prior"   drawn with the others (study 3)
##    x1 = "null"    0 (half of study 4's datasets)
##    x1 = <number>  a target t-value: x1's partial f is t / sqrt(n) (study 6)
## =====================================================================

jzs_scale <- sqrt(2) / 4

sample_covariance <- function(predictors) crossprod(predictors) / nrow(predictors)
sd_given_others   <- function(covariance, j = 1) sqrt(1 / solve(covariance)[j, j])

## Coefficients from the generating prior, in units of the noise SD
draw_coefficients_B <- function(generating_prior, covariance) {
  if (generating_prior != "JZS")
    return(drop(draw_coefficients_from_prior(1, covariance, prior = generating_prior)))
  ## JZS: b | g ~ N(0, g Sigma^-1), g ~ Inverse-Gamma(1/2, s^2 / 2), so each
  ## coefficient's partial f is Cauchy(0, s)
  g <- 1 / rgamma(1, shape = 1 / 2, rate = jzs_scale^2 / 2)
  MASS::mvrnorm(1, rep(0, ncol(covariance)), g * solve(covariance))
}

dataset_id_B <- function(generating_prior, k, n, x1, purpose, replication)
  sprintf("B-%s-%s-k%d-n%d-x1_%s-r%d", purpose, generating_prior, k, n, x1, replication)

## One dataset of set B. purpose names the study's use ("sbc", "calibration",
## "coherence") and is part of the ID, so the studies do not share datasets
## by accident; the same arguments always give the same dataset.
generate_B <- function(generating_prior, k, n, x1 = "prior", purpose, replication,
                       predictor_correlation = 0.3) {
  with_stream_of(dataset_id_B(generating_prior, k, n, x1, purpose, replication), {
    covariance_population <- matrix(predictor_correlation, k, k); diag(covariance_population) <- 1
    predictors <- matrix(rnorm(n * k), n) %*% chol(covariance_population)
    predictors <- scale(predictors, center = TRUE, scale = FALSE)
    attributes(predictors)[c("scaled:center")] <- NULL
    colnames(predictors) <- paste0("x", seq_len(k))
    covariance <- sample_covariance(predictors)
    intercept  <- rnorm(1, 0, 10)
    noise_sd   <- rlnorm(1, 0, 1)
    coefficients <- draw_coefficients_B(generating_prior, covariance)   # noise units
    if (identical(x1, "null")) coefficients[1] <- 0
    if (is.numeric(x1)) coefficients[1] <- x1 / sqrt(n) / sd_given_others(covariance)
    coefficients <- noise_sd * coefficients                            # data units
    outcome <- intercept + drop(predictors %*% coefficients) + rnorm(n, 0, noise_sd)
    f_squared <- sum((coefficients %*% covariance) * coefficients) / noise_sd^2
    partial_f <- coefficients[1] * sd_given_others(covariance) / noise_sd
    list(predictors = predictors, outcome = outcome, intercept = intercept, noise_sd = noise_sd,
         coefficients = coefficients, r_squared = f_squared / (1 + f_squared),
         partial_correlation_x1 = partial_f / sqrt(1 + partial_f^2))
  })
}

## x1's sample partial correlation, from its t-value
sample_partial_r <- function(predictors, outcome) {
  fit <- lm.fit(cbind(1, predictors), outcome)
  df  <- length(outcome) - ncol(predictors) - 1
  sigma_hat <- sqrt(sum(fit$residuals^2) / df)
  r_inverse <- chol2inv(qr.R(fit$qr))
  t_value   <- fit$coefficients[2] / (sigma_hat * sqrt(r_inverse[2, 2]))
  unname(t_value / sqrt(t_value^2 + df))
}
