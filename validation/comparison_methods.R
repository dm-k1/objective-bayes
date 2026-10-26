## =====================================================================
##  The methods compared in the validation studies (03 on): one set, used
##  everywhere. The flat-R method itself lives in R/bayes_factor.R and
##  R/conditional.R.
##
##  Four priors. Two follow from flat on the correlation with one predictor:
##    flat-R                 R-squared ~ Beta(1/2, 1) for every k (this package)
##    uniform correlations   R-squared ~ Beta(k/2, 1): every correlation matrix of
##                           (y, x1 ... xk) equally likely; rstanarm stan_lm's
##                           R-squared prior with eta = 1
##  Two are the standard defaults (BayesFactor / JASP; Liang et al., 2008):
##    JZS default            Cauchy prior, scale sqrt(2)/4 ("medium")
##    Zellner-Siow           Cauchy prior, scale 1 (R-squared ~ Beta(k/2, 1/2))
##
##  Testing one coefficient, two ways:
##    joint        the model with against without it (all four priors)
##    conditional  only the coefficient's partial correlation, with the prior the
##                 joint prior implies for one coefficient (flat-R and uniform
##                 correlations). For uniform correlations that implied prior is
##                 a flat partial correlation, for any k (check C4 in study 05).
##  Comparing whole models (studies 08, 09): the four priors.
##
##  Also here, for the archive and check C3: the approximate per-coefficient
##  flat partial correlation and the R2-D2 Bayes factor.
## =====================================================================

## ---- Method names and colours (the same in every study) --------------------
one_coefficient_methods <- c("JZS default", "Zellner-Siow", "flat-R joint", "flat-R conditional",
                             "uniform correlations joint", "uniform correlations conditional")
model_methods <- c("JZS default", "Zellner-Siow", "flat-R", "uniform correlations")
method_colours <- c(`JZS default` = "grey40", `Zellner-Siow` = "#7570b3",
                    `flat-R` = "#d95f02", `flat-R joint` = "#d95f02",
                    `flat-R conditional` = "#e7298a",
                    `uniform correlations` = "#a6761d", `uniform correlations joint` = "#a6761d",
                    `uniform correlations conditional` = "#1b9e77",
                    BIC = "#66a61e", `ordinary (lm)` = "black")


## ---- R-squared ~ Beta(a, 1) with a uniform direction ------------------------
## The g-prior mixture g/n = B * W, B ~ Beta(a, k/2 - a), W ~ Inverse-Gamma(1, 1/2),
## for a <= k/2 (a = 1/2: flat-R; a = k/2: uniform correlations, B = 1).
## Log Bayes factor against the intercept-only model: an average over B of a
## 1-D integral over W, both on grids (checks U1 and U2 in study 03).
log_bayes_factor_beta_r_squared <- function(r_squared, n, k, a, points_B = 100, points_w = 2000) {
  if (k == 0) return(0)
  log_w <- seq(log(1e-9), log(1e9), length.out = points_w); w <- exp(log_w)
  log_w_weight <- log(0.5) - log_w - 1 / (2 * w) + log(log_w[2] - log_w[1])   # IG(1, 1/2) x w
  B <- if (k / 2 - a <= 1e-12) 1 else qbeta((seq_len(points_B) - 0.5) / points_B, a, k / 2 - a)
  g <- n * outer(B, w)
  log_integrand <- (n - 1 - k) / 2 * log1p(g) - (n - 1) / 2 * log1p(g * (1 - r_squared)) +
                   matrix(log_w_weight, length(B), points_w, byrow = TRUE)
  per_B <- apply(log_integrand, 1, function(v) { m <- max(v); m + log(sum(exp(v - m))) })
  largest <- max(per_B); largest + log(mean(exp(per_B - largest)))
}


## ---- Whole models: log Bayes factor against the intercept-only model --------
## simple = FALSE returns the log Bayes factor (simple = TRUE overflows to Inf)
log_bayes_factor_model <- function(method, r_squared, n, k) {
  if (k == 0) return(0)
  switch(method,
    `flat-R`               = log_bayes_factor_flat_r_vs_intercept_only(r_squared, n, k),
    `uniform correlations` = log_bayes_factor_beta_r_squared(r_squared, n, k, a = k / 2),
    `JZS default`  = BayesFactor::linearReg.R2stat(n, k, r_squared, rscale = sqrt(2) / 4, simple = FALSE)$bf,
    `Zellner-Siow` = BayesFactor::linearReg.R2stat(n, k, r_squared, rscale = 1, simple = FALSE)$bf,
    BIC            = -(n * log(1 - r_squared) + k * log(n)) / 2,
    stop("Unknown method: ", method))
}


## ---- One coefficient --------------------------------------------------------
## Sample partial correlation of one coefficient, from its t-statistic
partial_correlation_of <- function(lm_fit, coefficient_row) {
  t_value <- summary(lm_fit)$coefficients[coefficient_row, "t value"]
  t_value / sqrt(t_value^2 + lm_fit$df.residual)
}

## Conditional test under uniform correlations: the implied prior for one
## coefficient is a flat partial correlation whatever k is (f^2 = W chi^2_1,
## the one-predictor prior), so this is the conditional flat-R test with the
## share fixed at 1. Uses log_peak_integral() from R/conditional.R.
log_bayes_factor_uniform_correlations_conditional <- function(partial_r, n, k, points_R = 400) {
  m <- n - k - 1
  R <- (seq_len(points_R) - 0.5) / points_R
  f <- R / sqrt(1 - R^2)
  b <- sqrt(n) * partial_r * f
  log_null <- log_peak_integral(0, m)
  log_ratio <- c(-n * f^2 / 2 + log_peak_integral( b, m) - log_null,
                 -n * f^2 / 2 + log_peak_integral(-b, m) - log_null)
  largest <- max(log_ratio)
  largest + log(mean(exp(log_ratio - largest)))
}

## Log Bayes factor for one coefficient (effect against none), every method.
## Needs the R-squared of the model with and without it, its partial
## correlation, n and the number of predictors k in the full model.
log_bayes_factor_one_coefficient <- function(method, r_squared_full, r_squared_without,
                                             partial_r, n, k) {
  joint <- function(model_method)
    log_bayes_factor_model(model_method, r_squared_full, n, k) -
    log_bayes_factor_model(model_method, r_squared_without, n, k - 1)
  switch(method,
    `JZS default`                = joint("JZS default"),
    `Zellner-Siow`               = joint("Zellner-Siow"),
    `flat-R joint`               = joint("flat-R"),
    `uniform correlations joint` = joint("uniform correlations"),
    `flat-R conditional` = log_bayes_factor_flat_r_conditional(partial_r, n, k),
    `uniform correlations conditional` = log_bayes_factor_uniform_correlations_conditional(partial_r, n, k),
    stop("Unknown method: ", method))
}

## All one-coefficient methods for one tested coefficient of a fitted model
log_bayes_factors_one_coefficient <- function(full_fit, tested_column, outcome, predictors,
                                              methods = one_coefficient_methods) {
  n <- length(outcome); k <- ncol(predictors)
  r_squared_full    <- summary(full_fit)$r.squared
  r_squared_without <- if (k == 1) 0 else
    summary(lm(outcome ~ predictors[, -tested_column, drop = FALSE]))$r.squared
  partial_r <- partial_correlation_of(full_fit, tested_column + 1)
  data.frame(method = methods,
             log_BF10 = vapply(methods, log_bayes_factor_one_coefficient, numeric(1),
                               r_squared_full = r_squared_full, r_squared_without = r_squared_without,
                               partial_r = partial_r, n = n, k = k),
             row.names = NULL)
}


## ---- Approximate per-coefficient flat partial correlation ------------------
## Flat prior on one partial correlation, likelihood evaluated on the slope with
## the estimated noise plugged in (not exact). Used by check C3 in study 05,
## which compares it with the exact version above, and by the archive.
partial_correlation_grid <- seq(-1, 1, length.out = 8003)[-c(1, 8003)]
partial_correlation_step <- partial_correlation_grid[2] - partial_correlation_grid[1]

likelihood_over_partial_correlation_grid <- function(lm_fit, predictor_index, predictors) {
  coefficient_row <- summary(lm_fit)$coefficients[predictor_index + 1, ]
  estimate        <- coefficient_row[[1]]
  standard_error  <- coefficient_row[[2]]
  degrees_of_freedom <- lm_fit$df.residual
  spread_left_after_others <- if (ncol(predictors) > 1)
    sd(resid(lm(predictors[, predictor_index] ~ predictors[, -predictor_index]))) else
    sd(predictors[, predictor_index])
  conversion_to_slope <- spread_left_after_others / summary(lm_fit)$sigma
  slopes <- partial_correlation_grid /
            (conversion_to_slope * sqrt(1 - partial_correlation_grid^2))
  list(slopes     = slopes,
       likelihood = dt((estimate - slopes) / standard_error, degrees_of_freedom) / standard_error,
       likelihood_at_zero = dt(estimate / standard_error, degrees_of_freedom) / standard_error)
}

bayes_factor_flat_partial_correlation <- function(lm_fit, predictor_index, predictors) {
  grid <- likelihood_over_partial_correlation_grid(lm_fit, predictor_index, predictors)
  evidence_effect <- sum(grid$likelihood * 0.5) * partial_correlation_step   # flat: density 1/2
  evidence_effect / grid$likelihood_at_zero
}


## ---- R2-D2 with the same R-squared prior (archive) --------------------------
## coefficient_j | sigma ~ N(0, sigma^2 * omega * share_j / var(x_j)),
## omega = R2 / (1 - R2), R2 ~ Beta(1/2, 1), shares ~ Dirichlet(1/2).
## Bayes factor by Monte Carlo over shared prior draws. Used by
## archive/03_simulations_with_r2d2.R.
draw_from_r2d2_prior <- function(number_of_predictors, number_of_draws = 1000,
                                 r2_shape_1 = 0.5, r2_shape_2 = 1, dirichlet_concentration = 0.5) {
  explained_share <- rbeta(number_of_draws, r2_shape_1, r2_shape_2)
  shares <- matrix(rgamma(number_of_draws * number_of_predictors, dirichlet_concentration),
                   number_of_draws, number_of_predictors)
  list(signal_to_noise = explained_share / (1 - explained_share),
       shares = shares / rowSums(shares))
}

log_bayes_factor_r2d2_vs_intercept_only <- function(centered_predictors, centered_outcome,
                                                    prior_draws) {
  number_of_predictors <- ncol(centered_predictors)
  sample_size          <- length(centered_outcome)
  outcome_sum_of_squares <- sum(centered_outcome^2)
  cross_products   <- crossprod(centered_predictors)
  predictor_outcome <- drop(crossprod(centered_predictors, centered_outcome))
  predictor_variances <- colSums(centered_predictors^2) / sample_size

  log_evidence_per_draw <- vapply(seq_along(prior_draws$signal_to_noise), function(m) {
    prior_sd <- sqrt(prior_draws$signal_to_noise[m] * prior_draws$shares[m, ] / predictor_variances)
    root <- chol(diag(number_of_predictors) + (prior_sd %o% prior_sd) * cross_products)
    projected <- backsolve(root, prior_sd * predictor_outcome, transpose = TRUE)
    -sum(log(diag(root))) - (sample_size - 1) / 2 * log1p(-sum(projected^2) / outcome_sum_of_squares)
  }, numeric(1))
  log_sum_exp(log_evidence_per_draw) - log(length(log_evidence_per_draw))
}

## JZS for one coefficient as a plain Bayes factor (archive)
bayes_factor_jzs_one_coefficient <- function(explained_share_full, explained_share_without,
                                             sample_size, number_of_predictors,
                                             rscale = sqrt(2) / 4) {
  log_with <- BayesFactor::linearReg.R2stat(sample_size, number_of_predictors, explained_share_full,
                                            rscale = rscale, simple = FALSE)$bf
  log_without <- if (number_of_predictors == 1) 0 else
    BayesFactor::linearReg.R2stat(sample_size, number_of_predictors - 1, explained_share_without,
                                  rscale = rscale, simple = FALSE)$bf
  exp(log_with - log_without)
}
