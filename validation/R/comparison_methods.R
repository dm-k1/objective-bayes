## =====================================================================
##  The methods compared in the validation studies (03 on): one set, used
##  everywhere. The three priors themselves live in R/bayes_factor.R and
##  R/conditional.R.
##
##  The three contenders, all R-squared ~ Beta(a, 1) with a uniform direction
##  (R/priors.R; study 00 is where they come from):
##    flat partial r   a = k/2: every correlation matrix of (y, x1 ... xk)
##                     equally likely; rstanarm stan_lm's prior with eta = 1
##    flat R           a = 1/2: R flat for every k
##    flat R-squared   a = 1:   R-squared flat for every k
##  Two references, the standard defaults (BayesFactor / JASP; Liang et al., 2008):
##    JZS default      Cauchy prior, scale sqrt(2)/4 ("medium")
##    Zellner-Siow     Cauchy prior, scale 1 (R-squared ~ Beta(k/2, 1/2))
##
##  Testing one coefficient, two ways:
##    joint        the model with against without it (every prior)
##    conditional  only the coefficient's partial correlation, with the prior
##                 the joint prior implies for one coefficient (the three
##                 contenders). For flat partial r that implied prior is a flat
##                 partial correlation, for any k (check C4 in study 05).
##  Comparing whole models (studies 08, 09): the three contenders and the two
##  references.
##
##  Also here, for check C3 of study 05: the approximate per-coefficient flat
##  partial correlation.
## =====================================================================

## ---- Method names and colours (the same in every study) --------------------
contenders <- c(`flat partial r` = "flat_partial_r", `flat R` = "flat_r",
                `flat R-squared` = "flat_r_squared")
model_methods <- c(names(contenders), "JZS default", "Zellner-Siow")
one_coefficient_methods <- c(paste(names(contenders), "conditional"),
                             paste(names(contenders), "joint"), "JZS default", "Zellner-Siow")

## Colour and line type for the figures: one colour per prior, as in study 00
## (its conditional and joint tests share it: conditional solid, joint
## dashed); the references in grey, dashed too where they sit next to
## conditional tests (they test jointly). The legend lists each prior's tests together, then the references.
## Names may carry a suffix ("flat R conditional (block)", "JZS (lmBF)").
method_colour_of <- function(method)
  if (startsWith(method, "flat partial r")) "#0072B2" else
  if (startsWith(method, "flat R-squared")) "#009E73" else
  if (startsWith(method, "flat R"))         "#D55E00" else
  if (startsWith(method, "JZS"))            "grey25"  else
  if (startsWith(method, "Zellner"))        "grey60"  else
  if (startsWith(method, "BIC"))            "#CC79A7" else "black"

method_scales <- function(methods, linetype = TRUE) {
  methods <- unique(as.character(methods))
  group <- match(sapply(methods, method_colour_of),
                 c("#0072B2", "#D55E00", "#009E73", "grey25", "grey60", "#CC79A7", "black"))
  methods <- methods[order(group, !grepl("conditional", methods))]
  colours <- setNames(sapply(methods, method_colour_of), methods)
  has_conditional <- any(grepl("conditional", methods))
  dashed <- grepl("joint", methods) | (has_conditional & grepl("^(JZS|Zellner)", methods))
  rows <- if (has_conditional) 2 else 1
  legend <- guide_legend(nrow = rows, byrow = FALSE, override.aes = list(alpha = 1))
  c(list(scale_colour_manual(values = colours, breaks = methods, name = NULL),
         theme(legend.position = "bottom", legend.key.width = unit(2.5, "lines"))),
    if (linetype) list(scale_linetype_manual(values = setNames(ifelse(dashed, "22", "solid"), methods),
                                             breaks = methods, name = NULL),
                       guides(colour = legend,
                              linetype = guide_legend(nrow = rows, byrow = FALSE)))
    else list(guides(colour = legend)))
}


## ---- R-squared ~ Beta(a, 1) as a mixture of g-priors (a check) ---------------
## The g-prior mixture g/n = B * W, B ~ Beta(a, k/2 - a), W ~ Inverse-Gamma(1, 1/2),
## for a <= k/2 (a = 1/2: flat R; a = 1: flat R-squared, k >= 3; a = k/2: flat
## partial r, B = 1). Log Bayes factor against the intercept-only model: an
## average over B of a 1-D integral over W, both on grids. An independent route
## to the package's Bayes factors (checks A4 and B1 in study 02). B is placed at
## its quantiles, so the grid needs many points when a is close to k/2 - a.
log_bayes_factor_beta_r_squared <- function(r_squared, n, k, a, points_B = 2000, points_w = 2000) {
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


## ---- Correlation matrices from the LKJ distribution -------------------------
## Every valid correlation matrix has density proportional to det^(eta - 1);
## eta = 1 is uniform. C-vine sampler of Lewandowski, Kurowicka & Joe (2009):
## draw partial correlations (independent Betas on (-1, 1) whose parameters
## step down by 1/2 per level, ending at Beta(eta, eta)), then convert them
## into ordinary correlations. Used by 00, 05 (C7) and 09 (V4).
draw_lkj_correlation_matrix <- function(dimension, eta = 1) {
  partial <- matrix(0, dimension, dimension); correlation <- diag(dimension)
  beta_parameter <- eta + (dimension - 1) / 2
  for (level in seq_len(dimension - 1)) {
    beta_parameter <- beta_parameter - 0.5
    for (i in (level + 1):dimension) {
      partial[level, i] <- 2 * rbeta(1, beta_parameter, beta_parameter) - 1
      value <- partial[level, i]
      if (level > 1) for (l in (level - 1):1)
        value <- value * sqrt((1 - partial[l, i]^2) * (1 - partial[l, level]^2)) +
                 partial[l, i] * partial[l, level]
      correlation[level, i] <- correlation[i, level] <- value
    }
  }
  correlation
}


## ---- Whole models: log Bayes factor against the intercept-only model --------
## simple = FALSE returns the log Bayes factor (simple = TRUE overflows to Inf)
log_bayes_factor_model <- function(method, r_squared, n, k) {
  if (k == 0) return(0)
  if (method %in% names(contenders))
    return(log_bayes_factor_vs_intercept_only(r_squared, n, k, prior = contenders[[method]]))
  switch(method,
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

## The conditional tests are in the package (R/conditional.R).

## Log Bayes factor for one coefficient (effect against none), every method.
## Needs the R-squared of the model with and without it, its partial
## correlation, n and the number of predictors k in the full model.
log_bayes_factor_one_coefficient <- function(method, r_squared_full, r_squared_without,
                                             partial_r, n, k) {
  joint <- function(model_method)
    log_bayes_factor_model(model_method, r_squared_full, n, k) -
    log_bayes_factor_model(model_method, r_squared_without, n, k - 1)
  prior_name <- sub(" (joint|conditional)$", "", method)
  if (method %in% c("JZS default", "Zellner-Siow")) return(joint(method))
  if (!prior_name %in% names(contenders)) stop("Unknown method: ", method)
  if (endsWith(method, "joint")) return(joint(prior_name))
  log_bayes_factor_conditional(partial_r, n, k, prior = contenders[[prior_name]])
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
## which compares it with the exact conditional test.
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
