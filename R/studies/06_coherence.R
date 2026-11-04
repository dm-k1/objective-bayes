## =====================================================================
##  Study 6: is the test of one coefficient coherent with one model?
##  (docs/appendix/C_study_design.md, study 6; replaces parts 2 and 3 of the old study 10)
##
##  Aims. For each prior: does the conditional test of x1 equal the Bayes
##  factor of the full model against its restriction by conditioning on x1's
##  partial correlation being zero (the Savage-Dickey ratio), and is that
##  restriction a proper model?
##  Data. Part 1: set B (R/data/data_set_B.R), data from flat partial r or
##  flat R-squared with x1's t-value set to 0, 1, 1.5, 2 or 2.5 (cycling over
##  the replications); k = 3 and 8, n = 50 and 200. Part 2: prior draws only.
##  Estimands. Part 1: the mean log difference between the conditional test and
##  Savage-Dickey, per prior and k. Part 2: the density of x1's partial
##  correlation near zero, and the other predictors' R-squared given it.
##  Methods. flatcor's conditional test; Savage-Dickey from the Stan posterior
##  of fit_regression() (logspline density at 0, a kernel estimate as a
##  cross-check). Flat R is covered by part 2 only: its prior density of one
##  partial correlation at 0 is infinite, so the ratio cannot be formed.
##  Performance. The mean difference and the mean absolute difference, each
##  with its MCSE (SD / sqrt(estimable datasets)); a dataset is estimable when
##  the Monte Carlo SE of its Savage-Dickey estimate, from its four chains, is
##  below 0.1.
##  n_sim. Part 1: 50 datasets per prior, k and n (precision.R,
##  coherence_mean: 36 estimable needed for MCSE 0.05 at SD 0.30). Part 2:
##  400,000 prior draws per prior and k.
##
##  Checks:
##    6.1 (old P3) Savage-Dickey by logspline and by kernel agree within 0.25
##        on every estimable dataset.
##    6.2 (old P4, flat partial r) The mean difference is within 2 MCSE of zero
##        at k = 3 and at k = 8.
##    6.3 (old P4, flat R-squared) A predicted difference, not a test of
##        coherence: under flat R-squared the partial correlations are tied
##        together through R-squared, so the conditional test is not the joint
##        model's Bayes factor. At k = 8 its mean absolute difference from
##        Savage-Dickey exceeds that of flat partial r (coherent by
##        construction, so its difference is Monte Carlo error alone) by more
##        than 2 MCSE; k = 3 is reported. (The old P4 applied the coherence
##        test to flat R-squared as well, and failed: that failure is this
##        prediction.)
##    6.4 (old P5) On the datasets that are not estimable, the conditional test
##        lies within 3 Monte Carlo SE of Savage-Dickey.
##    6.5 (old P6) Under flat partial r and flat R-squared, x1's partial
##        correlation = 0 is a proper restriction: the density near 0 settles
##        (within 10% from eps = .02 to .005) and the others keep their R^2.
##    6.6 (old P7) Under flat R, the density keeps growing as eps shrinks and
##        the others' median R^2 falls below .01 at eps = .005: the restriction
##        collapses to the intercept-only model.
##    6.7 (new) Under flat partial r the restriction gives the others
##        R^2 ~ Beta((k - 1)/2, 3/2) (stan_lm's prior with eta = 3/2 on k - 1
##        predictors): KS distance below its .01 critical value.
##
##  Run from the project folder:  Rscript R/studies/06_coherence.R
##  Smoke test: VALIDATION_SMOKE=1 Rscript R/studies/06_coherence.R
##  Raw results are cached in cache/study_06/ (delete them to rerun); the
##  tables and checks are saved to results/study_06.rds.
## =====================================================================
library(SimDesign)
library(dplyr)
library(tidyr)
library(ggplot2)
source("R/shared/settings.R")             # flatcor, cores, record_check(), save_figure()
source("R/shared/comparison_methods.R")   # contenders
source("R/shared/precision.R")            # n_sim_of()
source("R/data/data_sets.R")            # with_stream_of()
source("R/data/data_set_B.R")           # generate_B(), sample_partial_r(), sd_given_others()

study_folder <- file.path(cache_folder, if (smoke_mode) "study_06_smoke" else "study_06")
dir.create(study_folder, recursive = TRUE, showWarnings = FALSE)


## =====================================================================
##  Part 1: the conditional test against Savage-Dickey
## =====================================================================
design <- createDesign(prior = c("flat_partial_r", "flat_r_squared"), k = c(3, 8), n = c(50, 200))
if (smoke_mode) design <- subset(design, k == 3 & n == 50)
n_sim <- n_sim_of("coherence_mean")
target_t_values <- c(0, 1, 1.5, 2, 2.5)
stan_iterations <- if (smoke_mode) c(warmup = 300, sampling = 600) else c(warmup = 1000, sampling = 5000)
workers <- max(1, number_of_cores %/% 4)    # fit_regression() runs its 4 chains in parallel

## the prior density of x1's partial correlation at 0 (both signs): 1/2 for
## flat partial r (flat on (-1, 1)); otherwise half the density of |f| at 0,
## the slope of the share mixture's CDF at a tiny f
log_prior_density_at_0 <- function(prior, k) {
  if (prior == "flat_partial_r") return(log(1 / 2))
  tiny <- 1e-6
  log(prior_cdf_of_f(tiny, k, points_share = 4000, shape = prior_shape(prior, k)) / tiny / 2)
}

simulate_data <- function(condition, fixed_objects) {
  r <- condition$REPLICATION
  target_t <- target_t_values[(r - 1) %% length(target_t_values) + 1]
  c(generate_B(condition$prior, condition$k, condition$n, x1 = target_t, purpose = "coherence",
               replication = r), target_t = target_t)
}

analyse_dataset <- function(condition, dat, fixed_objects) {
  k <- condition$k; n <- condition$n
  dataset_id <- dataset_id_B(condition$prior, k, n, dat$target_t, "coherence", condition$REPLICATION)
  data <- data.frame(dat$predictors, y = dat$outcome)
  formula <- reformulate(colnames(dat$predictors), response = "y")
  partial_r <- sample_partial_r(dat$predictors, dat$outcome)
  conditional <- log_bayes_factor_conditional(partial_r, n, k, prior = condition$prior)
  fit <- suppressWarnings(fit_regression(formula, data, prior = condition$prior, chains = 4,
                                         iter_warmup = stan_iterations[["warmup"]],
                                         iter_sampling = stan_iterations[["sampling"]],
                                         seed = abs(digest::digest2int(dataset_id)) %% .Machine$integer.max))
  draws <- as.data.frame(fit$draws)
  partial_f   <- draws$x1 * sd_given_others(sample_covariance(dat$predictors)) / draws$sigma
  partial_rho <- partial_f / sqrt(1 + partial_f^2)
  log_density_at_0 <- function(values)
    log(logspline::dlogspline(0, logspline::logspline(values, lbound = -1, ubound = 1)))
  by_chain <- vapply(split(partial_rho, fit$chain_ids), log_density_at_0, numeric(1))
  kernel   <- density(partial_rho, bw = "SJ", n = 2^12, from = -1, to = 1)
  prior_at_0 <- log_prior_density_at_0(condition$prior, k)
  c(target_t = dat$target_t, conditional = conditional,
    savage_dickey        = prior_at_0 - log_density_at_0(partial_rho),
    savage_dickey_se     = sd(by_chain) / sqrt(length(by_chain)),
    savage_dickey_kernel = prior_at_0 - log(approx(kernel$x, kernel$y, 0)$y),
    divergent = fit$divergent)
}

summarise_cell <- function(condition, results, fixed_objects)
  c(mean_difference = mean(results[, "conditional"] - results[, "savage_dickey"]))

results_file <- file.path(study_folder, "simulation.rds")
if (file.exists(results_file)) {
  simulation <- readRDS(results_file)
  message("Using cached results: ", results_file)
} else {
  simulation <- runSimulation(design, replications = n_sim, generate = simulate_data,
                              analyse = analyse_dataset, summarise = summarise_cell,
                              parallel = workers > 1, ncores = workers,
                              packages = c("flatcor", "logspline"), store_results = TRUE,
                              control = list(include_replication_index = TRUE), progress = FALSE,
                              save_details = list(out_rootdir = study_folder,
                                                  tmpfilename = "resume_simulation.qs2"))
  saveRDS(simulation, results_file)
}

per_dataset <- SimExtract(simulation, what = "results") |>
  mutate(difference = conditional - savage_dickey, estimable = savage_dickey_se < 0.1)

coherence <- per_dataset |>
  filter(estimable) |>
  group_by(prior, k) |>
  summarise(estimable = n(), mean_difference = mean(difference),
            mcse = sd(difference) / sqrt(n()),
            mean_absolute_difference = mean(abs(difference)),
            mean_absolute_difference_mcse = sd(abs(difference)) / sqrt(n()),
            largest_difference = max(abs(difference)),
            .groups = "drop")
cat("\nConditional test minus Savage-Dickey (log scale), estimable datasets, per prior and k:\n")
print(as.data.frame(coherence), digits = 3, row.names = FALSE)

estimable <- filter(per_dataset, estimable)
record_check("6.1 Savage-Dickey by logspline and by kernel agree (estimable datasets)",
             all(abs(estimable$savage_dickey - estimable$savage_dickey_kernel) < 0.25),
             sprintf("%d estimable of %d; largest difference %.3f", nrow(estimable), nrow(per_dataset),
                     max(abs(estimable$savage_dickey - estimable$savage_dickey_kernel))))
coherent_at <- function(prior_name) {
  rows <- filter(coherence, prior == prior_name)
  list(passed = nrow(rows) > 0 && all(abs(rows$mean_difference) <= 2 * rows$mcse),
       details = paste(sprintf("k = %d: %.3f (MCSE %.3f, %d estimable)", rows$k, rows$mean_difference,
                               rows$mcse, rows$estimable), collapse = "; "))
}
result <- coherent_at("flat_partial_r")
record_check("6.2 flat_partial_r: conditional test within 2 MCSE of Savage-Dickey",
             result$passed, result$details)
## 6.3: flat R-squared is predicted to differ; flat partial r (coherent) gives
## the size of a difference that is Monte Carlo error alone
excess_over_coherent <- function(k_value) {
  rows <- filter(coherence, k == k_value)
  r2 <- filter(rows, prior == "flat_r_squared"); pr <- filter(rows, prior == "flat_partial_r")
  if (nrow(r2) == 0 || nrow(pr) == 0) return(c(excess = NA, mcse = NA))
  c(excess = r2$mean_absolute_difference - pr$mean_absolute_difference,
    mcse = sqrt(r2$mean_absolute_difference_mcse^2 + pr$mean_absolute_difference_mcse^2))
}
excess <- lapply(c(3, 8), excess_over_coherent)
record_check("6.3 flat_r_squared: conditional test differs from Savage-Dickey at k = 8 (predicted)",
             isTRUE(excess[[2]][["excess"]] > 2 * excess[[2]][["mcse"]]),
             sprintf("mean |difference| above flat partial r: k = 8: %.3f (MCSE %.3f); k = 3: %.3f (MCSE %.3f)",
                     excess[[2]][["excess"]], excess[[2]][["mcse"]], excess[[1]][["excess"]], excess[[1]][["mcse"]]))
not_estimable <- filter(per_dataset, !estimable)
record_check("6.4 not estimable: conditional within 3 Monte Carlo SE of Savage-Dickey",
             all(abs(not_estimable$difference) < 3 * not_estimable$savage_dickey_se),
             sprintf("%d datasets", nrow(not_estimable)))


## =====================================================================
##  Part 2: what "x1's partial correlation = 0" means under each joint prior
## =====================================================================
prior_draws <- if (smoke_mode) 40000 else 400000
widths <- c(0.1, 0.02, 0.005)
near_zero <- bind_rows(lapply(c(3, 8), function(k) bind_rows(lapply(unname(contenders), function(prior)
  with_stream_of(sprintf("prior-draws-%s-k%d", prior, k), {
    b <- draw_coefficients_from_prior(prior_draws, diag(k), prior = prior)
    rho1 <- b[, 1] / sqrt(1 + b[, 1]^2)    # x1's partial correlation (identity covariance, sigma = 1)
    others <- rowSums(b[, -1, drop = FALSE]^2)
    others_r_squared <- others / (1 + others)
    bind_rows(lapply(widths, function(eps) {
      keep <- abs(rho1) < eps
      tibble(prior = prior, k = k, eps = eps, kept = sum(keep), density_near_0 = mean(keep) / (2 * eps),
             median_r_squared_others = median(others_r_squared[keep]),
             ks_beta_3_2 = if (prior == "flat_partial_r")
               unname(suppressWarnings(ks.test(others_r_squared[keep], "pbeta", (k - 1) / 2, 1.5)$statistic))
             else NA_real_)
    }))
  })))))
cat("\nThe joint prior conditioned on |rho1| < eps: density near 0 and the others' R^2:\n")
print(as.data.frame(near_zero), digits = 3, row.names = FALSE)

for (k_value in c(3, 8)) {
  for (prior_name in c("flat_partial_r", "flat_r_squared")) {
    rows <- filter(near_zero, prior == prior_name, k == k_value)
    record_check(sprintf("6.5 %s, k = %d: rho1 = 0 is a proper restriction", prior_name, k_value),
                 abs(rows$density_near_0[3] / rows$density_near_0[2] - 1) < 0.1 &&
                   rows$median_r_squared_others[3] > 0.7 * rows$median_r_squared_others[1],
                 sprintf("density near 0: %s; others' median R^2: %s",
                         paste(sprintf("%.3f", rows$density_near_0), collapse = " "),
                         paste(sprintf("%.3f", rows$median_r_squared_others), collapse = " ")))
  }
  rows <- filter(near_zero, prior == "flat_r", k == k_value)
  record_check(sprintf("6.6 flat_r, k = %d: rho1 = 0 drags the others to 0", k_value),
               all(diff(rows$density_near_0) > 0) && rows$median_r_squared_others[3] < 0.01,
               sprintf("density near 0: %s; others' median R^2: %s",
                       paste(sprintf("%.3f", rows$density_near_0), collapse = " "),
                       paste(sprintf("%.3f", rows$median_r_squared_others), collapse = " ")))
  row <- filter(near_zero, prior == "flat_partial_r", k == k_value, eps == 0.02)
  critical <- 1.63 / sqrt(row$kept)
  record_check(sprintf("6.7 flat_partial_r, k = %d: the restriction gives the others Beta((k-1)/2, 3/2)", k_value),
               row$ks_beta_3_2 < critical,
               sprintf("KS distance %.4f, critical value %.4f (%d draws with |rho1| < .02)",
                       row$ks_beta_3_2, critical, row$kept))
}


## ---- Figure: the conditional and joint tests against Savage-Dickey -------------
save_figure(
  ggplot(per_dataset, aes(savage_dickey, conditional, colour = factor(n), shape = estimable)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_point(size = 2.2) +
    facet_grid(prior ~ k, labeller = labeller(k = function(v) paste("k =", v))) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "estimable") +
    labs(x = "Savage-Dickey log BF10 under the joint model", y = "conditional test, log BF10",
         colour = "n", title = "The conditional test against one prior for the whole model"),
  "study_06_coherence", width = 10, height = 8)

save_study_results("06", list(per_dataset = per_dataset, coherence = coherence, near_zero = near_zero))
report_checks()
