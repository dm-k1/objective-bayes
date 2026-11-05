## =====================================================================
##  Study 3: is the posterior calibrated?
##  (docs/appendix/C_study_design.md, study 3; replaces the SBC part of the old study 01)
##
##  Aims. Simulation-based calibration (Talts et al., 2018) of flatcor's Stan
##  model under each prior, with the SBC package: if the sampler is right, the
##  rank of each true value among its posterior draws is uniform.
##  Data. Set B (R/data/data_set_B.R): parameters and data drawn from each
##  prior, with the model's proper stand-ins for the intercept (N(0, 10)) and
##  the noise SD (lognormal(0, 1)); k = 8, n = 60, and k = 3, n = 50.
##  Estimands. The ranks of the true R-squared, the first coefficient b[1] and
##  the noise SD sigma among the posterior draws.
##  Methods. flatcor's Stan model (system.file("stan", "flatcor.stan")) under
##  each prior; as a contrast expected to fail, flat R sampled directly in the
##  coefficients (its density has a peak at zero the sampler misses).
##  Performance. The SBC package computes the ranks (thinned so that each fit
##  contributes 100 ranks) and draws the ECDF-difference plots. Uniformity is
##  judged per quantity by a chi-square over 20 rank bins (Talts et al.).
##  n_sim. 1,000 fits per prior and setting (precision.R, sbc_fits).
##
##  Checks:
##    3.1 (old S) Under each prior and setting, the ranks of R-squared, b[1] and
##        sigma are uniform: chi-square over 20 bins, p > .01.
##    The coefficient-coordinates contrast is reported, not checked.
##
##  Run from the project folder:  Rscript R/studies/03_posterior_calibration.R
##  Smoke test: VALIDATION_SMOKE=1 Rscript R/studies/03_posterior_calibration.R
##  Raw results are cached in cache/study_03/flatcor_<build>/ (delete it to rerun); the
##  tables and checks are saved to results/studies/study_03.rds.
## =====================================================================
library(SBC)
library(patchwork)
library(dplyr)
library(ggplot2)
source("R/shared/settings.R")             # flatcor, cores, record_check(), save_figure()
source("R/shared/comparison_methods.R")   # contenders
source("R/shared/precision.R")            # n_sim_of()
source("R/data/data_sets.R")            # with_stream_of()
source("R/data/data_set_B.R")           # generate_B()

study_folder <- study_cache_folder("03")
future::plan(future::multisession, workers = number_of_cores)

n_fits <- n_sim_of("sbc_fits")
settings <- data.frame(k = c(8, 3), n = c(60, 50))
if (smoke_mode) settings <- settings[1, , drop = FALSE]
iterations <- if (smoke_mode) c(warmup = 300, sampling = 300) else c(warmup = 1000, sampling = 1000)


## ---- The models -------------------------------------------------------------------
## the executables go to this study's cache, never into the installed package
## (read-only in an Apptainer image)
direction_model <- cmdstanr::cmdstan_model(system.file("stan", "flatcor.stan", package = "flatcor"),
                                           dir = study_folder)
coefficients_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file("
functions {
  // the flat R prior as a density on the coefficients: P(total) x P(combination | total) x conversion
  real flat_r_lpdf(vector b, real sigma, matrix predictor_covariance) {
    int k = rows(b);
    real f_squared = quad_form(predictor_covariance, b) / square(sigma);
    real f = sqrt(f_squared);
    return -1.5 * log1p(f_squared)
           - (log(2) + 0.5 * k * log(pi()) - lgamma(0.5 * k) + (k - 1) * log(f))
           + 0.5 * log_determinant(predictor_covariance) - k * log(sigma);
  }
}
data {
  int<lower=1> N; int<lower=1> K; matrix[N, K] X; vector[N] y;
  int<lower=0, upper=1> prior_only; int<lower=0, upper=1> proper_stand_ins;
  real<lower=0> r_squared_shape;
}
transformed data {
  matrix[N, K] Xc;
  for (j in 1:K) Xc[, j] = X[, j] - mean(X[, j]);
  matrix[K, K] predictor_covariance = crossprod(Xc) / N;
}
parameters { real alpha; real<lower=0> sigma; vector[K] b; }
model {
  alpha ~ normal(0, 10);
  sigma ~ lognormal(0, 1);
  target += flat_r_lpdf(b | sigma, predictor_covariance);
  if (!prior_only) y ~ normal(alpha + Xc * b, sigma);
}
generated quantities {
  real R_squared = quad_form(predictor_covariance, b) / (quad_form(predictor_covariance, b) + square(sigma));
}"))


## ---- SBC per prior and setting ----------------------------------------------------
## One SBC dataset: the true values (R-squared, b[1], sigma) and the Stan data
sbc_dataset <- function(prior, k, n, replication) {
  dat <- generate_B(prior, k, n, x1 = "prior", purpose = "sbc", replication = replication)
  list(variables = list(R_squared = dat$r_squared, `b[1]` = dat$coefficients[1], sigma = dat$noise_sd),
       generated = list(N = n, K = k, X = dat$predictors, y = dat$outcome, prior_only = 0L,
                        proper_stand_ins = 1L, r_squared_shape = prior_shape(prior, k)))
}

run_sbc <- function(model, prior, k, n, label) {
  cache_file <- file.path(study_folder, sprintf("sbc_%s_k%d_n%d.rds", label, k, n))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  datasets <- do.call(bind_datasets, lapply(seq_len(n_fits), function(replication) {
    one <- sbc_dataset(prior, k, n, replication)
    SBC_datasets(variables = posterior::as_draws_matrix(posterior::as_draws_rvars(
                   lapply(one$variables, posterior::rvar))), generated = list(one$generated))
  }))
  backend <- SBC_backend_cmdstan_sample(model, chains = 1, iter_warmup = iterations[["warmup"]],
                                        iter_sampling = iterations[["sampling"]], seed = 3,
                                        refresh = 0, show_messages = FALSE)
  result <- compute_SBC(datasets, backend, keep_fits = FALSE,
                        thin_ranks = iterations[["sampling"]] / 100)
  stats <- result$stats |> mutate(prior = label, k = k, n = n)
  saveRDS(stats, cache_file)
  stats
}

stats <- bind_rows(c(
  lapply(seq_len(nrow(settings)), function(row) bind_rows(lapply(unname(contenders), function(prior)
    run_sbc(direction_model, prior, settings$k[row], settings$n[row], prior)))),
  list(run_sbc(coefficients_model, "flat_r", 8, 60, "flat_r in coefficient coordinates"))))


## ---- Uniformity of the ranks, and the checks ----------------------------------------
uniformity <- stats |>
  group_by(prior, k, n, variable) |>
  summarise(fits = n(), p_value = {
              bins <- table(cut(rank / max_rank, breaks = seq(0, 1, length.out = 21), include.lowest = TRUE))
              suppressWarnings(chisq.test(bins)$p.value)
            }, .groups = "drop")
cat("\nRank uniformity (chi-square over 20 bins):\n")
print(as.data.frame(uniformity), digits = 3, row.names = FALSE)

for (row in seq_len(nrow(uniformity))) with(uniformity[row, ], {
  if (prior == "flat_r in coefficient coordinates")
    message(sprintf("contrast (not a check): %s, k = %d, %s: p = %.3g", prior, k, variable, p_value))
  else
    record_check(sprintf("3.1 %s, k = %d, n = %d: ranks of %s uniform", prior, k, n, variable),
                 p_value > 0.01, sprintf("chi-square p = %.3f over %d fits", p_value, fits))
})

## one panel row per prior (plot_ecdf_diff() needs one set of datasets at a time)
save_figure(
  patchwork::wrap_plots(lapply(unique(stats$prior), function(prior_name)
    plot_ecdf_diff(as.data.frame(filter(stats, prior == prior_name, k == settings$k[1]))) +
      labs(title = prior_name)), ncol = 1) +
    patchwork::plot_annotation(title = sprintf(
      "Simulation-based calibration, k = %d, n = %d: ECDF difference of the ranks",
      settings$k[1], settings$n[1])),
  "study_03_ecdf_difference", width = 10, height = 14)

save_study_results("03", list(stats = stats, uniformity = uniformity))
report_checks()
