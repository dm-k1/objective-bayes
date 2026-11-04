## =====================================================================
##  Precision targets: how many simulated datasets each study needs
##
##  Fixed before anything is run (validation/ADEMP.md; Morris, White &
##  Crowther, 2019, Table 6). For a proportion p the Monte Carlo standard
##  error is sqrt(p (1 - p) / n_sim), so n_sim = p (1 - p) / MCSE^2, taken at
##  the worst case p = .5 when p is not known in advance.
## =====================================================================

## Monte Carlo standard error of an estimated proportion. It uses the
## Agresti-Coull adjusted proportion (x + 2) / (n_sim + 4), so a share of 0 or 1
## does not get an SE of zero.
mcse_of_proportion <- function(p, n_sim) {
  adjusted <- (p * n_sim + 2) / (n_sim + 4)
  sqrt(adjusted * (1 - adjusted) / (n_sim + 4))
}

## Datasets needed for a proportion near p to reach a target MCSE
n_sim_for_proportion <- function(target_mcse, p = 0.5) ceiling(p * (1 - p) / target_mcse^2)

## Datasets needed for a mean with standard deviation sd to reach a target MCSE
n_sim_for_mean <- function(target_mcse, sd) ceiling((sd / target_mcse)^2)

## Effective draws needed for a Kolmogorov-Smirnov check whose critical
## distance at the .01 level (1.63 / sqrt(n)) is at most a given distance
n_sim_for_ks <- function(critical_distance) ceiling((1.63 / critical_distance)^2)

## The targets of validation/ADEMP.md, and the numbers of datasets that follow.
## n_sim is what the studies run; where it exceeds the formula, the reason is
## given.
precision_targets <- list(
  ## a share (rate of an answer, inclusion) reported per design cell: MCSE .02
  rate_per_cell     = list(mcse = 0.02, n_sim = n_sim_for_proportion(0.02)),               # 625
  ## interval coverage near .95, per cell: MCSE .01
  coverage_per_cell = list(mcse = 0.01, n_sim = n_sim_for_proportion(0.01, p = 0.95)),     # 475
  ## KS checks of a distribution at the .01 level with critical distance .01
  ks_draws          = list(mcse = NA, n_sim = n_sim_for_ks(0.01)),                          # 26,569
  ## a mean log BF difference (coherence), SD about 0.30 in the last run:
  ## (0.30 / 0.05)^2 = 36 estimable datasets; 50 per cell because some fits
  ## are not estimable
  coherence_mean    = list(mcse = 0.05, n_sim = 50),
  ## simulation-based calibration: fits per prior and setting (power against
  ## miscalibration of about 5% of the rank distribution)
  sbc_fits          = list(mcse = NA, n_sim = 1000),
  ## calibration of the Bayes factor: datasets per cell, so that bins of at
  ## least 625 datasets have shares with MCSE <= .02
  bf_calibration    = list(mcse = 0.02, n_sim = 20000),
  ## slopes of the median log BF against log n: R-squared draws per setting
  ## (slope MCSE about 0.02)
  rate_draws        = list(mcse = 0.02, n_sim = 2000),
  ## growing k with a small true model: datasets per setting (bounded by the
  ## cost of the k = 500 regressions at n = 2,000)
  growing_k_data    = list(mcse = NA, n_sim = 400))

## Smoke mode (VALIDATION_SMOKE=1): every study runs end to end on a handful
## of datasets and a corner of its design, to show the code works before the
## full run. Results in smoke mode are not results.
smoke_mode <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
if (smoke_mode) {   # smoke-test figures stay out of the tracked results
  figures_folder <- file.path(cache_folder, "smoke_figures")
  dir.create(figures_folder, recursive = TRUE, showWarnings = FALSE)
}
n_sim_of <- function(target) if (smoke_mode) 10L else precision_targets[[target]]$n_sim
