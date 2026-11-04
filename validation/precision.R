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

## The targets of validation/ADEMP.md, and the numbers of datasets that follow
precision_targets <- list(
  rate_per_cell     = list(mcse = 0.02, n_sim = n_sim_for_proportion(0.02)),             # 625
  coverage_per_cell = list(mcse = 0.01, n_sim = n_sim_for_proportion(0.01, p = 0.95)),   # 475
  coherence_mean    = list(mcse = 0.05, n_sim = n_sim_for_mean(0.05, sd = 0.30)))       # 36

## Smoke mode (VALIDATION_SMOKE=1): every study runs end to end on a handful
## of datasets and a corner of its design, to show the code works before the
## full run. Results in smoke mode are not results.
smoke_mode <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
n_sim_of <- function(target) if (smoke_mode) 10L else precision_targets[[target]]$n_sim
