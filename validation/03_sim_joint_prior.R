## =====================================================================
##  Validation 03: the joint possible-range prior vs alternatives
##  4 correlated predictors; x1's effect varies, x3 is always null.
##  Methods: JZS default, joint possible-range prior (this project),
##  per-coefficient flat partial r, R2-D2 with the same R-squared prior.
##
##  Run from the project folder. The simulation runs in parallel and is
##  cached: it only reruns when its settings change (or RERUN_SIMULATION = TRUE).
## =====================================================================
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
source("R/possible_range.R")
source("validation/comparison_methods.R")

source("validation/settings.R")   # parallel cores + cached simulation runner

results_file     <- file.path(results_folder, "sim_joint_prior.rds")
RERUN_SIMULATION <- FALSE   # TRUE forces a rerun even if cached results match

## ---- Settings --------------------------------------------------------
K_PRED  <- 4
EFFECTS <- c(0, 0.1, 0.2, 0.3, 0.5)
OTHER   <- c(0.3, 0, 0)
TEST    <- c(1, 3)
NS      <- c(20, 50, 100, 200, 500)
RHO     <- c(0, 0.5)
REPS    <- 200
RSCALE  <- sqrt(2) / 4

set.seed(99)
draws_k  <- draw_r2d2_prior(K_PRED)       # shared prior draws (full model)
draws_k1 <- draw_r2d2_prior(K_PRED - 1)   # shared prior draws (reduced models)

## ---- Monte Carlo check: two independent draw sets, same dataset -----
set.seed(5)
Xt <- matrix(rnorm(100 * 4), 100); yt <- drop(Xt %*% c(0.2, 0.3, 0, 0)) + rnorm(100)
Xtc <- scale(Xt, scale = FALSE); ytc <- yt - mean(yt)
c(draws_A = log_bayes_factor_r2d2(Xtc, ytc, draw_r2d2_prior(4)),
  draws_B = log_bayes_factor_r2d2(Xtc, ytc, draw_r2d2_prior(4)))   # should differ by < ~0.1

## ---- One dataset -----------------------------------------------------
sim_one <- function(b1, n, rho) {
  Sig  <- matrix(rho, K_PRED, K_PRED); diag(Sig) <- 1
  beta <- c(b1, OTHER)
  X    <- matrix(rnorm(n * K_PRED), n) %*% chol(Sig)
  y    <- drop(X %*% beta) + rnorm(n, sd = sqrt(1 - drop(t(beta) %*% Sig %*% beta)))
  Xc   <- scale(X, scale = FALSE); yc <- y - mean(y)

  fit    <- lm(y ~ X)
  R2f    <- summary(fit)$r.squared
  lbf_jf <- log_bayes_factor_vs_intercept_only(R2f, n, K_PRED)
  lbf_rf <- log_bayes_factor_r2d2(Xc, yc, draws_k)

  map_dfr(TEST, function(j) {
    R2r <- summary(lm(y ~ X[, -j]))$r.squared
    tibble(coef  = paste0("x", j), truth = beta[j],
           method = c("JZS default", "joint possible",
                      "per-coefficient flat r", "R2-D2 matched"),
           BF10  = c(bayes_factor_jzs(R2f, R2r, n, K_PRED, rscale = RSCALE),
                     exp(lbf_jf - log_bayes_factor_vs_intercept_only(R2r, n, K_PRED - 1)),
                     bayes_factor_flat_partial_r(fit, j, X),
                     exp(lbf_rf - log_bayes_factor_r2d2(Xc[, -j, drop = FALSE], yc, draws_k1))))
  })
}

## ---- Run (parallel; cached until the settings change) -----------------
sims <- run_simulation_grid(
  expand_grid(b1 = EFFECTS, n = NS, rho = RHO),
  results_file = results_file, rerun = RERUN_SIMULATION,
  settings = list(REPS = REPS, K_PRED = K_PRED, OTHER = OTHER, TEST = TEST, RSCALE = RSCALE),
  simulate_cell = function(cell) {
    map_dfr(seq_len(REPS), function(i)
      mutate(sim_one(cell$b1, cell$n, cell$rho), rep = i)) |>
      mutate(b1 = cell$b1, n = cell$n, rho = cell$rho, .before = 1)
  })

## ---- Summaries ---------------------------------------------------------
summ <- sims |>
  mutate(correct        = ifelse(truth == 0, BF10 < 1/3, BF10 > 3),
         mislead        = ifelse(truth == 0, BF10 > 3,  BF10 < 1/3),
         mislead_strong = ifelse(truth == 0, BF10 > 10, BF10 < 1/10)) |>
  group_by(coef, rho, b1, n, method) |>
  summarise(correct_3  = mean(correct),
            mislead_3  = mean(mislead),
            mislead_10 = mean(mislead_strong),
            n_runs     = n(), .groups = "drop")

## x1: effect varies
summ |> filter(coef == "x1") |>
  select(rho, b1, n, method, correct_3, mislead_3, mislead_10) |>
  arrange(rho, b1, n, method) |> print(n = Inf)

## x3: always null (pooled over x1 effect sizes)
sims |> filter(coef == "x3") |>
  group_by(rho, n, method) |>
  summarise(null_found   = mean(BF10 < 1/3),
            false_effect = mean(BF10 > 3),
            false_strong = mean(BF10 > 10), .groups = "drop") |>
  arrange(rho, n, method) |> print(n = Inf)

cols <- c(`JZS default` = "grey40", `joint possible` = "#d95f02",
          `per-coefficient flat r` = "#1b9e77", `R2-D2 matched` = "#7570b3")

ggplot(filter(summ, coef == "x1"), aes(n, correct_3, colour = method)) +
  geom_line(linewidth = 0.9) + geom_point() +
  scale_x_log10(breaks = NS) + scale_colour_manual(values = cols) +
  facet_grid(rho ~ b1, labeller = label_both) +
  labs(y = "P(BF > 3 for the true hypothesis)", title = "x1: correct evidence")

ggplot(filter(summ, coef == "x1"), aes(n, mislead_3, colour = method)) +
  geom_line(linewidth = 0.9) + geom_point() +
  scale_x_log10(breaks = NS) + scale_colour_manual(values = cols) +
  facet_grid(rho ~ b1, labeller = label_both) +
  labs(y = "P(BF > 3 for the wrong hypothesis)", title = "x1: misleading evidence")
