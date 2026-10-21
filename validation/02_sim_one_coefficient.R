## =====================================================================
##  Validation 02: one coefficient at a time
##  JZS default vs a flat prior on the partial correlation (the
##  one-predictor case of the possible-range idea), two predictors
##  (x and a covariate z).
##
##   A. simulation   true effect x sample size x predictor correlation
##   B. summaries    misleading / correct evidence, coverage, bias
##   C. confusion and calibration in different "worlds" of true effects
##   D. precision matching: which null threshold makes the flat-r null
##      calls as reliable as JZS's?
##
##  Run from the project folder. Section A runs in parallel and is cached:
##  it only reruns when its settings change (or RERUN_SIMULATION = TRUE).
## =====================================================================
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(BayesFactor)

source("validation/settings.R")   # parallel cores + cached simulation runner

results_file     <- file.path(results_folder, "sim_one_coefficient.rds")
RERUN_SIMULATION <- FALSE   # TRUE forces a rerun even if cached results match


## ===== A. SIMULATION =====

## ---- Settings ---------------------------------------------------
EFFECTS <- c(0, 0.1, 0.2, 0.3, 0.5)   # true beta of x (population SDs all 1)
NS      <- c(20, 50, 100, 200, 500)
RHO_XZ  <- c(0, 0.5)                  # correlation between x and covariate z
REPS    <- 200
GAMMA   <- 0.5                        # effect of covariate z
RSCALE  <- sqrt(2) / 4                # BayesFactor 'medium' for regression

G_B <- seq(-2, 2, length.out = 4001);              H_B <- G_B[2] - G_B[1]   # beta grid (JZS)
G_R <- seq(-1, 1, length.out = 8003)[-c(1, 8003)]; H_R <- G_R[2] - G_R[1]   # rho grid

## True partial correlation of x and y given z
true_partial_r <- function(beta, rxz) {
  se2 <- 1 - (beta^2 + GAMMA^2 + 2 * beta * GAMMA * rxz)
  beta * sqrt(1 - rxz^2) / sqrt(beta^2 * (1 - rxz^2) + se2)
}

## ---- One dataset, all methods -----------------------------------
sim_one <- function(beta, n, rxz) {
  z    <- rnorm(n)
  x    <- rxz * z + sqrt(1 - rxz^2) * rnorm(n)
  se_e <- sqrt(1 - (beta^2 + GAMMA^2 + 2 * beta * GAMMA * rxz))   # Var(y) = 1
  y    <- beta * x + GAMMA * z + rnorm(n, sd = se_e)

  fit <- lm(y ~ x + z)
  cf  <- summary(fit)$coefficients
  bh  <- cf["x", 1]; se <- cf["x", 2]; dfr <- fit$df.residual
  R2f <- summary(fit)$r.squared
  R2r <- summary(lm(y ~ z))$r.squared
  q   <- qt(0.975, dfr)

  ## 1. JZS default: exact BF; posterior approximated by Cauchy on beta
  BF_jzs <- linearReg.R2stat(n, 2, R2f, rscale = RSCALE, simple = TRUE) /
            linearReg.R2stat(n, 1, R2r, rscale = RSCALE, simple = TRUE)
  Lb   <- dt((bh - G_B) / se, dfr) / se
  post <- Lb * dcauchy(G_B, 0, RSCALE)
  cdf  <- cumsum(post) / sum(post)
  jzs  <- tibble(method = "JZS default", BF10 = BF_jzs,
                 pmean = sum(G_B * post) / sum(post),
                 lo = G_B[which.max(cdf >= 0.025)],
                 hi = G_B[which.max(cdf >= 0.975)])

  ## 2. Flat prior on the partial correlation rho in (-1, 1)
  ##    Likelihood evaluated on beta (exact t form); rho <-> beta mapping
  c_xz   <- sd(resid(lm(x ~ z))) / summary(fit)$sigma      # SD(x | z) / residual SD
  b_of_r <- G_R / (c_xz * sqrt(1 - G_R^2))                  # beta for each rho
  Lr     <- dt((bh - b_of_r) / se, dfr) / se                # beta-likelihood
  m1     <- sum(Lr * 0.5) * H_R                             # flat prior on rho: density 1/2
  m0     <- dt(bh / se, dfr) / se                           # beta = 0
  pr     <- Lr / sum(Lr); cr <- cumsum(pr)
  i_lo   <- which.max(cr >= 0.025); i_hi <- which.max(cr >= 0.975)
  prf    <- tibble(method = "partial r flat", BF10 = m1 / m0,
                   pmean  = sum(b_of_r * pr),
                   lo     = b_of_r[i_lo], hi = b_of_r[i_hi],
                   r_mean = sum(G_R * pr), r_lo = G_R[i_lo], r_hi = G_R[i_hi])

  ## 3. Frequentist reference
  freq <- tibble(method = "frequentist", BF10 = NA_real_,
                 pmean = bh, lo = bh - q * se, hi = bh + q * se)

  bind_rows(jzs, prf, freq) |> mutate(p = cf["x", 4])
}

## ---- Run the grid (parallel; cached until the settings change) -----
design <- expand_grid(beta = EFFECTS, n = NS, rho_xz = RHO_XZ)

sims <- run_simulation_grid(
  design, results_file = results_file, rerun = RERUN_SIMULATION,
  settings = list(REPS = REPS, GAMMA = GAMMA, RSCALE = RSCALE),
  simulate_cell = function(cell) {
    map_dfr(seq_len(REPS), function(i)
      mutate(sim_one(cell$beta, cell$n, cell$rho_xz), rep = i)) |>
      mutate(beta = cell$beta, n = cell$n, rho_xz = cell$rho_xz, .before = 1)
  }) |>
  mutate(true_r = true_partial_r(beta, rho_xz))

## ===== B. SUMMARIES =====

## ---- Summaries --------------------------------------------------
summ <- sims |>
  group_by(beta, n, rho_xz, method) |>
  summarise(
    BF10_gt3   = mean(BF10 > 3),
    BF10_gt10  = mean(BF10 > 10),
    BF01_gt3   = mean(BF10 < 1/3),
    BF01_gt10  = mean(BF10 < 1/10),
    p_sig      = mean(p < .05),
    coverage   = mean(lo <= beta & beta <= hi),
    r_coverage = mean(r_lo <= true_r & true_r <= r_hi),   # partial r method only
    width      = mean(hi - lo),
    bias       = mean(pmean - beta),
    rmse       = sqrt(mean((pmean - beta)^2)),
    .groups = "drop") |>
  mutate(
    mislead_3  = ifelse(beta == 0, BF10_gt3,  BF01_gt3),
    mislead_10 = ifelse(beta == 0, BF10_gt10, BF01_gt10),
    correct_3  = ifelse(beta == 0, BF01_gt3,  BF10_gt3),
    correct_10 = ifelse(beta == 0, BF01_gt10, BF10_gt10))

## Check the fix: at n = 500 the partial-r interval should match the frequentist one
summ |>
  filter(n == 500, method != "JZS default") |>
  select(rho_xz, beta, method, coverage, r_coverage, width) |>
  arrange(rho_xz, beta, method) |>
  print(n = Inf)

## Testing
summ |>
  filter(method != "frequentist") |>
  select(rho_xz, beta, n, method, mislead_3, mislead_10, correct_3, correct_10) |>
  arrange(rho_xz, beta, n, method) |>
  print(n = Inf)

## Estimation
summ |>
  select(rho_xz, beta, n, method, coverage, r_coverage, width, bias, rmse, p_sig) |>
  arrange(rho_xz, beta, n, method) |>
  print(n = Inf)

## ---- Plots ------------------------------------------------------
cols  <- c(`JZS default` = "grey40", `partial r flat` = "#1b9e77", frequentist = "black")
bayes <- filter(summ, method != "frequentist")
facet <- facet_grid(rho_xz ~ beta, labeller = label_both)
xax   <- scale_x_log10(breaks = NS)

# 1. Misleading evidence
ggplot(bayes, aes(n, mislead_3, colour = method)) +
  geom_line(linewidth = 0.9) + geom_point() + xax + facet +
  scale_colour_manual(values = cols) +
  labs(y = "P(BF > 3 for the wrong hypothesis)", title = "Misleading evidence")

# 2. Correct evidence, with p < .05 for reference
ggplot(bayes, aes(n, correct_3, colour = method)) +
  geom_line(linewidth = 0.9) + geom_point() +
  geom_line(data = filter(summ, method == "frequentist"), aes(y = p_sig),
            colour = "black", linetype = "dashed") +
  xax + facet + scale_colour_manual(values = cols) +
  labs(y = "P(BF > 3 for the true hypothesis)", title = "Correct evidence",
       subtitle = "dashed = P(p < .05)")

# 3. Coverage
ggplot(summ, aes(n, coverage, colour = method)) +
  geom_hline(yintercept = 0.95, linetype = "dotted") +
  geom_line(linewidth = 0.9) + geom_point() + xax + facet +
  scale_colour_manual(values = cols) +
  labs(y = "coverage of beta", title = "Coverage of 95% intervals")

# 4. Bias
ggplot(summ, aes(n, bias, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 0.9) + geom_point() + xax + facet +
  scale_colour_manual(values = cols) +
  labs(y = "posterior mean - true beta", title = "Bias")


## ===== C. CONFUSION AND CALIBRATION =====

sims <- readRDS(results_file)   # cached results from section A

## How often each true effect occurs in the world you study (rows sum to 1)
worlds <- tribble(
  ~world,               ~`0`, ~`0.1`, ~`0.2`, ~`0.3`, ~`0.5`,
  "balanced",            0.5,  0.125,  0.125,  0.125,  0.125,
  "small effects common",0.5,  0.25,   0.15,   0.07,   0.03,
  "large effects common",0.5,  0.03,   0.07,   0.15,   0.25,
  "few true effects",    0.8,  0.05,   0.05,   0.05,   0.05) |>
  pivot_longer(-world, names_to = "beta", values_to = "w") |>
  mutate(beta = as.numeric(beta))

bayes <- sims |>
  filter(method != "frequentist") |>
  mutate(truth = ifelse(beta == 0, "null", "effect"),
         call  = case_when(BF10 > 3   ~ "effect",
                           BF10 < 1/3 ~ "null",
                           TRUE       ~ "undecided")) |>
  group_by(beta, n, rho_xz, method) |>
  mutate(w_run = 1 / n()) |>          # each cell sums to 1
  ungroup()

weighted <- bayes |>
  inner_join(worlds, by = "beta", relationship = "many-to-many") |>
  mutate(wt = w * w_run)

## ---- 1. Confusion table per world, method, n --------------------
confusion <- weighted |>
  group_by(world, rho_xz, n, method) |>
  summarise(
    says_null       = sum(wt[call == "null"]),
    says_effect     = sum(wt[call == "effect"]),
    null_precision  = sum(wt[call == "null"   & truth == "null"])   / says_null,
    effect_precision= sum(wt[call == "effect" & truth == "effect"]) / says_effect,
    decisive        = says_null + says_effect,
    accuracy_decided= sum(wt[call == truth]) / decisive,
    wrong_null      = sum(wt[call == "null"   & truth == "effect"]),   # missed effects
    wrong_effect    = sum(wt[call == "effect" & truth == "null"]),     # false alarms
    .groups = "drop") |>
  mutate(null_lean = log(wrong_null / pmax(wrong_effect, 1e-4)))       # >0: errors lean to null

confusion |>
  select(world, rho_xz, n, method, decisive, null_precision,
         effect_precision, accuracy_decided, wrong_null, wrong_effect) |>
  arrange(world, rho_xz, n, method) |>
  print(n = Inf)

## ---- 2. Calibration: when a method says BF = k, are the odds k:1? -
## With 50% nulls the prior odds are 1, so observed odds of effect vs
## null should equal the stated BF if the method is calibrated
## for that world.
calib <- weighted |>
  filter(world != "few true effects") |>
  mutate(logBF = log(BF10),
         bin   = cut(logBF, breaks = c(-Inf, seq(-4, 4, by = 1), Inf))) |>
  group_by(world, method, bin) |>
  summarise(stated   = weighted.mean(pmin(pmax(logBF, -5), 5), wt),
            observed = log(sum(wt[truth == "effect"]) / sum(wt[truth == "null"])),
            mass     = sum(wt), .groups = "drop") |>
  filter(is.finite(observed), mass > 0.002)

ggplot(calib, aes(stated, observed, colour = method)) +
  geom_abline(linetype = "dashed") +
  geom_line(linewidth = 0.9) + geom_point(aes(size = mass)) +
  scale_colour_manual(values = c(`JZS default` = "grey40", `partial r flat` = "#1b9e77")) +
  facet_wrap(~ world) +
  labs(x = "stated log BF10", y = "observed log odds (effect : null)",
       size = "share of runs",
       title = "Calibration: does the Bayes factor mean what it says?",
       subtitle = "on the dashed line = calibrated; above = understates evidence for effects; below = overstates it")


## ===== D. PRECISION MATCHING =====

K_GRID <- 10^seq(0, 2, by = 0.01)   # candidate BF01 thresholds, 1 to 100

## Null precision and number of calls for one method at null threshold k
## (effect calls stay at BF10 > 3)
prec_at <- function(df, k) {
  null_call <- df$BF10 < 1 / k
  eff_call  <- df$BF10 > 3
  tibble(k         = k,
         says_null = sum(df$wt[null_call]),
         null_prec = sum(df$wt[null_call & df$truth == "null"]) / sum(df$wt[null_call]),
         decisive  = sum(df$wt[null_call | eff_call]))
}

matched <- weighted |>
  group_by(world, rho_xz, n) |>
  group_modify(function(d, key) {
    jzs    <- filter(d, method == "JZS default")
    pr     <- filter(d, method == "partial r flat")
    target <- prec_at(jzs, 3)
    empty  <- tibble(k = NA_real_, says_null = NA_real_,
                     null_prec = NA_real_, decisive = NA_real_)
    best <- if (is.na(target$null_prec)) empty else {
      curve <- map_dfr(K_GRID, ~ prec_at(pr, .x))
      hit   <- filter(curve, !is.na(null_prec), null_prec >= target$null_prec)
      if (nrow(hit)) slice_min(hit, k, n = 1, with_ties = FALSE) else empty
    }
    tibble(jzs_null_prec = target$null_prec,
           jzs_says_null = target$says_null,
           jzs_decisive  = target$decisive,
           k_needed      = best$k,
           pr_null_prec  = best$null_prec,
           pr_says_null  = best$says_null,
           pr_decisive   = best$decisive)
  }) |>
  ungroup() |>
  mutate(extra_null_calls = pr_says_null - jzs_says_null,
         extra_decisive   = pr_decisive  - jzs_decisive)

matched |>
  mutate(across(where(is.double), ~ round(.x, 3))) |>
  select(world, rho_xz, n, k_needed, jzs_null_prec, pr_null_prec,
         jzs_says_null, pr_says_null, extra_null_calls, extra_decisive) |>
  arrange(world, rho_xz, n) |>
  print(n = Inf)

## Plot: extra null calls at equal reliability
ggplot(matched, aes(n, extra_null_calls, colour = factor(rho_xz))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 0.9) + geom_point() +
  scale_x_log10(breaks = c(20, 50, 100, 200, 500)) +
  facet_wrap(~ world) +
  labs(x = "n", y = "extra null calls (share of all studies)",
       colour = "rho_xz",
       title = "Your null calls, made as reliable as JZS's",
       subtitle = "above 0: same reliability, more conclusions; below 0: fewer")

## Plot: the threshold needed
ggplot(matched, aes(n, k_needed, colour = factor(rho_xz))) +
  geom_hline(yintercept = 3, linetype = "dotted") +
  geom_line(linewidth = 0.9) + geom_point() +
  scale_x_log10(breaks = c(20, 50, 100, 200, 500)) +
  scale_y_log10() +
  facet_wrap(~ world) +
  labs(x = "n", y = "BF01 threshold needed", colour = "rho_xz",
       title = "Null threshold that matches JZS's null reliability",
       subtitle = "dotted = the usual 3")
