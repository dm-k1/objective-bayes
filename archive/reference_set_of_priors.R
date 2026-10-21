## ARCHIVED: earlier approach (reference set of Cauchy priors scored at
## one-observation resolution). Abandoned: the one-observation sigma region
## did not work as a definition of the possible range. Kept for history only.

## =================================================================
##  Reference set of priors at one-observation resolution
##  1. Pilot fit           -> b_hat, SE, sigma = SE * sqrt(N)
##  2. Grid of Cauchy priors (mu, s), in sigma units
##  3. Score each prior on predicting ONE observation
##  4. Keep priors with support >= 1/k  (k = 8 main, 32 comparison)
##  5. Full-data posterior for every kept prior (MCMC + integration)
##  6. Combine with equal weight; compare with the frequentist interval
## =================================================================

library(cmdstanr)
library(posterior)
library(loo)
library(dplyr)
library(purrr)
library(ggplot2)

vec <- function(v) array(v, dim = length(v))   # safe for K = 1

## ---- Settings ---------------------------------------------------
NU        <- 1       # prior family: 1 = Cauchy
S_WIDE    <- 100     # pilot prior scale (effectively flat)
K_MAIN    <- 8       # support threshold 1/8 (Royall: fairly strong)
K_ALT     <- 32      # comparison threshold 1/32 (strong)
D_HALF    <- 6       # grid: mu = b_hat +/- D_HALF * sigma
R_MIN     <- 0.01    # grid: s from R_MIN * sigma ...
R_MAX     <- 10      # ... to R_MAX * sigma (linear spacing)
N_D       <- 121
N_R       <- 121
KHAT_MAX  <- 0.7     # PSIS reliability threshold
AGREE_MAX <- 0.25    # max MCMC vs integration gap (post. SD, and log BF)
MAX_REFIT <- 100     # cap on MCMC refits; integration fills in beyond it

## ---- 1. Data: regression with K predictors ----------------------
set.seed(1)
N <- 100
X_raw <- cbind(x1 = rnorm(N), x2 = rnorm(N))
y_raw <- as.vector(X_raw %*% c(0.4, 0) + rnorm(N))   # x1 real, x2 null

terms <- colnames(X_raw)
K     <- ncol(X_raw)
X     <- matrix(scale(X_raw), N, K)
y     <- as.vector(scale(y_raw))

## ---- 2. Stan model: prior per coefficient passed in as data ------
stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
  vector[K] m;
  vector<lower=0>[K] s;
  vector<lower=0>[K] nu;
}
parameters {
  real a;
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  a ~ normal(0, 10);
  sigma ~ exponential(1);
  beta ~ student_t(nu, m, s);
  y ~ normal(a + X * beta, sigma);
}"
mod <- cmdstan_model(write_stan_file(stan_code))

run_fit <- function(m, s, iter, seed = 1) {
  mod$sample(
    data = list(N = N, K = K, X = X, y = y,
                m = vec(m), s = vec(s), nu = vec(rep(NU, K))),
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = iter,
    adapt_delta = 0.95, max_treedepth = 12,
    init = function() list(a = 0, beta = vec(rep(0, K)), sigma = 1),
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE,
    seed = seed
  )
}

get_draws <- function(fit) {
  dr <- as_draws_matrix(fit$draws(c("a", "beta", "sigma")))
  list(A   = as.numeric(dr[, "a"]),
       B   = as.matrix(dr[, paste0("beta[", 1:K, "]")]),
       Sig = as.numeric(dr[, "sigma"]))
}

lprior         <- function(b, m, s) dt((b - m) / s, df = NU, log = TRUE) - log(s)
logsumexp_mean <- function(x) max(x) + log(mean(exp(x - max(x))))
lse            <- function(x) max(x) + log(sum(exp(x - max(x))))
trap           <- function(x, f) sum(diff(x) * (head(f, -1) + tail(f, -1)) / 2)

## ---- 3. Pilot fit: b_hat, SE, and the one-observation sigma -----
pilot <- run_fit(m = rep(0, K), s = rep(S_WIDE, K), iter = 5000)
print(pilot$summary(c("a", "beta", "sigma")))   # check rhat ~ 1, ess > 400

pd     <- get_draws(pilot)
b_hat  <- colMeans(pd$B)
SE     <- apply(pd$B, 2, sd)
sigma1 <- SE * sqrt(N)                 # uncertainty from ONE observation

print(tibble(term = terms, b_hat = b_hat, SE = SE, sigma_one_obs = sigma1))

## Pilot Bayes factors: Rao-Blackwellized Savage-Dickey
rb_density0 <- function(k, d) {
  xk   <- X[, k]; sxx <- sum(xk^2); D <- length(d$A)
  rest <- if (K > 1) X[, -k, drop = FALSE] %*% t(d$B[, -k, drop = FALSE]) else 0
  R    <- sweep(matrix(y, N, D) - rest, 2, d$A)
  mean(dnorm(0, colSums(xk * R) / sxx, d$Sig / sqrt(sxx)))
}
BF_pilot <- sapply(1:K, function(k) exp(lprior(0, 0, S_WIDE)) / rb_density0(k, pd))

## Full-data marginal of beta_k (flat prior): Student-t, N - K - 1 df
df_t  <- N - K - 1
lmarg <- function(b, k) {
  sc <- SE[k] * sqrt((df_t - 2) / df_t)
  dt((b - b_hat[k]) / sc, df = df_t, log = TRUE) - log(sc)
}

## Full-data posterior mean, SD and log BF10 for one prior (integration)
post_int <- function(k, m, s) {
  b  <- sort(c(seq(b_hat[k] - 12 * SE[k], b_hat[k] + 12 * SE[k], length.out = 401),
               m + s * sinh(seq(-12, 12, length.out = 801))))
  lm <- lmarg(b, k)
  lp <- lm + lprior(b, m, s)
  f  <- exp(lp - max(lp))
  Z  <- trap(b, f)
  pm <- trap(b, b * f) / Z
  lr <- lp - lprior(b, 0, S_WIDE)
  lratio <- max(lr) + log(trap(b, exp(lr - max(lr))))
  c(mean = pm, sd = sqrt(trap(b, (b - pm)^2 * f) / Z),
    logBF = log(BF_pilot[k]) + lratio)
}

## =================================================================
##  STEP 3-4. SCORE EVERY PRIOR ON ONE OBSERVATION, KEEP THE GOOD ONES
## =================================================================

## How much probability the prior's prediction for ONE observation
## puts on the estimate: Cauchy(m, s) convolved with Normal(0, sigma^2)
log_score1 <- function(k, m, s) {
  v  <- sort(c(m + s * sinh(seq(-12, 12, length.out = 801)),
               b_hat[k] + sigma1[k] * seq(-10, 10, length.out = 401)))
  lf <- lprior(v, m, s) + dnorm(b_hat[k], v, sigma1[k], log = TRUE)
  max(lf) + log(trap(v, exp(lf - max(lf))))
}
log_best1 <- function(k) dnorm(0, 0, sigma1[k], log = TRUE)   # point prior at b_hat

grid <- map_dfr(1:K, function(k) {
  g <- expand.grid(d = seq(-D_HALF, D_HALF, length.out = N_D),
                   r = seq(R_MIN, R_MAX, length.out = N_R))
  m <- b_hat[k] + g$d * sigma1[k]
  s <- g$r * sigma1[k]
  tibble(term = terms[k], d = g$d, r = g$r, m = m, s = s,
         log_support = map2_dbl(m, s, function(m, s) log_score1(k, m, s)) - log_best1(k))
}) |>
  mutate(keep_main = log_support >= -log(K_MAIN),
         keep_alt  = log_support >= -log(K_ALT))

## The region must close inside the grid
grid |>
  group_by(term) |>
  summarise(hits_side_alt = any(keep_alt & abs(d) == max(abs(d))),
            hits_top_alt  = any(keep_alt & r == max(r)),
            n_main = sum(keep_main), n_alt = sum(keep_alt)) |>
  print()   # hits_* FALSE = region closed; otherwise raise D_HALF / R_MAX

## =================================================================
##  STEP 5. FULL-DATA POSTERIOR FOR EVERY KEPT PRIOR
## =================================================================

kept <- filter(grid, keep_alt)        # the 1/32 set contains the 1/8 set

## 5a. Integration (exact up to grid resolution)
int <- t(mapply(function(term, m, s) post_int(match(term, terms), m, s),
                kept$term, kept$m, kept$s))
kept <- kept |>
  mutate(int_mean = int[, "mean"], int_sd = int[, "sd"], int_logBF = int[, "logBF"])

## 5b. MCMC: reweight pilot draws, chain from neighbours, refit if needed
reweight_from <- function(k, m, s, base, label) {
  bk   <- base$bk
  logw <- lprior(bk, m, s) - lprior(bk, base$m, base$s)
  ps   <- suppressWarnings(psis(logw, r_eff = 1))
  wn   <- weights(ps, log = FALSE, normalize = TRUE)
  pm   <- sum(wn * bk)
  tibble(term = terms[k], m = m, s = s,
         post_mean = pm,
         post_sd   = sqrt(sum(wn * (bk - pm)^2)),
         BF10      = base$BF * exp(logsumexp_mean(logw)),
         khat      = ps$diagnostics$pareto_k,
         rhat      = NA_real_,
         method    = label)
}

refit_base <- function(k, m, s) {
  mm <- rep(0, K); ss <- rep(S_WIDE, K); mm[k] <- m; ss[k] <- s
  fit  <- run_fit(mm, ss, iter = 2000, seed = sample.int(1e6, 1))
  bk   <- get_draws(fit)$B[, k]
  logr <- lprior(bk, 0, S_WIDE) - lprior(bk, m, s)
  list(m = m, s = s, bk = bk, rhat = fit$summary("beta")$rhat[k],
       BF = BF_pilot[k] / exp(logsumexp_mean(logr)))
}

bases <- map(1:K, function(k)
  list(list(m = 0, s = S_WIDE, bk = pd$B[, k], BF = BF_pilot[k], rhat = 1)))

mc <- pmap_dfr(select(kept, term, m, s), function(term, m, s) {
  k <- match(term, terms)
  reweight_from(k, m, s, bases[[k]][[1]], "reweight")
})

bad <- mc |> filter(khat > KHAT_MAX) |> arrange(term, m, desc(s))
message(nrow(bad), " priors need more than the pilot")

n_refit  <- 0
new_rows <- vector("list", nrow(bad))
for (i in seq_len(nrow(bad))) {
  k <- match(bad$term[i], terms); m <- bad$m[i]; s <- bad$s[i]
  good_bases <- keep(bases[[k]], ~ .x$rhat <= 1.01)
  dist <- map_dbl(good_bases, ~ ((.x$m - m) / SE[k])^2 + (log(.x$s) - log(s))^2)
  row  <- NULL
  for (j in head(order(dist), 3)) {
    r <- reweight_from(k, m, s, good_bases[[j]], "chain")
    if (r$khat <= KHAT_MAX) { row <- r; break }
  }
  if (is.null(row)) {
    if (n_refit >= MAX_REFIT) next
    b <- refit_base(k, m, s); n_refit <- n_refit + 1
    bases[[k]] <- c(bases[[k]], list(b))
    row <- tibble(term = terms[k], m = m, s = s,
                  post_mean = mean(b$bk), post_sd = sd(b$bk),
                  BF10 = b$BF, khat = NA_real_, rhat = b$rhat,
                  method = "refit")
  }
  new_rows[[i]] <- row
}
updates <- bind_rows(new_rows)
message(n_refit, " refits; ", nrow(updates), " of ", nrow(bad), " priors resolved")
if (nrow(updates) > 0)
  mc <- rows_update(mc, updates, by = c("term", "m", "s"))

## 5c. Final values: MCMC where reliable and in agreement, else integration
final <- kept |>
  left_join(mc, by = c("term", "m", "s")) |>
  mutate(d_mean  = abs(post_mean - int_mean) / int_sd,
         d_logBF = abs(log(BF10) - int_logBF),
         mcmc_ok = !is.na(post_mean) &
                   (is.na(khat) | khat <= KHAT_MAX) &
                   (is.na(rhat) | rhat <= 1.01) &
                   d_mean <= AGREE_MAX & d_logBF <= AGREE_MAX,
         source  = ifelse(mcmc_ok, method, "integration"),
         pmean   = ifelse(mcmc_ok, post_mean, int_mean),
         psd     = ifelse(mcmc_ok, post_sd,   int_sd),
         logBF   = ifelse(mcmc_ok, log(BF10), int_logBF))

final |>
  filter(!is.na(post_mean), is.na(khat) | khat <= KHAT_MAX, is.na(rhat) | rhat <= 1.01) |>
  group_by(term) |>
  summarise(median_d_mean = median(d_mean), max_d_mean = max(d_mean),
            max_d_logBF = max(d_logBF)) |>
  print()
count(final, term, source) |> print()

## =================================================================
##  STEP 6-7. EQUAL-WEIGHT COMBINATION AND COMPARISON
## =================================================================

## Exact mixture density of the kept priors' full-data posteriors
mixture_density <- function(k, ck) {
  b   <- seq(b_hat[k] - 15 * SE[k], b_hat[k] + 15 * SE[k], length.out = 6001)
  lm  <- lmarg(b, k)
  acc <- numeric(length(b))
  for (i in seq_len(nrow(ck))) {
    lp  <- lm + lprior(b, ck$m[i], ck$s[i])
    f   <- exp(lp - max(lp))
    acc <- acc + f / trap(b, f)
  }
  tibble(beta = b, dens = acc / nrow(ck))
}
mix_quant <- function(b, dens, p) {
  cdf <- cumsum(c(0, diff(b) * (head(dens, -1) + tail(dens, -1)) / 2))
  approx(cdf / max(cdf), b, xout = p, ties = mean)$y
}

combine <- function(flag, k_val) {
  map(1:K, function(k) {
    ck <- filter(final, term == terms[k], .data[[flag]])
    md <- mixture_density(k, ck)
    list(
      summary = tibble(
        term       = terms[k],
        support    = paste0("1/", k_val),
        n_priors   = nrow(ck),
        mix_mean   = mean(ck$pmean),
        within_sd  = sqrt(mean(ck$psd^2)),
        between_sd = sqrt(mean((ck$pmean - mean(ck$pmean))^2)),
        lo95       = mix_quant(md$beta, md$dens, 0.025),
        hi95       = mix_quant(md$beta, md$dens, 0.975),
        BF10_min   = exp(min(ck$logBF)),
        BF10_max   = exp(max(ck$logBF)),
        BF10_avg   = exp(lse(ck$logBF) - log(nrow(ck)))),
      density = mutate(md, term = terms[k], support = paste0("1/", k_val)))
  })
}

res_main <- combine("keep_main", K_MAIN)
res_alt  <- combine("keep_alt",  K_ALT)

results <- bind_rows(map(c(res_main, res_alt), "summary")) |>
  left_join(tibble(term = terms, b_hat = b_hat, SE = SE), by = "term") |>
  mutate(mix_sd     = sqrt(within_sd^2 + between_sd^2),
         price      = mix_sd / SE,                     # 1 = frequentist width
         width_ratio = (hi95 - lo95) / (2 * 1.96 * SE),
         freq_lo95  = b_hat - 1.96 * SE,
         freq_hi95  = b_hat + 1.96 * SE,
         P_effect   = BF10_avg / (1 + BF10_avg))

print(select(results, term, support, n_priors, mix_mean, within_sd, between_sd,
             mix_sd, SE, price))
print(select(results, term, support, lo95, hi95, freq_lo95, freq_hi95, width_ratio))
print(select(results, term, support, BF10_min, BF10_max, BF10_avg, P_effect))

## =================================================================
##  PLOTS
## =================================================================
ref_lines <- tibble(term = terms, b_hat = b_hat)
lin_y <- scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02)))

# 1. The reference set: support at one-observation resolution
ggplot(grid, aes(m, s)) +
  geom_raster(aes(fill = pmax(log_support, -log(1000)))) +
  geom_contour(aes(z = log_support, linetype = after_stat(factor(level))),
               breaks = -log(c(K_ALT, K_MAIN)), colour = "black") +
  geom_vline(data = ref_lines, aes(xintercept = b_hat), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_viridis_c(labels = function(x) paste0("1/", round(exp(-x)))) +
  scale_linetype_manual(values = c("dashed", "solid"),
                        labels = paste0("1/", c(K_ALT, K_MAIN)),
                        name = "support\nboundary") +
  lin_y +
  facet_wrap(~ term, scales = "free") +
  labs(x = "prior location (mu)", y = "prior scale (s)",
       fill = "predicts one\nobservation as\nwell as best x",
       title = "The reference set of priors",
       subtitle = "dashed = estimate, dotted = 0")

# 2. Where each kept prior's full-data posterior ends up
ggplot(filter(final, keep_main), aes(m, s)) +
  geom_raster(aes(fill = (pmean - b_hat[match(term, terms)]) / SE[match(term, terms)])) +
  geom_vline(data = ref_lines, aes(xintercept = b_hat), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_gradient2(midpoint = 0) +
  lin_y +
  facet_wrap(~ term, scales = "free") +
  labs(x = "prior location (mu)", y = "prior scale (s)",
       fill = "posterior mean\n(SE from b_hat)",
       title = "How far each prior pulls the full-data posterior")

# 3. Evidence across the reference set
ggplot(filter(final, keep_main), aes(m, s)) +
  geom_raster(aes(fill = logBF)) +
  geom_vline(data = ref_lines, aes(xintercept = b_hat), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_gradient2(midpoint = 0) +
  lin_y +
  facet_wrap(~ term, scales = "free") +
  labs(x = "prior location (mu)", y = "prior scale (s)", fill = "log BF10",
       title = "Evidence for an effect, per prior in the reference set")

# 4. The combined posterior vs the data alone
dens_all <- bind_rows(map(c(res_main, res_alt), "density"))
freq_dens <- map_dfr(1:K, function(k) {
  b <- seq(b_hat[k] - 6 * SE[k], b_hat[k] + 6 * SE[k], length.out = 400)
  tibble(term = terms[k], beta = b, dens = dnorm(b, b_hat[k], SE[k]))
})
ggplot() +
  geom_line(data = freq_dens, aes(beta, dens), colour = "grey40",
            linetype = "dashed", linewidth = 0.9) +
  geom_line(data = dens_all, aes(beta, dens, colour = support), linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  coord_cartesian(xlim = range(freq_dens$beta)) +
  facet_wrap(~ term, scales = "free") +
  labs(x = "effect (beta)", y = "density", colour = "reference set",
       title = "Combined posterior over all reasonable priors",
       subtitle = "grey dashed = the data alone (frequentist)")


## =================================================================
##  ADD-ON: the same reference set, centred on 0 instead of b_hat
##  (the estimate is not used to build it; sigma is, and is n-free)
## =================================================================
D_HALF0 <- 10      # wider grid, so the 1/32 region also fits
R_MAX0  <- 25
N_D0    <- 121
N_R0    <- 121

## One-observation score, with the prediction evaluated at 'center'
log_score1_c <- function(k, m, s, center) {
  v  <- sort(c(m + s * sinh(seq(-12, 12, length.out = 801)),
               center + sigma1[k] * seq(-10, 10, length.out = 401)))
  lf <- lprior(v, m, s) + dnorm(center, v, sigma1[k], log = TRUE)
  max(lf) + log(trap(v, exp(lf - max(lf))))
}

grid0 <- map_dfr(1:K, function(k) {
  g <- expand.grid(d = seq(-D_HALF0, D_HALF0, length.out = N_D0),
                   r = seq(R_MIN, R_MAX0, length.out = N_R0))
  m <- g$d * sigma1[k]
  s <- g$r * sigma1[k]
  tibble(term = terms[k], d = g$d, r = g$r, m = m, s = s,
         log_support = map2_dbl(m, s, function(m, s) log_score1_c(k, m, s, 0)) -
                       log_best1(k))
}) |>
  mutate(keep_main = log_support >= -log(K_MAIN),
         keep_alt  = log_support >= -log(K_ALT))

grid0 |>
  group_by(term) |>
  summarise(hits_side = any(keep_alt & abs(d) == max(abs(d))),
            hits_top  = any(keep_alt & r == max(r)),
            n_main = sum(keep_main), n_alt = sum(keep_alt)) |>
  print()   # should be FALSE

## Full-data posterior and BF for every prior in the 0-centred set
kept0 <- filter(grid0, keep_alt)
int0  <- t(mapply(function(term, m, s) post_int(match(term, terms), m, s),
                  kept0$term, kept0$m, kept0$s))
kept0 <- kept0 |>
  mutate(pmean = int0[, "mean"], psd = int0[, "sd"], logBF = int0[, "logBF"])

## Equal-weight combination for any reference set
combine_set <- function(df, flag, label) {
  map_dfr(1:K, function(k) {
    ck <- filter(df, term == terms[k], .data[[flag]])
    md <- mixture_density(k, ck)
    tibble(term      = terms[k],
           reference = label,
           n_priors  = nrow(ck),
           mix_mean  = mean(ck$pmean),
           mix_sd    = sqrt(mean(ck$psd^2) + mean((ck$pmean - mean(ck$pmean))^2)),
           lo95      = mix_quant(md$beta, md$dens, 0.025),
           hi95      = mix_quant(md$beta, md$dens, 0.975),
           BF10_min  = exp(min(ck$logBF)),
           BF10_max  = exp(max(ck$logBF)),
           BF10_avg  = exp(lse(ck$logBF) - log(nrow(ck))))
  })
}

## The JZS default Bayes factor for each coefficient (BayesFactor package)
library(BayesFactor)
df_bf <- setNames(data.frame(y, X), c("y", terms))
jzs <- map_dfr(terms, function(tm) {
  f_full  <- as.formula(paste("y ~", paste(terms, collapse = " + ")))
  bf_full <- lmBF(f_full, data = df_bf, progress = FALSE)
  others  <- setdiff(terms, tm)
  bf <- if (length(others) == 0) {
    extractBF(bf_full)$bf
  } else {
    f_red <- as.formula(paste("y ~", paste(others, collapse = " + ")))
    extractBF(bf_full / lmBF(f_red, data = df_bf, progress = FALSE))$bf
  }
  tibble(term = tm, reference = "JZS default", BF10_avg = bf)
})

## Side-by-side comparison
comparison <- bind_rows(
  combine_set(final, "keep_main", "centred on b_hat, 1/8"),
  combine_set(final, "keep_alt",  "centred on b_hat, 1/32"),
  combine_set(kept0, "keep_main", "centred on 0, 1/8"),
  combine_set(kept0, "keep_alt",  "centred on 0, 1/32"),
  jzs) |>
  left_join(tibble(term = terms, b_hat = b_hat, SE = SE), by = "term") |>
  mutate(price    = mix_sd / SE,
         BF01_avg = 1 / BF10_avg,
         P_effect = BF10_avg / (1 + BF10_avg)) |>
  arrange(term, reference)

print(select(comparison, term, reference, BF10_avg, BF01_avg, P_effect,
             BF10_min, BF10_max))
print(select(comparison, term, reference, mix_mean, lo95, hi95, price))

## Evidence across the 0-centred reference set
ggplot(filter(kept0, keep_main), aes(m, s)) +
  geom_tile(aes(fill = logBF)) +
  geom_vline(data = tibble(term = terms, b_hat = b_hat),
             aes(xintercept = b_hat), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_gradient2(midpoint = 0) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~ term, scales = "free") +
  labs(x = "prior location (mu)", y = "prior scale (s)", fill = "log BF10",
       title = "Evidence across the 0-centred reference set (1/8)",
       subtitle = "dashed = estimate, dotted = 0")


comparison_sw <- map_dfr(list(list(final, "centred on b_hat"),
                              list(kept0, "centred on 0")), function(x) {
  x[[1]] |>
    filter(keep_main) |>
    group_by(term) |>
    summarise(reference = paste(x[[2]], "1/8, support-weighted"),
              BF10_avg  = exp(lse(logBF + log_support) - lse(log_support)),
              .groups = "drop")
}) |>
  mutate(BF01_avg = 1 / BF10_avg)
print(comparison_sw)


## =================================================================
##  DATA-BOUNDED REFERENCE SET
##  Prior predictive check on SD(y): no effect can explain more than
##  all of y's variance  ->  |beta_k| <= B_k  (partial-slope bound).
##  Candidate priors: Cauchy(mu, s) truncated to [-B_k, B_k].
##  Scale axis: tau = prior SD / SD of the uniform on [-B_k, B_k]
##  (0 = point prior, 1 = flat over every possible effect).
##  Equal weight over mu in [-B, B] and tau in (0, 1).
## =================================================================
N_MU_B  <- 81
N_TAU_B <- 41
SLACK   <- 1 + 2 / sqrt(2 * N)    # allowance for sampling error in the SDs

## ---- 1. The bound per coefficient, from the spread of y ---------
## Partial slope = partial correlation * SD(resid y) / SD(resid x),
## and |partial correlation| <= 1.
B_k <- sapply(1:K, function(k) {
  if (K == 1) return(sd(y) / sd(X[, 1]) * SLACK)
  ry <- resid(lm(y ~ X[, -k]))
  rx <- resid(lm(X[, k] ~ X[, -k]))
  sd(ry) / sd(rx) * SLACK
})
print(tibble(term = terms, bound = B_k, b_hat = b_hat,
             jzs_mass_beyond_bound = 1 - 2 / pi * atan(B_k / 0.354)))

## ---- 2. Truncated Cauchy prior --------------------------------------
lprior_tr <- function(b, m, s, B) {
  lz <- log(pt((B - m) / s, NU) - pt((-B - m) / s, NU))
  ifelse(abs(b) <= B, lprior(b, m, s) - lz, -Inf)
}

clip <- function(x, B) pmin(pmax(x, -B), B)

## SD of a truncated prior (for the tau axis)
tr_sd <- function(m, s, B) {
  b <- sort(unique(c(seq(-B, B, length.out = 2001),
                     clip(m + s * sinh(seq(-10, 10, length.out = 401)), B))))
  f <- exp(lprior_tr(b, m, s, B))
  Z <- trap(b, f); mu1 <- trap(b, b * f) / Z
  sqrt(trap(b, (b - mu1)^2 * f) / Z)
}

## ---- 3. Full-data posterior and BF under one truncated prior ------
post_int_tr <- function(k, m, s, B) {
  b  <- sort(unique(c(seq(-B, B, length.out = 2001),
                      clip(m + s * sinh(seq(-10, 10, length.out = 401)), B),
                      clip(b_hat[k] + SE[k] * seq(-10, 10, length.out = 401), B))))
  lmg <- lmarg(b, k)
  lp  <- lmg + lprior_tr(b, m, s, B)
  f   <- exp(lp - max(lp))
  Z   <- trap(b, f)
  pm  <- trap(b, b * f) / Z
  lr  <- lp - lprior(b, 0, S_WIDE)          # relative to the pilot prior
  c(mean  = pm,
    sd    = sqrt(trap(b, (b - pm)^2 * f) / Z),
    logBF = log(BF_pilot[k]) + max(lr) + log(trap(b, exp(lr - max(lr)))))
}

## ---- 4. The reference set: a square in (mu, tau) -------------------
taus  <- seq(0.01, 0.99, length.out = N_TAU_B)
s_tab <- 10^seq(-4, 4, length.out = 200)      # multiplied by B below

bounded <- map_dfr(1:K, function(k) {
  B <- B_k[k]
  map_dfr(seq(-B, B, length.out = N_MU_B), function(m) {
    # translate tau into the Cauchy scale s for this location
    tau_tab <- sapply(B * s_tab, function(s) tr_sd(m, s, B)) / (B / sqrt(3))
    s_vals  <- approx(tau_tab, B * s_tab, xout = taus, rule = 2, ties = mean)$y
    res     <- t(sapply(s_vals, function(s) post_int_tr(k, m, s, B)))
    tibble(term = terms[k], B = B, m = m, tau = taus, s = s_vals,
           pmean = res[, "mean"], psd = res[, "sd"], logBF = res[, "logBF"])
  })
})

## ---- 5. Combine with equal weight ----------------------------------
mix_tr <- function(k, ck) {
  B   <- B_k[k]
  b   <- seq(-B, B, length.out = 4001)
  lmg <- lmarg(b, k)
  acc <- numeric(length(b))
  for (i in seq_len(nrow(ck))) {
    lp  <- lmg + lprior_tr(b, ck$m[i], ck$s[i], B)
    f   <- exp(lp - max(lp))
    acc <- acc + f / trap(b, f)
  }
  tibble(beta = b, dens = acc / nrow(ck))
}

bounded_res <- map(1:K, function(k) {
  ck <- filter(bounded, term == terms[k])
  md <- mix_tr(k, ck)
  list(
    summary = tibble(
      term       = terms[k],
      reference  = "data-bounded (SD of y)",
      n_priors   = nrow(ck),
      mix_mean   = mean(ck$pmean),
      within_sd  = sqrt(mean(ck$psd^2)),
      between_sd = sqrt(mean((ck$pmean - mean(ck$pmean))^2)),
      lo95       = mix_quant(md$beta, md$dens, 0.025),
      hi95       = mix_quant(md$beta, md$dens, 0.975),
      BF10_min   = exp(min(ck$logBF)),
      BF10_max   = exp(max(ck$logBF)),
      BF10_avg   = exp(lse(ck$logBF) - log(nrow(ck)))),
    density = mutate(md, term = terms[k]))
})

bounded_summary <- bind_rows(map(bounded_res, "summary")) |>
  left_join(tibble(term = terms, b_hat = b_hat, SE = SE), by = "term") |>
  mutate(mix_sd    = sqrt(within_sd^2 + between_sd^2),
         price     = mix_sd / SE,
         freq_lo95 = b_hat - 1.96 * SE,
         freq_hi95 = b_hat + 1.96 * SE,
         BF01_avg  = 1 / BF10_avg,
         P_effect  = BF10_avg / (1 + BF10_avg))

print(select(bounded_summary, term, mix_mean, lo95, hi95, freq_lo95, freq_hi95,
             within_sd, between_sd, price))

## ---- 6. Evidence compared with the other constructions -------------
bind_rows(
  select(bounded_summary, term, reference, BF10_avg, BF10_min, BF10_max),
  filter(comparison, reference %in% c("centred on 0, 1/8", "JZS default")) |>
    select(term, reference, BF10_avg, BF10_min, BF10_max)) |>
  mutate(BF01_avg = 1 / BF10_avg,
         P_effect = BF10_avg / (1 + BF10_avg)) |>
  arrange(term, reference) |>
  print()

## =================================================================
##  PLOTS
## =================================================================
ref_b <- tibble(term = terms, b_hat = b_hat)

# A. Evidence across the reference square
ggplot(bounded, aes(m, tau)) +
  geom_tile(aes(fill = logBF)) +
  geom_vline(data = ref_b, aes(xintercept = b_hat), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_gradient2(midpoint = 0) +
  facet_wrap(~ term, scales = "free_x") +
  labs(x = "prior location (mu), within the possible range",
       y = "tau: 0 = point prior, 1 = flat over all possible effects",
       fill = "log BF10",
       title = "Evidence across all priors the data allow",
       subtitle = "dashed = estimate, dotted = 0")

# B. How far each prior pulls the posterior
ggplot(bounded, aes(m, tau)) +
  geom_tile(aes(fill = (pmean - b_hat[match(term, terms)]) / SE[match(term, terms)])) +
  geom_vline(data = ref_b, aes(xintercept = b_hat), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_gradient2(midpoint = 0) +
  facet_wrap(~ term, scales = "free_x") +
  labs(x = "prior location (mu)", y = "tau",
       fill = "posterior mean\n(SE from b_hat)",
       title = "How much each prior moves the posterior")

# C. Combined posterior vs the data alone
freq_d <- map_dfr(1:K, function(k) {
  b <- seq(b_hat[k] - 6 * SE[k], b_hat[k] + 6 * SE[k], length.out = 400)
  tibble(term = terms[k], beta = b, dens = dnorm(b, b_hat[k], SE[k]))
})
ggplot() +
  geom_line(data = freq_d, aes(beta, dens), colour = "grey40",
            linetype = "dashed", linewidth = 0.9) +
  geom_line(data = bind_rows(map(bounded_res, "density")),
            aes(beta, dens), colour = "steelblue", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  coord_cartesian(xlim = range(freq_d$beta)) +
  facet_wrap(~ term, scales = "free") +
  labs(x = "effect (beta)", y = "density",
       title = "Combined posterior over all data-allowed priors",
       subtitle = "grey dashed = the data alone (frequentist)")


flat_bf <- map_dfr(1:K, function(k) {
  B  <- B_k[k]
  b  <- seq(-B, B, length.out = 4001)
  lr <- lmarg(b, k) + log(1 / (2 * B)) - lprior(b, 0, S_WIDE)   # flat prior on [-B, B]
  tibble(term = terms[k],
         BF10_flat = BF_pilot[k] * exp(max(lr)) * trap(b, exp(lr - max(lr))))
}) |>
  left_join(select(bounded_summary, term, BF10_square = BF10_avg), by = "term") |>
  mutate(BF01_flat = 1 / BF10_flat, BF01_square = 1 / BF10_square)
print(flat_bf)

## Bayes factor for a t-test with a flat prior on the point-biserial
## correlation r (every possible group-outcome correlation equally plausible).
## Integrates over r in (-1, 1) directly: the prior there is simply 1/2.
bf_possible_t <- function(t, n1, n2, eps = 1e-6) {
  df   <- n1 + n2 - 2
  neff <- sqrt(n1 * n2 / (n1 + n2))
  a    <- (n1 / (n1 + n2)) * (n2 / (n1 + n2))       # p * q
  lik_r <- function(r) {
    d <- r / sqrt(a * (1 - r^2))                     # r -> Cohen's d
    v <- suppressWarnings(dt(t, df, ncp = d * neff))
    v[!is.finite(v)] <- 0                            # extreme r: no contribution
    v * 0.5                                          # flat prior on r
  }
  m1 <- integrate(lik_r, -1 + eps, 1 - eps,
                  subdivisions = 1000, rel.tol = 1e-8)$value
  m0 <- dt(t, df)
  c(BF10 = m1 / m0, BF01 = m0 / m1)
}

## Compare with the JASP / BayesFactor default (Cauchy, scale 0.707)
compare_t <- function(t, n1, n2) {
  library(BayesFactor)
  jzs <- exp(ttest.tstat(t, n1, n2, rscale = "medium", simple = FALSE)$bf)
  rbind(possible_range = bf_possible_t(t, n1, n2),
        jzs_default    = c(BF10 = jzs, BF01 = 1 / jzs))
}

compare_t(t = 2.5, n1 = 40, n2 = 40)
compare_t(t = 0.5, n1 = 40, n2 = 40)

hyp_chart <- function(t, n1, n2) {
  df   <- n1 + n2 - 2
  neff <- sqrt(n1 * n2 / (n1 + n2))
  d    <- seq(0, 1, by = 0.05)
  tibble::tibble(d  = d,
                 LR = dt(t, df, ncp = d * neff) / dt(t, df))   # vs d = 0
}
print(hyp_chart(0.5, 40, 40), n = 21)
