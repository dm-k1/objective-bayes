## ARCHIVED: direct cmdstanr version of the flat-R model, using a
## unit-vector parameterization (explained share + direction). Superseded
## by fit_flat_r_model() in R/flat_r.R. Kept as a fallback: if brms
## shows persistent divergent transitions (the prior has a sharp peak at
## all-zero coefficients), this parameterization avoids that peak.

## ---- C. Stan model for estimation on a single dataset -------------
## (K >= 2; for one predictor use the flat-r method)
library(cmdstanr)
joint_code <- "
data {
  int<lower=1> N;
  int<lower=2> K;
  matrix[N, K] X;
  vector[N] y;
}
transformed data {
  matrix[N, K] Xc;
  for (j in 1:K) Xc[, j] = X[, j] - mean(X[, j]);
  matrix[K, K] L      = cholesky_decompose(crossprod(Xc) / N);
  matrix[K, K] Linv_t = mdivide_left_tri_low(L, diag_matrix(rep_vector(1, K)))';
}
parameters {
  real alpha;
  real<lower=0> sigma;
  real<lower=0, upper=1> R2;      // population R-squared
  unit_vector[K] u;               // direction, uniform on the sphere
}
transformed parameters {
  vector[K] beta = sigma * sqrt(R2 / (1 - R2)) * (Linv_t * u);
}
model {
  target += -log(sigma);          // Jeffreys prior on sigma
  R2 ~ beta(0.5, 1);              // = flat on r when K = 1
  y ~ normal(alpha + Xc * beta, sigma);
}
"
joint_mod <- cmdstan_model(write_stan_file(joint_code))

## Example on one simulated dataset
set.seed(7)
Sig <- matrix(0.5, 4, 4); diag(Sig) <- 1
X   <- matrix(rnorm(100 * 4), 100) %*% chol(Sig)
y   <- drop(X %*% c(0.3, 0.3, 0, 0)) + rnorm(100, sd = 0.75)
fit_j <- joint_mod$sample(data = list(N = 100, K = 4, X = X, y = y),
                          chains = 4, parallel_chains = 4,
                          adapt_delta = 0.95, refresh = 0)
fit_j$summary(c("beta", "R2", "sigma"))
