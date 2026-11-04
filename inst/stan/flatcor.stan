// One of the three priors, R^2 ~ Beta(a, 1) with a uniform direction, sampled
// in its own coordinates:
//   R          the overall correlation (factor 1): R^2 ~ Beta(a, 1), with
//              a = K/2 (flat_partial_r: every correlation matrix of
//              (y, x1 ... xK) equally likely), 1/2 (flat_r) or 1 (flat_r_squared)
//   direction  how the total is divided, uniform on the sphere (factor 2)
// and the coefficients follow from them (factor 3 is the conversion):
//   b = sigma * f * (direction expressed in predictor units),  f = R / sqrt(1 - R^2)
// In these coordinates the prior has no peak at zero, so the sampler
// reaches small effects for any number of predictors.
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;                 // predictors (centered here)
  vector[N] y;
  int<lower=0, upper=1> prior_only;      // 1: ignore y (sample the prior)
  int<lower=0, upper=1> proper_stand_ins; // 1: proper priors on alpha and sigma (for checks)
  real<lower=0> r_squared_shape;         // a in R^2 ~ Beta(a, 1)
}
transformed data {
  matrix[N, K] Xc;
  for (j in 1:K) Xc[, j] = X[, j] - mean(X[, j]);
  // Any square root of the predictor covariance works; Cholesky is the cheapest.
  // It turns a direction in "explained variance" units into predictor units.
  matrix[K, K] root = cholesky_decompose(crossprod(Xc) / N);
}
parameters {
  real alpha;
  real<lower=0> sigma;
  // with one predictor the sign lives in R; with more, in the direction
  real<lower=(K == 1 ? -1 : 0), upper=1> R;
  // the direction is z / |z|. Any density for z that depends only on its length
  // makes the direction uniform; this one keeps |z| near 1, away from 0, where
  // the direction would be undefined (Stan's unit_vector lets |z| reach 0,
  // which gives divergent transitions)
  vector[K > 1 ? K : 0] z;
}
transformed parameters {
  real f = R / sqrt(1 - square(R));          // Cohen's f (signed when K = 1)
  vector[K] direction;
  if (K == 1) direction = rep_vector(1, 1);
  else        direction = z / norm2(z);
  vector[K] b = sigma * f * mdivide_right_tri_low(direction', root)';   // b' Sigma b = sigma^2 f^2
}
model {
  // factor 1: R^2 ~ Beta(a, 1), which is density 2a |R|^(2a - 1) on R (flat
  // for a = 1/2; with K = 1, flat_partial_r and flat_r both have a = 1/2)
  if (r_squared_shape != 0.5) target += (2 * r_squared_shape - 1) * log(abs(R));
  // factor 2: direction uniform
  if (K > 1) target += -0.5 * square((norm2(z) - 1) / 0.1);
  if (proper_stand_ins) {
    // Checks that sample the prior, or simulate data from it, need proper priors
    // on alpha and sigma. R, f and the direction don't depend on these.
    alpha ~ normal(0, 10);
    sigma ~ lognormal(0, 1);
  } else {
    target += -log(sigma);                   // Jeffreys prior on the noise; alpha flat
  }
  if (!prior_only)
    y ~ normal(alpha + Xc * b, sigma);
}
generated quantities {
  // per observation: log-likelihood (for loo) and a predicted value (for pp_check)
  vector[N] log_lik;
  array[N] real y_rep;
  real R_squared = square(R);               // population explained share
  {
    vector[N] mu = alpha + Xc * b;
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
    y_rep = normal_rng(mu, sigma);
  }
}
