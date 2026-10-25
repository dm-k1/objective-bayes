// The flat-R prior, sampled in its own coordinates:
//   R          the overall correlation, flat (factor 1)
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
}
transformed data {
  matrix[N, K] Xc;
  for (j in 1:K) Xc[, j] = X[, j] - mean(X[, j]);
  // Any square root of the predictor covariance works; Cholesky is the cheapest.
  // It turns a direction in "explained variance" units into predictor units.
  matrix[K, K] root = cholesky_decompose(crossprod(Xc) / N);
  int K_sphere = K > 1 ? K : 2;   // a unit vector needs at least 2 dimensions
}
parameters {
  real alpha;
  real<lower=0> sigma;
  // with one predictor the sign lives in R; with more, in the direction
  real<lower=(K == 1 ? -1 : 0), upper=1> R;
  unit_vector[K_sphere] direction_raw;
}
transformed parameters {
  real f = R / sqrt(1 - square(R));          // Cohen's f (signed when K = 1)
  vector[K] direction;
  if (K == 1) direction = rep_vector(1, 1);
  else        direction = direction_raw;
  vector[K] b = sigma * f * mdivide_right_tri_low(direction', root)';   // b' Sigma b = sigma^2 f^2
}
model {
  // factor 1: R flat (the declared bounds); factor 2: direction uniform (unit_vector)
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
