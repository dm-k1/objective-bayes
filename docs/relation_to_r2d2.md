# Is this just R2-D2? Relation to the R2-D2 prior

Short answer: **no.** Both priors start from a Beta(½, 1) prior on the explained
share of variance (R²), but they turn it into coefficients differently. They
are not the same prior, not even with one predictor.

## Two differences

### 1. How R² becomes coefficients

| | Possible-range prior (this project) | R2-D2 (brms `R2D2()`, Zhang et al. 2022) |
|---|---|---|
| Given a drawn R² | the size of the effect is **fixed** by it; only the direction is random | R² sets the **variance** of a normal distribution for the coefficients |
| Consequence | the R² implied by the coefficients follows Beta(½, 1) exactly | the implied R² does **not** follow Beta(½, 1): the normal draw can land well below or above it |

### 2. How the total is split over predictors

| | Possible-range prior | R2-D2 |
|---|---|---|
| Split | uniform direction, measured after whitening the predictors (accounts for their correlations) | Dirichlet split per predictor, scaled by each predictor's own SD |
| Correlated predictors | handled in the prior | not taken into account in the split |

## What it means: one predictor

With one predictor there is nothing to split, so only difference 1 matters.
Implied prior on the correlation r (1,000,000 draws each):

| | Possible-range (flat r) | brms R2D2 with Beta(½, 1) |
|---|---|---|
| P(\|r\| < 0.1) | 0.10 | **0.25** |
| P(\|r\| < 0.5) | 0.50 | **0.65** |
| P(\|r\| > 0.9) | 0.10 | 0.08 |
| median \|r\| | 0.50 | **0.31** |

The possible-range prior is exactly flat on the correlation ("every value
equally plausible"). R2-D2 with "the same" R² prior puts 2.5 times as much
weight on near-zero correlations: it leans toward small effects.

With several predictors the two diverge further, because of difference 2.

## What it means in practice

In the joint simulation (`validation/03_sim_joint_prior.R`, 4 predictors),
R2-D2 with Beta(½, 1) and Dirichlet(½):

- found true nulls slowest: P(BF₀₁ > 3) for the null predictor at n = 200 was
  0.16, against 0.41 for the possible-range prior and 0.61 for JZS;
- almost never called a small real effect null;
- had slightly more power with correlated predictors.

This fits difference 1: R2-D2's alternative already leans toward small
effects, so data must look clearly different from small effects before it
favours "no effect".

## Why the custom brms prior is needed

brms has no built-in way to express "flat on the correlation" or a uniform
direction after whitening. The closest built-in option,

```r
prior(R2D2(mean_R2 = 1/3, prec_R2 = 1.5), class = b)   # Beta(1/2, 1) on R2
```

is a different prior (see above). That is why `fit_possible_range_model()` adds the
possible-range prior as a custom density (`docs/derivation.md`, section 3).

## Reproduce the table

```r
set.seed(1)
number_of_draws <- 1e6

## Possible-range prior, one predictor: flat on the correlation
correlation_possible <- runif(number_of_draws, -1, 1)

## brms-style R2D2, one predictor: R2 ~ Beta(1/2, 1) sets the variance
explained_share  <- rbeta(number_of_draws, 0.5, 1)
signal_to_noise  <- explained_share / (1 - explained_share)
effect_in_noise_units <- rnorm(number_of_draws, 0, sqrt(signal_to_noise))
correlation_r2d2 <- effect_in_noise_units / sqrt(1 + effect_in_noise_units^2)

summarise_correlation <- function(r) c(
  abs_below_0.1 = mean(abs(r) < 0.1),
  abs_below_0.5 = mean(abs(r) < 0.5),
  abs_above_0.9 = mean(abs(r) > 0.9),
  median_abs    = median(abs(r)))

round(rbind(possible_range = summarise_correlation(correlation_possible),
            r2d2           = summarise_correlation(correlation_r2d2)), 3)
```

## References

- Zhang, Y. D., Naughton, B. P., Bondell, H. D., & Reich, B. J. (2022).
  Bayesian regression using a prior on the model fit: The R2-D2 shrinkage prior.
  *Journal of the American Statistical Association*.
- Aguilar, J. E., & Bürkner, P.-C. (2023). Intuitive joint priors for Bayesian
  linear multilevel models: The R2D2M2 prior. *Electronic Journal of Statistics*.
