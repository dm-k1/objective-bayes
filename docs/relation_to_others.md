# Relation to other priors

**In one sentence:** most default priors put a prior on each coefficient, so the
explained variance they claim grows with every predictor; the priors built to
avoid that put a prior on the total, and the flat-R prior is one of those, with
the total flat on $R$ and an even split.

## Which priors drift with the number of predictors

A prior drifts when the $R^2$ it implies before seeing data rises as
predictors are added ([`1_from_behaviour.md`](1_from_behaviour.md), step 2).

| Prior | Drifts? | Why |
|---|---|---|
| independent normal / Student-t per coefficient (usual brms/Stan setup; rstanarm's autoscaled normal) | yes | each coefficient adds its own expected share |
| ridge with a fixed penalty | yes | = independent normals |
| uniform box per standardized coefficient | yes | mass sits in the corners |
| uniform over all coefficient vectors with $R^2 \le 1$ | yes | implies $R^2 \sim$ Beta($k/2$, 1) |
| Zellner's g-prior with fixed $g$ (e.g. $g = n$) | yes | $R^2$ rises with $k$ for fixed $g$ |
| JZS / Zellner–Siow (BayesFactor, JASP) | yes, milder | heavy-tailed mixture over $g$, but still per-dimension |
| spike-and-slab with a fixed slab | yes | with the number of included predictors |
| independent flat partial correlations, combined | yes, strongly | $1 - R^2$ is a product of one factor per predictor |
| rstanarm `stan_lm` | yes, by design | $R^2 \sim$ Beta($k/2$, η) |
| hyper-g (Liang et al., 2008), Maruyama & George (2011) | partly | prior on the shrinkage $g/(1+g)$ |
| regularized horseshoe | no | global scale tied to the expected number of relevant predictors |
| R2-D2 (Zhang et al., 2022), R2D2M2 (Aguilar & Bürkner, 2023) | no | fixed prior on $R^2$, Dirichlet split |
| **flat-R** | no | flat on $R$ for every $k$, even split |

## Is this just R2-D2?

No. Both can start from Beta(½, 1) on $R^2$, but they turn it into
coefficients differently, so they are different priors, even with one predictor.

| | Flat-R | R2-D2 (brms `R2D2()`) |
|---|---|---|
| Given a drawn $R^2$ | the length of the effect is **fixed** by it; only the direction is random | $R^2$ sets the **variance** of a normal distribution for the coefficients |
| Implied $R^2$ | follows Beta(½, 1) exactly | does **not** follow Beta(½, 1): the normal draw lands below or above it |
| Split over predictors | uniform direction in explained variance (accounts for correlated predictors) | Dirichlet split per predictor, scaled by each predictor's own SD |
| Correlated predictors | handled in the prior | not taken into account in the split |

**One predictor** (only the first difference matters). Implied prior on the
correlation, 1,000,000 draws each:

| | Flat-R | brms R2D2 with Beta(½, 1) |
|---|---|---|
| P(\|r\| < 0.1) | 0.10 | **0.25** |
| P(\|r\| < 0.5) | 0.50 | **0.65** |
| P(\|r\| > 0.9) | 0.10 | 0.08 |
| median \|r\| | 0.50 | **0.31** |

R2-D2 with "the same" $R^2$ prior puts 2.5 times as much weight on near-zero
correlations: it leans toward small effects. So it finds true nulls more
slowly and rarely calls a small real effect null (`validation/03_simulations.R`).

brms's closest built-in option,

```r
prior(R2D2(mean_R2 = 1/3, prec_R2 = 1.5), class = b)   # Beta(1/2, 1) on R2
```

is the different prior above. A custom brms prior can state the flat-R density,
but brms samples the coefficients directly and misses the prior's peak at zero
for larger $k$. That is why `fit_flat_r_model()` uses its own Stan model,
sampled in $R$ and a direction ([`2_how_it_calculates.md`](2_how_it_calculates.md), section 6).

Reproduce the one-predictor table:

```r
set.seed(1)
number_of_draws <- 1e6

## Flat-R, one predictor: flat on the correlation
correlation_flat_r <- runif(number_of_draws, min = -1, max = 1)

## brms-style R2D2, one predictor: R2 ~ Beta(1/2, 1) sets the variance
r_squared_r2d2   <- rbeta(number_of_draws, shape1 = 0.5, shape2 = 1)
f_squared_r2d2   <- r_squared_r2d2 / (1 - r_squared_r2d2)          # Cohen's f^2
effect_r2d2      <- rnorm(number_of_draws, mean = 0, sd = sqrt(f_squared_r2d2))
correlation_r2d2 <- effect_r2d2 / sqrt(1 + effect_r2d2^2)

summarise_correlation <- function(r) c(
  abs_below_0.1 = mean(abs(r) < 0.1),
  abs_below_0.5 = mean(abs(r) < 0.5),
  abs_above_0.9 = mean(abs(r) > 0.9),
  median_abs    = median(abs(r)))

round(rbind(flat_r = summarise_correlation(correlation_flat_r),
            r2d2   = summarise_correlation(correlation_r2d2)), 3)
```

## Closest relatives

- Jeffreys' flat prior on a correlation; Ly, Verhagen & Wagenmakers on
  Jeffreys' default tests: flat-R is this prior with one predictor.
- Kucharský, Wagenmakers, van den Bergh & Ly (2026), *Test*: exact Bayes factor
  for one partial correlation with a stretched-beta prior. Flat-R puts the
  same kind of prior on the *total* instead of on each partial correlation.

## References

- Aguilar, J. E., & Bürkner, P.-C. (2023). Intuitive joint priors for Bayesian
  linear multilevel models: The R2D2M2 prior. *Electronic Journal of Statistics*.
- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012). Criteria
  for Bayesian model choice with application to variable selection. *Annals of Statistics*.
- Liang, F., Paulo, R., Molina, G., Clyde, M. A., & Berger, J. O. (2008).
  Mixtures of g priors for Bayesian variable selection. *JASA*.
- Maruyama, Y., & George, E. I. (2011). Fully Bayes factors with a generalized
  g-prior. *Annals of Statistics*.
- Zhang, Y. D., Naughton, B. P., Bondell, H. D., & Reich, B. J. (2022).
  Bayesian regression using a prior on the model fit: The R2-D2 shrinkage prior.
  *JASA*.
