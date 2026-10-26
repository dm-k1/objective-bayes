# Relation to other priors

**In one sentence:** the flat-R prior is a mixture of Zellner g-priors, the
family of the JZS default, and the only member whose implied $R^2$ stays
the same for every number of predictors *and* is Jeffreys' flat correlation for
one predictor (constant in $k$ alone is not enough: $R^2 \sim$ Beta(½, $b$) for
every $k$ is a g-mixture for every $b$, checked in `validation/09`, V5; see
[`2_how_it_calculates.md`](2_how_it_calculates.md), section 8). Its closest neighbour is the other
natural way to extend Jeffreys' flat correlation to several predictors: flat
over all correlation matrices, which is rstanarm `stan_lm`'s $R^2$ prior with
η = 1.

## Its family: mixtures of g-priors

Zellner–Siow (JZS, the BayesFactor/JASP default), hyper-g, hyper-g/n,
beta-prime, robust and intrinsic priors all put a prior on Zellner's $g$. The
flat-R prior does too, with $g/n = B \cdot W$, $B \sim$ Beta(½, (k−1)/2),
$W \sim$ Inverse-Gamma(1, ½) ([`2_how_it_calculates.md`](2_how_it_calculates.md),
section 8). Median implied $R^2$ before any data ($n = 100$):

| Prior | $k = 1$ | $k = 3$ | $k = 10$ | $k = 30$ |
|---|---|---|---|---|
| Zellner–Siow (scale 1) | 0.50 | 0.84 | 0.95 | 0.98 |
| hyper-g/n ($a = 3$) | 0.54 | 0.87 | 0.96 | 0.99 |
| robust (Bayarri et al., 2012) | 0.52 | 0.72 | 0.77 | 0.78 |
| intrinsic | 0.39 | 0.59 | 0.65 | 0.66 |
| **flat-R** | 0.25 | 0.25 | 0.25 | 0.25 |

The flat-R hyperprior involves Kummer's confluent hypergeometric function, so
it is not a special case of the tCCH family of Li & Clyde (2018), which unifies
hyper-g, beta-prime, robust and intrinsic priors. A search of the g-prior
literature did not turn it up; that is a finding of a search, not a proof.

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
| rstanarm `stan_lm` | yes, by design | $R^2 \sim$ Beta($k/2$, η); η = 1 is "every correlation matrix equally likely" (next section) |
| hyper-g (Liang et al., 2008), Maruyama & George (2011) | partly | prior on the shrinkage $g/(1+g)$ |
| regularized horseshoe | no | global scale tied to the expected number of relevant predictors |
| R2-D2 (Zhang et al., 2022), R2D2M2 (Aguilar & Bürkner, 2023) | no | fixed prior on $R^2$, Dirichlet split |
| **flat-R** | no | flat on $R$ for every $k$, even split |

## The other extension: flat over all correlation matrices

With one predictor, flat-R is Jeffreys' flat correlation. There are two natural
ways to extend that to $k$ predictors, and they differ in one claim:

| | what is flat | implied $R^2$ | median $R^2$, $k$ = 1 / 3 / 10 / 30 |
|---|---|---|---|
| **flat-R** | the total $R$, whatever $k$ is | Beta(½, 1) | 0.25 / 0.25 / 0.25 / 0.25 |
| **uniform correlations** | every correlation matrix of $(y, x_1 \ldots x_k)$ | Beta($k/2$, 1) | 0.25 / 0.63 / 0.87 / 0.95 |

- **Where Beta($k/2$, 1) comes from.** Under the LKJ distribution with shape η
  (Lewandowski, Kurowicka & Joe, 2009; η = 1 is uniform), the squared multiple
  correlation of one variable on the others is Beta($k/2$, η), the classical
  null law of the sample $R^2$ (Muirhead, 1982, Theorem 5.2.2). How it was
  found and checked is in the next subsection.
- **It is rstanarm's `stan_lm` prior.** `stan_lm` uses exactly this
  construction: $\sqrt{R^2}$ times a uniformly random direction in orthogonalised
  predictor coordinates (its QR step), with $R^2 \sim$ Beta($k/2$, η) (Gabry &
  Goodrich, rstanarm vignette). So flat-R and `stan_lm` share everything except
  the prior on $R^2$. But `stan_lm` has no default η: the user states where
  $R^2$ is expected, and η follows; the rstanarm documentation does not
  recommend η = 1. Fixing η = 1, because it is the one value that gives Jeffreys'
  flat correlation at $k = 1$, and using it for Bayes factors, is not found there.
- **The same growth is forced by consistency between model sizes.** As a
  g-mixture (section 8 of [`2_how_it_calculates.md`](2_how_it_calculates.md)) it
  has $B \equiv 1$: a fixed hyperprior $g/n \sim$ Inverse-Gamma(1, ½), so the
  prior is the same scale mixture of normals for every $k$, and models of
  different sizes agree with each other. Flat-R gives that up to keep the total
  fixed.
- **Zellner–Siow is the η = ½ member.** With scale 1, $f^2 = \chi^2_k/\chi^2_1$, so
  $R^2 \sim$ Beta($k/2$, ½) exactly (the medians in the table above, 0.50 / 0.84 /
  0.95 / 0.98, are those of Beta($k/2$, ½)).
- **Keeping every single correlation flat leads to the same place.** Under
  LKJ(1) single correlations narrow as $k$ grows (median $|\rho|$ of $y$ with one
  predictor 0.51, 0.35, 0.23 for $k$ = 1, 3, 8). The marginally uniform prior of
  Barnard, McCulloch & Meng (2000), the correlations of an inverse-Wishart with
  $k + 2$ degrees of freedom, keeps every single correlation flat for any $k$, and
  implies the same $R^2 \sim$ Beta($k/2$, 1) (simulation, KS ≤ 0.007 for $k$ =
  2–8, 40,000 draws). So "flat on the correlations", in either sense, gives
  $a = k/2$; only "flat on the total" gives flat-R.

### How Beta($k/2$, 1) was found: simulate, then derive, then check

**In one sentence:** draw whole correlation matrices with no preference,
compute the $R^2$ each one implies, and compare the spread of those $R^2$
values with the Beta distribution that a short argument predicts.

1. **Make "the correlations, flat, as $k$ grows" well defined.** The $k$
   correlations of $y$ with the predictors can't be drawn one by one: together
   with the correlations among the predictors they must form a valid
   correlation matrix, or $R^2$ can exceed 1. So the whole matrix of
   $(y, x_1 \ldots x_k)$ is drawn at once, with every valid matrix equally
   likely (LKJ, η = 1). At $k = 1$ that is a flat correlation, the starting
   point of [`1_from_behaviour.md`](1_from_behaviour.md), step 1.
2. **Draw such matrices.** The LKJ algorithm does not draw correlations
   directly. It draws **partial correlations** (= the correlation left after
   the earlier variables are taken out). These can be chosen independently
   and always give a valid matrix. Each is a Beta on (−1, 1), and they are
   converted back into ordinary correlations: `draw_uniform_correlation_matrix()`
   in `validation/09_variable_selection.R`.
3. **Compute each matrix's population $R^2$.** With $\rho$ the correlations of
   $y$ with the predictors and $\Sigma$ the correlations among the predictors:

   ```r
   r_squared <- rho %*% solve(Sigma, rho)   # population R-squared of y on the x's
   ```

   Thousands of matrices at $k$ = 1, 2, 3, 5 and 8 give the spread of $R^2$:
   the $R^2$ prior that "every correlation matrix equally likely" implies.
4. **Predict the formula, then test it.** $R^2$ builds up one predictor at a time:

   $$
   1 - R^2 = (1 - p_1^2)(1 - p_2^2) \cdots (1 - p_k^2),
   $$

   with $p_j$ the partial correlations of step 2: independent Betas whose
   parameters step up by ½ each time. A product of Betas like that collapses
   into one Beta, which predicts $R^2 \sim$ Beta($k/2$, η). The simulated $R^2$
   matches that Beta to within simulation noise (KS ≈ 0.01; check V4 in 09),
   and η = 2 matches too, so the formula holds beyond the flat case.

Neither half is enough alone: the simulation without the formula is only
numbers, and the formula without the simulation could hide a slip. The literal
reading "every single correlation flat" was checked the same way, with the
correlations of an inverse-Wishart with $k + 2$ degrees of freedom (Barnard et
al., 2000): the same Beta($k/2$, 1).

This is the $R^2$ prior of the **uniform correlations** extension. The flat-R
share per predictor, Beta(½, (k−1)/2), is a different result found the same way:
[`1_from_behaviour.md`](1_from_behaviour.md), step 4.

**In variable selection** (`validation/README.md`, 09) neither dominates:
flat-R separates real from null predictors better in small samples with
uncorrelated predictors, uniform correlations does better when predictors are
correlated. Uniform correlations is at least as good as JZS everywhere; flat-R
beats JZS in small uncorrelated samples and trails it slightly (AUC ≤ 0.007)
with correlated predictors. Which extension describes
real research better is an empirical question: how fast explained variance
grows with the number of predictors in published regressions. Uniform
correlations predicts a median $R^2$ of 0.76 with 5 predictors and 0.84 with 8;
flat-R predicts 0.25 at any size.

## The conditional Lindley paradox

Priors with one scale shared by all coefficients (Zellner–Siow, hyper-g, and the
joint flat-R prior) let large coefficients drive the test of a small one toward
"no effect" (Som, Hans & MacEachern, 2016). Their remedy, block hyper-g priors,
gives groups of coefficients their own scale. The conditional flat-R test does
the same for one coefficient: its own share of the flat total
([`2_how_it_calculates.md`](2_how_it_calculates.md), section 9).

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
  same kind of prior on the *total* instead of on each partial correlation. The
  conditional flat-R test uses the same inputs (the sample partial correlation,
  $n$ and $k$), with the prior the flat total implies for one coefficient.
- rstanarm `stan_lm` (Gabry & Goodrich): the same construction with
  $R^2 \sim$ Beta($k/2$, η), for estimation only; with η = 1 it is the
  "uniform correlations" extension above.

## References

- Li, Y., & Clyde, M. A. (2018). Mixtures of g-priors in generalized linear
  models. *JASA*.
- Som, A., Hans, C. M., & MacEachern, S. N. (2016). A conditional Lindley paradox
  in Bayesian linear models. *Biometrika*.
- Wetzels, R., & Wagenmakers, E.-J. (2012). A default Bayesian hypothesis test
  for correlations and partial correlations. *Psychonomic Bulletin & Review*.
- Kucharský, Š., Wagenmakers, E.-J., van den Bergh, D., & Ly, A. (2026).
  Analytic posterior distribution and Bayes factor for Pearson partial
  correlations. *TEST*.

- Aguilar, J. E., & Bürkner, P.-C. (2023). Intuitive joint priors for Bayesian
  linear multilevel models: The R2D2M2 prior. *Electronic Journal of Statistics*.
- Barnard, J., McCulloch, R., & Meng, X.-L. (2000). Modeling covariance
  matrices in terms of standard deviations and correlations, with application
  to shrinkage. *Statistica Sinica*.
- Gabry, J., & Goodrich, B. Estimating regularized linear models with
  rstanarm (`stan_lm` vignette). https://mc-stan.org/rstanarm/articles/lm.html
- Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
  correlation matrices based on vines and extended onion method. *Journal of
  Multivariate Analysis*.
- Muirhead, R. J. (1982). *Aspects of Multivariate Statistical Theory*. Wiley.
- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012). Criteria
  for Bayesian model choice with application to variable selection. *Annals of Statistics*.
- Liang, F., Paulo, R., Molina, G., Clyde, M. A., & Berger, J. O. (2008).
  Mixtures of g priors for Bayesian variable selection. *JASA*.
- Maruyama, Y., & George, E. I. (2011). Fully Bayes factors with a generalized
  g-prior. *Annals of Statistics*.
- Zhang, Y. D., Naughton, B. P., Bondell, H. D., & Reich, B. J. (2022).
  Bayesian regression using a prior on the model fit: The R2-D2 shrinkage prior.
  *JASA*.
