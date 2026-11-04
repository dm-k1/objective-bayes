# Relation to other priors

**In one sentence:** the three priors, $R^2 \sim$ Beta($a$, 1) with a uniform
direction, sit in the family of the JZS default (mixtures of Zellner g-priors)
and next to rstanarm `stan_lm`: flat partial r ($a = k/2$) *is* `stan_lm`'s
$R^2$ prior with η = 1, flat R ($a = ½$) is the only g-mixture whose implied
$R^2$ stays the same for every $k$ *and* is Jeffreys' flat correlation for one
predictor, and flat R-squared ($a = 1$) keeps the total flat on the $R^2$
ruler.

Constant in $k$ alone does not single out flat R: $R^2 \sim$ Beta(½, $b$) for
every $k$ is a g-mixture for every $b$ (checked in `validation/09`, V5; see
[`2_how_it_calculates.md`](2_how_it_calculates.md), section 8). Where the three
come from: [`1_from_behaviour.md`](1_from_behaviour.md), step 0 (study 00).

## Its family: mixtures of g-priors

Zellner–Siow (JZS, the BayesFactor/JASP default), hyper-g, hyper-g/n,
beta-prime, robust and intrinsic priors all put a prior on Zellner's $g$. The
three priors do too, with $g/n = B \cdot W$, $W \sim$ Inverse-Gamma(1, ½), and
$B \sim$ Beta($a$, $k/2 - a$): Beta(½, (k−1)/2) for flat R, Beta(1, $k/2 - 1$)
for flat R-squared ($k \ge 3$; at $k = 2$ it is flat partial r, at $k = 1$ it
is not a g-mixture), $B \equiv 1$ for flat partial r
([`2_how_it_calculates.md`](2_how_it_calculates.md), section 8). How much $R^2$
each member claims before any data, as $k$ grows, is the table in that section:
from 0.25 at every $k$ (flat R) and 0.5 at every $k$ (flat R-squared) to 0.98
at $k = 30$ (Zellner–Siow).

The flat R hyperprior involves Kummer's confluent hypergeometric function, so
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
| rstanarm `stan_lm`, and **flat partial r** (η = 1) | yes, by design | $R^2 \sim$ Beta($k/2$, η); η = 1 is "every correlation matrix equally likely" (next section) |
| hyper-g (Liang et al., 2008), Maruyama & George (2011) | partly | prior on the shrinkage $g/(1+g)$ |
| regularized horseshoe | no | global scale tied to the expected number of relevant predictors |
| R2-D2 (Zhang et al., 2022), R2D2M2 (Aguilar & Bürkner, 2023) | no | fixed prior on $R^2$, Dirichlet split |
| **flat R** | no | flat on $R$ for every $k$, even split |
| **flat R-squared** | no | flat on $R^2$ for every $k$, even split |

## Three extensions: flat on the correlations, on R, or on R²

With one predictor, Jeffreys' flat correlation is the starting point. Valid
correlation matrices allow three ways to extend it to $k$ predictors (study 00),
and they differ in one claim:

| | what is flat | implied $R^2$ | median $R^2$, $k$ = 1 / 3 / 10 / 30 | Jeffreys at $k = 1$? |
|---|---|---|---|---|
| **flat partial r** | every correlation matrix of $(y, x_1 \ldots x_k)$ | Beta($k/2$, 1) | 0.25 / 0.63 / 0.87 / 0.95 | yes |
| **flat R** | the total $R$, whatever $k$ is | Beta(½, 1) | 0.25 / 0.25 / 0.25 / 0.25 | yes |
| **flat R-squared** | the total $R^2$, whatever $k$ is | Beta(1, 1) | 0.50 / 0.50 / 0.50 / 0.50 | no: $r$ has density $\lvert r\rvert$ |

- **Where Beta($k/2$, 1) comes from.** Under the LKJ distribution with shape η
  (Lewandowski, Kurowicka & Joe, 2009; η = 1 is uniform), the squared multiple
  correlation of one variable on the others is Beta($k/2$, η), the classical
  null law of the sample $R^2$ (Muirhead, 1982, Theorem 5.2.2). How it was
  found and checked is in the next subsection.
- **It is rstanarm's `stan_lm` prior.** `stan_lm` uses exactly this
  construction: $\sqrt{R^2}$ times a uniformly random direction in orthogonalised
  predictor coordinates (its QR step), with $R^2 \sim$ Beta($k/2$, η) (Gabry &
  Goodrich, rstanarm vignette). So flat R and `stan_lm` share everything except
  the prior on $R^2$. But `stan_lm` has no default η: the user states where
  $R^2$ is expected, and η follows; the rstanarm documentation does not
  recommend η = 1. Why η = 1 in particular: the subsection "Why η = 1" below.
  Using it for Bayes factors is not found there. In this package:
  `prior = "flat_partial_r"` in every function.
- **The same growth is forced by consistency between model sizes.** As a
  g-mixture (section 8 of [`2_how_it_calculates.md`](2_how_it_calculates.md)) it
  has $B \equiv 1$: a fixed hyperprior $g/n \sim$ Inverse-Gamma(1, ½), so the
  prior is the same scale mixture of normals for every $k$, and models of
  different sizes agree with each other. Flat R gives that up to keep the total
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
  $a = k/2$; only "flat on the total" gives flat R.

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
   converted back into ordinary correlations: `draw_lkj_correlation_matrix()`
   in `validation/comparison_methods.R` (also used by study 00).
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

This is the $R^2$ prior of the **flat partial r** extension. The share per
predictor, Beta(½, (k−1)/2) for every member of the family, is a different
result found the same way: [`1_from_behaviour.md`](1_from_behaviour.md), step 4.

### Why η = 1

**In one sentence:** η = 1 is the only shape for which **every coefficient's
partial correlation is Jeffreys' flat correlation, whatever the number of other
predictors**.

| | What it says | Unique to η = 1? | Shown by |
|---|---|---|---|
| 1 | Under LKJ(η), the partial correlation of $y$ with one predictor given all the others is Beta(η, η) on (−1, 1), **for every $k$**. Flat only at η = 1 | yes | `validation/05`, C7 |
| 2 | At $k = 1$ the prior is Jeffreys' flat correlation | yes (LKJ(η) at $k = 1$ is Beta(η, η) on $r$) | the same, $k = 1$ |
| 3 | As a g-mixture, one hyperprior for every $k$: $g/n \sim$ Inverse-Gamma(1, ½), so models of every size share one scale | no (every η gives a fixed hyperprior) | doc 2, section 8 |
| 4 | "Every single correlation flat" (Barnard et al., 2000) leads to the same $R^2 \sim$ Beta($k/2$, 1) | — | simulation above |

Point 1 is why it matters for testing one coefficient: the test of a coefficient
is Jeffreys' correlation test on its partial correlation, with or without
covariates (in its fixed-X form: see "What is not new" below), and it follows from one prior on the whole model, not from a
separate rule per test. It holds for correlated predictors too, because the
partial correlation already takes the other predictors out.

Why point 1 holds: the LKJ sampler builds the matrix from independent partial
correlations, level by level, and the partials at the last level (two variables,
given all the others) are Beta(η, η) on (−1, 1), whatever the number of
variables. LKJ gives every ordering of the variables the same distribution, so
put $y$ and the tested predictor last: their partial correlation given all the
others is that last-level Beta(η, η).

Simulated (4,000 matrices per cell, KS distance):

| | $k = 1$ | $k = 3$ | $k = 8$ |
|---|---|---|---|
| η = 1, to flat | 0.017 | 0.010 | 0.011 |
| η = 2, to flat | 0.107 | 0.099 | 0.110 |
| η = 2, to Beta(2, 2) | 0.014 | 0.012 | 0.027 |

**What is not new.** A flat prior on each partial correlation is the default
partial-correlation test (Wetzels & Wagenmakers, 2012; Kucharský et al., 2026,
with a flat prior, α = 1). Ours is not the same Bayes factor, but it differs
from the exact one of Kucharský et al. in one assumption only (`validation/05`,
C8). Theirs is for multivariate normal data, so the predictors are random and
the tested predictor's residual SD is an unknown parameter; ours conditions on
X, as regression does, so the effect is in units of the observed residual SD.
Give ours that extra uncertainty and the two agree to 5e-7; the priors on the
noise and the other coefficients are the same in both. As it is, ours is lower
on the log scale by up to 0.26 (n = 20, k = 8, r = 0), 0.10 at n = 50 and
0.005 at n = 1000, largest at r = 0 (where ours does not depend on k).
What flat partial r adds is: the fixed-X version of their test, that this
test follows from one prior on the whole model, and that the same prior gives
Bayes factors between whole models and estimates (`stan_lm` with η = 1).

**In variable selection** ([`5_validation.md`](5_validation.md), 09) the three are compared
with JZS, Zellner–Siow and BIC. They rank real above null predictors almost
equally (flat partial r and flat R-squared within 0.025 of flat R in AUC in
every cell; +0.003 and −0.002 on average). They differ in where they draw the line: at
n = 500, uncorrelated, beta-binomial model prior, flat partial r includes 0.13
of 8 null predictors and misses 0.84 of 4 real ones, close to Zellner–Siow
(0.12, 0.89) and BIC (0.05, 1.01); flat R (0.91, 0.45) and flat R-squared
(0.92, 0.44) include more and miss fewer, more than JZS (0.42, 0.58). Which extension describes real
research better is partly an empirical question: how fast explained variance
grows with the number of predictors in published regressions. Flat partial r
predicts a median $R^2$ of 0.76 with 5 predictors and 0.84 with 8; flat R
predicts 0.25 and flat R-squared 0.5 at any size.

### Flat R-squared, and the earlier "growing" experiment

An earlier version of study 09 tried $R^2 \sim$ Beta($a$, 1) with
$a = k/(k+1)$: flat R at one predictor, moving toward flat R-squared as $k$
grows. It barely differed from flat R and was dropped. Flat R-squared is now a
full member of the comparison instead ($a = 1$ for every $k$), built and
checked like the other two, so that experiment is superseded. Its one
structural difference: at $k = 1$ it is flat on $r^2$, not on $r$, so it does
not start from Jeffreys' correlation (study 00, version C).

## The conditional Lindley paradox

Priors with one scale shared by all coefficients (Zellner–Siow, hyper-g, and the
joint tests of all three priors) let large coefficients drive the test of a
small one toward "no effect" (Som, Hans & MacEachern, 2016). Their remedy,
block hyper-g priors, gives groups of coefficients their own scale. The
conditional tests do the same for one coefficient: under flat R and flat
R-squared its own share of the total, under flat partial r its own flat partial
correlation
([`2_how_it_calculates.md`](2_how_it_calculates.md), section 10).

## Is this just R2-D2?

No. Both can start from Beta(½, 1) on $R^2$ (flat R), but they turn it into
coefficients differently, so they are different priors, even with one
predictor. The same holds for flat R-squared against R2-D2 with Beta(1, 1).

| | Flat R | R2-D2 (brms `R2D2()`) |
|---|---|---|
| Given a drawn $R^2$ | the length of the effect is **fixed** by it; only the direction is random | $R^2$ sets the **variance** of a normal distribution for the coefficients |
| Implied $R^2$ | follows Beta(½, 1) exactly | does **not** follow Beta(½, 1): the normal draw lands below or above it |
| Split over predictors | uniform direction in explained variance (accounts for correlated predictors) | Dirichlet split per predictor, scaled by each predictor's own SD |
| Correlated predictors | handled in the prior | not taken into account in the split |

**One predictor** (only the first difference matters). Implied prior on the
correlation, 1,000,000 draws each:

| | Flat R | brms R2D2 with Beta(½, 1) |
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

is the different prior above. A custom brms prior can state the flat R density,
but brms samples the coefficients directly and misses the prior's peak at zero
for larger $k$. That is why `fit_regression()` uses its own Stan model,
sampled in $R$ and a direction ([`2_how_it_calculates.md`](2_how_it_calculates.md), section 6).

Reproduce the one-predictor table:

```r
set.seed(1)
number_of_draws <- 1e6

## Flat R, one predictor: flat on the correlation
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
  Jeffreys' default tests: flat R is this prior with one predictor.
- Kucharský, Wagenmakers, van den Bergh & Ly (2026), *Test*: exact Bayes factor
  for one partial correlation with a stretched-beta prior. With α = 1 it is
  the flat partial r conditional test with random instead of fixed
  predictors (`validation/05`, C8). Flat R puts the
  same kind of prior on the *total* instead of on each partial correlation. The
  conditional flat R test uses the same inputs (the sample partial correlation,
  $n$ and $k$), with the prior the flat total implies for one coefficient.
- rstanarm `stan_lm` (Gabry & Goodrich): the same construction with
  $R^2 \sim$ Beta($k/2$, η), for estimation only; with η = 1 it is the
  "flat partial r" extension above.

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
