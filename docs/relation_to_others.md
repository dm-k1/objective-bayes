# Relation to other priors

The three priors in this project all put $R^2 \sim$ Beta($a$, 1) on the
explained variance and spread it over the predictors in a uniformly random
direction. They belong to the same family as the JZS default (mixtures of
Zellner *g*-priors) and sit next to rstanarm's `stan_lm`:

- flat partial r ($a = k/2$) is `stan_lm`'s $R^2$ prior with η = 1;
- flat R ($a = ½$) is the only *g*-mixture whose implied $R^2$ is the same for
  every $k$ and that is Jeffreys' flat correlation for one predictor;
- flat R-squared ($a = 1$) keeps the total flat on the $R^2$ scale.

Keeping the implied $R^2$ constant in $k$ does not by itself single out flat
R: $R^2 \sim$ Beta(½, $b$) is a *g*-mixture for every $k$ and every $b$
(checked in `validation/09`, V5; see
[`2_how_it_calculates.md`](2_how_it_calculates.md), section 8). Where the
three priors come from is described in
[`1_from_behaviour.md`](1_from_behaviour.md), step 0 (study 00).

## The family: mixtures of *g*-priors

Zellner–Siow (JZS, the BayesFactor and JASP default), hyper-*g*,
hyper-*g*/*n*, beta-prime, robust and intrinsic priors all put a prior on
Zellner's *g*. The three priors here do too, with $g/n = B \cdot W$,
$W \sim$ Inverse-Gamma(1, ½) and $B \sim$ Beta($a$, $k/2 - a$). This gives
Beta(½, $(k-1)/2$) for flat R, Beta(1, $k/2 - 1$) for flat R-squared
($k \ge 3$; at $k = 2$ it equals flat partial r, and at $k = 1$ it is not a
*g*-mixture) and $B \equiv 1$ for flat partial r
([`2_how_it_calculates.md`](2_how_it_calculates.md), section 8). The table in
that section shows how much $R^2$ each member expects before any data as $k$
grows: from a median of .25 at every $k$ (flat R) and .50 at every $k$ (flat
R-squared) to .98 at $k = 30$ (Zellner–Siow).

The flat R hyperprior involves Kummer's confluent hypergeometric function. We
did not find it to be a special case of the tCCH family of Li and Clyde
(2018), which unifies the hyper-*g*, beta-prime, robust and intrinsic priors,
but we have not shown that it is not one. A search of the *g*-prior
literature did not turn it up, which is the result of a search and not a
proof that it is new.

## Which priors drift with the number of predictors

A prior **drifts** (= the $R^2$ it implies before seeing data rises as
predictors are added; [`1_from_behaviour.md`](1_from_behaviour.md), step 2).
This happens when each added predictor brings its own expected share of
explained variance.

| Prior | Drifts? | Why |
|---|---|---|
| independent normal / Student-t per coefficient (usual brms/Stan setup; rstanarm's autoscaled normal) | yes | each coefficient adds its own expected share |
| ridge with a fixed penalty | yes | = independent normals |
| uniform box per standardised coefficient | yes | mass sits in the corners |
| uniform over all coefficient vectors with $R^2 \le 1$ | yes | implies $R^2 \sim$ Beta($k/2$, 1) |
| Zellner's g-prior with fixed $g$ (e.g. $g = n$) | yes | $R^2$ rises with $k$ for fixed $g$ |
| JZS / Zellner–Siow (BayesFactor, JASP) | yes | with scale 1, $R^2 \sim$ Beta($k/2$, ½), which rises faster with $k$ than flat partial r (median .98 against .95 at $k = 30$) |
| spike-and-slab with a fixed slab | yes | with the number of included predictors |
| independent flat partial correlations, combined | yes, strongly | $1 - R^2$ is a product of one factor per predictor |
| rstanarm `stan_lm`, and **flat partial r** (η = 1) | yes, by design | $R^2 \sim$ Beta($k/2$, η); η = 1 is "every correlation matrix equally likely" (next section) |
| hyper-g (Liang et al., 2008), Maruyama & George (2011) | partly | prior on the shrinkage $g/(1+g)$ |
| regularised horseshoe | no | global scale tied to the expected number of relevant predictors |
| R2-D2 (Zhang et al., 2022), R2D2M2 (Aguilar & Bürkner, 2023) | no | fixed prior on $R^2$, Dirichlet split |
| **flat R** | no | flat on $R$ for every $k$, even split |
| **flat R-squared** | no | flat on $R^2$ for every $k$, even split |

## Three extensions: flat on the correlations, on R, or on R²

With one predictor, Jeffreys' flat correlation is the starting point. Valid
correlation matrices allow three ways to extend it to $k$ predictors (study
00), and they differ in one claim:

| | what is flat | implied $R^2$ | median $R^2$, $k$ = 1 / 3 / 10 / 30 | Jeffreys at $k = 1$? |
|---|---|---|---|---|
| **flat partial r** | every correlation matrix of $(y, x_1 \ldots x_k)$ | Beta($k/2$, 1) | .25 / .63 / .87 / .95 | yes |
| **flat R** | the total $R$, whatever $k$ is | Beta(½, 1) | .25 / .25 / .25 / .25 | yes |
| **flat R-squared** | the total $R^2$, whatever $k$ is | Beta(1, 1) | .50 / .50 / .50 / .50 | no: $r$ has density $\lvert r\rvert$ |

- **Beta($k/2$, 1) is the distribution of $R^2$ under a uniform correlation
  matrix.** Under the LKJ distribution with shape η (Lewandowski, Kurowicka &
  Joe, 2009; η = 1 is uniform), the squared multiple correlation of one
  variable on the others is Beta($k/2$, η). This is the form of the classical
  null distribution of the sample $R^2$, with η in place of $(n - k - 1)/2$
  (Muirhead, 1982, Theorem 5.2.2). The next subsection
  describes how this was found and checked.
- **It is rstanarm's `stan_lm` prior.** `stan_lm` uses exactly this
  construction: $\sqrt{R^2}$ times a uniformly random direction in
  orthogonalised predictor coordinates (its QR step), with
  $R^2 \sim$ Beta($k/2$, η) (Gabry & Goodrich, rstanarm vignette). All three
  priors here share that construction with `stan_lm` and differ only in the
  prior on $R^2$. `stan_lm` has no default η: the user states where $R^2$ is
  expected and η follows from it, and the rstanarm documentation does not
  recommend η = 1 or use it for Bayes factors. Why η = 1 is the natural choice
  is the subject of the subsection "Why η = 1" below. In this package the
  prior is `prior = "flat_partial_r"` in every function.
- **Consistency between model sizes forces the same growth.** As a
  *g*-mixture (section 8 of [`2_how_it_calculates.md`](2_how_it_calculates.md))
  flat partial r has $B \equiv 1$: a fixed hyperprior $g/n \sim$
  Inverse-Gamma(1, ½), so the prior is the same scale mixture of normals for
  every $k$, and models of different sizes agree with each other. Flat R gives
  this up to keep the total fixed.
- **Zellner–Siow is the η = ½ member.** With scale 1,
  $f^2 = \chi^2_k/\chi^2_1$, so $R^2 \sim$ Beta($k/2$, ½) exactly; the medians
  of .50, .84, .95 and .98 in the table in section 8 of doc 2 are those of
  Beta($k/2$, ½).
- **Keeping every single correlation flat leads to the same $R^2$ prior.**
  Under LKJ(1), single correlations narrow as $k$ grows (the median $|\rho|$ of
  $y$ with one predictor is .50, .35 and .23 for $k$ = 1, 3 and 8). The
  marginally uniform prior of Barnard, McCulloch and Meng (2000), which takes
  the correlations of an inverse-Wishart with $k + 2$ degrees of freedom (the
  dimension plus one), keeps every single correlation flat for any $k$; that
  is their result. We found no source for the further fact that it implies
  the same $R^2 \sim$ Beta($k/2$, 1); it is checked by simulation in
  `docs/figures/make_figures.R` (KS = .0025, .0027 and .0038 for $k$ = 1, 3
  and 8, with 40,000 draws each). Flat correlations in either sense therefore
  give $a = k/2$, and only a flat total gives flat R.

### How Beta($k/2$, 1) was found: simulate, derive, check

The approach was to draw whole correlation matrices with no preference,
compute the $R^2$ each one implies, and compare the spread of those values
with the Beta distribution that a short argument predicts.

1. **Make "the correlations, flat, as $k$ grows" well defined.** The $k$
   correlations of $y$ with the predictors cannot be drawn one at a time:
   together with the correlations among the predictors they must form a valid
   correlation matrix, or $R^2$ can exceed 1. The whole matrix of
   $(y, x_1 \ldots x_k)$ is therefore drawn at once, with every valid matrix
   equally likely (LKJ with η = 1). At $k = 1$ this is a flat correlation, the
   starting point of [`1_from_behaviour.md`](1_from_behaviour.md), step 1.
2. **Draw such matrices.** The LKJ algorithm does not draw the correlations
   directly. It draws **partial correlations** (= the correlation left after
   the earlier variables are taken out), which can be chosen independently
   and always give a valid matrix. Each is a Beta variable on (−1, 1), and they
   are converted back into ordinary correlations by
   `draw_lkj_correlation_matrix()` in `validation/comparison_methods.R` (also
   used by study 00).
3. **Compute each matrix's population $R^2$.** With $\rho$ the correlations of
   $y$ with the predictors and $\Sigma$ the correlations among the predictors:

   ```r
   r_squared <- rho %*% solve(Sigma, rho)   # population R-squared of y on the x's
   ```

   Thousands of matrices at $k$ = 1, 2 and 5 give the distribution of $R^2$
   that "every correlation matrix equally likely" implies.
4. **Predict the formula, then test it.** $R^2$ builds up one predictor at a
   time:

   $$
   1 - R^2 = (1 - p_1^2)(1 - p_2^2) \cdots (1 - p_k^2),
   $$

   with $p_j$ the partial correlations of step 2: independent Beta variables
   whose parameters increase by ½ at each step. A product of such Betas
   collapses into a single Beta, which predicts $R^2 \sim$ Beta($k/2$, η). For
   η = 1 the simulated $R^2$ matches that Beta to within simulation noise at
   $k$ = 1, 2 and 5 (KS = .010–.013; check V4 in study 09). Other values of η
   are not checked for $R^2$.

Each half needs the other: the simulation without the formula gives only
numbers, and the formula without the simulation could hide a slip. The
literal reading "every single correlation flat" was checked in the same way,
with the correlations of an inverse-Wishart with $k + 2$ degrees of freedom
(Barnard et al., 2000), and gives the same Beta($k/2$, 1) by simulation
(`docs/figures/make_figures.R`).

This is the $R^2$ prior of flat partial r. The share per predictor,
Beta(½, $(k-1)/2$) for every member of the family, is a separate result found
in the same way ([`1_from_behaviour.md`](1_from_behaviour.md), step 4).

### Why η = 1

η = 1 is the only shape for which every coefficient's partial correlation is
Jeffreys' flat correlation, whatever the number of other predictors.

| | What it says | Unique to η = 1? | Shown by |
|---|---|---|---|
| 1 | Under LKJ(η), the partial correlation of $y$ with one predictor given all the others is Beta(η, η) on (−1, 1), **for every $k$**. Flat only at η = 1 | yes | `validation/05`, C7 |
| 2 | At $k = 1$ the prior is Jeffreys' flat correlation | yes (LKJ(η) at $k = 1$ is Beta(η, η) on $r$) | the same, $k = 1$ |
| 3 | As a g-mixture, one hyperprior for every $k$: $g/n \sim$ Inverse-Gamma(1, ½), so models of every size share one scale | no (every η gives a fixed hyperprior) | doc 2, section 8 |
| 4 | "Every single correlation flat" (Barnard et al., 2000) leads to the same $R^2 \sim$ Beta($k/2$, 1) | — | simulation (`docs/figures/make_figures.R`); no source found for the $R^2$ part |

Point 1 is what matters for testing one coefficient. Under η = 1 the test of
a coefficient is Jeffreys' correlation test on its partial correlation, with
or without covariates (in its fixed-*X* form; see "What is not new" below),
and it follows from one prior on the whole model rather than from a separate
rule for each test. It also holds for correlated predictors, because the
partial correlation already removes the other predictors.

Point 1 holds for the following reason. The LKJ sampler builds the matrix
from independent partial correlations, level by level, and the partial
correlations at the last level (two variables, given all the others) are
Beta(η, η) on (−1, 1) whatever the number of variables. LKJ gives every
ordering of the variables the same distribution, so $y$ and the tested
predictor can be placed last; their partial correlation given all the others
is then that last-level Beta(η, η).

Simulated KS distances (4,000 matrices per cell; validation run of
4 November 2026, commit 21d42dc):

| | $k = 1$ | $k = 3$ | $k = 8$ |
|---|---|---|---|
| η = 1, to flat | .016 | .013 | .013 |
| η = 2, to flat | .097 | .104 | .101 |
| η = 2, to Beta(2, 2) | .012 | .018 | .014 |

**What is not new.** A Bayes factor for one partial correlation with a flat
prior on it is the test of Kucharský et al. (2026) with α = 1. Our Bayes factor differs from the exact one
of Kucharský et al. in one assumption only (`validation/05`, C8). Theirs is
derived for multivariate normal data, so the predictors are random and the
tested predictor's residual SD is an unknown parameter. Ours conditions on
*X*, as regression does, so the effect is expressed in units of the observed
residual SD. When ours is given that extra uncertainty, the two agree to
within 5e-7, and the priors on the noise and the other coefficients are the
same in both. Without it, ours is lower on the log scale by up to 0.26
(*n* = 20, *k* = 8, *r* = 0), by 0.10 at *n* = 50 and by 0.005 at
*n* = 1,000, with the largest difference at *r* = 0 (where ours does not
depend on *k*). Three things were not found in the literature searched (see
"What is new here, and what is not"): the fixed-*X* version of their test,
the fact that this test agrees with one prior on the whole model, and the use
of the same prior for Bayes factors between whole models and for estimation
(`stan_lm` with η = 1).

**Variable selection** ([`4_validation.md`](4_validation.md), study 09)
compares the three priors with JZS, Zellner–Siow and BIC. They rank real
above null predictors almost equally: flat partial r and flat R-squared are
within .025 of flat R in AUC in every cell, and +.003 and −.002 on average.
They differ in where they draw the line. At *n* = 500, with uncorrelated
predictors and the beta-binomial model prior, flat partial r includes 0.13
of 8 null predictors and misses 0.84 of 4 real ones, close to Zellner–Siow
(0.12 and 0.89) and BIC (0.05 and 1.01). Flat R (0.91 and 0.45) and flat
R-squared (0.92 and 0.44) include more and miss fewer, more so than JZS
(0.42 and 0.58). Which extension describes real research better is partly an
empirical question: how fast explained variance grows with the number of
predictors in published regressions. Flat partial r predicts a median $R^2$
of .76 with 5 predictors and .84 with 8; flat R predicts .25 and flat
R-squared .50 at any number of predictors.

### Flat R-squared at one predictor

Flat R-squared ($a = 1$ for every $k$) differs structurally from the other
two in one respect: at $k = 1$ it is flat on $r^2$ rather than on $r$, so it
does not start from Jeffreys' correlation (study 00, version C).

## The conditional Lindley paradox

Priors with one scale shared by all coefficients (Zellner–Siow, hyper-*g*,
and the joint tests of all three priors here) let large coefficients push the
test of a small one towards "no effect" (Som, Hans & MacEachern, 2016). Their
remedy, block hyper-*g* priors, gives groups of coefficients their own scale.
The conditional tests do the same for one coefficient: under flat R and flat
R-squared it gets its own share of the total, and under flat partial r its
own flat partial correlation
([`2_how_it_calculates.md`](2_how_it_calculates.md), section 10).

## How flat R differs from R2-D2

Both can start from Beta(½, 1) on $R^2$, but they turn it into coefficients
differently, so they are different priors even with one predictor. The same
holds for flat R-squared against R2-D2 with Beta(1, 1).

| | Flat R | R2-D2 (brms `R2D2()`) |
|---|---|---|
| Given a drawn $R^2$ | the length of the effect is **fixed** by it; only the direction is random | $R^2$ sets the **variance** of a normal distribution for the coefficients |
| Implied $R^2$ | follows Beta(½, 1) exactly | does **not** follow Beta(½, 1): the normal draw lands below or above it |
| Split over predictors | uniform direction in explained variance (accounts for correlated predictors) | Dirichlet split per predictor, scaled by each predictor's own SD |
| Correlated predictors | handled in the prior | not taken into account in the split |

With one predictor only the first difference matters. The implied prior on
the correlation, from 1,000,000 draws each:

| | Flat R | brms R2D2 with Beta(½, 1) |
|---|---|---|
| P(\|r\| < .10) | .10 | **.25** |
| P(\|r\| < .50) | .50 | **.65** |
| P(\|r\| > .90) | .10 | .08 |
| median \|r\| | .50 | **.31** |

With nominally the same $R^2$ prior, R2-D2 puts 2.5 times as much weight on
near-zero correlations, so it leans towards small effects. R2-D2 is not part
of this project's comparison.

brms's closest built-in option,

```r
prior(R2D2(mean_R2 = 1/3, prec_R2 = 1.5), class = b)   # Beta(1/2, 1) on R2
```

is the different prior described above. A custom brms prior can state the
flat R density, but brms samples the coefficients directly and misses the
prior's peak at zero for larger $k$. For that reason `fit_regression()` uses
its own Stan model, sampled in $R$ and a direction
([`2_how_it_calculates.md`](2_how_it_calculates.md), section 6).

The one-predictor table can be reproduced with:

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

## How fast evidence for a true null grows

Study 02 (B3b) measures the rate at which the evidence for a true null grows
with the sample size: against the intercept-only model, log BF₀₁ grows like
$a \log n$, so BF₀₁ grows like $n^{a}$. That is the polynomial rate of every
*local* prior (= a prior on the effect whose density is positive at zero), and
Johnson and Rossell (2010) show it is unavoidable for such priors: the
evidence for a true null grows much more slowly than the evidence against a
false one. Their remedy is *non-local* priors (= alternatives with zero
density at the null), which require choosing an effect size of interest; in
psychology, Pramanik and Johnson (2024) argue for the same reason that default
tests confirm true nulls too slowly, and Datta and Johnson (2026) build the
non-local version for partial correlations.

Within the local priors the three branches sit at different points on this
scale. A model with $k$ predictors under a standard local prior such as a
fixed *g*-prior reaches the rate $n^{k/2}$; flat partial r matches it
($a = k/2$), while flat R ($a = ½$) and flat R-squared ($a = 1$) stay at the
rate of a single parameter whatever $k$ is. This is the same behaviour that
study 06 shows for one coefficient: with many predictors, flat R and flat
R-squared rarely confirm a true null. Flat partial r is as fast as a local
prior can be without choosing an effect size; a non-local prior would be
faster, at the cost of that choice.

## What is new here, and what is not

This section is based on a literature search of November 2026. It is not a
systematic review, and the papers listed at the end of it still need a full
reading before anything is written up.

The following points are not new, and are cited as such:

- **The test of one coefficient under flat partial r.** A flat prior on one
  partial correlation, with a closed-form Bayes factor in $n$, $k$ and the
  sample partial correlation, is Kucharský et al. (2026) with α = 1. Ours
  differs in one assumption, fixed instead of random predictors
  (`validation/05`, C8). Other default tests for one partial correlation are
  Wetzels and Wagenmakers (2012) and the Bayes factor functions for partial
  correlations of Datta and Johnson (2026), which use non-local alternatives.
- **Flat partial r as a joint prior.** It is LKJ(1), uniform over correlation
  matrices (Lewandowski et al., 2009), and `stan_lm`'s $R^2 \sim$
  Beta($k/2$, 1) (Gabry & Goodrich). That each partial correlation given all
  the others is then flat follows from the vine construction.
- **A shared scale pushing small coefficients towards zero.** This is the
  conditional Lindley paradox (Som, Hans & MacEachern, 2016); block
  *g*-priors and their Dirichlet-process mixtures (Porwal & Rodríguez, 2024)
  are the proposed remedies.
- **The criteria.** Bayarri et al. (2012); the information inconsistency of
  standard priors in normal linear models is studied by Mulder, Berger, Peña
  and Bayarri (2021).

The following points were not found in the literature searched:

- **The test of one coefficient as the Bayes factor of one joint model.**
  Kucharský et al. justify their test by Jeffreys's desiderata, with separate
  improper priors for the other parameters. In the parts read we found no
  relation to a prior on the whole model (no mention of LKJ, `stan_lm` or
  Savage–Dickey). Here the flat partial r conditional test agrees numerically
  with the Bayes factor of the joint LKJ(1) model (within a factor 1.10 of
  Savage–Dickey on the 16 of 20 datasets where it is estimable;
  `validation/10`). This does not hold for flat R-squared (a factor 2.40 at
  $k = 8$), and under flat R the restricted model is not a proper alternative:
  it collapses to the intercept-only model (study 10, P7). In the terms of
  the literature on **compatible priors** (= priors for nested models derived
  from one prior on the full model; Dawid & Lauritzen, 2001; Consonni &
  Veronese, 2008), this is compatibility by *conditioning*, the notion that
  the Savage–Dickey ratio assumes. It is one of three: the submodel prior can
  also be obtained by marginalisation or by Kullback–Leibler projection, and
  Consonni and Veronese compare the three for *g*-priors without favouring
  conditioning (an earlier version of their paper calls it "not particularly
  advisable"). Conditioning on a single value also depends on the variable
  conditioned on (the Borel–Kolmogorov paradox); here it is the partial
  correlation of the tested predictor. The coherence claim is therefore
  specific: under flat partial r the conditional test equals the Bayes factor
  of the full model against its conditional restriction, and under flat R
  that restriction does not exist. Under projection, flat R would have a
  proper restricted model.
- **How $R^2$-based defaults drift with $k$, and that for the conditional
  tests the drift moves the threshold and leaves the ranking unchanged.** No
  source found shows that, under a fixed prior on $R$ or $R^2$, the
  conditional test of one coefficient changes as other predictors are added
  while the ranking of datasets stays identical (AUC = .701 for all three at
  $k = 20$) and only the share left undecided changes (92% of true nulls under
  flat R; `validation/06`). The joint tests of flat R and flat R-squared do
  separate worse than JZS when predictors are correlated (AUC −.064 at
  $k = 2$, correlation .70). In the parts of Kucharský et al., Som et al. and
  Liang et al. read so far we found no simulations of this kind.
- **Under flat partial r, the conditional restriction is not flat partial r
  with one predictor fewer.** Conditioning the LKJ(1) prior on the tested
  partial correlation being zero leaves the other predictors with
  $R^2 \sim$ Beta($(k-1)/2$, 3/2), which is `stan_lm`'s prior with η = 3/2
  on $k - 1$ predictors, not the flat partial r prior Beta($(k-1)/2$, 1)
  (prior draws with |partial correlation| < .01: KS distance .006 at $k = 3$
  and .008 at $k = 8$ to Beta($(k-1)/2$, 3/2), against .15 and .19 to
  Beta($(k-1)/2$, 1); the medians .37 and .73 are those of study 10, P6).
  This is why flat partial r's joint test, which uses the flat partial r prior
  for the smaller model, differs from its conditional test: the two tests use
  different priors for the model without the tested predictor. It follows
  from the vine construction of LKJ (the submatrix of the other variables is
  LKJ(3/2)); it was checked here by simulation only.
- **The three priors as one family**, $R^2 \sim$ Beta($a$, 1) with
  $a = k/2$, ½ or 1, compared side by side on the Bayarri et al. criteria
  ([`3_comparison.md`](3_comparison.md)).

The ranked list of what to read is in [`reading_list.md`](reading_list.md). Still to be read in full: Kucharský et al. (2026), on whether $k$ is
discussed beyond its role in the formula; Som et al. (2016), for its theorems
on tests of single coefficients; Mulder et al. (2021); the Bayes factor
functions paper; and, most relevant to the coherence claim, Consonni and
Veronese (2008) and Dawid and Lauritzen (2001) on compatible priors across
nested models, which may already contain the conditioning result for flat R
or for LKJ priors. Of lower priority: Steel (2020), a review of Bayesian
model averaging whose discussion of *g*-priors and priors over model size
bears on study 09; and Pramanik and Johnson (2024), for how a psychology
audience is likely to read the slow growth of evidence for nulls.

The review of objective priors by Consonni, Fouskakis, Liseo and Ntzoufras
(2018) has been read (its model-comparison part, sections 3.1 to 3.6). It
does not discuss how the test of one coefficient depends on the number of
other predictors, priors on $R^2$ or on (partial) correlations, LKJ, or
operating characteristics such as how often true nulls are confirmed. It
does cover three points that bear on this project: compatibility of priors
across models (above); the unbalanced learning rates of local priors, with
evidence for a true null growing only as a power of $n$ (Dawid, 2011, as cited there;
Johnson & Rossell, 2010; see "How fast evidence for a true null grows"); and
block hyper-*g* priors for the case where one coefficient is much larger than
the rest (Som et al.), as well as the tendency of hyper-*g* priors to inflate
the inclusion probabilities of unimportant covariates (Dellaportas et al.,
2012, as cited there).

## Closest relatives

- Jeffreys' flat prior on a correlation, and the work of Ly, Verhagen and
  Wagenmakers on Jeffreys' default tests: flat R and flat partial r reduce to
  this prior with one predictor.
- Kucharský, Wagenmakers, van den Bergh and Ly (2026), *TEST*: the exact Bayes
  factor for one partial correlation with a stretched-beta prior. With α = 1
  it is the flat partial r conditional test with random instead of fixed
  predictors (`validation/05`, C8). Flat R puts the same kind of prior on the
  total instead of on each partial correlation. The conditional flat R test
  uses the same inputs (the sample partial correlation, $n$ and $k$), with the
  prior that the flat total implies for one coefficient.
- rstanarm `stan_lm` (Gabry & Goodrich): the same construction with
  $R^2 \sim$ Beta($k/2$, η), used for estimation only; with η = 1 it is flat
  partial r.

## References

- Aguilar, J. E., & Bürkner, P.-C. (2023). Intuitive joint priors for Bayesian
  linear multilevel models: The R2D2M2 prior. *Electronic Journal of Statistics*.
- Barnard, J., McCulloch, R., & Meng, X.-L. (2000). Modeling covariance
  matrices in terms of standard deviations and correlations, with application
  to shrinkage. *Statistica Sinica*.
- Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012). Criteria
  for Bayesian model choice with application to variable selection. *Annals of Statistics*.
- Consonni, G., & Veronese, P. (2008). Compatibility of prior specifications
  across linear models. *Statistical Science*, *23*(3), 332–353.
  https://doi.org/10.1214/08-STS258
- Consonni, G., Fouskakis, D., Liseo, B., & Ntzoufras, I. (2018). Prior
  distributions for objective Bayesian analysis. *Bayesian Analysis*, *13*(2),
  627–679. https://doi.org/10.1214/18-BA1103
- Datta, S., & Johnson, V. E. (2026). Bayes factor functions for testing
  partial correlation coefficients. *Journal of Mathematical Psychology*.
  arXiv:2503.10787.
- Dawid, A. P., & Lauritzen, S. L. (2001). Compatible prior distributions. In
  E. I. George (Ed.), *Bayesian methods with applications to science, policy
  and official statistics* (pp. 109–118). Office for Official Publications of
  the European Communities.
- Gabry, J., & Goodrich, B. Estimating regularized linear models with
  rstanarm (`stan_lm` vignette). https://mc-stan.org/rstanarm/articles/lm.html
- Johnson, V. E., & Rossell, D. (2010). On the use of non-local prior
  densities in Bayesian hypothesis tests. *Journal of the Royal Statistical
  Society: Series B*, *72*(2), 143–170. https://doi.org/10.1111/j.1467-9868.2009.00730.x
- Kucharský, Š., Wagenmakers, E.-J., van den Bergh, D., & Ly, A. (2026).
  Analytic posterior distribution and Bayes factor for Pearson partial
  correlations. *TEST*.
- Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
  correlation matrices based on vines and extended onion method. *Journal of
  Multivariate Analysis*.
- Li, Y., & Clyde, M. A. (2018). Mixtures of g-priors in generalized linear
  models. *JASA*.
- Liang, F., Paulo, R., Molina, G., Clyde, M. A., & Berger, J. O. (2008).
  Mixtures of g priors for Bayesian variable selection. *JASA*.
- Maruyama, Y., & George, E. I. (2011). Fully Bayes factors with a generalized
  g-prior. *Annals of Statistics*.
- Muirhead, R. J. (1982). *Aspects of Multivariate Statistical Theory*. Wiley.
- Mulder, J., Berger, J. O., Peña, V., & Bayarri, M. J. (2021). On the
  prevalence of information inconsistency in normal linear models. *TEST*.
- Porwal, A., & Rodríguez, A. (2024). Dirichlet process mixtures of block g
  priors for model selection and prediction in linear models. arXiv 2411.00471.
- Pramanik, S., & Johnson, V. E. (2024). Efficient alternatives for Bayesian
  hypothesis tests in psychology. *Psychological Methods*, *29*(2), 243–261.
  https://doi.org/10.1037/met0000482
- Som, A., Hans, C. M., & MacEachern, S. N. (2016). A conditional Lindley paradox
  in Bayesian linear models. *Biometrika*.
- Steel, M. F. J. (2020). Model averaging and its use in economics. *Journal
  of Economic Literature*, *58*(3), 644–719. https://doi.org/10.1257/jel.20191385
- Wetzels, R., & Wagenmakers, E.-J. (2012). A default Bayesian hypothesis test
  for correlations and partial correlations. *Psychonomic Bulletin & Review*.
- Zhang, Y. D., Naughton, B. P., Bondell, H. D., & Reich, B. J. (2022).
  Bayesian regression using a prior on the model fit: The R2-D2 shrinkage prior.
  *JASA*.
