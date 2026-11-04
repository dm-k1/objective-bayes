# What this project shows

**The claim.** For testing one regression coefficient with a default
("objective") Bayes factor, it matters what the default keeps flat as
predictors are added. Keep each partial correlation flat and the test means
the same at every number of predictors $k$. Keep the total ($R$ or $R^2$)
flat and the test drifts with $k$, and the drift is not a gain in accuracy:
the test sees the same, it only decides less.

All numbers come from the validation run of 4 November 2026 (commit
643c0b8; 131 of 132 checks pass), with references to the study that
produced them. The full results are in [`5_validation.md`](5_validation.md)
and [`../validation/report.html`](../validation/report.html).

---

## 1. Not everything can stay flat (study 00)

With one predictor, Jeffreys' flat correlation is the natural default. With
several, draw valid correlation matrices of $(y, x_1 \ldots x_k)$ and the
choice is forced: keep each partial correlation flat and $R^2$ grows with $k$
(P($R^2$ > 0.5) from 0.29 at $k = 1$ to 0.97 at $k = 10$); keep $R$ or $R^2$
flat and each partial correlation shrinks (P(|partial r| > 0.3) from 0.70 to
0.22 under flat $R$). That leaves three defaults, one family
$R^2 \sim$ Beta($a$, 1):

| Prior | $a$ | What stays flat |
|---|---|---|
| flat partial r | $k/2$ | each partial correlation (LKJ(1); `stan_lm` with η = 1) |
| flat R | ½ | the multiple correlation $R$ |
| flat R-squared | 1 | $R^2$ |

![Study 00: what valid correlation matrices force](../validation/results/figures/00_partial_r_R_R2_by_k.png)

## 2. All three are correct and meet the formal criteria (studies 01, 02, 05)

This is not a contest of a right against a wrong implementation. All three
pass every correctness check (Bayes factors against numerical integration,
brute-force sampling and a second derivation; simulation-based calibration of
the Stan model) and every criterion of Bayarri, Berger, Forte &
García-Donato (2012): predictive matching, invariance, information
consistency, model-selection consistency, intrinsic consistency. The
difference between them is behaviour.

## 3. The test drifts with k under flat R and flat R-squared, but sees the same (studies 04, 06)

Across 216 designs ($k$ = 2–20, $n$ = 50–1000), the three conditional tests
**rank datasets identically**: the AUC for separating a true null from a
small effect (slope 0.1) is the same for all three (0.701 at $k = 20$). None
is more accurate. What differs is what they *say* at BF 3:

| Truth, $k = 20$ | flat partial r | flat R | flat R-squared |
|---|---|---|---|
| null: says null / undecided / says effect | **0.88** / 0.11 / 0.01 | **0.06** / 0.92 / 0.02 | **0.26** / 0.72 / 0.02 |
| slope 0.1: says null / undecided / says effect | 0.60 / 0.26 / **0.15** | 0.005 / 0.78 / **0.21** | 0.04 / 0.75 / **0.21** |

Under flat partial r the null row is the same at $k = 2$ (0.86) as at
$k = 20$ (0.88). Under flat R it falls from 0.55 to 0.06: with 20 predictors
it leaves 92% of true nulls undecided. Its low rate of calling small effects
null comes from rarely saying anything, not from seeing effects better.

![Study 06: what each test says, as k grows](../validation/results/figures/06_outcomes.png)

## 4. Only flat partial r is one coherent model (study 10)

A test of one coefficient should be the Bayes factor of the joint model
against the same model with that coefficient zero (the Savage–Dickey ratio).
For flat partial r it is, within a factor 1.10. For flat R-squared it drifts
from it as $k$ grows (a factor 2.40 at $k = 8$; the one failed check of the
run, a property of that prior: its partial correlations are tied together
through $R^2$). For flat R the ratio is not defined.

![Study 10: the conditional and joint tests against the joint model](../validation/results/figures/10_coherence.png)

## 5. What flat partial r costs

- **It leans to the null for small effects with little data**: at $k = 20$ it
  calls 60% of 0.1 effects null. A BF₀₁ > 3 from it does not rule out a small
  effect; report it with the support curve (`R/support_curve.R`).
- **It shrinks less, so it predicts somewhat worse** at small $n$ with many
  predictors (prediction error 0.72–0.81 of OLS at $k = 8$, $n = 30$, against
  0.37–0.77 for flat R; study 07). In return its estimates are less biased
  (−8% of a real effect, against −28% for flat R) and its intervals keep
  coverage at 0.916 or above (flat R: down to 0.760).
- **In variable selection it is conservative**: few null predictors in,
  more real ones missed, much like Zellner–Siow (study 09).

## 6. What is new, and what is not

The test itself is not new: a flat prior on one partial correlation is
Kucharský et al. (2026) with α = 1, and the joint prior is LKJ(1), `stan_lm`
with η = 1. Not found in the literature searched: that this test is the
Bayes factor of one joint model, and that the drift of the $R$- and
$R^2$-based defaults with $k$ is a shift in where they decide, not in what
they see. Details and the papers still to read:
[`relation_to_others.md`](relation_to_others.md), "What is new here, and
what is not".

---

**Recommendation** ([`4_comparison.md`](4_comparison.md), section 8): flat
partial r, with the costs above stated when it is used. The package keeps no
default; every function asks for `prior =`.
