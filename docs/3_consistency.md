# 3. Consistency: arguments behind the default-prior criteria

Section 10 of [`2_how_it_calculates.md`](2_how_it_calculates.md) lists the seven
criteria of Bayarri et al. (2012). This document gives the arguments for the
three that need one. Each is labelled **proof** (complete, short) or **sketch**
(the steps are right, but a careful write-up with all conditions is still to do).
Numerical checks back each one up.

Throughout: the flat-R prior is the mixture of Zellner g-priors with

$$
t = \frac{g}{n} = B \cdot W, \qquad B \sim \mathrm{Beta}\!\left(\tfrac12, \tfrac{k-1}{2}\right),
\qquad W \sim \text{Inverse-Gamma}\!\left(1, \tfrac12\right)
$$

(section 8 of document 2), so the density of $t$ is

$$
p(t) = \frac{1}{2t^2}\, E\!\left[B\, e^{-B/(2t)}\right]
$$

(take the expectation over $B$ of the Inverse-Gamma density of $W = t/B$, times
the change of variables $1/B$). Two bounds on it do all the work below:

- **large $t$:** since $0 < B \le 1$ and $E[B] = 1/k$,
  $\;\dfrac{e^{-1/(2t)}}{2k\,t^2} \;\le\; p(t) \;\le\; \dfrac{1}{2k\,t^2}$;
- **small $t$:** only small $B$ contribute, where the Beta density behaves like
  $b^{-1/2}/\mathrm{B}(\tfrac12, \tfrac{k-1}{2})$, so
  $p(t) \approx c_k\, t^{-1/2}$ as $t \to 0$ for $k \ge 2$, with
  $c_k = \sqrt2\,\Gamma(\tfrac32) / \mathrm{B}(\tfrac12, \tfrac{k-1}{2})$ (checked numerically: for $k = 4$, $\sqrt t\,p(t) = 0.7978$ at $t = 10^{-4}$, against $c_4 = 0.7979$).
  (For $k = 1$, $B = 1$ and $p(t)$ vanishes faster than any power as $t \to 0$.)

The Bayes factor against the intercept-only model is (Liang et al., 2008)

$$
\mathrm{BF}_k = \int_0^\infty (1+nt)^{(n-1-k)/2}\,\bigl(1 + nt(1-R^2)\bigr)^{-(n-1)/2}\, p(t)\, dt .
$$

---

## Information consistency: proof

**Claim.** $\mathrm{BF}_k \to \infty$ as $R^2 \to 1$ if and only if $n \ge k + 3$.

As $R^2 \to 1$ the integrand increases to $(1+nt)^{(n-1-k)/2} p(t)$ (monotone
convergence), so $\mathrm{BF}_k$ tends to $\int_0^\infty (1+nt)^{(n-1-k)/2} p(t)\,dt$.
Near 0 this integral is finite. For large $t$, by the bounds above, the
integrand lies between two constants times $t^{(n-1-k)/2 - 2}$, so the integral
diverges if and only if $(n-1-k)/2 - 2 \ge -1$, that is $n \ge k + 3$. ∎

Checked in `validation/02_bayes_factor.R` (B2b): at $n = k + 2$ the log Bayes
factor levels off (0.91 for $k = 1$, 0.63 for $k = 3$), at $n = k + 3$ it keeps rising.

---

## Intrinsic consistency: proof

**Claim.** As $n \to \infty$, the prior for the coefficients converges to a fixed
proper prior.

The prior for $\beta$ given $\sigma$ depends on the data only through
$\Sigma_n = X_c^\top X_c / n$ (the hyperprior on $t = g/n$ does not involve $n$).
If $\Sigma_n \to \Sigma$, positive definite, as for predictors sampled from a
population, the flat-R density (section 5 of document 2) is continuous in
$\Sigma$ for every $\beta \ne 0$, so it converges pointwise to the flat-R density
with $\Sigma$, which is proper. By Scheffé's lemma the priors converge in total
variation. ∎

---

## Model selection consistency, fixed k: sketch

**Claim.** For nested models with a fixed number of predictors: if the smaller
model is true, $\log \mathrm{BF}_{01}$ grows like $\tfrac12 \log n$; if the larger
model is true, $\log \mathrm{BF}_{10}$ grows linearly in $n$.

*The smaller model is true* (take it to be the intercept-only model; the
nested case reduces to it). Then $nR^2$ stays bounded in probability
($\approx \chi^2_k$). Put $t = s/n$:

$$
\mathrm{BF}_k = \int_0^\infty (1+s)^{-k/2}
\underbrace{\left(\frac{1+s}{1+s(1-R^2)}\right)^{(n-1)/2}}_{\to\ \exp\left(\frac{nR^2}{2}\frac{s}{1+s}\right),\ \text{bounded}}
\; p\!\left(\tfrac{s}{n}\right) \frac{ds}{n} .
$$

With $p(s/n) \approx c_k (s/n)^{-1/2}$ for small $s/n$:

$$
\mathrm{BF}_k \approx n^{-1/2}\, c_k \int_0^\infty (1+s)^{-k/2} s^{-1/2}\,
e^{\frac{nR^2}{2}\frac{s}{1+s}}\, ds ,
$$

and the integral is finite and bounded in probability for $k \ge 2$
(for $k = 1$ the same rate follows from $E[t^{-1/2}] < \infty$). So
$\mathrm{BF}_{10} = O_p(n^{-1/2})$: the evidence for the true null grows like
$\tfrac12 \log n$, **the same rate for every fixed k**. Checked numerically: the
slope of the median $\log \mathrm{BF}_{01}$ against $\log n$ is 0.53
(`validation/02`, B3).

*The larger model is true.* Then $R^2_{\text{large}} \to \rho_1^2 > \rho_0^2 \leftarrow R^2_{\text{small}}$,
and for $t$ bounded away from 0 the ratio of the two integrands behaves like
$\bigl((1-R^2_{\text{small}})/(1-R^2_{\text{large}})\bigr)^{(n-1)/2}$, which grows
exponentially; the prior contributes only polynomial factors in $n$. So
$\log \mathrm{BF}_{10}$ grows linearly in $n$ (`validation/02`, B4: 0.9, 9.3, 54.1
at $n$ = 200, 1000, 5000).

**What a full proof still needs:** the dominated-convergence steps written out
with explicit bounds, and the non-nested case.

---

## When k grows with n: numerical only

The sketch above also shows why this is harder. Under the null,
$nR^2 \approx \chi^2_k$ grows with $k$, so the factor
$\exp\bigl(\tfrac{nR^2}{2}\tfrac{s}{1+s}\bigr)$ is no longer bounded: the evidence
for the null still grows, but more slowly. `validation/08_growing_k.R` checks
this for $k = n/4$ and $k = \sqrt n$. Under the null, with 200 null datasets
per $n$ (median $\log\mathrm{BF}_{01}$):

| $n$ | 50 | 200 | 1000 | 5000 | 20000 |
|---|---|---|---|---|---|
| $k$ fixed at 5 | 1.19 | 1.86 | 2.59 | 3.43 | 4.12 |
| $k = \sqrt n$ | 1.18 | 1.64 | 2.25 | 2.98 | 3.49 |
| $k = n/4$ | 0.99 | 1.34 | 1.80 | 2.15 | 2.49 |

Consistent, but with $k = n/4$ at about half the rate of fixed $k$. Zellner–Siow
accumulates evidence for the null much faster here, because its prior expects
$R^2$ near 1 when $k$ is large, which the data easily rule out; flat-R's
alternative stays at "$R$ flat", which is harder to tell apart from the null.
That is the same budget property, seen from the null side.
