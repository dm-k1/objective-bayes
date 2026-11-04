# Reading list

What to read before presenting or writing up this project, in order of
priority. Each entry says why it matters here and what to look for. Status:
**not read**, **partly read** (a summary or some sections), or **read**. All
citations were checked against the publisher or a library record.

The project makes two claims that look new (see
[appendix B](../appendix/B_relation_to_other_priors.md), "What is new here, and what
is not"):

- **Claim A.** Under defaults that keep the total (*R* or *R*²) flat, the
  conditional test of one coefficient confirms true nulls less often as other
  predictors are added, while it ranks datasets exactly as flat partial r
  does: the drift moves the threshold, not the accuracy (study 06).
- **Claim B.** Flat partial r's conditional test is the Bayes factor of the
  full model against its restriction by conditioning, and under flat R that
  restriction collapses to the intercept-only model (study 10).

---

## 1. Read before claiming anything is new

| # | Reference | Why, and what to look for | Status |
|---|---|---|---|
| 1 | Consonni, G., & Veronese, P. (2008). Compatibility of prior specifications across linear models. *Statistical Science*, *23*(3), 332–353. https://doi.org/10.1214/08-STS258 | The closest literature to claim B. Compares marginalisation, conditioning and Kullback–Leibler projection for deriving priors of nested linear models from one full-model prior (with *g*-priors). Look for: whether the degenerate conditional restriction under a prior on the total, or the conditional restriction of LKJ-type priors, is already shown; and their argument against conditioning, which a reader may raise against claim B. Preprint: arXiv:1102.2981. | not read |
| 2 | Dawid, A. P., & Lauritzen, S. L. (2001). Compatible prior distributions. In E. I. George (Ed.), *Bayesian methods with applications to science, policy and official statistics* (pp. 109–118). Office for Official Publications of the European Communities. | The origin of "compatible priors" and of the conditioning-versus-projection distinction, including the dependence of conditioning on the parametrisation (Borel–Kolmogorov). Look for: the definitions to use when stating claim B. Library copy needed. | not read |
| 3 | Kucharský, Š., Wagenmakers, E.-J., van den Bergh, D., & Ly, A. (2026). Analytic posterior distribution and Bayes factor for Pearson partial correlations. *TEST*, *35*(2), 535–558. | The flat partial r conditional test with random instead of fixed predictors (study 05, C8). Look for: any discussion of how the number of control variables affects the evidence beyond its role in the formula (claim A), and of a joint prior for all variables (claim B). | partly read |
| 4 | Rouder, J. N., & Morey, R. D. (2012). Default Bayes factors for model selection in regression. *Multivariate Behavioral Research*, *47*(6), 877–903. https://doi.org/10.1080/00273171.2012.734737 | The source of the JZS regression default in BayesFactor and JASP, the main reference method here. Look for: whether they discuss how the test of one covariate depends on the number or strength of the other covariates (study 04, M4; study 06), which would bear on claim A. | not read |
| 5 | Som, A., Hans, C. M., & MacEachern, S. N. (2016). A conditional Lindley paradox in Bayesian linear models. *Biometrika*, *103*(4), 993–999. | Shows that a shared scale makes the evidence for one coefficient depend on the others. Look for: whether its theorems already cover the conditional tests' independence of what the others explain, and how it relates to claim A. | not read |
| 6 | Pramanik, S., & Johnson, V. E. (2024). Efficient alternatives for Bayesian hypothesis tests in psychology. *Psychological Methods*, *29*(2), 243–261. https://doi.org/10.1037/met0000482 | Argues that default tests confirm true nulls too slowly, for a psychology audience. Look for: whether the dependence on the number of covariates is discussed (claim A), and how to position flat partial r against their non-local alternatives. | not read |

## 2. Know well enough to defend the project

| # | Reference | Why, and what to look for | Status |
|---|---|---|---|
| 7 | Bayarri, M. J., Berger, J. O., Forte, A., & García-Donato, G. (2012). Criteria for Bayesian model choice with application to variable selection. *Annals of Statistics*, *40*(3), 1550–1577. | The seven criteria checked in study 02, and the robust prior. Look for: the exact definitions of each criterion, so that "met" in the docs means what the paper means. | not read |
| 8 | Liang, F., Paulo, R., Molina, G., Clyde, M. A., & Berger, J. O. (2008). Mixtures of *g* priors for Bayesian variable selection. *Journal of the American Statistical Association*, *103*(481), 410–423. | The *g*-prior mixtures the three priors belong to, Zellner–Siow, and information consistency. Look for: their conditions for information consistency (compare with *n* ≥ *k* + 3 here) and model-selection consistency. | not read |
| 9 | Johnson, V. E., & Rossell, D. (2010). On the use of non-local prior densities in Bayesian hypothesis tests. *Journal of the Royal Statistical Society: Series B*, *72*(2), 143–170. https://doi.org/10.1111/j.1467-9868.2009.00730.x | Why evidence for a true null grows only as a power of *n* under local priors; the frame for study 02's rates (*a* log *n*). Look for: the rate result and its conditions. | not read |
| 10 | Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random correlation matrices based on vines and extended onion method. *Journal of Multivariate Analysis*, *100*(9), 1989–2001. | The LKJ distribution and its vine construction: why each partial correlation is flat under LKJ(1), and why its conditional restriction gives the other predictors *R*² ~ Beta((*k*−1)/2, 3/2). Look for: the marginal distribution of a submatrix under LKJ(η). | not read |
| 11 | Wagenmakers, E.-J., Lodewyckx, T., Kuriyal, H., & Grasman, R. (2010). Bayesian hypothesis testing for psychologists: A tutorial on the Savage–Dickey method. *Cognitive Psychology*, *60*(3), 158–189. | The Savage–Dickey ratio used in study 10, and the condition under which it equals the Bayes factor. Look for: that condition, which is exactly compatibility by conditioning. | not read |
| 12 | Ly, A., Verhagen, J., & Wagenmakers, E.-J. (2016). Harold Jeffreys's default Bayes factor hypothesis tests: Explanation, extension, and application in psychology. *Journal of Mathematical Psychology*, *72*, 19–32. | Jeffreys' flat prior on a correlation, the starting point of all three priors at *k* = 1. | not read |

## 3. Check the sources already cited

These are cited in the docs; read enough of each to confirm the docs
describe them correctly.

| # | Reference | What to check | Status |
|---|---|---|---|
| 13 | Mulder, J., Berger, J. O., Peña, V., & Bayarri, M. J. (2021). On the prevalence of information inconsistency in normal linear models. *TEST*, *30*(1), 103–132. | That the information-consistency statements in doc 2 agree with it. | not read |
| 14 | Datta, S., & Johnson, V. E. (2026). Bayes factor functions for testing partial correlation coefficients. *Journal of Mathematical Psychology*. arXiv:2503.10787. | The non-local alternative for partial correlations; how it compares with the flat partial r test. | not read |
| 15 | Li, Y., & Clyde, M. A. (2018). Mixtures of *g*-priors in generalized linear models. *Journal of the American Statistical Association*, *113*(524), 1828–1845. https://doi.org/10.1080/01621459.2018.1469992 | The tCCH family: whether flat R's hyperprior (a ₁F₁ in −*n*/(2*g*)) is or is not a member. The docs only say it was not found to be one. | not read |
| 16 | Zhang, Y. D., Naughton, B. P., Bondell, H. D., & Reich, B. J. (2022). Bayesian regression using a prior on the model fit: The R2-D2 shrinkage prior. *Journal of the American Statistical Association*, *117*(538), 862–874. https://doi.org/10.1080/01621459.2020.1825449 | The comparison with R2-D2 in `relation_to_others.md`. | not read |
| 17 | Barnard, J., McCulloch, R., & Meng, X.-L. (2000). Modeling covariance matrices in terms of standard deviations and correlations, with application to shrinkage. *Statistica Sinica*, *10*(4), 1281–1311. | The marginally uniform prior; the docs cite it for flat single correlations only. | not read |
| 18 | Wetzels, R., & Wagenmakers, E.-J. (2012). A default Bayesian hypothesis test for correlations and partial correlations. *Psychonomic Bulletin & Review*, *19*(6), 1057–1064. | That it is described as a JZS-based test, not a flat prior on the partial correlation. | not read |
| 19 | Schad, D. J., Nicenboim, B., Bürkner, P.-C., Betancourt, M., & Vasishth, S. (2023). Workflow techniques for the robust use of Bayes factors. *Psychological Methods*, *28*(6), 1404–1426. https://doi.org/10.1037/met0000472 | The calibration check of study 10, P1–P2. | not read |
| 20 | Schönbrodt, F. D., & Wagenmakers, E.-J. (2018). Bayes factor design analysis: Planning for compelling evidence. *Psychonomic Bulletin & Review*, *25*(1), 128–142. | The operating characteristics of studies 03 and 06 (how often each answer occurs). | not read |
| 21 | Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018). Validating Bayesian inference algorithms with simulation-based calibration. arXiv:1804.06788. | The calibration check of study 01. | not read |
| 22 | Gabry, J., & Goodrich, B. Estimating regularized linear models with rstanarm (`stan_lm` vignette). https://mc-stan.org/rstanarm/articles/lm.html | That flat partial r is `stan_lm` with η = 1, and how they derive *R*² ~ Beta(*K*/2, η). | not read |

## 4. Background, lower priority

| # | Reference | Why | Status |
|---|---|---|---|
| 23 | Consonni, G., Fouskakis, D., Liseo, B., & Ntzoufras, I. (2018). Prior distributions for objective Bayesian analysis. *Bayesian Analysis*, *13*(2), 627–679. https://doi.org/10.1214/18-BA1103 | The survey of objective priors; it led to entries 1 and 2. Its model-comparison part does not discuss claims A or B. | read (sections 3.1–3.6) |
| 24 | Steel, M. F. J. (2020). Model averaging and its use in economics. *Journal of Economic Literature*, *58*(3), 644–719. https://doi.org/10.1257/jel.20191385 | *g*-priors and priors over model size in model averaging; bears on study 09. | not read |

---

**Suggested order:** 1, 4 and 3 first (they decide whether claims A and B are
new and how to state them), then 11 and 10 (the two tools claim B rests on),
then 7, 8 and 9 (to defend the criteria and rates), then the rest as needed.
