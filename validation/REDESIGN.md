# Validation redesign

Working plan for the branch `validation-redesign`. It records the decisions
taken so far; the ADEMP descriptions (next step) fill in each study's design
before any code changes. When the redesign is done, this file is replaced by
the finished validation document.

## Why

The current validation (studies 00–11, run of 4 November 2026) is correct,
but it was built study by study:

- The number of simulated datasets per study was chosen by run time, not
  from a target Monte Carlo error, and the docs report results without their
  Monte Carlo standard errors (Morris, White & Crowther, 2019; Siepe et al.,
  2024).
- Several checks compare a maximum or a fixed number against a threshold
  that does not take the Monte Carlo error into account (for example study
  10, P4).
- Each study generates its own data, also where designs overlap, so
  comparisons across studies are not paired.
- The order follows how the studies were built: correctness checks (01, 02,
  05) are interleaved with behaviour studies (03, 04), the default-prior
  criteria are split over 02 and 08, and study 10 mixes calibration with
  coherence.

## Decisions

1. **Established frameworks where they exist, all in R.**
   - `SimDesign` (Chalmers & Adkins, 2020) for every simulation study:
     design grids, generate–analyse–summarise, L'Ecuyer-CMRG random streams,
     and SLURM array jobs (`runArraySimulation()`).
   - `rsimsum` (Gasparini, 2018) for performance measures with their Monte
     Carlo standard errors, including paired differences between methods.
   - `SBC` (Modrák, Kim and others) for simulation-based calibration.
   - `bridgesampling` as an independent route to the joint-model Bayes
     factor.
   - `Rmpfr` for a high-precision ground truth of the closed-form Bayes
     factors.
   - `BAS` (Clyde) for additional reference priors where `bas.lm` offers them
     (to check: whether the robust and intrinsic priors are available for
     linear models).
   - No Python.
2. **The number of datasets follows from a target Monte Carlo standard
   error, fixed in advance** (Morris et al., 2019, Table 6), using the
   worst-case variance where the value is unknown. No stopping rule that
   depends on the results.
3. **Every reported number carries its Monte Carlo standard error**, and every
   check's threshold is set relative to that error.
4. **Shared data-generating mechanisms**, defined by seed (a dataset ID maps to
   its own random stream), so that all methods, and where useful several
   studies, analyse the same datasets:
   - **A. Fixed true effects**: a grid of *k*, *n*, predictor correlation,
     the pattern of the other slopes and the slope of *x*₁. Used by the
     studies on one coefficient, on factors and on estimation.
   - **B. Drawn from each prior**: used by the calibration of the posterior
     and of the Bayes factor.
   - **C. Variable selection**: 12 candidates, 4 real.
   - **D. Growing *k***: *k* rising with *n*.
5. **The motivating simulation (current study 00) moves out of the
   validation**, to `docs/motivating_simulation.R`, with its figures and
   table in `docs/figures/`. It shows why there are three priors; it does
   not validate anything. *Done.*
6. **Order: correctness, then formal properties, then behaviour.**
7. **The package lives in its own repository**,
   [dm-k1/flatcor](https://github.com/dm-k1/flatcor), and this repository
   validates a fixed release of it (v1.0.0, pinned in `renv.lock`). The
   validation loads it with `library(flatcor)`; the internal functions some
   checks need are taken from its namespace in `validation/R/settings.R`.

## The new structure

| New | Question | From the current studies | Data |
|---|---|---|---|
| (motivation) | Why three priors? | 00 | own |
| 1 | Is the prior sampled correctly? | 01, P-checks | none |
| 2 | Are the Bayes factors computed correctly? | 02, A1–A4; 05; new: Rmpfr, bridgesampling | none / few |
| 3 | Is the posterior calibrated? | 01, S (with `SBC`) | B |
| 4 | Is the Bayes factor calibrated? | 10, part 1 | B |
| 5 | Does each prior meet the default-prior criteria? | 02, B1–B5; 08 | none / D |
| 6 | Is the test of one coefficient coherent with one model? | 10, parts 2–3 | from the prior |
| 7 | How does the test of one coefficient behave? | 03, 04, 06 (merged) | A |
| 8 | How does it behave for factors? | 11 | A, extended |
| 9 | How good are the estimates? | 07 | A (the same datasets as 7) |
| 10 | How does it select variables? | 09 | C |

## Steps

1. Install the frameworks and record them in `renv.lock`. *Done (30641ed).*
2. Write the ADEMP description per study, with the target Monte Carlo SE and
   the *n*<sub>sim</sub> that follows (`ADEMP.md`). *Done; approved as drafted.*
3. Shared data sets and study 7 as the worked sample (`precision.R`,
   `data_sets.R`, `07_one_coefficient.R`). *Done.*
4. Port the other studies in the new order. *Done; every port confirmed.*
5. Restructure and present (see "Layout" below), and make it
   cluster-ready: SLURM array scripts, a note on Apptainer, tests for seeding
   and for resuming. *Done, except automated tests for seeding and resuming
   (seeding is confirmed by hand: every study gives the same check lines on 1
   and 4 cores).*
6. Run the full validation on the cluster, then rewrite the docs and the
   report from the new results.

## Status (handoff, 4 November 2026)

Every study is written to the conventions of the sample (ADEMP header,
SimDesign for repeated simulation, `n_sim_of()` targets, Monte Carlo SE on
every estimate, checks relative to it, per-dataset streams from
`with_stream_of()`, cache in `cache/study_NN/`, own SimDesign resume file),
and every one runs end to end in smoke mode, through `run_all.R` and through
the notebook. Smoke mode: `VALIDATION_SMOKE=1 Rscript validation/R/NN_*.R`
(a corner of the design, 10 datasets; checks are shown but do not stop the
run, and results and figures go to `cache/`).

| Study | File | 1 core = 4 cores | Notes |
|---|---|---|---|
| 1 | `01_prior_sampling.R` | yes | full run done: all 10 checks pass (66 s) |
| 2 | `02_bayes_factors.R` | not needed (deterministic) | full run done (6 min): 42 of 45 checks pass; see below |
| 3 | `03_posterior_calibration.R` | yes | SBC package, `future` multisession; the ECDF figure is one panel per prior |
| 4 | `04_bayes_factor_calibration.R` | yes | |
| 5 | `05_default_prior_criteria.R` | yes | |
| 6 | `06_coherence.R` | yes | Stan inside SimDesign, one worker per 4 cores |
| 7 | `07_one_coefficient.R` | yes | |
| 8 | `08_factors.R` | yes | |
| 9 | `09_estimation.R` | yes | about 6 h on 7 cores; a cluster job |
| 10 | `10_variable_selection.R` | yes | |

**Open after the full run of study 2** (decisions for the author):

- *2.6, bridge sampling:* 8 of 9 cases pass. Flat R-squared at k = 1, n = 80
  fails: package -1.7301, bridge sampling -1.7348, a difference of 0.0047
  against 3 x its error of 0.0012; it is the only fit with divergent
  transitions (13). bridgesampling's error treats the draws as independent,
  so it probably understates the error of MCMC draws. Options: a higher
  `adapt_delta` for that case, or an error that allows for autocorrelation.
- *2.2, the 200-bit ground truth:* the closed forms of flat partial r (and of
  flat R-squared at k = 2) are off by up to 3.1 on the log scale where the
  evidence is astronomical (R-squared = .99999, n = 10^6, log BF about
  5.8 million, a relative error of 5e-7); `integrate()` agrees with the
  ground truth to 1e-9. 8 of 90 settings exceed the absolute tolerance of
  1e-6, all with log BF > 1,000. Options: a relative tolerance, a tolerance
  for |log BF| < 1,000 only with the rest reported, or a fix in flatcor.

Known pitfalls, already handled in the ports: SimDesign masks a function named
`summarise` and reserves all-capital names; give every `runSimulation()` its
own resume file (`save_details`); `SimExtract()` needs a `summarise`
function; `rsimsum` needs an `se` column for coverage; Stan models are
compiled into `cache/` (or a temporary folder), never into the installed
package, which is read-only in an Apptainer image.

## Layout (built)

```
validation/
  R/                all code: settings, precision, comparison methods, data sets,
                    studies 01-10, and present.R (the notebook's helpers)
  results/          one study_NN.rds per study: its tables, checks, figures and
                    the run's setup (tracked); run_summary.md from run_all.R
  figures/          every figure as .png (tracked)
  validation.ipynb  the notebook (R kernel): per study, its question, a PASS/FAIL
                    table of its checks, its main tables (gt) and figures; a
                    setting at the top chooses "results" (read the .rds only),
                    "cached" (rerun, reuse cache/) or "fresh" (clear the cache
                    and rerun), and smoke (the smoke test instead)
  run_all.R         the same without the notebook, one study after another
  slurm_array.sh    one study per SLURM array task (see CLUSTER.md)
  CLUSTER.md        running on a cluster, with the Docker image under Apptainer
cache/              raw per-dataset results, logs, smoke results (gitignored)
```

Each study script separates computing from presenting: it ends with
`save_study_results()`, which writes `validation/results/study_NN.rds`; the
notebook only presents. The notebook needs `gt` and `IRkernel` (both in
`renv.lock`; `IRkernel::installspec()` registers the `ir` kernel for R 4.6;
Positron's own R kernel also runs it).

The old figures in `validation/results/figures/` belong to the old studies:
the README and the docs still show them, so they stay until the docs are
rewritten from the new results (step 6).

`slurm_array.sh` is tested here only with mock SLURM variables (one smoke
task); the Apptainer commands of `CLUSTER.md` are untested, as there is no
Apptainer on this machine.

## References

- Chalmers, R. P., & Adkins, M. C. (2020). Writing effective and reliable
  Monte Carlo simulations with the SimDesign package. *The Quantitative
  Methods for Psychology*, *16*(4), 248–280.
- Gasparini, A. (2018). rsimsum: Summarise results from Monte Carlo simulation
  studies. *Journal of Open Source Software*, *3*(26), 739.
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
  studies to evaluate statistical methods. *Statistics in Medicine*, *38*(11),
  2074–2102. https://doi.org/10.1002/sim.8086
- Siepe, B. S., Bartoš, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W., &
  Pawel, S. (2024). Simulation studies for methodological research in
  psychology: A standardized template for planning, preregistration, and
  reporting. *Psychological Methods*. https://doi.org/10.1037/met0000695
