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
   checks need are taken from its namespace in `validation/settings.R`.

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

1. Install the frameworks and record them in `renv.lock`. *Done on `main`
   (commit 30641ed).*
2. Write one ADEMP description per new study (aims, data-generating
   mechanism, estimands, methods, performance measures), with the target Monte
   Carlo standard error and the number of datasets that follows from it.
   **For review before any code changes.**
3. Build the shared data-generating mechanisms and port study 7 as a worked
   sample with `SimDesign` and `rsimsum`; test it with a small run.
4. Port the other studies in the new order.
5. Make it cluster-ready: a smoke mode that runs every study end to end in
   minutes, a SLURM array script per study, and a note on running the Docker
   image with Apptainer. Tests for the seeding (same results on any number of
   cores) and for resuming interrupted runs.
6. Run the full validation on the cluster, then rewrite the docs and the
   report from the new results.

## Open questions

- Whether the explaining docs (docs 0–4, `findings.md`) are still needed
  next to the report, or can be folded into it.

- The target Monte Carlo standard error per quantity (step 2).
- Whether `BAS` offers the robust and intrinsic priors for linear models.

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
