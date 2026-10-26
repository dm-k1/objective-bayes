# Validation run

- **started:** 2026-10-26 16:00 CET
- **commit:** 606c955 (with uncommitted changes)
- **R:** R version 4.6.1 (2026-06-24)
- **CmdStan:** 2.40.0
- **machine:** Darwin arm64
- **cores:** 7
- **studies:** 01 02 05

| Study | Status | Minutes | Checks passed | Checks failed |
|---|---|---|---|---|
| 01_prior_and_posterior | pass | 7.1 | 13 | 0 |
| 02_bayes_factor | pass | 0.7 | 25 | 0 |
| 05_conditional_test | pass | 0.1 | 7 | 0 |

Studies 03, 04, 06, 07 and 09 are comparisons (results in their logs and
figures); they have few or no PASS/FAIL checks. Full output: `logs/`.
