## =====================================================================
##  Run the validation studies, each in its own R process, and summarise
##
##  In one sentence: run the chosen studies in order, keep a log per study,
##  and write one summary of every PASS/FAIL check with the run's setup
##  (commit, R, CmdStan, cores), so a run anywhere can be compared with a run
##  here.
##
##  Run from the project folder:
##    Rscript run_all.R            all studies, 01 to 10 (hours: 07, 09)
##    Rscript run_all.R quick      the correctness checks 01 and 02
##  (The motivating simulation is docs/motivating_simulation.R; it is not part
##  of the validation.)
##    Rscript run_all.R 03 04 09   a selection
##  Cores: all but one by default; set VALIDATION_CORES to choose
##  (e.g. VALIDATION_CORES=32 on a cloud machine).
##
##  The studies are R/studies/NN_*.R; each saves its tables and checks to
##  results/study_NN.rds and its figures to figures/.
##  This script writes the full output of each study to cache/logs/<study>.log
##  (gitignored), and to results/:
##    run_summary.md      the setup, one row per study, and every failed check
##    session_info.txt    R and package versions of this run
##  The notebook validation.ipynb presents the same results.
##  A study that fails does not stop the others. VALIDATION_SMOKE=1 runs every
##  study on a corner of its design and 10 datasets, to test the code; its
##  results, figures and summary then go to cache/ instead.
##  The exit status is 1 if any study failed.
## =====================================================================

## In the order of docs/ADEMP.md: correctness (01-04), formal
## properties (05-06), behaviour (07-10)
all_studies <- c("01_prior_sampling", "02_bayes_factors", "03_posterior_calibration",
                 "04_bayes_factor_calibration", "05_default_prior_criteria", "06_coherence",
                 "07_one_coefficient", "08_factors", "09_estimation", "10_variable_selection")
quick_studies <- c("01_prior_sampling", "02_bayes_factors")

## ---- Which studies -------------------------------------------------------------
requested <- commandArgs(trailingOnly = TRUE)
studies <- if (length(requested) == 0) all_studies else
  if (identical(requested, "quick")) quick_studies else {
    chosen <- all_studies[substr(all_studies, 1, 2) %in% requested]
    unknown <- setdiff(requested, substr(all_studies, 1, 2))
    if (length(unknown) > 0)
      stop("Unknown study number(s): ", paste(unknown, collapse = ", "),
           ". Use 01 to 10, or 'quick'.")
    chosen
  }

smoke_mode     <- identical(Sys.getenv("VALIDATION_SMOKE"), "1")
results_folder <- if (smoke_mode) file.path("cache", "smoke_results") else "results"
logs_folder    <- file.path("cache", "logs")
for (folder in c(results_folder, logs_folder)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)

## ---- The setup of this run ----------------------------------------------------
git_commit <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
                       error = function(e) character(0), warning = function(w) character(0))
git_dirty  <- tryCatch(length(system2("git", c("status", "--porcelain", "--", "R"),
                                      stdout = TRUE, stderr = FALSE)) > 0,
                       error = function(e) NA, warning = function(w) NA)
cmdstan_version <- if (requireNamespace("cmdstanr", quietly = TRUE))
  tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) "not installed") else
  "cmdstanr not installed"
cores <- Sys.getenv("VALIDATION_CORES", as.character(max(1, parallel::detectCores() - 1)))
setup <- c(
  started  = format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
  commit   = if (length(git_commit) == 1) paste0(git_commit, if (isTRUE(git_dirty)) " (with uncommitted changes)") else "unknown",
  R        = R.version.string,
  CmdStan  = cmdstan_version,
  machine  = paste(Sys.info()[["sysname"]], Sys.info()[["machine"]]),
  cores    = cores,
  studies  = paste(substr(studies, 1, 2), collapse = " "),
  smoke    = if (smoke_mode) "yes (results are not results)" else "no")
writeLines(capture.output(sessionInfo()), file.path(results_folder, "session_info.txt"))

## ---- Run each study in its own process -----------------------------------------
## A fresh process per study: no objects carried over, as when run by hand.
run_study <- function(study) {
  log_file <- file.path(logs_folder, paste0(study, ".log"))
  message(sprintf("[%s] %s ...", format(Sys.time(), "%H:%M"), study))
  started <- Sys.time()
  status  <- system2(file.path(R.home("bin"), "Rscript"), file.path("R", "studies", paste0(study, ".R")),
                     stdout = log_file, stderr = log_file,
                     env = paste0("VALIDATION_CORES=", cores))
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  log_lines <- readLines(log_file, warn = FALSE)
  checks <- grep("^\\[(PASS|FAIL)\\]", log_lines, value = TRUE)
  message(sprintf("        exit %d, %.1f minutes, %d checks passed, %d failed", status, minutes,
                  sum(startsWith(checks, "[PASS]")), sum(startsWith(checks, "[FAIL]"))))
  list(study = study, exit_status = status, minutes = minutes,
       passed = sum(startsWith(checks, "[PASS]")), failed = grep("^\\[FAIL\\]", checks, value = TRUE))
}

outcomes <- list()
for (study in studies) outcomes[[study]] <- run_study(study)

## ---- Summary --------------------------------------------------------------------
study_status <- function(outcome)
  if (is.na(outcome$exit_status)) "skipped" else
  if (outcome$exit_status == 0 && length(outcome$failed) == 0) "pass" else "FAIL"
summary_lines <- c(
  "# Validation run",
  "",
  paste0("- **", names(setup), ":** ", setup),
  "",
  "| Study | Status | Minutes | Checks passed | Checks failed |",
  "|---|---|---|---|---|",
  vapply(outcomes, function(o) sprintf("| %s | %s | %.1f | %d | %d |", o$study, study_status(o),
                                       o$minutes, o$passed, length(o$failed)), character(1)),
  "",
  "Every check is judged relative to its Monte Carlo error (docs/ADEMP.md).",
  "Full output: `cache/logs/`.")
all_failed <- unlist(lapply(outcomes, function(o) if (length(o$failed)) paste0(o$study, ": ", o$failed)))
if (length(all_failed) > 0)
  summary_lines <- c(summary_lines, "", "## Failed", "", paste0("- ", all_failed))
writeLines(summary_lines, file.path(results_folder, "run_summary.md"))
cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")

any_failed <- any(vapply(outcomes, study_status, character(1)) != "pass")
quit(status = if (any_failed) 1 else 0)
