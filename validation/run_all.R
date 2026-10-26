## =====================================================================
##  Run the validation studies, each in its own R process, and summarise
##
##  In one sentence: run the chosen studies in order, keep a log per study,
##  and write one summary of every PASS/FAIL check with the run's setup
##  (commit, R, CmdStan, cores), so a run anywhere can be compared with a run
##  here.
##
##  Run from the project folder:
##    Rscript validation/run_all.R            all studies, 01 to 09 (hours: 07)
##    Rscript validation/run_all.R quick      the correctness checks 01, 02, 05
##    Rscript validation/run_all.R 03 04 09   a selection
##  Cores: all but one by default; set VALIDATION_CORES to choose
##  (e.g. VALIDATION_CORES=32 on a cloud machine).
##
##  Writes to validation/results/:
##    logs/<study>.log    the full output of each study
##    run_summary.md      the setup, one row per study, and every failed check
##    session_info.txt    R and package versions of this run
##  A study that fails does not stop the others (except 04, which needs 03).
##  The exit status is 1 if any study failed.
## =====================================================================

all_studies <- c("01_prior_and_posterior", "02_bayes_factor", "03_simulations",
                 "04_matched_comparison", "05_conditional_test", "06_broad_simulation",
                 "07_estimation", "08_growing_k", "09_variable_selection")
quick_studies <- c("01_prior_and_posterior", "02_bayes_factor", "05_conditional_test")

## ---- Which studies -------------------------------------------------------------
requested <- commandArgs(trailingOnly = TRUE)
studies <- if (length(requested) == 0) all_studies else
  if (identical(requested, "quick")) quick_studies else {
    chosen <- all_studies[substr(all_studies, 1, 2) %in% requested]
    unknown <- setdiff(requested, substr(all_studies, 1, 2))
    if (length(unknown) > 0)
      stop("Unknown study number(s): ", paste(unknown, collapse = ", "),
           ". Use 01 to 09, or 'quick'.")
    chosen
  }

results_folder <- file.path("validation", "results")
logs_folder    <- file.path(results_folder, "logs")
dir.create(logs_folder, recursive = TRUE, showWarnings = FALSE)

## ---- The setup of this run ----------------------------------------------------
git_commit <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
                       error = function(e) character(0), warning = function(w) character(0))
git_dirty  <- tryCatch(length(system2("git", c("status", "--porcelain", "--", "R", "inst", "validation/*.R"),
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
  studies  = paste(substr(studies, 1, 2), collapse = " "))
writeLines(capture.output(sessionInfo()), file.path(results_folder, "session_info.txt"))

## ---- Run each study in its own process -----------------------------------------
## A fresh process per study: no objects carried over, as when run by hand.
run_study <- function(study) {
  log_file <- file.path(logs_folder, paste0(study, ".log"))
  message(sprintf("[%s] %s ...", format(Sys.time(), "%H:%M"), study))
  started <- Sys.time()
  status  <- system2(file.path(R.home("bin"), "Rscript"), file.path("validation", paste0(study, ".R")),
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
for (study in studies) {
  if (study == "04_matched_comparison" && "03_simulations" %in% names(outcomes) &&
      outcomes[["03_simulations"]]$exit_status != 0) {
    message("Skipping 04: it needs the results of 03, which failed.")
    outcomes[[study]] <- list(study = study, exit_status = NA, minutes = 0, passed = 0,
                              failed = "skipped: 03 failed")
    next
  }
  outcomes[[study]] <- run_study(study)
}

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
  "Studies 03, 04, 06, 07 and 09 are comparisons (results in their logs and",
  "figures); they have few or no PASS/FAIL checks. Full output: `logs/`.")
all_failed <- unlist(lapply(outcomes, function(o) if (length(o$failed)) paste0(o$study, ": ", o$failed)))
if (length(all_failed) > 0)
  summary_lines <- c(summary_lines, "", "## Failed", "", paste0("- ", all_failed))
writeLines(summary_lines, file.path(results_folder, "run_summary.md"))
cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")

any_failed <- any(vapply(outcomes, study_status, character(1)) != "pass")
quit(status = if (any_failed) 1 else 0)
