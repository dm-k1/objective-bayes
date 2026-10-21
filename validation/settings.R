## =====================================================================
##  Shared settings for the validation scripts
##   - parallel: use all cores but one (leaves one free for the computer)
##   - caching:  simulation results and brms fits are saved in
##               validation/results/ and reused when nothing changed
## =====================================================================

## ---- Parallel ----------------------------------------------------------------
## Forking (mclapply) is not available on Windows, so it runs serially there.
number_of_cores <- if (.Platform$OS.type == "windows") 1 else
  max(1, parallel::detectCores() - 1)
options(mc.cores = number_of_cores)   # brms also uses this: chains run in parallel

## ---- Cache folders -------------------------------------------------------------
results_folder <- file.path("validation", "results")
fits_folder    <- file.path(results_folder, "fits")
dir.create(fits_folder, recursive = TRUE, showWarnings = FALSE)

## ---- Cached, parallel simulation runner ---------------------------------------
## Runs `simulate_cell(cell)` for every row of `design`, in parallel, and saves
## the combined results. The cache is linked to the settings: the saved file
## stores `settings` (and the design), and is only reused when they are
## identical. Change REPS, the effect sizes, the sample sizes, ... and the
## simulation reruns automatically. Set rerun = TRUE to force a rerun.
##
## Each design cell gets its own seed (base_seed + row number), so results are
## reproducible and do not depend on the number of cores.
run_simulation_grid <- function(design, simulate_cell, results_file, settings,
                                base_seed = 2026, rerun = FALSE) {
  fingerprint <- list(settings = settings, design = as.data.frame(design))
  if (!rerun && file.exists(results_file)) {
    cached <- readRDS(results_file)
    if (identical(attr(cached, "fingerprint"), fingerprint)) {
      message("Using cached results: ", results_file)
      return(cached)
    }
    message("Settings changed since the cached run: rerunning.")
  }
  message(sprintf("Running %d design cells on %d cores ...", nrow(design), number_of_cores))
  started <- Sys.time()
  cell_results <- parallel::mclapply(seq_len(nrow(design)), function(row) {
    set.seed(base_seed + row)
    simulate_cell(design[row, , drop = FALSE])
  }, mc.cores = number_of_cores)
  failed <- vapply(cell_results, inherits, logical(1), what = "try-error")
  if (any(failed)) stop("Simulation failed in design rows: ", paste(which(failed), collapse = ", "),
                        "\nFirst error: ", cell_results[[which(failed)[1]]])
  results <- dplyr::bind_rows(cell_results)
  attr(results, "fingerprint") <- fingerprint
  saveRDS(results, results_file)
  message(sprintf("Done in %.1f minutes; saved to %s",
                  as.numeric(difftime(Sys.time(), started, units = "mins")), results_file))
  results
}

## ---- Pass/fail checks ------------------------------------------------------------
## check() records one check and prints PASS or FAIL; report_checks() prints the
## summary and stops with an error if anything failed.
check_results <- data.frame(check = character(), passed = logical(), details = character())
check <- function(name, passed, details = "") {
  passed <- isTRUE(passed)
  check_results[nrow(check_results) + 1, ] <<- list(name, passed, details)
  cat(sprintf("[%s] %-55s %s\n", if (passed) "PASS" else "FAIL", name, details))
}
report_checks <- function() {
  cat("\n")
  print(check_results, right = FALSE, row.names = FALSE)
  number_failed <- sum(!check_results$passed)
  if (number_failed > 0) stop(sprintf("%d of %d checks FAILED.", number_failed, nrow(check_results)))
  cat(sprintf("\nAll %d checks passed.\n", nrow(check_results)))
}
