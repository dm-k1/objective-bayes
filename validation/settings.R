## =====================================================================
##  Shared settings for the validation scripts
##   - parallel: use all cores but one (leaves one free for the computer)
##   - caching:  simulation results are saved in
##               validation/results/ and reused when nothing changed
##   - figures:  saved to validation/results/figures/
##   - checks:   record_check() and report_checks() for pass/fail summaries
## =====================================================================

## ---- Parallel ----------------------------------------------------------------
## parallel_map(items, work) runs work(item) for every item on all cores but one.
## It forks the R session (fast, shares everything) where that is allowed: under
## Rscript on macOS and Linux. Where forking is not allowed (Windows; Positron
## and RStudio sessions), it starts fresh R processes instead and copies the
## global objects and attached packages into them.
number_of_cores <- max(1, parallel::detectCores() - 1)
options(mc.cores = number_of_cores)

forking_allowed <- .Platform$OS.type != "windows" &&
  !inherits(try(parallel::mclapply(1:2, identity, mc.cores = 2), silent = TRUE), "try-error")

parallel_map <- function(items, work) {
  if (forking_allowed) return(parallel::mclapply(items, work, mc.cores = number_of_cores))
  cluster <- parallel::makeCluster(number_of_cores)
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster, setwd, getwd())
  attached_packages <- .packages()
  parallel::clusterCall(cluster, function(packages)
    for (package in rev(packages)) library(package, character.only = TRUE), attached_packages)
  parallel::clusterExport(cluster, ls(globalenv()), envir = globalenv())
  parallel::parLapply(cluster, items, work)
}

## ---- Cache folders -------------------------------------------------------------
results_folder <- file.path("validation", "results")
fits_folder    <- file.path(results_folder, "fits")
figures_folder <- file.path(results_folder, "figures")
dir.create(fits_folder,    recursive = TRUE, showWarnings = FALSE)
dir.create(figures_folder, recursive = TRUE, showWarnings = FALSE)

## ---- Figures -------------------------------------------------------------------
## Saves a ggplot to validation/results/figures/<name>.png, and also shows it
## when running interactively. (Under Rscript a printed plot would otherwise
## end up in Rplots.pdf in the project folder.)
save_figure <- function(figure, name, width = 10, height = 6) {
  ggplot2::ggsave(file.path(figures_folder, paste0(name, ".png")), figure,
                  width = width, height = height, dpi = 150)
  if (interactive()) print(figure)
  invisible(figure)
}

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
  cell_results <- parallel_map(seq_len(nrow(design)), function(row) {
    set.seed(base_seed + row)
    simulate_cell(design[row, , drop = FALSE])
  })
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
## record_check() records one check and prints PASS or FAIL; report_checks()
## prints the summary and stops with an error if anything failed.
check_results <- data.frame(check = character(), passed = logical(), details = character())
record_check <- function(name, passed, details = "") {
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
