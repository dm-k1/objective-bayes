## =====================================================================
##  Shared settings for the validation scripts
##   - parallel: use all cores but one (leaves one free for the computer),
##               or VALIDATION_CORES if set
##   - caching:  simulation results and Stan fits are saved in cache/
##               (gitignored) and reused when nothing changed
##   - figures:  saved to validation/results/figures/ (tracked)
##   - checks:   record_check() and report_checks() for pass/fail summaries
## =====================================================================

## ---- The package under test ------------------------------------------------------
## flatcor lives in its own repository (https://github.com/dm-k1/flatcor) and is
## installed at the version pinned in renv.lock. Some checks also call internal
## functions; they are taken from the package namespace here, in one place.
library(flatcor)
flatcor_internals <- c("sample_r_squared", "prior_shape", "prior_cdf_of_f", "log_sum_exp",
                       "log_sum_exp_per_row", "log_likelihood_ratio_block",
                       "log_bayes_factor_from_prior_cdf_of_f",
                       "log_bayes_factor_shape_vs_intercept_only",
                       "log_bayes_factor_flat_r_vs_intercept_only")
for (internal in flatcor_internals)
  assign(internal, utils::getFromNamespace(internal, "flatcor"), envir = globalenv())

## ---- Parallel ----------------------------------------------------------------
## parallel_map(items, work) runs work(item) for every item on all cores but one.
## It forks the R session (fast, shares everything) where that is allowed: under
## Rscript on macOS and Linux. Where forking is not allowed (Windows; Positron
## and RStudio sessions), it starts fresh R processes instead and copies the
## global objects and attached packages into them.
## VALIDATION_CORES overrides it (a cloud machine or a cluster job that should
## use exactly the cores it was given).
number_of_cores <- as.integer(Sys.getenv("VALIDATION_CORES",
                                         max(1, parallel::detectCores() - 1)))
options(mc.cores = number_of_cores)

forking_allowed <- .Platform$OS.type != "windows" &&
  !inherits(try(parallel::mclapply(1:2, identity, mc.cores = 2), silent = TRUE), "try-error")

## Each item gets its own independent random stream (L'Ecuyer-CMRG, as in
## parallel::nextRNGStream), started from the caller's stream: a set.seed()
## before the call fixes them all, whatever the number of cores, so the same
## script gives the same numbers on any machine. (Seeding each item with
## set.seed(base + i) instead is not safe: consecutive seeds can give draws
## that are not quite uniform across items.)
parallel_map <- function(items, work) {
  start <- sample.int(.Machine$integer.max, 1)
  caller_kind <- RNGkind()[1]
  caller_seed <- get(".Random.seed", envir = globalenv())
  RNGkind("L'Ecuyer-CMRG")
  set.seed(start)
  streams <- vector("list", length(items))
  stream <- get(".Random.seed", envir = globalenv())
  for (index in seq_along(items)) {
    streams[[index]] <- stream
    stream <- parallel::nextRNGStream(stream)
  }
  RNGkind(caller_kind)
  assign(".Random.seed", caller_seed, envir = globalenv())
  seeded_work <- function(index) {
    assign(".Random.seed", streams[[index]], envir = globalenv())
    work(items[[index]])
  }
  if (forking_allowed)
    return(parallel::mclapply(seq_along(items), seeded_work, mc.cores = number_of_cores))
  cluster <- parallel::makeCluster(number_of_cores)
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster, setwd, getwd())
  ## attach the same packages (flatcor among them); the internal functions the
  ## checks use travel with the global objects below
  attached_packages <- .packages()
  parallel::clusterCall(cluster, function(packages) {
    for (package in rev(packages)) library(package, character.only = TRUE)
    NULL
  }, attached_packages)
  parallel::clusterExport(cluster, ls(globalenv()), envir = globalenv())
  parallel::parLapply(cluster, seq_along(items), seeded_work)
}

## ---- Folders -------------------------------------------------------------------
## validation/results/ holds only what is tracked (the run summary, the first
## result, the figures); everything that can be recomputed goes to cache/
results_folder <- file.path("validation", "results")
figures_folder <- file.path(results_folder, "figures")
cache_folder   <- "cache"
fits_folder    <- file.path(cache_folder, "fits")
dir.create(fits_folder,    recursive = TRUE, showWarnings = FALSE)
dir.create(figures_folder, recursive = TRUE, showWarnings = FALSE)

## ---- Figures -------------------------------------------------------------------
## A white background with a light grid, for every study's figures
ggplot2::theme_set(ggplot2::theme_bw())
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
