## =====================================================================
##  Times one study on a fraction of its datasets
##
##    VALIDATION_TIMING=0.02 Rscript cluster/time_study.R 07
##
##  runs R/studies/07_*.R in timing mode (R/shared/settings.R): the full
##  design and full-length chains, but 2% of the datasets (at least 10).
##  It writes cache/timing/time_07_at_0.02.rds with the elapsed time, the CPU
##  time of R and of every process it started (forked workers, CmdStan), and
##  the datasets the study used against its full number.
##  cluster/estimate_run_time.R scales these up to the full run; the
##  "Run-time estimate" workflow (.github/workflows/timing.yml) runs both.
##
##  Studies 01 and 02 do not depend on the number of datasets
##  (R/shared/precision.R): timing them is running them in full.
## =====================================================================

.timing_study <- commandArgs(trailingOnly = TRUE)[1]
.timing_file  <- Sys.glob(sprintf("R/studies/%s_*.R", .timing_study))
stopifnot("give a study number, as in Rscript cluster/time_study.R 07" = length(.timing_file) == 1,
          "set VALIDATION_TIMING to the fraction of datasets, as in 0.02" =
            !is.na(suppressWarnings(as.numeric(Sys.getenv("VALIDATION_TIMING", NA)))))

.timing_start <- proc.time()
.timing_error <- tryCatch({ source(.timing_file); NA_character_ },
                          error = function(e) conditionMessage(e))
.timing_used  <- proc.time() - .timing_start

.timing_result <- list(
  study    = .timing_study,
  fraction = as.numeric(Sys.getenv("VALIDATION_TIMING")),
  cores    = as.integer(Sys.getenv("VALIDATION_CORES", max(1, parallel::detectCores() - 1))),
  cpu_model = tryCatch(sub(".*:\\s*", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]),
                       error = function(e) NA_character_, warning = function(w) NA_character_),
  elapsed_seconds = unname(.timing_used[["elapsed"]]),
  cpu_seconds = unname(sum(.timing_used[c("user.self", "sys.self", "user.child", "sys.child")], na.rm = TRUE)),
  n_sim_used = get0("n_sim_used", ifnotfound = list()),
  error = .timing_error)

dir.create(file.path("cache", "timing"), recursive = TRUE, showWarnings = FALSE)
saveRDS(.timing_result, file.path("cache", "timing",
                                  sprintf("time_%s_at_%s.rds", .timing_study, .timing_result$fraction)))
cat(sprintf("\nStudy %s at fraction %s: %.1f minutes elapsed, %.2f CPU hours%s\n",
            .timing_study, .timing_result$fraction, .timing_result$elapsed_seconds / 60,
            .timing_result$cpu_seconds / 3600,
            if (is.na(.timing_error)) "" else paste0("; stopped with an error: ", .timing_error)))
if (!is.na(.timing_error)) quit(status = 1)
