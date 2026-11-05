## =====================================================================
##  Scales the timing runs up to the full run
##
##    Rscript cluster/estimate_run_time.R [folder]      (default cache/timing)
##
##  reads the time_NN_at_<fraction>.rds files of cluster/time_study.R and
##  prints a Markdown table: per study, the CPU hours of the full run and the
##  hours it takes on 4, 8, 16 and 32 cores.
##
##  How: a study's cost is a fixed part (compiling, reference draws,
##  bootstraps) plus a part proportional to its datasets. Two timing runs at
##  different fractions give both, by a straight line through them, read off
##  at fraction 1. With one run (or two that used the same number of
##  datasets) the whole time is scaled up, which overestimates. Studies 01
##  and 02 do not depend on the number of datasets, so their run is the full
##  run. The hours on more cores assume the parallel efficiency measured in
##  the timing run; more cores rarely do better, so read them as lower bounds.
## =====================================================================

folder <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(folder)) folder <- file.path("cache", "timing")
files <- list.files(folder, pattern = "^time_.*[.]rds$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No timing runs in ", folder)

runs <- do.call(rbind, lapply(files, function(file) {
  run <- readRDS(file)
  ratios <- vapply(run$n_sim_used, function(n) n[["used"]] / n[["full"]], numeric(1))
  data.frame(study = run$study, fraction = run$fraction, cores = run$cores, cpu_model = run$cpu_model,
             ## the share of the full run's datasets this run used (1 when the
             ## study does not depend on them)
             scale = if (length(ratios) == 0) 1 else mean(ratios),
             elapsed_hours = run$elapsed_seconds / 3600, cpu_hours = run$cpu_seconds / 3600,
             failed = !is.na(run$error))
}))

estimates <- do.call(rbind, lapply(split(runs, runs$study), function(study_runs) {
  study <- study_runs$study[1]
  study_runs <- study_runs[!study_runs$failed, ]
  if (nrow(study_runs) == 0)
    return(data.frame(study = study, method = "failed", cpu_hours = NA, efficiency = NA,
                      timed_minutes = NA))
  largest <- study_runs[which.max(study_runs$scale), ]
  full_at <- function(hours) {
    if (length(unique(study_runs$scale)) >= 2) {
      line <- lm(hours ~ study_runs$scale)
      slope <- coef(line)[[2]]
      if (slope > 0) return(c(coef(line)[[1]] + slope, 0))
    }
    c(hours[which.max(study_runs$scale)] / largest$scale, 1)
  }
  cpu <- full_at(study_runs$cpu_hours)
  data.frame(study = largest$study,
             method = if (largest$scale == 1) "measured in full" else
                      if (cpu[2] == 0) "fixed + per-dataset line" else "scaled up (overestimate)",
             cpu_hours = cpu[1],
             efficiency = largest$cpu_hours / (largest$elapsed_hours * largest$cores),
             timed_minutes = sum(study_runs$elapsed_hours) * 60)
}))

hours_on <- function(cores) estimates$cpu_hours / (cores * pmin(1, estimates$efficiency))
table <- data.frame(Study = estimates$study, Estimate = estimates$method,
                    `Timing runs (min)` = sprintf("%.0f", estimates$timed_minutes),
                    `CPU hours` = sprintf("%.1f", estimates$cpu_hours),
                    `Parallel efficiency` = sprintf("%.0f%%", 100 * estimates$efficiency),
                    `Hours on 4 cores` = sprintf("%.1f", hours_on(4)),
                    `8 cores` = sprintf("%.1f", hours_on(8)),
                    `16 cores` = sprintf("%.1f", hours_on(16)),
                    `32 cores` = sprintf("%.1f", hours_on(32)), check.names = FALSE)
total <- function(cores) sprintf("**%.1f**", sum(hours_on(cores), na.rm = TRUE))
table <- rbind(table, data.frame(Study = "**all**", Estimate = "", `Timing runs (min)` = "",
                                 `CPU hours` = sprintf("**%.0f**", sum(estimates$cpu_hours, na.rm = TRUE)),
                                 `Parallel efficiency` = "", `Hours on 4 cores` = total(4),
                                 `8 cores` = total(8), `16 cores` = total(16), `32 cores` = total(32),
                                 check.names = FALSE))

cat("## Estimated run time of the full validation\n\n")
cat(sprintf("Timed on: %s (%d cores). Fractions of the datasets: %s.\n\n",
            paste(unique(na.omit(runs$cpu_model)), collapse = "; "), runs$cores[1],
            paste(sort(unique(runs$fraction[runs$fraction < 1])), collapse = ", ")))
cat(paste0("| ", paste(names(table), collapse = " | "), " |\n"))
cat(paste0("|", paste(rep("---", ncol(table)), collapse = "|"), "|\n"))
cat(paste0("| ", apply(table, 1, paste, collapse = " | "), " |", collapse = "\n"), "\n\n")
if (any(runs$failed))
  cat("Failed timing runs (left out): ", paste(unique(runs$study[runs$failed]), collapse = ", "), "\n\n")
cat("The hours on more cores are lower bounds: they assume the parallel efficiency of the",
    "timing run. A study runs on all its cores at once, so the whole validation on one",
    "machine takes the sum of the column; as a SLURM array (cluster/slurm_array.sh) it takes",
    "about as long as the slowest study. Apple-silicon cores are faster than these, often by",
    "a factor of 1.5 to 2.\n")
