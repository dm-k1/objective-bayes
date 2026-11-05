## Estimated run time of the full validation

Timed on: AMD EPYC 7763 64-Core Processor; AMD EPYC 9V45 96-Core Processor; AMD EPYC 9V74 80-Core Processor; Intel(R) Xeon(R) 6975P-C; INTEL(R) XEON(R) PLATINUM 8573C (4 cores). Fractions of the datasets: 0.02, 0.05.

| Study | Estimate | Timing runs (min) | Core hours | Hours on 4 cores | 8 cores | 16 cores | 32 cores |
|---|---|---|---|---|---|---|---|
| 01 | measured in full | 1 | 0.1 | 0.0 | 0.0 | 0.0 | 0.0 |
| 02 | measured in full | 10 | 0.7 | 0.2 | 0.1 | 0.0 | 0.0 |
| 03 | fixed + per-dataset line | 3 | 0.7 | 0.2 | 0.1 | 0.0 | 0.0 |
| 04 | fixed + per-dataset line | 23 | 22.3 | 5.6 | 2.8 | 1.4 | 0.7 |
| 05 | fixed + per-dataset line | 2 | 1.7 | 0.4 | 0.2 | 0.1 | 0.1 |
| 06 | scaled up (overestimate) | 11 | 1.8 | 0.4 | 0.2 | 0.1 | 0.1 |
| 07 | fixed + per-dataset line | 24 | 20.4 | 5.1 | 2.6 | 1.3 | 0.6 |
| 08 | fixed + per-dataset line | 44 | 39.1 | 9.8 | 4.9 | 2.4 | 1.2 |
| 09 | fixed + per-dataset line | 35 | 32.7 | 8.2 | 4.1 | 2.0 | 1.0 |
| 10 | fixed + per-dataset line | 6 | 1.6 | 0.4 | 0.2 | 0.1 | 0.0 |
| **all** |  |  | **121** | **30.3** | **15.1** | **7.6** | **3.8** | 

The hours on more cores are lower bounds: they assume every core is used as well as in the timing run. A study runs on all its cores at once, so the whole validation on one machine takes the sum of the column; as a SLURM array (cluster/slurm_array.sh) it takes about as long as the slowest study. Apple-silicon cores are faster than these, often by a factor of 1.5 to 2.

Source: the timing files of workflow run 42030390433 (5 November 2026, the workflow of commit f2a972f), scaled with `cluster/estimate_run_time.R` as of commit 3f9c104 (elapsed hours times cores). Rerun with `gh run download 42030390433 --pattern 'time-*' --dir timing` and `Rscript cluster/estimate_run_time.R timing`.
