#!/bin/bash
## =====================================================================
##  The validation on a SLURM cluster: one study per array task
##
##  Submit from the project folder (the logs folder must exist first):
##    mkdir -p cache/logs
##    sbatch cluster/slurm_array.sh                 all ten studies
##    sbatch --array=7,9 cluster/slurm_array.sh     a selection
##    VALIDATION_IMAGE=$PWD/flatcor-validation.sif sbatch cluster/slurm_array.sh
##                                                     inside the container
##  Without VALIDATION_IMAGE the studies run with the cluster's own R (which
##  then needs the packages of renv.lock and CmdStan). With it, they run in
##  the Apptainer image built from cluster/Dockerfile (cluster/README.md).
##
##  Each task runs R/studies/NN_*.R with VALIDATION_CORES set to the cores
##  of its allocation, so the study uses exactly those. The numbers do not
##  depend on the number of cores (every dataset has its own random stream).
##  A task writes cache/logs/NN_*.log; the study writes its tables and checks
##  to results/studies/study_NN.rds and its figures to figures/.
##  When all tasks have finished, the notebook (setting "results") shows them.
##
##  The resources below are one setting for every task. Study 9 is the
##  longest (about 6 hours on 7 cores); adjust --time and --cpus-per-task to
##  the cluster, or submit study 9 on its own with more cores.
##  VALIDATION_SMOKE=1 sbatch ... runs the smoke test on the cluster first.
## =====================================================================
#SBATCH --job-name=flatcor-validation
#SBATCH --array=1-10
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=cache/logs/slurm_%A_%a.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

studies=(01_prior_sampling 02_bayes_factors 03_posterior_calibration
         04_bayes_factor_calibration 05_default_prior_criteria 06_coherence
         07_one_coefficient 08_factors 09_estimation 10_variable_selection)
study=${studies[$((SLURM_ARRAY_TASK_ID - 1))]}

export VALIDATION_CORES=$SLURM_CPUS_PER_TASK
export VALIDATION_SMOKE=${VALIDATION_SMOKE:-0}
## Compiled Stan models of flatcor go to R's user cache; one per task, so two
## tasks never compile into the same file at once
export R_USER_CACHE_DIR=$PWD/cache/r_user_cache/task_$SLURM_ARRAY_TASK_ID
mkdir -p cache/logs "$R_USER_CACHE_DIR"
log=cache/logs/$study.log
status=0

echo "$(date '+%F %T') $study on $(hostname), $VALIDATION_CORES cores, smoke=$VALIDATION_SMOKE"
if [[ -n "${VALIDATION_IMAGE:-}" ]]; then
  ## --no-home: R inside sees only the image's packages, not a user library
  ## in the home folder; the project is mounted at /project
  apptainer exec --no-home --bind "$PWD":/project --pwd /project \
    --env VALIDATION_CORES="$VALIDATION_CORES",VALIDATION_SMOKE="$VALIDATION_SMOKE" \
    --env R_USER_CACHE_DIR=/project/cache/r_user_cache/task_$SLURM_ARRAY_TASK_ID \
    "$VALIDATION_IMAGE" Rscript "R/studies/$study.R" > "$log" 2>&1 || status=$?
else
  Rscript "R/studies/$study.R" > "$log" 2>&1 || status=$?
fi
echo "$(date '+%F %T') $study finished: exit $status;" \
     "$(grep -c '^\[PASS\]' "$log" || true) checks passed, $(grep -c '^\[FAIL\]' "$log" || true) failed"
exit $status
