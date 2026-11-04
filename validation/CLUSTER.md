# Running the validation on a cluster

The full validation takes hours (study 9 alone about 6 hours on 7 cores), so
it is meant for a cluster. Two pieces make that work:

- `validation/slurm_array.sh` runs one study per SLURM array task, with
  `VALIDATION_CORES` set from the task's allocation;
- the Docker image of `validation/Dockerfile` (R 4.6.1, the packages of
  `renv.lock`, CmdStan 2.40.0), converted to an Apptainer image, gives every
  task the same environment as a run anywhere else.

The results do not depend on the number of cores: every simulated dataset has
its own random stream, derived from its ID (`validation/ADEMP.md`).

## 1. The Apptainer image

Clusters run Apptainer (formerly Singularity) rather than Docker. Apptainer
can convert the Docker image, but it cannot build from the Dockerfile itself,
so build the Docker image first on a machine with Docker:

```bash
docker build -f validation/Dockerfile -t flatcor-validation .
docker save flatcor-validation -o flatcor-validation.tar
```

Copy `flatcor-validation.tar` to the cluster and convert it there:

```bash
apptainer build flatcor-validation.sif docker-archive://flatcor-validation.tar
```

(If the image is pushed to a registry instead, `apptainer pull
flatcor-validation.sif docker://<registry>/<name>:<tag>` does the same.)
Build the image on an `x86_64` machine if the cluster is `x86_64`: an image built
on Apple silicon is `arm64` and will not run there.

To run a single study in the image by hand, from the project folder:

```bash
apptainer exec --no-home --bind "$PWD":/project --pwd /project \
  --env VALIDATION_CORES=16,R_USER_CACHE_DIR=/project/cache/r_user_cache \
  flatcor-validation.sif Rscript validation/R/07_one_coefficient.R
```

Unlike Docker, Apptainer runs as the user who calls it, and keeps the image
read-only. That is why the command has three flags:

- `--bind "$PWD":/project --pwd /project` mounts the project, where every
  result, figure and cache file is written;
- `--no-home` stops the home folder from being mounted, so R uses only the
  packages in the image and not a personal R library on the cluster;
- `R_USER_CACHE_DIR` puts flatcor's compiled Stan models in the project's
  `cache/` instead of the (now unmounted) home folder.

The studies compile their Stan models into `cache/` or a temporary folder,
never into the read-only image. If a compile still fails with a
"Read-only file system" error, add `--writable-tmpfs`: Apptainer then keeps a
small writable layer in memory for the run, and discards it afterwards.

## 2. The SLURM array

From the project folder on the cluster:

```bash
mkdir -p cache/logs
VALIDATION_SMOKE=1 VALIDATION_IMAGE=$PWD/flatcor-validation.sif sbatch validation/slurm_array.sh
```

runs the smoke test of all ten studies first (minutes per study, with checks
that are shown but mean nothing). Then the full run:

```bash
VALIDATION_IMAGE=$PWD/flatcor-validation.sif sbatch validation/slurm_array.sh
```

or a selection, for example studies 7 and 9:

```bash
VALIDATION_IMAGE=$PWD/flatcor-validation.sif sbatch --array=7,9 validation/slurm_array.sh
```

Without `VALIDATION_IMAGE` the tasks use the cluster's own `Rscript`, which
then needs the packages of `renv.lock` and CmdStan installed. `sbatch` passes
the environment of the submitting shell to the tasks, which is how
`VALIDATION_IMAGE` and `VALIDATION_SMOKE` reach them.

The script asks for 16 cores, 32 GB and 24 hours per task, the same for every
study. Change the `#SBATCH` lines to suit the cluster, or override them when
submitting (`sbatch --cpus-per-task=32 --time=12:00:00 --array=9 ...`).

Each task writes:

- `cache/logs/slurm_<job>_<task>.out`: the start, the end and the number of
  checks passed and failed;
- `cache/logs/NN_<study>.log`: the study's full output;
- `validation/results/study_NN.rds` and `validation/figures/`: the study's
  tables, checks and figures, which are tracked.

A study that is interrupted (a time limit, a node failure) resumes where it
stopped when it is submitted again: SimDesign keeps a resume file per study in
`cache/study_NN/`, and the finished parts are cached there.

## 3. After the run

Copy `validation/results/` and `validation/figures/` back (or commit them on
the cluster), and open `validation/validation.ipynb` with the setting
`"results"`: it reads the `.rds` files only and shows every study's checks,
tables and figures. `Rscript validation/run_all.R` on a single large machine
does the same as the array, one study after the other, and also writes
`validation/results/run_summary.md`.
