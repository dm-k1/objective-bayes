# Running the validation on a cluster

The full validation takes hours (study 9 alone about 6 hours on 7 cores), so
it is meant for a cluster. Two pieces make that work:

- `cluster/slurm_array.sh` runs one study per SLURM array task, with
  `VALIDATION_CORES` set from the task's allocation;
- the container image of `cluster/Dockerfile` (R 4.6.1, the packages of
  `renv.lock`, CmdStan 2.40.0), built and smoke-tested by GitHub Actions and
  pulled with Apptainer, gives every task the same environment as a run
  anywhere else.

The results do not depend on the number of cores: every simulated dataset has
its own random stream, derived from its ID (`docs/appendix/C_study_design.md`).

## How long it takes

The workflow "Run-time estimate" (`.github/workflows/timing.yml`) measures this
on GitHub's runners, for free on a public repository. It starts when that file
changes on `validation-redesign`, and from Actions > Run-time estimate > Run
workflow once it is on `main`. Each study runs in timing mode: its full design and
full-length chains, on 2% and 5% of its datasets, each on its own 4-core runner.
`cluster/estimate_run_time.R` separates the fixed part of a study's cost from
the part per dataset, and reads off the full run. The run's summary page shows
a table with the CPU hours per study and the hours on 4, 8, 16 and 32 cores.
That is enough to choose a machine, or the time limit of the SLURM tasks.

The same can be run anywhere, a study at a time:

```bash
VALIDATION_TIMING=0.02 Rscript cluster/time_study.R 07
VALIDATION_TIMING=0.05 Rscript cluster/time_study.R 07
Rscript cluster/estimate_run_time.R
```

Timing mode writes only to `cache/timing/`, so it never mixes with a real run.

## 1. The Apptainer image

Clusters run Apptainer (formerly Singularity) rather than Docker. The image
is built by GitHub Actions ([`.github/workflows/image.yml`](../.github/workflows/image.yml))
whenever the environment changes: it builds `cluster/Dockerfile` on an
`x86_64` runner, runs every study in smoke mode inside it, and only then
publishes it to the GitHub container registry. On the cluster, pull it once:

```bash
apptainer pull flatcor-validation.sif docker://ghcr.io/dm-k1/flatcor-validation:validation-redesign
```

The tag is the branch the image was built from (`validation-redesign` now;
`latest` once the branch is merged into `main`). The workflow's run on GitHub
(Actions > "Validation image") shows whether the image for the current
`renv.lock` passed its smoke test. If the pull asks for credentials, the
package is still private: make it public in its settings on GitHub
(Packages > flatcor-validation > Package settings), or set
`APPTAINER_DOCKER_USERNAME` and `APPTAINER_DOCKER_PASSWORD` (a token with
`read:packages`) before pulling.

To build the image yourself instead, on an `x86_64` machine with Docker (an
image built on Apple silicon is `arm64` and will not run on an `x86_64`
cluster):

```bash
docker build -f cluster/Dockerfile -t flatcor-validation .
docker save flatcor-validation -o flatcor-validation.tar
```

then copy `flatcor-validation.tar` to the cluster and convert it there:

```bash
apptainer build flatcor-validation.sif docker-archive://flatcor-validation.tar
```

To run a single study in the image by hand, from the project folder:

```bash
apptainer exec --no-home --bind "$PWD":/project --pwd /project \
  --env VALIDATION_CORES=16,R_USER_CACHE_DIR=/project/cache/r_user_cache \
  flatcor-validation.sif Rscript R/studies/07_one_coefficient.R
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
VALIDATION_SMOKE=1 VALIDATION_IMAGE=$PWD/flatcor-validation.sif sbatch cluster/slurm_array.sh
```

runs the smoke test of all ten studies first (minutes per study, with checks
that are shown but mean nothing). Then the full run:

```bash
VALIDATION_IMAGE=$PWD/flatcor-validation.sif sbatch cluster/slurm_array.sh
```

or a selection, for example studies 7 and 9:

```bash
VALIDATION_IMAGE=$PWD/flatcor-validation.sif sbatch --array=7,9 cluster/slurm_array.sh
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
- `results/studies/study_NN.rds` and `figures/`: the study's
  tables, checks and figures, which are tracked.

A study that is interrupted (a time limit, a node failure) resumes where it
stopped when it is submitted again: SimDesign keeps a resume file per study in
`cache/study_NN/`, and the finished parts are cached there.

## 3. After the run

Copy `results/` and `figures/` back (or commit them on
the cluster), and open `validation.ipynb` with the setting
`"results"`: it reads the `.rds` files only and shows every study's checks,
tables and figures. `Rscript run_all.R` on a single large machine
does the same as the array, one study after the other, and also writes
`results/run_summary.md`.

The studies need not run in one go or on one machine: run them in order, each
where it fits (`Rscript run_all.R 03` on a laptop, a cloud machine for 07,
`sbatch --array=8` here), and copy each study's `results/studies/study_NN.rds`,
`results/session_info/study_NN.txt` and figures into the repository. Every
study records its own processor, cores and minutes, and `run_summary.md`
covers every study present.
