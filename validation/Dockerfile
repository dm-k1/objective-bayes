# The environment for running the package and its validation anywhere:
# the same R, the same packages (renv.lock) and the same CmdStan.
#
# The image holds only the environment; the project folder is mounted when
# the container runs, so results are written straight into it.
#
# Build once (20-40 minutes: packages and CmdStan compile from source):
#   docker build -f validation/Dockerfile -t flatcor-validation .
# (from the project folder; Docker reads validation/Dockerfile.dockerignore
# next to this file, which needs BuildKit, the default since Docker 23)
# Run every study, or a selection (see validation/run_all.R):
#   docker run --rm -v "$PWD":/project flatcor-validation
#   docker run --rm -v "$PWD":/project flatcor-validation Rscript validation/run_all.R quick
#   docker run --rm -v "$PWD":/project -e VALIDATION_CORES=32 flatcor-validation

# R version as recorded in renv.lock
FROM rocker/r-ver:4.6.1

# System libraries that the R packages compile against (graphics, fonts,
# networking) and git (run_all.R records the commit)
RUN apt-get update && apt-get install -y --no-install-recommends \
      git make g++ \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
      libpng-dev libjpeg-dev libtiff5-dev \
    && rm -rf /var/lib/apt/lists/*

# R packages at the versions in renv.lock, into the image's site library
# (cmdstanr comes from the Stan r-universe, the rest from CRAN)
COPY renv.lock /tmp/renv.lock
RUN Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")' \
 && Rscript -e 'renv::restore(lockfile = "/tmp/renv.lock", library = .libPaths()[1], prompt = FALSE)'

# CmdStan, the version the validation was last run with
ENV CMDSTAN_VERSION=2.40.0
RUN Rscript -e 'cmdstanr::install_cmdstan(dir = "/opt/cmdstan", version = Sys.getenv("CMDSTAN_VERSION"), cores = parallel::detectCores())'
ENV CMDSTAN=/opt/cmdstan/cmdstan-${CMDSTAN_VERSION}

# git refuses to read a mounted repository owned by another user unless told
RUN git config --system --add safe.directory /project

WORKDIR /project
CMD ["Rscript", "validation/run_all.R"]
