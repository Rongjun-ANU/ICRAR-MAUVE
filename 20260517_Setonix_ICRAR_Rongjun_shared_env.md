# 20260517 Setonix shared ICRAR_Rongjun environment guide

Goal: create a shared Pawsey/Setonix Conda/Mamba environment named `ICRAR_Rongjun`, based on the local `/opt/miniconda3/envs/ICRAR` environment, under the project shared software area:

```text
/software/projects/pawsey1308
```

This guide follows the pattern in:

```text
/Users/Igniz/Desktop/ICRAR/MAUVE/20260322 Setonix Shared Set up.md
```

The important Pawsey constraint from that note is that Conda/Mamba environments should live inside a Singularity overlay, because Conda creates many small files and should not be installed directly as normal files under `/software/projects/pawsey1308`.

## Local source environment

The local environment check on 2026-05-17 found:

- Local Conda: `/opt/miniconda3/bin/conda`, version `26.1.1`
- Local env: `/opt/miniconda3/envs/ICRAR`
- Local Python: `3.13.3`
- Local Mamba: only available after `conda activate ICRAR`, version `2.1.1`
- Local machine: macOS arm64
- Full package source split: `404` total packages, `388` Conda/Mamba-sourced, `16` Pip-sourced
- Full local package list: `/Users/Igniz/Desktop/ICRAR/20260517_ICRAR_local_env_full_package_list.md`
- Machine-readable local reference: `/Users/Igniz/Desktop/ICRAR/20260517_ICRAR_local_conda_list_reference.json`

Because the local environment is macOS arm64 and Setonix is Linux, a literal full Conda export is not a safe Setonix recipe. A test candidate that pinned every local Conda package failed for Linux because it included macOS/runtime packages such as `libcxx` and other platform-level dependencies.

The provided recipe is stricter than a normal `--from-history` export. It pins all Python-facing Conda packages and all Pip packages with exact versions, while letting Mamba choose the Linux-compatible low-level runtime libraries. This is the closest correct match for Setonix: exact at the Python/science package layer, Linux-native below that.

This stricter recipe was dry-run solved for `linux-64` with Conda before writing this guide.

Local recipe file:

```text
/Users/Igniz/Desktop/ICRAR/20260517_ICRAR_Rongjun_environment.yml
```

Verification helper:

```text
/Users/Igniz/Desktop/ICRAR/20260517_compare_ICRAR_env.py
```

## Final Setonix layout

After installation, the relevant shared files should look like:

```text
/software/projects/pawsey1308/
├── bin/
│   └── ICRAR_Rongjun
├── containers/
│   ├── mambaforge_latest.sif
│   └── ICRAR_Rongjun_overlay.img
└── docs/
    ├── 20260517_ICRAR_Rongjun_environment.yml
    ├── 20260517_ICRAR_local_conda_list_reference.json
    └── 20260517_compare_ICRAR_env.py
```

Inside the overlay, the Conda prefix will be:

```text
/ICRAR_Rongjun
```

## 1. Copy the recipe to Setonix

From the local Mac, copy the recipe and verification files to Setonix. Replace `<pawsey_user>` if needed:

```bash
scp \
  /Users/Igniz/Desktop/ICRAR/20260517_ICRAR_Rongjun_environment.yml \
  /Users/Igniz/Desktop/ICRAR/20260517_ICRAR_local_conda_list_reference.json \
  /Users/Igniz/Desktop/ICRAR/20260517_compare_ICRAR_env.py \
  <pawsey_user>@setonix-03.pawsey.org.au:/software/projects/pawsey1308/docs/
```

If direct `scp` into `/software/projects/pawsey1308/docs` is blocked, copy it to your home directory first, then move it after logging in.

## 2. Log in to Setonix

```bash
ssh <pawsey_user>@setonix-03.pawsey.org.au
```

Then set the common variables:

```bash
module load singularity/4.1.0-slurm

BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ICRAR_Rongjun_overlay.img
ENVYML=$BASE/docs/20260517_ICRAR_Rongjun_environment.yml

mkdir -p $BASE/bin
mkdir -p $BASE/containers
mkdir -p $BASE/docs
mkdir -p $SCR/fake_home_ICRAR_Rongjun
mkdir -p $SCR/pip_cache_ICRAR_Rongjun
mkdir -p $SCR/conda_pkgs_ICRAR_Rongjun
mkdir -p $SCR/singularity_cache_ICRAR_Rongjun
mkdir -p $SCR/singularity_tmp_ICRAR_Rongjun
```

## 3. Pull or reuse the shared Mambaforge container

If this file already exists, reuse it:

```text
/software/projects/pawsey1308/containers/mambaforge_latest.sif
```

If it does not exist, create it:

```bash
cd $BASE/containers
singularity pull docker://condaforge/mambaforge:latest
```

## 4. Create the overlay

Use a large overlay because this astronomy Python stack has heavy compiled dependencies.

```bash
cd $BASE/containers
singularity overlay create --size 51200 "$OVERLAY"
```

This creates a 50 GB overlay:

```text
/software/projects/pawsey1308/containers/ICRAR_Rongjun_overlay.img
```

## 5. Build the environment inside the overlay

```bash
module load singularity/4.1.0-slurm

BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ICRAR_Rongjun_overlay.img
ENVYML=$BASE/docs/20260517_ICRAR_Rongjun_environment.yml

singularity exec --cleanenv \
  --bind $BASE:/work \
  --bind $SCR:$SCR \
  --overlay $OVERLAY \
  $IMG \
  bash -lc "
    export HOME=$SCR/fake_home_ICRAR_Rongjun
    export XDG_CACHE_HOME=$SCR/pip_cache_ICRAR_Rongjun
    export PIP_CACHE_DIR=$SCR/pip_cache_ICRAR_Rongjun
    export CONDA_PKGS_DIRS=$SCR/conda_pkgs_ICRAR_Rongjun
    export PYTHONNOUSERSITE=1

    source /opt/conda/etc/profile.d/conda.sh

    mamba env create -p /ICRAR_Rongjun -f /work/docs/20260517_ICRAR_Rongjun_environment.yml
    /ICRAR_Rongjun/bin/python -m pip check
  "
```

Do not silently relax package versions if this fails. The recipe has already solved for `linux-64` in a dry run, so a Setonix failure should be treated as a specific conflict to inspect and document.

## 6. Test the environment

```bash
module load singularity/4.1.0-slurm

BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ICRAR_Rongjun_overlay.img

singularity exec --cleanenv \
  --bind $BASE:/work \
  --bind $SCR:$SCR \
  --overlay $OVERLAY \
  $IMG \
  bash -lc "
    export PYTHONNOUSERSITE=1
    /ICRAR_Rongjun/bin/python -V
    /ICRAR_Rongjun/bin/mamba --version
    /ICRAR_Rongjun/bin/python -m pip check
    /ICRAR_Rongjun/bin/python -c 'import numpy, pandas, astropy, matplotlib, reproject, spectral_cube, ppxf, mpdaf, specutils; print(\"ICRAR_Rongjun import test OK\")'
  "
```

## 6.1. Verify against the local reference

After the environment is created, export the Setonix package list and compare it against the local reference:

```bash
module load singularity/4.1.0-slurm

BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ICRAR_Rongjun_overlay.img

singularity exec --cleanenv \
  --bind $BASE:/work \
  --bind $SCR:$SCR \
  --overlay $OVERLAY \
  $IMG \
  bash -lc "
    /opt/conda/bin/conda list -p /ICRAR_Rongjun --json > /work/docs/20260517_ICRAR_Rongjun_setonix_conda_list.json
    /ICRAR_Rongjun/bin/python /work/docs/20260517_compare_ICRAR_env.py \
      /work/docs/20260517_ICRAR_local_conda_list_reference.json \
      /work/docs/20260517_ICRAR_Rongjun_setonix_conda_list.json
  "
```

Expected result:

```text
OK: Python-facing Conda packages and Pip packages match the local reference.
Low-level system libraries were not required to match because macOS and linux-64 differ.
```

## 7. Create the shared wrapper command

Create:

```text
/software/projects/pawsey1308/bin/ICRAR_Rongjun
```

with:

```bash
cat > /software/projects/pawsey1308/bin/ICRAR_Rongjun <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

module load singularity/4.1.0-slurm >/dev/null 2>&1 || true

BASE=/software/projects/pawsey1308
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ICRAR_Rongjun_overlay.img

SCR=/scratch/pawsey1308/${USER}
FAKE_HOME=$SCR/fake_home_ICRAR_Rongjun
CACHE=$SCR/singularity_cache_ICRAR_Rongjun
TMP=$SCR/singularity_tmp_ICRAR_Rongjun
PIP_CACHE=$SCR/pip_cache_ICRAR_Rongjun
PKGS=$SCR/conda_pkgs_ICRAR_Rongjun

mkdir -p "$FAKE_HOME" "$CACHE" "$TMP" "$PIP_CACHE" "$PKGS"

export SINGULARITY_TMPDIR="$TMP"
export SINGULARITY_CACHEDIR="$CACHE"

exec singularity exec --cleanenv \
  --bind "$BASE:/work" \
  --bind "$SCR:$SCR" \
  --overlay "$OVERLAY" \
  "$IMG" \
  bash -lc '
    export HOME='"$FAKE_HOME"'
    export XDG_CACHE_HOME='"$PIP_CACHE"'
    export PIP_CACHE_DIR='"$PIP_CACHE"'
    export CONDA_PKGS_DIRS='"$PKGS"'
    export PYTHONNOUSERSITE=1
    export PATH=/ICRAR_Rongjun/bin:/opt/conda/bin:$PATH
    exec "$@"
  ' bash "$@"
EOF

chmod 755 /software/projects/pawsey1308/bin/ICRAR_Rongjun
```

## 8. Use the environment interactively

Each user should put the shared bin directory on `PATH`:

```bash
export PATH=/software/projects/pawsey1308/bin:$PATH
```

Then run commands through the wrapper:

```bash
ICRAR_Rongjun python -V
ICRAR_Rongjun mamba --version
ICRAR_Rongjun python -c "import numpy, pandas, astropy; print('ok')"
ICRAR_Rongjun jupyter --version
```

The wrapper is intentionally command-oriented. Use:

```bash
ICRAR_Rongjun python my_script.py
ICRAR_Rongjun ipython
ICRAR_Rongjun jupyter lab --no-browser
```

Do not try to `conda activate /ICRAR_Rongjun` directly from the Setonix login shell, because the env lives inside the overlay.

## 9. Use the environment in Slurm

Example Slurm script:

```bash
#!/bin/bash -l
#SBATCH --account=pawsey1308
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=ICRAR_Rongjun_test
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

module load singularity/4.1.0-slurm
export PATH=/software/projects/pawsey1308/bin:$PATH

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

srun -n 1 -c $SLURM_CPUS_PER_TASK \
  ICRAR_Rongjun python -c "import numpy, pandas, astropy, ppxf, mpdaf, specutils; print('Slurm import test OK')"
```

## 10. Maintainer update procedure

To update the environment later, edit:

```text
/software/projects/pawsey1308/docs/20260517_ICRAR_Rongjun_environment.yml
```

Then run:

```bash
module load singularity/4.1.0-slurm

BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ICRAR_Rongjun_overlay.img

singularity exec --cleanenv \
  --bind $BASE:/work \
  --bind $SCR:$SCR \
  --overlay $OVERLAY \
  $IMG \
  bash -lc "
    export HOME=$SCR/fake_home_ICRAR_Rongjun
    export XDG_CACHE_HOME=$SCR/pip_cache_ICRAR_Rongjun
    export PIP_CACHE_DIR=$SCR/pip_cache_ICRAR_Rongjun
    export CONDA_PKGS_DIRS=$SCR/conda_pkgs_ICRAR_Rongjun
    export PYTHONNOUSERSITE=1

    source /opt/conda/etc/profile.d/conda.sh
    mamba env update -p /ICRAR_Rongjun -f /work/docs/20260517_ICRAR_Rongjun_environment.yml
    /ICRAR_Rongjun/bin/python -m pip check
  "
```

## Practical notes

- The local Mac environment included macOS-only Conda packages, so do not use a literal full `conda env export --no-builds` as the Setonix recipe.
- The Setonix recipe pins the Python/science package layer much more strictly than `--from-history`: all Python-facing Conda packages plus all Pip packages use exact versions.
- The low-level Linux libraries will not be byte-for-byte identical to macOS packages; they must be native Linux packages.
- The verification script checks the part that should match exactly and reports any missing or version-mismatched Python-facing package.
- The environment name for humans is `ICRAR_Rongjun`; the actual Conda prefix inside the overlay is `/ICRAR_Rongjun`.
- The shared command users run is `/software/projects/pawsey1308/bin/ICRAR_Rongjun`.
