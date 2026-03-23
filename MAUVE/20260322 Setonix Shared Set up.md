# 20260322 Setonix Shared Set up

Due to pawsey policy, we have to use containers to avoid installing large numbers of small files (common with Conda/Mamba) directly. 

Below is a **clean, corrected project-wide installation guide** for reinstalling the shared stack under:

```text
/software/projects/pawsey1308
```

Clarifications:

- shared/common things go directly under `/software/projects/pawsey1308`
- private CADC credentials should live under each user’s own existing folder, e.g.
  - `/software/projects/pawsey1308/$USER/cadc_home`
- shared wrapper commands should live in:
  - `/software/projects/pawsey1308/bin`

------

# Project-wide shared software setup for pawsey1308

## 0. Final layout

After setup, the shared/common structure should look like this:

```text
/software/projects/pawsey1308/
├── bin/
│   ├── conda1308
│   ├── mamba1308
│   ├── ngistenv1308
│   ├── mapviewer1308
│   ├── cadc-get-cert
│   ├── vcp
│   └── vls
├── containers/
│   ├── mambaforge_latest.sif
│   ├── ngistenv1308_overlay.img
│   └── cadc_overlay.img
├── docs/
├── ngist/
├── ngist_supplementary_public/
├── setonix/                 # optional helper scripts / slurm templates
├── rhuang/                  # personal folder, already exists
```

### Meaning of each area

- **`/software/projects/pawsey1308/ngist`**
  shared nGIST repo
- **`/software/projects/pawsey1308/ngist_supplementary_public`**
  shared supplementary repo
- **`/software/projects/pawsey1308/containers`**
  shared container image and shared overlays
- **`/software/projects/pawsey1308/bin`**
  shared wrapper commands for everyone
- **`/software/projects/pawsey1308/$USER/cadc_home`**
  each user’s private CADC certificate store
- **`/scratch/pawsey1308/$USER/...`**
  per-user caches, logs, personal softwares

------

# Part A — Maintainer installation steps (For Rongjun only)

This section is for the person setting up the shared project-wide stack.

------

## 1. Create the shared directories

Run:

```bash
BASE=/software/projects/pawsey1308

mkdir -p $BASE/bin
mkdir -p $BASE/containers
mkdir -p $BASE/docs
mkdir -p $BASE/setonix
```

------

## 2. Clone the shared repositories

## 2.1 nGIST repo

Choose whether the project should track `main`, `dev-branch`, or a specific tag.
For a shared install, I recommend using **`main`** unless the project explicitly wants `dev-branch`.

```bash
cd /software/projects/pawsey1308
git clone https://github.com/geckos-survey/ngist.git
cd ngist
git checkout main
git pull --ff-only origin main
```

------

## 2.2 supplementary repo

```bash
cd /software/projects/pawsey1308
git clone https://github.com/geckos-survey/ngist_supplementary_public.git
cd ngist_supplementary_public
git checkout main
git pull --ff-only origin main
```

------

## 3. Pull the shared base container

We will use one shared base Singularity image for both the nGIST environment and the CADC tools.

```bash
module load singularity/4.1.0-slurm

cd /software/projects/pawsey1308/containers
singularity pull docker://condaforge/mambaforge:latest
```

This should create:

```text
/software/projects/pawsey1308/containers/mambaforge_latest.sif
```

------

## 4. Create the shared nGIST overlay

This overlay will contain the shared environment `/ngistenv1308`.

```bash
BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ngistenv1308_overlay.img

mkdir -p $SCR/fake_home_1308
mkdir -p $SCR/pip_cache_1308
mkdir -p $SCR/conda_pkgs_1308

cd $BASE/containers
singularity overlay create --size 51200 "$OVERLAY"
```

The overlay is used because the environment contains lots of small files, and we do not want them to count as normal files on `/software`.

------

## 5. Build the shared nGIST environment inside the overlay

Install the environment as:

```text
/ngistenv1308
```

inside the overlay.

```bash
module load singularity/4.1.0-slurm

BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ngistenv1308_overlay.img

singularity exec --cleanenv \
  --bind $BASE:/work \
  --bind $SCR:$SCR \
  --overlay $OVERLAY \
  $IMG \
  bash -lc "
    export HOME=$SCR/fake_home_1308
    export XDG_CACHE_HOME=$SCR/pip_cache_1308
    export PIP_CACHE_DIR=$SCR/pip_cache_1308
    export CONDA_PKGS_DIRS=$SCR/conda_pkgs_1308
    export PYTHONNOUSERSITE=1

    source /opt/conda/etc/profile.d/conda.sh

    mamba env create -f /work/ngist_supplementary_public/ngistenv.yaml -p /ngistenv1308
    /ngistenv1308/bin/python -m pip install -e /work/ngist --no-deps
  "
```

### Test the environment

```bash
module load singularity/4.1.0-slurm

BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ngistenv1308_overlay.img

singularity exec --cleanenv \
  --bind $BASE:/work \
  --bind $SCR:$SCR \
  --overlay $OVERLAY \
  $IMG \
  bash -lc "
    /ngistenv1308/bin/python -V
    /ngistenv1308/bin/python -m pip show ngistPipeline | egrep 'Name|Version|Location'
    /ngistenv1308/bin/ngistPipeline --help | head
  "
```

------

## 6. Create the shared CADC overlay

This overlay contains the CADC/VOS tools, but **not** shared credentials.

```bash
BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
CADC_OVERLAY=$BASE/containers/cadc_overlay.img

mkdir -p $SCR/fake_home_cadc
mkdir -p $SCR/pip_cache_cadc
mkdir -p $SCR/conda_pkgs_cadc

cd $BASE/containers
singularity overlay create --size 4096 "$CADC_OVERLAY"
```

------

## 7. Install CADC/VOS tools into the CADC overlay

Install them into `/cadcenv` inside the overlay.

```bash
module load singularity/4.1.0-slurm

BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
CADC_OVERLAY=$BASE/containers/cadc_overlay.img

singularity exec --cleanenv \
  --bind $BASE:/work \
  --bind $SCR:$SCR \
  --overlay $CADC_OVERLAY \
  $IMG \
  bash -lc "
    export HOME=$SCR/fake_home_cadc
    export XDG_CACHE_HOME=$SCR/pip_cache_cadc
    export PIP_CACHE_DIR=$SCR/pip_cache_cadc
    export CONDA_PKGS_DIRS=$SCR/conda_pkgs_cadc
    export PYTHONNOUSERSITE=1

    source /opt/conda/etc/profile.d/conda.sh

    mamba create -y -p /cadcenv python=3.11 pip
    conda run -p /cadcenv pip install --no-cache-dir cadcutils vos
  "
```

------

# Part B — Create the shared wrapper commands (For Rongjun only)

All of the following scripts go into:

```text
/software/projects/pawsey1308/bin
```

These are the commands that everyone in the project will use.

------

## 1. Shared `conda1308`

Create:

```bash
cat > /software/projects/pawsey1308/bin/conda1308 <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

module load singularity/4.1.0-slurm >/dev/null 2>&1 || true

BASE=/software/projects/pawsey1308
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ngistenv1308_overlay.img

SCR=/scratch/pawsey1308/${USER}
FAKE_HOME=$SCR/fake_home_1308
CACHE=$SCR/conda_cache_1308
TMP=$SCR/singularity_tmp_1308
PKGS=$SCR/conda_pkgs_1308

mkdir -p "$FAKE_HOME" "$CACHE" "$TMP" "$PKGS"

export SINGULARITY_TMPDIR="$TMP"
export SINGULARITY_CACHEDIR="$CACHE"

exec singularity exec --cleanenv \
  --bind "$BASE:/work" \
  --bind "$SCR:$SCR" \
  --overlay "$OVERLAY" \
  "$IMG" \
  bash -lc "
    export HOME=$FAKE_HOME
    export XDG_CACHE_HOME=$CACHE
    export PIP_CACHE_DIR=$CACHE
    export CONDA_PKGS_DIRS=$PKGS
    export PYTHONNOUSERSITE=1
    source /opt/conda/etc/profile.d/conda.sh
    conda \"\$@\"
  " bash "$@"
EOF

chmod 755 /software/projects/pawsey1308/bin/conda1308
```

------

## 2. Shared `mamba1308`

```bash
cat > /software/projects/pawsey1308/bin/mamba1308 <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

module load singularity/4.1.0-slurm >/dev/null 2>&1 || true

BASE=/software/projects/pawsey1308
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ngistenv1308_overlay.img

SCR=/scratch/pawsey1308/${USER}
FAKE_HOME=$SCR/fake_home_1308
CACHE=$SCR/conda_cache_1308
TMP=$SCR/singularity_tmp_1308
PKGS=$SCR/conda_pkgs_1308

mkdir -p "$FAKE_HOME" "$CACHE" "$TMP" "$PKGS"

export SINGULARITY_TMPDIR="$TMP"
export SINGULARITY_CACHEDIR="$CACHE"

exec singularity exec --cleanenv \
  --bind "$BASE:/work" \
  --bind "$SCR:$SCR" \
  --overlay "$OVERLAY" \
  "$IMG" \
  bash -lc "
    export HOME=$FAKE_HOME
    export XDG_CACHE_HOME=$CACHE
    export PIP_CACHE_DIR=$CACHE
    export CONDA_PKGS_DIRS=$PKGS
    export PYTHONNOUSERSITE=1
    source /opt/conda/etc/profile.d/conda.sh
    mamba \"\$@\"
  " bash "$@"
EOF

chmod 755 /software/projects/pawsey1308/bin/mamba1308
```

------

## 3. Shared `ngistenv1308`

```bash
cat > /software/projects/pawsey1308/bin/ngistenv1308 <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
exec /software/projects/pawsey1308/bin/conda1308 run -p /ngistenv1308 "$@"
EOF

chmod 755 /software/projects/pawsey1308/bin/ngistenv1308
```

Usage example later will be:

```bash
ngistenv1308 python -V
ngistenv1308 ngistPipeline --help
```

------

## 4. Shared `mapviewer1308`

```bash
cat > /software/projects/pawsey1308/bin/mapviewer1308 <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

module load singularity/4.1.0-slurm >/dev/null 2>&1 || true

BASE=/software/projects/pawsey1308
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ngistenv1308_overlay.img

export SINGULARITYENV_DISPLAY="${DISPLAY:-}"
export SINGULARITYENV_XAUTHORITY="$HOME/.Xauthority"
export SINGULARITYENV_QT_X11_NO_MITSHM=1

exec singularity exec --cleanenv \
  --bind /tmp/.X11-unix:/tmp/.X11-unix \
  --bind "$HOME/.Xauthority:$HOME/.Xauthority" \
  --overlay "$OVERLAY" \
  "$IMG" \
  /ngistenv1308/bin/Mapviewer "$@"
EOF

chmod 755 /software/projects/pawsey1308/bin/mapviewer1308
```

------

## 5. Shared `cadc-get-cert`

Each user’s CADC certificate is stored in:

```text
/software/projects/pawsey1308/$USER/cadc_home
```

Create:

```bash
cat > /software/projects/pawsey1308/bin/cadc-get-cert <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

module load singularity/4.1.0-slurm >/dev/null 2>&1 || true

BASE=/software/projects/pawsey1308
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/cadc_overlay.img

SCR=/scratch/pawsey1308/${USER}
CADC_HOME=$BASE/${USER}/cadc_home

mkdir -p "$SCR/singularity_tmp_cadc" "$SCR/singularity_cache_cadc" "$CADC_HOME"
chmod 700 "$CADC_HOME" 2>/dev/null || true

export SINGULARITY_TMPDIR="$SCR/singularity_tmp_cadc"
export SINGULARITY_CACHEDIR="$SCR/singularity_cache_cadc"

exec singularity exec --cleanenv \
  --bind "$BASE:/work" \
  --bind "$SCR:$SCR" \
  --overlay "$OVERLAY" \
  "$IMG" \
  bash -lc '
    export HOME=/work/'"${USER}"'/cadc_home
    /cadcenv/bin/cadc-get-cert "$@"
  ' bash "$@"
EOF

chmod 755 /software/projects/pawsey1308/bin/cadc-get-cert
```

------

## 6. Shared `vcp`

```bash
cat > /software/projects/pawsey1308/bin/vcp <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

module load singularity/4.1.0-slurm >/dev/null 2>&1 || true

BASE=/software/projects/pawsey1308
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/cadc_overlay.img

SCR=/scratch/pawsey1308/${USER}
CADC_HOME=$BASE/${USER}/cadc_home

mkdir -p "$SCR/singularity_tmp_cadc" "$SCR/singularity_cache_cadc" "$CADC_HOME"
chmod 700 "$CADC_HOME" 2>/dev/null || true

export SINGULARITY_TMPDIR="$SCR/singularity_tmp_cadc"
export SINGULARITY_CACHEDIR="$SCR/singularity_cache_cadc"

exec singularity exec --cleanenv \
  --bind "$BASE:/work" \
  --bind "$SCR:$SCR" \
  --overlay "$OVERLAY" \
  "$IMG" \
  bash -lc '
    export HOME=/work/'"${USER}"'/cadc_home
    /cadcenv/bin/vcp "$@"
  ' bash "$@"
EOF

chmod 755 /software/projects/pawsey1308/bin/vcp
```

------

## 7. Shared `vls`

```bash
cat > /software/projects/pawsey1308/bin/vls <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

module load singularity/4.1.0-slurm >/dev/null 2>&1 || true

BASE=/software/projects/pawsey1308
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/cadc_overlay.img

SCR=/scratch/pawsey1308/${USER}
CADC_HOME=$BASE/${USER}/cadc_home

mkdir -p "$SCR/singularity_tmp_cadc" "$SCR/singularity_cache_cadc" "$CADC_HOME"
chmod 700 "$CADC_HOME" 2>/dev/null || true

export SINGULARITY_TMPDIR="$SCR/singularity_tmp_cadc"
export SINGULARITY_CACHEDIR="$SCR/singularity_cache_cadc"

exec singularity exec --cleanenv \
  --bind "$BASE:/work" \
  --bind "$SCR:$SCR" \
  --overlay "$OVERLAY" \
  "$IMG" \
  bash -lc '
    export HOME=/work/'"${USER}"'/cadc_home
    /cadcenv/bin/vls "$@"
  ' bash "$@"
EOF

chmod 755 /software/projects/pawsey1308/bin/vls
```

------

# Part C — what every team member needs to do (For Everyone)

This is the user-side setup.

------

## 1. Add the shared project `bin` to `PATH`

Each user adds this line to their own `~/.bashrc`:

```bash
export PATH=/software/projects/pawsey1308/bin:$PATH
```

Then reload:

```bash
source ~/.bashrc
```

That is enough.
Users do **not** need to copy these scripts into `~/bin`.

------

## 2. Check that the shared stack works

Run:

```bash
conda1308 --version
mamba1308 --version
ngistenv1308 python -V
ngistenv1308 ngistPipeline --help | head
```

------

# Part D — how everyone should use the shared stack (For Everyone)

------

## 1. Using nGIST from the command line

### Basic usage

```bash
ngistenv1308 ngistPipeline --config configFiles/MasterConfig.yaml --default-dir configFiles/defaultDir
```

### Running Python inside the shared environment

```bash
ngistenv1308 python -V
ngistenv1308 python -c "import ngistPipeline; print('ok')"
```

------

## 2. Using Mapviewer

Users must connect with X forwarding:

```bash
ssh -Y <username>@setonix-03.pawsey.org.au
```

Then:

```bash
mapviewer1308
```

If they are on macOS, XQuartz must be running locally.

------

## 3. Using CADC tools

### First-time setup for each user

Each user must obtain their own CADC proxy certificate:

```bash
cadc-get-cert -u <their_CADC_username>
```

This stores the credential under:

```text
/software/projects/pawsey1308/$USER/cadc_home/.ssl/cadcproxy.pem
```

### List ARC files

```bash
vls arc:projects/mauve/cubes/v3tk
```

### Copy an ARC file to scratch

```bash
mkdir -p /scratch/pawsey1308/$USER/mauve/cubes/v3tk

vcp \
  arc:projects/mauve/cubes/v3tk/IC3392_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits.gz \
  /scratch/pawsey1308/$USER/mauve/cubes/v3tk/
```

------

## 4. Using the shared nGIST stack in Slurm jobs

Use the shared wrapper directly.

Example job script under `/software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial` using the shared environment:

```bash
#!/bin/bash -l
#SBATCH --account=pawsey1308
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --job-name=IC3392_v3tk_v7.6.8
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --mail-user=rongjun.huang@research.uwa.edu.au,astro@rongjun-huang.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module load singularity/4.1.0-slurm

cd /software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial

# (Optional but usually wise) prevent oversubscription from BLAS/OpenMP libs:
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# run + stream output to both screen (slurm out) and logfile
srun -n 1 -c $SLURM_CPUS_PER_TASK stdbuf -oL -eL \
  ngistenv1308 ngistPipeline \
  --config configFiles/IC3392_MAUVE_MasterConfig_v7.6.8_setonix.yaml \
  --default-dir configFiles/defaultDir \
  2>&1 | tee IC3392_v3tk_v7.6.8.log

exit ${PIPESTATUS[0]}


```

------

# Part E — update procedure for the maintainer (For Rongjun only)

When the shared stack needs updating:

------

## 1. Update the repos

```bash
cd /software/projects/pawsey1308/ngist
git fetch --all --tags
git checkout main
git pull --ff-only

cd /software/projects/pawsey1308/ngist_supplementary_public
git fetch --all --tags
git checkout main
git pull --ff-only
```

------

## 2. Update the shared environment

```bash
module load singularity/4.1.0-slurm

BASE=/software/projects/pawsey1308
SCR=/scratch/pawsey1308/$USER
IMG=$BASE/containers/mambaforge_latest.sif
OVERLAY=$BASE/containers/ngistenv1308_overlay.img

mkdir -p $SCR/fake_home_1308
mkdir -p $SCR/pip_cache_1308
mkdir -p $SCR/conda_pkgs_1308

singularity exec --cleanenv \
  --bind $BASE:/work \
  --bind $SCR:$SCR \
  --overlay $OVERLAY \
  $IMG \
  bash -lc "
    export HOME=$SCR/fake_home_1308
    export XDG_CACHE_HOME=$SCR/pip_cache_1308
    export PIP_CACHE_DIR=$SCR/pip_cache_1308
    export CONDA_PKGS_DIRS=$SCR/conda_pkgs_1308
    export PYTHONNOUSERSITE=1
    source /opt/conda/etc/profile.d/conda.sh

    mamba env update -p /ngistenv1308 -f /work/ngist_supplementary_public/ngistenv.yaml
    /ngistenv1308/bin/python -m pip install -e /work/ngist --no-deps
  "
```

### Re-test

```bash
ngistenv1308 python -m pip show ngistPipeline | egrep 'Name|Version|Location'
ngistenv1308 ngistPipeline --help | head
```

------

# Final operational summary

## Shared/common

- `/software/projects/pawsey1308/bin`
- `/software/projects/pawsey1308/containers`
- `/software/projects/pawsey1308/ngist`
- `/software/projects/pawsey1308/ngist_supplementary_public`
- `/software/projects/pawsey1308/docs`
- `/software/projects/pawsey1308/setonix`

## Per-user private

- `/software/projects/pawsey1308/$USER/cadc_home`

## Command names for the team

- `conda1308`
- `mamba1308`
- `ngistenv1308`
- `mapviewer1308`
- `cadc-get-cert`
- `vcp`
- `vls`

