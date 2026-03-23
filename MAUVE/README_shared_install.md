# README_shared_install

This document explains how the shared nGIST stack for project `pawsey1308` is organised, how it was installed, and what each team member needs to do to use it.

The **main purpose** of this shared setup is to provide a common, project-wide **nGIST environment** and some convenient **wrapper commands** that everyone in the project can use. The **CADC / VOS tools** are included as an **optional extra** for users who need to transfer data from CADC / CANFAR.

The shared documentation and example files live in the project docs folder:

- `/software/projects/pawsey1308/docs/README_shared_install.md`
- `/software/projects/pawsey1308/docs/IC3392_MAUVE_MasterConfig_v7.6.8_setonix.yaml`
- `/software/projects/pawsey1308/docs/IC3392_v3tk_v7.6.8_setonix.slurm`

The two IC3392 files above are provided as a concrete example of how to run nGIST on Setonix with the shared setup.

------

## 1. Overview

The shared setup follows this model:

- **Shared project-wide software** lives under:

  ```bash
  /software/projects/pawsey1308
  ```

- **Large working data, caches, logs, and outputs** should live under

  ```bash
  /scratch/pawsey1308/mauve
  ```

- **User-specific personal files** that should not be shared, for example CADC proxy certificates, should live under each user’s **existing personal folder** inside the project software area, for example:

  ```bash
  /software/projects/pawsey1308/$USER/cadc_home
  ```

------

## 2. Directory layout

The shared project tree is organised like this:

```text
/software/projects/pawsey1308/
├── bin/
│ ├── conda1308 # if you like conda
│ ├── mamba1308 # if you like mamba
│ ├── ngistenv1308
│ ├── mapviewer1308
│ ├── cadc-get-cert # optional
│ ├── vcp # optional
│ └── vls # optional
├── containers/
│ ├── mambaforge_latest.sif
│ ├── ngistenv1308_overlay.img
│ └── cadc_overlay.img
├── docs/
│ ├── README_shared_install.md # PLEASE READ ME
│ ├── cp_cubes_from_acacia_to_scratch.sh
│ ├── cp_products_from_acacia_to_scratch.sh
│ ├── cp_cubes_from_scratch_to_acacia.sh
│ ├── cp_products_from_scratch_to_acacia.sh
│ ├── IC3392_MAUVE_MasterConfig_v7.6.8_setonix.yaml # example of yaml file
│ └── IC3392_v3tk_v7.6.8_setonix.slurm # example of slurm file
├── ngist/ # ngist repo
├── ngist_supplementary_public/ # ngist working repo
├── setonix/ 
├── rhuang/
└── ...
```

### What each part is for

- `ngist/`
  Shared nGIST repository.
- `ngist_supplementary_public/`
  Shared supplementary repository, including example configs and `ngistenv.yaml`.
- `containers/mambaforge_latest.sif`
  Shared base container image.
- `containers/ngistenv1308_overlay.img`
  Shared writable overlay that contains the project-wide nGIST environment at `/ngistenv1308`.
- `bin/`
  Shared wrapper commands that all project members can use.
- `<user>/cadc_home`
  Per-user CADC credential store. This is **not shared** even though it lives under the project tree.
- `docs/`
  Shared documentation and working examples for the team.

------

## 3. Design principles

The setup is based on the following rules.

### 3.1 Repositories stay outside the container

The shared repositories:

```bash
/software/projects/pawsey1308/ngist
/software/projects/pawsey1308/ngist_supplementary_public
```

remain normal directories so they can be:

- updated with `git pull`
- switched between branches or tags
- edited directly
- shared easily across the project

### 3.2 Environments live inside overlays

The Conda or Mamba style environment is **not** installed as a plain directory tree on `/software`, because that creates too many small files and exceeds the quota limit.

Instead, the shared environment is installed **inside a Singularity overlay**:

```bash
/software/projects/pawsey1308/containers/ngistenv1308_overlay.img
```

and mounted at runtime as:

```bash
/ngistenv1308
```

inside the container.

### 3.3 Shared commands are provided via wrappers

Users do **not** need to type long `singularity exec ...` commands. Instead, the shared wrappers in:

```bash
/software/projects/pawsey1308/bin
```

provide convenient commands like:

- `conda1308`
- `mamba1308`
- `ngistenv1308`
- `mapviewer1308`

### 3.4 Data and outputs belong on scratch

Large input cubes, outputs, caches, and logs should not be kept under the software tree.

Use:

```bash
/scratch/pawsey1308/mauve
```

for:

- input FITS cubes
- working directories
- logs
- output products
- container caches
- Conda package caches

------

## 4. Shared commands provided

The following commands are available project-wide via:

```bash
/software/projects/pawsey1308/bin
```

### Main nGIST commands

- `conda1308`
- `mamba1308`
- `ngistenv1308`
- `mapviewer1308`

### Optional CADC / ARC / VOS commands

- `cadc-get-cert`
- `vcp`
- `vls`

The **main recommended entry point for normal nGIST usage** is:

```bash
ngistenv1308
```

------

## 5. What every team member should do

### 5.1 Add shared bin to `PATH`

Each user should add this line to their own `~/.bashrc`:

```bash
export PATH=/software/projects/pawsey1308/bin:$PATH
```

Then reload:

```bash
source ~/.bashrc
```

------

### 5.2 Recommended modules in `~/.bashrc`

It is also recommended to add these to `~/.bashrc`:

```bash
module load singularity/4.1.0-slurm
module load rclone/1.68.1
```

### Notes

- `singularity/4.1.0-slurm` is recommended because the shared wrappers and Mapviewer use Singularity.

- `rclone/1.68.1` is recommended because many users may want to move data to or from object storage.

- If either module is unavailable in the future, use:

  ```bash
  module spider singularity
  module spider rclone
  ```

  and adjust the version.

------

### 5.3 Quick checks

After updating `~/.bashrc`, open a new shell or run:

```bash
source ~/.bashrc
```

Then check:

```bash
conda1308 --version
mamba1308 --version
ngistenv1308 python -V
ngistenv1308 ngistPipeline --help | head
```

If these commands work, the shared nGIST stack is available.

------

## 6. Normal nGIST usage

### 6.1 Running nGIST from the command line

Example:

```bash
ngistenv1308 ngistPipeline --config configFiles/MasterConfig.yaml --default-dir configFiles/defaultDir
```

This runs `ngistPipeline` inside the shared environment `/ngistenv1308`.

------

### 6.2 Running arbitrary Python commands inside the shared environment

Example:

```bash
ngistenv1308 python -V
ngistenv1308 python -c "import ngistPipeline; print('ok')"
```

------

### 6.3 Checking the installed nGIST version

```bash
ngistenv1308 python -m pip show ngistPipeline | egrep 'Name|Version|Location'
```

This should returns

```bash
Name: ngistPipeline
Version: 7.6.8
Location: /ngistenv1308/lib/python3.11/site-packages
```

------

### 6.4 Running Mapviewer

To use the GUI Mapviewer, users must connect with X forwarding, for example:

```bash
ssh -Y username@setonix.pawsey.org.au
```

Then run:

```bash
mapviewer1308
```

### Notes	

- On macOS, **XQuartz** must be running locally.
- `mapviewer1308` is just a wrapper around the shared containerised Mapviewer for quick visualization.

------

## 7. Recommended project workflow

A good default workflow is:

### Shared under `/software/projects/pawsey1308`

- repos
- wrappers
- containers
- overlays
- documentation
- helper scripts

Example scratch structure:

```text
/scratch/pawsey1308/mauve/
├── v3tk/
│ ├── GalaxyID_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits
│ └── GalaxyID_mask.fits
├── v3tk_v7.6.8/
│ └── GalaxyID/
```

------

## 8. Example Slurm usage

A typical Slurm script under `/software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial` using the shared environment looks like this.

**Important:** users should update the email notification address before submitting.

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

### Notes

- The concrete example files
  - `IC3392_MAUVE_MasterConfig_v7.6.8_setonix.yaml`
  - `IC3392_v3tk_v7.6.8_setonix.slurm`
    are stored in `/software/projects/pawsey1308/docs/`.
- `stdbuf -oL -eL ... | tee ...` gives live line-buffered logging.
- `exit ${PIPESTATUS[0]}` ensures Slurm sees the correct exit code from `ngistPipeline`.
- `set -euo pipefail` means:
  - `-e`: stop if a command fails
  - `-u`: stop on undefined variables
  - `pipefail`: a pipeline fails if any command in it fails

------

## 9. Updating the shared nGIST stack, maintainer only

This section is for the maintainer of the shared setup.

### 9.1 Update the repositories

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

### 9.2 Update the shared environment

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

# 10. Optional CADC / ARC / VOS tools

This section is optional. It is **not required for normal nGIST usage**.

The shared setup also provides:

- `cadc-get-cert`
- `vcp`
- `vls`

These are intended for users who need to access CANFAR / ARC data.

### Important

Each user must keep their own credentials in their own folder:

```bash
/software/projects/pawsey1308/$USER/cadc_home
```

### First-time setup for a user

Run:

```bash
cadc-get-cert -u <your_CADC_username>
```

This stores the CADC proxy in your own private location.

### Acacia reminder

The project also uses **Acacia** for permanent storage. Before using `rclone`, users must follow the Acacia instructions from Pawsey and from the project to create and save their Acacia access key and secret key correctly.

Once `rclone` is configured, some useful quick inspection commands are:

```bash
rclone lsd pawsey1308:mauve
rclone ls pawsey1308:mauve
rclone lsd pawsey1308:mauve/cubes
rclone lsd pawsey1308:mauve/products
```

The docs folder also provides helper scripts for copying between Acacia and scratch:

```bash
/software/projects/pawsey1308/docs/cp_cubes_from_acacia_to_scratch.sh
/software/projects/pawsey1308/docs/cp_products_from_acacia_to_scratch.sh
/software/projects/pawsey1308/docs/cp_cubes_from_scratch_to_acacia.sh
/software/projects/pawsey1308/docs/cp_products_from_scratch_to_acacia.sh
```

### Example: list ARC files

```bash
vls arc:projects/mauve/cubes/v3tk
```

### Example: copy a file from ARC to scratch

```bash
mkdir -p /scratch/pawsey1308/$USER/mauve/cubes/v3tk

vcp \
  arc:projects/mauve/cubes/v3tk/IC3392_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits.gz \
  /scratch/pawsey1308/$USER/mauve/cubes/v3tk/
```

------

## 11. Summary for users

### One-time setup

Add to `~/.bashrc`:

```bash
export PATH=/software/projects/pawsey1308/bin:$PATH
module load singularity/4.1.0-slurm
module load rclone/1.68.1
```

Then reload:

```bash
source ~/.bashrc
```

### Main commands

```bash
conda1308 --version
mamba1308 --version
ngistenv1308 python -V
ngistenv1308 ngistPipeline --help | head
mapviewer1308
```

### Optional CADC commands

```bash
cadc-get-cert -u <your_CADC_username>
vls arc:projects/mauve/cubes/v3tk
vcp arc:projects/mauve/cubes/v3tk/ /scratch/pawsey1308/$USER/
```

------

## 12. Summary for the maintainer

### Shared or common under `/software/projects/pawsey1308`

- `bin/`
- `containers/`
- `docs/`
- `ngist/`
- `ngist_supplementary_public/`
- example YAML, Slurm files, and Acacia helper scripts in `docs/`

### Per-user personal under `/software/projects/pawsey1308/$USER`

- `cadc_home/`

### Commands for the team

- `conda1308`
- `mamba1308`
- `ngistenv1308`
- `mapviewer1308`
- `cadc-get-cert`
- `vcp`
- `vls`
