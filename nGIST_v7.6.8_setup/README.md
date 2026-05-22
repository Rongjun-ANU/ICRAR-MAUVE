# nGIST v7.6.8 MAUVE Setonix setup

> WARNING: use `make_gist_config_try.py` to generate the current YAML files.
> The older `make_gist_config.py` is kept for reference/provenance and should
> not be the default script for this setup.

This repository contains the MAUVE configuration setup files used to generate
nGIST/GIST YAML configuration files for the v7.6.8 Setonix runs.

The repository is intended to keep the master configuration, target input
tables, generated galaxy YAML files, and older setup references together in one
version-controlled place.

## Repository layout

```text
config_setup/
  MAUVE_MasterConfig_v7.6.8_setonix.yaml       # master nGIST/GIST config
  make_gist_config_try.py                      # actual generator to run
  make_gist_config.py                          # older/reference generator
  v3tk_v7.6.8_setonix.slurm                    # slurm template with GALID placeholder
  27_creation.sh                               # create 26 YAMLs and 26 slurm scripts
  27_send.sh                                   # send generated files to Setonix
  27_setonix.sh                                # submit the 26 slurm scripts on Setonix
  GIST_setupinput_v1.fits                      # pointer to setup FITS table
  cube_centers_v3tk.csv                        # cube centers for v3tk cubes
  *_MAUVE_MasterConfig_v7.6.8_setonix.yaml     # generated galaxy configs
  input_tables/                                # MAUVE input tables and redshifts
  old/                                         # older config/script versions

config_setup.zip                               # archived copy of setup folder
```

## Main configuration

The active master config is:

```text
config_setup/MAUVE_MasterConfig_v7.6.8_setonix.yaml
```

Important current settings include:

- Wavelength range: `4800-8900 Angstrom`
- Voronoi target S/N: `40`
- Velocity scale: `41 km/s`
- Stellar/SFH mask: `specMask_KIN_narrow10`
- Gas mask: `specMask_GAS_narrow10`
- Gas emission-line config: `emissionLines_ppxf_8900.config`
- Gas fitting level: `BOTH`
- Setonix cube path: `/scratch/pawsey1308/mauve/cubes/v3tk/`
- Setonix product path: `/scratch/pawsey1308/mauve/products/v3tk_v7.6.8/`

## Creating a galaxy YAML file

Run `make_gist_config_try.py` from inside `config_setup`:

```bash
cd config_setup
python make_gist_config_try.py IC3392
```

This creates:

```text
IC3392_MAUVE_MasterConfig_v7.6.8_setonix.yaml
```

The script fills galaxy-specific values from the setup input table:

- `GENERAL.RUN_ID`
- `GENERAL.INPUT`
- `GENERAL.OUTPUT`
- `GENERAL.REDSHIFT`
- `GENERAL.NCPU`
- `READ_DATA.ORIGIN`
- `READ_DATA.EBmV`
- `SPATIAL_MASKING.MASK`
- `KIN.SIGMA`
- `CONT.SIGMA`

Manual center override:

```bash
python make_gist_config_try.py IC3392 -cpu 128 -center 219,219
```

## Creating the 27-galaxy Setonix batch

The 27-galaxy MAUVE batch uses 26 cube IDs because `NGC4567` and `NGC4568`
are in the combined cube `NGC4567_8`.

From inside `config_setup`, run:

```bash
./27_creation.sh
```

This runs the active generator once per cube ID, for example:

```bash
./make_gist_config_try.py NGC4567_8
```

For each cube ID, it creates:

```text
{GALID}_MAUVE_MasterConfig_v7.6.8_setonix.yaml
{GALID}_v3tk_v7.6.8_setonix.slurm
```

The slurm scripts are copied from `v3tk_v7.6.8_setonix.slurm` by replacing the
literal `GALID` placeholder with the current cube ID.

`NGC4192` and `NGC4501` are large, multi-pointing galaxies. Their generated
slurm scripts are special-cased for Setonix `highmem` instead of `work`:

```text
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=128
#SBATCH --mem=980G
#SBATCH --time=96:00:00
```

The 26 cube IDs are:

```text
IC3392
NGC4064
NGC4189
NGC4192
NGC4293
NGC4294
NGC4298
NGC4302
NGC4330
NGC4351
NGC4383
NGC4388
NGC4394
NGC4396
NGC4402
NGC4405
NGC4419
NGC4457
NGC4501
NGC4522
NGC4567_8
NGC4580
NGC4606
NGC4607
NGC4694
NGC4698
```

To send the generated files to Setonix from your local machine, run:

```bash
./27_send.sh rhuang
```

By default, this connects to `rhuang@setonix.pawsey.org.au`. You can override
the host if needed:

```bash
./27_send.sh rhuang setonix.pawsey.org.au
```

The send script copies:

- The 26 generated YAML files to
  `/software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial/configFiles/`
- The 26 generated slurm scripts plus `27_setonix.sh` to
  `/software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial/`

It uses `tar` streamed through `ssh`, not `scp`, because Setonix prints a login
notice that can break the `scp`/SFTP protocol in non-interactive sessions. The
script also disables macOS AppleDouble metadata during tar creation and removes
any matching `._*` sidecar files from the two target directories.

On Setonix, submit all jobs from the tutorial directory with:

```bash
./27_setonix.sh
```

That script sequentially runs:

```bash
sbatch {GALID}_v3tk_v7.6.8_setonix.slurm
```

for all 26 cube IDs. The `sbatch` commands are issued one by one, but the jobs
can then run together according to the Setonix scheduler.

## Why `make_gist_config_try.py`

`make_gist_config_try.py` is the script used for this setup. It adds two
important conveniences compared with the older `make_gist_config.py`:

- It can resolve `GIST_setupinput_v1.fits` when that file is a text pointer to a
  FITS file under `input_tables/`.
- If `-center` is not supplied, it reads the cube center from
  `cube_centers_v3tk.csv`.
- It handles the combined `NGC4567_8` cube, whose cube-center row is combined
  but whose setup-table rows are still listed separately as `NGC4567` and
  `NGC4568`.

The original `make_gist_config.py` remains in the repository as a reference to
the older workflow.

For `NGC4567_8`, the generated config keeps `NGC4567_8` as the `RUN_ID`,
input cube name, and mask name. The script uses the mean `z`, `EBmV`, and
initial `SIGMA` from the separate `NGC4567` and `NGC4568` setup-table rows.

## Current generated configs

The repository currently includes generated v7.6.8 Setonix configs for:

- `IC3392`
- `NGC4383`
- `NGC4396`
- `NGC4419`
- `NGC4501`
- `NGC4567_8`
- `NGC4698`

## Required Python packages

The generator scripts require:

```bash
python
astropy
numpy
pyyaml
```

On an existing science environment, the usual check is:

```bash
python -c "import astropy, numpy, yaml"
```

## Notes before running nGIST/GIST

This repository does not contain the full MUSE datacubes, masks, nGIST
`configFiles`, stellar templates, or output products. The generated YAML files
expect those resources to exist on Setonix at the configured paths.

Before launching a production run, check:

- The generated galaxy YAML has the intended `REDSHIFT`, `EBmV`, `SIGMA`, and
  `ORIGIN`.
- The input cube exists under `/scratch/pawsey1308/mauve/cubes/v3tk/`.
- The galaxy mask file, for example `IC3392_mask.fits`, is available where
  nGIST/GIST expects it.
- The nGIST `configFiles` directory contains the referenced masks, LSF files,
  emission-line config, and template definitions.
- The required stellar population templates are available under
  `spectralTemplates`.

## Version-control notes

The authoritative editable setup files are in `config_setup/`. Older files are
kept under `config_setup/old/` for provenance.

When changing the master config or setup tables, make a new commit that clearly
states what changed, for example:

```bash
git add README.md config_setup/
git commit -m "Document nGIST v7.6.8 Setonix setup"
```
