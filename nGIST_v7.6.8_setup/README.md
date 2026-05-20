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
