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
  27_galaxies.sh                               # shared cube-ID list and argument filtering
  27_creation.sh                               # create all/selected YAMLs and slurm scripts
  27_send.sh                                   # send all/selected generated files to Setonix
  27_setonix.sh                                # submit all/selected slurm scripts on Setonix
  27_status.sh                                 # check all/selected completion/running/timeout status
  vcp_from_v3tk_to_scratch.sh                  # stage local v3tk and PHANGS native public cubes to Setonix scratch
  QC_ngist_v3tk_v768.py                        # post-run QC PDF generator
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
- PHANGS-native AO masks: `specMask_KIN_narrow10_AO` and
  `specMask_GAS_narrow10_AO` for `NGC4254`, `NGC4321`, and `NGC4535`
- PHANGS-native AO read-data method: `MUSE_WFMAON`
- Gas emission-line config: `emissionLines_ppxf_8900.config`
- Gas fitting level: `BOTH`
- Setonix cube path: `/scratch/pawsey1308/mauve/cubes/v3tk/`
- Setonix product path: `/scratch/pawsey1308/mauve/products/v3tk_v7.6.8/`
- PHANGS-native public cube source:
  `vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/`

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

For local validation only, the generator also accepts `MAUVE_CUBE_DIR` and
`MAUVE_PRODUCTS_DIR` environment overrides. Production Setonix configs should
use the default scratch/product paths unless there is a deliberate run-location
change.

## Staging input cubes to Setonix scratch

The generated YAML files expect input cubes under:

```text
/scratch/pawsey1308/mauve/cubes/v3tk/
```

Use `config_setup/vcp_from_v3tk_to_scratch.sh` on Setonix/CADC-capable
environments to stage the current cube inputs. With no galaxy arguments, it
builds a manifest from the existing MAUVE v3tk cube release and appends these
public PHANGS-MUSE native cubes:

```text
vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4254_PHANGS_DATACUBE_native.fits
vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4321_PHANGS_DATACUBE_native.fits
vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4535_PHANGS_DATACUBE_native.fits
```

These are full 18-42 GB cubes, so production nGIST should use the staged local
scratch files. VOSspace is still useful for lightweight checks such as `vls`,
`vcat --head`, or small `vcp` cutouts, but the full nGIST workflow should not
treat `vos:` paths as normal local FITS inputs.

To transfer only selected galaxies, pass the GALID list. For example, this
stages only the public PHANGS-native `NGC4254` cube:

```bash
./vcp_from_v3tk_to_scratch.sh NGC4254
```

Preview the selected sources and worker count without contacting CADC or
copying files:

```bash
./vcp_from_v3tk_to_scratch.sh --dry-run NGC4254
```

The transfer script defaults to `JOBS=5` and caps any larger exported `JOBS`
value back to 5 for Setonix overlay-quota safety. It also shrinks the effective
worker count to the number of selected sources, so a one-galaxy transfer creates
one worker overlay, not five.

## Creating the Setonix batch

The `27_*.sh` script family keeps its historical name, but the current selectable
MAUVE batch uses 30 cube IDs because `NGC4567` and `NGC4568` are in the combined
cube `NGC4567_8` and the three PHANGS-native galaxies are included.

From inside `config_setup`, run:

```bash
./27_creation.sh
```

With no galaxy arguments, this processes all 30 cube IDs. To regenerate only
selected cube IDs, pass them as positional arguments:

```bash
./27_creation.sh NGC4383 NGC4419
```

The script rejects unknown galaxy IDs before generating files.

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

The default `work` script requests:

```text
#SBATCH --partition=work
#SBATCH --cpus-per-task=128
#SBATCH --mem=230G
#SBATCH --time=24:00:00
```

Large, multi-pointing galaxies are special-cased for Setonix `highmem` instead
of `work`:

```text
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=128
#SBATCH --mem=980G
#SBATCH --time=96:00:00
```

The current `highmem` queue list is:

```text
NGC4192
NGC4254
NGC4321
NGC4501
NGC4535
NGC4569
```

Several galaxies reached the gas module more than once without completing it
within the normal 24-hour `work` walltime. Their generated slurm scripts are
special-cased for Setonix `long`:

```text
#SBATCH --partition=long
#SBATCH --cpus-per-task=128
#SBATCH --mem=230G
#SBATCH --time=96:00:00
```

The current `long` queue list is:

```text
NGC4293
NGC4298
NGC4302
NGC4330
NGC4383
NGC4396
NGC4419
NGC4457
NGC4567_8
NGC4698
```

`NGC4580` is intentionally not in this list: its second run completed the gas
module and advanced to the SFH module, so the next restart should skip gas.

The 30 cube IDs are:

```text
IC3392
NGC4064
NGC4189
NGC4192
NGC4254
NGC4293
NGC4294
NGC4298
NGC4302
NGC4321
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
NGC4535
NGC4567_8
NGC4569
NGC4580
NGC4606
NGC4607
NGC4694
NGC4698
```

`make_gist_config_try.py` still needs a pixel origin for each generated YAML.
The origin normally comes from `cube_centers_v3tk.csv`. If the selected galaxy
has no center row and its cube already exists under the configured cube path,
the script falls back to the cube midpoint; otherwise pass `-center x,y` or add
a verified center row before generating its YAML.

On Setonix, regenerate the cube-center and cube-size CSV snapshots from the
staged scratch cube directory with:

```bash
./check_v3tk_cube_centers.sh
./check_v3tk_cube_sizes.sh
```

Both helpers default to `/scratch/pawsey1308/mauve/cubes/v3tk` and now scan
`*_DATACUBE*.fits,*_DATACUBE*.fits.gz`. This covers ordinary MAUVE v3tk cubes
and staged PHANGS-native cubes such as
`NGC4254_PHANGS_DATACUBE_native.fits`. PHANGS-native rows are written with the
MAUVE galaxy IDs (`NGC4254`, `NGC4321`, `NGC4535`), not `NGC4254_PHANGS`.
Each helper creates a temporary per-run copy of the nGIST container overlay
before calling `${HOME}/bin/conda`, then removes it on exit. This avoids the
shared `ngist_overlay.img` write lock collision when another nGIST container
process is already using the base overlay.

To send the generated files to Setonix from your local machine, run:

```bash
./27_send.sh rhuang
```

By default, this connects to `rhuang@setonix.pawsey.org.au`. You can override
the host if needed:

```bash
./27_send.sh rhuang setonix.pawsey.org.au
```

To send only selected cube IDs, pass the galaxy IDs after the username. If the
second argument is a known galaxy ID, the script keeps the default Setonix host:

```bash
./27_send.sh rhuang NGC4383 NGC4419
```

For the default `rhuang@setonix.pawsey.org.au` login, you can also pass galaxy
IDs directly:

```bash
./27_send.sh NGC4569
```

With an explicit host, put the host before the galaxy IDs:

```bash
./27_send.sh rhuang setonix.pawsey.org.au NGC4383 NGC4419
```

The send script copies:

- The selected generated YAML files, all 30 by default, to
  `/software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial/configFiles/`
- The selected generated slurm scripts, all 30 by default, plus
  `27_galaxies.sh`, `27_setonix.sh`, and `27_status.sh` to
  `/software/projects/pawsey1308/ngist_supplementary_public/ngistTutorial/`

It uses `tar` streamed through `ssh`, not `scp`, because Setonix prints a login
notice that can break the `scp`/SFTP protocol in non-interactive sessions. The
script also disables macOS AppleDouble metadata during tar creation and removes
any matching `._*` sidecar files from the two target directories.

On Setonix, submit all jobs from the tutorial directory with:

```bash
./27_setonix.sh
```

To submit only selected cube IDs:

```bash
./27_setonix.sh NGC4383 NGC4419
```

That script sequentially runs:

```bash
sbatch {GALID}_v3tk_v7.6.8_setonix.slurm
```

for the selected cube IDs, or all 30 cube IDs by default. The `sbatch` commands
are issued one by one, but the jobs can then run together according to the
Setonix scheduler.

To check job completion and timeout status on Setonix, run from the same
tutorial directory:

```bash
./27_status.sh
```

To check only selected cube IDs:

```bash
./27_status.sh NGC4383 NGC4419
```

The status script reads each product `LOGFILE` under:

```text
/scratch/pawsey1308/mauve/products/v3tk_v7.6.8/{GALID}/LOGFILE
```

Finished jobs are identified by:

```text
MainPipeline: nGIST completed successfully.
```

For unfinished jobs, it checks the run log in the tutorial directory, for
example:

```text
NGC4383_v3tk_v7.6.8.log
```

If the run log contains `DUE TO TIME LIMIT`, the status is reported as
`TIMEOUT_RESBATCH` and the script prints the `sbatch` command to resubmit that
galaxy. If the job is still visible in `squeue`, it is reported as `RUNNING`.
If the run log exists but is still empty, it is reported as `RUNNING_EMPTY_LOG`.

The report is printed to screen and saved in the current directory as:

```text
27_status_log_YYYYmmdd_HHMMSS.txt
```

For unfinished jobs, the script also gives a conservative rough remaining-time
estimate using the product `CONFIG` beside the `LOGFILE`. It checks which
modules are enabled and which enabled modules have not completed yet, then
sums the estimated remaining module times.

For gas fitting, this setup runs `GAS.LEVEL: BOTH`, so the full gas workload is:

```text
GAS_WORK = SPECTRA + BINS
```

Here `SPECTRA` is the cube-read count from `Read a total of ... spectra`, used
as the spaxel-level workload proxy, and `BINS` is the Voronoi-bin count.
If BIN-level gas has already completed, the gas estimate switches to the
SPAXEL-only resume case and scales by `SPECTRA`. SFH estimates scale by `BINS`.
If LS or user modules are enabled but no estimator is implemented, they are
listed as `NA` in the remaining-module breakdown.

Each module estimate uses the maximum scaled module time from comparable
finished jobs, not the median. These estimates are approximate because nGIST
can skip completed modules after restart, but it does not checkpoint every
internal step within a module.

The status script also prints `long` queue warnings when either:

- the same module has restarted more than once without completing, which means
  another 24-hour `work` job may repeat the same failure mode;
- at least one remaining module estimate is longer than 22 hours;
- the summed `EST_REMAIN` is longer than 22 hours, leaving too little margin for
  the 24-hour `work` walltime.

The warning text reads the local `{GALID}_v3tk_v7.6.8_setonix.slurm` file. If
the job is still on `work`, it recommends switching to `long`; if it is already
on `long` or `highmem`, it warns not to resubmit that galaxy on `work`.

## QC PDF script

The QC plotting script is:

```text
config_setup/QC_ngist_v3tk_v768.py
```

It is a post-run diagnostic script for the `v3tk_v7.6.8` products. It reads the
generated nGIST maps for one galaxy and writes a multi-page PDF named:

```text
{GALID}_v3tk_v7.6.8_QC.pdf
```

The script keeps `IC3392` as the visible default galaxy:

```python
galaxyid = 'IC3392'
```

but also supports a command-line override:

```bash
cd config_setup
python QC_ngist_v3tk_v768.py NGC4419
```

Input products expected for each galaxy are:

```text
{GALID}/{GALID}_gas_bin_maps.fits
{GALID}/{GALID}_gas_spaxel_maps.fits
{GALID}/{GALID}_kin_maps.fits
{GALID}/{GALID}_sfh_maps.fits
```

The current local data path in the script is:

```text
/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8
```

The corresponding Pawsey product path is noted in the script as:

```text
/scratch/pawsey1308/mauve/products/v3tk_v7.6.8
```

Change `data_folder` before running on Setonix or another filesystem.

The PDF keeps the original five-subplot page structure. For each of the main
gas lines, it plots both BIN-level and SPAXEL-level maps:

```text
HB4861
OIII5006
OI6300
OI6363
HA6562
NII6583
SII6716
```

For each line it shows:

```text
FLUX
FLUX_ERR
VEL
SIGMA
SIGMA_ERR
```

The script then adds derived gas diagnostic pages, again for both BIN and
SPAXEL levels where available:

```text
Ha/Hb
NII/Ha
SII/NII
Ha-Hb vel
Ha/Hb sigma
```

It also plots BIN-level stellar kinematics from `{GALID}_kin_maps.fits`:

```text
V
SIGMA
H3
H4
FORM_ERR_SIGMA
```

and gas-vs-stellar or gas-vs-gas consistency checks:

```text
Vs-Vha
Sigma_s - Sigma_ha
SII6716/SII6730
VHa-VNII
SHa-SNII
```

For SPAXEL-level gas-vs-stellar comparisons, the script uses the stellar
reference maps stored inside the gas SPAXEL product:

```text
V_STARS2
SIGMA_STARS2
```

The final maps come from `{GALID}_sfh_maps.fits`:

```text
AGE
METAL
EBV
```

Gas velocity, sigma, and sigma-error color limits are shared between BIN and
SPAXEL maps for the same quantity. This makes the two levels easier to compare
directly. Some ratio and residual maps keep fixed limits, for example H-alpha
minus H-beta velocity, H-alpha over H-beta sigma, SII6716/SII6730, and H-alpha
minus NII kinematics.

The script was kept intentionally close to the older `QC_ngist.py` workflow:
the page structure, repeated `AGE`/`METAL` style, and five-map layout are
preserved so the new v7.6.8 outputs can be checked against the familiar QC
format with minimal changes.

The QC script requires:

```bash
python
astropy
numpy
matplotlib
```

## Why `make_gist_config_try.py`

`make_gist_config_try.py` is the script used for this setup. It adds two
important conveniences compared with the older `make_gist_config.py`:

- It can resolve `GIST_setupinput_v1.fits` when that file is a text pointer to a
  FITS file under `input_tables/`.
- If `-center` is not supplied, it reads the cube center from
  `cube_centers_v3tk.csv`, then falls back to a local staged cube midpoint when
  the CSV has no row for the selected galaxy.
- It handles the combined `NGC4567_8` cube, whose cube-center row is combined
  but whose setup-table rows are still listed separately as `NGC4567` and
  `NGC4568`.
- It uses the PHANGS-native public cube filenames for `NGC4254`, `NGC4321`, and
  `NGC4535`, for example
  `/scratch/pawsey1308/mauve/cubes/v3tk/NGC4254_PHANGS_DATACUBE_native.fits`.
- It switches those three PHANGS-native AO-enabled galaxies to the
  `MUSE_WFMAON` read-data method, `specMask_KIN_narrow10_AO` for KIN/CONT/SFH,
  and `specMask_GAS_narrow10_AO` for GAS. The AO read-data method is needed
  before spatial masking/binning, while these AO masks copy the narrow10 masks
  and add the observed-frame global-NaN range `5806.25-5963.75 Angstrom` as a
  `sky` mask entry.

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
- The input cube exists under `/scratch/pawsey1308/mauve/cubes/v3tk/`; for
  `NGC4254`, `NGC4321`, and `NGC4535`, this should be the staged
  `*_PHANGS_DATACUBE_native.fits` file.
- PHANGS-native AO-enabled configs use `READ_DATA: METHOD: MUSE_WFMAON` and
  reference the `_AO` spectral masks, which include the observed-frame
  `5806.25-5963.75 Angstrom` global-NaN mask.
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
