# SFH16 7000-run mass and QC batch

Run from the `further` directory on CANFAR, with the science Python environment
active. The wrapper covers the same 39 cube IDs as the current Setonix
`config_setup/27_galaxies.sh` (40 galaxies; `NGC4567_8` is one combined cube).
It includes NGC4579 so it can be processed when its products become available.
It does not depend on the Setonix setup repository being installed on CANFAR.

```bash
# Print the selection and routing without writing or running anything.
bash QC_ngist_v3tk_v768_7000_SFH16.sh --dry-run

# Process all available SFH16 runs sequentially.
bash QC_ngist_v3tk_v768_7000_SFH16.sh

# Process one or several explicitly selected runs; either ID spelling works.
bash QC_ngist_v3tk_v768_7000_SFH16.sh NGC4567_8
bash QC_ngist_v3tk_v768_7000_SFH16.sh NGC4189_SFH16 NGC4380

# Explicit paths, including an alternative location for PDF output.
bash QC_ngist_v3tk_v768_7000_SFH16.sh \
  --root /arc/projects/mauve/products/v3tk_v7.6.8_7000 \
  --cube-root /arc/projects/mauve/cubes/v3tk \
  --output-dir ./QC NGC4567_8
```

The default product root is `/arc/projects/mauve/products/v3tk_v7.6.8_7000`,
falling back to `$PWD/v3tk_v7.6.8_7000` only when the former directory is absent.
An explicit `--root` selects the parent containing both `G/` and `G_SFH16/`.
The wrapper routes both Mass and QC through one private input view, including
on CANFAR where the original QC normally prefers the `/arc` root. An alternative
explicit root therefore works even when `/arc` is present.
Python is selected from `PYTHON_BIN`, the active `CONDA_PREFIX`, the local ICRAR
interpreter if present, then `python` on PATH. It must support both existing
Python scripts and their dependencies. The existing photometry and catalogue
tables beside `Mass.py` are still required.

For each physical cube ID `G`, the wrapper uses the four-file SFH16 transfer
layout. Only the following inputs come from the new SFH16 run:

```text
<root>/G_SFH16/CONFIG
<root>/G_SFH16/LOGFILE
<root>/G_SFH16/G_SFH16_sfh_maps.fits
<root>/G_SFH16/G_SFH16_sfh_weights.fits
```

The unchanged inputs are read from the neighbouring original 7000 directory:

```text
<root>/G/G_spatial_binning_maps.fits
<root>/G/G_kin_maps.fits
<root>/G/G_gas_bin_maps.fits
<root>/G/G_gas_spaxel_maps.fits
```

The original 7000 run must be the one whose spatial bins and fixed kinematics
were retained for that SFH16 fit. The wrapper never uses inherited-file
duplicates in `G_SFH16/`, even if they remain from an older transfer. It leaves
those duplicates in place; removing them is a separate storage-cleanup task.
The mask is not read by either step and need not be copied into `G_SFH16/`.
FITS discovery also accepts `.fits.gz`, `SPATIAL_BINNING_maps`, `SFH_maps`,
`sfh-weights`, `KIN_maps`, and gas `BIN`/`SPAXEL` variants (including `GAS`).

Only run after the SFH16 job and product transfer finish. Presence checks are
not a completion, checksum, FITS-integrity, or scientific-readiness audit.
The wrapper requires the SFH16 CONFIG and LOGFILE but does not parse the log's
completion marker or compare bin assignments. It does not fetch missing files.

The mass step creates a private temporary directory containing a combined input
view. It invokes `Mass.py -g G --root <view> --product-subdir .
--fallback-root <view> --ngist-config <root>/G_SFH16/CONFIG
--cube-root <cube-root>`, preserving
the physical galaxy ID for cube lookup, distance, inclination, and catalogue
matching. Both product search roots point to the private view, and the exact
SFH16 CONFIG is supplied explicitly. Existing Mass checks of that CONFIG and
the SFH weights remain active.

The binning input has an SFH16 basename in the temporary view so that Mass's
shared-keyword discovery and output naming work. It is hardlinked when possible;
otherwise only that binning map is copied temporarily. The SFH, KIN and gas
inputs use symlinks to their selected source files. No persistent inherited
FITS copies are created. Temporary inputs are removed on exit, including after
Mass or QC failure.

Mass must exit successfully and create a nonempty, fresh output in the private
view before the wrapper moves it to:

```text
<root>/G_SFH16/G_SFH16_spatial_binning_maps_further.fits
```

A previous mass product cannot satisfy this check. The ordinary-run mass
product is never selected as an output.

QC then combines the new SFH and mass maps with the original KIN and gas maps.
The shell wrapper's embedded Python launcher parses the existing
`QC_ngist_v3tk_v768_7000.py`, replaces only its top-level `data_folder` assignment
in memory with the private view, and executes the existing plotting code with
`G_SFH16` as its sole galaxy argument. This handles the original `/arc` priority
without modifying the QC source, plotting calculations, or original shell file.
If the expected assignment changes in a future QC version, the launcher fails
explicitly rather than selecting another product tree. Other outputs are:

```text
<output-dir>/G_SFH16_v3tk_v7.6.8_7000_QC.pdf
<output-dir>/logs_SFH16/G_SFH16_mass.log
<output-dir>/logs_SFH16/G_SFH16_qc.log
```

The default output directory is `$PWD/QC`. Rerunning regenerates the SFH16 mass
product and QC PDF and replaces those per-run logs. The ordinary-run files
have different names and are not selected. Avoid concurrent invocations on the
same run because they would write the same products and logs.

Without galaxy arguments, missing product sets are reported and skipped, while
available runs proceed. With explicit galaxy arguments, missing product sets
count as failures. A Mass failure prevents QC for that run; handled Mass/QC
failures do not stop the remaining selected runs. Unexpected filesystem errors
while creating the view or publishing outputs abort the batch. The summary
reports passed, skipped, and failed counts. Exit status is nonzero if any run
fails or no run succeeds. A zero exit
with skips means the available subset succeeded, not the complete sample.
`--dry-run` prints plans only; it does not establish input readiness.

Optional environment variables:

- `PYTHON_BIN`: Python executable.
- `CUBE_ROOT`: default cube directory, overridden by `--cube-root`.
- `PHANGS_CUBE_ROOT`: first location checked for the PHANGS-native cubes for
  NGC4254, NGC4321, and NGC4535, followed by the selected cube root and staging
  root. Existing nonempty `.fits` or `.fits.gz` cubes are reused and retained.
- `PHANGS_STAGING_ROOT`: temporary download parent (default `/scratch`). Missing
  native cubes are downloaded with `vcp` into a private directory here, passed
  to Mass with `--phangs-cube-root`, and removed immediately after Mass succeeds
  or fails. Exit cleanup also removes this private directory. Only one cube is
  staged at a time. The directory must have room for the native cube.
- `PHANGS_NATIVE_VOS_DIR`: download source, default
  `vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES`, matching ordinary `Mass.sh`.
  Downloads require `vcp` and working CANFAR/VOS access. A download failure
  prevents Mass/QC for that galaxy, removes partial temporary files, and allows
  the batch to continue.
- `TEMPLATE_LIBRARY_DIR`: passed as `--template-library-dir` for Mass's
  template-based M/L conversion. Existing Mass validity checks remain active.

Local verification:

```bash
bash -n QC_ngist_v3tk_v768_7000_SFH16.sh
/opt/miniconda3/envs/ICRAR/bin/python tests/test_qc_sfh16_wrapper.py
```

After deploying the updated shell wrapper to CANFAR, retry only the three
missing-cube failures:

```bash
bash QC_ngist_v3tk_v768_7000_SFH16.sh NGC4254 NGC4321 NGC4535
```

The wrapper tests use temporary product placeholders, recording Mass/QC
stand-ins, and the actual Mass filename-discovery/output-naming helpers.
They check split-source routing, the combined ID, preserved ordinary products
and old duplicates, compressed/uppercase filenames, temporary-copy fallback,
fresh-output gating, failure handling, and cleanup. They verify routing and
control flow, not stellar masses or real QC map values. A production CANFAR
Mass/QC run is still needed to validate real data, bin compatibility, available
templates/cubes, and runtime dependencies.
