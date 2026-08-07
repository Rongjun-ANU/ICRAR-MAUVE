# 2026-05-21 Notes: `v3tk_VRI_image.sh` Pipeline

## Purpose

`v3tk_VRI_image.sh` is the top-level shell wrapper for making the final VRI visual products in this directory.

It starts from existing VRI FITS products:

```text
*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits
```

and produces, for each galaxy:

```text
<GALAXY>_observed_VRI.png
<GALAXY>_observed_VRI.pdf
<GALAXY>_legacy_reprojected.jpg
<GALAXY>_combined_VRI.png
```

The shell script does not create the VRI FITS files. Those are upstream prerequisites produced by `v3tk_to_VRI.py`, normally via `v3tk_to_VRI.sh`.

## Main Command

Run the full VRI image pipeline from:

```bash
cd /Users/Igniz/Desktop/ICRAR/v3tk_to_VRI
./v3tk_VRI_image.sh
```

The script writes a timestamped log:

```text
v3tk_VRI_image_YYYYMMDD_HHMMSS.log
```

Example:

```text
v3tk_VRI_image_20260521_131131.log
```

## Pipeline Summary

`v3tk_VRI_image.sh` runs three Python scripts:

```text
Step 1: v3tk_get_legacy.py
Step 2: v3tk_observed_VRI_image.py
Step 3: v3tk_combined_VRI_image.py
```

Data flow:

```text
*_v3tk_VRI.fits
  -> *_legacy_reprojected.jpg
  -> *_observed_VRI.png / *_observed_VRI.pdf
  -> *_combined_VRI.png
```

Then `auto_arrange_and_combine.py` can be run separately to arrange the per-galaxy `*_combined_VRI.png` panels into a single mosaic.

## Top-Level Shell Script

File:

```text
v3tk_VRI_image.sh
```

Important behavior:

- `set -euo pipefail` makes the script stop on unhandled errors, unset variables, or failed pipe components.
- It changes directory to the directory containing the script, so relative paths are stable.
- It creates a timestamped log file and tees all output to that log.
- It detects CPU count using Python, `getconf`, `nproc`, or `sysctl`.
- It tries to activate the `ICRAR` conda environment.
- If `ICRAR` is not active, it falls back to `conda run -n ICRAR python`.
- It runs all three stages with `--overwrite`.

The script-level controls are:

```bash
MAX_DL=1 ./v3tk_VRI_image.sh
```

`MAX_DL` limits concurrent Legacy Survey downloads. The default is:

```text
MAX_DL=2
```

This throttle matters because the Legacy Survey service can rate-limit simultaneous HTTP requests.

The stage order is intentional: the standalone observed renderer now uses the
reprojected Legacy image to match its inner footprint boundary. Therefore the
Legacy product must exist before the observed PNG/PDF is rendered.

## Upstream Prerequisite: `v3tk_to_VRI.sh` and `v3tk_to_VRI.py`

Although `v3tk_VRI_image.sh` does not call the converter, its input files come from:

```text
v3tk_to_VRI.sh
v3tk_to_VRI.py
```

`v3tk_to_VRI.sh`:

- logs to `v3tk_to_VRI.log`;
- reads compressed MUSE cube files from:

```text
/arc/projects/mauve/cubes/v3tk
```

- copies each `*_v3tk.fits.gz` into the working directory;
- unzips it;
- runs:

```bash
python v3tk_to_VRI.py <cube.fits>
```

- removes the local uncompressed input cube after conversion.

`v3tk_to_VRI.py` collapses a 3D MUSE cube into V/R/I flux and AB magnitude maps. It writes a VRI FITS file with these HDUs:

```text
V_FLUX
V_MAG
R_FLUX
R_MAG
I_FLUX
I_MAG
```

Important converter defaults:

- input cube HDU: `DATA`;
- V filter: `bessell-V`;
- R filter: `bessell-R`;
- I filter: `bessell-I`;
- raw cube flux scale: `1e-20 erg/s/cm^2/Angstrom`;
- default NaN policy: replace non-finite cube samples with zero for visualization;
- partial filter overlap is enabled by default, so filters can be truncated and renormalized to the cube wavelength range.

The converter preserves the spatial WCS in the 2D output HDUs.

## Observed VRI Rendering (Pipeline Step 2)

File:

```text
v3tk_observed_VRI_image.py
```

Command used by the shell script:

```bash
python v3tk_observed_VRI_image.py --input-dir . --overwrite --workers "$NCPU"
```

Required inputs for each galaxy are the VRI FITS product and the matching:

```text
<GALAXY>_legacy_reprojected.jpg
```

Default input pattern:

```text
*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits
```

For each matched FITS file, it writes:

```text
<GALAXY>_observed_VRI.png
<GALAXY>_observed_VRI.pdf
```

The galaxy ID is parsed from the FITS filename by removing:

```text
_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits
```

Rendering details:

- reads `V_FLUX`, `R_FLUX`, and `I_FLUX`;
- prepares non-finite values as zero;
- clips negative values to zero for display;
- maps channels as:

```text
I -> red
R -> green
V -> blue
```

- applies a gray-world white balance;
- scales channels using luminance percentiles;
- uses Astropy Lupton/asinh RGB rendering;
- applies a color-preserving, highlight-protected gamma lift through one shared
  multiplier per RGB pixel;
- builds a circular galaxy-centre reference aperture around the valid-footprint
  bounding-box centre, with radius `0.22` times the shorter footprint dimension;
- uses the maximum unquantized Rec.709 luminance in that aperture to set one
  global VRI multiplier, targeting `255/255` by default, so foreground stars
  outside the aperture do not determine the galaxy scale;
- applies any required RGB-gamut protection through an additional shared
  multiplier per pixel, so one highly colored pixel cannot limit the dynamic
  range of the whole galaxy or clip channels independently; a colored central
  reference may therefore reach channel 255 with luminance slightly below 255
  rather than being forced to white;
- calculates a smooth Legacy color/brightness baseline only from adjacent
  pixels outside the VRI footprint and uses it to adjust the inner 20-pixel VRI
  edge region by default;
- ignores all Legacy pixels covered by the VRI footprint and retains the
  complete high-frequency VRI residual so the edge stays sharp;
- never copies or alpha-blends raw Legacy RGB into a valid VRI pixel;
- leaves pixels outside the standalone observed footprint black;
- writes PNG at native image dimensions;
- writes a PDF with the same visual appearance.

Important defaults:

```text
percentile-low: 1.0
percentile-high: 99.9
stretch: 1.0
Q: 8.0
gamma: 0.65
post-boost: 2.75
target-max-luminance: 255.0 (on the 0-255 display scale)
center-radius-fraction: 0.22 (of the shorter valid-footprint dimension)
legacy-transition-width: 20 pixels
workers: os.cpu_count()
overwrite: true
```

Useful dry run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python v3tk_observed_VRI_image.py --dry-run
```

Potential warning:

If Matplotlib cannot write to `~/.matplotlib`, it creates a temporary cache directory. This is not a pipeline failure, but setting `MPLCONFIGDIR` to a writable directory can reduce startup overhead in multiprocessing runs.

Example:

```bash
MPLCONFIGDIR=/private/tmp/mplconfig ./v3tk_VRI_image.sh
```

## Legacy Survey Download and Reprojection (Pipeline Step 1)

File:

```text
v3tk_get_legacy.py
```

Command used by the shell script:

```bash
python v3tk_get_legacy.py --overwrite --workers "$NCPU" --max-concurrent-downloads "$MAX_DL"
```

For each matched VRI FITS file, it writes:

```text
<GALAXY>_legacy_reprojected.jpg
```

Temporary cutout files:

```text
<GALAXY>_Legacy.jpg
<GALAXY>_Legacy.fits
```

By default, the temporary JPEG is deleted after successful reprojection, and the temporary FITS is deleted immediately after WCS use. Use `--keep-cutout` if the downloaded JPEG should be retained.

What the script does:

- selects a 2D flux HDU from the VRI FITS, preferring `R_FLUX`;
- gets the target MUSE WCS and shape;
- computes the image-center sky coordinate;
- estimates the MUSE footprint size from the WCS corners;
- requests a Legacy Survey cutout centered on the MUSE field;
- downloads both a JPEG cutout and a FITS cutout;
- uses the FITS cutout header for authoritative Legacy WCS;
- verifies that the cutout covers the target MUSE footprint;
- enlarges the requested cutout up to three attempts if coverage is insufficient;
- reprojects each RGB channel to the MUSE grid using `reproject_adaptive(conserve_flux=True)`;
- writes the reprojected RGB image as JPEG.

Important defaults:

```text
layer: ls-dr9
pixscale: 0.262 arcsec/pixel
margin-factor: 1.02
timeout: 120 seconds
workers: min(2, os.cpu_count()) when run directly
max-concurrent-downloads: 2
overwrite: true
```

The shell script passes `--workers "$NCPU"`, but download concurrency is separately limited by `MAX_DL`.

Useful dry run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python v3tk_get_legacy.py --dry-run
```

If HTTP 429 or 503 occurs, the downloader retries with exponential backoff and random jitter. If rate limits continue, run:

```bash
MAX_DL=1 ./v3tk_VRI_image.sh
```

## Combine Observed MUSE and Legacy Background (Pipeline Step 3)

File:

```text
v3tk_combined_VRI_image.py
```

Command used by the shell script:

```bash
python v3tk_combined_VRI_image.py --overwrite --workers "$NCPU"
```

For each matched VRI FITS file, it requires:

```text
<GALAXY>_observed_VRI.png
<GALAXY>_legacy_reprojected.jpg
```

and writes:

```text
<GALAXY>_combined_VRI.png
```

Combination logic:

- reads a 2D flux HDU from the VRI FITS, preferring `R_FLUX`;
- builds a valid-pixel mask:

```text
finite and > 0
```

- flips the mask vertically to match the PNG orientation written by `v3tk_observed_VRI_image.py`;
- removes the Legacy reprojected background wherever the valid mask is true;
- inserts observed VRI RGB pixels opaquely at those valid mask pixels, with no
  alpha blend or overlapping Legacy contribution;
- retains Legacy RGB only where the valid mask is false;
- writes the result as PNG.

The combined image therefore shows:

- MUSE observed VRI rendering where the VRI flux map is valid;
- Legacy Survey background where the MUSE rendering is invalid, blank, NaN, or non-positive.

Size checks:

- observed PNG must match the VRI FITS 2D shape;
- legacy reprojected JPG must match the VRI FITS 2D shape;
- mismatch raises an error rather than silently combining misaligned images.

Useful dry run:

```bash
/opt/miniconda3/envs/ICRAR/bin/python v3tk_combined_VRI_image.py --dry-run
```

## Current Workspace Snapshot

On 2026-05-21, this directory contains 26 VRI FITS inputs matching:

```text
*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits
```

The same directory also contains existing per-galaxy products such as:

```text
IC3392_observed_VRI.png
IC3392_observed_VRI.pdf
IC3392_legacy_reprojected.jpg
IC3392_combined_VRI.png
```

and equivalent products for other galaxies.

## Relevant Files

Top-level runner:

```text
v3tk_VRI_image.sh
```

Upstream VRI FITS generation:

```text
v3tk_to_VRI.sh
v3tk_to_VRI.py
```

Pipeline stages called by `v3tk_VRI_image.sh`:

```text
v3tk_observed_VRI_image.py
v3tk_get_legacy.py
v3tk_combined_VRI_image.py
```

Optional mosaic arranger after the per-galaxy combined images are made:

```text
auto_arrange_and_combine.py
```

Related R-band analogues in the repository:

```text
v3tk_R_image.sh
v3tk_observed_R_image.py
v3tk_combined_R_image.py
```

These are not run by `v3tk_VRI_image.sh`, but they follow a similar pattern for R-only products.

## Dependency State

The validated Python environment is:

```text
/opt/miniconda3/envs/ICRAR/bin/python
```

Relevant package versions checked on 2026-05-21:

```text
numpy 2.2.6
astropy 7.0.2
speclite 0.20
matplotlib 3.10.3
Pillow 11.3.0
reproject 0.14.1
```

Stage-specific dependencies:

- `v3tk_to_VRI.py`: `numpy`, `astropy`, `speclite`;
- `v3tk_observed_VRI_image.py`: `numpy`, `astropy`, `matplotlib`, `Pillow`;
- `v3tk_get_legacy.py`: `numpy`, `astropy`, `reproject`, `Pillow`, Python `urllib`;
- `v3tk_combined_VRI_image.py`: `numpy`, `astropy`, `Pillow`.

## Failure Modes and Checks

No VRI FITS matched:

```text
No files matched pattern '*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits'
```

This means the upstream conversion has not been run, the script is in the wrong directory, or the pattern differs.

Missing Legacy image during observed step:

```text
Missing Legacy image required by observed renderer: <GALAXY>_legacy_reprojected.jpg.
Run v3tk_get_legacy.py first.
```

This happens if pipeline Step 2 is run independently before Step 1.

Size mismatch:

The scripts intentionally fail if FITS dimensions, observed PNG size, and reprojected Legacy image size do not agree. This protects against combining products from different runs or stale files.

Download rate limiting:

The Legacy Survey server may return rate-limit or temporary service errors. The downloader retries, but `MAX_DL=1` is safer if repeated HTTP errors occur.

Network required:

Step 1 requires internet access to:

```text
https://www.legacysurvey.org/viewer/jpeg-cutout
https://www.legacysurvey.org/viewer/fits-cutout
```

## Recommended Routine

If VRI FITS products already exist:

```bash
cd /Users/Igniz/Desktop/ICRAR/v3tk_to_VRI
MAX_DL=2 ./v3tk_VRI_image.sh
```

If Legacy Survey rate-limits the run:

```bash
MAX_DL=1 ./v3tk_VRI_image.sh
```

After per-galaxy combined images exist, make a combined mosaic:

```bash
/opt/miniconda3/envs/ICRAR/bin/python auto_arrange_and_combine.py *combined_VRI.png 16 9
```

For a quick mosaic without proof:

```bash
/opt/miniconda3/envs/ICRAR/bin/python auto_arrange_and_combine.py --fast *combined_VRI.png 16 9
```

## Important Scientific and Display Assumptions

The VRI FITS files are visualization products, not calibrated standard broadband photometry in every case. Partial filter overlap is enabled by default in `v3tk_to_VRI.py`, so if the MUSE cube does not fully cover a Bessell filter, the filter response is truncated and renormalized to the cube wavelength range.

The observed PNG/PDF products are display renderings:

- non-finite samples are handled for visualization;
- negative display values are clipped;
- gray-world white balance and Lupton/asinh stretch affect visual color;
- `gamma=0.65` lifts faint and mid-tone detail after the Lupton mapping;
- `post_boost=2.75` controls the strength of that lift;
- for maximum channel value `m`, the curve is
  `m_out = m + post_boost * (m**gamma - m) * (1 - m)`, clipped to the display
  range only as a numerical safeguard;
- the `(1 - m)` shoulder fades the lift toward white, protecting bright-region
  contrast, while one shared multiplier per pixel preserves mapped RGB ratios.
- one global VRI multiplier is set by the maximum unquantized Rec.709 luminance
  inside a centered circular aperture whose radius is `0.22` times the shorter
  valid-footprint dimension, targeting `255/255` by default;
- foreground stars outside that aperture may also saturate, but do not set the
  galaxy scale; shared-channel gamut protection can leave a colored central
  reference slightly below Rec.709 luminance 255 rather than destroy its color;
- only adjacent Legacy pixels outside the VRI footprint contribute to the
  low-frequency RGB baseline used through the inner 20-pixel footprint region;
- high-frequency VRI detail is retained, and covered or raw Legacy pixels are
  never copied or alpha-blended into the VRI footprint;
- the edge treatment is intentionally display-oriented and does not alter the
  VRI FITS measurements.

The combined image is a visual overlay product:

- MUSE VRI rendering is used where `R_FLUX` is finite and positive;
- Legacy Survey is used elsewhere;
- the observed VRI edge uses a 20-pixel low-frequency baseline derived only
  from adjacent exterior Legacy pixels, without a raw seam or image overlap;
- the mask is based on the selected flux map, not on all three V/R/I channels simultaneously.

These choices are appropriate for visually inspecting the MUSE footprint against a Legacy Survey background, but they should not be interpreted as a direct calibrated RGB photometric product without checking the underlying FITS maps.
