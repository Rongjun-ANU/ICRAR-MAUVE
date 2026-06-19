# v3tk_to_VRI

This repository contains the scripts and generated products for converting MAUVE/MUSE `v3tk` cubes and the supported PHANGS-MUSE native public cubes into VRI visual products, reprojecting Legacy Survey RGB cutouts onto the MUSE footprint, combining the observed and Legacy images, and arranging the per-galaxy outputs into fixed-aspect-ratio mosaics.

The repository intentionally includes the current data products as well as the code:

- VRI FITS products: `*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits` and supported `*_PHANGS_DATACUBE_native_VRI.fits`; downstream image scripts also accept `.fits.gz` products directly.
- observed VRI renderings: `*_observed_VRI.png` and `*_observed_VRI.pdf`
- reprojected Legacy Survey backgrounds: `*_legacy_reprojected.jpg`
- combined observed-on-Legacy products: `*_combined_VRI.png`
- all-galaxy mosaics: `All_combined_VRI*.png` and `ALL_combined_VRI.png`
- arrangement reports: `*.proof.txt`
- run logs and dated documentation notes

At the time this README was written, the folder contains 26 VRI FITS products and no individual file above GitHub's normal 100 MB file limit. The local checkout is about 825 MB including `.git`, so pushes and clones may take some time.

## Repository Layout

| File pattern | Purpose |
| --- | --- |
| `v3tk_to_VRI.py` | Converts one 3D MUSE cube into V, R, and I flux/magnitude maps. It accepts both `.fits` and `.fits.gz` inputs directly, or a GALID when the matching v3tk cube is in the current directory. |
| `v3tk_to_VRI.sh` | Batch wrapper that first checks the current directory for selected local `*_v3tk.fits` or `*_v3tk.fits.gz` cubes, then falls back to `/arc/projects/mauve/cubes/v3tk`, runs `v3tk_to_VRI.py` directly on local filesystem inputs without copying, unzipping, or deleting them, and stages only supported PHANGS native public `vos:` inputs with `vcp` when needed. |
| `v3tk_observed_VRI_image.py` | Renders the VRI FITS products into native-size PNG and PDF images. It discovers both `.fits` and `.fits.gz` products when the pattern ends in `.fits`. |
| `v3tk_get_legacy.py` | Downloads Legacy Survey cutouts, obtains matching WCS information, and reprojects the RGB image to the MUSE grid. It uses the same `.fits`/`.fits.gz` adaptive product discovery. |
| `v3tk_combined_VRI_image.py` | Combines the observed VRI rendering with the reprojected Legacy image, using valid MUSE pixels where available and Legacy pixels outside the MUSE footprint. It also accepts compressed VRI FITS products. |
| `v3tk_VRI_image.sh` | Top-level image pipeline wrapper for observed rendering, Legacy reprojection, and combined image generation. |
| `auto_arrange_and_combine.py` | Packs many per-galaxy image panels into one fixed-ratio mosaic. Uses OR-Tools when available for proof-aware layouts and can fall back to fast heuristic layouts. |
| `v3tk_observed_R_image.py`, `v3tk_combined_R_image.py`, `v3tk_R_image.sh` | Older or parallel R-band-only workflow files. |
| `20260521_*_documentation.md` | Detailed notes for the VRI image pipeline and the mosaic arranger. |

## Current Data Products

The current checked-in galaxy set is:

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

For each galaxy, the usual complete VRI product set is:

```text
<GALAXY>_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits
<GALAXY>_observed_VRI.png
<GALAXY>_observed_VRI.pdf
<GALAXY>_legacy_reprojected.jpg
<GALAXY>_combined_VRI.png
```

The main all-galaxy mosaic products currently include:

```text
ALL_combined_VRI.png
All_combined_VRI_16_9.png
All_combined_VRI_99_62.png
All_combined_VRI_152_102.png
```

Each `All_combined_VRI*.png` product has a matching `*.proof.txt` report describing the arrangement status, canvas size, input image sizes, placements, and whether a mathematical proof of optimality was obtained.

## Pipeline Summary

The full VRI workflow is:

```text
raw upstream cube: *_v3tk.fits or *_v3tk.fits.gz
  -> v3tk_to_VRI.py
  -> *_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits
  -> v3tk_observed_VRI_image.py
  -> *_observed_VRI.png / *_observed_VRI.pdf
  -> v3tk_get_legacy.py
  -> *_legacy_reprojected.jpg
  -> v3tk_combined_VRI_image.py
  -> *_combined_VRI.png
  -> auto_arrange_and_combine.py
  -> All_combined_VRI*.png + *.proof.txt
```

The supported PHANGS-native path is the same after conversion, with the input and first product names changed:

```text
public PHANGS native cube: NGC4254_PHANGS_DATACUBE_native.fits
  -> v3tk_to_VRI.py
  -> NGC4254_PHANGS_DATACUBE_native_VRI.fits
  -> image, Legacy, combined-image, and mosaic stages
```

The three supported public PHANGS-MUSE native cubes are:

```text
vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4254_PHANGS_DATACUBE_native.fits
vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4321_PHANGS_DATACUBE_native.fits
vos:phangs/RELEASES/PHANGS-MUSE/DR1.0/DATACUBES/NGC4535_PHANGS_DATACUBE_native.fits
```

Use the `native` PHANGS products for consistency with the MAUVE native-grid workflow. Full VRI conversion reads the whole cube. Local filesystem cubes, including `.fits.gz`, are passed directly to Astropy and are not copied, unzipped, or deleted by the wrapper. Only public `vos:` inputs are staged into the working directory with `vcp`. Header-only checks or small cutouts can use CANFAR VOSspace helpers such as `vcat --head` or `vcp` cutout syntax, but the full-cube workflow should not treat a `vos:` URI as a normal local FITS path.

The converter writes these FITS HDUs:

```text
V_FLUX
V_MAG
R_FLUX
R_MAG
I_FLUX
I_MAG
```

The image rendering maps channels as:

```text
I -> red
R -> green
V -> blue
```

The combined image shows observed MUSE VRI pixels where the VRI flux map is valid and reprojected Legacy Survey RGB pixels outside the observed footprint.

## Environment

The scripts are designed to run in the local `ICRAR` conda environment. The Python dependencies used across the workflow are:

```text
numpy
astropy
speclite
matplotlib
pillow
reproject
ortools
```

Install equivalent dependencies with conda-forge if setting up a new machine:

```bash
conda create -n ICRAR -c conda-forge python numpy astropy speclite matplotlib pillow reproject ortools
conda activate ICRAR
```

`ortools` is only required for certified mosaic optimization in `auto_arrange_and_combine.py`. If OR-Tools is unavailable, the arranger can still run in heuristic-only mode with `--fast`.

## Main Commands

Run commands from the repository root:

```bash
cd /Users/Igniz/Desktop/ICRAR/v3tk_to_VRI
```

### 1. Convert upstream cubes to VRI FITS

Use this only on a system where the upstream cube directory exists:

```bash
./v3tk_to_VRI.sh
```

To restrict conversion to selected galaxies, pass the GALIDs explicitly. This
is the recommended check before processing the three PHANGS-native MAUVE
galaxies:

```bash
./v3tk_to_VRI.sh --dry-run NGC4254 NGC4321 NGC4535
./v3tk_to_VRI.sh NGC4254 NGC4321 NGC4535
```

With positional GALID arguments, the wrapper selects only those galaxies. It
first checks the current working directory for matching local v3tk cube names
and uses those inputs in place without deleting them, then falls back to
`/arc/projects/mauve/cubes/v3tk`. For
`NGC4254`, `NGC4321`, and `NGC4535`, it selects
`*_PHANGS_DATACUBE_native.fits` from the local target directory when present, or
from public PHANGS VOSspace via `vcp` when not present. It does not fall through
to the all-galaxy local discovery list.

The wrapper expects:

```text
/arc/projects/mauve/cubes/v3tk/*_v3tk.fits
/arc/projects/mauve/cubes/v3tk/*_v3tk.fits.gz
```

It also stages the three supported PHANGS-native public cubes listed above with `vcp` if they are not already present in the local target directory. If both compressed and uncompressed copies exist for the same local v3tk cube, the wrapper processes the first matching stem once. For `NGC4254`, `NGC4321`, and `NGC4535`, the PHANGS-native source is preferred over an old local v3tk-style source. It writes `v3tk_to_VRI.log`; during each per-cube run it passes local `.fits` or `.fits.gz` inputs directly to `v3tk_to_VRI.py` and writes the `_VRI.fits` output in the working directory. It does not unzip, copy, or delete local filesystem input cubes.

To convert a single local cube directly:

```bash
python v3tk_to_VRI.py path/to/cube.fits
python v3tk_to_VRI.py path/to/cube.fits.gz
python v3tk_to_VRI.py NGC4380
```

The GALID form expects the matching `NGC4380_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits`
or `.fits.gz` cube to be in the current directory.

Useful converter options:

```bash
python v3tk_to_VRI.py path/to/cube.fits --row-chunk 8
python v3tk_to_VRI.py path/to/cube.fits --output MY_OUTPUT_VRI.fits
python v3tk_to_VRI.py path/to/cube.fits --no-allow-partial-overlap
```

When no `--output` is supplied, `v3tk_to_VRI.py` strips either `.fits` or `.fits.gz`
before appending `_VRI.fits`, so `cube.fits.gz` writes `cube_VRI.fits`.

### 2. Generate observed, Legacy, and combined VRI images

The main image pipeline is:

```bash
./v3tk_VRI_image.sh
```

The wrapper runs:

```bash
python v3tk_observed_VRI_image.py --input-dir . --overwrite --workers "$NCPU"
python v3tk_get_legacy.py --overwrite --workers "$NCPU" --max-concurrent-downloads "$MAX_DL"
python v3tk_combined_VRI_image.py --overwrite --workers "$NCPU"
```

`MAX_DL` controls the number of simultaneous Legacy Survey downloads. The default is 2. Use a lower value if the Legacy Survey service rate-limits requests:

```bash
MAX_DL=1 ./v3tk_VRI_image.sh
```

Dry-run the individual stages before writing files:

```bash
python v3tk_observed_VRI_image.py --dry-run
python v3tk_get_legacy.py --dry-run
python v3tk_combined_VRI_image.py --dry-run
```

The image-stage scripts default to `*_DATACUBE*_VRI.fits`, and `.fits` patterns also discover matching `.fits.gz` counterparts. This lets both the older v3tk VRI FITS products and the PHANGS-native VRI products be rendered and combined without manual unzip.

### 3. Build an all-galaxy mosaic

Certified/default mode:

```bash
python auto_arrange_and_combine.py *combined_VRI.png 16 9
```

Quick heuristic mode without a mathematical proof:

```bash
python auto_arrange_and_combine.py --fast *combined_VRI.png 16 9
```

Longer certified attempt:

```bash
python auto_arrange_and_combine.py --time-limit 600 *combined_VRI.png 16 9
```

Custom proof/report path:

```bash
python auto_arrange_and_combine.py --proof-file VRI_16_9_arrangement_report.txt *combined_VRI.png 16 9
```

The arranger skips existing all-galaxy outputs beginning with `ALL_`, `All_`, or `All`, so previous mosaics are not accidentally used as input panels.

## Mosaic Proof Reports

`auto_arrange_and_combine.py` writes a report next to the output image, usually:

```text
All_combined_VRI_16_9.proof.txt
```

Important fields:

| Field | Meaning |
| --- | --- |
| `ratio` | Requested fixed output aspect ratio, such as `16:9`. |
| `status` | `OPTIMAL`, `FEASIBLE`, or `HEURISTIC_ONLY`. |
| `final_k` | Integer scale factor for the output canvas. |
| `final_canvas` | Final mosaic size in pixels. |
| `density` | Fraction of the canvas area occupied by input panels. |
| `placement table` | Input file names, dimensions, and `(x, y)` paste positions. |

Interpretation:

- `OPTIMAL` means OR-Tools certified that no smaller integer canvas exists for the requested ratio under the script's rectangle-packing assumptions.
- `FEASIBLE` means the saved layout is valid and non-overlapping, but global optimality was not proven within the time limit.
- `HEURISTIC_ONLY` means `--fast` was used and no mathematical proof was attempted.

The arranger can also reuse compatible existing proof reports as warm starts and synchronize compatible sibling mosaics so that related outputs share the same layout.

## GitHub Notes

This repository is configured with:

```text
origin  https://github.com/Rongjun-ANU/v3tk_to_VRI.git
branch  main
```

Because this repository intentionally includes generated FITS/images/logs, it is larger than a code-only repository. This is acceptable for the current contents because no individual file is over 100 MB. If a future raw cube, mosaic, archive, or derived product exceeds 100 MB, GitHub will reject a normal push for that file; use Git LFS or keep that file outside the repository.

If FITS products are compressed before sharing, track both FITS extensions:

```bash
git lfs track '*.fits'
git lfs track '*.fits.gz'
```

Typical Git commands:

```bash
git status
git add README.md
git commit -m "Add project README"
git push origin main
```

## Existing Detailed Notes

The dated documentation files contain deeper operational notes:

- `20260521_v3tk_VRI_image_pipeline_documentation.md`
- `20260521_auto_arrange_and_combine_documentation.md`

Use those when debugging the pipeline or checking details such as default render parameters, Legacy reprojection behavior, OR-Tools runtime behavior, and proof-report interpretation.
