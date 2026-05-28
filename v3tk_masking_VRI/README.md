# v3tk Masking VRI

This repository contains the VRI-band v3tk masking products and helper scripts used to
build nGIST-compatible spatial masks for a Virgo/MUSE galaxy sample. It includes the
input VRI datacubes, the generated binary FITS masks, diagnostic overlay images, log
output, and combined overview mosaics.

The current snapshot is intended to be a self-contained record of the masking run:

- 26 input `*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits` files; gzip-compressed `.fits.gz` inputs are also accepted by the masking script.
- 26 output `*_mask.fits` files.
- 26 original `*_combined_VRI.png` reference images.
- 26 diagnostic `*_combined_VRI_mask.png` overlays.
- Combined all-galaxy overview mosaics and proof reports.
- Two Python scripts used for mask generation and image arrangement.

## Repository Contents

| Path or pattern | Count | Description |
| --- | ---: | --- |
| `make_ngist_masks_from_catalogs_VRI.py` | 1 | Main masking script. Generates binary spatial masks and diagnostic overlays from VRI FITS datacubes. |
| `auto_arrange_and_combine.py` | 1 | Mosaic arranger. Packs per-galaxy PNGs into dense fixed-ratio overview images, with optional OR-Tools proof reports. |
| `*_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits` or `.fits.gz` | 26 | Input VRI FITS datacubes used for WCS, pixel geometry, and science-data footprint. Compressed inputs are read directly with Astropy and do not need to be unzipped first. |
| `*_mask.fits` | 26 | Binary nGIST-compatible spatial masks. `0` means unmasked; `1` means masked. |
| `*_combined_VRI.png` | 26 | Per-galaxy VRI reference images. |
| `*_combined_VRI_mask.png` | 26 | Per-galaxy diagnostic overlays showing masked foreground/background objects. |
| `All_combined_VRI_16_9.png` | 1 | 16:9 combined mosaic of original VRI images. |
| `All_combined_VRI_mask_16_9.png` | 1 | 16:9 combined mosaic of masked diagnostic overlays. |
| `ALL_combined_VRI_mask.png` | 1 | Earlier combined masked mosaic. |
| `*.proof.txt` | 3 | Layout proof or fallback reports from the arranger. |
| `v3tk_masking.log` | 1 | Full log from the masking run. |

The folder size is roughly 537 MB. At the time this README was written, no file was
larger than GitHub's 100 MB hard per-file limit. The largest files are the combined
mosaic PNGs at about 58 MB each, which GitHub may warn about but should accept.

## Upstream Input Image Preparation

This repository assumes that the per-galaxy VRI reference images already exist as
`*_combined_VRI.png` files. If starting from earlier MAUVE/MUSE products, first check
out the companion company repository:

[Rongjun-ANU/v3tk_to_VRI](https://github.com/Rongjun-ANU/v3tk_to_VRI)

That repository is used to create the input images for this masking repository: RGB
images from the VRI bands of available MAUVE-MUSE galaxies, combined with rescaled
Legacy Survey images as the background. After those `*_combined_VRI.png` images and
matching VRI FITS products are available, use this repository to generate the spatial
masks, diagnostic mask overlays, and all-galaxy overview mosaics.

## Galaxy Sample

The included galaxies are:

`IC3392`, `NGC4064`, `NGC4189`, `NGC4192`, `NGC4293`, `NGC4294`, `NGC4298`,
`NGC4302`, `NGC4330`, `NGC4351`, `NGC4383`, `NGC4388`, `NGC4394`, `NGC4396`,
`NGC4402`, `NGC4405`, `NGC4419`, `NGC4457`, `NGC4501`, `NGC4522`,
`NGC4567_8`, `NGC4580`, `NGC4606`, `NGC4607`, `NGC4694`, and `NGC4698`.

## Masking Workflow

The masking script is:

```bash
python make_ngist_masks_from_catalogs_VRI.py
```

With no command-line arguments, it processes all files matching:

```text
*_DATACUBE*_VRI.fits
*_DATACUBE*_VRI.fits.gz
```

To process only selected galaxies, pass one or more shell patterns:

```bash
python make_ngist_masks_from_catalogs_VRI.py 'NGC4383*_DATACUBE*_VRI.fits'
python make_ngist_masks_from_catalogs_VRI.py 'NGC4383*_DATACUBE*_VRI.fits.gz'
python make_ngist_masks_from_catalogs_VRI.py 'NGC4383*_DATACUBE*_VRI.fits' 'NGC4419*_DATACUBE*_VRI.fits'
```

When a command-line pattern ends in `.fits`, the script also checks the corresponding `.fits.gz` pattern. This means the common `.fits` commands above work when only compressed inputs are present.

For each input `XXX_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk_VRI.fits` or `.fits.gz`, the script writes:

- `XXX_mask.fits`: binary mask in the FITS spatial grid.
- `XXX_combined_VRI_mask.png`: diagnostic overlay image.
- `v3tk_masking.log`: run log, overwritten on each full script run.

The mask convention is:

- `0`: unmasked spaxel, usually target galaxy or sky.
- `1`: masked spaxel, usually foreground star or background object.

## Masking Logic

`make_ngist_masks_from_catalogs_VRI.py` is designed as a conservative pre-processing
tool for nGIST. In broad terms it does the following:

1. Loads the VRI FITS datacube or image and derives a 2D WCS.
2. Uses the `R_FLUX` extension, when available, to define the valid MUSE field of view.
3. Queries Gaia DR3 for foreground-star candidates.
4. Builds an optional target-galaxy footprint from `R_MAG` to avoid masking target-galaxy structure as background.
5. Queries Legacy Surveys DR9 through NOIRLab Data Lab TAP when available.
6. Uses spectroscopic or catalog evidence from SDSS and NED for background-object checks when needed.
7. Uses PS1, VizieR, SkyMapper, or SDSS fallbacks when the higher-quality catalog path is unavailable.
8. Rasterizes the selected foreground/background objects into a binary FITS mask.
9. Writes a diagnostic PNG overlay with outlines on the galaxy image or FITS-derived background.

Important implementation details:

- The script prefers pixel-locked overlays so the mask can be visually checked against the FITS grid.
- Gaia masking defaults to a foreground-selection mode based on parallax and proper motion evidence.
- Legacy DR9 background-galaxy masking uses morphology and photo-z information where available.
- The target-footprint veto is intended to reduce accidental masking of target-galaxy HII regions or internal structure.
- The script continues with warnings if optional catalog services are unavailable, but mask completeness may change.

## Mosaic Workflow

The mosaic script is:

```bash
python auto_arrange_and_combine.py '*_combined_VRI.png' 16 9
python auto_arrange_and_combine.py '*_combined_VRI_mask.png' 16 9
```

It loads the individual PNGs at native pixel size and packs them into a fixed-ratio
canvas. By default it uses a deterministic heuristic warm start and OR-Tools CP-SAT
optimization. It writes both an output PNG and a `.proof.txt` report.

Useful options:

```bash
python auto_arrange_and_combine.py '*_combined_VRI_mask.png' 16 9 --time-limit 300
python auto_arrange_and_combine.py '*_combined_VRI_mask.png' 16 9 --fast
python auto_arrange_and_combine.py '*_combined_VRI_mask.png' 16 9 --no-reuse-existing
```

The script skips files beginning with `ALL_` or `All_` so previously generated mosaics
are not accidentally packed into new mosaics.

## Python Environment

Recommended Python version: Python 3.10 or newer.

Core dependencies for mask generation:

```bash
python -m pip install numpy astropy matplotlib
```

Recommended optional dependencies:

```bash
python -m pip install pillow astroquery pyvo scipy
```

Dependency for OR-Tools-backed mosaic proof reports:

```bash
python -m pip install ortools
```

One combined install command is:

```bash
python -m pip install numpy astropy matplotlib pillow astroquery pyvo scipy ortools
```

## Reproducing the Current Products

From the repository root:

If the `*_combined_VRI.png` inputs are missing, first generate them using the companion
`Rongjun-ANU/v3tk_to_VRI` workflow described above. Then run:

```bash
python make_ngist_masks_from_catalogs_VRI.py
python auto_arrange_and_combine.py '*_combined_VRI.png' 16 9
python auto_arrange_and_combine.py '*_combined_VRI_mask.png' 16 9
```

After regeneration, inspect:

- `v3tk_masking.log` for catalog-query warnings or failed galaxies.
- Each `*_combined_VRI_mask.png` for visually incorrect masks.
- Each `*.proof.txt` for the mosaic layout status and geometry.

## Git Notes

This repository intentionally tracks the FITS and PNG data products for reproducibility.
The virtual environment and local cache files are ignored through `.gitignore`.

Ignored local-only files:

- `.venv/`
- `.DS_Store`
- `__pycache__/`
- `.pytest_cache/`
- `.ipynb_checkpoints/`

If future products exceed GitHub's 100 MB per-file limit, use Git LFS for those files
before committing:

```bash
git lfs install
git lfs track '*.fits'
git lfs track '*.fits.gz'
git lfs track '*.png'
git add .gitattributes
```

For the current snapshot, plain Git is sufficient because all files are below 100 MB.
