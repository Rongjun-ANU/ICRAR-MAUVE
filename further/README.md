# MAUVE further post-processing

This folder contains the active MAUVE second-project post-processing scripts
that run on nGIST products under `v3tk_v7.6.8/{GALID}`.

## PHANGS-native cube inputs

`Mass.py` and `proxy_EWHa.py` read full datacubes from a local filesystem path.
By default the shell wrappers point `--cube-root` at:

```text
/arc/projects/mauve/cubes/v3tk
```

The three MAUVE galaxies that overlap the PHANGS-MUSE DR1 native public release
use these cube filenames:

```text
NGC4254_PHANGS_DATACUBE_native.fits
NGC4321_PHANGS_DATACUBE_native.fits
NGC4535_PHANGS_DATACUBE_native.fits
```

These are expected to be staged locally, for example on Setonix scratch, before
full Mass or proxy-EW processing. Do not point these scripts at a `vos:` URI for
full-cube processing. Use CANFAR tools such as `vcat --head` or cutouts only for
header-only or small inspection jobs; full FITS integration should use a local
staged file.

## Running selected galaxies

Both wrappers accept positional GALID arguments. With no arguments they run the
default sample, including the three PHANGS-native galaxies. With arguments they
run only the requested subset:

```bash
./Mass.sh NGC4254
./Mass.sh NGC4254 NGC4321 NGC4535

./proxy_EWHa.sh NGC4254
./proxy_EWHa.sh NGC4254 NGC4321 NGC4535
```

The Python scripts can also be run directly with explicit roots:

```bash
python Mass.py -g NGC4254 --root . --fallback-root . --cube-root /arc/projects/mauve/cubes/v3tk
python proxy_EWHa.py -g NGC4254 --root . --fallback-root . --cube-root /arc/projects/mauve/cubes/v3tk
```
