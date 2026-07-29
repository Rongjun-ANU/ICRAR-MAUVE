#!/usr/bin/env python3
"""Check PHANGS native cube STAT HDUs for fillable variance gaps.

By default this scans NGC4254, NGC4321, and NGC4535 under
/scratch/pawsey1308/mauve/cubes/v3tk, or one requested GALID if supplied.
The checker reads each {GALID}_PHANGS_DATACUBE_native.fits STAT HDU, limits the
search to the open 4700-9350 A interval, skips the AO/LGS gap, skips
{GALID}_mask.fits masked spaxels, and requires positive DATA for spectral
STAT-gap selection. Joint DATA/STAT NaNs are grouped into 8-connected components
on each wavelength plane. A component is spatially fillable only when its full
adjacent boundary has finite DATA and positive finite STAT. The log is replaced
at the start of each run after any previous log is rotated to a timestamped
filename, and it lists only fillable targets.
"""

from phangs_variance_tools import check_main


if __name__ == "__main__":
    raise SystemExit(check_main())
