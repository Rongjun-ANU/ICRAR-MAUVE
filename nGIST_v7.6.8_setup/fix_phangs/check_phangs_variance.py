#!/usr/bin/env python3
"""Check PHANGS native cube STAT HDUs for fillable variance gaps.

By default this scans NGC4254, NGC4321, and NGC4535 under
/scratch/pawsey1308/mauve/cubes/v3tk, or one requested GALID if supplied.
The checker reads each {GALID}_PHANGS_DATACUBE_native.fits STAT HDU, limits the
search to the open 4700-9350 A interval, skips the AO/LGS gap, skips
{GALID}_mask.fits masked spaxels, and excludes only wavelength samples where the
corresponding DATA value is non-finite from spectral STAT-gap selection. Samples
where both DATA and STAT are NaN are checked separately for spatial filling from
a 1 arcsec-diameter circle. Only fillable targets are printed and logged, with
coordinates, wavelength, method, and reason.
"""

from phangs_variance_tools import check_main


if __name__ == "__main__":
    raise SystemExit(check_main())
