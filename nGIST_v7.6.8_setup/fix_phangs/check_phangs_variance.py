#!/usr/bin/env python3
"""Check PHANGS native cube STAT HDUs for fillable variance gaps.

By default this scans NGC4254, NGC4321, and NGC4535 under
/scratch/pawsey1308/mauve/cubes/v3tk, or one requested GALID if supplied.
The checker reads each {GALID}_PHANGS_DATACUBE_native.fits STAT HDU, limits the
search to 4750-9100 A, skips the AO/LGS gap, skips {GALID}_mask.fits masked
spaxels, and skips spaxels whose DATA HDU is non-finite in the checked range.
Only fillable target gaps are printed and logged, with (x,y), z-index range,
wavelength range, length, and reason.
"""

from phangs_variance_tools import check_main


if __name__ == "__main__":
    raise SystemExit(check_main())
