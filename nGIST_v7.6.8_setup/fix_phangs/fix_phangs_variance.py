#!/usr/bin/env python3
"""Write PHANGS native fixed cubes with fillable STAT gaps patched.

This uses the same target-gap selection as check_phangs_variance.py. For every
fillable STAT gap, it copies the original cube to
{GALID}_PHANGS_DATACUBE_native_fixed.fits and fills the bad STAT samples with
the mean of the available immediate positive bounds. Either one or both sides
may provide a bound. Where DATA and STAT are both NaN, it fills both from their
separate spatial means within a 1 arcsec-diameter circle at the same wavelength,
but only when the complete 8-connected NaN component is surrounded by finite
DATA and positive finite STAT pixels. The input cube is not modified, and no
fixed cube is written when no fillable gaps are found.
"""

from phangs_variance_tools import fix_main


if __name__ == "__main__":
    raise SystemExit(fix_main())
