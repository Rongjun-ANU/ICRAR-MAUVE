#!/usr/bin/env python3
"""Write PHANGS native fixed cubes with fillable STAT gaps patched.

This uses the same target-gap selection as check_phangs_variance.py. For every
fillable STAT gap, it copies the original cube to
{GALID}_PHANGS_DATACUBE_native_fixed.fits and fills the bad STAT samples with
the nanmean of the immediate bracketing STAT values. The input cube is not
modified, and no fixed cube is written when no fillable gaps are found.
"""

from phangs_variance_tools import fix_main


if __name__ == "__main__":
    raise SystemExit(fix_main())
