"""Read small spectral samples; emit compressed maps to stdout, no remote writes."""
import sys
import numpy as np
from astropy.io import fits

root = '/scratch/pawsey1308/mauve/cubes/v3tk/'
out = {}
for gal in ['NGC4689', 'IC3392', 'NGC4606']:
    path = root + gal + '_DATACUBE_FINAL_WCS_Pall_mad_red_v3tk.fits'
    with fits.open(path, memmap=True) as h:
        hdr = h['DATA'].header
        step = hdr.get('CD3_3', hdr.get('CDELT3'))
        wave = hdr['CRVAL3'] + (np.arange(h['DATA'].shape[0]) + 1 - hdr['CRPIX3']) * step
        out[gal+'_header'] = np.array(hdr.tostring())
        out[gal+'_primary'] = np.array(h[0].header.tostring())
        for center in [5100, 5500, 6000, 6800]:
            idx = np.argsort(abs(wave-center))[:9]
            flux = np.array(h['DATA'].data[idx], dtype=float)
            var = np.array(h['STAT'].data[idx], dtype=float)
            good = np.isfinite(flux) & np.isfinite(var) & (var > 0)
            if 'DQ' in h:
                good &= h['DQ'].data[idx] == 0
            flux[~good] = np.nan
            var[~good] = np.nan
            out[f'{gal}_flux_{center}'] = np.nanmedian(flux, axis=0).astype('f4')
            out[f'{gal}_var_{center}'] = np.nanmedian(var, axis=0).astype('f4')
        print(gal, h['DATA'].shape, 'NFILES', h[0].header.get('NFILES'), file=sys.stderr)
np.savez_compressed(sys.stdout.buffer, **out)
