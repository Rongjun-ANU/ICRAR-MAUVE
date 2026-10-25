from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.spatial import ConvexHull
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks
import json
from astropy.io import fits

ROOT = Path(__file__).parent
d = np.load(ROOT/'sampled_maps.npz')
fig, axs = plt.subplots(3, 4, figsize=(15, 12), constrained_layout=True)
for row, gal in enumerate(['NGC4689', 'IC3392', 'NGC4606']):
    for col, w in enumerate([5100, 5500, 6000, 6800]):
        v = d[f'{gal}_var_{w}']
        z = np.log10(np.where(v > 0, v, np.nan))
        valid = np.isfinite(z)
        print(gal, w, v.shape, 'finite', valid.sum(), 'variance_percentiles', np.nanpercentile(v, [5, 50, 95]))
        ax = axs[row, col]
        lo, hi = np.nanpercentile(z, [5, 95])
        im = ax.imshow(z, origin='lower', vmin=lo, vmax=hi, cmap='viridis')
        ax.set_title(f'{gal}: log variance, {w} A')
        fig.colorbar(im, ax=ax, shrink=.7)
fig.savefig(ROOT/'variance_windows.png', dpi=150)
plt.close(fig)

# Fit the orientation from the external valid-data hull, independently of
# the internal variance troughs. Coordinates here are zero-based pixels.
gal = 'NGC4689'
z = np.log10(d[f'{gal}_var_5500'])
valid = np.isfinite(z)
y, x = np.indices(z.shape)
points = np.column_stack([x[valid], y[valid]])
hull = points[ConvexHull(points).vertices]
def transform(xx, yy, angle):
    t = np.deg2rad(angle)
    return xx*np.cos(t)+yy*np.sin(t), -xx*np.sin(t)+yy*np.cos(t)
def area(angle):
    a, b = transform(hull[:, 0], hull[:, 1], angle)
    return np.ptp(a)*np.ptp(b)
angle = minimize_scalar(area, bounds=(-20, 0), method='bounded').x
u, v = transform(x, y, angle)
bounds = [float(np.min(u[valid])), float(np.max(u[valid])), float(np.min(v[valid])), float(np.max(v[valid]))]
profiles = {}
troughs = {}
for w in [5100, 5500, 6000, 6800]:
    zz = np.log10(d[f'{gal}_var_{w}'])
    profiles[w] = []
    troughs[w] = []
    for axis, coord, other, blo, bhi, olo, ohi, expected in [
        ('u', u, v, *bounds[:2], *bounds[2:], 1),
        ('v', v, u, *bounds[2:], *bounds[:2], 2),
    ]:
        bins = np.arange(np.floor(blo), np.ceil(bhi)+1)
        prof = np.full(len(bins)-1, np.nan)
        select = np.isfinite(zz) & (other>olo+30) & (other<ohi-30)
        ii = np.floor(coord).astype(int)
        for k, lo in enumerate(bins[:-1]):
            vals = zz[select & (ii==int(lo))]
            if len(vals)>20:
                prof[k] = np.median(vals)
        finite = np.isfinite(prof)
        prof = np.interp(np.arange(len(prof)), np.where(finite)[0], prof[finite])
        smooth = gaussian_filter(prof, 2)
        # Subtract a broad local trend to isolate narrow overlap depressions.
        residual = smooth-gaussian_filter(prof, 35)
        peaks, props = find_peaks(-residual, prominence=.07, distance=100)
        peaks = [p for p in peaks if 40<p<len(prof)-40]
        peaks = sorted(sorted(peaks, key=lambda p: residual[p])[:expected])
        assert len(peaks)==expected, (w, axis, peaks)
        centers = [float(bins[p]+.5) for p in peaks]
        profiles[w].append((bins[:-1]+.5, smooth, residual))
        troughs[w].append(centers)

uc = np.median([troughs[w][0][0] for w in troughs])
vc = np.median([troughs[w][1] for w in troughs], axis=0)
# Estimate strip extent from half-depth crossings in the detrended profiles.
# Use the union over wavelength to avoid assigning overlaps exclusively.
band_samples = [[], [], []]
for w in profiles:
    for axis in range(2):
        coords, smooth, residual = profiles[w][axis]
        for j, center in enumerate(troughs[w][axis]):
            k = int(np.argmin(abs(coords-center)))
            threshold = residual[k]/2
            left = right = k
            while left>0 and residual[left-1]<threshold: left-=1
            while right<len(coords)-1 and residual[right+1]<threshold: right+=1
            band_samples[0 if axis==0 else j+1].append([float(coords[left]-.5),float(coords[right]+.5)])
bands = [[min(a[0] for a in s),max(a[1] for a in s)] for s in band_samples]
footprints = np.zeros((6,)+z.shape, dtype=np.uint8)
for row in range(3):
    vlo = bounds[2] if row==0 else bands[row][0]
    vhi = bounds[3] if row==2 else bands[row+1][1]
    for col in range(2):
        ulo = bounds[0] if col==0 else bands[0][0]
        uhi = bands[0][1] if col==0 else bounds[1]
        footprints[row*2+col] = valid & (u>=ulo) & (u<=uhi) & (v>=vlo) & (v<=vhi)
# Conservative 20-pixel exclusions on either side of each trough center;
# these are a chosen safety margin, not measured footprint uncertainties.
margin = 20
labels = np.zeros(z.shape, dtype=np.int16)
for row in range(3):
    vlo, vhi = [bounds[2], *vc, bounds[3]][row:row+2]
    for col in range(2):
        ulo, uhi = [bounds[0], uc, bounds[1]][col:col+2]
        sel = valid & (u>ulo+margin) & (u<uhi-margin) & (v>vlo+margin) & (v<vhi-margin)
        labels[sel] = row*2+col+1
hdr = fits.Header.fromstring(str(d[gal+'_header']))
# Retain only celestial WCS using Astropy's celestial projection.
from astropy.wcs import WCS
hdr2 = WCS(hdr).celestial.to_header()
hdr2['BUNIT'] = 'region ID'
hdr2['COMMENT'] = 'Approximate interior regions inferred from variance troughs; not OB IDs.'
hdr2['COMMENT'] = '0 excludes overlaps, edges and invalid pixels; 1-6 bottom to top, left to right.'
fits.PrimaryHDU(labels, header=hdr2).writeto(ROOT/'NGC4689_approximate_interiors.fits', overwrite=True)
hdus = [fits.PrimaryHDU()]
for i in range(6):
    fh = WCS(hdr).celestial.to_header()
    fh['COMMENT'] = 'Inferred footprint; overlaps belong to multiple regions. Not verified OB IDs.'
    hdus.append(fits.ImageHDU(footprints[i], header=fh, name=f'R{i+1}'))
hdus.append(fits.ImageHDU(footprints.sum(axis=0).astype('u1'), header=WCS(hdr).celestial.to_header(), name='MEMBERS'))
fits.HDUList(hdus).writeto(ROOT/'NGC4689_approximate_shared_footprints.fits', overwrite=True)
summary = dict(angle_degrees=float(angle), uv_bounds=bounds,
               transform='u=x*cos(angle)+y*sin(angle); v=-x*sin(angle)+y*cos(angle); zero-based x,y',
               trough_centers_by_wavelength=troughs, exclusion_halfwidth_pixels=margin,
               region_counts={str(i):int((labels==i).sum()) for i in range(1,7)},
               caution='Regions are geometric labels, not verified OB001-OB006 identities or unblended exposure spectra.')
summary['approximate_overlap_bands_uv'] = bands
summary['overlap_method'] = 'Union of half-depth detrended variance trough widths across four windows; approximate, not measured exposure support.'
summary['single_pointing_hull_sizes_pixels'] = {}
for control in ['IC3392', 'NGC4606']:
    cv = np.isfinite(d[control+'_var_5500'])
    yy, xx = np.where(cv)
    pp = np.column_stack([xx, yy])
    hh = pp[ConvexHull(pp).vertices]
    def control_area(a):
        aa, bb = transform(hh[:,0], hh[:,1], a)
        return np.ptp(aa)*np.ptp(bb)
    ca = minimize_scalar(control_area, bounds=(20,60), method='bounded').x
    aa, bb = transform(hh[:,0], hh[:,1], ca)
    summary['single_pointing_hull_sizes_pixels'][control] = [float(np.ptp(aa)),float(np.ptp(bb))]
(ROOT/'geometry.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))

fig, axs = plt.subplots(2, 2, figsize=(12, 11), constrained_layout=True)
im=axs[0,0].imshow(z, origin='lower', cmap='viridis', vmin=np.nanpercentile(z,5), vmax=np.nanpercentile(z,95))
fig.colorbar(im, ax=axs[0,0], shrink=.7, label='log10 propagated variance')
axs[0,0].contour(u, levels=[uc], colors='white', linewidths=1)
axs[0,0].contour(v, levels=vc, colors='white', linewidths=1)
axs[0,0].set_title('NGC4689: 5500 A variance + inferred strip centers')
axs[0,1].imshow(np.ma.masked_equal(labels,0), origin='lower', cmap='tab10', vmin=1, vmax=10)
for i in range(1,7):
    yy,xx=np.where(labels==i)
    axs[0,1].text(np.median(xx), np.median(yy), f'R{i}', ha='center', va='center', fontsize=15)
axs[0,1].set_title('Approximate interiors; overlaps and edges excluded')
for w in profiles:
    for j in range(2):
        coords, smooth, residual=profiles[w][j]
        axs[1,j].plot(coords, residual, label=f'{w} A', lw=1)
for j, centers in enumerate([[uc],vc]):
    for c in centers: axs[1,j].axvline(c,color='k',ls='--',lw=.8)
    axs[1,j].set_xlabel(('u' if j==0 else 'v')+' coordinate (pixels)')
    axs[1,j].set_ylabel('log variance minus broad spatial trend')
    axs[1,j].set_ylim(-.25,.12)
    axs[1,j].set_xlim((bounds[:2] if j==0 else bounds[2:])[0]+25, (bounds[:2] if j==0 else bounds[2:])[1]-25)
    axs[1,j].legend()
    axs[1,j].set_title('Repeated overlap troughs across wavelength')
fig.savefig(ROOT/'NGC4689_inferred_pointings.png', dpi=160)
plt.close(fig)

fig, axs = plt.subplots(1,2,figsize=(11,7),constrained_layout=True)
axs[0].imshow(z,origin='lower',cmap='gray',vmin=np.nanpercentile(z,5),vmax=np.nanpercentile(z,95))
for i in range(6):
    axs[0].contour(footprints[i],levels=[.5],colors=[plt.get_cmap('tab10')(i)],linewidths=1.2)
    yy,xx=np.where(labels==i+1)
    axs[0].text(np.median(xx),np.median(yy),f'R{i+1}',ha='center',color='white')
axs[0].set_title('Approximate footprints: shared overlaps')
im=axs[1].imshow(np.ma.masked_equal(footprints.sum(axis=0),0),origin='lower',cmap='viridis',vmin=1,vmax=4)
fig.colorbar(im,ax=axs[1],ticks=[1,2,4],shrink=.6,label='Number of candidate pointing memberships')
axs[1].set_title('Overlap strips belong to 2; crossings to 4')
fig.savefig(ROOT/'NGC4689_shared_footprints.png',dpi=150)
plt.close(fig)
