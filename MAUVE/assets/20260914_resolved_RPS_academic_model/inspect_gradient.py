"""Execute data/plane cells of the exact requested notebook; save report extracts."""
from pathlib import Path
import contextlib
import hashlib
import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
ROOT = Path('/Users/Igniz/Desktop/ICRAR/further')
PATH = ROOT / '20260914_check_SF_gradient_scan_VIVA.ipynb'
nb = json.loads(PATH.read_text())
os.chdir(ROOT)
ns = {'hashlib': hashlib}
cells = [3, 5, 9, 11, 13, 16, 20, 28, 30]
with (OUT / 'gradient_execution.log').open('w') as log, contextlib.redirect_stdout(log):
    for cell in cells:
        print('EXECUTING', cell, flush=True)
        exec(compile(''.join(nb['cells'][cell-1]['source']), f'{PATH.name}:cell{cell}', 'exec'), ns)
for key in ['gradient_summary', 'preferred_sides', 'alignment', 'map_support',
            'prepeak_reference', 'control_scatter', 'galaxy_metrics']:
    ns[key].to_csv(OUT / (key + '.csv'), index=False)
g = 'NGC4654'
m = ns['spatial'][g]
d = ns['direction_inputs'][g]
side = ns['preferred_sides'].set_index('GALID').loc[g]
fixed = ns['sides_at_angle'](d['pa'], d['values'], 315.0)
result = dict(notebook=str(PATH), sha256=hashlib.sha256(PATH.read_bytes()).hexdigest(),
              executed_cells=cells, named_galaxy=g, inferred_sides=side.to_dict(),
              fixed_VIVA_PA_sides=fixed,
              formal_PA_uncertainty_computed=False,
              continuum_SNR25_cut_in_gradient_notebook=False,
              scope='Data cells and plane estimator; original plotting/gallery cells not executed.')
(OUT / 'ngc4654_audit.json').write_text(json.dumps(result, indent=2, default=lambda x: x.item()))
fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.1), layout='constrained')
scale = float(m['components'].iloc[0].kpc_per_arcsec)
mask = m['joint_mask']
for ax, key, title in zip(axes[:2], ['delta_local', 'delta_env'],
                          ['Within-galaxy residual', 'Relative to pre-peak control']):
    arr = np.where(mask, m[key], np.nan)
    im = ax.pcolormesh(m['east']*scale, m['north']*scale, arr, shading='nearest',
                      cmap='RdBu_r', vmin=-.65, vmax=.65, rasterized=True)
    ax.set_aspect('equal'); ax.invert_xaxis()
    ax.set(xlabel='East offset (kpc)', ylabel='North offset (kpc)', title=title)
    for theta, color in [(float(side.theta_SF), '#00b8c4'), (315., '#db9d21')]:
        length = 3.0
        ax.annotate('', xy=(length*np.sin(np.deg2rad(theta)),length*np.cos(np.deg2rad(theta))),
                    xytext=(0,0), arrowprops=dict(arrowstyle='->', color=color, lw=2))
    fig.colorbar(im, ax=ax, shrink=.65, label='Residual (dex)')
ax = axes[2]
sep = ns['angle_distance'](d['pa'], 315.)
for select, label, color in [(sep<90, 'NW / VIVA-facing', '#2166ac'),
                             (sep>90, 'Opposite hemisphere', '#b2182b')]:
    values = np.sort(d['values'][select,1])
    ax.plot(values, np.arange(1,len(values)+1)/len(values), color=color, label=label)
ax.axvline(0, color='0.3', ls=':')
ax.set(xlabel='Environmental residual (dex)', ylabel='Cumulative area fraction',
       xlim=(-1.2,1.2), title='Fixed 315-degree hemispheres')
ax.legend(fontsize=8, loc='lower right')
fig.suptitle('NGC4654: descriptive SF asymmetry on the notebook joint SF mask', fontsize=12)
fig.savefig(OUT/'figure_04_ngc4654.png',dpi=200,bbox_inches='tight')
fig.savefig(OUT/'figure_04_ngc4654.pdf',bbox_inches='tight')
plt.close(fig)
print('NGC4654_LIVE_SELECTED_CELLS_PASS', json.dumps(result['inferred_sides'], default=str))
