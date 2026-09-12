import json, pathlib, os, contextlib, ast, hashlib
import numpy as np
import pandas as pd
from astropy.io import fits

root = pathlib.Path('/Users/Igniz/Desktop/ICRAR/further')
out = pathlib.Path('/private/tmp/mauve_20260911')
os.chdir(root)
p = root / '20260909_check_corrected_Halpha_surface_density_SF_NSF_ND_by_stage.ipynb'
nb = json.loads(p.read_text())
ns = {}
checks = []
with (out / 'live_halpha_execution.log').open('w') as log, contextlib.redirect_stdout(log):
    for cell in [3, 5, 7, 9]:
        exec(compile(''.join(nb['cells'][cell-1]['source']), f'{p.name}:cell{cell}', 'exec'), ns)
    original_loader = ns['load_halpha_category_maps']
    def audited_loader(row, geometry):
        maps = original_loader(row, geometry)
        with fits.open(row.sfr_path, memmap=True) as hd:
            ha = np.asarray(hd['HA6562_FLUX_CORR'].data, dtype=float)
            hb = np.asarray(hd['HB4861_FLUX_CORR'].data, dtype=float)
            raw_ha = np.asarray(hd['HA6562_FLUX'].data, dtype=float)
            raw_hb = np.asarray(hd['HB4861_FLUX'].data, dtype=float)
            with np.errstate(divide='ignore', invalid='ignore'):
                raw_bd = raw_ha / raw_hb
                log_hb = maps['log_lha_surface'] + np.log10(hb / ha)
            common = maps['bpt_common_valid']
            maps['log_bpt_line_surface']['HB4861'] = np.where(common, log_hb, np.nan)
            keep = maps['valid_disc_mask'] & np.isfinite(maps['snr_postfit']) & (maps['snr_postfit'] > 25)
            for cat in ['sf', 'nsf']:
                for footprint in ['all_retained', 'mass_window', 'radial_window', 'common_bpt_mass', 'common_bpt_radius']:
                    use = keep & maps[cat + '_mask']
                    if 'mass' in footprint:
                        use &= (maps['log_sigma_star'] >= 7) & (maps['log_sigma_star'] < 9.5)
                    if 'radi' in footprint:
                        use &= (maps['radius_re'] >= 0) & (maps['radius_re'] < 2.5)
                    if 'common_bpt' in footprint:
                        use &= common
                    checks.append({'GALID': row.GALID, 'stage': row.stage, 'category': cat,
                                   'footprint': footprint, 'n': int(use.sum()),
                                   'n_raw_below_286': int((use & (raw_bd < 2.86)).sum())})
        return maps
    ns['load_halpha_category_maps'] = audited_loader
    ns['BPT_LINE_ORDER'] = ns['BPT_LINE_ORDER'] + ('HB4861',)
    for cell in [11, 13]:
        exec(compile(''.join(nb['cells'][cell-1]['source']), f'{p.name}:cell{cell}', 'exec'), ns)
    # Extract only the unchanged central-estimator functions, omitting plot/bootstrap calls.
    for cell, names in [(22, ['stage_bpt_line_estimator']),
                        (25, ['add_occupancy_weighted_halpha', 'stage_j_estimator']),
                        (28, ['build_halpha_share_table', 'stage_q_estimator'])]:
        tree = ast.parse(''.join(nb['cells'][cell-1]['source']))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in names:
                exec(compile(ast.Module(body=[node], type_ignores=[]), f'{p.name}:cell{cell}:function', 'exec'), ns)
    frame = ns['HALPHA_BY_GALAXY_BIN']
    scalar = frame[[c for c in frame if not c.startswith('values_') and not c.startswith('log_lha_values')]]
    scalar.to_csv(out / 'halpha_by_galaxy_bin.csv', index=False)
    ns['STAGE_HALPHA_PROFILES'].to_csv(out / 'stage_halpha_profiles.csv', index=False)
    ns['stage_bpt_line_estimator'](frame).to_csv(out / 'stage_bpt_line_profiles_with_hbeta.csv', index=False)
    j = ns['add_occupancy_weighted_halpha'](frame)
    ns['stage_j_estimator'](j).to_csv(out / 'stage_j_profiles.csv', index=False)
    pd.DataFrame(checks).to_csv(out / 'balmer_common_support_check.csv', index=False)
    # An independent check of the actual per-galaxy category partition.
    eligible = scalar[scalar.eligible_usable]
    sums = eligible.groupby(['GALID', 'variable', 'bin_left']).f_category.sum()
    assert np.allclose(sums, 1)
    assert ns['all_retained_strict']
    checks_out = {'strict_snr': True, 'category_partition_max_residual': float(np.max(abs(sums-1))),
                  'sample_by_stage': ns['sample'].groupby('stage').size().to_dict(),
                  'notebook_sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                  'scope': 'Selected data cells and central estimator functions; extra direct Hbeta diagnostic; no bootstrap rerun or notebook edits.'}
    (out / 'live_halpha_checks.json').write_text(json.dumps(checks_out, indent=2))
print('HALPHA_SELECTED_CELLS_AND_DIRECT_HBETA_PASS', flush=True)
