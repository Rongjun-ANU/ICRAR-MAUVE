"""Read-only recalculation of central estimators, without bootstrap or fitting."""
from pathlib import Path
import json, os
import numpy as np
import pandas as pd

ROOT = Path('/Users/Igniz/Desktop/ICRAR/further')
OUT = Path(__file__).resolve().parent
os.chdir(ROOT)
nb = json.loads((ROOT/'20260909_check_Sigma_SFR_intensity_by_stage.ipynb').read_text())
ns = {'__name__': '__main__'}
for ci in (3, 5, 7):
    exec(compile(''.join(nb['cells'][ci-1]['source']), f'N103-cell-{ci}', 'exec'), ns)
    ns['display'] = lambda *args, **kwargs: None

original_loader = ns['load_galaxy_maps']
fraction_rows = []
def audited_loader(row, geometry):
    maps = original_loader(row, geometry)
    keep = maps['valid_disc_mask'] & np.isfinite(maps['snr_postfit']) & (maps['snr_postfit'] > 25)
    counts = {cat: int((keep & maps[key]).sum()) for cat,key in [('SF','sf_mask'),('ND','nd_mask'),('NSF','nsf_mask')]}
    assert sum(counts.values()) == int(keep.sum())
    fraction_rows.append(dict(GALID=row.GALID,stage=row.stage,N_valid=int(keep.sum()),**{'N_'+k:v for k,v in counts.items()},**{'f_'+k:v/keep.sum() for k,v in counts.items()}))
    return maps
ns['load_galaxy_maps'] = audited_loader
for ci in (9,11):
    exec(compile(''.join(nb['cells'][ci-1]['source']), f'N103-cell-{ci}', 'exec'), ns)
for name in ['SAMPLE_QC','SFR_BY_GALAXY_BIN','STAGE_SFR_PROFILES','SFR_SNR_SUMMARY']:
    ns[name].to_csv(OUT/('live_'+name+'.csv'),index=False)
fractions = pd.DataFrame(fraction_rows)
fractions.to_csv(OUT/'live_category_fractions_per_product.csv',index=False)
stage = fractions.groupby('stage')[['f_SF','f_ND','f_NSF']].mean().reindex(ns['STAGE_ORDER'])
stage.to_csv(OUT/'live_category_fractions_by_stage.csv')
print('LIVE EQUAL-GALAXY FULL-FOOTPRINT FRACTIONS\n',stage.to_string())
stored = pd.read_csv(OUT/'N104_c15_o4_t1.csv').set_index('stage').loc[list(ns['STAGE_ORDER']),['f_SF','f_ND','f_NSF']]
diff = float(np.max(np.abs(stage.to_numpy()-stored.to_numpy())))
print('MAX ABS DIFFERENCE FROM SAVED FOUR-DECIMAL FRACTIONS:',diff)
assert diff < 0.000051
profiles = ns['STAGE_SFR_PROFILES']
res=[]
for axis in ['log_sigma_star','radius_re']:
    ref=profiles[profiles.stage.eq('pre-peak') & profiles.variable.eq(axis)].set_index('bin_left')
    for _,row in profiles[profiles.variable.eq(axis) & ~profiles.stage.eq('pre-peak') & profiles.in_fixed_reliable_range].iterrows():
        r=dict(stage=row.stage,variable=axis,bin_left=row.bin_left,bin_right=row.bin_right,N_gal_usable=row.N_gal_usable,N_gal_SF_support=row.N_gal_SF_support)
        for metric in ['occupancy','survivor_intensity','total_hii']:
            a,b=row[metric],ref.loc[row.bin_left,metric]
            r['delta_log_'+metric]=float(np.log10(a/b)) if a>0 and b>0 else np.nan
        res.append(r)
contrasts=pd.DataFrame(res)
contrasts.to_csv(OUT/'live_sfr_decomposition.csv',index=False)
old=pd.read_csv(OUT/'N103_c13_o2_t1.csv')
join=contrasts.merge(old,on=['stage','variable','bin_left','bin_right'],suffixes=('_new','_saved'),validate='one_to_one')
diffs={}
for metric in ['occupancy','survivor_intensity','total_hii']:
    key='delta_log_'+metric
    diffs[key]=float(np.nanmax(np.abs(join[key+'_new']-join[key+'_saved'])))
    assert diffs[key]<0.000006
    assert np.array_equal(np.isfinite(join[key+'_new']),np.isfinite(join[key+'_saved']))
print('DECOMPOSITION MAX ABS DIFFERENCES:',diffs)
print('PASS: source setup/loader/central-estimator cells 3,5,7,9,11 only; no bootstrap, plots, raw fits, or product writes.')
(OUT/'central_result_check.json').write_text(json.dumps(dict(status='PASS',fraction_max_abs_difference=diff,decomposition_max_abs_differences=diffs,cells=[3,5,7,9,11],products=len(fractions)),indent=2))
