"""Read-only full-template comparison to unmodified nGIST 7.6.8 functions."""
import ast
import copy
import json
from pathlib import Path
import sys
import logging
import argparse
import hashlib
import re
import yaml
import scipy
import astropy
import ppxf as ppxf_package
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d
from ppxf.ppxf_util import gaussian_filter1d, log_rebin

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--weights', type=Path, required=True)
parser.add_argument('--ngist-source-dir', type=Path, required=True,
    help='Directory containing miles.py and _auxiliary.py from the pinned nGIST revision')
parser.add_argument('--output-dir', type=Path, required=True)
args=parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)
# Exercise the exact functions now implemented inside Mass.py, without running
# its module-level CLI/file-loading workflow. No runtime normalization module.
mass_path=ROOT/'Mass.py'
mass_tree=ast.parse(mass_path.read_text())
function_names={'sha256_file','parse_template_filename_key','_resource_path',
                'load_run_setup','load_lsf','build_template_norm_flux_lookup',
                'r_light_mass_to_light'}
constant_names={'NGIST_REVISION','TEMPLATE_VELSCALE_RATIO','_TEMPLATE_RE'}
mass_nodes=[node for node in mass_tree.body
    if isinstance(node,ast.FunctionDef) and node.name in function_names
    or isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id in constant_names for t in node.targets)]
mass_namespace=dict(globals(),__file__=str(mass_path))
exec(compile(ast.Module(body=mass_nodes,type_ignores=[]),str(mass_path),'exec'),mass_namespace)
globals().update({name:mass_namespace[name] for name in function_names})

class QuietStatus:
    @staticmethod
    def running(*args): pass
    @staticmethod
    def updateDone(*args): pass
    @staticmethod
    def updateFailed(*args): pass

ns=dict(np=np,fits=fits,interp1d=interp1d,gaussian_filter1d=gaussian_filter1d,
        log_rebin=log_rebin,printStatus=QuietStatus,logging=logging)
import os,glob
ns.update(os=os,glob=glob)
for file,names in [('miles.py',{'age_metal_alpha','prepareSpectralTemplateLibrary'}),('_auxiliary.py',{'getLSF'})]:
    path=args.ngist_source_dir/file
    tree=ast.parse(path.read_text())
    functions=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in names]
    assert len(functions)==len(names)
    exec(compile(ast.Module(body=functions,type_ignores=[]),str(path),'exec'),ns)

weights_path=args.weights.resolve()
with fits.open(weights_path) as h:
    w=h[1].data['WEIGHTS'].astype(float); grid=h[2].data.copy(); hdr=h[1].header.copy()
keys=[(round(10**float(row['LOGAGE']),2),round(float(row['METAL']),2),round(float(row['ALPHA']),2)) for row in grid]
setup=load_run_setup(weights_path,hdr,ROOT)
data,temp=load_lsf(setup['data_lsf_path'],setup['template_lsf_path'],setup['redshift'])
lookup,meta=build_template_norm_flux_lookup(ROOT/'EMILES_baseFe_safe',setup['lmin'],setup['lmax'],lsf_data=data,lsf_templates=temp,template_velscale=setup['template_velscale'],grid_keys=keys)
factors=np.array([lookup[key] for key in keys])
cfg=copy.deepcopy(setup['config'])
cfg['GENERAL']['TEMPLATE_DIR']=str(ROOT)
cfg['GENERAL']['CONFIG_DIR']=str(ROOT)
refdata,reftemp=ns['getLSF'](cfg,'SFH')
sample=np.linspace(4600,9300,100)
np.testing.assert_array_equal(data(sample),refdata(sample))
np.testing.assert_array_equal(temp(sample),reftemp(sample))
reference=ns['prepareSpectralTemplateLibrary'](cfg,setup['lmin'],setup['lmax'],setup['template_velscale'],refdata,reftemp,'SFH',True)
light=reference[0].reshape(reference[0].shape[0],-1)
refkeys=[(round(10**float(a),2),round(float(m),2),round(float(z),2)) for a,m,z in zip(reference[4].ravel(),reference[5].ravel(),reference[6].ravel())]
order=[refkeys.index(k) for k in keys]
cfg['SFH']['NORM_TEMP']='UNNORMALIZED_VALIDATION_ONLY'
raw_result=ns['prepareSpectralTemplateLibrary'](cfg,setup['lmin'],setup['lmax'],setup['template_velscale'],refdata,reftemp,'SFH',True)
raw=raw_result[0].reshape(raw_result[0].shape[0],-1)[:,order]
reference_means=raw.mean(axis=0)
np.testing.assert_allclose(factors,reference_means,rtol=2e-13,atol=0)
np.testing.assert_allclose(raw/factors[None,:],light[:,order],rtol=2e-13,atol=2e-13)

names='IMF slope MH Age U B V R I J H K UminusV BminusV VminusR VminusI VminusJ VminusH VminusK ML_U ML_B ML_V ML_R ML_I ML_J ML_H ML_K F439W F555W F675W F814W C439_555 C555_675 C555_814'.split()
tbl=np.genfromtxt(ROOT/'BaSTI+Chabrier.dat',dtype=None,encoding=None,names=names,comments='#',skip_header=2,autostrip=True)
phot={(round(float(r['Age']),2),round(float(r['MH']),2)):r for r in tbl if r['IMF']=='Ch'}
ml=np.array([phot[k[:2]]['ML_R'] for k in keys]); lum=np.array([10**(-.4*phot[k[:2]]['R']) for k in keys])
new,frac=r_light_mass_to_light(w,ml,lum,factors)
light_r=np.sum((w/w.sum(axis=1)[:,None])*ml[None,:],axis=1)
# Positive values mean that the unconverted light_r assumption overestimates
# M/L relative to the corrected interpretation for this same saved mixture.
light_r_error_dex=np.log10(light_r/new)
light_r_error_percent=100.*(light_r/new-1.)
old_lookup={}
for p in sorted((ROOT/'EMILES_baseFe_safe').glob('*.fits')):
    with fits.open(p) as h:
        wave=h[0].header['CRVAL1']+h[0].header['CDELT1']*np.arange(h[0].header['NAXIS1'])
        old_lookup[parse_template_filename_key(p)]=np.nanmean(h[0].data[(wave>=setup['lmin'])&(wave<=setup['lmax'])])
oldf=np.array([old_lookup[k] for k in keys])
old,_=r_light_mass_to_light(w,ml,lum,oldf)
delta=np.log10(new/old)
summary=dict(weights=str(weights_path), weights_sha256=sha256_file(weights_path),
    photometry_sha256=sha256_file(ROOT/'BaSTI+Chabrier.dat'),
    reference_source_sha256={name:sha256_file(args.ngist_source_dir/name) for name in ('miles.py','_auxiliary.py')},bins=w.shape[0],ssps=len(keys),setup={k:v for k,v in setup.items() if k!='config'},metadata=meta,
    max_factor_relative_error=float(np.max(np.abs(factors/reference_means-1))),
    max_normalized_spectrum_absolute_error=float(np.max(np.abs(raw/factors[None,:]-light[:,order]))),
    factor_new_over_old_percentiles=np.percentile(factors/oldf,[0,16,50,84,100]).tolist(),
    delta_log10_ml_percentiles=np.percentile(delta,[0,16,50,84,100]).tolist(),
    light_r_error_definition='log10(light_r / corrected_light_norm_to_r)',
    light_r_error_dex_percentiles=np.percentile(light_r_error_dex,[0,16,50,84,100]).tolist(),
    light_r_error_percent_percentiles=np.percentile(light_r_error_percent,[0,16,50,84,100]).tolist(),
    light_r_overestimate_bins=int(np.count_nonzero(light_r_error_dex>0)),
    light_r_underestimate_bins=int(np.count_nonzero(light_r_error_dex<0)),
    median_light_r_ml=float(np.median(light_r)),
    max_fraction_sum_error=float(np.max(np.abs(frac.sum(axis=1)-1))),
    median_old_ml=float(np.median(old)),median_new_ml=float(np.median(new)))
out=args.output_dir/f'validation_{int(setup["lmax"])}.json'
out.write_text(json.dumps(summary,indent=2)+'\n')
np.savez(args.output_dir/f'comparison_{int(setup["lmax"])}.npz',new=new,old=old,delta=delta,
         factors=factors,old_factors=oldf,light_r=light_r,light_r_error_dex=light_r_error_dex,
         light_r_error_percent=light_r_error_percent)
print(json.dumps(summary,indent=2))

