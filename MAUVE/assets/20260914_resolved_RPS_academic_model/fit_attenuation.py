"""Partial fit of the observed survivor-intensity contrast; galaxies are resampled."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
df = pd.read_csv(OUT/'sfr_by_galaxy_bin.csv')
df = df[df.variable.eq('log_sigma_star') & df.center.between(7.625,9.125)]
centres = np.sort(df.center.unique())
stages = ['pre-peak','close-to-peak','post-peak']
NBOOT = 10000
rng = np.random.default_rng(20260914)
matrices, boot, central = {}, {}, {}
for stage in stages:
    sub = df[df.stage.eq(stage)].copy()
    ids = sorted(sub.GALID.unique())
    sub.loc[~sub.eligible_usable, ['f_SF','mean_Sigma_SFR_all_HII_galaxy']] = np.nan
    f = sub.pivot(index='GALID',columns='center',values='f_SF').reindex(index=ids,columns=centres).to_numpy()
    t = sub.pivot(index='GALID',columns='center',values='mean_Sigma_SFR_all_HII_galaxy').reindex(index=ids,columns=centres).to_numpy()
    assert np.array_equal(np.isfinite(f),np.isfinite(t))
    matrices[stage] = (ids,f,t)
    # The same complete galaxy draw is used for every bin and numerator/denominator.
    draw = rng.integers(0,len(ids),size=(NBOOT,len(ids)))
    central[stage] = np.nansum(t,axis=0)/np.nansum(f,axis=0)
    numerator = np.nansum(t[draw],axis=1)
    denominator = np.nansum(f[draw],axis=1)
    boot[stage] = np.divide(numerator,denominator,
                           out=np.full(numerator.shape,np.nan),where=denominator>0)

rows, results, loo_rows = [], [], []
fig, axes = plt.subplots(1,2,figsize=(10.8,4.1),layout='constrained')
for ax,stage,color in zip(axes,stages[1:],['#1a9850','#d73027']):
    y = np.log10(central[stage]/central['pre-peak'])
    yb = np.log10(boot[stage]/boot['pre-peak'])
    valid = np.isfinite(yb).all(axis=1)
    yb = yb[valid]
    errors = np.std(yb,axis=0,ddof=1)
    weights = 1/errors**2
    x = centres-8.5
    X = np.column_stack([np.ones(len(x)),x])
    transform = np.linalg.solve(X.T@(weights[:,None]*X), X.T*weights)
    pars = transform@y
    draws = np.sum(yb[:,:,None]*transform.T[None,:,:],axis=1)
    constant = float(np.sum(weights*y)/weights.sum())
    cb = np.sum(yb*weights[None,:],axis=1)/weights.sum()
    assert np.isfinite(draws).all() and np.isfinite(cb).all()
    q = np.percentile(cb,[16,50,84])
    linear_q = np.percentile(draws,[16,50,84],axis=0)
    fit = X@pars
    attenuation = 10**constant
    result = dict(stage=stage,bins=centres.tolist(),n_bootstrap=NBOOT,
        finite_complete_draws=int(valid.sum()),constant_dex=constant,
        constant_p16=float(q[0]),constant_p84=float(q[2]),
        attenuation_factor=attenuation,
        factor_p16=float(10**q[0]),factor_p84=float(10**q[2]),
        effective_exposure=float(-np.log(10)*constant),
        exposure_p16=float(-np.log(10)*q[2]),exposure_p84=float(-np.log(10)*q[0]),
        slope_dex_per_dex=float(pars[1]),slope_p16=float(linear_q[0,1]),slope_p84=float(linear_q[2,1]),
        intercept_linear=float(pars[0]),
        rms_constant_dex=float(np.sqrt(np.mean((y-constant)**2))),
        rms_linear_dex=float(np.sqrt(np.mean((y-fit)**2))),
        equal_bin_constant_dex=float(np.mean(y)),
        max_abs_constant_residual_dex=float(np.max(abs(y-constant))),
        bootstrap_correlation_min=float(np.min(np.corrcoef(yb.T)[np.triu_indices(len(x),1)])),
        bootstrap_correlation_max=float(np.max(np.corrcoef(yb.T)[np.triu_indices(len(x),1)])))
    # Leave one whole galaxy out, recomputing the profile with fixed full-fit weights.
    for omit_stage in ['pre-peak',stage]:
        ids,f,t = matrices[omit_stage]
        for j,gid in enumerate(ids):
            keep = np.arange(len(ids))!=j
            estimate = np.nansum(t[keep],axis=0)/np.nansum(f[keep],axis=0)
            yp = np.log10(estimate/central['pre-peak']) if omit_stage==stage else np.log10(central[stage]/estimate)
            a = float(np.sum(weights*yp)/weights.sum())
            loo_rows.append(dict(target_stage=stage,omitted_stage=omit_stage,GALID=gid,constant_dex=a))
    loo = [r['constant_dex'] for r in loo_rows if r['target_stage']==stage]
    result.update(leave_one_galaxy_min_dex=min(loo),leave_one_galaxy_max_dex=max(loo))
    results.append(result)
    lo,hi=np.percentile(yb,[16,84],axis=0)
    ax.errorbar(centres,y,yerr=[y-lo,hi-y],fmt='o',color=color,capsize=3,label='Observed contrast')
    ax.axhline(constant,color=color,label='Constant attenuation fit')
    ax.fill_between(centres,q[0],q[2],color=color,alpha=.14,label='Galaxy-bootstrap 16-84%')
    ax.plot(centres,fit,color='0.25',ls='--',label='Linear-density sensitivity')
    ax.axhline(0,color='0.5',lw=.8)
    ax.set(xlabel=r'$\log_{10}(\Sigma_*/M_\odot\,{\rm kpc}^{-2})$',
           ylabel=r'$\log_{10}(I_{\rm SF,stage}/I_{\rm SF,pre})$',title=stage,ylim=(-.85,.55))
    ax.legend(fontsize=7.5)
    for j,c in enumerate(centres):
        rows.append(dict(stage=stage,center=c,observed_delta=y[j],p16=lo[j],p84=hi[j],
                         sigma_boot=errors[j],constant_prediction=constant,linear_prediction=fit[j],
                         n_pre_usable=int(np.isfinite(matrices['pre-peak'][1][:,j]).sum()),
                         n_target_usable=int(np.isfinite(matrices[stage][1][:,j]).sum())))
pd.DataFrame(rows).to_csv(OUT/'attenuation_fit_profiles.csv',index=False)
pd.DataFrame(loo_rows).to_csv(OUT/'attenuation_leave_one_galaxy.csv',index=False)
(OUT/'attenuation_fit_results.json').write_text(json.dumps(results,indent=2))
fig.savefig(OUT/'figure_07_attenuation_fit.png',dpi=200,bbox_inches='tight')
fig.savefig(OUT/'figure_07_attenuation_fit.pdf',bbox_inches='tight')
plt.close(fig)
assert all(r['finite_complete_draws']>9900 for r in results)
assert all(np.isfinite(r['constant_dex']) for r in results)
print(json.dumps(results,indent=2))
