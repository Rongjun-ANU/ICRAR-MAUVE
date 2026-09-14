"""Report-only calculations from frozen scalar extracts; no FITS or notebook writes."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import ndtr
from scipy.integrate import quad, solve_ivp

OUT = Path(__file__).resolve().parent
STAGES = ['pre-peak', 'close-to-peak', 'post-peak']
COLORS = dict(zip(STAGES, ['#2166ac', '#1a9850', '#d73027']))
CC = {'SF': '#2166ac', 'NSF': '#b2182b', 'ND': '#e69f00'}
NBOOT = 10000
RNG = np.random.default_rng(20260914)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 200,
                     'font.family': 'DejaVu Sans'})
sfr = pd.read_csv(OUT / 'sfr_by_galaxy_bin.csv')
ha = pd.read_csv(OUT / 'halpha_by_galaxy_bin.csv')
ids = {s: sorted(sfr.loc[sfr.stage.eq(s), 'GALID'].unique()) for s in STAGES}
idx = {s: RNG.integers(0, len(ids[s]), (NBOOT, len(ids[s]))) for s in STAGES}

def mean_boot(v, s):
    a = np.asarray(v, float)[idx[s]]
    n = np.isfinite(a).sum(axis=1)
    return np.divide(np.nansum(a, axis=1), n, out=np.full(NBOOT, np.nan), where=n>0)

def avg(v):
    v = np.asarray(v, float)
    return float(np.mean(v[np.isfinite(v)])) if np.isfinite(v).any() else np.nan

def log(v):
    v = np.asarray(v, float)
    return np.log10(v, where=v>0, out=np.full(v.shape, np.nan))

def interval(v):
    a = v[np.isfinite(v)]
    return np.percentile(a, [16,84]) if len(a) else [np.nan,np.nan]

rows, draws = [], {}
for stage in STAGES:
    for variable in ['log_sigma_star', 'radius_re']:
        low, high = (7,9.5) if variable == 'log_sigma_star' else (0,2.5)
        centers = sorted(sfr.loc[sfr.variable.eq(variable) & sfr.center.ge(low) & sfr.center.lt(high), 'center'].unique())
        for center in centers:
            g = sfr[sfr.stage.eq(stage) & sfr.variable.eq(variable) & sfr.center.eq(center)].set_index('GALID').reindex(ids[stage])
            h = ha[ha.stage.eq(stage) & ha.variable.eq(variable) & ha.center.eq(center)]
            v = g.eligible_usable.fillna(False).astype(bool)
            o = g.f_SF.where(v).to_numpy(float)
            t = g.mean_Sigma_SFR_all_HII_galaxy.where(v).to_numpy(float)
            ob, tb = mean_boot(o,stage), mean_boot(t,stage)
            ib = np.divide(tb,ob,out=np.full(NBOOT,np.nan),where=ob>0)
            vals = {'O': (avg(o), ob, int(v.sum())),
                    'I': (avg(t)/avg(o) if avg(o)>0 else np.nan, ib, int((v & g.N_SF.ge(20)).sum())),
                    'T': (avg(t), tb, int(v.sum()))}
            for category in ['SF','NSF','ND']:
                c = h[h.category.eq(category)].set_index('GALID').reindex(ids[stage])
                use = c.eligible_usable.fillna(False).astype(bool)
                supported = c.category_support.fillna(False).astype(bool)
                f = c.f_category.where(use).to_numpy(float)
                vals['F_'+category] = (avg(f), mean_boot(f,stage), int(use.sum()))
                if category == 'NSF':
                    i = c.lha_mean_within_category.where(supported).to_numpy(float)
                    j = (c.sum_lha_within_category / c.N_usable).where(use).to_numpy(float)
                    vals['I_NSF'] = (avg(i), mean_boot(i,stage), int(supported.sum()))
                    vals['J_NSF'] = (avg(j), mean_boot(j,stage), int(use.sum()))
                    bp = c.bpt_common_support.fillna(False).astype(bool)
                    for line in ['HA6562','HB4861','OIII5006','NII6583','SII_SUM']:
                        lv = c['mean_'+line+'_within_category'].where(bp).to_numpy(float)
                        vals[line] = (avg(lv), mean_boot(lv,stage), int(bp.sum()))
            for metric,(value,boot,n) in vals.items():
                lo,hi = interval(boot)
                rows.append({'stage':stage,'variable':variable,'center':center,'metric':metric,
                             'value':value,'lo':lo,'hi':hi,'n_support':n,
                             'zero_draw_fraction':float(np.mean(boot==0)),
                             'nonfinite_draw_fraction':float(np.mean(~np.isfinite(boot)))})
                draws[(stage,variable,center,metric)] = boot
profile = pd.DataFrame(rows)
profile.to_csv(OUT/'fresh_summary_bootstrap.csv',index=False)

contrasts=[]
for stage in STAGES[1:]:
    target=profile[profile.stage.eq(stage)]
    for r in target.itertuples():
        ref=profile[profile.stage.eq('pre-peak') & profile.variable.eq(r.variable) & profile.center.eq(r.center) & profile.metric.eq(r.metric)].iloc[0]
        d=log(draws[(stage,r.variable,r.center,r.metric)])-log(draws[('pre-peak',r.variable,r.center,r.metric)])
        lo,hi=interval(d)
        contrasts.append({'stage':stage,'variable':r.variable,'center':r.center,'metric':r.metric,
                          'delta':float(log(r.value)-log(ref.value)), 'lo':lo,'hi':hi,
                          'n_target':r.n_support,'n_pre':ref.n_support,
                          'finite_draw_fraction':float(np.isfinite(d).mean())})
contrast=pd.DataFrame(contrasts)
contrast.to_csv(OUT/'fresh_two_stage_contrasts.csv',index=False)

def save(fig,name):
    fig.savefig(OUT/(name+'.png'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)

labels={'log_sigma_star':r'$\log_{10}(\Sigma_\star/[M_\odot\,{\rm kpc}^{-2}])$', 'radius_re':r'$R/R_e$'}
fig,axs=plt.subplots(2,3,figsize=(12,7),layout='constrained')
for k,var in enumerate(labels):
    for j,cat in enumerate(['SF','NSF','ND']):
        ax=axs[k,j]
        for s in STAGES:
            g=profile[profile.stage.eq(s)&profile.variable.eq(var)&profile.metric.eq('F_'+cat)].sort_values('center')
            ax.plot(g.center,g.value,c=COLORS[s],label=s,lw=1.7)
            poor = g.n_support < 3
            ax.scatter(g.loc[poor,'center'],g.loc[poor,'value'],facecolors='white',edgecolors=COLORS[s],zorder=4,s=38)
            ax.fill_between(g.center,g.lo,g.hi,color=COLORS[s],alpha=.14)
        ax.set(xlabel=labels[var],ylabel='Usable-area fraction',ylim=(-.03,1.03),title=cat)
        ax.grid(alpha=.15)
axs[0,0].legend(fontsize=8)
fig.suptitle('Live MAUVE fractions: strict post-fit S/N > 25\nEqual-galaxy means; 10,000 whole-system bootstrap draws',fontsize=12)
save(fig,'figure_01_live_occupancy')

fig,axs=plt.subplots(2,3,figsize=(12,7),layout='constrained')
specs=[('O','Selected SF occupancy'),('I','Surviving-SF intensity'),('T','HII-traced SFR / usable area'),('F_NSF','NSF occupancy'),('I_NSF','Halpha within NSF'),('J_NSF','NSF Halpha / usable area')]
for ax,(metric,title) in zip(axs.ravel(),specs):
    for s in STAGES[1:]:
        g=contrast[contrast.stage.eq(s)&contrast.variable.eq('log_sigma_star')&contrast.metric.eq(metric)].sort_values('center')
        ax.plot(g.center,g.delta,c=COLORS[s],label=s,lw=1.7)
        ax.fill_between(g.center,g.lo,g.hi,color=COLORS[s],alpha=.14)
        poor=(g.n_target<3)|(g.n_pre<3)
        ax.scatter(g.loc[poor,'center'],g.loc[poor,'delta'],facecolors='white',edgecolors=COLORS[s],zorder=4,s=38)
    ax.axhline(0,color='.4',lw=.8,ls='--');ax.grid(alpha=.15)
    ax.set(title=title,xlabel=labels['log_sigma_star'],ylabel='Stage minus pre-peak (dex)',ylim=(-2.8,.5))
axs[0,0].legend(fontsize=8)
fig.suptitle('Both stages resampled; open points have fewer than 3 supporting systems\nUndefined logarithms are omitted, not replaced by zero',fontsize=11)
save(fig,'figure_02_live_decomposition')

ratios=[]
fig,axs=plt.subplots(2,3,figsize=(12,7),layout='constrained')
for k,var in enumerate(labels):
    for j,(num,den,title) in enumerate([('OIII5006','HB4861','[O III] / Hbeta'),('NII6583','HA6562','[N II] / Halpha'),('SII_SUM','HA6562','[S II] / Halpha')]):
        ax=axs[k,j]
        for s in STAGES[1:]:
            g=profile[profile.stage.eq(s)&profile.variable.eq(var)&profile.metric.eq(num)].sort_values('center')
            values=[];los=[];his=[];poor=[]
            for r in g.itertuples():
                def val(st,m):
                    return profile[profile.stage.eq(st)&profile.variable.eq(var)&profile.center.eq(r.center)&profile.metric.eq(m)].value.iloc[0]
                v=float(log(val(s,num)/val(s,den))-log(val('pre-peak',num)/val('pre-peak',den)))
                d=log(draws[(s,var,r.center,num)])-log(draws[(s,var,r.center,den)])-log(draws[('pre-peak',var,r.center,num)])+log(draws[('pre-peak',var,r.center,den)])
                lo,hi=interval(d);values.append(v);los.append(lo);his.append(hi)
                np_=int(profile[profile.stage.eq('pre-peak')&profile.variable.eq(var)&profile.center.eq(r.center)&profile.metric.eq(num)].n_support.iloc[0])
                poor.append(min(r.n_support,np_)<3)
                ratios.append({'stage':s,'variable':var,'center':r.center,'ratio':title,'delta':v,'lo':lo,'hi':hi,'n_target':r.n_support,'n_pre':np_})
            ax.plot(g.center,values,c=COLORS[s],lw=1.7,label=s)
            ax.fill_between(g.center,los,his,color=COLORS[s],alpha=.14)
            ax.scatter(g.center.to_numpy()[poor],np.array(values)[poor],facecolors='white',edgecolors=COLORS[s],s=38,zorder=4)
        ax.axhline(0,c='.4',lw=.8,ls='--');ax.grid(alpha=.15)
        ax.set(xlabel=labels[var],ylabel='NSF ratio change (dex)',title=title,ylim=(-.8,1.15))
axs[0,0].legend(fontsize=8)
fig.suptitle('Common BPT-line support; direct Hbeta denominator\nRatios of equal-galaxy linear line means; both-stage bootstrap',fontsize=12)
save(fig,'figure_03_live_line_ratios')
pd.DataFrame(ratios).to_csv(OUT/'fresh_true_bpt_contrasts.csv',index=False)

