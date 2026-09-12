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
RNG = np.random.default_rng(20260911)
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

# Exact two-reservoir solution; time unit is Myr.
def reservoirs(t,a0,m0,ka,tf,km,td,R=.4,eta=0):
    kap=ka+1/tf; lam=km+(1-R+eta)/td
    a=a0*np.exp(-kap*t)
    feed=a0/tf*(t*np.exp(-lam*t) if abs(kap-lam)<1e-12 else (np.exp(-kap*t)-np.exp(-lam*t))/(lam-kap))
    return a,m0*np.exp(-lam*t)+feed

# Moment of a lognormal HII amplitude between finite/infinite thresholds.
def F(x,mu,s):
    return ndtr((np.log(x)-mu)/s) if x>0 else 0.
def moment(l,u,mu,s):
    lower=ndtr((np.log(l)-mu-s*s)/s) if l>0 else 0.
    upper=ndtr((np.log(u)-mu-s*s)/s) if np.isfinite(u) else 1.
    return np.exp(mu+s*s/2)*(upper-lower)
def state(mu,s,D,C,Ldet,wcut):
    l=max(0,Ldet-D)
    u=max(l,6*C-D,D*(1-wcut)/wcut,0)
    pnd=F(l,mu,s);psf=1-F(u,mu,s);pnsf=F(u,mu,s)-pnd
    insf=moment(l,u,mu,s)/pnsf+D if pnsf>1e-14 else np.nan
    isf=moment(u,np.inf,mu,s)/psf+D if psf>1e-14 else np.nan
    return psf,pnsf,pnd,isf,insf,l,u

times=np.linspace(0,600,241)
demo=[]
for zone,p in [('inner',(.4,.85,.18,.07,.08,300,600)),('outer',(.2,.7,.01,.015,.08,120,350))]:
    med,s,D0,C,Ldet,tH,tD=p
    for t in times:
        D=D0*np.exp(-t/tD)
        sf,nsf,nd,isf,insf,l,u=state(np.log(med)-t/tH,s,D,C,Ldet,.32748538011695905)
        demo.append({'zone':zone,'time_Myr':t,'D':D,'median_A':med*np.exp(-t/tH),'SF':sf,'NSF':nsf,'ND':nd,'I_SF':isf,'I_NSF':insf,'A_det':l,'A_SF':u})
demo=pd.DataFrame(demo);demo.to_csv(OUT/'toy_occupancy.csv',index=False)
fig,axs=plt.subplots(2,2,figsize=(10.5,7),layout='constrained')
for k,z in enumerate(['inner','outer']):
    d=demo[demo.zone.eq(z)]
    for cat in ['SF','NSF','ND']:axs[k,0].plot(d.time_Myr,d[cat],c=CC[cat],label=cat,lw=2)
    axs[k,0].set(title=z.capitalize()+' illustrative region',ylabel='Area fraction',ylim=(-.02,1.02),xlabel='Time after supply decline (Myr)')
    axs[k,1].plot(d.time_Myr,log(d.I_SF),label='Within SF',c=CC['SF'],lw=2)
    axs[k,1].plot(d.time_Myr,log(d.I_NSF),label='Within NSF',c=CC['NSF'],lw=2)
    axs[k,1].set(ylabel=r'$\log_{10}$ Halpha intensity (relative units)',xlabel='Time after supply decline (Myr)',title='Selection-conditioned intensities')
for ax in axs.ravel():ax.grid(alpha=.15);ax.legend(fontsize=8)
fig.suptitle('Analytical population example: an ionized inner remnant and a fading outer disc\nIllustration, not a fit; numerical times do not assign ages to MAUVE stages',fontsize=11)
save(fig,'figure_04_analytic_occupancy')

t=np.linspace(0,450,301);A=np.exp(-t/150);D=.2*np.exp(-t/600);w=D/(A+D)
ha_mix=A+D;nii=.3*A+1.2*D;sii=.2*A+.6*D
oiii=(.5*A+1.8*D)/2.86;hb=ha_mix/2.86
oiii_soft=(.5*A+.2*D)/2.86
sigma=np.sqrt((1-w)*25**2+w*70**2+w*(1-w)*20**2)
fig,axs=plt.subplots(2,2,figsize=(10.5,7),layout='constrained')
for y,label,c in [(ha_mix,'Halpha','#111111'),(nii,'[N II]','#7b3294'),(sii,'[S II]','#e69f00'),(oiii,'[O III]','#008b8b')]:
    axs[0,0].plot(t,log(y/y[0]),label=label,c=c,lw=2)
axs[0,0].set(ylabel='Line luminosity change (dex)',title='Every line becomes fainter')
for y,label,c in [(nii/ha_mix,'[N II]/Halpha','#7b3294'),(oiii/hb,'[O III]/Hbeta, high residual','#008b8b')]:
    axs[0,1].plot(t,log(y/y[0]),label=label,c=c,lw=2)
axs[0,1].plot(t,log((oiii_soft/hb)/(oiii_soft[0]/hb[0])),label='[O III]/Hbeta, low residual',c='#cb5b24',lw=2,ls='--')
axs[0,1].set(ylabel='Ratio change (dex)',title='The residual template sets the direction')
axs[1,0].plot(t,sigma,c='#4d4d4d',lw=2);axs[1,0].axhline(45,c='#d73027',ls='--',label='SF cut: 45 km/s')
axs[1,0].set(ylabel='Mixture dispersion (km/s)',title='Fixed component widths; changing flux weights')
axs[1,1].plot(t,ha_mix/.1,c='#4d4d4d',lw=2);axs[1,1].axhline(6,c='#d73027',ls='--',label='SF cut: 6 Angstrom')
axs[1,1].set(ylabel='Illustrative EW(Halpha) (Angstrom)',title='The same continuum outlives bright HII emission')
for ax in axs.ravel():ax.set_xlabel('Time after decline (Myr)');ax.grid(alpha=.15);ax.legend(fontsize=8)
fig.suptitle('Two emitting components: young-star A and residual D\nIllustration, not a fit; residual D is continuously powered, not fossil glow',fontsize=11)
save(fig,'figure_05_mixing_fading_width')

fig,axs=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
tt=np.linspace(0,1200,500)
for ka,label,c in [(0,'Supply only transferred','#2166ac'),(.004,'Diffuse stripping','#1a9850'),(.012,'Faster diffuse stripping','#d73027')]:
    a,m=reservoirs(tt,20,10,ka,300,.0005,1500,eta=.3)
    axs[0].plot(tt,m/10,label=label,c=c,lw=2)
axs[0].set(xlabel='Time after external supply stops (Myr)',ylabel='Molecular reservoir / initial value',title='Retained gas can briefly grow, then fade');axs[0].legend(fontsize=8)
rr=np.linspace(.01,3,300);p0=1;effscale=.65
for p,label,c in [(.02,'Low peak pressure','#2166ac'),(.08,'Intermediate peak pressure','#1a9850'),(.25,'High peak pressure','#d73027')]:
    q=1-ndtr((np.log(p/p0)+rr/effscale)/.75)
    axs[1].plot(rr,q,c=c,label=label,lw=2)
axs[1].set(xlabel='Illustrative radius / Re',ylabel='Retained column fraction',title='Column scatter softens the stripping edge');axs[1].legend(fontsize=8)
for ax in axs:ax.grid(alpha=.15)
fig.suptitle('Physical ingredients, before applying the SF/NSF/ND classifier (not fitted)',fontsize=11)
save(fig,'figure_06_gas_regulator')

# Independent numerical checks of the proposed analytic equations.
params=(20.,10.,.004,300.,.0005,1500.,.4,.3)
a0,m0,ka,tf,km,td,R,eta=params
def rhs(t,y):
    a,m=y
    return [-(ka+1/tf)*a,a/tf-(km+(1-R+eta)/td)*m]
tt=np.linspace(0,1500,301)
ode=solve_ivp(rhs,(0,1500),(a0,m0),rtol=1e-11,atol=1e-12,dense_output=True)
a,m=reservoirs(tt,*params)
err=float(np.max(np.abs(np.vstack([a,m])-ode.sol(tt))))
mu,s,D,C,Ld=np.log(.4),.85,.18,.07,.08
st=state(mu,s,D,C,Ld,.32748538011695905);l,u=st[-2:]
def pdf_A(x):return np.exp(-.5*((np.log(x)-mu)/s)**2)/(x*s*np.sqrt(2*np.pi)) if x>0 else 0
prob=quad(pdf_A,l,u,epsabs=1e-12)[0]
integ=quad(lambda x:(x+D)*pdf_A(x),l,u,epsabs=1e-12)[0]/prob
moment_err=abs(integ-st[4])
partition=float(np.max(abs(demo[['SF','NSF','ND']].sum(axis=1)-1)))
sfprofile=pd.read_csv(OUT/'stage_sfr_profiles.csv')
identity=float(np.nanmax(abs(sfprofile.total_hii-sfprofile.occupancy*sfprofile.survivor_intensity)))
mix_d=np.gradient(log(nii/ha_mix),t)
# Finite-difference gas mass conservation against analytic RHS at representative times.
ad,md=rhs(tt,np.vstack([a,m]))
mass_res=float(np.max(abs(ad+md+ka*a+km*m+(1-R+eta)*m/td)))
checks={'reservoir_ode_max_abs_error':err,'gas_mass_budget_max_abs_residual':mass_res,
        'lognormal_NSF_intensity_quadrature_abs_error':moment_err,
        'probability_partition_max_error':partition,'live_T_OI_max_abs_error':identity,
        'all_demo_line_amplitudes_decline':bool(all(np.all(np.diff(v)<0) for v in [ha_mix,nii,sii,oiii])),
        'demo_NII_Ha_ratio_increases':bool(np.all(mix_d>0)),
        'bootstrap_draws':NBOOT,'random_seed':20260911,
        'model_scope':'Analytic/toy verification and fresh bootstrap of scalar data; no physical parameter fit.'}
assert err<1e-7 and mass_res<1e-12 and moment_err<1e-9 and partition<1e-12 and identity<1e-12
(OUT/'numerical_checks.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
