"""Constant-coefficient HI-only gas model and photon-budget-limited HOLMES model.

Illustrative calculations, not a fit or a new extraction of MAUVE maps.
Every observational anchor is rebuilt from the existing five-line export.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parent / '20260914_resolved_RPS_academic_model'
records = json.loads((SOURCE / 'input_fingerprints.json').read_text())
for record in records:
    assert hashlib.sha256(Path(record['path']).read_bytes()).hexdigest() == record['sha256'], record['path']
(OUT / 'source_fingerprints.json').write_text(json.dumps(records, indent=2)+'\n')
src = SOURCE / 'stage_bpt_line_profiles_with_hbeta.csv'
data = pd.read_csv(src)
data = data[(data.variable == 'log_sigma_star') & (data.center == 8.625)]
rows = []
for (stage, category), group in data.groupby(['stage','category']):
    lines = group.set_index('line')
    lum = lines.mean_line_surface.to_dict()
    rows.append(dict(stage=stage, category=category, n_gal_common=int(lines.N_gal_support.min()),
        ha_surface=lum['HA6562'], hb_surface=lum['HB4861'],
        nii_ha=lum['NII6583']/lum['HA6562'], sii_ha=lum['SII_SUM']/lum['HA6562'],
        oiii_hb=lum['OIII5006']/lum['HB4861']))
anchors = pd.DataFrame(rows)
anchors.to_csv(OUT/'mauve_anchor_table.csv', index=False)
idx = anchors.set_index(['stage','category'])
pre, post = idx.loc[('pre-peak','NSF')], idx.loc[('post-peak','NSF')]

C_ALPHA = 4.983582089552239e-42
H_PLANCK, C_LIGHT, WAVELENGTH_ALPHA = 6.62607015e-27, 2.99792458e10, 6562.8e-8
P_ALPHA, Q_HOLMES, SIGMA_OLD = 1/2.206, 7e40, 10**8.625
EPS_ALPHA = H_PLANCK*C_LIGHT/WAVELENGTH_ALPHA*P_ALPHA
L_HOLMES = EPS_ALPHA*Q_HOLMES*SIGMA_OLD  # upper benchmark: f_abs=1, all stars old
RETURN, LOADING, TAU_DEP, TAU_ION = .4, 0., 2., .003
GAMMA_H2 = (1-RETURN+LOADING)/TAU_DEP
SFR0 = C_ALPHA*(pre.ha_surface-L_HOLMES)  # total absorbed young component
H20, HI0 = 1000*TAU_DEP*SFR0, 10.
TAU_CONV = HI0/(GAMMA_H2*H20)
GAMMA_STRIP = 3.

def normalized_sfr(t, gamma_hi, gamma_h2, tau_phi0):
    t=np.asarray(t)
    if abs(gamma_hi-gamma_h2)<1e-9:
        return (1+t/tau_phi0)*np.exp(-gamma_h2*t)
    return (np.exp(-gamma_h2*t)+(np.exp(-gamma_hi*t)-np.exp(-gamma_h2*t))
            /(tau_phi0*(gamma_h2-gamma_hi)))

def filtered_exponential(t, gamma, tau=TAU_ION):
    t=np.asarray(t)
    if abs(1-gamma*tau)<1e-9:
        return (1+t/tau)*np.exp(-t/tau)
    return (np.exp(-gamma*t)-gamma*tau*np.exp(-t/tau))/(1-gamma*tau)

def filtered_sfr(t, gamma_hi, gamma_h2, tau_phi0):
    # Equal gas-rate branch is checked by independent numerical convolution below.
    assert abs(gamma_hi-gamma_h2)>1e-9
    coefficient=1/(tau_phi0*(gamma_h2-gamma_hi))
    return ((1-coefficient)*filtered_exponential(t,gamma_h2)
            +coefficient*filtered_exponential(t,gamma_hi))

def peak_time(gamma_hi, gamma_h2, tau_phi0):
    if 1/tau_phi0<=gamma_h2: return 0.
    if abs(gamma_hi-gamma_h2)<1e-9: return 1/gamma_h2-tau_phi0
    return np.log(gamma_hi/(gamma_h2*(1+(gamma_hi-gamma_h2)*tau_phi0)))/(gamma_hi-gamma_h2)

t = np.linspace(0,2,2001)
gas_rows, checks, cases = [], {}, {}
for label, compression, strip in [('balanced',1.,3.),('closed_no_RPS',1.,0.),('compression',.25,3.)]:
    tauconv=TAU_CONV*compression
    gammahi=1/tauconv+strip
    tauphi0=H20*tauconv/HI0
    frac=normalized_sfr(t,gammahi,GAMMA_H2,tauphi0)
    filtered=filtered_sfr(t,gammahi,GAMMA_H2,tauphi0)
    atomic=HI0*np.exp(-gammahi*t)
    molecular=H20*frac
    numerical=solve_ivp(lambda tt,y:[-gammahi*y[0],y[0]/tauconv-GAMMA_H2*y[1],
        (y[1]/H20-y[2])/TAU_ION],(0,2),[HI0,H20,1.],t_eval=t,
        rtol=1e-10,atol=1e-12,max_step=.001)
    assert numerical.success
    gas_error=float(np.max(np.abs(numerical.y[:2]-np.array([atomic,molecular]))/np.maximum(np.array([atomic,molecular]),1e-8)))
    kernel_error=float(np.max(np.abs(numerical.y[2]-filtered)))
    assert gas_error<1e-7 and kernel_error<1e-8
    assert np.all(frac>=np.exp(-GAMMA_H2*t)-1e-12)
    assert np.max(molecular)<=H20+HI0+1e-10
    peak=peak_time(gammahi,GAMMA_H2,tauphi0)
    peak_fraction=float(normalized_sfr(peak,gammahi,GAMMA_H2,tauphi0))
    if peak:
        residual=(HI0*np.exp(-gammahi*peak)/tauconv-GAMMA_H2*H20*peak_fraction)
        assert abs(residual)<1e-10
    cases[label]=dict(tau_conv=tauconv,gamma_hi=gammahi,tau_phi0=tauphi0,
        sfr_1gyr=float(frac[1000]),young_ha_1gyr=float(filtered[1000]),
        t_peak_gyr=peak,peak_sfr_fraction=peak_fraction,peak_dex=float(np.log10(peak_fraction)))
    checks[label]=dict(gas_ode_max_rel=gas_error,filtered_ode_max_abs=kernel_error)
    gas_rows.append(pd.DataFrame(dict(case=label,t_gyr=t,hi=atomic,h2=molecular,
        sfr=SFR0*frac,sfr_fraction=frac,young_ha_fraction=filtered,
        total_ha=(pre.ha_surface-L_HOLMES)*filtered+L_HOLMES,
        holmes_ha_weight=L_HOLMES/((pre.ha_surface-L_HOLMES)*filtered+L_HOLMES))))
gas=pd.concat(gas_rows,ignore_index=True)
gas.to_csv(OUT/'gas_and_connected_Halpha.csv',index=False)

# Non-balanced equal-rate branch, checked independently, including its peak.
deg_gamma, deg_tauphi = .3, 1.
deg=solve_ivp(lambda tt,y:[-deg_gamma*y[0],y[0]-deg_gamma*y[1]],
    (0,2),[1/deg_tauphi,1.],t_eval=t,rtol=1e-11,atol=1e-12,max_step=.002)
deg_error=float(np.max(np.abs(deg.y[1]-normalized_sfr(t,deg_gamma,deg_gamma,deg_tauphi))))
assert deg_error<1e-9
checks['equal_rate_max_abs']=deg_error

# Photon-budget-constrained initial normalization, never fitting endpoint spectra.
endpoints={'nii_ha':1.5,'sii_ha':1.,'oiii_hb':3.}
comparison=[]
for line,rholmes in endpoints.items():
    # Preserve the exact exported-denominator audit for Appendix D.
    denom='hb_surface' if line=='oiii_hb' else 'ha_surface'
    old=L_HOLMES/2.86 if line=='oiii_hb' else L_HOLMES
    w0,w1=old/pre[denom],old/post[denom]
    ryoung=(pre[line]-w0*rholmes)/(1-w0)
    predicted=(1-w1)*ryoung+w1*rholmes
    slope=(post[line]-pre[line])/(1/post[denom]-1/pre[denom])
    young_required=pre[line]-slope/pre[denom]
    holmes_required=young_required+slope/old
    old_needed=slope/(rholmes-young_required)
    comparison.append(dict(line=line,r_holmes_assumed=rholmes,r_young_initial=ryoung,
        w_pre=w0,w_post=w1,ratio_pre=pre[line],ratio_post=post[line],
        ratio_post_predicted=predicted,residual_dex=np.log10(predicted/post[line]),
        r_young_exact_two_point=young_required,r_holmes_required_at_fiducial=holmes_required,
        old_balmer_needed=old_needed,photon_budget_multiple=old_needed/old))
comparison=pd.DataFrame(comparison)
comparison.to_csv(OUT/'line_budget_actual_denominator.csv',index=False)
actual_comparison=comparison.copy()
# Model 0 uses the common Halpha weight for both Balmer denominators.
comparison=[]
for line,rholmes in endpoints.items():
    w0,w1=L_HOLMES/pre.ha_surface,L_HOLMES/post.ha_surface
    ryoung=(pre[line]-w0*rholmes)/(1-w0)
    predicted=(1-w1)*ryoung+w1*rholmes
    contrast=(post[line]-pre[line])/(w1-w0)
    young_required=pre[line]-w0*contrast
    comparison.append(dict(line=line,r_holmes_assumed=rholmes,r_young_initial=ryoung,
        w_pre=w0,w_post=w1,ratio_pre=pre[line],ratio_post=post[line],
        ratio_post_predicted=predicted,residual_dex=np.log10(predicted/post[line]),
        r_young_exact_two_point=young_required,r_holmes_required_at_fiducial=young_required+contrast,
        old_balmer_needed=(post[line]-pre[line])/((rholmes-young_required)*(1/post.ha_surface-1/pre.ha_surface)),
        photon_budget_multiple=contrast/(rholmes-young_required)))
comparison=pd.DataFrame(comparison)
comparison.to_csv(OUT/'line_budget_constraints.csv',index=False)

# Main-text effective young + HOLMES examples, without an unnecessary split.
fyoung=np.geomspace(1,.01,501)
line_rows=[]
for name,initial in [('bright_NSF_anchor',pre.ha_surface),('faint_example',2e38)]:
    young=(initial-L_HOLMES)*fyoung
    old=np.full_like(young,L_HOLMES)
    total=young+old
    for line,ends in {'nii_ha':(.35,1.5),'sii_ha':(.30,1.),'oiii_hb':(.30,3.)}.items():
        lineflux=ends[0]*young+ends[1]*old
        ratio=lineflux/total # equal decrements for O3 demonstration
        assert np.all(np.diff(ratio)>0) and np.all(np.diff(lineflux)<0)
        line_rows.append(pd.DataFrame(dict(case=name,line=line,young_fraction=fyoung,
            total_ha=total,holmes_weight=old/total,ratio=ratio,
            line_surface=lineflux/(2.86 if line=='oiii_hb' else 1.))))
demonstration=pd.concat(line_rows,ignore_index=True)
demonstration.to_csv(OUT/'young_HOLMES_fading.csv',index=False)

# Luminosity addition and weights agree even for unequal decrements.
hal=np.array([7.,3.,1.]); decrements=np.array([2.86,3.1,3.4]); ratios=np.array([.2,.5,3.])
hb=hal/decrements
direct=np.sum(ratios*hb)/np.sum(hb)
weighted=np.sum((hal/hal.sum()/decrements)/np.sum(hal/hal.sum()/decrements)*ratios)
checks['unequal_decrement_identity_abs']=float(abs(direct-weighted))
assert abs(direct-weighted)<1e-14

required_young_fraction=(post.ha_surface-L_HOLMES)/(pre.ha_surface-L_HOLMES)
balanced=cases['balanced']
target_time=brentq(lambda time: filtered_sfr(time,balanced['gamma_hi'],GAMMA_H2,balanced['tau_phi0'])-required_young_fraction,0,20)
instant_target=brentq(lambda time: normalized_sfr(time,balanced['gamma_hi'],GAMMA_H2,balanced['tau_phi0'])-required_young_fraction,0,20)
# Spatially elevated initial state at unchanged coefficients, not faster conversion.
spatial_factor=1.5
base=gas[gas.case=='balanced'].reset_index(drop=True)
control=gas[gas.case=='closed_no_RPS'].reset_index(drop=True)
spatial= pd.DataFrame(dict(t_gyr=t,reference_sfr=control.sfr,
    accumulated_sfr=spatial_factor*base.sfr,
    spatial_excess_dex=np.log10(spatial_factor*base.sfr/control.sfr)))
spatial.to_csv(OUT/'spatial_excess_at_fixed_efficiency.csv',index=False)
spatial_ode=solve_ivp(lambda tt,y:[-balanced['gamma_hi']*y[0],y[0]/TAU_CONV-GAMMA_H2*y[1]],
    (0,2),[spatial_factor*HI0,spatial_factor*H20],t_eval=t,rtol=1e-10,atol=1e-12,max_step=.002)
checks['spatial_initial_state_max_abs']=float(np.max(np.abs(spatial_ode.y[1]/H20-spatial_factor*base.sfr_fraction)))
assert checks['spatial_initial_state_max_abs']<1e-9
assert np.all(np.diff(spatial.accumulated_sfr)<0) and np.all(spatial.spatial_excess_dex>0)
# Actual line fading factors from unmodified observed luminosity denominators.
fading={'ha':float(post.ha_surface/pre.ha_surface),'hb':float(post.hb_surface/pre.hb_surface)}
for line,denom in [('nii_ha','ha_surface'),('sii_ha','ha_surface'),('oiii_hb','hb_surface')]:
    fading[line]=float(post[line]*post[denom]/(pre[line]*pre[denom]))
summary=dict(epsilon_alpha=EPS_ALPHA,p_alpha=P_ALPHA,q_holmes=Q_HOLMES,sigma_old=SIGMA_OLD,
    holmes_ha_ceiling_fiducial=L_HOLMES,holmes_ha_FSPS=L_HOLMES*5/7,
    sfr0=SFR0,h2_0=H20,hi_0=HI0,tau_dep=TAU_DEP,gamma_h2=GAMMA_H2,
    tau_conv=TAU_CONV,cases=cases,mass_budget_peak_bound=1+HI0/H20,
    mass_budget_peak_dex=float(np.log10(1+HI0/H20)),
    required_young_fraction=required_young_fraction,
    minimum_decline_time_gyr=float(-np.log(required_young_fraction)/GAMMA_H2),
    conditional_balanced_target_time_gyr=target_time,
    instantaneous_balanced_target_time_gyr=instant_target,
    response_relative_correction_at_1gyr=float(balanced['young_ha_1gyr']/balanced['sfr_1gyr']-1),
    spatial_factor=spatial_factor,spatial_excess_dex_at_1gyr=float(spatial.spatial_excess_dex.iloc[1000]),
    observed_line_fading=fading,actual_denominator_budget=actual_comparison.to_dict('records'),
    equivalent_tau_dep_upper_at_1gyr=float((1-RETURN)/(-np.log(required_young_fraction))),
    line_budget=comparison.to_dict('records'),checks=checks,
    input_csv_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
    classification='Illustrative model and conditional mean-value consistency checks; not a fit')
(OUT/'numerical_audit.json').write_text(json.dumps(summary,indent=2)+'\n')

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':160})
colors={'balanced':'#24628b','closed_no_RPS':'#757575','compression':'#b75522'}
fig,axs=plt.subplots(1,2,figsize=(9,3.6),constrained_layout=True)
for name,grp in gas[gas.case!='compression'].groupby('case',sort=False):
    label={'balanced':'Initial balance + stripping','closed_no_RPS':'Closed, no stripping','compression':'Faster conversion + stripping'}[name]
    axs[0].plot(grp.t_gyr,grp.sfr_fraction,color=colors[name],label=label)
    axs[1].plot(grp.t_gyr,grp.hi/HI0,color=colors[name])
axs[0].plot(t,spatial_factor*base.sfr_fraction,color='#b75522',label='1.5 x initial gas + stripping')
axs[0].axhline(1,color='k',lw=.6,ls=':');axs[0].legend(fontsize=8)
axs[0].set(xlabel='Time since model onset [Gyr]',ylabel=r'$\Sigma_{\rm SFR}(t)/\Sigma_{\rm SFR,reference}(0)$')
axs[1].set(xlabel='Time since model onset [Gyr]',ylabel=r'$\Sigma_{\rm HI}(t)/\Sigma_{\rm HI}(0)$',yscale='log',ylim=(1e-7,1.2))
for ext in ['png','pdf']:fig.savefig(OUT/f'figure01_SFR_response.{ext}')
plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(10,3.3),constrained_layout=True)
for name,style in [('bright_NSF_anchor','-'),('faint_example','--')]:
    sub=demonstration[(demonstration.case==name)&(demonstration.line=='nii_ha')]
    label='Bright NSF scale' if name.startswith('bright') else 'Faint illustration'
    axs[0].plot(sub.young_fraction,sub.holmes_weight,style,label=label)
    axs[1].plot(sub.total_ha,sub.ratio,style,label=label)
    axs[2].plot(sub.young_fraction,sub.line_surface/sub.line_surface.iloc[0],style,label=label)
axs[0].set(xscale='log',xlabel='Remaining young Hα fraction',ylabel='HOLMES Hα weight',xlim=(1,.01));axs[0].legend(fontsize=8)
axs[1].set(xscale='log',xlabel='Total Hα [erg s$^{-1}$ kpc$^{-2}$]',ylabel='[N II]/Hα',xlim=(4e39,3e37))
axs[2].set(xscale='log',xlabel='Remaining young Hα fraction',ylabel='Remaining [N II] fraction',xlim=(1,.01))
for ext in ['png','pdf']:fig.savefig(OUT/f'figure02_fading_and_ratio.{ext}')
plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(9.5,3.3),constrained_layout=True)
for ax,(_,row) in zip(axs,comparison.iterrows()):
    ax.plot([0,1],[row.ratio_pre,row.ratio_post_predicted],'o-',label='Constant HOLMES model')
    ax.plot([0,1],[row.ratio_pre,row.ratio_post],'s--',color='#b75522',label='Exported NSF means')
    ax.set(xticks=[0,1],xticklabels=['Pre-peak','Post-peak'],ylabel={'nii_ha':'[N II]/Hα','sii_ha':'[S II]/Hα','oiii_hb':'[O III]/Hβ'}[row.line])
axs[0].legend(fontsize=8)
for ext in ['png','pdf']:fig.savefig(OUT/f'figure03_MAUVE_budget_test.{ext}')
plt.close(fig)
print(json.dumps(summary,indent=2))
