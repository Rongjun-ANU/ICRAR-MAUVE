"""HI-only analytical predictions and conditional HII/non-HII mixing audit.

No notebooks/FITS are modified or executed. Scalar anchors are regenerated from
the 14 September export only after verifying its recorded notebook/pipeline
fingerprints. The large source FITS maps are not individually fingerprinted.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
OLD = OUT.parent / '20260914_resolved_RPS_academic_model'
records = json.loads((OLD / 'input_fingerprints.json').read_text())
for record in records:
    path = Path(record['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
        raise RuntimeError(f'Source changed since scalar extraction: {path}')
(OUT / 'source_fingerprints.json').write_text(json.dumps(records, indent=2)+'\n')
source = OLD / 'stage_bpt_line_profiles_with_hbeta.csv'
data = pd.read_csv(source)
data = data.loc[(data.variable == 'log_sigma_star') & (data.center == 8.625)]
anchors = []
for (stage, category), group in data.groupby(['stage', 'category']):
    lines = group.set_index('line')
    luminosities = lines.mean_line_surface.to_dict()
    anchors.append(dict(stage=stage, category=category,
        n_gal_common=int(lines.N_gal_support.min()),
        ha_surface=luminosities['HA6562'], hb_surface=luminosities['HB4861'],
        nii_surface=luminosities['NII6583'], sii_surface=luminosities['SII_SUM'],
        oiii_surface=luminosities['OIII5006'],
        nii_ha=luminosities['NII6583']/luminosities['HA6562'],
        sii_ha=luminosities['SII_SUM']/luminosities['HA6562'],
        oiii_hb=luminosities['OIII5006']/luminosities['HB4861'],
        corrected_balmer=luminosities['HA6562']/luminosities['HB4861']))
adf = pd.DataFrame(anchors)
adf.to_csv(OUT/'mauve_anchor_table.csv',index=False)
indexed = adf.set_index(['stage','category'])

C_ALPHA = 4.983582089552239e-42
TAU_DEP = 2.0
RETURN = .4
ETA = 0.0
HI0 = 10.0
GAMMA_STRIP_HI = 3.0
F_CAPTURE = .95
W_NON_NULL = .10
TAU_ION = .003
B_HII = B_NON = 2.86
L_PRE_SF = float(indexed.loc[('pre-peak','SF'),'ha_surface'])
SFR0 = C_ALPHA * L_PRE_SF/F_CAPTURE
H20 = 1000.0*TAU_DEP*SFR0
GAMMA_CONS = (1-RETURN+ETA)/TAU_DEP
TAU_CONV = HI0/(GAMMA_CONS*H20)
GAMMA_HI = 1/TAU_CONV+GAMMA_STRIP_HI

def normalized_sfr(t, gamma_hi=GAMMA_HI, gamma_cons=GAMMA_CONS):
    t=np.asarray(t)
    if abs(gamma_hi-gamma_cons)<1e-10:
        return (1+gamma_cons*t)*np.exp(-gamma_cons*t)
    return (gamma_hi*np.exp(-gamma_cons*t)
            -gamma_cons*np.exp(-gamma_hi*t))/(gamma_hi-gamma_cons)

def reservoirs(t, gamma_hi=GAMMA_HI):
    return HI0*np.exp(-gamma_hi*np.asarray(t)), H20*normalized_sfr(t,gamma_hi)

t=np.linspace(0,3,1501)
hi,h2=reservoirs(t)
f_sfr=normalized_sfr(t)
f_control=normalized_sfr(t,1/TAU_CONV)
sol=solve_ivp(lambda tt,y:[-GAMMA_HI*y[0],y[0]/TAU_CONV-GAMMA_CONS*y[1]],
    (0,3),[HI0,H20],t_eval=t,rtol=1e-11,atol=1e-12,max_step=.002)
assert sol.success,sol.message
ode_error=float(np.max(np.abs(sol.y-np.vstack([hi,h2]))/np.maximum(np.vstack([hi,h2]),1e-12)))
# Balance has a zero first derivative; second derivative is exact analytically.
initial_sfr_derivative=(HI0/TAU_CONV-GAMMA_CONS*H20)/(1000*TAU_DEP)
initial_sfr_second_derivative=-GAMMA_HI*GAMMA_CONS*SFR0
# Stable degenerate branch checked independently with matched reservoir balance.
deg_cons=.3
deg_hi0=deg_cons*H20*TAU_CONV
deg=solve_ivp(lambda tt,y:[-deg_cons*y[0],y[0]/TAU_CONV-deg_cons*y[1]],
    (0,3),[deg_hi0,H20],t_eval=t,rtol=1e-11,atol=1e-12,max_step=.002)
deg_exact=H20*normalized_sfr(t,deg_cons,deg_cons)
deg_error=float(np.max(np.abs(deg.y[1]-deg_exact)/deg_exact))

ion=solve_ivp(lambda tt,y:[(SFR0*normalized_sfr(tt)-y[0])/TAU_ION],
    (0,3),[SFR0],t_eval=t,rtol=1e-10,atol=1e-12,max_step=.001)
assert ion.success,ion.message
l_full=ion.y[0]/C_ALPHA
l_ha=F_CAPTURE*l_full
l_hii=(1-W_NON_NULL)*l_ha
l_non=W_NON_NULL*l_ha
templates={'nii_ha':(.22,.60),'sii_ha':(.20,.45),'oiii_hb':(.57,.77)}
ratios={key: np.full_like(t,(1-W_NON_NULL)*ends[0]+W_NON_NULL*ends[1])
        for key,ends in templates.items()}
null=pd.DataFrame(dict(time_gyr=t,hi=hi,h2=h2,sigma_sfr=SFR0*f_sfr,
    sfr_fraction=f_sfr,closed_no_rps_fraction=f_control,
    ha_surface=l_ha,hii_ha_surface=l_hii,non_hii_ha_surface=l_non,
    w_non_hii=np.full_like(t,W_NON_NULL),**ratios))
null.to_csv(OUT/'gas_and_fixed_partition_predictions.csv',index=False)

diagnostics=adf[['stage','category','n_gal_common']].copy()
for key,(rh,rd) in templates.items():
    diagnostics['w_'+key]=(adf[key]-rh)/(rd-rh)
    diagnostics['admissible_'+key]=diagnostics['w_'+key].between(0,1)
diagnostics.to_csv(OUT/'single_ratio_weight_diagnostics.csv',index=False)

# Conditional line-only interpolation: fixes observed NSF Halpha and NII/Ha
# at two stages. It does NOT fit a temporal trajectory or independent spectra.
pre=indexed.loc[('pre-peak','NSF')]
post=indexed.loc[('post-peak','NSF')]
rh,rd=templates['nii_ha']
w0=(pre.nii_ha-rh)/(rd-rh)
w1=(post.nii_ha-rh)/(rd-rh)
l0,l1=float(pre.ha_surface),float(post.ha_surface)
gamma_hii=-np.log((1-w1)*l1/((1-w0)*l0))
gamma_non=-np.log(w1*l1/(w0*l0))
tl=np.linspace(0,1,501)
lh=(1-w0)*l0*np.exp(-gamma_hii*tl)
ld=w0*l0*np.exp(-gamma_non*tl)
lt=lh+ld
wl=ld/lt
line_ratios={key:(1-wl)*ends[0]+wl*ends[1] for key,ends in templates.items()}
line_predictions=pd.DataFrame(dict(time_coordinate_gyr=tl,ha_surface=lt,
    hii_ha_surface=lh,non_hii_ha_surface=ld,w_non_hii=wl,**line_ratios))
line_predictions['nii_surface']=lt*line_predictions.nii_ha
line_predictions['sii_surface']=lt*line_predictions.sii_ha
line_predictions['oiii_surface']=lt/B_HII*line_predictions.oiii_hb
# A finite-interval photon-budget realization, conditional on assigning the
# NSF normalization the same normalized gas-shaped young-source history.
# This determines required allocation fractions; it does not predict them.
source_fraction=np.interp(tl,t,ion.y[0]/SFR0)
nsf_young_max=(l0/F_CAPTURE)*source_fraction
budget_hii=lh/nsf_young_max
budget_leak=ld/nsf_young_max
budget_unabs=1-budget_hii-budget_leak
line_predictions['required_f_hii']=budget_hii
line_predictions['required_f_leak_absorbed']=budget_leak
line_predictions['required_f_unabs']=budget_unabs
line_predictions.to_csv(OUT/'conditional_NSF_line_interpolation.csv',index=False)
endpoint_rows=[]
for i,stage in [(0,'pre-peak'),(-1,'post-peak')]:
    for key in templates:
        observed=float(indexed.loc[(stage,'NSF'),key])
        predicted=float(line_predictions.iloc[i][key])
        endpoint_rows.append(dict(stage=stage,ratio=key,observed=observed,
            predicted=predicted,fractional_residual=predicted/observed-1))
pd.DataFrame(endpoint_rows).to_csv(OUT/'line_interpolation_residuals.csv',index=False)

# Algebraic recovery test, no invented MAUVE measurement covariance.
h=np.array([e[0] for e in templates.values()])
d=np.array([e[1]-e[0] for e in templates.values()])
synthetic_w=.37
synthetic_cov=np.diag([.02,.02,.03])**2
cinv=np.linalg.inv(synthetic_cov)
recovered_w=float(d@cinv@(h+synthetic_w*d-h)/(d@cinv@d))
synthetic_recovery_error=abs(recovered_w-synthetic_w)
# Unequal-Balmer conversion and direct luminosity recomputation.
test_wa=.4; test_bh=2.86; test_bd=3.10
test_wb=(test_wa/test_bd)/((1-test_wa)/test_bh+test_wa/test_bd)
test_ro_direct=((1-test_wa)*.57/test_bh+test_wa*.77/test_bd)/((1-test_wa)/test_bh+test_wa/test_bd)
test_ro_mix=(1-test_wb)*.57+test_wb*.77
balmer_identity_error=abs(test_ro_direct-test_ro_mix)

plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans',
    'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,
    'grid.alpha':.16,'savefig.bbox':'tight'})
fig,axs=plt.subplots(1,2,figsize=(10.4,4.0))
axs[0].semilogy(t,hi/HI0,label='H I / initial',color='#5470a5')
axs[0].semilogy(t,f_sfr,label='H2 and SFR / initial',color='#b44b40')
axs[0].semilogy(t,f_control,'--',label='Closed, no RPS control',color='#666666')
axs[0].semilogy(t,np.exp(-GAMMA_CONS*t),':',label='Zero-supply lower bound',color='#333333')
axs[0].set(xlabel='Elapsed model time (Gyr)',ylabel='Fraction of initial state',ylim=(1e-4,1.1),title='H I-only removal; depletion time 2 Gyr')
axs[0].legend(fontsize=8,loc='lower left')
for tau,color in [(2,'#b44b40'),(1,'#5470a5'),(.5,'#52866a')]:
    bc=(1-RETURN)/tau
    tauconv=HI0/(bc*H20)
    axs[1].plot(t,normalized_sfr(t,1/tauconv+GAMMA_STRIP_HI,bc),label=f'Depletion time {tau:g} Gyr',color=color)
axs[1].axhline(.407,color='#666666',ls=':',label='Surviving-SF factor 0.407')
axs[1].set(xlabel='Elapsed model time (Gyr)',ylabel='SFR / initial SFR',ylim=(0,1.03),title='Sensitivity; gas columns held fixed')
axs[1].legend(fontsize=8)
fig.tight_layout();fig.savefig(OUT/'figure01_HI_only_SFR.png',dpi=220);fig.savefig(OUT/'figure01_HI_only_SFR.pdf');plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(10.4,4.0))
axs[0].plot(t[t<=1],l_ha[t<=1]/l_ha[0],label='Total Halpha',color='#333333')
axs[0].plot(t[t<=1],l_hii[t<=1]/l_hii[0],'--',label='H II Halpha',color='#b44b40')
axs[0].plot(t[t<=1],l_non[t<=1]/l_non[0],':',label='Non-H II Halpha',color='#5470a5')
axs[0].set(xlabel='Elapsed gas-model time (Gyr)',ylabel='Luminosity / own initial value',title='Fixed photon allocation: common fading')
axs[0].legend(fontsize=8)
axs[1].semilogy(tl,lt/l0,label='Total Halpha',color='#333333')
axs[1].semilogy(tl,lh/lh[0],label='H II Halpha',color='#b44b40')
axs[1].semilogy(tl,ld/ld[0],label='Non-H II Halpha',color='#5470a5')
axs[1].plot(tl,wl,'--',label='Non-H II Halpha weight',color='#52866a')
axs[1].set(xlabel='Assumed interpolation coordinate (Gyr)',ylabel='Fraction or normalized luminosity',ylim=(.03,1.1),title='NSF Halpha + NII-calibrated example')
axs[1].legend(fontsize=8)
fig.tight_layout();fig.savefig(OUT/'figure02_photon_partition_and_fading.png',dpi=220);fig.savefig(OUT/'figure02_photon_partition_and_fading.pdf');plt.close(fig)

fig,axs=plt.subplots(1,3,figsize=(11.0,3.8),sharex=True)
for ax,(key,(rh,rd)) in zip(axs,templates.items()):
    ax.axhspan(rh,rd,color='#e4e8ed',label='Fixed endpoint interval')
    ax.plot(tl,line_ratios[key],color='#b44b40',label='Conditional interpolation')
    ax.scatter([0,1],[pre[key],post[key]],marker='x',s=65,color='#111111',label='MAUVE NSF means',zorder=4)
    ax.set(xlabel='Interpolation coordinate (Gyr)',ylabel={'nii_ha':'[N II] / Halpha','sii_ha':'[S II] / Halpha','oiii_hb':'[O III] / Hbeta'}[key])
axs[0].legend(fontsize=7,loc='upper left')
fig.tight_layout();fig.savefig(OUT/'figure03_multi_ratio_falsification.png',dpi=220);fig.savefig(OUT/'figure03_multi_ratio_falsification.pdf');plt.close(fig)

summary=dict(c_alpha=C_ALPHA,hi0=HI0,h20=H20,sfr0=SFR0,tau_dep_gyr=TAU_DEP,
    tau_conv_gyr=TAU_CONV,return_fraction=RETURN,eta=ETA,gamma_strip_hi=GAMMA_STRIP_HI,
    gamma_strip_h2=0.0,gamma_loss_hi=GAMMA_HI,gamma_cons_h2=GAMMA_CONS,
    tau_ion_gyr=TAU_ION,f_capture=F_CAPTURE,w_non_null=W_NON_NULL,templates=templates,
    sfr_fraction_1gyr=float(normalized_sfr(1)),sfr_dex_1gyr=float(np.log10(normalized_sfr(1))),
    closed_control_fraction_1gyr=float(normalized_sfr(1,1/TAU_CONV)),
    rps_to_closed_control_1gyr=float(normalized_sfr(1)/normalized_sfr(1,1/TAU_CONV)),
    zero_supply_bound_1gyr=float(np.exp(-GAMMA_CONS)),
    lower_time_for_survivor_factor=float(-np.log(.407)/GAMMA_CONS),
    upper_tau_dep_for_survivor_factor_1gyr=float((1-RETURN)/(-np.log(.407))),
    ha_null_fraction_1gyr=float(np.interp(1,t,l_ha)/l_ha[0]),
    nsf_ha_factor=float(l1/l0),nsf_w0=float(w0),nsf_w1=float(w1),
    effective_gamma_ha_hii=gamma_hii,effective_gamma_ha_non_hii=gamma_non,
    nsf_hii_ha_factor=float(lh[-1]/lh[0]),nsf_non_hii_ha_factor=float(ld[-1]/ld[0]),
    nsf_sii_line_factor=float(line_predictions.sii_surface.iloc[-1]/line_predictions.sii_surface.iloc[0]),
    nsf_oiii_line_factor=float(line_predictions.oiii_surface.iloc[-1]/line_predictions.oiii_surface.iloc[0]),
    required_photon_fractions_initial=[float(budget_hii[0]),float(budget_leak[0]),float(budget_unabs[0])],
    required_photon_fractions_final=[float(budget_hii[-1]),float(budget_leak[-1]),float(budget_unabs[-1])],
    required_photon_fraction_minimum=float(min(budget_hii.min(),budget_leak.min(),budget_unabs.min())),
    ode_relative_error=ode_error,degenerate_ode_relative_error=deg_error,
    initial_sfr_derivative=initial_sfr_derivative,
    initial_sfr_second_derivative=initial_sfr_second_derivative,
    synthetic_weight_recovery_error=synthetic_recovery_error,
    unequal_balmer_identity_error=balmer_identity_error,
    scalar_export_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
assert ode_error<1e-8 and deg_error<1e-8
assert abs(initial_sfr_derivative)<1e-12
assert np.all(np.diff(f_sfr)<=1e-12)
assert np.all(f_sfr>=np.exp(-GAMMA_CONS*t)-1e-12)
assert synthetic_recovery_error<1e-12 and balmer_identity_error<1e-12
assert np.all(np.isfinite(lt)) and np.all((wl>=0)&(wl<=1))
assert np.all(budget_unabs>=0) and np.all(budget_hii>=0) and np.all(budget_leak>=0)
assert np.max(np.abs(budget_hii+budget_leak+budget_unabs-1))<1e-12
assert abs(line_predictions.nii_ha.iloc[0]-pre.nii_ha)<1e-12
assert abs(line_predictions.nii_ha.iloc[-1]-post.nii_ha)<1e-12
(OUT/'numerical_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
