"""New illustrative gas, emission and class calculations; no fitted physical clocks."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, quad
from scipy.special import ndtr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
C_HA=4.983582e-42
LREF=1e38
TAU_ION=5.0
RETURN=.4
PARAMS={
    'inner':dict(HI0=2.4,H20=8.,tau_conv=1000.,tau_HI=160.,tau_H2=650.,
                 tau_dep=2000.,Lcont0=2.,tau_cont=1200.,continuum=1.2),
    'outer':dict(HI0=8.,H20=4.,tau_conv=6666.666666666667,tau_HI=80.,tau_H2=180.,
                 tau_dep=2000.,Lcont0=.05,tau_cont=100.,continuum=.05)}
t=np.linspace(0,1000,501)
s_log=1.1
w_bpt=.45
w_sigma=(45**2-20**2)/(65**2-20**2)
w_limit=min(w_bpt,w_sigma)

def gas(t,p):
    k=1/p['tau_HI']+1/p['tau_conv']
    lam=1/p['tau_H2']+(1-RETURN)/p['tau_dep']
    HI=p['HI0']*np.exp(-k*t)
    if np.isclose(lam,k,rtol=0,atol=1e-13):
        H2=np.exp(-lam*t)*(p['H20']+p['HI0']*t/p['tau_conv'])
    else:
        H2=p['H20']*np.exp(-lam*t)+p['HI0']/p['tau_conv']*(np.exp(-k*t)-np.exp(-lam*t))/(lam-k)
    return HI,H2

def response(t,rate):
    if np.isclose(rate*TAU_ION,1):
        return (1+t/TAU_ION)*np.exp(-t/TAU_ION)
    return (np.exp(-rate*t)-rate*TAU_ION*np.exp(-t/TAU_ION))/(1-rate*TAU_ION)

def young_lum(t,p):
    k=1/p['tau_HI']+1/p['tau_conv']
    lam=1/p['tau_H2']+(1-RETURN)/p['tau_dep']
    c2=p['HI0']/p['tau_conv']/(lam-k)
    c1=p['H20']-c2
    # Gas in Msun/pc2; time Myr -> SFR in Msun/yr/kpc2: factors cancel.
    return (c1*response(t,lam)+c2*response(t,k))/p['tau_dep']/C_HA/LREF

def cdf(lum,mu):
    return ndtr((np.log(lum)-mu)/s_log) if lum>0 else 0.

def moment_above(lum,mu):
    if lum<=0:return np.exp(mu+s_log**2/2)
    return np.exp(mu+s_log**2/2)*ndtr((mu+s_log**2-np.log(lum))/s_log)

rows=[];max_rel=0.
fig,axes=plt.subplots(2,3,figsize=(11.7,7.0),layout='constrained')
for j,(zone,p) in enumerate(PARAMS.items()):
    HI,H2=gas(t,p); lum=young_lum(t,p)
    def rhs(_,y):
        hi,h2,L=y
        return [-hi/p['tau_conv']-hi/p['tau_HI'],
                hi/p['tau_conv']-((1-RETURN)/p['tau_dep']+1/p['tau_H2'])*h2,
                (h2/p['tau_dep']/C_HA/LREF-L)/TAU_ION]
    sol=solve_ivp(rhs,[0,1000],[p['HI0'],p['H20'],lum[0]],t_eval=t,rtol=1e-10,atol=1e-12)
    max_rel=max(max_rel,float(np.max(abs(sol.y-np.array([HI,H2,lum]))/np.maximum(abs(sol.y),1e-8))))
    zone_rows=[]
    for ti,hi,h2,ly in zip(t,HI,H2,lum):
        lc=p['Lcont0']*np.exp(-ti/p['tau_cont'])
        mu=np.log(ly)
        det=max(.8-lc,0.)
        sf=max(det,6*p['continuum']-lc,lc*(1-w_limit)/w_limit,0.)
        nd=cdf(det,mu);fsf=1-cdf(sf,mu);fnsf=1-nd-fsf
        m_sf=moment_above(sf,mu)
        m_nd=np.exp(mu+s_log**2/2)-moment_above(det,mu)
        m_nsf=np.exp(mu+s_log**2/2)-m_sf-m_nd
        I_sf=(m_sf+lc*fsf)/fsf*C_HA*LREF if fsf>0 else np.nan
        I_nsf=(m_nsf+lc*fnsf)/fnsf*LREF if fnsf>1e-10 else np.nan
        zone_rows.append(dict(zone=zone,time_Myr=ti,Sigma_HI=hi,Sigma_H2=h2,
            young_Halpha_Lref=ly,continuing_Halpha_Lref=lc,
            F_SF=fsf,F_NSF=fnsf,F_ND=nd,I_SF_apparent=I_sf,I_NSF_Halpha=I_nsf,
            J_NSF_Halpha=(m_nsf+lc*fnsf)*LREF))
    z=pd.DataFrame(zone_rows);rows.extend(zone_rows)
    axes[j,0].plot(t,HI/p['HI0'],label='Atomic gas')
    axes[j,0].plot(t,H2/p['H20'],label='Molecular gas / true SFR')
    axes[j,0].plot(t,lum/lum[0],ls='--',label='Young-star H-alpha')
    axes[j,0].set(title=f'{zone.capitalize()}: regulator response',ylabel='Fraction of initial value',ylim=(0,1.05))
    axes[j,0].legend(fontsize=7.5)
    for c,color in [('SF','#2166ac'),('NSF','#b2182b'),('ND','#e69f00')]:
        axes[j,1].plot(t,z['F_'+c],label=c,color=color)
    axes[j,1].set(title='Illustrative class occupancy',ylabel='Area fraction',ylim=(0,1))
    axes[j,1].legend(fontsize=8)
    axes[j,2].plot(t,np.log10(z.I_SF_apparent/z.I_SF_apparent.iloc[0]),label='Selected SF: apparent SFR',color='#2166ac')
    axes[j,2].plot(t,np.log10(z.I_NSF_Halpha/z.I_NSF_Halpha.iloc[0]),label='Selected NSF: H-alpha',color='#b2182b')
    axes[j,2].set(title='Selection-conditioned intensity',ylabel='Change from initial value (dex)')
    axes[j,2].legend(fontsize=7.5)
for ax in axes.flat:ax.set_xlabel('Illustrative elapsed time (Myr)')
fig.savefig(OUT/'figure_05_gas_to_classes.png',dpi=200,bbox_inches='tight')
fig.savefig(OUT/'figure_05_gas_to_classes.pdf',bbox_inches='tight');plt.close(fig)
df=pd.DataFrame(rows);df.to_csv(OUT/'analytical_predictions.csv',index=False)
df[df.time_Myr.isin([0,200,600,1000])].to_csv(OUT/'analytical_prediction_landmarks.csv',index=False)
assert np.allclose(df[['F_SF','F_NSF','F_ND']].sum(axis=1),1)
assert (df[['F_SF','F_NSF','F_ND']]>=-1e-12).all().all()

# Distinct experiment: exact line mixing with both physical components fading.
tm=np.linspace(0,600,301);ly=8*np.exp(-tm/200);lc=2*np.exp(-tm/700)
w=lc/(ly+lc)
mix=pd.DataFrame(dict(time_Myr=tm,Halpha=ly+lc,w_cont=w,
                     NII_Halpha=.25*(1-w)+1.*w,
                     SII_Halpha=.18*(1-w)+.65*w,
                     OIII_Hbeta=.6*(1-w)+.3*w))
mix.to_csv(OUT/'line_mixing_predictions.csv',index=False)
fig,axs=plt.subplots(1,3,figsize=(11,3.6),layout='constrained')
axs[0].plot(tm,ly+lc,label='Total H-alpha');axs[0].plot(tm,ly,label='Young-star component')
axs[0].plot(tm,lc,label='Continuing component')
axs[0].set(ylabel='Luminosity / reference luminosity');axs[0].legend(fontsize=8)
for name in ['NII_Halpha','SII_Halpha','OIII_Hbeta']:axs[1].plot(tm,mix[name],label=name.replace('_',' / '))
axs[1].set(ylabel='Linear line ratio');axs[1].legend(fontsize=8)
axs[2].plot(np.log10(mix.NII_Halpha),np.log10(mix.OIII_Hbeta),color='#6a3d9a')
axs[2].scatter(np.log10(mix.NII_Halpha.iloc[[0,-1]]),np.log10(mix.OIII_Hbeta.iloc[[0,-1]]),c=['#2166ac','#b2182b'])
axs[2].set(xlabel='log10 [N II] / H-alpha',ylabel='log10 [O III] / H-beta',
           title='A rightward, downward mixing path')
for ax in axs[:2]:ax.set_xlabel('Illustrative elapsed time (Myr)')
fig.savefig(OUT/'figure_06_line_mixing.png',dpi=200,bbox_inches='tight')
fig.savefig(OUT/'figure_06_line_mixing.pdf',bbox_inches='tight');plt.close(fig)

# Independent quadrature validates the truncated lognormal moment.
mu=np.log(3.);threshold=2.7
pdf=lambda u:np.exp(-.5*((u-mu)/s_log)**2)/(s_log*np.sqrt(2*np.pi))
num=quad(lambda u:np.exp(u)*pdf(u),np.log(threshold),mu+14*s_log,epsabs=1e-10)[0]
moment_error=abs(num-moment_above(threshold,mu))/num
assert max_rel<1e-7 and moment_error<1e-9
checks=dict(max_relative_ODE_solution_error=max_rel,truncated_moment_relative_error=moment_error,
            category_partition_max_error=float(abs(df[['F_SF','F_NSF','F_ND']].sum(axis=1)-1).max()),
            parameters=PARAMS,Lref_erg_s_kpc2=LREF,C_Halpha=C_HA,tau_ion_Myr=TAU_ION,
            lognormal_sigma_ln=s_log,w_BPT_proxy=w_bpt,w_sigma=w_sigma,w_limit=w_limit,
            Halpha_detection_Lref=.8,return_fraction=RETURN,physical_parameters_fitted=False)
(OUT/'analytical_checks.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
