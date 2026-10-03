"""Atomic-only initial-column family for report sections 3.5 and 8.3.

3 October 2026: replace the old simultaneous HI/H2 scaling experiment.
Read the saved reference normalization; do not re-extract observational data
or regenerate any spectral results. Factors are illustrative, not fitted.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp, quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parent / '20261002_Model0_Derivation/numerical_audit.json'
source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
saved = json.loads(SOURCE.read_text())
HI0, H20 = saved['hi_0'], saved['h2_0']
CONV, DEP = saved['tau_conv'], saved['tau_dep']
GHI, GH2 = saved['cases']['balanced']['gamma_hi'], saved['gamma_h2']
SFR0 = saved['sfr0']
assert np.isclose(HI0/CONV, GH2*H20, rtol=1e-13)
assert np.isclose(1e-3*H20/DEP, SFR0, rtol=1e-13)

def molecular(t, factor, hi0=HI0, h20=H20, conv=CONV, ghi=GHI, gh2=GH2):
    """Molecular column including the common surviving initial reservoir."""
    t = np.asarray(t)
    if ghi == gh2:
        supplied = hi0/conv*t*np.exp(-gh2*t)
    else:
        supplied = hi0/conv*(np.exp(-ghi*t)-np.exp(-gh2*t))/(gh2-ghi)
    return h20*np.exp(-gh2*t)+factor*supplied

def peak_time(factor):
    phi_time = H20*CONV/(factor*HI0)
    if 1/phi_time <= GH2:
        return 0.0
    if GHI == GH2:
        return 1/GH2-phi_time
    return float(np.log(GHI/(GH2*(1+(GHI-GH2)*phi_time)))/(GHI-GH2))

t = np.linspace(0, 2, 2001)
reference = molecular(t, 1.)
cases = [('leading',1.5), ('reference',1.), ('trailing',.5)]
curves, summaries, rows, errors = {}, [], [], []
for region, factor in cases:
    h2 = molecular(t, factor)
    hi = factor*HI0*np.exp(-GHI*t)
    numerical = solve_ivp(lambda time,y:[-GHI*y[0], y[0]/CONV-GH2*y[1]],
        (0,2),[factor*HI0,H20],t_eval=t,method='DOP853',rtol=2e-12,atol=2e-13)
    assert numerical.success
    error = float(np.max(np.abs(numerical.y[1]-h2))/H20)
    assert error < 1e-10
    errors.append(error)
    frac = h2/H20
    offset = np.log10(h2/reference)
    peak = peak_time(factor)
    if peak>0:
        assert abs(factor*HI0/CONV*np.exp(-GHI*peak)-GH2*molecular(peak,factor))<1e-11
    curves[region] = (factor,frac,offset)
    summaries.append(dict(region=region,c_HI=factor,initial_HI=factor*HI0,
        initial_H2=H20,initial_SFR=SFR0,tau_phi_initial=H20*CONV/(factor*HI0),
        initial_fractional_SFR_slope=(factor-1)*GH2,
        t_peak_gyr=peak,peak_fraction=float(molecular(peak,factor)/H20),
        sfr_fraction_1gyr=float(frac[1000]),offset_dex_1gyr=float(offset[1000]),
        sfr_1gyr=float(1e-3*h2[1000]/DEP),h2_1gyr=float(h2[1000]),
        offset_dex_2gyr=float(offset[-1])))
    for time,atomic,mol,f,delta in zip(t,hi,h2,frac,offset):
        rows.append(dict(region=region,c_HI=factor,t_gyr=time,hi=atomic,h2=mol,
            sfr=1e-3*mol/DEP,sfr_fraction=f,offset_dex=delta))

assert all(curves[r][1][0]==1 for r,_ in cases)
assert np.all(curves['leading'][1][1:]>curves['reference'][1][1:])
assert np.all(curves['reference'][1][1:]>curves['trailing'][1][1:])
assert np.all(np.diff(curves['reference'][1])<0)
assert np.all(np.diff(curves['trailing'][1])<0)
# A direct counterexample to multiplying the full reference H2: the true
# difference is zero at onset, whereas that incorrect expression is not.
assert molecular(0,1.5)-molecular(0,.5)==0
assert (1.5-.5)*molecular(0,1.)>0
expected=(1.5-.5)*(reference-H20*np.exp(-GH2*t))
assert np.allclose(molecular(t,1.5)-molecular(t,.5),expected,rtol=1e-12,atol=1e-13)

other_errors=[]
for ghi,gh2,conv in [(0.2,0.8,10.),(.5,.5,4.)]:
    for factor in [.5,1.,1.5]:
        sol=solve_ivp(lambda time,y:[-ghi*y[0],y[0]/conv-gh2*y[1]],
            (0,2),[factor*HI0,H20],t_eval=t,method='DOP853',rtol=2e-12,atol=2e-13)
        analytic=molecular(t,factor,conv=conv,ghi=ghi,gh2=gh2)
        error=float(np.max(np.abs(sol.y[1]-analytic))/H20)
        assert error<1e-10
        other_errors.append(error)
        integral=quad(lambda u:HI0/conv*np.exp(-ghi*u)*np.exp(-gh2*(1-u)),0,1)[0]
        assert np.isclose(molecular(1,factor,conv=conv,ghi=ghi,gh2=gh2),
                          H20*np.exp(-gh2)+factor*integral,rtol=1e-13)

for filename,records in [('spatial_curves.csv',rows),('spatial_summary.csv',summaries)]:
    with (OUT/filename).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(records[0]))
        writer.writeheader();writer.writerows(records)

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
colors={'leading':'#b75522','reference':'#333333','trailing':'#24628b'}
fig,axes=plt.subplots(1,2,figsize=(9,3.8),layout='constrained')
for region,_ in cases:
    factor,frac,offset=curves[region]
    label=region.capitalize()+rf': $c_{{\rm HI}}={factor:g}$'
    axes[0].plot(t,frac,color=colors[region],label=label,lw=1.7)
    axes[1].plot(t,offset,color=colors[region],lw=1.7)
axes[0].axhline(1,color='#999999',ls=':',lw=.8)
axes[0].set(xlabel='Time since model onset [Gyr]',
    ylabel=r'$\Sigma_{\rm SFR}(x,t)/\Sigma_{\rm SFR,0}$',xlim=(0,2),ylim=(.49,1.08))
axes[0].legend(frameon=False,fontsize=8,loc='lower left')
axes[1].set(xlabel='Time since model onset [Gyr]',
    ylabel=r'$\log_{10}[\Sigma_{\rm SFR}(x,t)/\Sigma_{\rm SFR}^{\rm reference}(t)]$',
    xlim=(0,2),ylim=(-.020,.020))
axes[1].text(.06,.96,'Matched reference has the same HI stripping',
    transform=axes[1].transAxes,va='top',fontsize=8)
# Show the small initial rise without implying a large starburst.
inset=axes[0].inset_axes([.53,.60,.44,.33])
for region,_ in cases:
    inset.plot(t,curves[region][1],color=colors[region],lw=1.1)
inset.axhline(1,color='#999999',ls=':',lw=.6)
inset.set(xlim=(0,.22),ylim=(.967,1.010),xticks=[0,.1,.2],yticks=[.98,1.00],title='Early response [Gyr]')
inset.tick_params(labelsize=7);inset.title.set_fontsize(8)
for extension in ['png','pdf']:
    fig.savefig(OUT/f'figure01_atomic_spatial_response.{extension}',dpi=190)
plt.close(fig)

audit=dict(source=str(SOURCE),source_sha256=source_hash,
    constants=dict(hi0=HI0,h20=H20,tau_conv=CONV,tau_dep=DEP,
                   gamma_HI=GHI,gamma_H2=GH2,sfr0=SFR0),cases=summaries,
    maximum_reference_normalized_H2_ODE_error=max(errors+other_errors),
    checks=dict(initial_equality=True,strict_spatial_ordering=True,
                corrected_difference_identity=True,equal_rate_limit=True,
                reversed_rate_ordering=True,peak=True,positive_supply_integral=True),
    classification='Illustrative atomic factors; saved MAUVE luminosity normalization; not a fit')
(OUT/'spatial_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
