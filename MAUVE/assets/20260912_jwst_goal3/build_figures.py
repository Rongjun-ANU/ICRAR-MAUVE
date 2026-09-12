"""Rebuild illustrative figures and equation images with the ICRAR Python environment."""
from pathlib import Path
import json
import os
import re
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/mauve_goal3_20260912/matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ASSETS = Path(__file__).resolve().parent
ROOT = ASSETS.parent.parent
SCRATCH = Path('/private/tmp/mauve_goal3_20260912')
SCRATCH.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.titleweight': 'bold', 'savefig.dpi': 240})
MP = 1.67262192369e-24
KB = 1.380649e-16
G = 6.67430e-8
PC = 3.085677581491367e18
MSUN = 1.988409870698051e33
D = 16.2e6 * PC
scale = 16.2e6 / 206264.806247
pram = .6 * MP * 1e-4 * 1e8**2 / KB
grav = np.pi * G * (100 * MSUN / PC**2)**2 / (2*KB)
flux = 1e37 / (4*np.pi*D**2)
assert abs(pram/7270-1)<.002
assert abs(grav/3.31e5-1)<.005
assert abs(scale/78.54-1)<.001
assert abs(flux/3.18e-16-1)<.002

fig, ax = plt.subplots(figsize=(8.1,3.8), layout='constrained')
labels=['ICM thermal pressure\n(chosen n and T)', 'ICM ram stress\n(chosen n and v)',
        'Ionized hydrogen\n(n_e = 10-100 cm$^{-3}$)', 'Self-gravitating slab\n(Sigma = 100 solar masses/pc$^2$)']
intervals=[(2e3,2e4),(pram,pram*10*1.5**2),(2e5,2e6),(grav,grav)]
for i,(lo,hi) in enumerate(intervals):
 ax.plot([lo,hi],[i,i],lw=10,solid_capstyle='round',color=['#6a8f9c','#c2753e','#285e76','#596659'][i])
 ax.plot([lo,hi],[i,i],'o',ms=6,color='#243945')
ax.set_yticks(range(4),labels);ax.invert_yaxis();ax.set_xscale('log');ax.set_xlim(1e3,4e6)
ax.set_xlabel(r'Pressure / $k_B$ [K cm$^{-3}$]');ax.grid(axis='x',alpha=.2)
ax.set_title('Pressure comparison: illustrative inputs, not MAUVE measurements',fontsize=11,pad=14)
fig.savefig(ASSETS/'pressure_comparison.png',bbox_inches='tight');plt.close(fig)

labels=['Pa-alpha (reference)','F335M PAH','F405N Br-alpha','F770W','F1000W','F1130W','F2100W','MUSE (illustrative)']
fwhm=np.array([.060,.111,.136,.301,.356,.390,.685,1.0])*scale
fig,ax=plt.subplots(figsize=(8.2,4.35),layout='constrained')
y=np.arange(len(labels));ax.barh(y,fwhm,color='#285e76',height=.62,label='Reference FWHM')
ax.scatter(2*fwhm,y,color='#c2753e',marker='|',s=130,label='Two-FWHM diameter')
for i,x in enumerate(fwhm):ax.text(x+1.8,i,f'{x:.1f} pc',va='center',fontsize=9)
ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlim(0,175);ax.set_xlabel('Physical scale at 16.2 Mpc [pc]')
ax.legend(loc='lower right',frameon=False,fontsize=9);ax.grid(axis='x',alpha=.17)
ax.set_title('The experiment has different resolutions for different diagnostics',fontsize=11,pad=12)
fig.savefig(ASSETS/'resolution_limits.png',bbox_inches='tight');plt.close(fig)

n=np.arange(3,31);delta=2.8*.15*np.sqrt(2/n)
fig,ax=plt.subplots(figsize=(7.7,3.6),layout='constrained')
ax.plot(n,delta,color='#285e76',lw=2.3)
ax.axhline(.20,color='#c2753e',ls='--',label='Illustrative target: 0.20 fraction difference')
for ng in [5,9,15,20]:
 val=2.8*.15*np.sqrt(2/ng);ax.plot(ng,val,'o',color='#285e76');ax.annotate(f'{val:.2f}',(ng,val),xytext=(4,8),textcoords='offset points',fontsize=9)
ax.set(xlabel='Independent galaxies per comparison group',ylabel='Detectable absolute fraction difference',ylim=(.08,.40),xlim=(3,30))
ax.grid(alpha=.18);ax.legend(frameon=False,fontsize=9,loc='upper right')
ax.set_title('80% power, 5% two-sided test; assumed galaxy scatter = 0.15',fontsize=11,pad=12)
fig.savefig(ASSETS/'power_fraction.png',bbox_inches='tight');plt.close(fig)

equations={}
for md in ROOT.glob('20260912_MAUVE_JWST_*.md'):
 text=md.read_text()
 text=text.translate(str.maketrans({'\u201c':'"','\u201d':'"','\u2018':"'",'\u2019':"'",'\u2013':'-','\u2014':'-','\u00a0':' '}))
 md.write_text(text)
 eqs=re.findall(r'\$\$\s*\n(.*?)\n\$\$',text,re.S)
 out=[]
 for i,eq in enumerate(eqs):
  path=SCRATCH/f'{md.stem}_eq_{i}.png'
  fig=plt.figure(figsize=(10,1));fig.text(.01,.5,'$'+eq.strip()+'$',fontsize=17,va='center')
  fig.savefig(path,bbox_inches='tight',pad_inches=.10,dpi=250,transparent=True);plt.close(fig)
  out.append(str(path))
 equations[md.stem]=out
record={'distance_Mpc':16.2,'pc_per_arcsec':scale,'ram_pressure_over_k':pram,'slab_pressure_over_k':grav,
        'line_flux_benchmark':flux,'ratio_fractional_error_snr10':np.sqrt(2)/10,
        'differential_attenuation_error_mag':2.5/np.log(10)*np.sqrt(2)/10,
        'fraction_difference_n9':2.8*.15*np.sqrt(2/9),'equations':equations}
(SCRATCH/'build_data.json').write_text(json.dumps(record,indent=2))
print(json.dumps({k:v for k,v in record.items() if k!='equations'},indent=2))
