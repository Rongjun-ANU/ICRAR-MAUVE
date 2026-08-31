"""Original mathematical checks and schematic figures for the 2026-08-31 FMR report.

No observed or simulation catalogue is used. All curves are illustrative.
Run with MPLBACKEND=Agg and a writable MPLCONFIGDIR.
"""
from pathlib import Path
import json
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
mp.mp.dps = 55
R, a, y = 0.44, 0.56, 0.012
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 150,
                     'savefig.dpi': 220, 'font.family': 'DejaVu Sans'})
blue, orange, green, purple = '#23638c', '#c77832', '#35765a', '#835994'


def F(x):
    x = np.asarray(x)
    return 1 - x / np.expm1(x)


def G(x):
    return np.asarray(x) - F(x)


def B(x):
    x = np.asarray(x)
    fp = ((x - 1) * np.exp(x) + 1) / np.expm1(x)**2
    return G(x) * fp / (F(x) * (1-fp))


def invG(v):
    return brentq(lambda x: G(x)-v, 1e-8, max(4., v+2.), xtol=1e-13)


def save(fig, name):
    fig.savefig(OUT / (name+'.png'), bbox_inches='tight')
    fig.savefig(OUT / (name+'.pdf'), bbox_inches='tight')
    plt.close(fig)


checks = {}
# Independent integration of gas, stars, gas metals, stellar metals, and losses.
eta, eps, phi = 0.22, 0.38, 1.0
A = a + eta
tau = 1/(A*eps)
times = np.geomspace(1e-3, 15, 700)*tau


def rhs(t, v):
    gas, stars, metals, stellar_metals, out, outmetals = v
    sfr = eps*gas
    z = metals/gas if gas > 0 else 0.
    return [phi-A*sfr, a*sfr, y*sfr-A*z*sfr,
            a*z*sfr, eta*sfr, eta*z*sfr]


sol = solve_ivp(rhs, (0, times[-1]), np.zeros(6), t_eval=times,
                rtol=2e-11, atol=1e-13)
assert sol.success
x = times/tau
g = phi*tau*(-np.expm1(-x))
stars = a/A*phi*tau*(x+np.expm1(-x))
z = y/A*F(x)
checks['ideal_gas_max_scaled_error'] = float(np.max(abs(sol.y[0]-g)/(phi*tau)))
checks['ideal_stars_max_scaled_error'] = float(np.max(abs(sol.y[1]-stars)/(a/A*phi*tau)))
checks['ideal_Z_max_yield_scaled_error'] = float(np.max(abs(sol.y[2]/sol.y[0]-z)/(y/A)))
checks['mass_budget_max_relative_error'] = float(np.max(abs(sol.y[0]+sol.y[1]+sol.y[4]-phi*times)/(phi*times)))
checks['metal_budget_max_relative_error'] = float(np.max(abs(sol.y[2]+sol.y[3]+sol.y[5]-y/a*sol.y[1])/(y/a*sol.y[1])))
Zstar = y/A*(x-2+(x+2)*np.exp(-x))/(x-1+np.exp(-x))
checks['stellar_Z_max_yield_scaled_error_x_gt_0p01'] = float(np.max(abs(sol.y[3,x>.01]/sol.y[1,x>.01]-Zstar[x>.01])/(y/A)))

# Differential-loading exact solution compared against metal-mass integration.
errs = []
for eta_d, zeta_d in [(0.8, 0.2), (0.8, 0.8), (0.8, 2.4)]:
    Am, Cz = a+eta_d, a+zeta_d
    rr = Cz/Am
    td = np.linspace(.001, 12, 900)/(Am*eps)
    def rd(t, v):
        return [phi-Am*eps*v[0], y*eps*v[0]-Cz*eps*v[1]]
    sd = solve_ivp(rd, (0,td[-1]), [0,0], t_eval=td, rtol=2e-11, atol=1e-13)
    xd = Am*eps*td
    if abs(rr-1)<1e-10:
        zz = y/Am*F(xd)
    else:
        hh = -np.expm1(-rr*xd)/rr - (np.exp(-xd)-np.exp(-rr*xd))/(rr-1)
        zz = y/Am*hh/(-np.expm1(-xd))
    errs.append(float(np.max(abs(sd.y[1]/sd.y[0]-zz)/(y/Cz))))
checks['differential_loading_max_scaled_errors'] = errs

# High-precision verification of the slope identity, independent finite differences.
def fm(xx): return 1-xx/mp.expm1(xx)
def gm(xx): return xx-fm(xx)
slope_errors=[]
for xx in [.01,.1,1.,3.,10.,30.]:
    xx=mp.mpf(xx)
    beta0=gm(xx)*mp.diff(fm,xx)/(fm(xx)*mp.diff(gm,xx))
    for be, ca in [(0,0),(1.5,0),(.3,-.5),(0,.4)]:
        def z_of_s(ss):
            xx2=mp.findroot(lambda w: gm(w)-gm(xx)*mp.exp((ca+be-1)*ss), (xx*.99,xx*1.01))
            return mp.exp(-ca*ss)*fm(xx2)
        actual=mp.diff(lambda ss: mp.log(z_of_s(ss)),0)
        predicted=beta0*(be-1)+(beta0-1)*ca
        slope_errors.append(float(abs(actual-predicted)))
checks['slope_identity_max_absolute_error'] = max(slope_errors)
checks['B_test_range'] = [float(B(.001)),float(B(50.))]

# Counterexample to arbitrary-history uniqueness: equal Mg, Mstar, eps, eta, SFR.
u1,u2,u3=.2,1.,2.
w=(np.exp(-u2)-np.exp(-u3))/(np.exp(-u1)-np.exp(-u3))
gas_a=np.exp(-u2); gas_b=w*np.exp(-u1)+(1-w)*np.exp(-u3)
met_a=y*u2*np.exp(-u2)
met_b=y*(w*u1*np.exp(-u1)+(1-w)*u3*np.exp(-u3))
checks['history_counterexample']={'early_late_mixture_weight_at_age_0p2':float(w),
    'same_gas_error':float(abs(gas_a-gas_b)), 'Z_single_pulse':float(met_a/gas_a),
    'Z_two_pulses':float(met_b/gas_b),'delta_log10_Z':float(np.log10((met_b/gas_b)/(met_a/gas_a)))}

# Figure 1: full ideal solution, with analytical limiting behaviour.
x=np.geomspace(.01,25,700); ff=F(x); gg=G(x)
fig,ax=plt.subplots(2,2,figsize=(8.4,6.1),layout='constrained')
ax[0,0].loglog(x,-np.expm1(-x),color=blue,label='Gas mass')
ax[0,0].loglog(x,x+np.expm1(-x),color=orange,label='Stellar mass')
ax[0,0].set_ylabel('Mass / its stated normalisation');ax[0,0].legend()
ax[0,1].semilogx(x,ff,color=blue,label=r'$Z_g/(y/A)$')
ax[0,1].semilogx(x,(x-2+(x+2)*np.exp(-x))/(x-1+np.exp(-x)),color=orange,label=r'$Z_\star/(y/A)$')
ax[0,1].set_ylabel('Normalised metallicity');ax[0,1].legend()
ax[1,0].semilogx(x,ff/gg,color=blue,label=r'$K_1$')
ax[1,0].semilogx(x,(1-ff/gg)/gg,color=orange,label=r'$K_3$')
ax[1,0].set_ylabel('Metal-retention functions');ax[1,0].legend()
ax[1,1].semilogx(x,B(x),color=purple)
ax[1,1].set_ylabel(r'$B=d\ln F/d\ln G$');ax[1,1].set_ylim(-.03,1.05)
for aa in ax.flat:
    aa.set_xlabel(r'Evolutionary stage $x=t/\tau$');aa.grid(alpha=.18)
fig.suptitle('Ideal gas-flow solution: original calculation, no data fit',fontsize=12)
save(fig,'figure_01_ideal_solution')

# Figure 2: gas fraction and elasticity.
mus=np.geomspace(.02,30,600)
fig,ax=plt.subplots(1,2,figsize=(8.4,3.8),layout='constrained')
for et,color in [(0,blue),(.5,orange),(2,purple)]:
    AA=a+et
    xs=np.array([invG(AA/a/mu) for mu in mus])
    zz=y/AA*F(xs)
    ax[0].loglog(mus,zz/(y/a),color=color,label=rf'$\eta={et:g}$')
    ax[1].semilogx(mus,B(xs),color=color)
ax[0].set_ylabel(r'$Z_g/(y/a)$');ax[0].legend()
ax[1].set_ylabel(r'$B$ (baseline SFR slope is $-B$)')
for aa in ax:
    aa.set_xlabel(r'Gas-to-stellar mass ratio $\mu=M_g/M_\star$');aa.grid(alpha=.18)
save(fig,'figure_02_gfmr_and_slope')

# Figure 3: fixed-final-Mstar, exact ideal-model positive family.
Mstar=10**10.7; psi0=.5*Mstar/1e9; eps0=.5; et=.44; AA=a+et
rel=np.linspace(.8,1.6,160)
rows=[]
fig,ax=plt.subplots(1,2,figsize=(8.4,3.8),layout='constrained')
for be,col in [(0,blue),(1,green),(1.5,orange)]:
    sfr=psi0*rel; ee=eps0*rel**be
    mg=sfr*1e9/ee; mu=mg/Mstar
    xs=np.array([invG(AA/a/mm) for mm in mu])
    zz=y/AA*F(xs); age=xs/(AA*ee)
    inf=AA*sfr/(-np.expm1(-xs))
    ax[0].plot(np.log10(rel),np.log10(zz),color=col,label=rf'$b={be:g}$')
    ax[1].plot(np.log10(rel),mu,color=col)
    for i in range(len(rel)):rows.append([be,sfr[i],ee[i],mu[i],zz[i],age[i],inf[i]])
ax[0].set_ylabel(r'$\log_{10}Z_g$');ax[0].legend()
ax[1].set_ylabel(r'$M_g/M_\star$')
for aa in ax:
    aa.set_xlabel(r'$\log_{10}(\mathrm{SFR}/\mathrm{SFR}_0)$');aa.grid(alpha=.18)
fig.suptitle(r'Exact fixed-mass families: $M_\star=10^{10.7}M_\odot$, $\eta=0.44$',fontsize=11)
save(fig,'figure_03_fixed_mass_families')
np.savetxt(OUT/'fixed_mass_families.csv',np.asarray(rows),delimiter=',',comments='',
           header='b,SFR_Msun_per_yr,epsilon_per_Gyr,gas_to_stellar_mass,Z_mass_fraction,age_Gyr,inflow_Msun_per_yr')
positive=np.asarray(rows)[np.asarray(rows)[:,0]==1.5]
checks['positive_exact_family']={'SFR_range':positive[[0,-1],1].tolist(),
 'Z_range':positive[[0,-1],4].tolist(),'age_Gyr_range':[float(positive[:,5].min()),float(positive[:,5].max())],
 'Mstar_Msun':Mstar,'minimum_dZ_dSFR':float(np.min(np.diff(positive[:,4])/np.diff(positive[:,1])))}

# Figure 4: sign map from the exact derivative identity.
bs=np.linspace(-.1,2,220);cs=np.linspace(-1.2,.8,220)
bb,cc=np.meshgrid(bs,cs)
fig,axes=plt.subplots(1,2,figsize=(8.4,3.8),layout='constrained')
for xx,aa in zip([.5,4.],axes):
    bv=float(B(xx));slope=bv*(bb-1)+(bv-1)*cc
    im=aa.contourf(bb,cc,slope,levels=np.linspace(-1.3,1.3,27),cmap='RdBu_r',extend='both')
    aa.contour(bb,cc,slope,levels=[0],colors='black',linewidths=1.2)
    aa.set_title(rf'$x={xx:g}$, $B={bv:.2f}$')
    aa.set_xlabel(r'$b=\partial\ln\epsilon/\partial\ln\mathrm{SFR}$')
axes[0].set_ylabel(r'$c=\partial\ln A/\partial\ln\mathrm{SFR}$')
fig.colorbar(im,ax=axes,label=r'Conditional slope $\beta$')
save(fig,'figure_04_sign_conditions')

# Figure 5: forced non-equilibrium ODE, same continuum equations.
# Dimensionless time t/tau; gas and metal mass are scaled to equilibrium.
omega=1.;amp=.03;tt=np.linspace(0,20+2*np.pi,6000)
fig,axes=plt.subplots(1,2,figsize=(8.4,3.8),layout='constrained')
dynamic={}
for mode,col in [('inflow',blue),('efficiency',orange)]:
    def forc(t,v):
        p=amp*np.cos(omega*t) if mode=='inflow' else 0.
        e=amp*np.cos(omega*t) if mode=='efficiency' else 0.
        return [1+p-(1+e)*v[0], (1+e)*(v[0]-v[1])]
    sd=solve_ivp(forc,(0,tt[-1]),[1,1],t_eval=tt,rtol=2e-11,atol=1e-13)
    sel=tt>=20
    ep=amp*np.cos(tt[sel]) if mode=='efficiency' else np.zeros(sel.sum())
    ss=sd.y[0,sel]*(1+ep);zz=sd.y[1,sel]/sd.y[0,sel]
    lx=np.log(ss);lz=np.log(zz)
    slope=np.cov(lx,lz,bias=True)[0,1]/np.var(lx)
    linear=(-omega**2 if mode=='inflow' else 1)/(1+omega**2)
    dynamic[mode]={'nonlinear_log_slope':float(slope),'linear_prediction':linear}
    axes[0].plot(lx/np.log(10),lz/np.log(10),color=col,label=mode)
    j=120
    axes[0].annotate('',xy=(lx[j+20]/np.log(10),lz[j+20]/np.log(10)),
                     xytext=(lx[j]/np.log(10),lz[j]/np.log(10)),arrowprops={'arrowstyle':'->','color':col})
oms=np.geomspace(.03,30,400)
axes[1].semilogx(oms,-oms**2/(1+oms**2),color=blue,label='inflow only')
axes[1].semilogx(oms,1/(1+oms**2),color=orange,label='efficiency only')
axes[0].set_xlabel(r'$\Delta\log_{10}\mathrm{SFR}$');axes[0].set_ylabel(r'$\Delta\log_{10}Z$')
axes[0].legend();axes[0].set_title(r'Nonlinear ODE, $\Omega=1$, 3% forcing')
axes[1].set_xlabel(r'$\Omega=\omega\tau$');axes[1].set_ylabel('Regression slope');axes[1].legend()
for aa in axes:aa.grid(alpha=.18)
save(fig,'figure_05_driver_response')
checks['dynamic_response']=dynamic

# Figure 6: exact mass-controlled statistical example of pooling confusion.
rng=np.random.default_rng(260831)
mm=np.repeat([9.5,10.5,11.0],300)
res=rng.normal(0,.12,len(mm))
ss=.8*(mm-10)+res
zz=8.65+.25*(mm-10)-.35*res
fig,ax=plt.subplots(figsize=(6.8,3.5),layout='constrained')
for mass,col in zip([9.5,10.5,11.0],[blue,orange,purple]):
    ii=mm==mass
    ax.scatter(ss[ii],zz[ii],s=9,alpha=.5,color=col,label=rf'$\log M_\star={mass:g}$')
ax.set_xlabel(r'$\log_{10}\mathrm{SFR}$ (arbitrary zero point)')
ax.set_ylabel(r'$12+\log_{10}(\mathrm{O/H})$ (synthetic)')
ax.legend();ax.grid(alpha=.18)
save(fig,'figure_06_conditioning_example')
checks['conditioning_example']={'pooled_slope':float(np.polyfit(ss,zz,1)[0]),
    'within_exact_mass_slope':float(np.polyfit(res,zz-(8.65+.25*(mm-10)),1)[0])}

assert checks['ideal_Z_max_yield_scaled_error']<1e-7
assert checks['mass_budget_max_relative_error']<1e-8
assert checks['metal_budget_max_relative_error']<1e-8
assert max(errs)<1e-7
assert checks['slope_identity_max_absolute_error']<1e-12
assert positive[-1,4]>positive[0,4]
assert all(abs(v['nonlinear_log_slope']-v['linear_prediction'])<.002 for v in dynamic.values())
(OUT/'numerical_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
