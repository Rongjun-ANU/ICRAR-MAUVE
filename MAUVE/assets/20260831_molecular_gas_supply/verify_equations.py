"""Dimensional/numerical checks and original analytical illustrations.

No observational data are read. Run with the ICRAR Python environment.
These checks establish arithmetic and conservation identities, not the
validity of a cloud-chemistry or stellar-history closure in Virgo.
"""
from pathlib import Path
from fractions import Fraction
import json
import math
import platform

import numpy as np
import astropy
from astropy import units as u
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


OUT = Path(__file__).resolve().parent
RATE = u.Msun / u.pc**2 / u.Myr
MH = 1.6735575e-24 * u.g
MU = 1.36
checks = []


def check(name, condition, detail):
    if not bool(condition):
        raise AssertionError(f"{name}: {detail}")
    checks.append({"name": name, "passed": True, "detail": detail})


def close(a, b, rtol=1e-11, atol=1e-12):
    return np.allclose(a, b, rtol=rtol, atol=atol)


unit_equal = (1 * u.Msun / u.yr / u.kpc**2).to_value(RATE)
check("surface_rate_units", unit_equal == 1, unit_equal)

for j in range(1, 257):
    # Rational arithmetic: exact identity, not a floating-point cancellation.
    F, D = Fraction(j, 17), Fraction(j + 7, 23)
    Bin, Bout = Fraction(j + 3, 31), Fraction(j + 1, 29)
    E, psi = Fraction(j + 2, 41), Fraction(j + 5, 37)
    R, Lpi = Fraction(2, 5), Fraction(j + 4, 43)
    mdot = F - D + Bin - Bout + E - psi
    phi = mdot + (1 - R) * psi + Lpi
    dictionary = F + Bin + E - R * psi + Lpi - D - Bout
    assert phi == dictionary
    assert mdot + psi - E == F - D + Bin - Bout
    q = Fraction(j, 11)
    assert (F + q) - (D + q) == F - D
check("exact_budget_dictionary", True, "256 rational-valued budgets; Eq7-9 and Eq21")

psi, R = Fraction(1, 50), Fraction(2, 5)
phi_no_loss = (1 - R) * psi
tau_ratio = psi / phi_no_loss
check("paper1_timescale_counterexample", tau_ratio == Fraction(5, 3), float(tau_ratio))

n = 100 / u.cm**3
Rd = 3e-17 * u.cm**3 / u.s
tform = (1 / (2 * Rd * n)).to_value(u.Myr)
check("chemical_factor_two", close(tform, 5.281347969004825, rtol=1e-7), tform)

t = np.linspace(0, 200, 501)
xeq, x0 = .8, .1
for tau in [4.0, 40.0]:
    a, b = xeq / tau, (1 - xeq) / tau
    x = xeq + (x0 - xeq) * np.exp(-t / tau)
    derivative = -(x0 - xeq) / tau * np.exp(-t / tau)
    assert close(derivative, a * (1 - x) - b * x)
    assert close(a / (a + b), xeq)
check("chemical_solution", True, "501 time samples at each of two rates with identical equilibrium")

sm, sa, C, sf = Fraction(11), Fraction(29), Fraction(7, 10), Fraction(1, 5)
sn = sm + sa
fdot = ((C - sf) * sn + sm * sf) / sn**2
check("neutral_fraction_budget", sn * fdot + (1 - sm / sn) * sf == C, "Eq18-20 exact")
check("neutral_phase_cancellation", (C - sf) + (-C) == -sf, "Chemistry cancels from total gas")

tau, psi0, a, b = 1000.0, .02, .01, -.001
deriv = tau * psi0 * (a + b)
check("gas_response_product_rule", close(deriv, tau * (a * psi0) + psi0 * (b * tau)), deriv)

Z, Zin, yield_p, psi, phi, R, loss = .008, .004, .006, .02, .04, .4, .007
mdot = phi - (1 - R) * psi - loss
metaldot = yield_p * psi + Zin * phi - Z * ((1 - R) * psi + loss)
check("metal_budget_elimination", close(metaldot - Z * mdot, yield_p * psi + (Zin - Z) * phi), "Eq2")
zrate = (metaldot - Z * mdot) / 20
check("metal_inversion", close((yield_p * psi - 20 * zrate) / (Z - Zin), phi), "Eq46")
delta = 1e-6
value = yield_p * psi / phi
numeric_derivative = (value * 10**delta - value * 10**(-delta)) / (2 * delta)
check("base10_differential", close(numeric_derivative, math.log(10) * value, rtol=1e-8), numeric_derivative)

fhe = MU / 1.4
F = .14 * fhe * .5 * 2**2.3
D = .30 * fhe * .2 * 2 / (-math.expm1(-3.8))
mean_local = .14 * fhe * .5 * (.5**2.3 + 3.5**2.3) / 2
clumping = mean_local / F
check("chemical_worked_example", close([F, D, F-D], [.33487128043, .11923889468, .21563238575]), [F, D, F-D])
check("beam_nonlinearity", close(clumping, 1.831780222936047), {"local_mean": mean_local, "mean_column_rate": F, "ratio": clumping})

physical_F_pref = (2 * 1.4 * MH * 2e-16 / u.s * 1e21 / u.cm**2).to_value(RATE)
check("formation_coefficient_rounding", close(physical_F_pref, .14, rtol=.025), physical_F_pref)
p = .15
IR_pref = (2 * 1.4 * MH * 4 * math.pi * p / ((1-p) * 3.6) * 1e5 / u.s / u.cm**2).to_value(RATE)
FUV_pref = (2 * 1.4 * MH * 4 * math.pi * p / (1-p) * 1e5 / u.s / u.cm**2 * 1.9).to_value(RATE)
check("IR_coefficient_rounding", close(IR_pref, .044, rtol=.025), IR_pref)
check("FUV_coefficient_rounding", close(FUV_pref, .30, rtol=.025), FUV_pref)

sigma_col = (MU * MH * 1e21 / u.cm**2).to_value(u.Msun / u.pc**2)
check("column_to_mass", close(sigma_col, 10.9, rtol=.005), sigma_col)
hi_coeff = (MU * MH * 1.823e18 / u.cm**2).to_value(u.Msun / u.pc**2)
check("HI_surface_density_coefficient", close(hi_coeff, .0199, rtol=.005), hi_coeff)

pc_arcsec = ((16.5 * u.Mpc) * (1 * u.arcsec).to_value(u.rad)).to_value(u.pc)
speed = (1 * u.km / u.s).to_value(u.pc / u.Myr)
ring_flux = (2 * math.pi * 50 * u.pc * 50 * u.Msun / u.pc**2 * 2 * u.km / u.s).to_value(u.Msun / u.yr)
check("angular_scale", close(pc_arcsec, 80.0, rtol=.001), pc_arcsec)
check("velocity_units", close(speed, 1.022712165, rtol=1e-8), speed)
check("ring_flux_example", close(ring_flux, .03212945024, rtol=1e-8), ring_flux)
rin, rout, k = 2.0, 3.0, .4
annulus = (k*rout**2-k*rin**2)/(math.pi*(rout**2-rin**2))
check("annulus_convergence_sign", close(annulus, k/math.pi) and annulus > 0, "Inward mass flux increasing outward supplies an annulus")

Tb = 1.222e6 * .000219 / (1.420405752**2 * 8.2 * 7.1)
nhi_3sigma = 3 * 1.823e18 * Tb * math.sqrt(16 * 1.4)
check("MHONGOOSE_sensitivity_reconstruction", close(nhi_3sigma, 10**19.77, rtol=.03), {"Tb_rms_K": Tb, "NHI_3sigma": nhi_3sigma})
check("beam_sensitivity_scaling", (8/1)**2 == 64 and ((8/1)**2)**2 == 4096, "Equal flux noise scaling, not an exposure forecast")

models = [(Fraction(3,100), Fraction(1,100)), (Fraction(3,10), Fraction(28,100))]
for ff, dd in models:
    assert ff-dd-Fraction(8,1000)+Fraction(8,1000)-Fraction(2,100) == 0
check("identical_snapshot_examples", True, "Both illustrative molecular budgets are stationary")
check("lifecycle_proxy", close(.7*20/20, .7), "0.7 Msun pc^-2 Myr^-1; 35 times the toy SFR")
check("net_rate_error", close(math.hypot(.03,.03), .04242640687119285), "Illustrative uncorrelated errors")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.labelcolor": "#273e47", "text.color": "#273e47",
                     "axes.titleweight": "bold", "axes.titlesize": 12})
teal, orange, dark = "#216c78", "#bc683e", "#273e47"

fig, ax = plt.subplots(figsize=(10.6, 4.3))
ax.set(xlim=(0,1), ylim=(0,1)); ax.axis("off")
def box(x,y,w,h,label):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.01,rounding_size=0.014",
                 fc="#eef5f5",ec="#8eaeb3",lw=1.2))
    ax.text(x+w/2,y+h/2,label,ha="center",va="center",fontsize=11)
def arrow(start,end,label,tx,ty,color=teal):
    ax.annotate("", xy=end, xytext=start,
                arrowprops=dict(arrowstyle="-|>",lw=1.8,color=color,mutation_scale=14))
    ax.text(tx,ty,label,ha="center",va="center",color=color,fontsize=11)
box(.015,.48,.235,.24,"Non-molecular ISM\ninside the aperture")
box(.38,.48,.24,.24,"Molecular reservoir\nCO-bright + CO-dark")
box(.75,.48,.235,.24,"Molecular gas\noutside the aperture")
box(.395,.085,.21,.18,"Stars")
arrow((.26,.655),(.37,.655),"Formation  F",.315,.785)
arrow((.37,.53),(.26,.53),"Destruction  D",.315,.415,orange)
arrow((.74,.655),(.63,.655),r"Inward  $B_m^+$",.685,.785)
arrow((.63,.53),(.74,.53),r"Outward  $B_m^-$",.685,.415,orange)
arrow((.55,.47),(.55,.275),r"Star formation  $\psi$",.72,.34,orange)
arrow((.45,.275),(.45,.47),r"Direct return  $E_m$",.265,.34)
ax.text(.5,.96,"Chemical conversion and molecular transport are separate processes",
        ha="center",va="center",fontsize=13,weight="bold")
ax.text(.5,.01,"Schematic bookkeeping only; arrow widths do not encode rates.",
        ha="center",va="bottom",fontsize=9,color="#596c72")
fig.savefig(OUT / "reservoir_budget.png",dpi=260,bbox_inches="tight",facecolor="white")
plt.close(fig)

fig, axes = plt.subplots(1,2,figsize=(10.6,4.35),gridspec_kw={"wspace":.38})
for tau,color in [(4,teal),(40,orange)]:
    x = xeq + (x0-xeq)*np.exp(-t/tau)
    axes[0].plot(t,x,color=color,lw=2.4,label=f"Relaxation time {tau} Myr")
axes[0].axhline(xeq,color="#82959a",ls="--",lw=1)
axes[0].set(xlabel="Time (Myr)",ylabel="Molecular fraction",xlim=(0,200),ylim=(0,1),
            title="A. Same equilibrium, different rates")
axes[0].legend(loc="lower right",frameon=False,fontsize=9)
axes[0].text(86,.83,"Same equilibrium fraction: 0.8",fontsize=9,color="#596c72")
axes[1].bar([0,1],[mean_local,F],color=[teal,orange],width=.58)
axes[1].set_xticks([0,1],["Average of\nlocal rates","Rate from\naveraged columns"],fontsize=10)
axes[1].set(ylabel=r"Formation rate ($M_\odot$ pc$^{-2}$ Myr$^{-1}$)",ylim=(0,.8),
            title="B. A nonlinear beam-averaging bias")
for j,val in enumerate([mean_local,F]):
    axes[1].text(j,val+.02,f"{val:.3f}",ha="center",fontsize=11,weight="bold")
axes[1].text(.5,.74,f"Ratio = {clumping:.2f}",ha="center",fontsize=10,color="#596c72")
fig.subplots_adjust(bottom=.24,top=.87)
fig.text(.5,.025,"Constructed analytical examples; no MAUVE observations or fitted rates.",
         ha="center",fontsize=9,color="#596c72")
fig.savefig(OUT / "identifiability_and_beam.png",dpi=260,bbox_inches="tight",facecolor="white")
plt.close(fig)

result = {"status": "passed", "check_count": len(checks), "checks": checks,
          "method": "Exact rational identities, dimensional checks with astropy, numerical analytic-solution checks; no SymPy dependency",
          "versions": {"python": platform.python_version(), "numpy": np.__version__,
                       "astropy": astropy.__version__, "matplotlib": matplotlib.__version__},
          "examples": {"F": F,"D_photo": D,"net_chemical": F-D,"beam_local_mean":mean_local,
                       "beam_mean_column_rate":F,"clumping_ratio":clumping,
                       "mass_column_Msun_pc2":sigma_col,"formation_Myr":tform,
                       "rate_change_timescale_Myr":sigma_col/(F-D),"ring_Msun_yr":ring_flux},
          "figures": ["reservoir_budget.png","identifiability_and_beam.png"],
          "limits": "Does not validate physical closures or make any MAUVE rate measurement."}
(OUT / "verification_math.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps(result,indent=2))
