"""Reproducible illustrations for the 2026-08-30 nebular metallicity report.

These are fixed-state emissivity calculations, NOT a photoionization model,
an abundance calibration, or an analysis of MAUVE observations.
"""

import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/nebular_metallicity_mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyneb as pn


OUT = Path(__file__).resolve().parent
O3 = pn.Atom("O", 3)
H1 = pn.RecAtom("H", 1)
NE = 100.0
ION_ABUNDANCE = 2.5e-4


def eps_o(temperature, wave):
    return np.asarray(O3.getEmissivity(tem=temperature, den=NE, wave=wave))


def eps_h(temperature):
    return np.asarray(H1.getEmissivity(tem=temperature, den=NE, wave=4861))


def solve_temperature(nebular_to_auroral):
    return float(O3.getTemDen(nebular_to_auroral, den=NE,
                 to_eval="(L(4959)+L(5007))/L(4363)"))


rows = []
for temperature in (5000., 6000., 8000., 10000., 12000., 15000., 20000.):
    neb = eps_o(temperature, 4959) + eps_o(temperature, 5007)
    aur = eps_o(temperature, 4363)
    strong = float(ION_ABUNDANCE * eps_o(temperature, 5007) / eps_h(temperature))
    rows.append({"Te_K": temperature, "OIII5007_over_Hbeta": strong,
                 "OIII4363_over_4959plus5007": float(aur / neb)})

at10k = rows[3]
temperature_recovered = solve_temperature(1. / at10k["OIII4363_over_4959plus5007"])
abundance_recovered = float(O3.getIonAbundance(
    int_ratio=at10k["OIII5007_over_Hbeta"], tem=temperature_recovered,
    den=NE, wave=5007, Hbeta=1.0))
assert abs(temperature_recovered / 10000. - 1.) < 0.003
assert abs(abundance_recovered / ION_ABUNDANCE - 1.) < 0.005

# Equal H+ emission measure, same density, same O++/H+ in the two components.
temps = np.array([8000., 12000.])
mix_neb = np.mean(eps_o(temps, 4959) + eps_o(temps, 5007))
mix_aur = np.mean(eps_o(temps, 4363))
mix_hb = np.mean(eps_h(temps))
mix_t = solve_temperature(float(mix_neb / mix_aur))
mix_5007_hb = float(ION_ABUNDANCE * np.mean(eps_o(temps, 5007)) / mix_hb)
mix_a = float(O3.getIonAbundance(mix_5007_hb, tem=mix_t, den=NE,
                               wave=5007, Hbeta=1.0))

critical = []
for element, ion, wave in (("S",2,6716),("S",2,6731),("N",2,6584),
                           ("N",2,5755),("O",3,5007),("O",3,4363),
                           ("S",3,9069),("S",3,6312)):
    atom = pn.Atom(element, ion)
    upper, lower = atom.getTransition(wave)
    critical.append({"ion": f"{element}{ion}", "wave_A_label":wave,
                     "upper_level":int(upper), "lower_level":int(lower),
                     "ncrit_cm-3_at_10000K":float(atom.getCritDensity(tem=10000,level=upper)),
                     "atom_file":atom.atomFile, "collision_file":atom.collFile})

lines = {"OII3727":3727.,"CII4267":4267.,"OIII4363":4363.,
         "OII_RL4650":4650.,"Hbeta":4861.,"OIII5007":5007.,
         "NII5755":5755.,"HeI5876":5876.,"SIII6312":6312.,
         "Halpha":6563.,"NII6584":6584.,"SII6731":6731.,
         "HeI7065":7065.,"ArIII7136":7136.,"OII7320":7320.,
         "OII7330":7330.,"SIII9069":9069.,"SIII9531":9531.}
coverage = {name:{str(upper):4800. <= wave <= upper for upper in (7000,8900,9300)}
            for name,wave in lines.items()}
assert coverage["OII7320"] == {"7000":False,"8900":True,"9300":True}
assert coverage["SIII9069"] == {"7000":False,"8900":False,"9300":True}
assert not any(coverage["OIII4363"].values())
assert not any(coverage["SIII9531"].values())

result = {"purpose":"Fixed-state emissivity and arithmetic illustration; no MAUVE data used",
          "pyneb_version":pn.__version__, "density_cm-3":NE,
          "O++_over_H+_input":ION_ABUNDANCE,
          "atomic_data":{"OIII_A":O3.atomFile,"OIII_collision":O3.collFile,
                         "HI_recombination":H1.recFitsFile},
          "table":rows,
          "isothermal_recovery":{"Te_K":temperature_recovered,
                                  "O++_over_H+":abundance_recovered},
          "equal_EM_8000_12000K":{"Te_EM_mean_K":10000.,"t2":0.04,
                                  "Te_inferred_K":mix_t,
                                  "OIII5007_over_Hbeta":mix_5007_hb,
                                  "O++_over_H+_inferred":mix_a,
                                  "bias_dex":math.log10(mix_a/ION_ABUNDANCE)},
          "critical_densities":critical,"MUSE_z0_coverage":coverage,
          "redshift_blue_entry":{"OIII4363":4800./4363.-1.,
                                  "OII3727":4800./3727.-1.},
          "N_enhancement_0p3dex_fixed_other_fluxes":{
              "D16_dOH_dex":1.264*0.3,"M13_N2_dOH_dex":0.462*0.3,
              "M13_O3N2_dOH_dex":0.214*0.3}}
S2 = pn.Atom("S", 2)
result["SII_density_ratio_at_10000K"] = [
    {"ne_cm-3":density,"I6716_over_I6731":float(
        S2.getEmissivity(tem=10000,den=density,wave=6716)/
        S2.getEmissivity(tem=10000,den=density,wave=6731))}
    for density in (10.,100.,1000.,10000.,100000.)]
(OUT / "physics_checks.json").write_text(json.dumps(result, indent=2)+"\n")

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                     "axes.spines.top":False,"axes.spines.right":False,
                     "svg.fonttype":"none"})
tgrid = np.linspace(5000.,20000.,180)
fig, axes = plt.subplots(1,2,figsize=(10.0,3.6),layout="constrained")
axes[0].plot(tgrid/1000, ION_ABUNDANCE*eps_o(tgrid,5007)/eps_h(tgrid),
             color="#19647e",lw=2.2)
axes[0].set(xlabel=r"Electron temperature ($10^3$ K)",
            ylabel=r"[O III] 5007 / H$\beta$",
            title="Strong-line response at fixed ionic abundance")
axes[1].plot(tgrid/1000,eps_o(tgrid,4363)/(eps_o(tgrid,4959)+eps_o(tgrid,5007)),
             color="#a54b32",lw=2.2)
axes[1].set(xlabel=r"Electron temperature ($10^3$ K)",
            ylabel="[O III] 4363 / (4959 + 5007)",
            title="Auroral ratio: the abundance cancels")
axes[1].set_yscale("log")
for ax in axes:
    ax.axvline(10,color="#7c8388",lw=1,ls="--")
    ax.grid(alpha=.17)
fig.suptitle(r"PyNeb illustration: O$^{++}$/H$^+$ = 2.5 $\times$ 10$^{-4}$; $n_e$ = 100 cm$^{-3}$",fontsize=11)
fig.savefig(OUT/"fixed_state_emissivities.svg",bbox_inches="tight")
plt.close(fig)

print(json.dumps(result,indent=2))
print("PASS: temperature/abundance round trip, wavelength limits, and algebra checks")
