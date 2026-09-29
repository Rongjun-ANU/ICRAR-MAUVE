"""Conditional resolved-gas and line predictions for the 29 September report.

Reads the 14 September scalar export only after checking that its source
notebook and SFR pipeline still match their recorded fingerprints. It does not
write to science notebooks or FITS products.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
OLD = OUT.parent / "20260914_resolved_RPS_academic_model"
FP = json.loads((OLD / "input_fingerprints.json").read_text())
for record in FP:
    path = Path(record["path"])
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != record["sha256"]:
        raise RuntimeError(f"Source changed since scalar extraction: {path}")
(OUT / "source_fingerprints.json").write_text(json.dumps(FP, indent=2) + "\n")

data = pd.read_csv(OLD / "stage_bpt_line_profiles_with_hbeta.csv")
bin_center = 8.625
data = data.loc[(data.variable == "log_sigma_star") & (data.center == bin_center)]
anchors = []
for (stage, category), group in data.groupby(["stage", "category"]):
    by_line = group.set_index("line")
    val = by_line.mean_line_surface.to_dict()
    counts = by_line.N_gal_support.to_dict()
    if all(np.isfinite(val.get(line, np.nan)) for line in
           ["HA6562", "HB4861", "NII6583", "SII_SUM", "OIII5006"]):
        anchors.append(dict(stage=stage, category=category,
                            n_gal_common=int(min(counts.values())),
                            ha_surface=val["HA6562"],
                            hb_surface=val["HB4861"],
                            nii_surface=val["NII6583"],
                            sii_surface=val["SII_SUM"],
                            oiii_surface=val["OIII5006"],
                            nii_ha=val["NII6583"] / val["HA6562"],
                            sii_ha=val["SII_SUM"] / val["HA6562"],
                            oiii_hb=val["OIII5006"] / val["HB4861"],
                            corrected_balmer=val["HA6562"] / val["HB4861"]))
anchor_df = pd.DataFrame(anchors)
anchor_df.to_csv(OUT / "mauve_anchor_table.csv", index=False)

def anchor(stage, category):
    return anchor_df.set_index(["stage", "category"]).loc[(stage, category)]

# All times below are Gyr, all neutral-gas columns are M_sun pc^-2, and
# luminosities are erg s^-1 kpc^-2. C_alpha follows the unchanged MAUVE script.
C_ALPHA = 4.983582089552239e-42  # M_sun yr^-1 / (erg s^-1)
TAU_DEP = 2.0                  # Gyr
RETURN_FRACTION = 0.4         # effective prompt return to modeled cold gas
ETA = 0.0                     # feedback mass loading in this isolated example
HI0 = 10.0                    # M_sun pc^-2: declared example, not MAUVE H I
GAMMA_STRIP_HI = 3.0          # Gyr^-1: effective interval coefficient
GAMMA_STRIP_H2 = 2.22         # Gyr^-1: effective interval coefficient
TAU_ION = 0.003               # Gyr = 3 Myr: illustrative response kernel
F_HII0 = 0.855
F_UNABS0 = 0.05
F_UNABS_MAX = 0.20
B_CASEB = 2.86
T_END = 1.0                   # Gyr, a demonstration coordinate, not a stage age

L_OBS_PRE_SF = float(anchor("pre-peak", "SF").ha_surface)
L_FULL_YOUNG0 = L_OBS_PRE_SF / (1.0 - F_UNABS0)
PSI0 = C_ALPHA * L_FULL_YOUNG0             # M_sun yr^-1 kpc^-2
H20 = PSI0 * TAU_DEP * 1000.0               # M_sun pc^-2
TAU_CONV = HI0 / (((1 - RETURN_FRACTION + ETA) / TAU_DEP) * H20)

GAMMA_LOSS_HI = 1 / TAU_CONV + GAMMA_STRIP_HI
GAMMA_LOSS_H2 = (1 - RETURN_FRACTION + ETA) / TAU_DEP + GAMMA_STRIP_H2

def hi_exact(t):
    return HI0 * np.exp(-GAMMA_LOSS_HI * np.asarray(t))

def h2_exact(t):
    t = np.asarray(t)
    if abs(GAMMA_LOSS_H2 - GAMMA_LOSS_HI) < 1e-11:
        return np.exp(-GAMMA_LOSS_H2 * t) * (H20 + HI0 * t / TAU_CONV)
    return (H20 * np.exp(-GAMMA_LOSS_H2 * t)
            + HI0 / TAU_CONV
            * (np.exp(-GAMMA_LOSS_HI * t) - np.exp(-GAMMA_LOSS_H2 * t))
            / (GAMMA_LOSS_H2 - GAMMA_LOSS_HI))

def psi_exact(t):
    return h2_exact(t) / (1000.0 * TAU_DEP)

t = np.linspace(0, T_END, 1001)
hi, h2, psi = hi_exact(t), h2_exact(t), psi_exact(t)

def rhs_constant(tt, state):
    a, m = state
    return [-GAMMA_LOSS_HI * a, a / TAU_CONV - GAMMA_LOSS_H2 * m]

sol = solve_ivp(rhs_constant, (0, T_END), [HI0, H20], t_eval=t,
                rtol=1e-11, atol=1e-12, max_step=0.002)
if not sol.success:
    raise RuntimeError(sol.message)
ode_rel = float(np.max(np.abs(sol.y - np.vstack([hi, h2])) /
                       np.maximum(np.vstack([hi, h2]), 1e-12)))

# K_alpha(a) = exp(-a/tau_ion)/tau_ion and pre-perturbation SFR = PSI0.
# This ODE is exactly the convolution for that specific illustrative kernel.
ion = solve_ivp(lambda tt, y: [(psi_exact(tt) - y[0]) / TAU_ION],
                (0, T_END), [PSI0], t_eval=t, rtol=1e-10, atol=1e-12,
                max_step=0.001)
if not ion.success:
    raise RuntimeError(ion.message)
psi_recent = ion.y[0]

# Conditional photon-partition closure. The young-photon fractions are mutually
# exclusive, and f_leak represents photons absorbed in diffuse, non-HII gas.
f_hii = F_HII0 * h2 / H20
f_unabs = F_UNABS0 + (F_UNABS_MAX - F_UNABS0) * (1 - hi / HI0)
f_leak = 1 - f_hii - f_unabs
if np.min(f_leak) < -1e-12 or np.max(f_hii + f_leak + f_unabs - 1) > 1e-12:
    raise RuntimeError("Young-photon partition is invalid")

l_full = psi_recent / C_ALPHA
l_hii = f_hii * l_full
l_non = f_leak * l_full
l_ha = l_hii + l_non
l_hb_hii, l_hb_non = l_hii / B_CASEB, l_non / B_CASEB
l_hb = l_hb_hii + l_hb_non
w_non = l_non / l_ha

# Rounded *envelopes* of the measured SF/NSF line-ratio range, not identified
# pure source spectra. Their dependence on observed endpoints is stated in the
# report: agreement at those endpoints is not an independent validation.
templates = {
    "NII_Ha": {"HII": 0.22, "non_HII": 0.60},
    "SII_Ha": {"HII": 0.20, "non_HII": 0.45},
    "OIII_Hb": {"HII": 0.57, "non_HII": 0.77},
}
l_nii = templates["NII_Ha"]["HII"] * l_hii + templates["NII_Ha"]["non_HII"] * l_non
l_sii = templates["SII_Ha"]["HII"] * l_hii + templates["SII_Ha"]["non_HII"] * l_non
l_oiii = templates["OIII_Hb"]["HII"] * l_hb_hii + templates["OIII_Hb"]["non_HII"] * l_hb_non
r_nii, r_sii, r_oiii = l_nii / l_ha, l_sii / l_ha, l_oiii / l_hb

# The luminosity ceiling from old-star photoionization assumes every such
# photon is absorbed. It is a bound within the cited HOLMES model, not a
# universal bound on shock-, AGN-, or young-leakage-powered non-HII emission.
MSTAR_SURFACE = 10**bin_center       # M_sun kpc^-2
Q_OLD_PER_MASS = 7e40              # photons s^-1 M_sun^-1 (Belfiore+ 2022)
H_NU_HA = 3.027e-12                # erg photon^-1
P_HA = 0.45                        # illustrative Case-B Halpha yield
L_OLD_MAX = H_NU_HA * P_HA * Q_OLD_PER_MASS * MSTAR_SURFACE
ALPHA_B = 2.6e-13                  # cm^3 s^-1, illustrative 1e4 K Case B
CM_PER_PC = 3.085677581491367e18
CM_PER_KPC = 1000 * CM_PER_PC
M_H = 1.6735575e-24              # g
MSUN_G = 1.98847e33
G_PER_MSUN_PC2 = MSUN_G / CM_PER_PC**2
EM_NON_END_PC_CM6 = (l_non[-1] / (H_NU_HA * P_HA * ALPHA_B)
                     / CM_PER_KPC**2 / CM_PER_PC)

# A local, short efficiency boost is a separate conditional experiment.
# It changes the molecular depletion time and is not labelled a measured
# pressure-to-efficiency law.
PULSE_AMP = 0.25
PULSE_PEAK = 0.060                 # Gyr
PULSE_SHAPE = 9                    # positive dimensionless shape exponent

def pulse(tval):
    x = np.asarray(tval) / PULSE_PEAK
    return PULSE_AMP * x**PULSE_SHAPE * np.exp(PULSE_SHAPE * (1 - x))

def tau_dep_pulse(tval):
    return TAU_DEP / (1 + pulse(tval))

def rhs_pulse(tt, y):
    a, m = y
    tau = tau_dep_pulse(tt)
    return [-GAMMA_LOSS_HI * a,
            a / TAU_CONV - ((1 - RETURN_FRACTION + ETA) / tau + GAMMA_STRIP_H2) * m]

p_sol = solve_ivp(rhs_pulse, (0, T_END), [HI0, H20], t_eval=t,
                  rtol=1e-10, atol=1e-11, max_step=0.0005)
if not p_sol.success:
    raise RuntimeError(p_sol.message)
p_h2 = p_sol.y[1]
p_psi = p_h2 / (1000.0 * tau_dep_pulse(t))

# Exact local logarithmic decline rate from molecular balance and the
# derivative of tau_dep(t). Negative Gamma means the local SFR is rising.
gamma_sfr = (GAMMA_LOSS_H2 - hi / (TAU_CONV * h2))
g = pulse(t)
gprime = (np.divide(PULSE_SHAPE * g, t, out=np.zeros_like(t), where=t > 0)
          - PULSE_SHAPE * g / PULSE_PEAK)
dln_tau_dt = -gprime / (1 + g)
p_gamma_sfr = ((1 - RETURN_FRACTION + ETA) / tau_dep_pulse(t)
               + GAMMA_STRIP_H2 - hi / (TAU_CONV * p_h2) + dln_tau_dt)
exposure = np.log(PSI0 / psi)
p_exposure = np.log(PSI0 / p_psi)

pred = pd.DataFrame(dict(time_Gyr=t, Sigma_HI=hi, Sigma_H2=h2,
                         Sigma_Phi=hi / TAU_CONV, Sigma_SFR=psi,
                         recent_SFR_equivalent=psi_recent,
                         tau_dep_Gyr=np.full(len(t), TAU_DEP),
                         gamma_loss_HI=GAMMA_LOSS_HI,
                         gamma_loss_H2=GAMMA_LOSS_H2,
                         Gamma_SFR=gamma_sfr, exposure=exposure,
                         Sigma_H2_pulse=p_h2, Sigma_SFR_pulse=p_psi,
                         tau_dep_pulse_Gyr=tau_dep_pulse(t),
                         Gamma_SFR_pulse=p_gamma_sfr,
                         exposure_pulse=p_exposure,
                         f_HII=f_hii, f_leak=f_leak, f_unabs=f_unabs,
                         L_Ha_HII=l_hii, L_Ha_non_HII=l_non,
                         L_Ha_total=l_ha, L_Hb_total=l_hb,
                         L_NII=l_nii, L_SII=l_sii, L_OIII=l_oiii,
                         w_non_HII=w_non, NII_Ha=r_nii,
                         SII_Ha=r_sii, OIII_Hb=r_oiii))
pred.to_csv(OUT / "numerical_predictions.csv", index=False)

peak = int(np.argmax(p_psi))
rising = np.where(p_gamma_sfr < 0)[0]
check = dict(
    analytic_vs_ODE_max_relative_error=ode_rel,
    gas_partition_max_error=float(np.max(np.abs(f_hii + f_leak + f_unabs - 1))),
    photon_fraction_minimum=float(np.min(np.vstack([f_hii, f_leak, f_unabs]))),
    photon_fraction_maximum=float(np.max(np.vstack([f_hii, f_leak, f_unabs]))),
    lha_additivity_max_relative_error=float(np.max(np.abs(l_hii + l_non - l_ha) / l_ha)),
    first_rising_time_Myr=float(t[rising[0]] * 1000) if len(rising) else None,
    pulse_peak_time_Myr=float(t[peak] * 1000),
    pulse_peak_relative_to_control=float(p_psi[peak] / PSI0),
    pulse_peak_relative_to_no_pulse=float(p_psi[peak] / psi[peak]),
    pulse_peak_Gamma_SFR_per_Gyr=float(p_gamma_sfr[peak]),
    pulse_peak_exposure=float(p_exposure[peak]),
    end_H2_ratio=float(h2[-1] / H20),
    end_SFR_ratio=float(psi[-1] / PSI0),
    end_Ha_ratio=float(l_ha[-1] / l_ha[0]),
    end_HII_Ha_ratio=float(l_hii[-1] / l_hii[0]),
    end_non_HII_Ha_ratio=float(l_non[-1] / l_non[0]),
    end_non_HII_weight=float(w_non[-1]),
    end_instantaneous_vs_kernel_relative=float((psi_recent[-1] - psi[-1]) / psi[-1]),
    old_star_caseB_Lha_ceiling=L_OLD_MAX,
    old_ceiling_over_post_NSF_Ha=float(L_OLD_MAX / anchor("post-peak", "NSF").ha_surface),
    final_non_HII_emission_measure_pc_cm6=float(EM_NON_END_PC_CM6),
    uniform_ionized_column_200pc_Msun_pc2=float(
        np.sqrt(EM_NON_END_PC_CM6 / 200) * 200 * CM_PER_PC * M_H / G_PER_MSUN_PC2),
    uniform_ionized_column_1000pc_Msun_pc2=float(
        np.sqrt(EM_NON_END_PC_CM6 / 1000) * 1000 * CM_PER_PC * M_H / G_PER_MSUN_PC2),
)
assert check["analytic_vs_ODE_max_relative_error"] < 1e-7
assert check["gas_partition_max_error"] < 1e-12
assert 0 <= check["photon_fraction_minimum"]
assert check["photon_fraction_maximum"] <= 1
assert check["lha_additivity_max_relative_error"] < 1e-12
assert check["end_non_HII_Ha_ratio"] < 1 and check["end_HII_Ha_ratio"] < 1
assert check["pulse_peak_relative_to_control"] > 1
(OUT / "numerical_checks.json").write_text(json.dumps(check, indent=2) + "\n")

parameters = dict(
    data_anchor=dict(log10_Sigma_star_center=bin_center,
                     source_csv=str(OLD / "stage_bpt_line_profiles_with_hbeta.csv"),
                     exact_source_hashes_match=True),
    gas=dict(Sigma_HI0_Msun_pc2=HI0, Sigma_H20_Msun_pc2=H20,
             Sigma_SFR0_Msun_yr_kpc2=PSI0,
             tau_conv_Gyr=TAU_CONV, tau_dep_Gyr=TAU_DEP,
             gamma_strip_HI_per_Gyr=GAMMA_STRIP_HI,
             gamma_strip_H2_per_Gyr=GAMMA_STRIP_H2,
             gamma_loss_HI_per_Gyr=GAMMA_LOSS_HI,
             gamma_loss_H2_per_Gyr=GAMMA_LOSS_H2,
             return_fraction=RETURN_FRACTION, feedback_loading=ETA),
    line_model=dict(C_Halpha=C_ALPHA, tau_ion_Myr=TAU_ION * 1000,
                    caseB_balmer=B_CASEB, f_HII_initial=F_HII0,
                    f_unabs_initial=F_UNABS0, f_unabs_limit=F_UNABS_MAX,
                    q_old_per_msun=Q_OLD_PER_MASS,
                    hnu_Halpha_erg=H_NU_HA, p_Halpha=P_HA,
                    alpha_B_caseB_cm3_s=ALPHA_B,
                    templates=templates),
    pulse=dict(peak_efficiency_amplitude=PULSE_AMP,
               center_Myr=PULSE_PEAK * 1000,
               shape_exponent=PULSE_SHAPE),
    illustrative_duration_Gyr=T_END,
    fitted_physical_parameters=False,
)
(OUT / "model_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")

plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 200})
tm = t * 1000
blue, orange, teal = "#236c94", "#b76533", "#168879"

fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.3), layout="constrained")
axes[0, 0].plot(tm, hi, color=orange, label=r"$\Sigma_{\rm HI}$")
axes[0, 0].plot(tm, h2, color=blue, label=r"$\Sigma_{\rm H_2}$")
axes[0, 0].set(ylabel=r"Gas surface density ($M_\odot\,\mathrm{pc}^{-2}$)",
               title="Exact two-reservoir response")
axes[0, 0].legend()
axes[0, 1].plot(tm, psi / PSI0, color=blue, label="Gas model")
axes[0, 1].plot(tm, p_psi / PSI0, color=teal, label="Local efficiency pulse")
axes[0, 1].axhline(1, color=".4", lw=.8, ls=":")
axes[0, 1].set(ylabel="SFR / pre-perturbation SFR", title="A short local response")
axes[0, 1].legend()
axes[1, 0].plot(tm, p_gamma_sfr, color=orange)
axes[1, 0].axhline(0, color=".4", lw=.8, ls=":")
axes[1, 0].set(xlabel="Illustrative elapsed time (Myr)",
               ylabel=r"$\Gamma_{\rm SFR}$ (Gyr$^{-1}$)",
               title="Negative = instantaneously rising SFR")
axes[1, 1].plot(tm, p_exposure, color=teal)
axes[1, 1].axhline(0, color=".4", lw=.8, ls=":")
axes[1, 1].set(xlabel="Illustrative elapsed time (Myr)",
               ylabel=r"$\mathcal{E}=\ln[\psi_0/\psi(t)]$",
               title="Negative = above initial SFR")
fig.savefig(OUT / "figure_01_reservoir_and_pulse.png")
fig.savefig(OUT / "figure_01_reservoir_and_pulse.pdf")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), layout="constrained")
for key, color, label in [(f_hii, blue, "Compact H II absorption"),
                          (f_leak, orange, "Diffuse non-H II absorption"),
                          (f_unabs, ".4", "Unabsorbed/modelled elsewhere")]:
    axes[0].plot(tm, key, color=color, label=label)
axes[0].set(xlabel="Illustrative elapsed time (Myr)", ylabel="Young-photon fraction",
            title="Mutually exclusive young-photon pathways", ylim=(0, 1))
axes[0].legend(fontsize=8)
for arr, color, label in [(l_hii, blue, "H II"), (l_non, orange, "Non-H II"),
                           (l_ha, "black", "Total H-alpha")]:
    axes[1].plot(tm, arr, color=color, label=label)
axes[1].axhline(L_OBS_PRE_SF, color=".55", ls="--", lw=.8,
                label="Pre-peak SF reference")
axes[1].axhline(anchor("post-peak", "NSF").ha_surface,
                color="#915d83", ls=":", lw=1,
                label="Post-peak NSF reference")
axes[1].set(xlabel="Illustrative elapsed time (Myr)",
            ylabel=r"$\mathcal{L}_{\mathrm{H}\alpha}$ (erg s$^{-1}$ kpc$^{-2}$)",
            yscale="log", title="Diffuse emission can peak before fading")
axes[1].legend(fontsize=7)
fig.savefig(OUT / "figure_02_photon_partition_and_fading.png")
fig.savefig(OUT / "figure_02_photon_partition_and_fading.pdf")
plt.close(fig)

fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
line_data = [(l_ha, "H-alpha", "#182c36"),
             (l_nii, "[N II]", "#b04734"), (l_sii, "[S II] doublet", "#b6802d"),
             (l_oiii, "[O III]", "#346eac")]
for arr, label, color in line_data:
    axes[0, 0].plot(tm, arr, label=label, color=color)
axes[0, 0].set(yscale="log", ylabel=r"Luminosity surface density (erg s$^{-1}$ kpc$^{-2}$)",
               title="Absolute line luminosities")
axes[0, 0].legend(fontsize=8)
for k, (ratio, name, line) in enumerate([(r_nii, "[N II] / H-alpha", "NII_Ha"),
                                          (r_sii, "[S II] / H-alpha", "SII_Ha"),
                                          (r_oiii, "[O III] / H-beta", "OIII_Hb")]):
    ax = axes.ravel()[k + 1]
    ax.plot(tm, ratio, color=["#b04734", "#b6802d", "#346eac"][k])
    ax.axhline(anchor("pre-peak", "SF")[line.lower()],
               color=".45", ls="--", lw=.8, label="Pre-peak SF")
    ax.axhline(anchor("post-peak", "NSF")[line.lower()],
               color="#915d83", ls=":", lw=.8, label="Post-peak NSF")
    ax.set(title=name, ylabel="Linear ratio")
    if k == 0:
        ax.legend(fontsize=7)
axes[1, 1].plot(tm, w_non, color=teal)
axes[1, 1].set(title="Non-H II fraction of H-alpha", ylabel=r"$w_\alpha$")
axes[1, 2].plot(np.log10(r_nii), np.log10(r_oiii), color="#6b4b91")
axes[1, 2].scatter(np.log10(r_nii[[0, -1]]), np.log10(r_oiii[[0, -1]]),
                   c=[blue, orange], s=35, zorder=3)
axes[1, 2].annotate("start", (np.log10(r_nii[0]), np.log10(r_oiii[0])),
                    xytext=(6, 4), textcoords="offset points", fontsize=8)
axes[1, 2].annotate("end", (np.log10(r_nii[-1]), np.log10(r_oiii[-1])),
                    xytext=(6, 4), textcoords="offset points", fontsize=8)
axes[1, 2].set(title="Conditional BPT trajectory",
               xlabel="log10([N II] / H-alpha)",
               ylabel="log10([O III] / H-beta)")
for ax in axes[1, :2]:
    ax.set_xlabel("Illustrative elapsed time (Myr)")
axes[0, 1].set_xlabel("Illustrative elapsed time (Myr)")
axes[0, 2].set_xlabel("Illustrative elapsed time (Myr)")
fig.savefig(OUT / "figure_03_lines_and_ratios.png")
fig.savefig(OUT / "figure_03_lines_and_ratios.pdf")
plt.close(fig)

print(json.dumps(check, indent=2))
print("MAUVE_ANCHORED_TWO_RESERVOIR_AND_LINES_PASS")
