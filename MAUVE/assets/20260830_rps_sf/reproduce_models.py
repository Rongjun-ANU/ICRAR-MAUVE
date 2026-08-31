#!/usr/bin/env python3
"""Illustrative RPS/SF models and independent numerical checks.

No MAUVE measurements are loaded and no parameters are fitted to observations.
Run with a Python environment containing numpy, scipy, and matplotlib.
All outputs are written beside this script. Units are documented in the report.
"""

from pathlib import Path
import csv
import json
import os
import platform
import tempfile

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rps_sf_20260830_mpl"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D


OUT = Path(__file__).resolve().parent
G = 6.67430e-8
KB = 1.380649e-16
MP = 1.67262192369e-24
MSUN = 1.98847e33
PC = 3.085677581491367e18
MYR = 365.25 * 86400 * 1e6
SURFACE = MSUN / PC**2
KMS_KPC_MYR = 1e5 / (1e3 * PC) * MYR

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 8.7, "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 150,
    "savefig.dpi": 190, "axes.axisbelow": True,
    "pdf.fonttype": 42,
})
COLORS = ["#176B91", "#BA4D32", "#72518B", "#41835B"]


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white", dpi=300)
    plt.close(fig)


def table(stem, columns):
    with (OUT / f"{stem}.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))


def weight(sigma_g=20.0, c_z=8.0, sigma_star=100.0, c_star=30.0):
    """Approximate two-component vertical weight, in dyn cm^-2."""
    return np.pi * G / 2 * (sigma_g * SURFACE) * (
        sigma_g * SURFACE + c_z / c_star * sigma_star * SURFACE)


checks = {}

# The isothermal Lane-Emden sphere provides an independent check of the
# Bonnor-Ebert coefficient; its turbulent use remains an approximation.
x0 = 1e-6
lane = solve_ivp(
    lambda x, y: [y[1], np.exp(-y[0]) - 2*y[1]/x],
    [x0, 20], [x0*x0/6 - x0**4/120, x0/3 - x0**3/30],
    rtol=2e-11, atol=2e-13, dense_output=True, max_step=0.025)
xcrit = brentq(lambda x: 2*np.exp(-lane.sol(x)[0]) - lane.sol(x)[1]**2, 5, 8)
psi, dpsi = lane.sol(xcrit)
be_coefficient = xcrit*xcrit*dpsi*np.exp(-psi/2)/np.sqrt(4*np.pi)
checks["bonnor_ebert"] = {
    "xi_critical": float(xcrit), "mass_coefficient": float(be_coefficient),
    "density_contrast": float(np.exp(psi)),
    "relative_error_vs_1p182": float(abs(be_coefficient/1.182-1)),
}
assert abs(be_coefficient/1.182 - 1) < 5e-4
assert abs(xcrit - 6.451) < 0.002

# A direct minimization of the razor-thin dispersion relation checks Q=1.
q_errors = []
for q in [0.55, 0.8, 1.0, 1.2, 2.0]:
    fit = minimize_scalar(lambda x: 1+x*x-2*x/q, bounds=(0, 10), method="bounded")
    q_errors.append(abs(fit.fun - (1-1/q**2)))
checks["thin_disk_minimum_absolute_error"] = float(max(q_errors))
assert max(q_errors) < 1e-9

# A frozen vertical-profile thickness approximation has a different threshold
# from the free-boundary modes of a pressure-confined slab.
eta = np.linspace(0, 2, 250)
xc = np.array([brentq(lambda x: 1-x*x-2*h*x**3, 0, 1.00001) for h in eta])
qcrit = 2*xc / ((1+xc**2)*(1+eta*xc))
thick_errors = []
for h, qc in zip(eta[::20], qcrit[::20]):
    fit = minimize_scalar(
        lambda x: 1+x*x-2*x/(qc*(1+h*x)), bounds=(0, 5), method="bounded")
    thick_errors.append(abs(fit.fun))
checks["thickness_threshold_absolute_residual"] = float(max(thick_errors))
assert max(thick_errors) < 1e-8
assert np.all(np.diff(qcrit) < 0)
table("thickness_threshold", {"eta_kappaH_over_cR": eta, "x_critical": xc, "Q_critical": qcrit})

# Figure 1: an algebraic diagnostic, not an empirical RPS calibration.
gas = np.linspace(0.07, 1.6, 300)
press = np.geomspace(0.2, 15, 350)
gg, pp = np.meshgrid(gas, press)
fig, axs = plt.subplots(1, 2, figsize=(8.3, 4.2), sharey=True, layout="constrained")
for ax, support in zip(axs, [1.0, 1.5]):
    enhancement = gg*np.sqrt(pp)/support
    im = ax.pcolormesh(gg, pp, np.log10(enhancement), shading="auto", cmap="RdBu_r", rasterized=True,
                       norm=TwoSlopeNorm(vmin=-1.1, vcenter=0, vmax=1.1))
    contour = ax.contour(gg, pp, enhancement, levels=[0.5, 1, 2], colors="#252525",
                        linewidths=[0.8, 1.7, 0.8])
    positions = [(0.8, (0.5*support/0.8)**2), (0.9, (support/0.9)**2),
                 (1.15, (2*support/1.15)**2)]
    ax.clabel(contour, inline=True, fmt={0.5: "E=0.5", 1: "E=1", 2: "E=2"},
              fontsize=9, manual=positions)
    ax.set(yscale="log", xlabel=r"Gas ratio $g=\Sigma_g/\Sigma_{g,0}$",
           title=rf"Vertical support ratio $s={support:g}$")
    ax.axvline(1, color="0.3", ls=":", lw=0.8)
axs[0].set_ylabel(r"Midplane pressure ratio $p=P_{\rm mid}/P_{{\rm mid},0}$")
fig.colorbar(im, ax=axs, label=r"$\log_{10}(\Sigma_{\rm SFR}/\Sigma_{{\rm SFR},0})$", shrink=0.92)
save(fig, "01_enhancement_regimes")

# Figure 2: the exact constant-coefficient regulator solution versus an ODE.
tau0 = 2000/0.6  # Myr: gas depletion time / (1-return+mass loading).
time = np.linspace(0, 1200, 2401)
step_cases = [
    ("Pressure only", 2.0, 0.0, 1.0),
    ("Pressure + slower stripping", 2.0, 1/500, 1.0),
    ("Weak boost + faster stripping", 1.1, 1/250, 1.0),
]
fig, axs = plt.subplots(2, 1, figsize=(7.5, 5.5), sharex=True, layout="constrained")
step_summary = []
step_columns = {"time_Myr": time}
for n, (label, b, k, supply) in enumerate(step_cases):
    K = k*tau0
    teq = tau0/(b+K)
    geq = supply/(b+K)
    g = geq+(1-geq)*np.exp(-time/teq)
    E = b*g
    numeric = solve_ivp(lambda t, y: [supply/tau0-(b/tau0+k)*y[0]],
                        [0, time[-1]], [1.0], t_eval=time, rtol=1e-11, atol=1e-13)
    error = np.max(np.abs(numeric.y[0]-g))
    assert error < 2e-10
    einf = b*geq
    crossing = float(teq*np.log((b-einf)/(1-einf))) if einf < 1 and b > 1 else None
    step_summary.append({"label": label, "b": b, "k_per_Myr": k, "supply_ratio": supply,
                         "K": K, "tau_eq_Myr": teq, "E_infinity": einf,
                         "t_cross_Myr": crossing, "ODE_max_absolute_error": float(error)})
    axs[0].plot(time, E, label=label, color=COLORS[n], lw=2)
    axs[1].plot(time, g, color=COLORS[n], lw=2)
    step_columns[f"case{n+1}_gas_ratio"] = g
    step_columns[f"case{n+1}_SFR_ratio"] = E
axs[0].axhline(1, color="0.3", ls="--", lw=1)
axs[0].set(ylabel=r"SFR ratio $\mathcal{E}$", ylim=(0, 2.15))
axs[0].legend(loc="upper right")
axs[1].set(xlabel="Time after the imposed change [Myr]", ylabel="Gas ratio g", ylim=(0, 1.05))
save(fig, "02_exact_regulator")
table("exact_regulator", step_columns)
checks["exact_regulator"] = step_summary

# Figure 3: a deliberately uncalibrated pressure pulse with gas-dependent W.
# The identical, constant supply in each zone isolates pressure and stripping.
t = np.linspace(-400, 650, 4201)
pbg = 1e4 * KB
pram_peak = 1.2e5 * KB
w0 = weight()
pmid0 = w0+pbg
a = 0.6
nu0 = 1/2000  # Myr^-1, SFR / total cold gas surface density.
tresp = 15.0
pulse_cases = [
    ("Windward disc", 0.6, 0.25, 1/350),
    ("Leeward disc", 0.12, 1.0, 1/160),
]
pulse_summary = []
pulse_columns = {"time_Myr": t, "ram_pulse_P_over_k": 1.2e5*np.exp(-t*t/(2*80**2))}
fig, axs = plt.subplots(3, 1, figsize=(7.5, 7.7), sharex=True, layout="constrained")
for n, (label, chi, u, kpeak) in enumerate(pulse_cases):
    def rhs(ti, y):
        pulse = np.exp(-ti*ti/(2*80**2))
        g, b, _, _, _, _ = y
        support = np.sqrt(1+u*pulse)
        p = (weight(20*g, 8*support)+pbg+chi*pram_peak*pulse)/pmid0
        target = np.sqrt(p)/support
        E = b*g
        return [a*nu0*(1-E)-kpeak*pulse*g, (target-b)/tresp,
                (E-y[2])/5, (E-y[3])/100, a*nu0*E, kpeak*pulse*g]

    sol = solve_ivp(rhs, [t[0], t[-1]], [1, 1, 1, 1, 0, 0],
                    t_eval=t, rtol=2e-10, atol=2e-12, max_step=0.8)
    assert sol.success
    g, b, ha, uv, consumed, stripped = sol.y
    E = g*b
    pulse = np.exp(-t*t/(2*80**2))
    s = np.sqrt(1+u*pulse)
    p = (weight(20*g, 8*s)+pbg+chi*pram_peak*pulse)/pmid0
    balance = g+consumed+stripped-1-a*nu0*(t-t[0])
    assert np.max(np.abs(balance)) < 2e-11
    assert np.min(g) > 0 and np.min(b) > 0
    idx = int(np.argmax(E))
    pulse_summary.append({
        "label": label, "chi": chi, "u_support_squared": u,
        "k_peak_per_Myr": kpeak, "pressure_sigma_Myr": 80, "response_time_Myr": tresp,
        "SFR_peak_ratio": float(E[idx]), "SFR_peak_time_Myr": float(t[idx]),
        "SFR_at_ram_peak": float(np.interp(0, t, E)),
        "gas_at_ram_peak": float(np.interp(0, t, g)),
        "SFR_at_200_Myr": float(np.interp(200, t, E)),
        "gas_at_200_Myr": float(np.interp(200, t, g)),
        "stripped_initial_gas_units_at_200_Myr": float(np.interp(200, t, stripped)),
        "max_mass_balance_residual": float(np.max(np.abs(balance))),
    })
    axs[0].plot(t, p, color=COLORS[n], lw=2, label=label)
    axs[1].plot(t, g, color=COLORS[n], lw=2)
    axs[2].plot(t, E, color=COLORS[n], lw=2, label=label+" intrinsic")
    axs[2].plot(t, uv, color=COLORS[n], ls=":", lw=1.5, label=label+" 100 Myr kernel")
    for key, vals in [("gas_ratio", g), ("efficiency_ratio", b), ("SFR_ratio", E),
                      ("midplane_pressure_ratio", p), ("support_ratio", s),
                      ("tracer_5Myr_ratio", ha), ("tracer_100Myr_ratio", uv),
                      ("locked_and_outflow_initial_gas_units", consumed),
                      ("stripped_initial_gas_units", stripped)]:
        pulse_columns[f"zone{n+1}_{key}"] = vals
for ax in axs:
    ax.axvline(0, color="0.5", ls="--", lw=0.9)
    ax.grid(alpha=0.17)
axs[0].legend()
axs[0].set(ylabel="Midplane pressure ratio p")
axs[1].set(ylabel="Gas ratio g", ylim=(0.1, 1.08))
axs[2].axhline(1, color="0.5", ls="--", lw=0.9)
axs[2].set(ylabel=r"SFR ratio $\mathcal{E}$", xlabel="Time relative to ram-pressure maximum [Myr]")
axs[2].legend(handles=[Line2D([0], [0], color="0.25", lw=2, label="Intrinsic SFR"),
                       Line2D([0], [0], color="0.25", lw=1.5, ls=":",
                              label="100 Myr exponential kernel")], loc="upper right")
save(fig, "03_pressure_pulse")
table("pressure_pulse", pulse_columns)
checks["pressure_pulse"] = pulse_summary

# Figure 4: pressure and turbulent support compete differently for cloud collapse.
pc_ratio = np.geomspace(0.4, 100, 400)
fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.9), layout="constrained")
for n, sc in enumerate([1.0, 1.2, 1.5]):
    axs[0].plot(pc_ratio, sc**4/np.sqrt(pc_ratio), color=COLORS[n], lw=2,
                label=f"Cloud support ratio = {sc:g}")
axs[0].axhline(1, color="0.4", ls="--", lw=0.9)
axs[0].set(xscale="log", yscale="log", xlabel="Cloud boundary-pressure ratio",
           ylabel=r"$M_{\rm BE}'/M_{\rm BE}$", title="Pressure-confined cloud threshold")
axs[0].legend(loc="lower left")
axs[1].plot(eta, qcrit, color=COLORS[0], lw=2)
axs[1].set(xlabel=r"$\eta_H=\kappa H/c_R$", ylabel=r"$Q_{\rm crit}$",
           title="Frozen-profile thickness approximation", ylim=(0.35, 1.03))
axs[1].annotate("Thinner layer: higher threshold", xy=(0.05, 0.95), xytext=(0.4, 0.88),
                arrowprops={"arrowstyle": "->", "color": "0.3"}, fontsize=9)
save(fig, "04_cloud_and_disk_thresholds")

# Figure 5: the free-boundary modal speed is not the turbulent linewidth.
zeta = np.geomspace(0.01, 100, 300)
qp_over_qg = 1/np.sqrt(1+zeta)
fig, ax = plt.subplots(figsize=(6.8, 3.4), layout="constrained")
ax.plot(zeta, qp_over_qg, color=COLORS[0], lw=2.3)
ax.axhline(1, color="0.4", ls="--", label="Classical Q at fixed surface density and linewidth")
ax.set(xscale="log", xlabel=r"$\zeta=P_{\rm ext}/(\pi G\Sigma_g^2)$",
       ylabel=r"$Q_P/Q_g=c_{\rm mode}/c_s$", ylim=(0, 1.08))
ax.text(0.02, 0.18, "Symmetric isothermal slab; free boundaries; long wavelengths.\n"
        "This curve is neither an RPS calibration nor an SFR law.", fontsize=9)
ax.legend(loc="upper right")
save(fig, "05_pressure_confined_mode")
table("free_boundary_mode", {"zeta_external_pressure": zeta, "Q_P_over_Q_g": qp_over_qg})

# Figure 6: gas rotation gives a measurable, conditional response-time lag.
response = np.linspace(0, 80, 300)
omega = (200/5)*KMS_KPC_MYR
lag = np.arctan(omega*response)
amplitude = 1/np.sqrt(1+(omega*response)**2)
phi = np.linspace(-np.pi, np.pi, 500)
tr = 15
lag15 = float(np.arctan(omega*tr))
amp15 = float(1/np.sqrt(1+(omega*tr)**2))
angular = amp15*np.cos(phi-lag15)
residual = tr*omega*(-amp15*np.sin(phi-lag15)) + angular - np.cos(phi)
assert max(abs(residual)) < 1e-14
fig, axs = plt.subplots(1, 2, figsize=(8.2, 3.8), layout="constrained")
axs[0].plot(response, np.rad2deg(lag), color=COLORS[0], lw=2)
axs[0].set(xlabel="Response time [Myr]", ylabel="Downstream phase lag [degrees]",
           title="Circular speed 200 km/s at 5 kpc")
axs[0].axvline(15, color="0.4", ls=":")
axs[1].plot(np.rad2deg(phi), np.cos(phi), color="0.4", ls="--", label="First pressure harmonic")
axs[1].plot(np.rad2deg(phi), angular, color=COLORS[0], lw=2, label="15 Myr response")
axs[1].set(xlabel="Azimuth from upstream direction [degrees]", ylabel="Perturbation / input amplitude",
           title="Rotation attenuates and shifts the response")
axs[1].legend(loc="lower center")
save(fig, "06_azimuthal_lag")
table("azimuthal_lag", {"response_time_Myr": response, "phase_lag_degrees": np.rad2deg(lag),
                       "amplitude_over_forcing": amplitude})
checks["angular_response"] = {"omega_per_Myr": omega, "lag_15Myr_degrees": float(np.rad2deg(lag15)),
                               "amplitude_15Myr": amp15,
                               "max_equation_residual": float(max(abs(residual)))}

# Worked dimensional examples: inputs are explicitly invented for illustration.
ne, mu_e, v, temp = 3e-4, 1.17, 1000e5, 2e7
rho_icm = mu_e*MP*ne
pram = rho_icm*v*v
pth = 1.95*ne*KB*temp
baseline_p = weight()+pth
perturbed_p = weight(20, 8*1.1)+pth+0.5*pram
p_example = perturbed_p/baseline_p
q_example = (40*KMS_KPC_MYR/MYR)*(8e5)/(np.pi*G*20*SURFACE)
be_mass = be_coefficient*(2e5)**4/(G**1.5*np.sqrt(1e5*KB))/MSUN
rho_mid = baseline_p/(8e5)**2
tff_mid = np.sqrt(3*np.pi/(32*G*rho_mid))/MYR
H_mid_pc = (20*SURFACE/(2*rho_mid))/PC
checks["worked_example"] = {
    "n_e_per_cm3": ne, "mu_e": mu_e, "v_kms": v/1e5, "T_K": temp,
    "P_ram_over_k_K_cm3": float(pram/KB), "P_thermal_over_k_K_cm3": float(pth/KB),
    "weight_baseline_over_k_K_cm3": float(weight()/KB),
    "weight_10pct_support_increase_over_k_K_cm3": float(weight(20, 8.8)/KB),
    "Pmid_baseline_over_k_K_cm3": float(baseline_p/KB),
    "Pmid_perturbed_over_k_K_cm3": float(perturbed_p/KB),
    "p_ratio_with_10pct_support_increase": float(p_example),
    "b_ratio_with_10pct_support_increase": float(np.sqrt(p_example)/1.1),
    "Q_g": float(q_example), "Q_g_after_g1p5_sR1p2": float(q_example*1.2/1.5),
    "M_BE_Msun_c2kms_Poverk1e5": float(be_mass),
    "pressure_ratio_needed_for_20pct_cloud_support_increase": 1.2**8,
    "baseline_midplane_tff_Myr": float(tff_mid), "baseline_H_eff_pc": float(H_mid_pc),
    "baseline_vertical_crossing_Myr": float(H_mid_pc*PC/(8e5)/MYR),
}

# Verify the pressure inverse over varying gas/support/phase/clumping factors.
rng = np.random.default_rng(27082026)
g, s, f, e, d, p = np.exp(rng.uniform(-1, 1, (6, 500)))
E = g*f*e*np.sqrt(d*p)/s
recovered = (E*s/(g*f*e))**2/d
inverse_error = np.max(np.abs(recovered/p-1))
assert inverse_error < 2e-15
checks["pressure_inverse_max_relative_error"] = float(inverse_error)
checks["runtime"] = {"python": platform.python_version(), "numpy": np.__version__,
                      "scipy": scipy.__version__, "matplotlib": matplotlib.__version__}
checks["scope"] = "Illustrative closures; algebra and implementation checked, not observational validation."
(OUT / "numerical_checks.json").write_text(json.dumps(checks, indent=2)+"\n")
print(json.dumps(checks, indent=2))
