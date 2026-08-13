#!/usr/bin/env python3
"""Reproducible two-driver gas-regulator experiment for the FMR review.

The calculation integrates

    dM_g/dt     = Phi - (1 - R + eta) epsilon M_g
    d(M_g Z)/dt = y epsilon M_g + Z_in Phi
                  - Z (1 - R + eta) epsilon M_g

for stochastic Ornstein-Uhlenbeck fluctuations in log(Phi) and log(epsilon).
Stellar mass is treated as fixed: this is a local/residual response experiment,
not a complete evolutionary population model.  The original inflow-only Forbes
(2014) limit is sigma_epsilon = 0; adding sigma_epsilon tests the efficiency
driver discussed by Wang & Lilly (2021) and Huang et al. (2026).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


R_RETURN = 0.4
ETA = 1.0
LAMBDA = 1.0 - R_RETURN + ETA
EPSILON_0 = 1.0
PHI_0 = 1.0
YIELD = 0.015
Z_IN = 0.002
TAU_EQ = 1.0 / (LAMBDA * EPSILON_0)


def _ou_step(x: np.ndarray, rho: float, rng: np.random.Generator) -> np.ndarray:
    """Exact one-step update for a unit-variance OU process."""
    return rho * x + np.sqrt(1.0 - rho * rho) * rng.standard_normal(x.shape)


def simulate_series(
    sigma_phi: float,
    sigma_epsilon: float,
    seed: int,
    *,
    n_steps: int = 32_000,
    burn_steps: int = 8_000,
    sample_stride: int = 8,
    dt: float | None = None,
    tau_driver: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return residual log10(SFR) and log10(Z) after burn-in."""
    dt = TAU_EQ / 80.0 if dt is None else dt
    tau_driver = 2.0 * TAU_EQ if tau_driver is None else tau_driver
    rho = float(np.exp(-dt / tau_driver))
    rng = np.random.default_rng(seed)

    # Start at the constant-driver equilibrium.
    mg = PHI_0 / (LAMBDA * EPSILON_0)
    z_eq = Z_IN + YIELD / LAMBDA
    mz = mg * z_eq
    x_phi = np.zeros(1)
    x_epsilon = np.zeros(1)
    log_sfr: list[float] = []
    log_z: list[float] = []

    for step in range(n_steps):
        x_phi = _ou_step(x_phi, rho, rng)
        x_epsilon = _ou_step(x_epsilon, rho, rng)
        # The -sigma^2/2 term keeps the arithmetic mean driver near its baseline.
        phi = PHI_0 * np.exp(sigma_phi * x_phi[0] - 0.5 * sigma_phi**2)
        epsilon = EPSILON_0 * np.exp(
            sigma_epsilon * x_epsilon[0] - 0.5 * sigma_epsilon**2
        )
        z = mz / mg
        dmg = phi - LAMBDA * epsilon * mg
        dmz = YIELD * epsilon * mg + Z_IN * phi - z * LAMBDA * epsilon * mg
        mg += dt * dmg
        mz += dt * dmz

        if step >= burn_steps and (step - burn_steps) % sample_stride == 0:
            log_sfr.append(float(np.log10(epsilon * mg)))
            log_z.append(float(np.log10(mz / mg)))

    sfr_residual = np.asarray(log_sfr) - np.mean(log_sfr)
    z_residual = np.asarray(log_z) - np.mean(log_z)
    return sfr_residual, z_residual


def correlation_and_slope(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Pearson correlation and OLS slope y(x)."""
    corr = float(np.corrcoef(x, y)[0, 1])
    slope = float(np.polyfit(x, y, 1)[0])
    return corr, slope


def scan_driver_ratio(
    ratios: np.ndarray,
    *,
    sigma_phi: float = 0.28,
    seeds: np.ndarray | None = None,
) -> list[dict[str, float]]:
    """Monte-Carlo scan of efficiency-to-inflow fluctuation amplitude."""
    if seeds is None:
        seeds = np.arange(12, dtype=int) + 7100
    rows: list[dict[str, float]] = []
    for ratio in ratios:
        correlations = []
        slopes = []
        for seed in seeds:
            x, y = simulate_series(
                sigma_phi=sigma_phi,
                sigma_epsilon=float(ratio * sigma_phi),
                seed=int(seed),
                n_steps=18_000,
                burn_steps=4_000,
                sample_stride=7,
            )
            corr, slope = correlation_and_slope(x, y)
            correlations.append(corr)
            slopes.append(slope)
        corr_values = np.asarray(correlations)
        slope_values = np.asarray(slopes)
        rows.append(
            {
                "sigma_epsilon_over_sigma_phi": float(ratio),
                "corr_median": float(np.median(corr_values)),
                "corr_p16": float(np.percentile(corr_values, 16)),
                "corr_p84": float(np.percentile(corr_values, 84)),
                "positive_fraction": float(np.mean(corr_values > 0.0)),
                "slope_median": float(np.median(slope_values)),
                "slope_p16": float(np.percentile(slope_values, 16)),
                "slope_p84": float(np.percentile(slope_values, 84)),
            }
        )
    return rows


def write_scan_csv(rows: list[dict[str, float]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_figure(rows: list[dict[str, float]], output_path: Path) -> dict[str, tuple[float, float]]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cases = {
        "Inflow only\n(Forbes-2014 limit)": (0.28, 0.0, 2026),
        "Efficiency only\n(two-driver extension)": (0.0, 0.28, 2027),
        "Both; efficiency dominates": (0.28, 0.42, 2028),
    }

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 4, figsize=(15.2, 3.8), constrained_layout=True)
    results: dict[str, tuple[float, float]] = {}
    colors = ["#3A6EA5", "#D1495B", "#2A9D8F"]
    for axis, (label, (sigma_phi, sigma_epsilon, seed)), color in zip(
        axes[:3], cases.items(), colors
    ):
        x, y = simulate_series(sigma_phi, sigma_epsilon, seed)
        corr, slope = correlation_and_slope(x, y)
        results[label.replace("\n", " ")] = (corr, slope)
        axis.hexbin(x, y, gridsize=38, mincnt=1, cmap="Greys", linewidths=0)
        x_line = np.linspace(np.percentile(x, 1), np.percentile(x, 99), 100)
        axis.plot(x_line, slope * x_line, color=color, lw=2.4)
        axis.axhline(0, color="#777777", lw=0.6)
        axis.axvline(0, color="#777777", lw=0.6)
        axis.set_title(label, fontsize=10.5, weight="bold")
        axis.set_xlabel(r"$\Delta\log_{10}{\rm SFR}$")
        axis.text(
            0.05,
            0.93,
            rf"$r={corr:+.2f}$" + "\n" + rf"$m={slope:+.3f}$",
            transform=axis.transAxes,
            va="top",
            fontsize=9.3,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82},
        )
    axes[0].set_ylabel(r"$\Delta\log_{10} Z_{\rm g}$")

    ratio = np.asarray([row["sigma_epsilon_over_sigma_phi"] for row in rows])
    median = np.asarray([row["corr_median"] for row in rows])
    p16 = np.asarray([row["corr_p16"] for row in rows])
    p84 = np.asarray([row["corr_p84"] for row in rows])
    axes[3].fill_between(ratio, p16, p84, color="#8AB17D", alpha=0.35, label="16–84%")
    axes[3].plot(ratio, median, color="#264653", lw=2.5, label="median")
    axes[3].axhline(0.0, color="#777777", lw=0.8)
    axes[3].axvline(1.0, color="#D1495B", lw=1.0, ls="--", label="equal rms")
    axes[3].set_title("Driver competition", fontsize=10.5, weight="bold")
    axes[3].set_xlabel(r"$\sigma_{\ln\epsilon}/\sigma_{\ln\Phi}$")
    axes[3].set_ylabel(r"Pearson $r(\Delta\log {\rm SFR},\Delta\log Z_{\rm g})$")
    axes[3].set_xlim(float(ratio.min()), float(ratio.max()))
    axes[3].set_ylim(-1.0, 1.0)
    axes[3].legend(loc="lower right", frameon=True, fontsize=8)

    fig.suptitle(
        "A fixed-mass regulator changes sign when efficiency fluctuations dominate inflow fluctuations",
        fontsize=12.5,
        weight="bold",
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    ratios = np.linspace(0.0, 2.5, 21)
    rows = scan_driver_ratio(ratios)
    csv_path = output_dir / "20260813_FMR_driver_competition_scan.csv"
    figure_path = output_dir / "20260813_FMR_driver_competition_toy_model.png"
    write_scan_csv(rows, csv_path)
    cases = make_figure(rows, figure_path)

    sign_change = next(
        (
            row["sigma_epsilon_over_sigma_phi"]
            for row in rows
            if row["corr_median"] > 0.0
        ),
        float("nan"),
    )
    print(f"tau_eq={TAU_EQ:.4f}; tau_driver={2 * TAU_EQ:.4f}")
    for label, (corr, slope) in cases.items():
        print(f"{label}: r={corr:+.4f}, slope={slope:+.5f}")
    print(f"first sampled positive median correlation ratio={sign_change:.3f}")
    print(figure_path)
    print(csv_path)


if __name__ == "__main__":
    main()
