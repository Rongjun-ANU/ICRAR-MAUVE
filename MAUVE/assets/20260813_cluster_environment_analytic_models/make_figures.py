#!/usr/bin/env python3
"""Create the analytical-environment review figures and provenance record."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ASSET_DIR = Path(__file__).resolve().parent
REGISTRY = ASSET_DIR / "source_registry.csv"
COLORS = {
    "navy": "#16324F",
    "blue": "#2F6690",
    "cyan": "#3A8D9D",
    "green": "#4A7C59",
    "gold": "#D5A021",
    "orange": "#C76D3A",
    "red": "#A23E48",
    "plum": "#6D597A",
    "ink": "#20252A",
    "paper": "#F7F4EE",
    "grid": "#D8D4CB",
}


def save(fig, stem):
    fig.savefig(ASSET_DIR / f"{stem}.png", dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(ASSET_DIR / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def broad_family(family):
    if family in {
        "transport", "tides", "cluster_potential", "cluster_dynamics",
        "starvation", "preprocessing",
    }:
        return "Foundations and environment"
    if family in {
        "cold_disk_rps", "cold_and_hot_rps", "hot_halo_rps", "cluster_icm",
        "model_validation", "star_formation_response",
    }:
        return "RPS equations and tests"
    if family in {"semi_analytic_population", "model_comparison"}:
        return "Population semi-analytics"
    if family in {
        "phase_space", "quenching_clock", "chemical_clock", "gas_regulator",
        "candidate_diagnostic",
    }:
        return "History and diagnostics"
    if family == "individual_inversion":
        return "Tailored inversions"
    return "Reviews / corrections"


def literature_timeline(rows):
    groups = [
        "Foundations and environment",
        "RPS equations and tests",
        "Population semi-analytics",
        "History and diagnostics",
        "Tailored inversions",
        "Reviews / corrections",
    ]
    palette = [COLORS["navy"], COLORS["blue"], COLORS["green"], COLORS["gold"], COLORS["red"], COLORS["plum"]]
    bins = list(range(1970, 2030, 5))
    centres = np.asarray(bins[:-1]) + 2.5
    counts = {group: np.zeros(len(centres), dtype=int) for group in groups}
    for row in rows:
        year = int(row["year"])
        idx = min(max((year - bins[0]) // 5, 0), len(centres) - 1)
        counts[broad_family(row["family"])][idx] += 1

    fig, ax = plt.subplots(figsize=(11.5, 5.8), constrained_layout=True)
    bottom = np.zeros(len(centres), dtype=int)
    for group, colour in zip(groups, palette):
        ax.bar(centres, counts[group], width=4.35, bottom=bottom, color=colour,
               edgecolor="white", linewidth=0.6, label=group)
        bottom += counts[group]
    ax.set_xlim(1970, 2028)
    ax.set_ylim(0, max(bottom) + 2.2)
    ax.set_xticks(np.arange(1970, 2030, 5))
    ax.set_xlabel("Publication year (five-year bins)")
    ax.set_ylabel("Included papers")
    ax.set_title("Analytical cluster-environment literature: from criteria to inference")
    ax.grid(axis="y", color=COLORS["grid"], linewidth=0.7)
    ax.set_axisbelow(True)
    ax.legend(ncol=3, frameon=False, loc="upper left", fontsize=8.6)
    ax.text(0.995, 0.98, f"{len(rows)} DOI-verified papers; search closed 2026-08-13",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.4, color="#555555")
    save(fig, "figure_01_literature_timeline")


def add_box(ax, xy, width, height, title, lines, colour):
    x, y = xy
    box = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        linewidth=1.25, edgecolor=colour, facecolor="white",
    )
    ax.add_patch(box)
    ax.text(x + 0.025, y + height - 0.055, title, ha="left", va="top",
            fontsize=11.3, fontweight="bold", color=colour)
    ax.text(x + 0.025, y + height - 0.105, "\n".join(lines), ha="left", va="top",
            fontsize=8.3, color=COLORS["ink"], linespacing=1.38)
    return box


def model_stack():
    fig, ax = plt.subplots(figsize=(12.2, 6.9), constrained_layout=True)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("A model stack for converting cluster observables into evolutionary constraints",
                 fontsize=15, loc="left", pad=13, color=COLORS["ink"])

    xvals = [0.035, 0.235, 0.435, 0.635, 0.835]
    width, height, y = 0.155, 0.58, 0.245
    specs = [
        ("Observed galaxy", ["Sigma_star(R)", "HI / CO / Halpha", "Z_gas(R), SFH", "inclination, PA"], COLORS["navy"]),
        ("Cluster and orbit", ["rho_ICM(r, theta)", "potential and R200", "R_proj, Delta v_los", "latent 3D trajectory"], COLORS["blue"]),
        ("Removal physics", ["restoring force", "pressure pulse", "viscous / thermal", "tidal terms"], COLORS["orange"]),
        ("Evolution response", ["gas continuity", "SF consumption", "metal enrichment", "fallback / mixing"], COLORS["green"]),
        ("Posterior outputs", ["R_strip, f_strip", "p_max, impulse", "time to/since peak", "mechanism mixture"], COLORS["red"]),
    ]
    for idx, (title, lines, colour) in enumerate(specs):
        add_box(ax, (xvals[idx], y), width, height, title, lines, colour)
        if idx < len(specs) - 1:
            arrow = FancyArrowPatch((xvals[idx] + width + 0.008, 0.535), (xvals[idx + 1] - 0.008, 0.535),
                                    arrowstyle="-|>", mutation_scale=14, lw=1.3, color="#666666")
            ax.add_patch(arrow)

    ax.text(0.035, 0.155, "Marginalize rather than fix:", fontsize=10.3, fontweight="bold", color=COLORS["plum"])
    ax.text(0.255, 0.155,
            "line-of-sight position, transverse speed, wind-disk angle, ICM substructure, initial gas profile, feedback, and pre-processing",
            fontsize=8.8, color=COLORS["ink"], va="center")
    ax.text(0.035, 0.085,
            "Key distinction: equations can predict susceptibility; a unique history requires multiple spatially resolved observables and explicit priors.",
            fontsize=9.4, color="#4B4B4B", style="italic")
    save(fig, "figure_02_model_stack")


def rps_regimes():
    x = np.linspace(-1.7, 1.7, 500)
    y = np.linspace(-1.7, 1.7, 500)
    xx, yy = np.meshgrid(x, y)
    quasi = xx
    impulse = xx + yy
    score = np.where(yy >= 0, quasi, impulse)
    cmap = LinearSegmentedColormap.from_list("strip", ["#E8EEF3", "#F4D7B4", "#B44B4E"])

    fig, ax = plt.subplots(figsize=(8.6, 7.1), constrained_layout=True)
    im = ax.contourf(xx, yy, score, levels=np.linspace(-2.2, 2.2, 17), cmap=cmap, extend="both")
    ax.axhline(0, color=COLORS["ink"], lw=1.4)
    ax.plot([0, 0], [0, 1.7], color=COLORS["ink"], lw=2.0)
    ax.plot(x[x <= 0], -x[x <= 0], color=COLORS["ink"], lw=2.0)
    ax.plot(x[x >= 0], -x[x >= 0], color=COLORS["ink"], lw=1.0, ls="--", alpha=0.6)
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.7, 1.7)
    ax.set_xlabel(r"log$_{10}$ peak-pressure ratio:  $P_{max}/P_{restore}$")
    ax.set_ylabel(r"log$_{10}$ pulse-duration ratio:  $\tau_{pulse}/\tau_z$")
    ax.set_title("Two limiting stripping regimes in a time-dependent ram-pressure pulse")
    ax.text(-1.55, 1.42, "Long pulse / quasi-static", fontsize=10.5, fontweight="bold")
    ax.text(-1.55, -1.48, "Short pulse / impulse-controlled", fontsize=10.5, fontweight="bold")
    ax.text(-1.45, 0.35, "bound or weakly affected", fontsize=9.5, color=COLORS["navy"])
    ax.text(0.58, 0.65, "pressure exceeds\nrestoring force", fontsize=9.5, color="#6F2227", ha="center")
    ax.text(0.85, -0.75, "integrated impulse\nunbinds gas", fontsize=9.5, color="#6F2227", ha="center")
    ax.text(-1.62, -1.62,
            "Schematic normalization only. Radius, gas column, vertical potential and pulse shape move the boundaries.",
            fontsize=8.2, color="#555555")
    cbar = fig.colorbar(im, ax=ax, shrink=0.82, pad=0.03)
    cbar.set_label("qualitative stripping susceptibility")
    save(fig, "figure_03_rps_regimes")


def applicability_matrix():
    rows = [
        "Gunn-Gott disk",
        "Finite-pulse / impulse",
        "Hot-halo stripping",
        "Transport + tides",
        "Projected phase space",
        "Quenching / chemical clocks",
        "Semi-analytic populations",
        "Tailored dynamical inversion",
    ]
    cols = [
        "Current\nR_strip", "Gas\nloss", "Peak\npressure", "Infall /\npeak time",
        "SFH", "Metallicity", "Population\nstats", "Individual\nreadiness",
    ]
    values = np.array([
        [3, 2, 1, 0, 0, 0, 1, 1],
        [3, 3, 2, 1, 0, 0, 1, 2],
        [1, 3, 2, 1, 2, 1, 3, 1],
        [0, 2, 0, 0, 1, 1, 2, 0],
        [0, 0, 0, 2, 1, 0, 3, 1],
        [0, 1, 0, 2, 3, 3, 3, 1],
        [1, 3, 2, 2, 3, 2, 3, 0],
        [3, 3, 3, 3, 3, 2, 0, 3],
    ], dtype=float)
    cmap = LinearSegmentedColormap.from_list("app", ["#F3F1EC", "#AFC8D6", "#4A7C59", "#A23E48"])
    fig, ax = plt.subplots(figsize=(11.2, 6.2), constrained_layout=True)
    im = ax.imshow(values, cmap=cmap, vmin=0, vmax=3, aspect="auto")
    ax.set_xticks(range(len(cols)), cols, fontsize=9)
    ax.set_yticks(range(len(rows)), rows, fontsize=9.5)
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, length=0)
    ax.set_title("What each model family can constrain without over-interpretation", loc="left", pad=18)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = int(values[i, j])
            ax.text(j, i, ["-", "limited", "useful", "strong"][value],
                    ha="center", va="center", fontsize=8.1,
                    color="white" if value >= 2 else COLORS["ink"],
                    fontweight="bold" if value == 3 else "normal")
    for edge in np.arange(-0.5, len(cols), 1):
        ax.axvline(edge, color="white", lw=1.1)
    for edge in np.arange(-0.5, len(rows), 1):
        ax.axhline(edge, color="white", lw=1.1)
    cbar = fig.colorbar(im, ax=ax, shrink=0.78, pad=0.025, ticks=[0, 1, 2, 3])
    cbar.ax.set_yticklabels(["not an output", "limited", "useful", "strong"])
    ax.text(0, 1.075,
            "Scores summarize the reviewed literature; they are not performance measurements and assume suitable input data.",
            transform=ax.transAxes, fontsize=8.4, color="#555555")
    save(fig, "figure_04_applicability_matrix")


def main():
    with REGISTRY.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    literature_timeline(rows)
    model_stack()
    rps_regimes()
    applicability_matrix()
    counts = Counter(broad_family(row["family"]) for row in rows)
    provenance = {
        "generated_on": "2026-08-13",
        "generator": str(Path(__file__).resolve()),
        "source_registry": str(REGISTRY.resolve()),
        "source_count": len(rows),
        "broad_family_counts": dict(counts),
        "figures": {
            "figure_01_literature_timeline": "Counts of DOI-verified sources by publication year and broad model family.",
            "figure_02_model_stack": "Author synthesis of the model layers needed for individual-galaxy inference.",
            "figure_03_rps_regimes": "Qualitative synthesis of quasi-static and impulse limits after Takeda et al. (1984) and Koppen et al. (2018).",
            "figure_04_applicability_matrix": "Author assessment of output scope; qualitative, not a benchmark.",
        },
    }
    (ASSET_DIR / "figure_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"figures": 4, "source_count": len(rows), "family_counts": dict(counts)}))


if __name__ == "__main__":
    main()
