#!/usr/bin/env python3
"""Assemble the narrative review and registry-derived appendices into Markdown."""

from __future__ import annotations

import csv
import html
import json
import re
from collections import Counter
from pathlib import Path


ASSET_DIR = Path(__file__).resolve().parent
REPORT_DIR = ASSET_DIR.parent.parent
OUTPUT = REPORT_DIR / (
    "20260813 Analytical Models of Galaxy Evolution in Cluster Environments "
    "and Ram Pressure Stripping.md"
)

FAMILY_LABELS = {
    "cold_disk_rps": "Cold-disk ram-pressure equations",
    "cold_and_hot_rps": "Cold- and hot-gas ram-pressure equations",
    "hot_halo_rps": "Hot-halo stripping and starvation",
    "star_formation_response": "Star-formation response to pressure",
    "transport": "Transport, evaporation, and viscous stripping",
    "tides": "Tides and harassment",
    "cluster_potential": "Cluster potential",
    "cluster_dynamics": "Cluster dynamics and caustics",
    "cluster_icm": "ICM and ram-pressure profiles",
    "preprocessing": "Preprocessing",
    "starvation": "Starvation and gas-supply cutoff",
    "model_validation": "Analytic-model validation",
    "semi_analytic_population": "Semi-analytic population models",
    "model_comparison": "Population-model comparisons",
    "phase_space": "Projected-phase-space inference",
    "quenching_clock": "Quenching clocks",
    "gas_regulator": "Gas regulators",
    "chemical_clock": "Chemical clocks",
    "candidate_diagnostic": "Candidate diagnostics",
    "individual_inversion": "Tailored individual-galaxy inversions",
    "review": "Reviews",
    "errata": "Published corrections",
}

FAMILY_ORDER = [
    "cold_disk_rps", "cold_and_hot_rps", "hot_halo_rps",
    "star_formation_response", "transport", "tides", "cluster_potential",
    "cluster_dynamics", "cluster_icm", "preprocessing", "starvation",
    "model_validation", "semi_analytic_population", "model_comparison",
    "phase_space", "quenching_clock", "gas_regulator", "chemical_clock",
    "candidate_diagnostic", "individual_inversion", "review", "errata",
]

READINESS = {
    "cold_disk_rps": "High for a susceptibility or stripping-radius map; low for a unique history.",
    "cold_and_hot_rps": "Useful for comparative susceptibility if galaxy and host profiles are constrained.",
    "hot_halo_rps": "Useful as a starvation prior; usually not directly observable for one galaxy.",
    "star_formation_response": "Useful as a parametric burst/fading prior, not a unique mechanism diagnostic.",
    "transport": "Order-of-magnitude use unless viscosity, magnetic suppression, and phase structure are constrained.",
    "tides": "Useful when old-stellar morphology and kinematics constrain the collisionless response.",
    "cluster_potential": "A required environmental prior, not an evolutionary solution by itself.",
    "cluster_dynamics": "Useful for membership and broad orbital priors; projection remains substantial.",
    "cluster_icm": "Useful as a mean pressure prior; local substructure must be marginalized.",
    "preprocessing": "Useful as a population or prior component; individual prior hosts are often uncertain.",
    "starvation": "Useful for response times after supply cutoff; does not identify what stopped the supply.",
    "model_validation": "Defines where simpler equations are reliable and where model discrepancy is required.",
    "semi_analytic_population": "Strong for populations and priors; weak as a direct inversion for one object.",
    "model_comparison": "Strong for discriminating population trends; not a unique individual history.",
    "phase_space": "Probabilistic only; report full infall/backsplash/interloper distributions.",
    "quenching_clock": "A statistical clock whose zero-point and star-formation history must be stated.",
    "gas_regulator": "Useful for gas-supply and consumption histories with resolved gas/SFR constraints.",
    "chemical_clock": "Complementary timing evidence; selection and yield/outflow assumptions matter.",
    "candidate_diagnostic": "Useful for prioritization and completeness tests, not causal proof.",
    "individual_inversion": "Highest individual readiness when several independent resolved observables agree.",
    "review": "Context and citation chaining; primary papers remain authoritative for numerical reuse.",
    "errata": "Mandatory companion for quantitative reuse of the corrected paper.",
}

TITLE_OVERRIDES = {
    "fossati_etal_2018": (
        "A Virgo Environmental Survey Tracing Ionised Gas Emission (VESTIGE). "
        "II. Constraining the quenching time in the stripped galaxy NGC 4330"
    ),
    "vollmer_etal_2018": (
        "Two uneven sisters. I. NGC 4388 - a strongly constrained ram-pressure "
        "stripping event"
    ),
    "vollmer_etal_2021": (
        "A Virgo Environmental Survey Tracing Ionised Gas Emission (VESTIGE). "
        "VIII. Modeling ram pressure stripping of diffuse gas in the Virgo cluster "
        "spiral galaxy NGC 4330"
    ),
}


def clean(value: str) -> str:
    value = html.unescape(value or "")
    value = re.sub(r"<[^>]+>", "", value)
    value = value.replace("\u00a0", " ")
    value = re.sub(r"\s+", " ", value).strip()
    return value


def authors_display(value: str, short: bool = False) -> str:
    names = []
    for raw in value.split("; "):
        raw = clean(raw)
        if not raw:
            continue
        if ", " in raw:
            family, given = raw.split(", ", 1)
            names.append(f"{given} {family}")
        else:
            names.append(raw)
    if short and len(names) > 3:
        return f"{names[0]} et al."
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return " & ".join(names)
    return ", ".join(names[:-1]) + ", & " + names[-1]


def bibliographic_line(row: dict, number: int, short_authors: bool = False) -> str:
    title = TITLE_OVERRIDES.get(row["key"], clean(row["title"]))
    authors = authors_display(row["authors"], short=short_authors)
    journal = clean(row["journal"])
    volume = clean(row["volume"])
    issue = clean(row["issue"])
    page = clean(row["page"])
    locator = volume
    if issue:
        locator += f"({issue})"
    if page:
        locator += f", {page}"
    return (
        f"**[{number:03d}] {authors} ({row['year']}).** "
        f"\"{title}.\" *{journal}* {locator}. "
        f"[DOI: {row['doi']}](https://doi.org/{row['doi']})."
    )


def main() -> None:
    body = (ASSET_DIR / "review_body.md").read_text(encoding="utf-8").rstrip()
    with (ASSET_DIR / "source_registry.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    with (ASSET_DIR / "web_spot_checks.tsv").open(newline="", encoding="utf-8") as handle:
        spot_checks = list(csv.DictReader(handle, delimiter="\t"))
    audit = json.loads((ASSET_DIR / "doi_resolution_audit.json").read_text(encoding="utf-8"))

    if len(rows) != 80 or audit["n_resolved"] != 80 or audit["n_failed"] != 0:
        raise RuntimeError("Expected a clean 80/80 DOI registry before assembly")

    numbered = {row["key"]: index for index, row in enumerate(rows, start=1)}
    appendix = [
        "",
        "# Appendix A. Source-by-source model register",
        "",
        (
            "This register is generated from the verified bibliography. The model class and "
            "inclusion role are review judgements; the bibliographic fields and DOI are from the "
            "2026-08-13 Crossref audit. The readiness statement says what the paper can safely "
            "contribute to an individual-galaxy analysis, not whether the paper is scientifically important."
        ),
        "",
    ]

    for family in FAMILY_ORDER:
        family_rows = [row for row in rows if row["family"] == family]
        if not family_rows:
            continue
        appendix.extend([f"## A.{FAMILY_ORDER.index(family) + 1} {FAMILY_LABELS[family]}", ""])
        for row in family_rows:
            number = numbered[row["key"]]
            appendix.extend([
                bibliographic_line(row, number, short_authors=True),
                "",
                f"- **Why included:** {clean(row['inclusion_role'])}.",
                f"- **Individual-galaxy use:** {READINESS[family]}",
                "",
            ])

    appendix.extend([
        "## A.23 Foundational DOI-less ICM model",
        "",
        (
            "**[081] A. Cavaliere & R. Fusco-Femiano (1976).** \"X-rays from hot "
            "plasma in clusters of galaxies.\" *Astronomy & Astrophysics* **49**, "
            "137-144. [ADS record](https://ui.adsabs.harvard.edu/abs/1976A%26A....49..137C/abstract)."
        ),
        "",
        "- **Why included:** foundational beta-profile representation of a cluster ICM.",
        "- **Individual-galaxy use:** a transparent radial density prior, not a local density measurement.",
        "",
        "# Appendix B. Complete verified references",
        "",
    ])

    for row in sorted(rows, key=lambda item: (int(item["year"]), clean(item["title"]).lower())):
        appendix.extend([bibliographic_line(row, numbered[row["key"]]), ""])
    appendix.extend([
        (
            "**[081] Cavaliere, A., & Fusco-Femiano, R. (1976).** \"X-rays from hot "
            "plasma in clusters of galaxies.\" *Astronomy & Astrophysics* **49**, 137-144. "
            "[ADS record](https://ui.adsabs.harvard.edu/abs/1976A%26A....49..137C/abstract)."
        ),
        "",
        "# Appendix C. Reproducibility and artifact manifest",
        "",
        f"- DOI seed records: **{audit['n_seeded_dois']}**.",
        f"- DOI records resolved through Crossref: **{audit['n_resolved']}**.",
        f"- Failed DOI resolutions: **{audit['n_failed']}**.",
        "- Foundational DOI-less records checked through ADS: **1**.",
        f"- Primary-source web spot checks recorded: **{len(spot_checks)}** "
        f"({100 * len(spot_checks) / (len(rows) + 1):.1f} percent of the 81-work corpus).",
        "- Search and verification closed: **2026-08-13 (Australia/Perth)**.",
        "",
        "## C.1 Machine-readable evidence",
        "",
        "- [DOI seeds](assets/20260813_cluster_environment_analytic_models/doi_seeds.tsv)",
        "- [Crossref DOI audit](assets/20260813_cluster_environment_analytic_models/doi_resolution_audit.json)",
        "- [Verified source registry](assets/20260813_cluster_environment_analytic_models/source_registry.csv)",
        "- [Primary-source spot checks](assets/20260813_cluster_environment_analytic_models/web_spot_checks.tsv)",
        "- [Figure provenance](assets/20260813_cluster_environment_analytic_models/figure_provenance.json)",
        "- [Acceptance audit](assets/20260813_cluster_environment_analytic_models/acceptance_audit.md)",
        "",
        "## C.2 Reproducible scripts",
        "",
        "- [Bibliography verifier](assets/20260813_cluster_environment_analytic_models/verify_bibliography.py)",
        "- [Figure generator](assets/20260813_cluster_environment_analytic_models/make_figures.py)",
        "- [Markdown assembler](assets/20260813_cluster_environment_analytic_models/assemble_review.py)",
        "- [PDF builder](assets/20260813_cluster_environment_analytic_models/build_pdf.py)",
        "",
        "## C.3 Figure files",
        "",
        "Each figure is supplied as PNG for the Markdown report and as PDF for reuse:",
        "",
        "1. Literature timeline (`figure_01_literature_timeline.*`).",
        "2. Individual-galaxy model stack (`figure_02_model_stack.*`).",
        "3. Peak-pressure versus impulse regimes (`figure_03_rps_regimes.*`).",
        "4. Model-family applicability matrix (`figure_04_applicability_matrix.*`).",
        "",
        "## C.4 Scope boundary",
        "",
        (
            "This is a deep systematic-style scoping review, not a formal claim of universal "
            "bibliographic completeness. Pure simulations and observational catalogues are represented "
            "only when they validate, calibrate, or operationalize an analytic or semi-analytic model. "
            "The registry is designed so that future papers can be added without changing the distinctions "
            "between analytic physics, semi-analytic populations, simulation-calibrated diagnostics, "
            "empirical clocks, and tailored inversions."
        ),
        "",
    ])

    OUTPUT.write_text(body + "\n" + "\n".join(appendix), encoding="utf-8")
    counts = Counter(row["family"] for row in rows)
    print(json.dumps({
        "output": str(OUTPUT),
        "doi_records": len(rows),
        "total_works": len(rows) + 1,
        "spot_checks": len(spot_checks),
        "family_count": len(counts),
    }))


if __name__ == "__main__":
    main()
