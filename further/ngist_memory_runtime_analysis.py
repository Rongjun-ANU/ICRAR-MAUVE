"""Helpers for the 2026-05-28 nGIST memory/runtime notebook."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import math
import re
from typing import Iterable

import numpy as np
import pandas as pd
from astropy.io import fits


EFFICIENCY_CUTOFF_JOBID = 43057897
DATETIME_FORMAT = "%m/%d/%y %H:%M:%S"

HIGHMEM_GALAXIES = {"NGC4192", "NGC4501"}
LONG_GALAXIES = {
    "NGC4293",
    "NGC4298",
    "NGC4302",
    "NGC4330",
    "NGC4383",
    "NGC4396",
    "NGC4419",
    "NGC4457",
    "NGC4567_8",
    "NGC4698",
}


def mem_to_gib(value: str) -> float:
    """Convert Slurm memory strings such as 220GiB or 21.1MiB to GiB."""
    if value == "N/A":
        return math.nan
    match = re.fullmatch(r"([0-9.]+)([A-Za-z]+)", value)
    if not match:
        return math.nan
    number = float(match.group(1))
    unit = match.group(2)
    if unit == "GiB":
        return number
    if unit == "MiB":
        return number / 1024.0
    if unit == "KiB":
        return number / 1024.0 / 1024.0
    if unit == "TiB":
        return number * 1024.0
    return math.nan


def _percent_to_float(value: str) -> float:
    if value == "N/A%":
        return math.nan
    return float(value.rstrip("%"))


def galaxy_from_jobname(jobname: str) -> str:
    return jobname.split("_v3tk", 1)[0]


def variant_from_jobname(jobname: str) -> str:
    if "bintest" in jobname:
        return "bin_test"
    if "spaxeltest" in jobname or "spax+" in jobname:
        return "spaxel_test"
    return "regular"


def parse_efficiency_table(path: Path, cutoff_jobid: int = EFFICIENCY_CUTOFF_JOBID) -> list[dict]:
    """Parse the recent Slurm efficiency rows down to and including cutoff_jobid."""
    row_re = re.compile(
        r"^(\d+)\s+(\S+)\s+(\S+)\s+(.+?)\s+"
        r"(\d+)\s+(\d+)\s+([0-9.]+)\s+([0-9.]+%)\s+([0-9.]+%)\s+"
        r"(\S+)\s+(\S+)\s+(\S+)\s+(\S+)"
    )
    rows: list[dict] = []
    for line in Path(path).read_text().splitlines():
        match = row_re.match(line)
        if not match:
            continue
        (
            jobid,
            jobname,
            state,
            end_text,
            reqcpu,
            alloc,
            effcores,
            effreq,
            effalloc,
            reqmem,
            averss,
            maxrss,
            mempct,
        ) = match.groups()
        jobid_int = int(jobid)
        galaxy = galaxy_from_jobname(jobname)
        reqmem_gib = mem_to_gib(reqmem)
        maxrss_gib = mem_to_gib(maxrss)
        rows.append(
            {
                "jobid": jobid_int,
                "jobname": jobname,
                "galaxy": galaxy,
                "job_variant": variant_from_jobname(jobname),
                "state": state,
                "end_text": end_text.strip(),
                "reqcpu": int(reqcpu),
                "alloc": int(alloc),
                "effcores": float(effcores),
                "effreq_pct": _percent_to_float(effreq),
                "effalloc_pct": _percent_to_float(effalloc),
                "reqmem_gib": reqmem_gib,
                "averss_gib": mem_to_gib(averss),
                "maxrss_gib": maxrss_gib,
                "mem_pct": _percent_to_float(mempct),
                "queue_inferred": infer_queue(galaxy, reqmem_gib),
                "standard_mem_cap": bool(
                    math.isfinite(reqmem_gib)
                    and math.isfinite(maxrss_gib)
                    and reqmem_gib <= 230
                    and maxrss_gib >= 0.98 * reqmem_gib
                ),
            }
        )
        if jobid_int == cutoff_jobid:
            break
    return rows


def infer_queue(galaxy: str, reqmem_gib: float) -> str:
    if math.isfinite(reqmem_gib) and reqmem_gib >= 900:
        return "highmem"
    if galaxy in LONG_GALAXIES:
        return "long/work_220GiB"
    if galaxy in HIGHMEM_GALAXIES:
        return "highmem_expected"
    return "work_220GiB"


def parse_config_value(path: Path, section: str, key: str) -> str | None:
    if not path.exists():
        return None
    current_section = None
    for raw_line in path.read_text(errors="replace").splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line:
            continue
        if not line.startswith(" ") and line.endswith(":"):
            current_section = line[:-1].strip()
            continue
        if current_section != section:
            continue
        match = re.match(rf"\s*{re.escape(key)}\s*:\s*(.+)\s*$", line)
        if match:
            return match.group(1).strip().strip("'\"")
    return None


def collect_cube_metrics(root: Path) -> pd.DataFrame:
    """Collect cube/map metrics from local v3tk_v7.6.8 product directories."""
    records: list[dict] = []
    for directory in sorted(Path(root).iterdir()):
        if not directory.is_dir():
            continue
        galaxy = directory.name
        binning_path = directory / f"{galaxy}_spatial_binning_maps.fits"
        mask_path = directory / f"{galaxy}_mask.fits"
        if not binning_path.exists() or not mask_path.exists():
            continue

        with fits.open(binning_path, memmap=True) as hdul:
            binid = hdul["BINID"].data
            ny, nx = binid.shape
            valid_binid = np.isfinite(binid) & (binid >= 0)
            available_spaxels = int(valid_binid.sum())
            n_bins = int(len(np.unique(binid[valid_binid]))) if available_spaxels else 0

        with fits.open(mask_path, memmap=True) as hdul:
            mask_table = hdul[1].data
            mask_zero = int((mask_table["MASK"] == 0).sum())
            mask_rows = int(len(mask_table))

        gas_bin_maps = directory / f"{galaxy}_gas_bin_maps.fits"
        gas_spaxel_maps = directory / f"{galaxy}_gas_spaxel_maps.fits"
        records.append(
            {
                "galaxy": galaxy,
                "ny": int(ny),
                "nx": int(nx),
                "shape": f"{ny}x{nx}",
                "total_cube_spaxels": int(ny * nx),
                "mask_rows": mask_rows,
                "available_spaxels": available_spaxels,
                "mask_zero_spaxels": mask_zero,
                "mask_matches_binid": bool(mask_zero == available_spaxels),
                "n_bins": n_bins,
                "gas_level_config": parse_config_value(directory / "CONFIG", "GAS", "LEVEL"),
                "input_cube": parse_config_value(directory / "CONFIG", "GENERAL", "INPUT"),
                "has_gas_bin_maps": gas_bin_maps.exists(),
                "has_gas_spaxel_maps": gas_spaxel_maps.exists(),
                "has_sfh_maps": (directory / f"{galaxy}_sfh_maps.fits").exists(),
                "has_kin_maps": (directory / f"{galaxy}_kin_maps.fits").exists(),
            }
        )
    if not records:
        return pd.DataFrame().set_index(pd.Index([], name="galaxy"))
    return pd.DataFrame.from_records(records).set_index("galaxy").sort_index()


def parse_log_time(line: str) -> datetime | None:
    match = re.match(r"(\d\d/\d\d/\d\d \d\d:\d\d:\d\d)", line)
    if not match:
        return None
    return datetime.strptime(match.group(1), DATETIME_FORMAT)


@dataclass
class AttemptState:
    galaxy: str
    attempt_index: int
    start_line: int
    start_time: datetime | None
    last_line: int
    last_time: datetime | None
    skipped_modules: int = 0
    current_gas_level: str | None = None
    read_start: datetime | None = None
    read_done: datetime | None = None
    read_total_spectra: int | None = None
    bin_level_start: datetime | None = None
    bin_fit_start: datetime | None = None
    bin_ppxf_done: datetime | None = None
    bin_ppxf_spectra: int | None = None
    bin_result_write: datetime | None = None
    bin_bestfit_write: datetime | None = None
    spaxel_level_start: datetime | None = None
    spaxel_fit_start: datetime | None = None
    spaxel_ppxf_done: datetime | None = None
    spaxel_ppxf_spectra: int | None = None
    spaxel_result_write: datetime | None = None
    spaxel_bestfit_write: datetime | None = None
    gas_maps_start: datetime | None = None
    sfh_start: datetime | None = None
    sfh_done: datetime | None = None

    def to_record(self) -> dict:
        record = {
            "galaxy": self.galaxy,
            "attempt_index": self.attempt_index,
            "start_line": self.start_line,
            "last_line": self.last_line,
            "skipped_modules": self.skipped_modules,
            "read_total_spectra": self.read_total_spectra,
            "bin_ppxf_spectra": self.bin_ppxf_spectra,
            "spaxel_ppxf_spectra": self.spaxel_ppxf_spectra,
            "has_bin_fit": self.bin_ppxf_done is not None,
            "has_spaxel_fit": self.spaxel_ppxf_done is not None,
            "has_spaxel_bestfit_write": self.spaxel_bestfit_write is not None,
            "has_gas_maps_start": self.gas_maps_start is not None,
            "has_sfh_done": self.sfh_done is not None,
            "is_continuation_only": (
                self.bin_level_start is None
                and (self.skipped_modules > 0 or self.sfh_start is not None)
            ),
        }
        for name in (
            "start_time",
            "last_time",
            "read_start",
            "read_done",
            "bin_level_start",
            "bin_fit_start",
            "bin_ppxf_done",
            "bin_result_write",
            "bin_bestfit_write",
            "spaxel_level_start",
            "spaxel_fit_start",
            "spaxel_ppxf_done",
            "spaxel_result_write",
            "spaxel_bestfit_write",
            "gas_maps_start",
            "sfh_start",
            "sfh_done",
        ):
            value = getattr(self, name)
            record[name] = value
        record.update(
            {
                "wall_hours": hours_between(self.start_time, self.last_time),
                "read_hours": hours_between(self.read_start, self.read_done),
                "bin_fit_hours": hours_between(self.bin_fit_start, self.bin_ppxf_done),
                "bin_bestfit_write_hours": hours_between(
                    self.bin_result_write, self.bin_bestfit_write
                ),
                "spaxel_fit_hours": hours_between(
                    self.spaxel_fit_start, self.spaxel_ppxf_done
                ),
                "spaxel_result_write_hours": hours_between(
                    self.spaxel_ppxf_done, self.spaxel_result_write
                ),
                "spaxel_bestfit_write_hours": hours_between(
                    self.spaxel_result_write, self.spaxel_bestfit_write
                ),
                "bestfit_to_maps_hours": hours_between(
                    self.spaxel_bestfit_write, self.gas_maps_start
                ),
                "sfh_hours": hours_between(self.sfh_start, self.sfh_done),
            }
        )
        return record


def hours_between(start: datetime | None, end: datetime | None) -> float:
    if start is None or end is None:
        return math.nan
    return (end - start).total_seconds() / 3600.0


def parse_logfile_attempts(path: Path, galaxy: str | None = None) -> list[dict]:
    """Parse LOGFILE into attempts and extract key read/fit/write stage times."""
    path = Path(path)
    galaxy = galaxy or path.parent.name
    attempts: list[AttemptState] = []
    current: AttemptState | None = None

    def finish_current() -> None:
        if current is not None:
            attempts.append(current)

    for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        timestamp = parse_log_time(line)
        if "_readData: Using the read-in routine" in line:
            finish_current()
            current = AttemptState(
                galaxy=galaxy,
                attempt_index=len(attempts) + 1,
                start_line=line_number,
                start_time=timestamp,
                last_line=line_number,
                last_time=timestamp,
            )
            continue
        if current is None:
            continue
        if timestamp is not None:
            current.last_time = timestamp
        current.last_line = line_number

        if "Results of the module are already in the output directory. Module is skipped." in line:
            current.skipped_modules += 1
        if "MUSE_WFM: Reading the MUSE-WFM cube" in line:
            current.read_start = timestamp
        elif "Finished reading the MUSE cube" in line:
            current.read_done = timestamp
            match = re.search(r"Read a total of (\d+) spectra", line)
            if match:
                current.read_total_spectra = int(match.group(1))
        elif "Using full spectral library for ppxf on BIN level" in line:
            current.current_gas_level = "BIN"
            current.bin_level_start = timestamp
        elif "Using full spectral library for ppxf on SPAXEL level" in line:
            current.current_gas_level = "SPAXEL"
            current.spaxel_level_start = timestamp
        elif "Running PPXF for emission lines analysis in parallel mode" in line:
            if current.current_gas_level == "BIN":
                current.bin_fit_start = timestamp
            elif current.current_gas_level == "SPAXEL":
                current.spaxel_fit_start = timestamp
        elif "ppxf_gas_wrapper: Running PPXF on" in line:
            match = re.search(r"on (\d+) spectra took ([0-9.]+)s", line)
            spectra = int(match.group(1)) if match else None
            if current.current_gas_level == "BIN":
                current.bin_ppxf_done = timestamp
                current.bin_ppxf_spectra = spectra
            elif current.current_gas_level == "SPAXEL":
                current.spaxel_ppxf_done = timestamp
                current.spaxel_ppxf_spectra = spectra
        elif f"{galaxy}_gas_bin.fits" in line and "Wrote:" in line:
            current.bin_result_write = timestamp
        elif f"{galaxy}_gas_bestfit_bin.fits" in line and "Wrote:" in line:
            current.bin_bestfit_write = timestamp
        elif f"{galaxy}_gas_spaxel.fits" in line and "Wrote:" in line:
            current.spaxel_result_write = timestamp
        elif f"{galaxy}_gas_bestfit_spaxel.fits" in line and "Wrote:" in line:
            current.spaxel_bestfit_write = timestamp
        elif "_writeFITS: Producing FITS maps from the emission-line analysis" in line:
            current.gas_maps_start = timestamp
        elif "ppxf_sfh_wrapper: Running pPXF in parallel mode" in line:
            current.sfh_start = timestamp
        elif "_writeFITS: Produced SFH maps in FITS format" in line:
            current.sfh_done = timestamp

    finish_current()
    return [attempt.to_record() for attempt in attempts]


def collect_logfile_attempts(root: Path) -> pd.DataFrame:
    records: list[dict] = []
    for logfile in sorted(Path(root).glob("*/LOGFILE")):
        records.extend(parse_logfile_attempts(logfile, logfile.parent.name))
    if not records:
        return pd.DataFrame()
    return pd.DataFrame.from_records(records).sort_values(
        ["galaxy", "attempt_index"], ignore_index=True
    )


def latest_attempt_summary(attempts: pd.DataFrame) -> pd.DataFrame:
    if attempts.empty:
        return pd.DataFrame()
    summary_records: list[dict] = []
    for galaxy, group in attempts.groupby("galaxy", sort=False):
        latest = group.sort_values("attempt_index").iloc[-1]
        completed_gas = group[group["has_gas_maps_start"]]
        max_spax_write = group["spaxel_bestfit_write_hours"].max(skipna=True)
        max_spax_fit = group["spaxel_fit_hours"].max(skipna=True)
        summary_records.append(
            {
                "galaxy": galaxy,
                "latest_attempt_index": int(latest["attempt_index"]),
                "latest_attempt_continuation_only": bool(latest["is_continuation_only"]),
                "latest_attempt_has_gas_maps": bool(latest["has_gas_maps_start"]),
                "latest_attempt_has_sfh_done": bool(latest["has_sfh_done"]),
                "n_log_attempts": int(len(group)),
                "n_gas_map_attempts": int(len(completed_gas)),
                "max_spaxel_bestfit_write_hours": float(max_spax_write)
                if pd.notna(max_spax_write)
                else math.nan,
                "max_spaxel_fit_hours": float(max_spax_fit) if pd.notna(max_spax_fit) else math.nan,
            }
        )
    return pd.DataFrame.from_records(summary_records).set_index("galaxy")


def build_job_dataframe(
    rows: Iterable[dict], metrics: pd.DataFrame, attempts: pd.DataFrame
) -> pd.DataFrame:
    jobs = pd.DataFrame.from_records(list(rows))
    if jobs.empty:
        return jobs
    jobs = jobs.merge(metrics.reset_index(), on="galaxy", how="left")
    attempt_summary = latest_attempt_summary(attempts)
    if not attempt_summary.empty:
        jobs = jobs.merge(attempt_summary.reset_index(), on="galaxy", how="left")
    jobs["memory_measurement"] = np.where(
        jobs["standard_mem_cap"], "censored_at_request", "measured"
    )
    jobs["completed_regular"] = (
        (jobs["state"] == "COMPLETED")
        & (jobs["job_variant"] == "regular")
        & jobs["maxrss_gib"].notna()
    )
    jobs["has_finished_gas_products"] = jobs["has_gas_bin_maps"].fillna(False) & jobs[
        "has_gas_spaxel_maps"
    ].fillna(False)
    jobs["likely_continuation_only_completion"] = (
        jobs["completed_regular"]
        & jobs["latest_attempt_continuation_only"].fillna(False)
        & ~jobs["latest_attempt_has_gas_maps"].fillna(False)
    )
    return jobs


def select_reliable_completed_gas_jobs(jobs: pd.DataFrame) -> pd.DataFrame:
    """Use one completed regular gas-stage job per galaxy for memory scaling."""
    if jobs.empty:
        return jobs
    candidates = jobs[
        jobs["completed_regular"]
        & jobs["has_finished_gas_products"].fillna(False)
        & ~jobs["likely_continuation_only_completion"].fillna(False)
        & ~jobs["standard_mem_cap"].fillna(False)
    ].copy()
    if candidates.empty:
        return candidates
    candidates = candidates.sort_values("jobid", ascending=False)
    return candidates.drop_duplicates("galaxy", keep="first").sort_values(
        "total_cube_spaxels"
    )


def representative_completed_galaxy_peaks(jobs: pd.DataFrame) -> pd.DataFrame:
    """Return one representative peak-memory row for each completed regular galaxy.

    The galaxy set is defined by galaxies with at least one completed regular job.
    The representative memory row is the highest observed regular MaxRSS for that
    galaxy, so resumed SFH-only completions do not hide earlier gas-stage peaks.
    """
    if jobs.empty:
        return jobs

    completed_galaxies = set(
        jobs.loc[
            (jobs["state"] == "COMPLETED")
            & (jobs["job_variant"] == "regular")
            & jobs["total_cube_spaxels"].notna(),
            "galaxy",
        ]
    )
    if not completed_galaxies:
        return jobs.iloc[0:0].copy()

    candidates = jobs[
        jobs["galaxy"].isin(completed_galaxies)
        & (jobs["job_variant"] == "regular")
        & jobs["maxrss_gib"].notna()
        & jobs["total_cube_spaxels"].notna()
        & (jobs["maxrss_gib"] > 1.0)
    ].copy()
    if candidates.empty:
        return candidates

    candidates = candidates.sort_values(["galaxy", "maxrss_gib", "jobid"], ascending=[True, False, False])
    representative = candidates.drop_duplicates("galaxy", keep="first").copy()
    representative["representative_peak_source"] = np.where(
        representative["state"] == "COMPLETED",
        "completed_job",
        "earlier_gas_attempt",
    )
    return representative.sort_values("total_cube_spaxels")


def fit_memory_scalings(reliable_jobs: pd.DataFrame) -> pd.DataFrame:
    variables = ["total_cube_spaxels", "available_spaxels", "n_bins"]
    records: list[dict] = []
    for variable in variables:
        subset = reliable_jobs[[variable, "maxrss_gib"]].dropna()
        x = subset[variable].to_numpy(dtype=float)
        y = subset["maxrss_gib"].to_numpy(dtype=float)
        if len(x) < 2:
            records.append(
                {
                    "metric": variable,
                    "n": len(x),
                    "slope_gib_per_unit": math.nan,
                    "intercept_gib": math.nan,
                    "r": math.nan,
                    "rmse_gib": math.nan,
                }
            )
            continue
        design = np.vstack([x, np.ones_like(x)]).T
        slope, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
        pred = slope * x + intercept
        rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
        records.append(
            {
                "metric": variable,
                "n": int(len(x)),
                "slope_gib_per_unit": float(slope),
                "intercept_gib": float(intercept),
                "r": float(np.corrcoef(x, y)[0, 1]),
                "rmse_gib": rmse,
            }
        )
    return pd.DataFrame.from_records(records).set_index("metric")


def build_queue_risk_table(jobs: pd.DataFrame, metrics: pd.DataFrame, fits: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty or fits.empty:
        return pd.DataFrame()
    total_fit = fits.loc["total_cube_spaxels"]
    slope = total_fit["slope_gib_per_unit"]
    intercept = total_fit["intercept_gib"]
    rmse = total_fit["rmse_gib"]
    records: list[dict] = []
    for galaxy, metric in metrics.iterrows():
        galaxy_jobs = jobs[jobs["galaxy"] == galaxy] if not jobs.empty else pd.DataFrame()
        observed_max = galaxy_jobs["maxrss_gib"].max(skipna=True) if not galaxy_jobs.empty else math.nan
        observed_completed_highmem = bool(
            (
                (galaxy_jobs["state"] == "COMPLETED")
                & (galaxy_jobs["reqmem_gib"] >= 900)
            ).any()
        ) if not galaxy_jobs.empty else False
        capped_standard = bool(galaxy_jobs["standard_mem_cap"].any()) if not galaxy_jobs.empty else False
        completed_standard = bool(
            (
                (galaxy_jobs["state"] == "COMPLETED")
                & (galaxy_jobs["reqmem_gib"] <= 230)
                & ~galaxy_jobs["standard_mem_cap"]
            ).any()
        ) if not galaxy_jobs.empty else False
        predicted = slope * metric["total_cube_spaxels"] + intercept
        upper = predicted + rmse
        if observed_completed_highmem:
            recommendation = "highmem confirmed"
        elif capped_standard and not completed_standard:
            recommendation = "highmem likely"
        elif upper >= 220:
            recommendation = "borderline or highmem"
        else:
            recommendation = "220GiB likely enough"
        records.append(
            {
                "galaxy": galaxy,
                "total_cube_spaxels": int(metric["total_cube_spaxels"]),
                "available_spaxels": int(metric["available_spaxels"]),
                "n_bins": int(metric["n_bins"]),
                "observed_maxrss_gib": float(observed_max) if pd.notna(observed_max) else math.nan,
                "standard_220gib_capped": capped_standard,
                "completed_standard_220gib": completed_standard,
                "completed_highmem": observed_completed_highmem,
                "predicted_maxrss_gib": float(predicted),
                "prediction_rmse_gib": float(rmse),
                "predicted_plus_rmse_gib": float(upper),
                "queue_risk": recommendation,
            }
        )
    return pd.DataFrame.from_records(records).sort_values(
        ["predicted_maxrss_gib", "total_cube_spaxels"], ascending=False
    )


def load_analysis(
    efficiency_path: Path = Path("finished_efficiency_20260528_111829.txt"),
    product_root: Path = Path("v3tk_v7.6.8"),
) -> dict[str, pd.DataFrame]:
    rows = parse_efficiency_table(efficiency_path)
    metrics = collect_cube_metrics(product_root)
    attempts = collect_logfile_attempts(product_root)
    jobs = build_job_dataframe(rows, metrics, attempts)
    reliable = select_reliable_completed_gas_jobs(jobs)
    fits = fit_memory_scalings(reliable)
    risk = build_queue_risk_table(jobs, metrics, fits)
    return {
        "efficiency_rows": pd.DataFrame.from_records(rows),
        "cube_metrics": metrics,
        "log_attempts": attempts,
        "jobs": jobs,
        "reliable_jobs": reliable,
        "memory_fits": fits,
        "queue_risk": risk,
    }


def add_memory_scaling_panels(axs, jobs: pd.DataFrame, reliable: pd.DataFrame, fits: pd.DataFrame) -> None:
    variables = [
        ("total_cube_spaxels", "Total cube spaxels"),
        ("available_spaxels", "Available spaxels"),
        ("n_bins", "Voronoi bins"),
    ]
    state_styles = {
        "COMPLETED": dict(marker="o", color="#2a9d8f"),
        "TIMEOUT": dict(marker="^", color="#e76f51"),
        "FAILED": dict(marker="x", color="#d62828"),
        "CANCELLED": dict(marker="s", color="#6c757d"),
    }
    for ax, (variable, label) in zip(axs, variables):
        for state, group in jobs.dropna(subset=[variable, "maxrss_gib"]).groupby("state"):
            style = state_styles.get(state, dict(marker=".", color="#555555"))
            ax.scatter(
                group[variable],
                group["maxrss_gib"],
                s=np.where(group["reqmem_gib"] >= 900, 80, 44),
                linewidth=1.0,
                facecolors="none" if state == "TIMEOUT" else style["color"],
                edgecolors=style["color"],
                marker=style["marker"],
                label=state,
                alpha=0.85,
            )
        capped = jobs[jobs["standard_mem_cap"].fillna(False) & jobs[variable].notna()]
        if not capped.empty:
            ax.scatter(
                capped[variable],
                capped["maxrss_gib"],
                marker="_",
                color="black",
                s=180,
                linewidth=2,
                label="220GiB cap",
            )
        if variable in fits.index and pd.notna(fits.loc[variable, "slope_gib_per_unit"]):
            x = reliable[variable].dropna().to_numpy(dtype=float)
            if len(x):
                xx = np.linspace(0, max(x.max(), jobs[variable].max()), 100)
                yy = fits.loc[variable, "slope_gib_per_unit"] * xx + fits.loc[variable, "intercept_gib"]
                ax.plot(xx, yy, color="#264653", lw=2, label=f"fit r={fits.loc[variable, 'r']:.2f}")
        ax.axhline(220, color="#222222", lw=1, ls="--")
        ax.axhline(980, color="#6c757d", lw=1, ls=":")
        ax.set_xlabel(label)
        ax.set_ylabel("MaxRSS [GiB]")
        ax.grid(alpha=0.25)
    handles, labels = axs[0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    axs[-1].legend(unique.values(), unique.keys(), loc="best", fontsize=8)


def case_study_table(attempts: pd.DataFrame, galaxies: Iterable[str]) -> pd.DataFrame:
    columns = [
        "galaxy",
        "attempt_index",
        "read_total_spectra",
        "wall_hours",
        "bin_fit_hours",
        "bin_bestfit_write_hours",
        "spaxel_fit_hours",
        "spaxel_bestfit_write_hours",
        "bestfit_to_maps_hours",
        "sfh_hours",
        "has_gas_maps_start",
        "has_sfh_done",
    ]
    subset = attempts[attempts["galaxy"].isin(list(galaxies))].copy()
    if subset.empty:
        return subset
    interesting = subset[
        subset["has_gas_maps_start"] | subset["has_spaxel_fit"] | subset["has_sfh_done"]
    ]
    return interesting[columns].sort_values(["galaxy", "attempt_index"])
