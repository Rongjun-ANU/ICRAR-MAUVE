from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_module():
    try:
        import ngist_memory_runtime_analysis as mod
    except ModuleNotFoundError as exc:
        raise AssertionError("ngist_memory_runtime_analysis module is missing") from exc
    return mod


def test_efficiency_parser_uses_recent_rows_and_maxrss():
    mod = load_module()
    rows = mod.parse_efficiency_table(ROOT / "finished_efficiency_20260528_111829.txt")

    assert len(rows) == 65
    assert rows[-1]["jobid"] == 43057897
    assert all(row["jobid"] >= 43057897 for row in rows)

    ic_both = next(row for row in rows if row["jobid"] == 43057897)
    ic_spaxel = next(row for row in rows if row["jobid"] == 43300739)
    ic_bin = next(row for row in rows if row["jobid"] == 43330802)

    assert ic_both["galaxy"] == "IC3392"
    assert ic_both["maxrss_gib"] == 113.0
    assert ic_both["averss_gib"] == 56.6
    assert ic_spaxel["maxrss_gib"] == 125.0
    assert ic_bin["maxrss_gib"] == 63.1


def test_cube_metrics_match_logfile_read_count_and_mask_for_ic3392():
    mod = load_module()
    metrics = mod.collect_cube_metrics(ROOT / "v3tk_v7.6.8")
    ic = metrics.loc["IC3392"]

    assert ic["shape"] == "438x437"
    assert int(ic["total_cube_spaxels"]) == 191406
    assert int(ic["available_spaxels"]) == 94575
    assert int(ic["mask_zero_spaxels"]) == int(ic["available_spaxels"])
    assert int(ic["n_bins"]) == 4062
    assert bool(ic["has_gas_bin_maps"])
    assert bool(ic["has_gas_spaxel_maps"])


def test_logfile_parser_captures_spaxel_bestfit_write_intervals():
    mod = load_module()
    attempts = mod.parse_logfile_attempts(ROOT / "v3tk_v7.6.8" / "NGC4192" / "LOGFILE", "NGC4192")
    finished = [attempt for attempt in attempts if attempt.get("gas_maps_start") is not None]

    assert len(finished) == 1
    attempt = finished[0]
    assert int(attempt["read_total_spectra"]) == 1778868
    assert 44.0 < attempt["spaxel_bestfit_write_hours"] < 45.0
    assert 0.65 < attempt["bin_fit_hours"] < 0.67
    assert 14.0 < attempt["spaxel_fit_hours"] < 14.1


def test_scaling_fit_prefers_total_cube_size_over_available_spaxels():
    mod = load_module()
    rows = mod.parse_efficiency_table(ROOT / "finished_efficiency_20260528_111829.txt")
    metrics = mod.collect_cube_metrics(ROOT / "v3tk_v7.6.8")
    attempts = mod.collect_logfile_attempts(ROOT / "v3tk_v7.6.8")
    jobs = mod.build_job_dataframe(rows, metrics, attempts)
    reliable = mod.select_reliable_completed_gas_jobs(jobs)
    fits = mod.fit_memory_scalings(reliable)

    assert fits.loc["total_cube_spaxels", "r"] > fits.loc["available_spaxels", "r"]
    assert fits.loc["total_cube_spaxels", "r"] > 0.95
    assert "NGC4580" not in set(reliable["galaxy"])


def test_representative_completed_peak_uses_ngc4580_gas_attempt():
    mod = load_module()
    rows = mod.parse_efficiency_table(ROOT / "finished_efficiency_20260528_111829.txt")
    metrics = mod.collect_cube_metrics(ROOT / "v3tk_v7.6.8")
    attempts = mod.collect_logfile_attempts(ROOT / "v3tk_v7.6.8")
    jobs = mod.build_job_dataframe(rows, metrics, attempts)
    representative = mod.representative_completed_galaxy_peaks(jobs)

    assert len(representative) == 18
    ngc4580 = representative.loc[representative["galaxy"] == "NGC4580"].iloc[0]
    assert int(ngc4580["jobid"]) == 43091812
    assert ngc4580["state"] == "TIMEOUT"
    assert ngc4580["maxrss_gib"] == 183.0
    assert ngc4580["representative_peak_source"] == "earlier_gas_attempt"


if __name__ == "__main__":
    for name, func in sorted(globals().items()):
        if name.startswith("test_") and callable(func):
            func()
            print(f"PASS {name}")
