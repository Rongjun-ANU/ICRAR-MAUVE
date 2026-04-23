from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr


LOG10 = np.log(10.0)


def _first_existing_hdu(hdul, candidates):
    available = {hdu.name for hdu in hdul}
    for name in candidates:
        if name in available:
            return hdul[name].data
    raise KeyError(f"None of the candidate HDUs were found: {candidates}")


def _propagate_log10_ratio_error(numerator, denominator, numerator_err, denominator_err):
    error = np.full_like(numerator, np.nan, dtype=float)
    valid = (
        np.isfinite(numerator)
        & np.isfinite(denominator)
        & np.isfinite(numerator_err)
        & np.isfinite(denominator_err)
        & (numerator > 0)
        & (denominator > 0)
        & (numerator_err >= 0)
        & (denominator_err >= 0)
    )
    if not np.any(valid):
        return error

    error[valid] = (
        np.sqrt(
            (numerator_err[valid] / numerator[valid]) ** 2
            + (denominator_err[valid] / denominator[valid]) ** 2
        )
        / LOG10
    )
    return error


def _propagate_o3n2_m13_error(oiii, hb, nii, ha, oiii_err, hb_err, nii_err, ha_err):
    o3_err = _propagate_log10_ratio_error(oiii, hb, oiii_err, hb_err)
    n2_err = _propagate_log10_ratio_error(nii, ha, nii_err, ha_err)

    total = np.full_like(o3_err, np.nan, dtype=float)
    valid = np.isfinite(o3_err) & np.isfinite(n2_err)
    if not np.any(valid):
        return total

    # Marino et al. (2013): 12 + log(O/H) = 8.533 - 0.214 * O3N2
    total[valid] = 0.214 * np.sqrt(o3_err[valid] ** 2 + n2_err[valid] ** 2)
    return total


def _in_bin_mask(values, left, right, is_last):
    if is_last:
        return (values >= left) & (values <= right)
    return (values >= left) & (values < right)


def _fill_nan_by_interpolation(values):
    filled = np.asarray(values, dtype=float).copy()
    finite = np.isfinite(filled)
    if np.all(finite):
        return filled
    if not np.any(finite):
        return np.zeros_like(filled)

    indices = np.arange(filled.size)
    filled[~finite] = np.interp(indices[~finite], indices[finite], filled[finite])
    return filled


def _digitize_to_bins(values, bin_edges):
    idx = np.digitize(values, bin_edges[1:-1], right=False)
    return np.clip(idx, 0, len(bin_edges) - 2)


def make_percentile_mass_bin_edges(
    values,
    n_bins=6,
    lower_percentile=2.5,
    upper_percentile=97.5,
):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise RuntimeError("Cannot define mass-bin edges from an empty sample.")

    vmin = np.nanpercentile(finite, lower_percentile)
    vmax = np.nanpercentile(finite, upper_percentile)
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        raise RuntimeError("Invalid mass range for percentile-based bin construction.")

    return np.linspace(vmin, vmax, n_bins + 1)


def calculate_binned_stats(x_data, y_data, bin_width=0.2, min_unique=20):
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)

    finite = np.isfinite(x_data) & np.isfinite(y_data)
    if np.sum(finite) == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])

    x_min = np.nanmin(x_data[finite])
    x_max = np.nanmax(x_data[finite])
    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min == x_max:
        return np.array([]), np.array([]), np.array([]), np.array([])

    bin_edges = np.arange(x_min, x_max + bin_width, bin_width)
    if bin_edges.size < 2:
        return np.array([]), np.array([]), np.array([]), np.array([])

    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    out_centers = []
    out_medians = []
    out_stds = []
    out_counts = []

    for i, (left, right) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
        in_bin = finite & (x_data >= left) & (x_data < right)
        if np.sum(in_bin) == 0:
            continue

        y_bin = y_data[in_bin]
        if len(np.unique(y_bin)) < min_unique:
            continue

        out_centers.append(centers[i])
        out_medians.append(np.nanmedian(y_bin))
        out_stds.append(np.nanstd(y_bin))
        out_counts.append(np.sum(in_bin))

    return (
        np.asarray(out_centers),
        np.asarray(out_medians),
        np.asarray(out_stds),
        np.asarray(out_counts),
    )


def load_o3n2_m13_hii_sample_with_errors(excluded_galaxies=("NGC4383",)):
    from astropy.io import fits

    excluded = set(excluded_galaxies)
    spatial_files = sorted(Path(".").glob("*_SPATIAL_BINNING_maps_extended.fits"))

    log_sigma_star = []
    log_sigma_sfr = []
    oh = []
    sigma_log_sfr = []
    sigma_oh = []
    galaxy_name = []

    for spatial_path in spatial_files:
        galaxy = spatial_path.name.split("_")[0]
        if galaxy in excluded:
            continue

        gas_path = Path(f"{galaxy}_gas_BIN_maps_extended.fits")
        if not gas_path.exists():
            continue

        with fits.open(spatial_path) as h_spatial, fits.open(gas_path) as h_gas:
            mass_map = h_spatial["LOGMASS_SURFACE_DENSITY"].data
            sfr_map = h_gas["LOGSFR_SURFACE_DENSITY_HII"].data
            oh_map = h_gas["O_H_O3N2_M13_HII"].data

            ha_corr = h_gas["HA6562_FLUX_corr_HII"].data
            hb_corr = h_gas["HB4861_FLUX_corr_HII"].data
            oiii_corr = h_gas["OIII5006_FLUX_corr_HII"].data
            nii_corr = h_gas["NII6583_FLUX_corr_HII"].data

            ha_err = _first_existing_hdu(
                h_gas,
                (
                    "HA6563_FLUX_ERR",
                    "HA6562_FLUX_ERR",
                    "HALPHA6563_FLUX_ERR",
                    "HALPHA_FLUX_ERR",
                ),
            )
            hb_err = h_gas["HB4861_FLUX_ERR"].data
            oiii_err = h_gas["OIII5006_FLUX_ERR"].data
            nii_err = _first_existing_hdu(h_gas, ("NII6584_FLUX_ERR", "NII6583_FLUX_ERR"))

        common_mask = (
            np.isfinite(mass_map)
            & np.isfinite(sfr_map)
            & np.isfinite(oh_map)
        )

        sfr_err_map = np.full_like(sfr_map, np.nan, dtype=float)
        valid_ha = np.isfinite(ha_corr) & np.isfinite(ha_err) & (ha_corr > 0) & (ha_err >= 0)
        sfr_err_map[valid_ha] = ha_err[valid_ha] / (ha_corr[valid_ha] * LOG10)

        oh_err_map = _propagate_o3n2_m13_error(
            oiii_corr,
            hb_corr,
            nii_corr,
            ha_corr,
            oiii_err,
            hb_err,
            nii_err,
            ha_err,
        )

        n_common = int(np.sum(common_mask))
        log_sigma_star.append(mass_map[common_mask].ravel())
        log_sigma_sfr.append(sfr_map[common_mask].ravel())
        oh.append(oh_map[common_mask].ravel())
        sigma_log_sfr.append(sfr_err_map[common_mask].ravel())
        sigma_oh.append(oh_err_map[common_mask].ravel())
        galaxy_name.append(np.full(n_common, galaxy, dtype=object))

    if not log_sigma_star:
        raise RuntimeError("No valid HII sample could be loaded.")

    return {
        "log_sigma_star": np.concatenate(log_sigma_star),
        "log_sigma_sfr": np.concatenate(log_sigma_sfr),
        "oh": np.concatenate(oh),
        "sigma_log_sfr": np.concatenate(sigma_log_sfr),
        "sigma_oh": np.concatenate(sigma_oh),
        "galaxy_name": np.concatenate(galaxy_name),
    }


def mask_sample(sample, mask):
    masked = {}
    for key, value in sample.items():
        if isinstance(value, np.ndarray) and value.shape == mask.shape:
            masked[key] = value[mask]
        else:
            masked[key] = value
    return masked


def compute_rho_curve_from_pipeline(
    log_sigma_star,
    log_sigma_sfr,
    oh,
    mass_bin_edges,
    bin_width=0.2,
    min_unique=20,
    min_trend_points=3,
):
    x = np.asarray(log_sigma_star)
    y = np.asarray(log_sigma_sfr)
    z = np.asarray(oh)

    n_bins = len(mass_bin_edges) - 1
    centers = 0.5 * (mass_bin_edges[:-1] + mass_bin_edges[1:])
    rho = np.full(n_bins, np.nan, dtype=float)
    p_value = np.full(n_bins, np.nan, dtype=float)
    counts = np.zeros(n_bins, dtype=int)
    clipped_counts = np.zeros(n_bins, dtype=int)

    delta_sfr_cloud = []
    delta_oh_cloud = []
    mass_bin_index_cloud = []

    for i in range(n_bins):
        left = mass_bin_edges[i]
        right = mass_bin_edges[i + 1]
        in_bin = _in_bin_mask(x, left, right, i == n_bins - 1)
        counts[i] = int(np.sum(in_bin))
        if counts[i] == 0:
            continue

        y_bin = y[in_bin]
        z_bin = z[in_bin]

        y_offset = y_bin - np.nanmean(y_bin)
        z_offset = z_bin - np.nanmean(z_bin)
        valid = np.isfinite(y_offset) & np.isfinite(z_offset)
        if np.sum(valid) <= 2:
            continue

        sigma_clip_mask = np.zeros_like(valid, dtype=bool)
        y_valid = y_offset[valid]
        z_valid = z_offset[valid]

        y_unique = np.unique(y_valid)
        z_unique = np.unique(z_valid)
        x_mean = np.mean(y_unique)
        x_std = np.std(y_unique)
        y_mean = np.mean(z_unique)
        y_std = np.std(z_unique)

        x_lower = x_mean - 3.0 * x_std
        x_upper = x_mean + 3.0 * x_std
        y_lower = y_mean - 3.0 * y_std
        y_upper = y_mean + 3.0 * y_std

        keep_valid = (
            (y_valid >= x_lower)
            & (y_valid <= x_upper)
            & (z_valid >= y_lower)
            & (z_valid <= y_upper)
        )
        sigma_clip_mask[np.where(valid)[0]] = keep_valid
        clipped_counts[i] = int(np.sum(sigma_clip_mask))
        if clipped_counts[i] <= 2:
            continue

        trend_x = y_offset[sigma_clip_mask]
        trend_y = z_offset[sigma_clip_mask]
        delta_sfr_cloud.append(trend_x)
        delta_oh_cloud.append(trend_y)
        mass_bin_index_cloud.append(np.full(trend_x.size, i, dtype=int))

        trend_centers, _, _, _ = calculate_binned_stats(
            trend_x,
            trend_y,
            bin_width=bin_width,
            min_unique=min_unique,
        )
        if trend_centers.size < min_trend_points:
            continue

        rho[i], p_value[i] = spearmanr(trend_x, trend_y)

    if delta_sfr_cloud:
        delta_sfr_cloud = np.concatenate(delta_sfr_cloud)
        delta_oh_cloud = np.concatenate(delta_oh_cloud)
        mass_bin_index_cloud = np.concatenate(mass_bin_index_cloud)
    else:
        delta_sfr_cloud = np.array([])
        delta_oh_cloud = np.array([])
        mass_bin_index_cloud = np.array([], dtype=int)

    return {
        "centers": centers,
        "rho": rho,
        "p_value": p_value,
        "counts": counts,
        "clipped_counts": clipped_counts,
        "delta_sfr": delta_sfr_cloud,
        "delta_oh": delta_oh_cloud,
        "mass_bin_index": mass_bin_index_cloud,
    }


def _build_primary_bin_model(sample, mass_bin_edges):
    x = sample["log_sigma_star"]
    y = sample["log_sigma_sfr"]
    z = sample["oh"]
    sigma_y = sample["sigma_log_sfr"]
    sigma_z = sample["sigma_oh"]

    n_bins = len(mass_bin_edges) - 1
    mean_sfr = np.full(n_bins, np.nan, dtype=float)
    mean_oh = np.full(n_bins, np.nan, dtype=float)
    intrinsic_sfr = np.full(n_bins, np.nan, dtype=float)
    intrinsic_oh = np.full(n_bins, np.nan, dtype=float)

    for i in range(n_bins):
        left = mass_bin_edges[i]
        right = mass_bin_edges[i + 1]
        in_bin = _in_bin_mask(x, left, right, i == n_bins - 1)
        if np.sum(in_bin) == 0:
            continue

        y_bin = y[in_bin]
        z_bin = z[in_bin]
        sigma_y_bin = sigma_y[in_bin]
        sigma_z_bin = sigma_z[in_bin]

        mean_sfr[i] = np.nanmean(y_bin)
        mean_oh[i] = np.nanmean(z_bin)

        y_std = np.nanstd(y_bin)
        z_std = np.nanstd(z_bin)
        y_err_rms = np.sqrt(np.nanmean(sigma_y_bin ** 2))
        z_err_rms = np.sqrt(np.nanmean(sigma_z_bin ** 2))

        intrinsic_sfr[i] = np.sqrt(max(y_std ** 2 - y_err_rms ** 2, 0.0))
        intrinsic_oh[i] = np.sqrt(max(z_std ** 2 - z_err_rms ** 2, 0.0))

    return {
        "mean_sfr": _fill_nan_by_interpolation(mean_sfr),
        "mean_oh": _fill_nan_by_interpolation(mean_oh),
        "intrinsic_sfr": _fill_nan_by_interpolation(intrinsic_sfr),
        "intrinsic_oh": _fill_nan_by_interpolation(intrinsic_oh),
    }


def _mass_error_proxy(log_sigma_star, base_error=0.01):
    x = np.asarray(log_sigma_star, dtype=float)
    x_min = np.nanmin(x)
    x_max = np.nanmax(x)
    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min == x_max:
        return np.full_like(x, base_error)

    return np.abs(2.0 * base_error - base_error * (x - x_min) / (x_max - x_min))


def run_sec34_style_null_test(
    sample,
    mass_bin_edges,
    n_realizations=250,
    random_seed=726,
    include_mass_error_proxy=True,
    mass_error_base=0.01,
    bin_width=0.2,
    min_unique=20,
    min_trend_points=3,
):
    x_true = np.asarray(sample["log_sigma_star"], dtype=float)
    y_obs = np.asarray(sample["log_sigma_sfr"], dtype=float)
    z_obs = np.asarray(sample["oh"], dtype=float)
    sigma_y = np.asarray(sample["sigma_log_sfr"], dtype=float)
    sigma_z = np.asarray(sample["sigma_oh"], dtype=float)

    finite = (
        np.isfinite(x_true)
        & np.isfinite(y_obs)
        & np.isfinite(z_obs)
        & np.isfinite(sigma_y)
        & np.isfinite(sigma_z)
    )
    if not np.all(finite):
        sample = mask_sample(sample, finite)
        x_true = sample["log_sigma_star"]
        y_obs = sample["log_sigma_sfr"]
        z_obs = sample["oh"]
        sigma_y = sample["sigma_log_sfr"]
        sigma_z = sample["sigma_oh"]

    observed = compute_rho_curve_from_pipeline(
        x_true,
        y_obs,
        z_obs,
        mass_bin_edges,
        bin_width=bin_width,
        min_unique=min_unique,
        min_trend_points=min_trend_points,
    )

    model = _build_primary_bin_model(sample, mass_bin_edges)
    bin_index = _digitize_to_bins(x_true, mass_bin_edges)
    n_bins = len(mass_bin_edges) - 1
    rho_matrix = np.full((n_realizations, n_bins), np.nan, dtype=float)

    sigma_x = None
    if include_mass_error_proxy:
        sigma_x = _mass_error_proxy(x_true, base_error=mass_error_base)

    rng = np.random.default_rng(random_seed)
    example_cloud = None

    for i in range(n_realizations):
        y_true = model["mean_sfr"][bin_index] + rng.normal(
            loc=0.0,
            scale=model["intrinsic_sfr"][bin_index],
        )
        z_true = model["mean_oh"][bin_index] + rng.normal(
            loc=0.0,
            scale=model["intrinsic_oh"][bin_index],
        )

        if include_mass_error_proxy:
            x_mock = x_true + rng.normal(loc=0.0, scale=sigma_x)
        else:
            x_mock = x_true.copy()

        y_mock = y_true + rng.normal(loc=0.0, scale=sigma_y)
        z_mock = z_true + rng.normal(loc=0.0, scale=sigma_z)

        mock_curve = compute_rho_curve_from_pipeline(
            x_mock,
            y_mock,
            z_mock,
            mass_bin_edges,
            bin_width=bin_width,
            min_unique=min_unique,
            min_trend_points=min_trend_points,
        )
        rho_matrix[i] = mock_curve["rho"]
        if example_cloud is None:
            example_cloud = mock_curve

    return {
        "sample_size": x_true.size,
        "mass_bin_edges": np.asarray(mass_bin_edges, dtype=float),
        "observed": observed,
        "primary_model": model,
        "null": {
            "rho_matrix": rho_matrix,
            "median_rho": np.nanmedian(rho_matrix, axis=0),
            "p16_rho": np.nanpercentile(rho_matrix, 16.0, axis=0),
            "p84_rho": np.nanpercentile(rho_matrix, 84.0, axis=0),
            "example_cloud": example_cloud,
        },
        "config": {
            "n_realizations": n_realizations,
            "random_seed": random_seed,
            "include_mass_error_proxy": include_mass_error_proxy,
            "mass_error_base": mass_error_base,
            "bin_width": bin_width,
            "min_unique": min_unique,
            "min_trend_points": min_trend_points,
        },
    }


def print_sec34_summary(result):
    edges = result["mass_bin_edges"]
    centers = result["observed"]["centers"]
    obs_rho = result["observed"]["rho"]
    obs_n = result["observed"]["counts"]
    null_p16 = result["null"]["p16_rho"]
    null_med = result["null"]["median_rho"]
    null_p84 = result["null"]["p84_rho"]

    print("=" * 132)
    print("Sec. 3.4-style null test summary")
    print("=" * 132)
    print(
        f"{'Bin':>3}  {'Mass range':>15}  {'center':>8}  {'N':>8}  "
        f"{'obs rho':>10}  {'null p16':>10}  {'null med':>10}  {'null p84':>10}"
    )
    for i in range(len(centers)):
        obs_txt = f"{obs_rho[i]:.4f}" if np.isfinite(obs_rho[i]) else "N/A"
        p16_txt = f"{null_p16[i]:.4f}" if np.isfinite(null_p16[i]) else "N/A"
        med_txt = f"{null_med[i]:.4f}" if np.isfinite(null_med[i]) else "N/A"
        p84_txt = f"{null_p84[i]:.4f}" if np.isfinite(null_p84[i]) else "N/A"
        print(
            f"{i+1:3d}  {edges[i]:6.2f}-{edges[i+1]:6.2f}  {centers[i]:8.3f}  {obs_n[i]:8d}  "
            f"{obs_txt:>10}  {p16_txt:>10}  {med_txt:>10}  {p84_txt:>10}"
        )
    print("=" * 132)


def plot_sec34_style_null_test(result, reference_curve=None, panel_mode="full"):
    observed = result["observed"]
    mock = result["null"]["example_cloud"]
    centers = observed["centers"]

    if panel_mode not in {"full", "rho_only"}:
        raise ValueError(f"Unsupported panel_mode={panel_mode!r}.")

    all_dx = [observed["delta_sfr"]]
    all_dz = [observed["delta_oh"]]
    if mock is not None:
        all_dx.append(mock["delta_sfr"])
        all_dz.append(mock["delta_oh"])

    x_lim = np.nanpercentile(np.abs(np.concatenate(all_dx)), 99.0)
    y_lim = np.nanpercentile(np.abs(np.concatenate(all_dz)), 99.0)
    x_lim = max(x_lim, 0.1)
    y_lim = max(y_lim, 0.05)

    if panel_mode == "rho_only":
        fig, rho_ax = plt.subplots(1, 1, figsize=(6.6, 5.2))
    else:
        fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))

        hb1 = axes[0].hexbin(
            observed["delta_sfr"],
            observed["delta_oh"],
            gridsize=55,
            mincnt=1,
            cmap="Greys",
            linewidths=0.0,
        )
        axes[0].axhline(0.0, color="gray", linestyle="--", linewidth=1.0, alpha=0.6)
        axes[0].axvline(0.0, color="gray", linestyle="--", linewidth=1.0, alpha=0.6)
        axes[0].set_xlim(-x_lim, x_lim)
        axes[0].set_ylim(-y_lim, y_lim)
        axes[0].set_title("Observed residual cloud", fontsize=12)
        axes[0].set_xlabel(r"$\Delta \log\,\Sigma_{\rm SFR}$", fontsize=11)
        axes[0].set_ylabel(r"$\Delta$[O/H]", fontsize=11)
        axes[0].grid(True, alpha=0.2)
        fig.colorbar(hb1, ax=axes[0], label="N")

        if mock is not None:
            hb2 = axes[1].hexbin(
                mock["delta_sfr"],
                mock["delta_oh"],
                gridsize=55,
                mincnt=1,
                cmap="Oranges",
                linewidths=0.0,
            )
            fig.colorbar(hb2, ax=axes[1], label="N")
        axes[1].axhline(0.0, color="gray", linestyle="--", linewidth=1.0, alpha=0.6)
        axes[1].axvline(0.0, color="gray", linestyle="--", linewidth=1.0, alpha=0.6)
        axes[1].set_xlim(-x_lim, x_lim)
        axes[1].set_ylim(-y_lim, y_lim)
        axes[1].set_title("Null cloud (one realization)", fontsize=12)
        axes[1].set_xlabel(r"$\Delta \log\,\Sigma_{\rm SFR}$", fontsize=11)
        axes[1].set_ylabel(r"$\Delta$[O/H]", fontsize=11)
        axes[1].grid(True, alpha=0.2)

        rho_ax = axes[2]

    rho_ax.fill_between(
        centers,
        result["null"]["p16_rho"],
        result["null"]["p84_rho"],
        color="tab:orange",
        alpha=0.25,
        label="Null 16-84%",
    )
    rho_ax.plot(
        centers,
        result["null"]["median_rho"],
        color="tab:orange",
        linewidth=2.0,
        label="Null median",
    )
    rho_ax.plot(
        centers,
        observed["rho"],
        color="black",
        linewidth=2.2,
        marker="o",
        label="Observed (error-qualified subset)",
    )
    if reference_curve is not None:
        ref_centers = np.asarray(reference_curve["centers"], dtype=float)
        ref_rho = np.asarray(reference_curve["rho"], dtype=float)
        ref_valid = np.isfinite(ref_centers) & np.isfinite(ref_rho)
        if np.any(ref_valid):
            rho_ax.plot(
                ref_centers[ref_valid],
                ref_rho[ref_valid],
                color="tab:blue",
                linewidth=1.8,
                linestyle="--",
                marker="s",
                markersize=4,
                label=reference_curve.get("label", "Reference observed curve"),
            )
    rho_ax.axhline(0.0, color="gray", linestyle="--", linewidth=1.0, alpha=0.6)
    rho_ax.set_xlabel(r"$\log\,\Sigma_*$ bin center", fontsize=11)
    rho_ax.set_ylabel(r"Spearman $\rho$", fontsize=11)
    rho_ax.set_title(r"Observed vs null $\rho(\log\,\Sigma_*)$", fontsize=12)
    rho_ax.grid(True, alpha=0.25)
    rho_ax.legend(loc="best", fontsize=10)

    if result["config"]["include_mass_error_proxy"]:
        mass_label = (
            "Includes Sanchez-style mass-error proxy "
            f"(base={result['config']['mass_error_base']:.3f} dex)"
        )
    else:
        mass_label = "Mass distribution fixed; no extra mass-error proxy"

    if panel_mode == "rho_only":
        fig.suptitle(
            "Null test: observed vs null Spearman $\\rho$",
            fontsize=13,
        )
        plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    else:
        fig.suptitle(
            "Null test: independent primary relations + propagated errors\n"
            # + mass_label
            ,
            fontsize=13,
        )
        plt.tight_layout()
    plt.show()
    return fig


def list_available_galaxies():
    return [path.name.split("_")[0] for path in sorted(Path(".").glob("*_SPATIAL_BINNING_maps_extended.fits"))]


def load_o3n2_m13_hii_galaxy_maps_with_errors(galaxy):
    from astropy.io import fits
    from astropy.wcs import WCS

    spatial_path = Path(f"{galaxy}_SPATIAL_BINNING_maps_extended.fits")
    gas_path = Path(f"{galaxy}_gas_BIN_maps_extended.fits")
    if not (spatial_path.exists() and gas_path.exists()):
        raise FileNotFoundError(f"Required FITS files for {galaxy} not found.")

    with fits.open(spatial_path) as h_spatial, fits.open(gas_path) as h_gas:
        mass_map = h_spatial["LOGMASS_SURFACE_DENSITY"].data
        sfr_map = h_gas["LOGSFR_SURFACE_DENSITY_HII"].data
        oh_map = h_gas["O_H_O3N2_M13_HII"].data
        gas_header = h_gas["LOGSFR_SURFACE_DENSITY_HII"].header

        ha_corr = h_gas["HA6562_FLUX_corr_HII"].data
        hb_corr = h_gas["HB4861_FLUX_corr_HII"].data
        oiii_corr = h_gas["OIII5006_FLUX_corr_HII"].data
        nii_corr = h_gas["NII6583_FLUX_corr_HII"].data

        ha_err = _first_existing_hdu(
            h_gas,
            (
                "HA6563_FLUX_ERR",
                "HA6562_FLUX_ERR",
                "HALPHA6563_FLUX_ERR",
                "HALPHA_FLUX_ERR",
            ),
        )
        hb_err = h_gas["HB4861_FLUX_ERR"].data
        oiii_err = h_gas["OIII5006_FLUX_ERR"].data
        nii_err = _first_existing_hdu(h_gas, ("NII6584_FLUX_ERR", "NII6583_FLUX_ERR"))

    common_mask = np.isfinite(mass_map) & np.isfinite(sfr_map) & np.isfinite(oh_map)

    sfr_err_map = np.full_like(sfr_map, np.nan, dtype=float)
    valid_ha = np.isfinite(ha_corr) & np.isfinite(ha_err) & (ha_corr > 0) & (ha_err >= 0)
    sfr_err_map[valid_ha] = ha_err[valid_ha] / (ha_corr[valid_ha] * LOG10)

    oh_err_map = _propagate_o3n2_m13_error(
        oiii_corr,
        hb_corr,
        nii_corr,
        ha_corr,
        oiii_err,
        hb_err,
        nii_err,
        ha_err,
    )

    valid_mask = common_mask & np.isfinite(sfr_err_map) & np.isfinite(oh_err_map)

    return {
        "galaxy": galaxy,
        "log_sigma_star_map": np.asarray(mass_map, dtype=float),
        "log_sigma_sfr_map": np.asarray(sfr_map, dtype=float),
        "oh_map": np.asarray(oh_map, dtype=float),
        "sigma_log_sfr_map": np.asarray(sfr_err_map, dtype=float),
        "sigma_oh_map": np.asarray(oh_err_map, dtype=float),
        "common_mask": np.asarray(common_mask, dtype=bool),
        "valid_mask": np.asarray(valid_mask, dtype=bool),
        "gas_header": gas_header,
        "wcs_celestial": WCS(gas_header).celestial,
    }


def galaxy_maps_to_sample(galaxy_maps):
    valid_mask = np.asarray(galaxy_maps["valid_mask"], dtype=bool)
    galaxy = galaxy_maps["galaxy"]
    n_valid = int(np.sum(valid_mask))
    if n_valid == 0:
        raise RuntimeError(f"No finite propagated-error spaxels are available for {galaxy}.")

    return {
        "log_sigma_star": galaxy_maps["log_sigma_star_map"][valid_mask].ravel(),
        "log_sigma_sfr": galaxy_maps["log_sigma_sfr_map"][valid_mask].ravel(),
        "oh": galaxy_maps["oh_map"][valid_mask].ravel(),
        "sigma_log_sfr": galaxy_maps["sigma_log_sfr_map"][valid_mask].ravel(),
        "sigma_oh": galaxy_maps["sigma_oh_map"][valid_mask].ravel(),
        "galaxy_name": np.full(n_valid, galaxy, dtype=object),
    }


def generate_sec34_mock_realization(
    sample,
    mass_bin_edges,
    random_seed=726,
    include_mass_error_proxy=True,
    mass_error_base=0.01,
):
    x_true = np.asarray(sample["log_sigma_star"], dtype=float)
    y_obs = np.asarray(sample["log_sigma_sfr"], dtype=float)
    z_obs = np.asarray(sample["oh"], dtype=float)
    sigma_y = np.asarray(sample["sigma_log_sfr"], dtype=float)
    sigma_z = np.asarray(sample["sigma_oh"], dtype=float)

    finite = (
        np.isfinite(x_true)
        & np.isfinite(y_obs)
        & np.isfinite(z_obs)
        & np.isfinite(sigma_y)
        & np.isfinite(sigma_z)
    )
    if not np.all(finite):
        sample = mask_sample(sample, finite)
        x_true = np.asarray(sample["log_sigma_star"], dtype=float)
        y_obs = np.asarray(sample["log_sigma_sfr"], dtype=float)
        z_obs = np.asarray(sample["oh"], dtype=float)
        sigma_y = np.asarray(sample["sigma_log_sfr"], dtype=float)
        sigma_z = np.asarray(sample["sigma_oh"], dtype=float)

    model = _build_primary_bin_model(sample, mass_bin_edges)
    bin_index = _digitize_to_bins(x_true, mass_bin_edges)
    rng = np.random.default_rng(random_seed)

    y_true = model["mean_sfr"][bin_index] + rng.normal(
        loc=0.0,
        scale=model["intrinsic_sfr"][bin_index],
    )
    z_true = model["mean_oh"][bin_index] + rng.normal(
        loc=0.0,
        scale=model["intrinsic_oh"][bin_index],
    )

    if include_mass_error_proxy:
        sigma_x = _mass_error_proxy(x_true, base_error=mass_error_base)
        x_mock = x_true + rng.normal(loc=0.0, scale=sigma_x)
    else:
        sigma_x = np.zeros_like(x_true)
        x_mock = x_true.copy()

    y_mock = y_true + rng.normal(loc=0.0, scale=sigma_y)
    z_mock = z_true + rng.normal(loc=0.0, scale=sigma_z)

    return {
        "log_sigma_star": x_mock,
        "log_sigma_sfr": y_mock,
        "oh": z_mock,
        "log_sigma_star_true": x_true,
        "log_sigma_sfr_true": y_true,
        "oh_true": z_true,
        "sigma_log_sigma_star": sigma_x,
        "sigma_log_sfr": sigma_y,
        "sigma_oh": sigma_z,
        "primary_model": model,
        "mass_bin_edges": np.asarray(mass_bin_edges, dtype=float),
    }


def run_sec34_style_null_test_by_galaxy(
    full_sample,
    galaxy_order=None,
    n_bins=6,
    lower_percentile=2.5,
    upper_percentile=97.5,
    n_realizations=120,
    random_seed=726,
    include_mass_error_proxy=True,
    mass_error_base=0.01,
    bin_width=0.2,
    min_unique=20,
    min_trend_points=3,
):
    if galaxy_order is None:
        galaxy_order = sorted(np.unique(full_sample["galaxy_name"]).tolist())

    galaxy_results = {}
    rng = np.random.default_rng(random_seed)

    for galaxy in galaxy_order:
        mask = np.asarray(full_sample["galaxy_name"] == galaxy)
        if not np.any(mask):
            continue

        galaxy_sample = mask_sample(full_sample, mask)
        mass_bin_edges = make_percentile_mass_bin_edges(
            galaxy_sample["log_sigma_star"],
            n_bins=n_bins,
            lower_percentile=lower_percentile,
            upper_percentile=upper_percentile,
        )
        galaxy_seed = int(rng.integers(0, np.iinfo(np.uint32).max))
        galaxy_results[galaxy] = run_sec34_style_null_test(
            galaxy_sample,
            mass_bin_edges=mass_bin_edges,
            n_realizations=n_realizations,
            random_seed=galaxy_seed,
            include_mass_error_proxy=include_mass_error_proxy,
            mass_error_base=mass_error_base,
            bin_width=bin_width,
            min_unique=min_unique,
            min_trend_points=min_trend_points,
        )

    return {
        "galaxy_order": [galaxy for galaxy in galaxy_order if galaxy in galaxy_results],
        "galaxy_results": galaxy_results,
        "config": {
            "n_bins": n_bins,
            "lower_percentile": lower_percentile,
            "upper_percentile": upper_percentile,
            "n_realizations": n_realizations,
            "random_seed": random_seed,
            "include_mass_error_proxy": include_mass_error_proxy,
            "mass_error_base": mass_error_base,
            "bin_width": bin_width,
            "min_unique": min_unique,
            "min_trend_points": min_trend_points,
        },
    }


def plot_sec34_galaxy_grid(
    galaxy_results,
    galaxy_order=None,
    combined_curve=None,
    excluded_galaxies=("NGC4383",),
    max_cols=6,
    mask_null_to_observed=True,
    panel_mode="full",
):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    if galaxy_order is None:
        galaxy_order = list(galaxy_results.keys())

    if panel_mode not in {"full", "rho_only"}:
        raise ValueError(f"Unsupported panel_mode={panel_mode!r}.")

    if panel_mode == "rho_only":
        galaxies_per_row = max(1, max_cols // 3)
        n_cols = galaxies_per_row
        n_rows = int(np.ceil(len(galaxy_order) / galaxies_per_row))
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(n_cols * 4.2, n_rows * 3.35),
            squeeze=False,
        )
    else:
        n_panels_per_galaxy = 3
        galaxies_per_row = max_cols // n_panels_per_galaxy
        if galaxies_per_row < 1:
            raise ValueError("max_cols must be at least 3.")

        n_rows = int(np.ceil(len(galaxy_order) / galaxies_per_row))
        fig, axes = plt.subplots(
            n_rows,
            max_cols,
            figsize=(max_cols * 3.15, n_rows * 3.35),
            squeeze=False,
        )

    for ax in axes.ravel():
        ax.set_visible(False)

    excluded = set(excluded_galaxies)

    for idx, galaxy in enumerate(galaxy_order):
        result = galaxy_results[galaxy]
        row = idx // galaxies_per_row
        block = idx % galaxies_per_row

        if panel_mode == "rho_only":
            ax_rho = axes[row, block]
            ax_rho.set_visible(True)
        else:
            col0 = block * n_panels_per_galaxy
            ax_obs = axes[row, col0]
            ax_null = axes[row, col0 + 1]
            ax_rho = axes[row, col0 + 2]
            ax_obs.set_visible(True)
            ax_null.set_visible(True)
            ax_rho.set_visible(True)

        observed = result["observed"]
        mock = result["null"]["example_cloud"]
        centers = observed["centers"]
        observed_valid_rho = np.isfinite(observed["rho"])

        if mask_null_to_observed:
            null_p16 = np.where(observed_valid_rho, result["null"]["p16_rho"], np.nan)
            null_p84 = np.where(observed_valid_rho, result["null"]["p84_rho"], np.nan)
            null_median = np.where(observed_valid_rho, result["null"]["median_rho"], np.nan)
        else:
            null_p16 = result["null"]["p16_rho"]
            null_p84 = result["null"]["p84_rho"]
            null_median = result["null"]["median_rho"]

        if panel_mode != "rho_only":
            dx_arrays = []
            dz_arrays = []
            if observed["delta_sfr"].size > 0 and observed["delta_oh"].size > 0:
                dx_arrays.append(observed["delta_sfr"])
                dz_arrays.append(observed["delta_oh"])
            if mock is not None and mock["delta_sfr"].size > 0 and mock["delta_oh"].size > 0:
                dx_arrays.append(mock["delta_sfr"])
                dz_arrays.append(mock["delta_oh"])

            if dx_arrays and dz_arrays:
                x_lim = np.nanpercentile(np.abs(np.concatenate(dx_arrays)), 99.0)
                y_lim = np.nanpercentile(np.abs(np.concatenate(dz_arrays)), 99.0)
            else:
                x_lim = 1.0
                y_lim = 0.1
            x_lim = max(float(x_lim), 0.1)
            y_lim = max(float(y_lim), 0.05)

            if observed["delta_sfr"].size > 0 and observed["delta_oh"].size > 0:
                ax_obs.hexbin(
                    observed["delta_sfr"],
                    observed["delta_oh"],
                    gridsize=28,
                    mincnt=1,
                    cmap="Greys",
                    linewidths=0.0,
                )
            ax_obs.axhline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
            ax_obs.axvline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
            ax_obs.set_xlim(-x_lim, x_lim)
            ax_obs.set_ylim(-y_lim, y_lim)
            ax_obs.grid(True, alpha=0.15)
            if galaxy in excluded:
                ax_obs.set_title(f"{galaxy} (excluded)\nObserved", fontsize=9)
            else:
                ax_obs.set_title(f"{galaxy}\nObserved", fontsize=9)
            ax_obs.set_xlabel(r"$\Delta \log\,\Sigma_{\rm SFR}$", fontsize=7)
            ax_obs.set_ylabel(r"$\Delta$[O/H]", fontsize=7)
            ax_obs.tick_params(labelsize=7)

            if mock is not None and mock["delta_sfr"].size > 0 and mock["delta_oh"].size > 0:
                ax_null.hexbin(
                    mock["delta_sfr"],
                    mock["delta_oh"],
                    gridsize=28,
                    mincnt=1,
                    cmap="Oranges",
                    linewidths=0.0,
                )
            ax_null.axhline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
            ax_null.axvline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
            ax_null.set_xlim(-x_lim, x_lim)
            ax_null.set_ylim(-y_lim, y_lim)
            ax_null.grid(True, alpha=0.15)
            ax_null.set_title("Null", fontsize=9)
            ax_null.set_xlabel(r"$\Delta \log\,\Sigma_{\rm SFR}$", fontsize=7)
            ax_null.set_ylabel(r"$\Delta$[O/H]", fontsize=7)
            ax_null.tick_params(labelsize=7)

        ax_rho.fill_between(
            centers,
            null_p16,
            null_p84,
            color="tab:orange",
            alpha=0.25,
        )
        ax_rho.plot(
            centers,
            null_median,
            color="tab:orange",
            linewidth=1.8,
        )
        if combined_curve is not None:
            ref_centers = np.asarray(combined_curve["centers"], dtype=float)
            ref_rho = np.asarray(combined_curve["rho"], dtype=float)
            ref_valid = np.isfinite(ref_centers) & np.isfinite(ref_rho)
            if np.any(ref_valid):
                ax_rho.plot(
                    ref_centers[ref_valid],
                    ref_rho[ref_valid],
                    color="tab:blue",
                    linewidth=1.3,
                    linestyle=":",
                    alpha=0.95,
                )
        obs_linestyle = "--" if galaxy in excluded else "-"
        ax_rho.plot(
            centers,
            observed["rho"],
            color="black",
            linewidth=1.8,
            linestyle=obs_linestyle,
            marker="o",
            markersize=3,
        )
        ax_rho.axhline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
        ax_rho.grid(True, alpha=0.2)
        ax_rho.set_xlabel(r"$\log\,\Sigma_*$ bin center", fontsize=7)
        ax_rho.set_ylabel(r"Spearman $\rho$", fontsize=7)
        ax_rho.tick_params(labelsize=7)
        if galaxy in excluded:
            ax_rho.set_title(f"{galaxy} (excluded)", fontsize=9)
        elif panel_mode == "rho_only":
            ax_rho.set_title(galaxy, fontsize=9)
        ax_rho.text(
            0.03,
            0.03,
            f"N={result['sample_size']}",
            transform=ax_rho.transAxes,
            ha="left",
            va="bottom",
            fontsize=7,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, linewidth=0.0),
        )

    if panel_mode == "rho_only":
        legend_handles = [
            Patch(facecolor="tab:orange", alpha=0.25, label="Null 16-84%"),
            Line2D([0], [0], color="tab:orange", linewidth=1.8, label="Null median"),
            Line2D([0], [0], color="black", linewidth=1.8, label="Galaxy observed"),
            Line2D([0], [0], color="black", linewidth=1.8, linestyle="--", label="Excluded galaxy observed"),
        ]
    else:
        legend_handles = [
            Patch(facecolor="0.5", alpha=0.6, label="Observed cloud"),
            Patch(facecolor="tab:orange", alpha=0.25, label="Null 16-84%"),
            Line2D([0], [0], color="tab:orange", linewidth=1.8, label="Null median"),
            Line2D([0], [0], color="black", linewidth=1.8, label="Galaxy observed"),
            Line2D([0], [0], color="black", linewidth=1.8, linestyle="--", label="Excluded galaxy observed"),
        ]
    if combined_curve is not None:
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color="tab:blue",
                linewidth=1.3,
                linestyle=":",
                label=combined_curve.get("label", "Combined observed"),
            )
        )

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=min(len(legend_handles), 6),
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, 0.995),
    )
    if panel_mode == "rho_only":
        fig.suptitle(
            "Per-galaxy null tests\nGalaxy-specific Spearman $\\rho$ trends",
            fontsize=14,
            y=1.015,
        )
        plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.955])
    else:
        fig.suptitle(
            "Per-galaxy null tests\nObserved residual cloud, one null realization, and galaxy-specific rho trend",
            fontsize=14,
            y=1.015,
        )
        plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.965])
    plt.show()
    return fig


def compute_offset_maps_from_mass_bins(
    log_sigma_star_map,
    log_sigma_sfr_map,
    oh_map,
    valid_mask,
    n_bins=12,
):
    x = np.asarray(log_sigma_star_map, dtype=float)
    y = np.asarray(log_sigma_sfr_map, dtype=float)
    z = np.asarray(oh_map, dtype=float)
    valid_mask = np.asarray(valid_mask, dtype=bool)

    if np.sum(valid_mask) == 0:
        raise RuntimeError("Cannot compute offset maps from an empty valid mask.")

    mass_min = np.nanmin(x[valid_mask])
    mass_max = np.nanmax(x[valid_mask])
    sigmaM_bin_edges = np.linspace(mass_min, mass_max, n_bins + 1)

    sfr_offset = np.full_like(y, np.nan, dtype=float)
    oh_offset = np.full_like(z, np.nan, dtype=float)
    sigM_binned_map = np.full_like(x, np.nan, dtype=float)

    for i in range(n_bins):
        bin_mask = _in_bin_mask(x, sigmaM_bin_edges[i], sigmaM_bin_edges[i + 1], i == n_bins - 1) & valid_mask
        if np.sum(bin_mask) == 0:
            continue

        bin_center = 0.5 * (sigmaM_bin_edges[i] + sigmaM_bin_edges[i + 1])
        sfr_bin_mean = np.nanmean(y[bin_mask])
        oh_bin_mean = np.nanmean(z[bin_mask])

        sfr_offset[bin_mask] = y[bin_mask] - sfr_bin_mean
        oh_offset[bin_mask] = z[bin_mask] - oh_bin_mean
        sigM_binned_map[bin_mask] = bin_center

    sfr_offset_max = np.nanmax(np.abs(sfr_offset[valid_mask]))
    oh_offset_max = np.nanmax(np.abs(oh_offset[valid_mask]))

    return {
        "sigmaM_bin_edges": sigmaM_bin_edges,
        "sfr_offset": sfr_offset,
        "oh_offset": oh_offset,
        "sigM_binned_map": sigM_binned_map,
        "sfr_offset_max": max(float(sfr_offset_max), 0.1),
        "oh_offset_max": max(float(oh_offset_max), 0.05),
    }


def local_bivariate_moran(A, B, valid_mask, neighbourhood="queen"):
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    valid_mask = np.asarray(valid_mask, dtype=bool)

    ny, nx = A.shape
    I_map = np.full_like(A, np.nan, dtype=float)
    quad_map = np.zeros_like(A, dtype=int)

    Aval = A[valid_mask]
    Bval = B[valid_mask]
    Astd = np.nanstd(Aval)
    Bstd = np.nanstd(Bval)
    if not np.isfinite(Astd) or Astd == 0:
        Astd = 1.0
    if not np.isfinite(Bstd) or Bstd == 0:
        Bstd = 1.0

    Az = (A - np.nanmean(Aval)) / Astd
    Bz = (B - np.nanmean(Bval)) / Bstd

    if neighbourhood == "queen":
        offsets = [
            (-1, -1), (-1, 0), (-1, 1),
            (0, -1),            (0, 1),
            (1, -1),  (1, 0),   (1, 1),
        ]
    elif neighbourhood == "rook":
        offsets = [(-1, 0), (0, -1), (0, 1), (1, 0)]
    else:
        raise ValueError("neighbourhood must be 'queen' or 'rook'")

    for i in range(ny):
        for j in range(nx):
            if not valid_mask[i, j]:
                continue

            neigh_vals = []
            for dy, dx in offsets:
                ii = i + dy
                jj = j + dx
                if 0 <= ii < ny and 0 <= jj < nx and valid_mask[ii, jj]:
                    neigh_vals.append(Bz[ii, jj])

            if not neigh_vals:
                continue

            lag_B = float(np.mean(neigh_vals))
            I_val = Az[i, j] * lag_B
            I_map[i, j] = I_val

            if Az[i, j] >= 0 and lag_B >= 0:
                quad_map[i, j] = 1
            elif Az[i, j] < 0 and lag_B >= 0:
                quad_map[i, j] = 2
            elif Az[i, j] < 0 and lag_B < 0:
                quad_map[i, j] = 3
            else:
                quad_map[i, j] = 4

    return I_map, quad_map


def _world_extent_from_header(gas_header, shape):
    from astropy.wcs import WCS

    wcs_celestial = WCS(gas_header).celestial
    y_size, x_size = shape
    x_coords = np.arange(x_size)
    y_coords = np.arange(y_size)
    xx, yy = np.meshgrid(x_coords, y_coords)
    ra, dec = wcs_celestial.pixel_to_world_values(xx, yy)
    return [float(ra.max()), float(ra.min()), float(dec.min()), float(dec.max())]


def _scalebar_length_deg(pix_scale_arcsec=0.2, distance_mpc=16.5):
    kpc_arcsec = 206265.0 / (distance_mpc * 1000.0)
    return kpc_arcsec / 3600.0


def build_simulated_moran_maps_for_galaxy(
    galaxy,
    random_seed=726,
    null_n_bins=6,
    offset_n_bins=12,
    include_mass_error_proxy=True,
    mass_error_base=0.01,
    neighbourhood="queen",
    threshold=1.0,
):
    galaxy_maps = load_o3n2_m13_hii_galaxy_maps_with_errors(galaxy)
    sample = galaxy_maps_to_sample(galaxy_maps)

    null_mass_bin_edges = make_percentile_mass_bin_edges(
        sample["log_sigma_star"],
        n_bins=null_n_bins,
        lower_percentile=2.5,
        upper_percentile=97.5,
    )
    mock = generate_sec34_mock_realization(
        sample,
        null_mass_bin_edges,
        random_seed=random_seed,
        include_mass_error_proxy=include_mass_error_proxy,
        mass_error_base=mass_error_base,
    )

    valid_mask = np.asarray(galaxy_maps["valid_mask"], dtype=bool)
    shape = galaxy_maps["log_sigma_star_map"].shape

    x_mock_map = np.full(shape, np.nan, dtype=float)
    y_mock_map = np.full(shape, np.nan, dtype=float)
    z_mock_map = np.full(shape, np.nan, dtype=float)
    x_mock_map[valid_mask] = mock["log_sigma_star"]
    y_mock_map[valid_mask] = mock["log_sigma_sfr"]
    z_mock_map[valid_mask] = mock["oh"]

    offset_data = compute_offset_maps_from_mass_bins(
        x_mock_map,
        y_mock_map,
        z_mock_map,
        valid_mask,
        n_bins=offset_n_bins,
    )
    I_map, quad_map = local_bivariate_moran(
        offset_data["sfr_offset"],
        offset_data["oh_offset"],
        valid_mask,
        neighbourhood=neighbourhood,
    )

    valid_I = np.isfinite(I_map)
    n_total = int(np.sum(valid_I))
    n_high = int(np.sum(I_map[valid_I] > threshold))
    n_mid = int(np.sum((I_map[valid_I] >= -threshold) & (I_map[valid_I] <= threshold)))
    n_low = int(np.sum(I_map[valid_I] < -threshold))

    return {
        "galaxy": galaxy,
        "threshold": float(threshold),
        "neighbourhood": neighbourhood,
        "valid_mask": valid_mask,
        "null_mass_bin_edges": null_mass_bin_edges,
        "mock_realization": mock,
        "x_mock_map": x_mock_map,
        "y_mock_map": y_mock_map,
        "z_mock_map": z_mock_map,
        "I_map": I_map,
        "quad_map": quad_map,
        "extent": _world_extent_from_header(galaxy_maps["gas_header"], shape),
        "scalebar_length_deg": _scalebar_length_deg(),
        "stats": {
            "n_total": n_total,
            "n_high": n_high,
            "n_mid": n_mid,
            "n_low": n_low,
            "pct_high": (100.0 * n_high / n_total) if n_total > 0 else 0.0,
            "pct_mid": (100.0 * n_mid / n_total) if n_total > 0 else 0.0,
            "pct_low": (100.0 * n_low / n_total) if n_total > 0 else 0.0,
        },
        **offset_data,
    }


def plot_simulated_moran_2x2(sim_result):
    import matplotlib.gridspec as gs_module  # noqa: F401
    import matplotlib.ticker as mticker
    from matplotlib.colors import BoundaryNorm, SymLogNorm
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    galaxy = sim_result["galaxy"]
    extent = sim_result["extent"]
    threshold = sim_result["threshold"]
    sfr_offset = sim_result["sfr_offset"]
    oh_offset = sim_result["oh_offset"]
    I_map = sim_result["I_map"]
    sigM_binned_map = sim_result["sigM_binned_map"]
    sigmaM_bin_edges = sim_result["sigmaM_bin_edges"]

    dx = abs(extent[1] - extent[0])
    dy = abs(extent[3] - extent[2])
    data_ratio = (dx / dy) if dy > 0 else 1.0
    panel_h = 4.2
    panel_w = panel_h * data_ratio
    cbar_w = 0.35
    W = 2 * panel_w + 2 * cbar_w + 1.2
    H = 2 * panel_h + 0.8

    fig, axs = plt.subplots(2, 2, figsize=(W, H), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.96, bottom=0.07, top=0.95, wspace=0.20, hspace=0.06)
    ax1, ax2 = axs[0, 0], axs[0, 1]
    ax3, ax4 = axs[1, 0], axs[1, 1]

    def add_cbar(ax, im, label, size="3.5%", pad=0.02, labelpad=6, labelsize=11):
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=size, pad=pad)
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label(label, fontsize=labelsize, labelpad=labelpad)
        return cbar

    def add_scalebar(ax, length_deg, label="1 kpc"):
        x0, x1 = ax.get_xlim()
        x_range = abs(x1 - x0)
        if x_range == 0:
            return
        frac_len = length_deg / x_range
        x_start = 0.05
        y_start = 0.95
        x_end = x_start + frac_len
        line, = ax.plot(
            [x_start, x_end],
            [y_start, y_start],
            transform=ax.transAxes,
            color="k",
            lw=3,
            clip_on=False,
        )
        line.set_in_layout(False)
        txt = ax.text(
            (x_start + x_end) / 2,
            y_start - 0.03,
            label,
            transform=ax.transAxes,
            ha="center",
            va="top",
            color="k",
        )
        txt.set_in_layout(False)

    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_major_locator(mticker.MaxNLocator(3))
        ax.yaxis.set_major_locator(mticker.MaxNLocator(5))
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)

    ax3.set_xlabel("RA (deg)", fontsize=16)
    ax4.set_xlabel("RA (deg)", fontsize=16)
    ax1.set_ylabel("Dec (deg)", fontsize=16)
    ax3.set_ylabel("Dec (deg)", fontsize=16)

    im_sfr = ax1.imshow(
        np.ma.masked_invalid(sfr_offset),
        origin="lower",
        cmap="coolwarm",
        extent=extent,
        vmin=-sim_result["sfr_offset_max"],
        vmax=sim_result["sfr_offset_max"],
        aspect="equal",
    )
    ax1.set_title(f"{galaxy}: simulated " + r"$\Delta\log\,\Sigma_{\mathrm{SFR}}$", fontsize=13)
    add_scalebar(ax1, sim_result["scalebar_length_deg"], label="1 kpc")
    add_cbar(ax1, im_sfr, r"$\Delta\log\,\Sigma_{\mathrm{SFR}}$ (dex)", labelpad=4)

    im_oh = ax2.imshow(
        np.ma.masked_invalid(oh_offset),
        origin="lower",
        cmap="coolwarm",
        extent=extent,
        vmin=-sim_result["oh_offset_max"],
        vmax=sim_result["oh_offset_max"],
        aspect="equal",
    )
    ax2.set_title(f"{galaxy}: simulated " + r"$\Delta$[O/H]", fontsize=13)
    add_cbar(ax2, im_oh, r"$\Delta$[O/H] (dex)")

    vI = np.nanmax(np.abs(I_map))
    if not np.isfinite(vI) or vI <= 0:
        vI = threshold
    norm = SymLogNorm(linthresh=threshold, linscale=1.0, vmin=-vI, vmax=vI, base=10)
    imI = ax3.imshow(
        np.ma.masked_invalid(I_map),
        origin="lower",
        cmap="coolwarm",
        extent=extent,
        norm=norm,
        aspect="equal",
    )
    ax3.set_title(f"{galaxy}: simulated Moran-like map", fontsize=13)
    cbarI = add_cbar(
        ax3,
        imI,
        r"$z(\Delta\log\Sigma_{\rm SFR}) \times z_{\rm lag}(\Delta{\rm [O/H]})$",
    )
    max_pow = int(np.floor(np.log10(vI))) if vI > 0 else 0
    pos_ticks = [10**k for k in range(0, max_pow + 1) if 10**k <= vI]
    pos_ticks_filtered = [p for p in pos_ticks if p != threshold]
    neg_ticks_filtered = [-p for p in reversed(pos_ticks_filtered)]
    ticks = sorted(set(neg_ticks_filtered + [-threshold, 0.0, threshold] + pos_ticks_filtered))
    cbarI.set_ticks(ticks)
    cbarI.set_ticklabels([f"{tick:g}" for tick in ticks])

    legend_text = (
        rf"$I_i > {threshold:g}$: {sim_result['stats']['pct_high']:.1f}% (N={sim_result['stats']['n_high']})"
        + "\n"
        + rf"$|I_i| \leq {threshold:g}$: {sim_result['stats']['pct_mid']:.1f}% (N={sim_result['stats']['n_mid']})"
        + "\n"
        + rf"$I_i < -{threshold:g}$: {sim_result['stats']['pct_low']:.1f}% (N={sim_result['stats']['n_low']})"
    )
    ax3.text(
        0.98,
        0.02,
        legend_text,
        transform=ax3.transAxes,
        fontsize=10,
        va="bottom",
        ha="right",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

    cmap_discrete = plt.cm.inferno
    norm_discrete = BoundaryNorm(sigmaM_bin_edges, cmap_discrete.N)
    im_sigM_binned = ax4.imshow(
        np.ma.masked_invalid(sigM_binned_map),
        origin="lower",
        cmap=cmap_discrete,
        norm=norm_discrete,
        extent=extent,
        aspect="equal",
    )
    ax4.set_title(f"{galaxy}: simulated binned " + r"$\log\,\Sigma_*$", fontsize=13)
    cbar_binned = add_cbar(
        ax4,
        im_sigM_binned,
        r"$\log\,\Sigma_*\;(M_\odot\,\mathrm{kpc}^{-2})$",
    )
    bin_centers = 0.5 * (sigmaM_bin_edges[:-1] + sigmaM_bin_edges[1:])
    cbar_binned.set_ticks(bin_centers)
    cbar_binned.set_ticklabels([f"{center:.2f}" for center in bin_centers])

    plt.show()
    return fig
