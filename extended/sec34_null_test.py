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
        trend_centers, _, _, _ = calculate_binned_stats(
            trend_x,
            trend_y,
            bin_width=bin_width,
            min_unique=min_unique,
        )
        if trend_centers.size < min_trend_points:
            continue

        rho[i], p_value[i] = spearmanr(trend_x, trend_y)
        delta_sfr_cloud.append(trend_x)
        delta_oh_cloud.append(trend_y)
        mass_bin_index_cloud.append(np.full(trend_x.size, i, dtype=int))

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


def _mass_error_proxy(log_sigma_star, base_error=0.05):
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
    mass_error_base=0.05,
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


def plot_sec34_style_null_test(result, reference_curve=None):
    observed = result["observed"]
    mock = result["null"]["example_cloud"]
    centers = observed["centers"]

    all_dx = [observed["delta_sfr"]]
    all_dz = [observed["delta_oh"]]
    if mock is not None:
        all_dx.append(mock["delta_sfr"])
        all_dz.append(mock["delta_oh"])

    x_lim = np.nanpercentile(np.abs(np.concatenate(all_dx)), 99.0)
    y_lim = np.nanpercentile(np.abs(np.concatenate(all_dz)), 99.0)
    x_lim = max(x_lim, 0.1)
    y_lim = max(y_lim, 0.05)

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

    axes[2].fill_between(
        centers,
        result["null"]["p16_rho"],
        result["null"]["p84_rho"],
        color="tab:orange",
        alpha=0.25,
        label="Null 16-84%",
    )
    axes[2].plot(
        centers,
        result["null"]["median_rho"],
        color="tab:orange",
        linewidth=2.0,
        label="Null median",
    )
    axes[2].plot(
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
            axes[2].plot(
                ref_centers[ref_valid],
                ref_rho[ref_valid],
                color="tab:blue",
                linewidth=1.8,
                linestyle="--",
                marker="s",
                markersize=4,
                label=reference_curve.get("label", "Reference observed curve"),
            )
    axes[2].axhline(0.0, color="gray", linestyle="--", linewidth=1.0, alpha=0.6)
    axes[2].set_xlabel(r"$\log\,\Sigma_*$ bin center", fontsize=11)
    axes[2].set_ylabel(r"Spearman $\rho$", fontsize=11)
    axes[2].set_title(r"Observed vs null $\rho(\log\,\Sigma_*)$", fontsize=12)
    axes[2].grid(True, alpha=0.25)
    axes[2].legend(loc="best", fontsize=10)

    if result["config"]["include_mass_error_proxy"]:
        mass_label = (
            "Includes Sanchez-style mass-error proxy "
            f"(base={result['config']['mass_error_base']:.3f} dex)"
        )
    else:
        mass_label = "Mass distribution fixed; no extra mass-error proxy"

    fig.suptitle(
        "Sec. 3.4-style null test: independent primary relations + propagated errors\n"
        + mass_label,
        fontsize=13,
    )
    plt.tight_layout()
    plt.show()
