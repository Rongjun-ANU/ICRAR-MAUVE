if old_paper_cells_ready(require_ew=True):
    # ------------------------------------------------------------------
    # Source cell 17 recreated with fixed galaxy colors: Δ[O/H] vs Δlog Σ_SFR in Σ_* bins
    # For all galaxies (excluding NGC4383) and Spearman rho vs Σ_* summary
    # Spearman rho uses only points inside 95% contour and only if trend exists
    # ------------------------------------------------------------------
    from pathlib import Path
    from astropy.io import fits

    EW_FILE_SUFFIX = '_proxy_EW_maps.fits'
    EW_HDU = 'PROXY_EWHA'
    EW_MIN = 6.0

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    import matplotlib.gridspec as gridspec
    from scipy.stats import spearmanr
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    # ------------------------------------------------------------------
    # Helper: read the required maps
    # ------------------------------------------------------------------
    def load_maps(gal):
        spatial_path = further_spatial_path(gal)
        gas_path = further_gas_path(gal)
        ew_path = further_ew_path(gal, EW_FILE_SUFFIX)

        if not (Path(spatial_path).exists() and Path(gas_path).exists() and Path(ew_path).exists()):
            raise FileNotFoundError(f'Required FITS files for {gal} not found.')

        with fits.open(spatial_path) as h_spatial:
            sigM = h_spatial['LOGMASS_SURFACE_DENSITY'].data
        with fits.open(gas_path) as h_gas:
            sigSFR = h_gas['LOGSFR_SURFACE_DENSITY_HII'].data
            ha_sigma_obs = h_gas['HA6562_SIGMA'].data
            ha_sigma_corr = h_gas['HA6562_SIGMA_CORR'].data if 'HA6562_SIGMA_CORR' in h_gas else None
            if ha_sigma_corr is None:
                raise KeyError(f'{gal}: missing HA6562_SIGMA_CORR in {gas_path}')
            if ha_sigma_corr.shape != ha_sigma_obs.shape:
                raise ValueError(f'{gal}: shape mismatch SIGMA_CORR {ha_sigma_corr.shape} vs SIGMA {ha_sigma_obs.shape}')
            sigma2_intrinsic = ha_sigma_obs**2 - ha_sigma_corr**2
            ha_sigma = np.full_like(ha_sigma_obs, np.nan, dtype=float)
            sigma_ok = np.isfinite(sigma2_intrinsic) & (sigma2_intrinsic > 0)
            ha_sigma[sigma_ok] = np.sqrt(sigma2_intrinsic[sigma_ok])
            if 'O_H_O3N2_M14_HII' in h_gas:
                oh_map = h_gas['O_H_O3N2_M14_HII'].data
                indicator_name = 'O3N2-M14'
            elif 'O_H_O3N2_M13_HII' in h_gas:
                print(f"Warning: O3N2-M14 not available for {gal}; using O3N2-M13 instead.")
                oh_map = h_gas['O_H_O3N2_M13_HII'].data
                indicator_name = 'O3N2-M13'
            else:
                available = [hdu.name for hdu in h_gas]
                raise KeyError(f"No O3N2 metallicity extension (M14/M13) found for {gal}. Available: {available}")
        with fits.open(ew_path) as h_ew:
            if EW_HDU not in h_ew:
                raise KeyError(f"Missing {EW_HDU} in {ew_path}")
            ew_ha = h_ew[EW_HDU].data
        ew_ha = align_map_to_shape(ew_ha, sigSFR.shape, label=f"{gal} {EW_HDU}")
        ew_mask = np.isfinite(ew_ha) & (ew_ha > EW_MIN)
        ew_mask = ew_mask & np.isfinite(ha_sigma) & (ha_sigma < 45.0)
        sigSFR = np.where(ew_mask, sigSFR, np.nan)
        oh_map = np.where(ew_mask, oh_map, np.nan)
        return sigM, sigSFR, oh_map, indicator_name

    # ------------------------------------------------------------------
    # Helper: calculate binned statistics with unique-count requirement
    # ------------------------------------------------------------------
    def calculate_binned_stats(x_data, y_data, bin_width=0.2, min_unique=20):
        """Calculate median and std in x-bins, requiring >= min_unique unique y values."""
        x_data = np.asarray(x_data)
        y_data = np.asarray(y_data)

        finite = np.isfinite(x_data) & np.isfinite(y_data)
        if np.sum(finite) == 0:
            return np.array([]), np.array([]), np.array([]), np.array([])

        x_min, x_max = np.nanmin(x_data[finite]), np.nanmax(x_data[finite])
        if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min == x_max:
            return np.array([]), np.array([]), np.array([]), np.array([])

        bin_edges = np.arange(x_min, x_max + bin_width, bin_width)
        if bin_edges.size < 2:
            return np.array([]), np.array([]), np.array([]), np.array([])

        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        medians = []
        stds = []
        counts = []
        valid_centers = []

        for i, (bin_start, bin_end) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
            in_bin = finite & (x_data >= bin_start) & (x_data < bin_end)
            if np.sum(in_bin) == 0:
                continue

            y_in_bin = y_data[in_bin]
            if len(np.unique(y_in_bin)) >= min_unique:
                medians.append(np.nanmedian(y_in_bin))
                stds.append(np.nanstd(y_in_bin))
                counts.append(len(y_in_bin))
                valid_centers.append(bin_centers[i])

        return (np.array(valid_centers), np.array(medians),
                np.array(stds), np.array(counts))

    # ------------------------------------------------------------------
    # Galaxy list and color map (reserve color for excluded galaxy)
    # ------------------------------------------------------------------
    all_bin_centers = []  # Σ_* bin centers for which rho is valid (per-galaxy)
    all_rho_values = []     # Spearman rho (per-galaxy)
    all_galaxy_indices = []

    # Also store all spaxels across galaxies to compute a global trend
    global_sigM_list = []
    global_sfr_off_list = []
    global_oh_off_list = []
    global_gal_index_list = []

    # Build full galaxy list from available FITS files (including excluded)
    galaxies_all = physics_ready_galaxies(require_ew=True)
    excluded_galaxies = {'NGC4383'}

    # Build galaxy -> fixed MAUVE color mapping over ALL galaxies (so excluded keeps its color)
    galaxy_color_dict = {gal: galaxy_color(gal) for gal in galaxies_all}
    print(f"Using fixed {FIXED_GALAXY_PALETTE_NAME} colors for source cell 17 recreation.")

    # Now define the working sample that EXCLUDES NGC4383 for analysis/plotting
    galaxies = [g for g in galaxies_all if g not in excluded_galaxies]

    print("=" * 80)
    print(f"Offset relations and Spearman rho vs Σ_* (excluding: {excluded_galaxies})")
    print("All galaxies (for color assignment):", ", ".join(galaxies_all))
    print("Galaxies used in this analysis:", ", ".join(galaxies))
    print("=" * 80)

    # ------------------------------------------------------------------
    # Loop over galaxies and compute Spearman rho per Σ_* bin (per-galaxy)
    # and collect global spaxels for combined trend
    # ------------------------------------------------------------------
    for gal_index, gal in enumerate(galaxies):
        print(f"\nProcessing {gal} ...")
        try:
            sigM, sigSFR, oh_o3n2_map, indicator_name = load_maps(gal)
        except Exception as e:
            print(f"  Skipping {gal}: {e}")
            continue

        valid_sfr = np.isfinite(sigSFR)
        valid_oh = np.isfinite(oh_o3n2_map)
        valid_sigM = np.isfinite(sigM)
        valid_all = valid_sfr & valid_oh & valid_sigM

        if np.sum(valid_all) == 0:
            print("  No valid pixels.")
            continue

        global_gal_index_list.append(np.full(np.sum(valid_all), gal_index, dtype=int))

        # 1. Define 12 Σ_* bins for this galaxy (spanning 95% population)
        n_bins = 6
        sigmaM_min = np.nanpercentile(sigM[valid_all], 2.5)
        sigmaM_max = np.nanpercentile(sigM[valid_all], 97.5)
        sigmaM_bin_edges = np.linspace(sigmaM_min, sigmaM_max, n_bins + 1)
        sigmaM_bin_centers = 0.5 * (sigmaM_bin_edges[:-1] + sigmaM_bin_edges[1:])

        # 2. Global symmetric limits for this galaxy (for reference / consistency)
        global_max_abs_sfr = 0
        global_max_abs_oh = 0
        for i in range(n_bins):
            if i == 0:
                bin_mask = (sigM >= sigmaM_bin_edges[i]) & (sigM < sigmaM_bin_edges[i+1]) & valid_all
            elif i == n_bins - 1:
                bin_mask = (sigM >= sigmaM_bin_edges[i]) & (sigM <= sigmaM_bin_edges[i+1]) & valid_all
            else:
                bin_mask = (sigM >= sigmaM_bin_edges[i]) & (sigM < sigmaM_bin_edges[i+1]) & valid_all

            if np.sum(bin_mask) > 0:
                sfr_bin = sigSFR[bin_mask]
                oh_bin = oh_o3n2_map[bin_mask]

                sfr_offset = sfr_bin - np.mean(sfr_bin)
                oh_offset = oh_bin - np.mean(oh_bin)

                if np.any(np.isfinite(sfr_offset)):
                    global_max_abs_sfr = max(global_max_abs_sfr, np.nanmax(np.abs(sfr_offset)))
                if np.any(np.isfinite(oh_offset)):
                    global_max_abs_oh = max(global_max_abs_oh, np.nanmax(np.abs(oh_offset)))

        print(f"  Σ* range (95%): {sigmaM_min:.3f} to {sigmaM_max:.3f}")
        print(f"  Global Δlog(Σ_SFR) axis limits: ±{global_max_abs_sfr:.3f}")
        print(f"  Global Δ[O/H] axis limits: ±{global_max_abs_oh:.3f}")

        # 3. Compute Spearman rho in each Σ_* bin using 3-sigma clipping on unique data (per galaxy)
        print("  Spearman rho in each Σ_* bin (3σ clip + trend, per galaxy):")
        for i in range(n_bins):
            if i == 0:
                bin_mask = (sigM >= sigmaM_bin_edges[i]) & (sigM < sigmaM_bin_edges[i+1]) & valid_all
            elif i == n_bins - 1:
                bin_mask = (sigM >= sigmaM_bin_edges[i]) & (sigM <= sigmaM_bin_edges[i+1]) & valid_all
            else:
                bin_mask = (sigM >= sigmaM_bin_edges[i]) & (sigM < sigmaM_bin_edges[i+1]) & valid_all

            n_pixels = np.sum(bin_mask)

            if n_pixels <= 0:
                print(f"    Bin {i+1:02d} ({sigmaM_bin_centers[i]:.2f}): n=0, r=N/A")
                continue

            # Values in this Σ_* bin
            sfr_bin = sigSFR[bin_mask]
            oh_bin = oh_o3n2_map[bin_mask]

            # Offsets within this bin
            sfr_bin_mean = np.mean(sfr_bin)
            oh_bin_mean = np.mean(oh_bin)

            sfr_offset = sfr_bin - sfr_bin_mean
            oh_offset = oh_bin - oh_bin_mean

            valid_corr = np.isfinite(sfr_offset) & np.isfinite(oh_offset)
            if np.sum(valid_corr) > 0:
                global_sigM_list.append(sigM[bin_mask][valid_corr])
                global_sfr_off_list.append(sfr_offset[valid_corr])
                global_oh_off_list.append(oh_offset[valid_corr])

            if np.sum(valid_corr) <= 2:
                spearman_rho = np.nan
                print(f"    Bin {i+1:02d} ({sigmaM_bin_centers[i]:.2f}): n={n_pixels}, r=N/A (insufficient valid points)")
                continue

            # ------------------------------------------------------------------
            # Calculate 3-sigma clipping mask on unique valid data
            # ------------------------------------------------------------------
            sigma_clip_mask = np.zeros_like(valid_corr, dtype=bool)
            if np.sum(valid_corr) > 0:
                # Get unique valid data using np.unique
                sfr_data = sfr_offset[valid_corr]
                oh_data = oh_offset[valid_corr]
            
                sfr_unique = np.unique(sfr_data)
                oh_unique = np.unique(oh_data)
            
                # Calculate mean and std on unique values
                x_mean = np.mean(sfr_unique)
                x_std = np.std(sfr_unique)
                y_mean = np.mean(oh_unique)
                y_std = np.std(oh_unique)
            
                # 3-sigma bounds
                x_lower = x_mean - 3 * x_std
                x_upper = x_mean + 3 * x_std
                y_lower = y_mean - 3 * y_std
                y_upper = y_mean + 3 * y_std
            
                # Apply 3-sigma clipping to all data points
                within_3sigma = (sfr_data >= x_lower) & (sfr_data <= x_upper) & \
                                (oh_data >= y_lower) & (oh_data <= y_upper)
                sigma_clip_mask[valid_corr] = within_3sigma

            # ------------------------------------------------------------------
            # Median trend using only points within 3-sigma
            # ------------------------------------------------------------------
            x_trend = sfr_offset[sigma_clip_mask]
            y_trend = oh_offset[sigma_clip_mask]

            centers, medians, stds, counts = calculate_binned_stats(x_trend, y_trend,
                                                                    bin_width=0.2,
                                                                    min_unique=20)

            # Require at least 2 median-trend points for a valid trend
            has_trend = centers.size >= 3

            # ------------------------------------------------------------------
            # Spearman rho: only for points within 3-sigma and only if trend exists
            # ------------------------------------------------------------------
            if has_trend and np.sum(sigma_clip_mask & valid_corr) > 2:
                try:
                    spearman_rho, spearman_p = spearmanr(sfr_offset[sigma_clip_mask & valid_corr],
                                                    oh_offset[sigma_clip_mask & valid_corr])
                    print(f"    Bin {i+1:02d} ({sigmaM_bin_centers[i]:.2f}): n={n_pixels}, r={spearman_rho:.3f} (p={spearman_p:.3g})")
                except Exception:
                    spearman_rho = np.nan
                    print(f"    Bin {i+1:02d} ({sigmaM_bin_centers[i]:.2f}): n={n_pixels}, r=N/A (error)")
            else:
                spearman_rho = np.nan
                print(f"    Bin {i+1:02d} ({sigmaM_bin_centers[i]:.2f}): n={n_pixels}, r=N/A (no trend or too few in 3σ)")

            # Store only valid rho values
            if np.isfinite(spearman_rho):
                all_bin_centers.append(sigmaM_bin_centers[i])
                all_rho_values.append(spearman_rho)
                all_galaxy_indices.append(gal_index)

    # ------------------------------------------------------------------
    # Prepare global arrays for combined Spearman vs Σ_* (all galaxies)
    # ------------------------------------------------------------------
    if global_sigM_list:
        global_sigM = np.concatenate(global_sigM_list)
        global_sfr_off = np.concatenate(global_sfr_off_list)
        global_oh_off = np.concatenate(global_oh_off_list)
    else:
        global_sigM = np.array([])
        global_sfr_off = np.array([])
        global_oh_off = np.array([])

    # ------------------------------------------------------------------
    # Convert collected per-galaxy Spearman rho data to arrays
    # ------------------------------------------------------------------
    all_bin_centers = np.array(all_bin_centers)
    all_rho_values = np.array(all_rho_values)
    all_galaxy_indices = np.array(all_galaxy_indices)

    print("\n" + "=" * 80)
    print("Summary: number of Σ_* bins with valid Spearman rho (per galaxy):", len(all_rho_values))
    print("=" * 80)

    # ------------------------------------------------------------------
    # Compute combined Spearman rho vs Σ_* over all galaxies (same method)
    # ------------------------------------------------------------------
    global_bin_centers = []
    global_rho_values = []

    if global_sigM.size > 0:
        # Define 12 Σ_* bins for the combined sample (spanning 95% population)
        n_bins_global = 12
        sigmaM_min_g = np.nanpercentile(global_sigM, 2.5)
        sigmaM_max_g = np.nanpercentile(global_sigM, 97.5)
        sigmaM_edges_g = np.linspace(sigmaM_min_g, sigmaM_max_g, n_bins_global + 1)
        sigmaM_centers_g = 0.5 * (sigmaM_edges_g[:-1] + sigmaM_edges_g[1:])

        print("\nComputing combined Spearman rho in Σ_* bins (all galaxies together)...")
        for i in range(n_bins_global):
            if i == 0:
                bin_mask_g = (global_sigM >= sigmaM_edges_g[i]) & (global_sigM < sigmaM_edges_g[i+1])
            elif i == n_bins_global - 1:
                bin_mask_g = (global_sigM >= sigmaM_edges_g[i]) & (global_sigM <= sigmaM_edges_g[i+1])
            else:
                bin_mask_g = (global_sigM >= sigmaM_edges_g[i]) & (global_sigM < sigmaM_edges_g[i+1])

            n_pix_g = np.sum(bin_mask_g)
            if n_pix_g <= 0:
                print(f"  Global bin {i+1:02d} ({sigmaM_centers_g[i]:.2f}): n=0, r=N/A")
                continue

            sfr_off_g = global_sfr_off[bin_mask_g]
            oh_off_g = global_oh_off[bin_mask_g]



            valid_corr_g = np.isfinite(sfr_off_g) & np.isfinite(oh_off_g)
            if np.sum(valid_corr_g) <= 2:
                print(f"  Global bin {i+1:02d} ({sigmaM_centers_g[i]:.2f}): n={n_pix_g}, r=N/A (insufficient valid points)")
                continue

            # ------------------------------------------------------------------
            # Calculate 3-sigma clipping mask on unique valid data (global)
            # ------------------------------------------------------------------
            sigma_clip_mask_global = np.zeros_like(valid_corr_g, dtype=bool)
            if np.sum(valid_corr_g) > 0:
                sfr_data_g = sfr_off_g[valid_corr_g]
                oh_data_g = oh_off_g[valid_corr_g]
            
                sfr_unique_g = np.unique(sfr_data_g)
                oh_unique_g = np.unique(oh_data_g)
            
                x_mean_g = np.mean(sfr_unique_g)
                x_std_g = np.std(sfr_unique_g)
                y_mean_g = np.mean(oh_unique_g)
                y_std_g = np.std(oh_unique_g)
            
                x_lower_g = x_mean_g - 3 * x_std_g
                x_upper_g = x_mean_g + 3 * x_std_g
                y_lower_g = y_mean_g - 3 * y_std_g
                y_upper_g = y_mean_g + 3 * y_std_g
            
                within_3sigma_g = (sfr_data_g >= x_lower_g) & (sfr_data_g <= x_upper_g) & \
                                  (oh_data_g >= y_lower_g) & (oh_data_g <= y_upper_g)
                sigma_clip_mask_global[valid_corr_g] = within_3sigma_g

            x_trend_g = sfr_off_g[sigma_clip_mask_global]
            y_trend_g = oh_off_g[sigma_clip_mask_global]

            centers_g, medians_g, stds_g, counts_g = calculate_binned_stats(x_trend_g, y_trend_g,
                                                                            bin_width=0.2,
                                                                            min_unique=20)

            has_trend_g = centers_g.size >= 3

            if has_trend_g and np.sum(sigma_clip_mask_global & valid_corr_g) > 2:
                try:
                    rho_g, p_g = spearmanr(sfr_off_g[sigma_clip_mask_global & valid_corr_g],
                                         oh_off_g[sigma_clip_mask_global & valid_corr_g])
                    global_bin_centers.append(sigmaM_centers_g[i])
                    global_rho_values.append(rho_g)
                    print(f"  Global bin {i+1:02d} ({sigmaM_centers_g[i]:.2f}): n={n_pix_g}, r={rho_g:.3f} (p={p_g:.3g})")
                except Exception:
                    print(f"  Global bin {i+1:02d} ({sigmaM_centers_g[i]:.2f}): n={n_pix_g}, r=N/A (error)")
            else:
                print(f"  Global bin {i+1:02d} ({sigmaM_centers_g[i]:.2f}): n={n_pix_g}, r=N/A (no trend or too few in 3σ)")

    global_bin_centers = np.array(global_bin_centers)
    global_rho_values = np.array(global_rho_values)

    # ------------------------------------------------------------------
    # Plot Spearman rho vs Σ_* for all galaxies, with galaxy color legend
    # and an additional combined (all galaxies) trend in black
    # ------------------------------------------------------------------
    if len(all_rho_values) > 0:
    
        galaxy_indices_unique = np.unique(all_galaxy_indices)
    
        # Create figure with main axis and a thin top axis for legend
        fig = plt.figure(figsize=(10, 7))
        gs = gridspec.GridSpec(2, 1, height_ratios=[0.18, 1.0], hspace=0.05)
        legend_ax = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[1, 0])
        legend_ax.axis('off')
    
        plotted_galaxies = [] # Keep track of which galaxies are actually plotted
    
        # Plot Spearman rho vs Σ_* for each galaxy
        for idx in galaxy_indices_unique:
            mask = all_galaxy_indices == idx
            gal_name = galaxies[idx]
            color = galaxy_color_dict[gal_name]
        
            # Sort by Σ_* so lines are not scrambled
            x_vals = all_bin_centers[mask]
            y_vals = all_rho_values[mask]
        
            # Filter: skip if only 1 data point
            if len(x_vals) <= 2:
                continue
            
            plotted_galaxies.append(gal_name)
        
            order = np.argsort(x_vals)
            x_sorted = x_vals[order]
            y_sorted = y_vals[order]
        
            ax.plot(x_sorted, y_sorted, color=color, linewidth=1.5, alpha=0.8)
            # No per-galaxy scatter points: line-only trend summary for this 1D diagnostic.
    
        # Add combined trend (all galaxies together) in black (no legend entry)
        if global_rho_values.size > 2:
            order_g = np.argsort(global_bin_centers)
            ax.plot(global_bin_centers[order_g], global_rho_values[order_g],
                    color='black', linewidth=2.0, alpha=0.9, linestyle='-')
            # No global scatter points: black line is the pooled trend summary.
    
        # Axis formatting
        ax.axhline(0.0, color='gray', linestyle='--', linewidth=1.0, alpha=0.5)
        ax.set_xlabel(r'$\log\,\Sigma_* \; (M_\odot\,\mathrm{kpc}^{-2})$', fontsize=12)
        ax.set_ylabel(r'Spearman $\rho$ [Δ[O/H] vs Δ$\log\,\Sigma_{\mathrm{SFR}}$]', fontsize=12)
        # ax.set_title(r'Spearman $\rho$ vs $\log\,\Sigma_*$ (all galaxies, excluding NGC4383)', fontsize=14)
        ax.grid(True, alpha=0.3)
    
        # ------------------------------------------------------------------
        # Galaxy legend at the top using dummy patches (excluding NGC4383)
        # ------------------------------------------------------------------
        legend_elements = [
            Patch(facecolor=galaxy_color_dict[gal], edgecolor='black', label=gal)
            for gal in plotted_galaxies
        ]
        if len(plotted_galaxies) > 0:
            legend_ax.legend(handles=legend_elements, loc='center', ncol=min(len(plotted_galaxies), 6),
                            frameon=False, fontsize=11, title='Galaxies', title_fontsize=13)
    
        plt.tight_layout()
        # plt.savefig('Plot_Paper_1/MAUVE_MUSE_SpearmanRho_vs_SigmaM_AllGalaxies_ExclNGC4383_EW.png',
                    # dpi=600, bbox_inches='tight')
        # plt.savefig('Plot_Paper_1/MAUVE_MUSE_SpearmanRho_vs_SigmaM_AllGalaxies_ExclNGC4383_EW.pdf',
                    # bbox_inches='tight')
        plt.show()
    else:
        print("No valid Spearman rho values to plot.")
else:
    print("Skipped copied previous-project analysis cell 6: no nested *_further mass/SFR/OH+EW products are ready yet.")
