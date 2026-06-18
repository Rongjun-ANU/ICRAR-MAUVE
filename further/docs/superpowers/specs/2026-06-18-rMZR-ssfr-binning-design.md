# Design Spec: sSFR-based Offset Calculation and Spearman Trends

Transitioning the resolved Mass-Metallicity Relation (rMZR) notebook [20260616_initial_check_rMZR_binned_by_EWHa.ipynb](file:///Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb) from binning offsets by $\Sigma_*$ (stellar mass surface density) to specific Star Formation Rate (sSFR).

## User Review Required

> [!IMPORTANT]
> Cell indexing: The user uses 1-based indexing. Cell 8 is Python-index 7, Cell 22 is Python-index 21, Cell 23 is Python-index 22, Cell 27 is Python-index 26, Cell 28 is Python-index 27, and Cell 29 is Python-index 28. All cell references in this design document use the user's 1-based indexing.

---

## Proposed Changes

### [20260616_initial_check_rMZR_binned_by_EWHa.ipynb](file:///Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb)

#### [MODIFY] Reusable Helpers (Cell 8)
- Define a global constant: `XRANGE_SSFR = (-11.2, -8.6)`.
- Update `plot_per_galaxy_spearman_by_mass_for_color` (renaming it to `plot_per_galaxy_spearman_by_ssfr_for_color` or updating its inline logic):
  - Calculate specific star formation rate `log_ssfr = data["logSigmaSFR"] - data["logSigmaM"]`.
  - Pass `log_ssfr` as the primary binning axis to `spearman_rho_by_primary_bin_for_color`.
  - Update the plot's x-limits to `ax.set_xlim(*XRANGE_SSFR)` and the x-label to:
    $$ax.set\_xlabel(r'\log\,\mathrm{sSFR} \ (\mathrm{yr}^{-1})', fontsize=12)$$

#### [MODIFY] Spearman Summary vs sSFR (Cell 22)
- For each galaxy, calculate specific star formation rate:
  $$log\_ssfr = sigSFR - sigM$$
- Replace stellar mass bin edges (`sigmaM_bin_edges`) with specific star formation rate bin edges:
  $$ssfr\_min = np.nanpercentile(log\_ssfr[valid\_all], 2.5)$$
  $$ssfr\_max = np.nanpercentile(log\_ssfr[valid\_all], 97.5)$$
  $$ssfr\_bin\_edges = np.linspace(ssfr\_min, ssfr\_max, n\_bins + 1)$$
- Update offset subtraction to use the mean values in each sSFR bin.
- For the combined global trend, calculate `ssfr_min_g` and `ssfr_max_g` on `global_ssfr` and use 12 bins.
- Update summary plot x-axis range to `(-11.2, -8.6)` and update all titles, labels, and legends to refer to sSFR.

#### [MODIFY] Offset Contours (Cell 23)
- Update offset calculations to bin spaxels by `log_ssfr` instead of `sigM` (using 6 dynamic percentile bins per galaxy).
- Ensure generated contours are plotted against these sSFR-derived offsets.

#### [MODIFY] Moran-like 2x2 Maps (Cells 27 & 28)
- Update offset map calculations to bin by specific star formation rate (`log_ssfr`) with 12 dynamic percentile bins.
- In the bottom-right panel, plot the binned sSFR map (`ssfr_binned_map`) instead of the binned $\Sigma_*$ map.
- Update the bottom-right colorbar and ticks to represent:
  $$\log\,\mathrm{sSFR} \ (\mathrm{yr}^{-1})$$

#### [MODIFY] Combined Figure with NGC4383 (Cell 29)
- Modify `compute_spearman_track` and `compute_spearman_track_from_offsets` to bin by sSFR (`log_ssfr`) instead of stellar mass.
- Update columns 3 and 4 of the combined figure (the Spearman tracks) to plot against sSFR with x-limits of `(-11.2, -8.6)`.
- Update corresponding titles and axis labels.

---

## Verification Plan

### Manual Verification
- We will execute the updated cells inside the IPython notebook against the real local products.
- We will check that the output figures display correctly without any `ValueError` or shape mismatch issues, and confirm that the bottom-right panel of the Moran-like maps correctly reflects sSFR.
