# Completed Setonix three-method R-band M/L experiment

Slurm job 51895534 completed successfully on 2026-10-31 at 14:12 UTC
(22:12 Australia/Perth), after approximately 4 hours 53 minutes.
All 51,609 fresh SFH fits succeeded; no bins failed.

Final notebook:
`../20261031_SETONIX_MAUVE_ML_Three_Methods_full_executed.ipynb`

Remote notebook:
`/home/rhuang/20261031_SETONIX_MAUVE_ML_Three_Methods_full_executed.ipynb`

Remote experiment:
`/scratch/pawsey1308/rhuang/mauve_ml_three_methods_20261031/d331b6665d204432/`

This local directory contains downloaded summaries, per-bin CSVs, manifests and
figures. Large fit HDF5 checkpoints, the SPS NPZ and output FITS maps remain in
the remote experiment directory.

## Plot refresh on 1 November 2026

The source and full executed three-method notebooks now include dedicated
raw-count M/L histograms and pairwise scatter panels with explicit 1:1 lines.
These were regenerated on Setonix from the completed full-fit CSVs; the SFH
fits, weights, numerical summaries and earlier notebook outputs were preserved.
The executed notebook has 31 cells, with separate plotting-refresh provenance.
Pre-refresh source/executed notebook copies are retained in
`notebook_archive/pre_plot_refresh_20261101/`.

For each of the four cases, `figures/` now also contains:

- `ml_histograms_raw_counts.png`: three overlaid distributions using identical
  60-bin edges, the full M/L range, and raw spatial-bin counts (no density normalization).
- `ml_pairwise_scatter_one_to_one.png`: all three method pairs with equal x/y
  limits and aspect, dashed 1:1 lines, and ratio-direction annotations; all bins
  are plotted without subsampling.
- `ml_histogram_counts.csv`: reproducible histogram counts and common edges.

Every method's histogram count sums to 2582, 2593, 22547 or 23887, respectively.
`plot_refresh_20261101.json` records input CSV hashes and plot provenance. The
updated Markdown/PDF report includes all eight new figures. The original
`comparisons.png` panels remain complementary difference/scatter diagnostics.

## Results

Define legacy error as `100*(light_r/light_norm_to_r - 1)`.
Negative values mean that the legacy interpretation underestimates R-band M/L.
Each fitted bin has equal weight in these statistics. These are not
luminosity-weighted galaxy-integrated M/L or mass corrections.

| Galaxy | Range (A) | Successful bins | Median legacy error (%) | P99 absolute error (%) | Maximum absolute error (%) | Worst BINID |
|---|---|---:|---:|---:|---:|---:|
| NGC4522 | 4800-7000 | 2582 | -4.1241 | 13.3597 | 17.5524 | 464 |
| NGC4522 | 4800-8900 | 2593 | -1.9173 | 5.2364 | 6.2577 | 207 |
| NGC4698 | 4800-7000 | 22547 | -1.1769 | 7.0311 | 11.7622 | 22223 |
| NGC4698 | 4800-8900 | 23887 | -0.7926 | 3.3881 | 5.4599 | 14455 |

Define corrected/native difference as
`100*(light_norm_to_r/ppxf_native - 1)`.

| Galaxy | Range (A) | Median difference (%) | P99 absolute difference (%) | Maximum absolute difference (%) | Worst BINID |
|---|---|---:|---:|---:|---:|
| NGC4522 | 4800-7000 | +0.9815 | 1.8248 | 2.3945 | 1069 |
| NGC4522 | 4800-8900 | +0.9608 | 1.6217 | 3.1494 | 2249 |
| NGC4698 | 4800-7000 | +0.4248 | 1.2253 | 1.8423 | 22460 |
| NGC4698 | 4800-8900 | +0.6029 | 1.1563 | 1.9650 | 23364 |

The largest direct light_r/native discrepancy was 16.5797%, in NGC4522's
7000-A BINID 464. All three pairwise comparisons are in
`all_cases_summary_full.csv` and the per-case `bin_results.csv` files.

## Interpretation and checks

- Original nGIST versus current Mass.py preparation passed for all 477 SSPs
  in each of the four cases. Prepared and LIGHT-normalized spectral samples
  matched exactly in the checks; divisor comparisons differed by less than
  1e-14 in relative terms.
- The corrected table-based M/L and reference-divisor table-control expression
  agreed across every bin, with maximum relative difference 5.2e-15 in the
  saved bin-result comparisons.
- A separate read-only audit of every saved fresh weight vector reconstructed
  both the table M/L and the native spectral-photometry M/L from the shared SSP
  masses and the respective R luminosities. Agreement was better than 3e-13
  relative after reading the diagnostic CSVs. See `verification_report.json`.
  This accounts for the corrected/native residual through table-versus-spectral
  R luminosities under the recorded conventions, rather than normalization
  reconstruction failure.
- The native route called the actual unmodified pPXF helper on fresh fits.
  SSP masses were derived from the matched BaSTI/Chabrier table. All routes
  therefore share that mass calibration; this is not an independent validation
  of the IMF, stellar masses or table calibration. The adopted Vega R,
  solar-magnitude and spectral-unit conventions remain explicit assumptions.
- The legacy bias is larger for the existing 7000-A configuration in both
  galaxies. MDEG also differs (12 versus 23), so the comparison does not isolate
  wavelength coverage alone. Corrected/native P99 differences are similar
  between configurations; their literal maxima are larger in the 8900-A cases.
- NGC4522 8900-A BINID 2249, the maximum corrected/native discrepancy of 3.15%,
  has post-fit S/N 10.6 and prominent spectral artifacts in its diagnostic plot.
  Retain this literal maximum but also report P99; a finite successful fit does
  not itself certify science-quality weights. Other maxima are retained without
  an additional post-hoc S/N selection.

## Verification and artifacts

Downloaded notebook schema was validated and contains no error outputs.
Per-bin CSV counts, positive finite M/L values, all pairwise medians/maxima and
worst-bin IDs were independently checked against the summaries. The range
comparison and representative worst-bin fit figures were visually inspected.
No refits or production products were changed during this final audit.

Useful files:

- `all_cases_summary_full.csv`: all pairwise statistics.
- `normalization_reproduction_summary.csv` and `.json`: four-case preparation test.
- `ssp_photometry_comparison.csv`: table/spectral R luminosities and shared masses.
- `verification_report.json`: independent saved-weight reconstruction audit.
- `range_summary_full.png`: P99 and maximum deviations by fitting configuration.
- `full/<case>/bin_results.csv`: all bin-level values and fit diagnostics.
- `full/<case>/figures/`: maps, comparisons and worst-bin spectra/SSP contributions.
