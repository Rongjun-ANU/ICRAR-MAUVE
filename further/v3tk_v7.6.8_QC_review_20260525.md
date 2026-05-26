# v3tk v7.6.8 Initial QC Review

Date: 2026-05-25

## Scope Checked

- Product root: `/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8/{GALID}`
- Expected QC PDFs: `/Users/Igniz/Desktop/ICRAR/further/{GALID}_v3tk_v7.6.8_QC.pdf`
- Found 26 galaxy product folders.
- Found 16 expected-location QC PDFs, each with 20 pages.
- Treating NGC4383 and NGC4419 as previous-run products and excluding them from the new-product science review.
- Active new-product review set: 16 galaxies with QC PDFs.

Supporting generated review files:

- Inventory: `.planning/qc_review_20260525/inventory.csv`
- HDU inventory: `.planning/qc_review_20260525/hdu_inventory.csv`
- Map summary metrics: `.planning/qc_review_20260525/map_metrics.csv`
- Derived QC metrics: `.planning/qc_review_20260525/derived_qc_metrics.csv`
- PDF page contact sheets: `.planning/qc_review_20260525/page_contact_sheets/`

## Output Completeness Flags

Expected QC PDFs are missing for active new-product folders:

- NGC4293, NGC4298, NGC4302, NGC4330, NGC4396, NGC4457, NGC4567_8, NGC4698.

These currently have only KIN/mask/spatial-binning products and lack gas/SFH products.

Excluded previous-run products:

- NGC4383 and NGC4419 have different 21-line gas products including ARIII7135, HEI5875, NI5197, NI5200, and SIII9068. I do not use them below.

## QC Script Issues

`QC_ngist_v3tk_v768.py` has two sign-label mismatches:

- The panel labelled `Vs-Vha` is computed as `HA6562_VEL - stellar V`.
- The panel labelled `Sigma_s - Sigma_ha` is computed as `HA6562_SIGMA - stellar SIGMA`.

The maps are still useful, but the labels currently imply the opposite sign.

## Scientific/QC Flags

H-beta is the weak link for extinction and BPT-style ratios. Median HB4861 Flux/Err is only 3.0 for NGC4606, 7.3 for NGC4694, 7.9 for NGC4580, 8.3 for NGC4064, 9.4 for IC3392, and 9.8 for NGC4405. H-alpha is generally much stronger.

Balmer decrement is suspicious in several products. Among pixels with Halpha and Hbeta Flux/Err > 5, the fraction with Halpha/Hbeta < 2.86 is 65.7% for NGC4064, 59.9% for NGC4405, 55.8% for NGC4580, 51.3% for NGC4606, and 48.5% for IC3392. Check stellar absorption/continuum subtraction/calibration before using these as extinction maps.

SII6716/SII6730 needs stronger filtering. Even at S/N > 5 in both lines, many pixels are outside the nominal 0.44-1.45 range: IC3392 67.3%, NGC4351 65.5%, NGC4064 56.7%, NGC4394 53.2%, NGC4694 52.0%, and NGC4405 51.5%.

Stellar sigma is often close to or below the practically useful resolution floor. Fraction of KIN sigma pixels below 50 km/s is 99.4% in NGC4405, 89.7% in NGC4351, 89.0% in IC3392, 86.7% in NGC4606, and 84.8% in NGC4522. Treat stellar dispersion maps cautiously unless resolution corrections are explicitly validated.

H-alpha versus stellar velocity offsets are scientifically interesting, but should be cut on H-alpha S/N. At H-alpha Flux/Err > 10, |V_stars - V_Halpha| > 50 km/s occurs in 68.4% of NGC4388, 35.2% of NGC4402, 22.1% of NGC4405, 18.8% of NGC4064, 17.8% of IC3392, and 15.7% of NGC4522.

SFH AGE/METAL maps often look much noisier/blockier than EBV. EBV has coherent dust morphology in many cases, while AGE/METAL can be dominated by binning/noise or template degeneracy. Interpret AGE/METAL only after checking S/N, chi2, and SFH-weight stability.

## Immediate Science Opportunities

1. Gas-stripping / extraplanar ionized-gas kinematics. Highest immediate targets: NGC4388, NGC4402, NGC4522, and NGC4607. The data already show coherent H-alpha structure plus large gas-star velocity residuals or edge-on disturbed morphology. A first notebook could rank asymmetry using H-alpha flux, Vstar-Vha, and sigma offsets with H-alpha S/N cuts.
2. A "kinematically disturbed Virgo disks" mini-sample. NGC4388 and NGC4402 are the strongest from the velocity residual metric; NGC4522 and NGC4607 are natural edge-on/stripping candidates. This is a compact science story: where does ionized gas stop following stellar rotation, and is the offset aligned with stripped/asymmetric gas?
3. Dust/continuum residual diagnostic from Balmer decrement failures. NGC4064, NGC4405, NGC4580, NGC4606, and IC3392 have large fractions of Halpha/Hbeta below 2.86 even with S/N>5. This is not a final dust result yet, but it is immediately useful as a project on H-beta absorption/continuum subtraction systematics versus EBV morphology.
4. Resolved shock/LIER candidate mapping in disturbed systems. Use S/N-masked NII/Ha, SII/Ha, OI/Ha, and gas velocity dispersion. Good first targets are NGC4388, NGC4694, NGC4064, NGC4189, and NGC4501, where morphology/line ratios look nontrivial enough to justify a focused map set.
5. Star-forming ring/spiral structure in NGC4501. This looks like the cleanest "pretty and immediate" disk case: coherent H-alpha and EBV spiral/ring morphology, good H-alpha S/N, and enough structure for radial/azimuthal profiles or comparison of gas, dust, and stellar-population maps.

## Next Analyses Worth Doing

1. Environmental gas-stripping kinematics: use H-alpha morphology plus Vstar-Vgas residuals for NGC4388, NGC4402, NGC4522, and NGC4607. These look like the highest-yield active new-product cases for extraplanar/asymmetric gas.
2. Resolved ionization classification: build S/N-masked BPT maps and compare NII/Ha, SII/Ha, OI/Ha where reliable. Separate star-forming disks, shocks/LIER-like regions, and AGN/outflow candidates.
3. Dust and continuum-systematics check: focus on galaxies with Halpha/Hbeta below intrinsic values. Compare Balmer decrement against EBV maps, residual spectra, and Hbeta absorption regions.
4. Stellar-population sanity workflow: do not use AGE/METAL maps alone. Compare SFH-weight grids, RED_CHI2, postfit S/N, and mass/light-weighted alternatives before interpreting gradients.
5. Sample-level product readiness dashboard: one row per galaxy with product completeness, PDF status, line-set family, median Halpha/Hbeta/OIII/SII S/N, and recommended science-safe products.
