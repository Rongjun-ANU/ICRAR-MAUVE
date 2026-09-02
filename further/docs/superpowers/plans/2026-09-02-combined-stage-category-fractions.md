# Combined Stage Category Fractions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create one dated notebook that measures SF, Balmer non-detection (ND), and Balmer-detected non-SF (NSF) fractions versus stellar-mass surface density and normalized radius for pre-peak, close-to-peak, post-peak, and all qualifying galaxies.

**Architecture:** Reuse the live category and valid-disc definitions from the three 2026-09-01 notebooks, but replace individual-galaxy plotting with a common-bin aggregation table. Build two estimators from the same per-galaxy counts: pooled spaxel counts and an equal-galaxy mean of within-bin fractions.

**Tech Stack:** Python, Jupyter notebook, NumPy, pandas, Matplotlib, Astropy FITS/WCS.

---

### Task 1: Build the dated notebook and live QC inventory

**Files:**
- Create: `20260902_check_combined_SF_fraction_categories_by_stage.ipynb`
- Reference only: the three `20260901_check_*peak_SF_fraction_and_rSFMS_categories.ipynb` notebooks

- [ ] **Step 1: Reuse the current map and mask definitions**

Preserve `valid_disc_mask`, the unchanged SF selection, and the Balmer-only partition:

```python
balmer_detected = detected["HA6562"] & detected["HB4861"]
nd = valid & ~balmer_detected
nsf = valid & balmer_detected & ~sf
assert np.array_equal(sf.astype(np.int8) + nd + nsf, valid.astype(np.int8))
```

- [ ] **Step 2: Discover every inclination-<80 degree product in stage classes 1, 2, and 3**

Load each product once, report pass/fail QC by galaxy and stage, and stop if any stage has no passing galaxy.

- [ ] **Step 3: Verify source parity**

Compare the reused constants and mask expressions against all three source notebooks and confirm that stage membership comes from the live catalogue and geometry files.

### Task 2: Compute common-bin counts and both estimators

**Files:**
- Modify: `20260902_check_combined_SF_fraction_categories_by_stage.ipynb`

- [ ] **Step 1: Construct common complete bin edges**

Use 1/3-dex bins for `log Sigma_star` and 0.2-wide bins for `R/Re`, aligned to the live finite valid-disc range so no valid pixel is silently dropped.

- [ ] **Step 2: Store per-galaxy counts in every bin**

```python
row = {
    "GALID": galid, "stage": stage, "variable": variable,
    "bin_left": lo, "bin_right": hi, "center": (lo + hi) / 2,
    "N_valid": int((valid & in_bin).sum()),
    "N_SF": int((sf & in_bin).sum()),
    "N_ND": int((nd & in_bin).sum()),
    "N_NSF": int((nsf & in_bin).sum()),
}
assert row["N_SF"] + row["N_ND"] + row["N_NSF"] == row["N_valid"]
```

- [ ] **Step 3: Compute pooled-spaxel fractions**

For each stage and the all-stage union, sum category and valid counts first, then calculate `f_category = sum(N_category) / sum(N_valid)`.

- [ ] **Step 4: Compute equal-galaxy fractions**

Within each variable bin, retain galaxies with at least 50 valid pixels, calculate `f_category,g = N_category,g / N_valid,g`, and average those galaxy fractions. This is equivalent to assigning each usable pixel in galaxy `g` weight `1/N_valid,g` within that bin. Store the contributing-galaxy count without adding percentile shading.

- [ ] **Step 5: Assert accounting and weighting**

Require the three category fractions to sum to one in every supported pooled and equal-galaxy bin; require every included galaxy's within-bin pixel weights to sum to one.

### Task 3: Plot, explain, and verify

**Files:**
- Modify: `20260902_check_combined_SF_fraction_categories_by_stage.ipynb`

- [ ] **Step 1: Plot the pooled estimator**

Create a 2-by-4 figure: rows are `log Sigma_star` and `R/Re`; columns are pre-peak, close-to-peak, post-peak, and all stages. Plot SF, ND, and NSF together with fixed category colours and show contributing valid-pixel counts.

- [ ] **Step 2: Plot the equal-galaxy estimator**

Use the same layout and axes with contributing-galaxy counts and no shaded regions. State explicitly that weights are recomputed within every displayed bin.

- [ ] **Step 3: Add compact output tables and interpretation limits**

Display stage-level galaxy/pixel totals plus the pooled and equal-galaxy result tables. Explain that pooled curves describe pixels, equal-galaxy curves describe the mean galaxy, NSF is a selection outcome rather than proof of absent star formation, and neither estimator establishes causality.

- [ ] **Step 4: Execute and verify the notebook**

Run the notebook against live files. Require zero execution errors, exact category accounting, all three stages represented, complete global-bin coverage, and per-bin equal-galaxy weight totals of one.

### Task 4: Add the all-galaxy radial-cut sensitivity figures and spatial maps

**Files:**
- Modify: `20260902_check_combined_SF_fraction_categories_by_stage.ipynb`

- [ ] **Step 1: Preserve the original 26-product results**

Keep `POOLED_FRACTIONS` and `EQUAL_GALAXY_FRACTIONS` unchanged for the original two figures.

- [ ] **Step 2: Recompute both estimators after the targeted cut**

For every galaxy, retain only valid spaxels with `R/Re < 4` before recounting both the stellar-density and radial bins. Independently rebuild the pooled-spaxel and equal-galaxy tables and verify stage counts remain 8 pre-peak, 5 close-to-peak, 13 post-peak, and 26 all-stage products.

- [ ] **Step 3: Plot both sensitivity versions**

Append one 2-by-4 pooled-spaxel figure and one 2-by-4 equal-galaxy figure with the radial threshold identified explicitly. Keep both figures free of shaded regions and use one common x-axis range across all four panels in each row.

- [ ] **Step 4: Re-execute and verify**

Require exact SF/ND/NSF accounting, unit within-bin galaxy weights, expected sample membership, four combined figures, and no stored errors or warnings.

- [ ] **Step 5: Append individual spatial category maps**

At the notebook end, reproduce the source WCS SF/ND/NSF category-map style for every passing galaxy. Overlay a dotted black `1 Re` contour and a solid black `4 Re` contour, retaining the unchanged valid-disc category footprint.
