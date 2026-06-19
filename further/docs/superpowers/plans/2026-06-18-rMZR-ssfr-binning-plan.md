# rMZR sSFR Binning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Transition the rMZR notebook from binning offsets and Spearman trends by $\Sigma_*$ (stellar mass surface density) to sSFR (specific star formation rate).

**Architecture:** Update cells 8, 22, 23, 27, 28, and 29 of the notebook to compute sSFR ($\log \Sigma_{\mathrm{SFR}} - \log \Sigma_*$), use it as the primary binning axis for offset mean subtraction and Spearman trends, update x-axis limits and labels, and replace the bottom-right Moran map helper plot with binned sSFR.

**Tech Stack:** Python 3, IPython notebook (`.ipynb`), `numpy`, `scipy`, `matplotlib`, `astropy.io.fits`.

## Global Constraints
- Target Path: `/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb`
- Cell Indexing: All references to cells are 1-based.
- Python environment: `/opt/miniconda3/envs/ICRAR/bin/python`.
- Limit changes strictly to the requested binning transition. Do not refactor unrelated code.

---

### Task 1: Update Reusable Helpers (Cell 8)

**Files:**
- Modify: `/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb` (Cell 8, index 7)

**Interfaces:**
- Consumes: `data["logSigmaSFR"]`, `data["logSigmaM"]`
- Produces: `plot_per_galaxy_spearman_by_ssfr_for_color` and `XRANGE_SSFR`

- [ ] **Step 1: Write a python check script to verify code cell content and structure**
  We will write a python script `/Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/test_cell8.py` that checks the current helper definitions.
  ```python
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
  cell = nb['cells'][7]
  src = "".join(cell['source'])
  assert "def plot_per_galaxy_spearman_by_mass_for_color" in src
  print("Cell 8 helper exists and contains expected function.")
  ```

- [ ] **Step 2: Run verification script**
  Run: `/opt/miniconda3/envs/ICRAR/bin/python /Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/test_cell8.py`
  Expected: Success output.

- [ ] **Step 3: Modify Cell 8 (index 7)**
  Write a script `/Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cell8.py` to replace the helper function:
  ```python
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
  cell = nb['cells'][7]
  src = "".join(cell['source'])
  
  old_code = """def plot_per_galaxy_spearman_by_mass_for_color(
    ax,
    per_galaxy,
    galaxies,
    *,
    linewidth=1.6,
    alpha=0.85,
    min_bins_to_plot=3,
):
    plotted = []
    for gal in galaxies:
        if gal not in per_galaxy:
            continue
        data = per_galaxy[gal]
        centers, rho_values, p_values, counts = spearman_rho_by_primary_bin_for_color(
            data["logSigmaM"],
            data["logSigmaSFR"],
            data["oh_m13"],
        )
        if rho_values.size < min_bins_to_plot:
            print(f"  {gal}: Spearman bins plotted=0 (valid bins={rho_values.size})")
            continue
        order = np.argsort(centers)
        ax.plot(
            centers[order],
            rho_values[order],
            color=galaxy_color(gal),
            linewidth=linewidth,
            alpha=alpha,
            label=gal,
        )
        plotted.append(gal)
        print(f"  {gal}: Spearman bins plotted={rho_values.size}")

    ax.axhline(0.0, color="0.45", linestyle="--", linewidth=1.0, alpha=0.65)
    ax.set_xlim(*XRANGE_RMZR)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel(r'$\log\,\Sigma_* \ (M_\odot\,\mathrm{kpc}^{-2})$', fontsize=12)
    ax.set_ylabel(r'Spearman $\\rho$ [Δ[O/H] vs Δ$\\log\\,\\Sigma_\\mathrm{SFR}$]', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.text(
        0.02,
        0.98,
        r"Per-galaxy Spearman $\\rho$ in $\\Sigma_*$ bins",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )
    return plotted"""

  new_code = """XRANGE_SSFR = (-11.2, -8.6)

def plot_per_galaxy_spearman_by_ssfr_for_color(
    ax,
    per_galaxy,
    galaxies,
    *,
    linewidth=1.6,
    alpha=0.85,
    min_bins_to_plot=3,
):
    plotted = []
    for gal in galaxies:
        if gal not in per_galaxy:
            continue
        data = per_galaxy[gal]
        log_ssfr = data["logSigmaSFR"] - data["logSigmaM"]
        centers, rho_values, p_values, counts = spearman_rho_by_primary_bin_for_color(
            log_ssfr,
            data["logSigmaSFR"],
            data["oh_m13"],
        )
        if rho_values.size < min_bins_to_plot:
            print(f"  {gal}: Spearman bins plotted=0 (valid bins={rho_values.size})")
            continue
        order = np.argsort(centers)
        ax.plot(
            centers[order],
            rho_values[order],
            color=galaxy_color(gal),
            linewidth=linewidth,
            alpha=alpha,
            label=gal,
        )
        plotted.append(gal)
        print(f"  {gal}: Spearman bins plotted={rho_values.size}")

    ax.axhline(0.0, color="0.45", linestyle="--", linewidth=1.0, alpha=0.65)
    ax.set_xlim(*XRANGE_SSFR)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel(r'$\\log\\,\\mathrm{sSFR} \\ (\\mathrm{yr}^{-1})$', fontsize=12)
    ax.set_ylabel(r'Spearman $\\rho$ [Δ[O/H] vs Δ$\\log\\,\\Sigma_\\mathrm{SFR}$]', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.text(
        0.02,
        0.98,
        r"Per-galaxy Spearman $\\rho$ in sSFR bins",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )
    return plotted"""

  if old_code in src:
      src = src.replace(old_code, new_code)
  else:
      # Try normalized spacing
      raise ValueError("Exact old code block not found in Cell 8.")
      
  cell['source'] = [line + '\n' for line in src.split('\n')][:-1]
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb', 'w') as f:
      json.dump(nb, f, indent=1)
  print("Cell 8 helper function updated successfully.")
  ```
  Run: `/opt/miniconda3/envs/ICRAR/bin/python /Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cell8.py`

- [ ] **Step 4: Verify cell 8 syntax**
  Verify that the modified notebook cell runs without syntax error.
  Run: `/opt/miniconda3/envs/ICRAR/bin/python -c "import json; nb=json.load(open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb')); exec(''.join(nb['cells'][7]['source']))"`
  Expected: Success output without errors.

- [ ] **Step 5: Commit**
  ```bash
  git add 20260616_initial_check_rMZR_binned_by_EWHa.ipynb
  git commit -m "feat: update Cell 8 helper function to bin by sSFR"
  ```

---

### Task 2: Update Spearman Summary vs sSFR (Cell 22)

**Files:**
- Modify: `/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb` (Cell 22, index 21)

- [ ] **Step 1: Write verification script to check Cell 22 contents**
  Write a script `/Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cell22.py` that loads, modifies, and saves Cell 22 code to use specific star formation rate `log_ssfr = sigSFR - sigM` for binning, updates the axis labels, and updates limits to `(-11.2, -8.6)`.
  *Note: The script will parse the cells JSON, find `cells[21]`, replace all instances of `sigmaM_bin_edges`, `sigmaM_bin_centers`, `sigmaM_min`, `sigmaM_max` with sSFR counterparts.*
  Let's write this script completely.
  ```python
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
  cell = nb['cells'][21]
  src = "".join(cell['source'])
  
  # Replacement 1: define log_ssfr
  src = src.replace(
      "        global_gal_index_list.append(np.full(np.sum(valid_all), gal_index, dtype=int))",
      "        global_gal_index_list.append(np.full(np.sum(valid_all), gal_index, dtype=int))\n        log_ssfr = sigSFR - sigM"
  )
  
  # Replacement 2: binning edges
  old_bin_def = """        # 1. Define 12 Σ_* bins for this galaxy (spanning 95% population)
        n_bins = 6
        sigmaM_min = np.nanpercentile(sigM[valid_all], 2.5)
        sigmaM_max = np.nanpercentile(sigM[valid_all], 97.5)
        sigmaM_bin_edges = np.linspace(sigmaM_min, sigmaM_max, n_bins + 1)
        sigmaM_bin_centers = 0.5 * (sigmaM_bin_edges[:-1] + sigmaM_bin_edges[1:])"""
        
  new_bin_def = """        # 1. Define 6 sSFR bins for this galaxy (spanning 95% population)
        n_bins = 6
        ssfr_min = np.nanpercentile(log_ssfr[valid_all], 2.5)
        ssfr_max = np.nanpercentile(log_ssfr[valid_all], 97.5)
        ssfr_bin_edges = np.linspace(ssfr_min, ssfr_max, n_bins + 1)
        ssfr_bin_centers = 0.5 * (ssfr_bin_edges[:-1] + ssfr_bin_edges[1:])"""
  src = src.replace(old_bin_def, new_bin_def)

  # Replacement 3: loop variables
  src = src.replace(
      "(sigM >= sigmaM_bin_edges[i])",
      "(log_ssfr >= ssfr_bin_edges[i])"
  ).replace(
      "< sigmaM_bin_edges[i+1])",
      "< ssfr_bin_edges[i+1])"
  ).replace(
      "<= sigmaM_bin_edges[i+1])",
      "<= ssfr_bin_edges[i+1])"
  )
  
  # Replacement 4: prints
  src = src.replace(
      'print(f"  Σ* range (95%): {sigmaM_min:.3f} to {sigmaM_max:.3f}")',
      'print(f"  sSFR range (95%): {ssfr_min:.3f} to {ssfr_max:.3f}")'
  ).replace(
      'print("  Spearman rho in each Σ_* bin (3σ clip + trend, per galaxy):")',
      'print("  Spearman rho in each sSFR bin (3σ clip + trend, per galaxy):")'
  ).replace(
      'sigmaM_bin_centers[i]',
      'ssfr_bin_centers[i]'
  ).replace(
      'all_bin_centers.append(sigmaM_bin_centers[i])',
      'all_bin_centers.append(ssfr_bin_centers[i])'
  )
  
  # Replacement 5: global list collection
  src = src.replace(
      'global_sigM_list.append(sigM[bin_mask][valid_corr])',
      'global_sigM_list.append(log_ssfr[bin_mask][valid_corr])'
  )
  
  # Replacement 6: global binning
  old_global_bin = """    if global_sigM_list:
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

    print("\\n" + "=" * 80)
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

        print("\\nComputing combined Spearman rho in Σ_* bins (all galaxies together)...")
        for i in range(n_bins_global):
            if i == 0:
                bin_mask_g = (global_sigM >= sigmaM_edges_g[i]) & (global_sigM < sigmaM_edges_g[i+1])
            elif i == n_bins_global - 1:
                bin_mask_g = (global_sigM >= sigmaM_edges_g[i]) & (global_sigM <= sigmaM_edges_g[i+1])
            else:
                bin_mask_g = (global_sigM >= sigmaM_edges_g[i]) & (global_sigM < sigmaM_edges_g[i+1])"""
                
  new_global_bin = """    if global_sigM_list:
        global_ssfr = np.concatenate(global_sigM_list)
        global_sfr_off = np.concatenate(global_sfr_off_list)
        global_oh_off = np.concatenate(global_oh_off_list)
    else:
        global_ssfr = np.array([])
        global_sfr_off = np.array([])
        global_oh_off = np.array([])

    # ------------------------------------------------------------------
    # Convert collected per-galaxy Spearman rho data to arrays
    # ------------------------------------------------------------------
    all_bin_centers = np.array(all_bin_centers)
    all_rho_values = np.array(all_rho_values)
    all_galaxy_indices = np.array(all_galaxy_indices)

    print("\\n" + "=" * 80)
    print("Summary: number of sSFR bins with valid Spearman rho (per galaxy):", len(all_rho_values))
    print("=" * 80)

    # ------------------------------------------------------------------
    # Compute combined Spearman rho vs sSFR over all galaxies (same method)
    # ------------------------------------------------------------------
    global_bin_centers = []
    global_rho_values = []

    if global_ssfr.size > 0:
        # Define 12 sSFR bins for the combined sample (spanning 95% population)
        n_bins_global = 12
        ssfr_min_g = np.nanpercentile(global_ssfr, 2.5)
        ssfr_max_g = np.nanpercentile(global_ssfr, 97.5)
        ssfr_edges_g = np.linspace(ssfr_min_g, ssfr_max_g, n_bins_global + 1)
        ssfr_centers_g = 0.5 * (ssfr_edges_g[:-1] + ssfr_edges_g[1:])

        print("\\nComputing combined Spearman rho in sSFR bins (all galaxies together)...")
        for i in range(n_bins_global):
            if i == 0:
                bin_mask_g = (global_ssfr >= ssfr_edges_g[i]) & (global_ssfr < ssfr_edges_g[i+1])
            elif i == n_bins_global - 1:
                bin_mask_g = (global_ssfr >= ssfr_edges_g[i]) & (global_ssfr <= ssfr_edges_g[i+1])
            else:
                bin_mask_g = (global_ssfr >= ssfr_edges_g[i]) & (global_ssfr < ssfr_edges_g[i+1])"""
  src = src.replace(old_global_bin, new_global_bin)

  # Replacement 7: more global references
  src = src.replace("sfr_off_g = global_sfr_off[bin_mask_g]", "sfr_off_g = global_sfr_off[bin_mask_g]")
  src = src.replace("sigmaM_centers_g[i]", "ssfr_centers_g[i]")
  src = src.replace("global_bin_centers.append(sigmaM_centers_g[i])", "global_bin_centers.append(ssfr_centers_g[i])")
  
  # Replacement 8: plot formatting
  src = src.replace(
      "ax.plot(x_sorted, y_sorted, color=color, linewidth=1.5, alpha=0.8)",
      "ax.plot(x_sorted, y_sorted, color=color, linewidth=1.5, alpha=0.8)"
  )
  src = src.replace(
      "ax.set_xlabel(r'$\\log\\,\\Sigma_* \\; (M_\\odot\\,\\mathrm{kpc}^{-2})$', fontsize=12)",
      "ax.set_xlabel(r'$\\log\\,\\mathrm{sSFR} \\ (\\mathrm{yr}^{-1})$', fontsize=12)"
  ).replace(
      "r\"Per-galaxy Spearman $\\rho$ in $\\Sigma_*$ bins\"",
      "r\"Per-galaxy Spearman $\\rho$ in sSFR bins\""
  ).replace(
      "ax.set_xlim(6, 10)", # wait, let's look at the original code's set_xlim
      "ax.set_xlim(-11.2, -8.6)"
  ).replace(
      "ax.axhline(0.0, color='gray', linestyle='--', linewidth=1.0, alpha=0.5)",
      "ax.axhline(0.0, color='gray', linestyle='--', linewidth=1.0, alpha=0.5)\n        ax.set_xlim(-11.2, -8.6)"
  )

  cell['source'] = [line + '\n' for line in src.split('\n')][:-1]
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb', 'w') as f:
      json.dump(nb, f, indent=1)
  print("Cell 22 updated.")
  ```
  Run: `/opt/miniconda3/envs/ICRAR/bin/python /Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cell22.py`

- [ ] **Step 2: Run verification check on Cell 22**
  Verify cell runs successfully against actual FITS data.
  Run:
  ```bash
  /opt/miniconda3/envs/ICRAR/bin/python -c "
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
  exec(''.join(nb['cells'][2]['source']))
  exec(''.join(nb['cells'][3]['source']))
  exec(''.join(nb['cells'][21]['source']))
  "
  ```
  Expected: Success output showing Spearman rho values for each sSFR bin.

- [ ] **Step 3: Commit**
  ```bash
  git add 20260616_initial_check_rMZR_binned_by_EWHa.ipynb
  git commit -m "feat: update Cell 22 to calculate offsets and Spearman summary in sSFR bins"
  ```

---

### Task 3: Update Offset Contours (Cell 23)

**Files:**
- Modify: `/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb` (Cell 23, index 22)

- [ ] **Step 1: Write apply script for Cell 23**
  Write a script `/Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cell23.py` to modify cell 23 to use sSFR dynamic binning.
  ```python
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
  cell = nb['cells'][22]
  src = "".join(cell['source'])
  
  # Replacement 1: define log_ssfr
  src = src.replace(
      "        valid_all = valid_sfr & valid_oh & valid_sigM",
      "        valid_all = valid_sfr & valid_oh & valid_sigM\n        log_ssfr = sigSFR - sigM"
  )
  
  # Replacement 2: binning definitions
  old_bin_def = """        # Define 12 Σ_* bins for this galaxy (spanning 95% population)
        n_bins = 6
        sigmaM_min = np.nanpercentile(sigM[valid_all], 2.5)
        sigmaM_max = np.nanpercentile(sigM[valid_all], 97.5)
        sigmaM_bin_edges = np.linspace(sigmaM_min, sigmaM_max, n_bins + 1)"""
        
  new_bin_def = """        # Define 6 sSFR bins for this galaxy (spanning 95% population)
        n_bins = 6
        ssfr_min = np.nanpercentile(log_ssfr[valid_all], 2.5)
        ssfr_max = np.nanpercentile(log_ssfr[valid_all], 97.5)
        ssfr_bin_edges = np.linspace(ssfr_min, ssfr_max, n_bins + 1)"""
  src = src.replace(old_bin_def, new_bin_def)

  # Replacement 3: loop variables
  src = src.replace(
      "(sigM >= sigmaM_bin_edges[i])",
      "(log_ssfr >= ssfr_bin_edges[i])"
  ).replace(
      "< sigmaM_bin_edges[i+1])",
      "< ssfr_bin_edges[i+1])"
  ).replace(
      "<= sigmaM_bin_edges[i+1])",
      "<= ssfr_bin_edges[i+1])"
  )

  cell['source'] = [line + '\n' for line in src.split('\n')][:-1]
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb', 'w') as f:
      json.dump(nb, f, indent=1)
  print("Cell 23 updated.")
  ```
  Run: `/opt/miniconda3/envs/ICRAR/bin/python /Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cell23.py`

- [ ] **Step 2: Run verification check on Cell 23**
  Verify cell runs successfully against actual FITS data.
  Run:
  ```bash
  /opt/miniconda3/envs/ICRAR/bin/python -c "
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
  exec(''.join(nb['cells'][2]['source']))
  exec(''.join(nb['cells'][3]['source']))
  exec(''.join(nb['cells'][7]['source']))
  exec(''.join(nb['cells'][22]['source']))
  "
  ```
  Expected: Success output showing processed galaxies and plotted contours.

- [ ] **Step 3: Commit**
  ```bash
  git add 20260616_initial_check_rMZR_binned_by_EWHa.ipynb
  git commit -m "feat: update Cell 23 contours to use sSFR binned offsets"
  ```

---

### Task 4: Update Moran-like 2x2 Maps (Cells 27 & 28)

**Files:**
- Modify: `/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb` (Cells 27 and 28, indices 26 and 27)

- [ ] **Step 1: Write apply script for Cell 27 and 28**
  Write a script `/Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cells27_28.py` to modify cells 27 and 28.
  ```python
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
      
  for idx in [26, 27]:
      cell = nb['cells'][idx]
      src = "".join(cell['source'])
      
      # 1. Define log_ssfr
      src = src.replace(
          "        valid_all = valid_sfr & valid_oh & valid_sigM",
          "        valid_all = valid_sfr & valid_oh & valid_sigM\n        log_ssfr = sigSFR - sigM"
      )
      
      # 2. Binning edges
      old_bin_def = """        # Define 12 Σ* bins
        n_bins = 12
        sigmaM_min = np.nanmin(sigM[valid_all])
        sigmaM_max = np.nanmax(sigM[valid_all])
        sigmaM_bin_edges = np.linspace(sigmaM_min, sigmaM_max, n_bins + 1)"""
        
      new_bin_def = """        # Define 12 sSFR bins
        n_bins = 12
        ssfr_min = np.nanmin(log_ssfr[valid_all])
        ssfr_max = np.nanmax(log_ssfr[valid_all])
        ssfr_bin_edges = np.linspace(ssfr_min, ssfr_max, n_bins + 1)"""
      src = src.replace(old_bin_def, new_bin_def)
      
      # 3. Replacements in loop
      src = src.replace(
          "(sigM >= sigmaM_bin_edges[i])",
          "(log_ssfr >= ssfr_bin_edges[i])"
      ).replace(
          "< sigmaM_bin_edges[i+1])",
          "< ssfr_bin_edges[i+1])"
      ).replace(
          "<= sigmaM_bin_edges[i+1])",
          "<= ssfr_bin_edges[i+1])"
      )
      
      # 4. Map names and plot definitions
      src = src.replace(
          "sigM_binned_map = np.full_like(sigM, np.nan)",
          "ssfr_binned_map = np.full_like(sigM, np.nan)"
      ).replace(
          "sigM_binned_map[bin_mask] = bin_center",
          "ssfr_binned_map[bin_mask] = bin_center"
      ).replace(
          "np.ma.masked_invalid(sigM_binned_map)",
          "np.ma.masked_invalid(ssfr_binned_map)"
      ).replace(
          "norm_discrete = BoundaryNorm(sigmaM_bin_edges, cmap_discrete.N)",
          "norm_discrete = BoundaryNorm(ssfr_bin_edges, cmap_discrete.N)"
      ).replace(
          "bin_centers = [(sigmaM_bin_edges[i] + sigmaM_bin_edges[i+1]) / 2 for i in range(n_bins)]",
          "bin_centers = [(ssfr_bin_edges[i] + ssfr_bin_edges[i+1]) / 2 for i in range(n_bins)]"
      )
      
      # 5. Label replacements
      src = src.replace(
          "r'$\\log\\,\\Sigma_*\\;(M_\\odot\\,\\mathrm{kpc}^{-2})$'",
          "r'$\\log\\,\\mathrm{sSFR}\\;(\\mathrm{yr}^{-1})$'"
      ).replace(
          "binned Σ* map",
          "binned sSFR map"
      )
      
      cell['source'] = [line + '\n' for line in src.split('\n')][:-1]
      
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb', 'w') as f:
      json.dump(nb, f, indent=1)
  print("Cells 27 & 28 updated.")
  ```
  Run: `/opt/miniconda3/envs/ICRAR/bin/python /Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cells27_28.py`

- [ ] **Step 2: Run verification checks on Cells 27 and 28**
  Run:
  ```bash
  /opt/miniconda3/envs/ICRAR/bin/python -c "
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
  exec(''.join(nb['cells'][2]['source']))
  exec(''.join(nb['cells'][3]['source']))
  exec(''.join(nb['cells'][26]['source']))
  exec(''.join(nb['cells'][27]['source']))
  "
  ```
  Expected: Success output (runs for NGC4189 and NGC4254).

- [ ] **Step 3: Commit**
  ```bash
  git add 20260616_initial_check_rMZR_binned_by_EWHa.ipynb
  git commit -m "feat: update Cells 27 and 28 2x2 maps to show binned sSFR map"
  ```

---

### Task 5: Update Combined Figure with NGC4383 (Cell 29)

**Files:**
- Modify: `/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb` (Cell 29, index 28)

- [ ] **Step 1: Write apply script for Cell 29**
  Write a script `/Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cell29.py` to modify cell 29.
  *Note: The script will update `compute_spearman_track` to compute sSFR dynamically per-galaxy or globally, and update `compute_spearman_track_from_offsets` to bin by sSFR.*
  ```python
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
  cell = nb['cells'][28]
  src = "".join(cell['source'])
  
  # 1. Update compute_spearman_track
  old_track_fn = """    def compute_spearman_track(sigM, sigSFR, oh_map, n_bins=6):
        valid_all = np.isfinite(sigM) & np.isfinite(sigSFR) & np.isfinite(oh_map)
        if np.sum(valid_all) == 0:
            return np.array([]), np.array([])

        sigM = sigM[valid_all]
        sigSFR = sigSFR[valid_all]
        oh_map = oh_map[valid_all]

        sigmaM_min = np.nanpercentile(sigM, 2.5)
        sigmaM_max = np.nanpercentile(sigM, 97.5)
        sigmaM_edges = np.linspace(sigmaM_min, sigmaM_max, n_bins + 1)
        sigmaM_centers = 0.5 * (sigmaM_edges[:-1] + sigmaM_edges[1:])

        valid_centers = []
        rho_values = []

        for i in range(n_bins):
            if i == 0:
                bin_mask = (sigM >= sigmaM_edges[i]) & (sigM < sigmaM_edges[i+1])
            elif i == n_bins - 1:
                bin_mask = (sigM >= sigmaM_edges[i]) & (sigM <= sigmaM_edges[i+1])
            else:
                bin_mask = (sigM >= sigmaM_edges[i]) & (sigM < sigmaM_edges[i+1])"""
                
  new_track_fn = """    def compute_spearman_track(sigM, sigSFR, oh_map, n_bins=6):
        valid_all = np.isfinite(sigM) & np.isfinite(sigSFR) & np.isfinite(oh_map)
        if np.sum(valid_all) == 0:
            return np.array([]), np.array([])

        sigM = sigM[valid_all]
        sigSFR = sigSFR[valid_all]
        oh_map = oh_map[valid_all]
        log_ssfr = sigSFR - sigM

        ssfr_min = np.nanpercentile(log_ssfr, 2.5)
        ssfr_max = np.nanpercentile(log_ssfr, 97.5)
        ssfr_edges = np.linspace(ssfr_min, ssfr_max, n_bins + 1)
        ssfr_centers = 0.5 * (ssfr_edges[:-1] + ssfr_edges[1:])

        valid_centers = []
        rho_values = []

        for i in range(n_bins):
            if i == 0:
                bin_mask = (log_ssfr >= ssfr_edges[i]) & (log_ssfr < ssfr_edges[i+1])
            elif i == n_bins - 1:
                bin_mask = (log_ssfr >= ssfr_edges[i]) & (log_ssfr <= ssfr_edges[i+1])
            else:
                bin_mask = (log_ssfr >= ssfr_edges[i]) & (log_ssfr < ssfr_edges[i+1])"""
  src = src.replace(old_track_fn, new_track_fn)

  # 2. Update compute_spearman_track_from_offsets
  old_track_offsets_fn = """    def compute_spearman_track_from_offsets(sigM, sfr_off, oh_off, n_bins=6):
        \"\"\"Like compute_spearman_track but takes pre-subtracted offsets (Option 1).\"\"\"
        valid_all = np.isfinite(sigM) & np.isfinite(sfr_off) & np.isfinite(oh_off)
        if np.sum(valid_all) == 0:
            return np.array([]), np.array([])

        sigM = sigM[valid_all]
        sfr_off = sfr_off[valid_all]
        oh_off = oh_off[valid_all]

        sigmaM_min = np.nanpercentile(sigM, 2.5)
        sigmaM_max = np.nanpercentile(sigM, 97.5)
        sigmaM_edges = np.linspace(sigmaM_min, sigmaM_max, n_bins + 1)
        sigmaM_centers = 0.5 * (sigmaM_edges[:-1] + sigmaM_edges[1:])

        valid_centers = []
        rho_values = []

        for i in range(n_bins):
            if i == 0:
                bin_mask = (sigM >= sigmaM_edges[i]) & (sigM < sigmaM_edges[i+1])
            elif i == n_bins - 1:
                bin_mask = (sigM >= sigmaM_edges[i]) & (sigM <= sigmaM_edges[i+1])
            else:
                bin_mask = (sigM >= sigmaM_edges[i]) & (sigM < sigmaM_edges[i+1])"""
                
  new_track_offsets_fn = """    def compute_spearman_track_from_offsets(sigM, sfr_off, oh_off, n_bins=6):
        \"\"\"Like compute_spearman_track but takes pre-subtracted offsets (Option 1).\"\"\"
        # Note: here sigM argument actually receives log_ssfr
        log_ssfr = sigM
        valid_all = np.isfinite(log_ssfr) & np.isfinite(sfr_off) & np.isfinite(oh_off)
        if np.sum(valid_all) == 0:
            return np.array([]), np.array([])

        log_ssfr = log_ssfr[valid_all]
        sfr_off = sfr_off[valid_all]
        oh_off = oh_off[valid_all]

        ssfr_min = np.nanpercentile(log_ssfr, 2.5)
        ssfr_max = np.nanpercentile(log_ssfr, 97.5)
        ssfr_edges = np.linspace(ssfr_min, ssfr_max, n_bins + 1)
        ssfr_centers = 0.5 * (ssfr_edges[:-1] + ssfr_edges[1:])

        valid_centers = []
        rho_values = []

        for i in range(n_bins):
            if i == 0:
                bin_mask = (log_ssfr >= ssfr_edges[i]) & (log_ssfr < ssfr_edges[i+1])
            elif i == n_bins - 1:
                bin_mask = (log_ssfr >= ssfr_edges[i]) & (log_ssfr <= ssfr_edges[i+1])
            else:
                bin_mask = (log_ssfr >= ssfr_edges[i]) & (log_ssfr < ssfr_edges[i+1])"""
  src = src.replace(old_track_offsets_fn, new_track_offsets_fn)

  # 3. Update global track offset collection
  old_offset_loop = """                            for bi in range(6):
                                liking = 1
                                if bi == 0:
                                    bm = (sM_v >= sM_edges[bi]) & (sM_v < sM_edges[bi+1])
                                elif bi == 5:
                                    bm = (sM_v >= sM_edges[bi]) & (sM_v <= sM_edges[bi+1])
                                else:
                                    bm = (sM_v >= sM_edges[bi]) & (sM_v < sM_edges[bi+1])"""
  # Wait, let's look at the exact code from cell_28.py lines 480-486:
  # for bi in range(6):
  #     if bi == 0:
  #         bm = (sM_v >= sM_edges[bi]) & (sM_v < sM_edges[bi+1])
  #     elif bi == 5:
  #         bm = (sM_v >= sM_edges[bi]) & (sM_v <= sM_edges[bi+1])
  #     else:
  #         bm = (sM_v >= sM_edges[bi]) & (sM_v < sM_edges[bi+1])
  old_offset_loop = """                            for bi in range(6):
                                if bi == 0:
                                    bm = (sM_v >= sM_edges[bi]) & (sM_v < sM_edges[bi+1])
                                elif bi == 5:
                                    bm = (sM_v >= sM_edges[bi]) & (sM_v <= sM_edges[bi+1])
                                else:
                                    bm = (sM_v >= sM_edges[bi]) & (sM_v < sM_edges[bi+1])"""

  new_offset_loop = """                            log_ssfr_v = sSFR_v - sM_v
                            ssfr_min = np.nanpercentile(log_ssfr_v, 2.5)
                            ssfr_max = np.nanpercentile(log_ssfr_v, 97.5)
                            ssfr_edges = np.linspace(ssfr_min, ssfr_max, 7)
                            for bi in range(6):
                                if bi == 0:
                                    bm = (log_ssfr_v >= ssfr_edges[bi]) & (log_ssfr_v < ssfr_edges[bi+1])
                                elif bi == 5:
                                    bm = (log_ssfr_v >= ssfr_edges[bi]) & (log_ssfr_v <= ssfr_edges[bi+1])
                                else:
                                    bm = (log_ssfr_v >= ssfr_edges[bi]) & (log_ssfr_v < ssfr_edges[bi+1])"""
  src = src.replace(old_offset_loop, new_offset_loop)

  # 4. Update the append to global_pixels (Option 1) to use log_ssfr_v instead of sM_v
  src = src.replace(
      "spearman_results[reg][ind]['global_pixels']['sigM'].extend(sM_v[bm][vc])",
      "spearman_results[reg][ind]['global_pixels']['sigM'].extend(log_ssfr_v[bm][vc])"
  )

  # 5. Label and Limits
  src = src.replace(
      "ax.set_xlabel(r'$\\log\\,\\Sigma_* \\; (M_\\odot\\,\\mathrm{kpc}^{-2})$', fontsize=22)",
      "ax.set_xlabel(r'$\\log\\,\\mathrm{sSFR} \\ (\\mathrm{yr}^{-1})$', fontsize=22)"
  ).replace(
      "ax.set_xlim(6, 10)",
      "ax.set_xlim(-11.2, -8.6)"
  ).replace(
      "sigmaM_centers[i]",
      "ssfr_centers[i]"
  )

  cell['source'] = [line + '\n' for line in src.split('\n')][:-1]
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb', 'w') as f:
      json.dump(nb, f, indent=1)
  print("Cell 29 updated.")
  ```
  Run: `/opt/miniconda3/envs/ICRAR/bin/python /Users/Igniz/.gemini/antigravity/brain/8957a08f-74fd-4cf6-9ed8-b8ca49c44796/scratch/apply_cell29.py`

- [ ] **Step 2: Run verification check on Cell 29**
  Verify cell runs successfully against actual FITS data.
  Run:
  ```bash
  /opt/miniconda3/envs/ICRAR/bin/python -c "
  import json
  with open('/Users/Igniz/Desktop/ICRAR/further/20260616_initial_check_rMZR_binned_by_EWHa.ipynb') as f:
      nb = json.load(f)
  exec(''.join(nb['cells'][2]['source']))
  exec(''.join(nb['cells'][3]['source']))
  exec(''.join(nb['cells'][28]['source']))
  "
  ```
  Expected: Success output.

- [ ] **Step 3: Commit**
  ```bash
  git add 20260616_initial_check_rMZR_binned_by_EWHa.ipynb
  git commit -m "feat: update Cell 29 combined figure to show Spearman tracks vs sSFR"
  ```
