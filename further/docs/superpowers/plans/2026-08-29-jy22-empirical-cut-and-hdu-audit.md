# JY22 Empirical Cut And HDU Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Change the JY22 gross-mismatch threshold from 9 to 25 with unambiguous empirical provenance, then audit the 423-extension gas FITS schema and recommend a smaller science/QC contract without deleting current products.

**Architecture:** Keep the threshold in the shared diagnostic constant so spatial and integrated inference remain identical. Propagate empirical wording through tests, the script changelog, runtime messages, README, primary-header provenance, and extension comments. Treat HDU reduction as a separate design review: measure the live schema and storage first, classify extensions by scientific role and reproducibility, and defer deletion until the user approves a target contract.

**Tech Stack:** Python 3.13, NumPy, Astropy FITS, pytest, Markdown.

---

### Task 1: Specify The New Empirical JY22 Cut

**Files:**
- Modify: `/Users/Igniz/Desktop/ICRAR/further/tests/test_model_grid_diagnostics.py`
- Test: `/Users/Igniz/Desktop/ICRAR/further/tests/test_model_grid_diagnostics.py`

- [ ] **Step 1: Change the threshold contract test before production code**

```python
assert JY22_EMPIRICAL_CHI2_MAX == pytest.approx(25.0)
```

- [ ] **Step 2: Add source-contract assertions for empirical wording**

```python
assert "empirical" in source.lower()
assert "chi2<=25" in source
assert "chi2>25" in source
```

- [ ] **Step 3: Run the focused tests and verify RED**

Run: `/opt/miniconda3/envs/ICRAR/bin/python -m pytest tests/test_model_grid_diagnostics.py -q`

Expected: failure because the implementation still exposes 9 and the old nominal wording.

### Task 2: Implement And Document The Empirical Cut

**Files:**
- Modify: `/Users/Igniz/Desktop/ICRAR/further/model_grid_diagnostics.py`
- Modify: `/Users/Igniz/Desktop/ICRAR/further/SFR+Z.py`
- Modify: `/Users/Igniz/Desktop/ICRAR/further/README.md`
- Test: `/Users/Igniz/Desktop/ICRAR/further/tests/test_model_grid_diagnostics.py`

- [ ] **Step 1: Change the shared constant**

```python
JY22_EMPIRICAL_CHI2_MAX = 25.0
JY22_FORMAL_CHI2_MAX = JY22_EMPIRICAL_CHI2_MAX
```

- [ ] **Step 2: Replace statistical/physical adequacy language with empirical gross-mismatch language**

Use `JY22_EMPIRICAL_CHI2_MAX` in production code and retain the former name only as a compatibility alias for existing notebooks. Use `empirical chi2<=25 gross-mismatch cut` in runtime output, FITS primary-header comments, extension comments, README QC documentation, and the script changelog. Preserve `JY22_VALID` as numerical validity and preserve existing flag bits.

- [ ] **Step 3: Run the focused test and verify GREEN**

Run: `/opt/miniconda3/envs/ICRAR/bin/python -m pytest tests/test_model_grid_diagnostics.py -q`

Expected: all tests pass.

- [ ] **Step 4: Run syntax and focused source checks**

Run: `/opt/miniconda3/envs/ICRAR/bin/python -m py_compile SFR+Z.py model_grid_diagnostics.py`

Run: `rg -n "chi2.?[<>]=?9|chi2>9|nominal.*chi2|1 dof" SFR+Z.py model_grid_diagnostics.py README.md tests/test_model_grid_diagnostics.py`

Expected: compilation succeeds and no active JY22 metadata retains the old threshold or interpretation.

### Task 3: Audit The Live FITS Extension Contract

**Files:**
- Inspect: `/Users/Igniz/Desktop/ICRAR/further/SFR+Z.py`
- Inspect: `/Users/Igniz/Desktop/ICRAR/further/model_grid_diagnostics.py`
- Inspect: `/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8/*/*_gas_bin_maps_further.fits`

- [ ] **Step 1: Measure the live products**

Use Astropy in read-only mode to record extension counts, shapes, dtypes, per-extension byte costs, duplicate HII/SF content, and the number of products with identical schemas.

- [ ] **Step 2: Classify every extension family**

Assign each family to core science, essential reproducibility/QC, optional deep diagnostics, or redundant/derivable candidate. Keep raw measurements needed to recompute science outputs and distinguish flags from duplicated convenience masks.

- [ ] **Step 3: Produce a concrete reduction proposal without changing FITS outputs**

Report exact extension families to retain, move to a separate QC FITS product, or consider dropping. Quantify the extension and approximate storage reduction, and identify downstream notebook/script references that must be checked before deletion.

### Task 4: Final Verification And Handoff

**Files:**
- Inspect: all modified files above
- Preserve: `/Users/Igniz/Desktop/ICRAR/further/20260829_check_surviving_SF_distribution.ipynb`

- [ ] **Step 1: Review the scoped diff and worktree**

Run: `git -c core.fsmonitor=false status --short`

Run: `git -c core.fsmonitor=false diff -- further/SFR+Z.py further/model_grid_diagnostics.py further/README.md further/tests/test_model_grid_diagnostics.py further/docs/superpowers/plans/2026-08-29-jy22-empirical-cut-and-hdu-audit.md`

Expected: only the requested threshold/provenance work and this plan are modified; the unrelated untracked notebook remains untouched.

- [ ] **Step 2: Append the required work-log entry**

Record files changed, focused tests, syntax checks, the 423-HDU audit result, unresolved risks, and that no production FITS files were rerun or deleted.

- [ ] **Step 3: Do not commit**

Leave all scoped changes uncommitted unless the user separately requests a commit or PR.
