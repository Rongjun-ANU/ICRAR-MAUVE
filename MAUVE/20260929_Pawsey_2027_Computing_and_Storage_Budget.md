---
title: "MAUVE 2027 Pawsey computing and storage budget"
subtitle: "One 8900-Angstrom reference pass and scalable SU accounting"
date: "29 September 2026"
author: "MAUVE-MUSE / ICRAR research note"
---

# MAUVE 2027 Pawsey computing and storage budget

**29 September 2026 | Setonix nGIST v7.6.8 | One 8900-Angstrom reference pass**

## Executive estimate

For **40 galaxies represented by 39 datacubes**, one complete 4800-8900 Angstrom pass is estimated at **1,143.018 node-hours**. Each job requests **128 physical CPU cores**. Pawsey's CPU accounting is **1 SU per requested physical core-hour**, so the pass costs approximately **1,143.018 hours x 128 cores = 146,306 SU**. One **20% allowance** for tests, retries, and related science is **29,261 SU**; the **buffered single-pass reference is 175,568 SU**.

The number of future science configurations is open. If **N** comparable 39-cube passes are ultimately planned, their modelled CPU budget is **`SU(N) = N x 175,567.53`**, with **`39 x N` principal jobs**. This is a scaling rule, not an assumption that N has already been chosen. Future wavelength and binning variants should be benchmarked before claiming they have exactly this cost.

The user's storage estimates imply **10.414 TB** to retain all input cubes plus **one** output family. With 15% headroom and upward rounding, this is a **15 TB one-pass Acacia planning value**. Active scratch projects to **20 TB** with working headroom if completed products are transferred away between passes. Multi-pass Acacia capacity depends on how many output families must coexist; Section 4 gives the formula.

## 1. Scope and live evidence

The [29 September status report](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/27_status_log_20260929_105522.txt) marks **35 regular 8900-Angstrom cube products finished**, representing 36 galaxies because NGC4567_8 contains two. Of those, **32 cube IDs** have a clean uninterrupted full-run interval in their product `LOGFILE`. NGC4222, NGC4298 and NGC4450 completed after skipped-module resumptions; the last successful segment is not a new complete-run timing. NGC4216, NGC4548, NGC4579 and NGC4689 have no finished regular product log. All 39 cube IDs still enter the budget, with seven missing full-run times imputed and labelled below.

The [29 September finished-job report](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/finished_efficiency_20260929_130155.txt) covers only the shorter 7000-Angstrom jobs. It has no start time, elapsed time, charged SU, or regular 8900-Angstrom rows. It cannot revise the 8900-Angstrom reference hours. The earlier [May](/Users/Igniz/Desktop/ICRAR/further/finished_efficiency_20260528_111829.txt) and [July](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/finished_efficiency_20260728_110558.txt) regular-job snapshots contain peak memory; **29 of the 32** clean full runs can be uniquely matched to them. The earlier figure's 26 memory points came from an overly strict same-minute join: NGC4189, NGC4254 and NGC4606 crossed that minute boundary by 3-6 seconds. The corrected join uses the same cube ID and a unique end-time match within 60 seconds.

The regular product `CONFIG` files use `READ_DATA` 4800-8900 Angstrom, `GENERAL.NCPU=128`, `GAS.LEVEL=BOTH`, and enabled SFH. `LOGFILE` intervals approximate wall time; they are not an exact Slurm charged-SU ledger. The 2025 [Computational Methodology](/Users/Igniz/Desktop/ICRAR/MAUVE/Computational-Methodology.pdf) used CANFAR-era per-realisation assumptions, so this estimate uses the live Setonix full runs instead.

## 2. One complete 8900-Angstrom pass

The 32 measured, uninterrupted full-run intervals sum to **633.4106 node-hours**. Their uncompressed input FITS sizes span **5.94-56.08 GiB**, and their elapsed intervals span **3.22-90.25 hours**. We sum the individual intervals rather than multiplying one summary statistic by 39.

For the seven untimed cubes, we use the plot-based **22 GiB** size boundary and the **maximum observed 8900-Angstrom full-run time in each group**. The **23 measured small cubes (<22 GiB)** have a maximum of **29.170 hours** (NGC4396); the **9 measured large cubes (>=22 GiB)** have a maximum of **90.253 hours** (NGC4654). The five known-size proxies are NGC4222 and NGC4450 (small) and NGC4216, NGC4298 and NGC4689 (large).

NGC4548 and NGC4579 have no row in the input-size inventory. Assigning each the user's revised **50 GiB** input size places both in the large group. Each receives the **90.253-hour observed large-group maximum**, a conservative planning proxy rather than a measurement. A sample maximum is not a guaranteed upper bound for an unmeasured cube, especially for NGC4216, whose 58.55-GiB input is larger than every directly timed cube.

| Contribution to one 39-cube pass | Cube IDs | Node-hours | Basis |
|:--|--:|--:|:--|
| Clean full 8900-Angstrom runs | 32 | 633.4106 | Sum of observed LOGFILE intervals |
| Known-size small proxies | 2 | 58.3406 | 2 x 29.170-hour observed maximum |
| Known-size large proxies | 3 | 270.7600 | 3 x 90.253-hour observed maximum |
| Assumed 50-GiB large proxies | 2 | 180.5067 | 2 x 90.253-hour observed maximum |
| **Complete reference pass** | **39** | **1,143.0178** | Unrounded values summed by the notebook |

### Hours to service units

[Pawsey's 2027 CPU accounting table](https://pawsey.atlassian.net/wiki/spaces/US/pages/51931328/The+Pawsey+Partner+Merit+Allocation+Scheme) assigns **one SU per requested physical CPU core-hour**. The full pipeline jobs request one node with **128 physical cores**. Thus one node-hour uses **128 SU**; multiply the total node-hours by 128. The efficiency report's `Alloc=256` counts logical threads and is not a 256-physical-core charge in this accounting.

```text
Reference node-hours = 633.4106 + 58.3406
                     + 270.7600 + 180.5067
                     approx 1,143.0178 node-hours
Reference SU         = unrounded reference node-hours x 128 physical cores
                     approx 146,306.28 SU
20% allowance        = unrounded reference SU x 0.20
                     approx 29,261.26 SU
Buffered one-pass SU = unrounded reference SU x 1.20
                     approx 175,567.53 SU
```

Displayed terms are rounded; all totals use the unrounded per-cube values in the notebook. The 20% is a planning allowance, not a measured retry rate. The earlier `pawsey1308` record shows **291,884 of 300,000 SU used in one quarter (97.3%)**; four 300,000-SU quarters are 1.2 million SU of annual allocation, but that single record does not establish annual usage.

## 3. Scale when the future run count is known

For **N positive integer** comparable passes, the notebook's `scale_reference_passes(N)` reports principal jobs, reference node-hours, reference SU, allowance SU, and total SU. The algebra is:

```text
Principal jobs(N) = 39 x N
Reference SU(N)   = 146,306.28 x N
Allowance SU(N)   = 29,261.26 x N
Buffered SU(N)    = 175,567.53 x N
```

N is deliberately **not assigned a proposal value here**. The 8900-Angstrom pass is the reference unit; different wavelength limits, separate gas binning, or reuse of existing modules may change future run costs. The [2027 Pawsey Partner scheme](https://pawsey.atlassian.net/wiki/spaces/US/pages/51931328/The+Pawsey+Partner+Merit+Allocation+Scheme) has a **1,000,000-SU full-year application minimum**, which is a policy floor separate from the **175,568-SU** buffered single-pass workload. A final annual application amount requires a chosen run count and justification of any work above the reference estimate.

## 4. Storage for one pass and future retention

The user reports **8 TB of output for one run type across 36 galaxies/35 completed cubes**, **1.5 TB for all input cubes**, and **14 TB now on scratch**. These decimal-TB figures are user-provided, not outputs of the job reports. Scaling by cube count gives one 39-cube output family of **`8 x 39/35 = 8.914 TB`**. Inputs plus one retained family use **`1.5 + 8 x 39/35 = 10.414 TB`**. A 15% planning margin is **11.976 TB**, rounded up in 5-TB steps to **15 TB** for one pass.

For **N output families retained together**, use the notebook's `scale_storage(N)`:

```text
Acacia minimum(N) = 1.5 + N x (8 x 39/35) TB
Acacia with margin(N) = 1.15 x Acacia minimum(N)
Planning request(N) = round upward to the next 5 TB
```

This assumes all N output families coexist on Acacia. The long-term capacity request should be set after deciding how many products will be retained and for how long. [Pawsey's Acacia guidance](https://pawsey.atlassian.net/wiki/spaces/US/pages/51924510/Acacia+-+Troubleshooting) directs requests above 10 TB to a separate managed-storage application; the [2027 scheme](https://pawsey.atlassian.net/wiki/spaces/US/pages/51931328/The+Pawsey+Partner+Merit+Allocation+Scheme) lists only 0.5 TB as the default project allocation.

Scratch is active working space, not another term in Acacia retention. The comparable 39-cube projection is **`14 x 39/35 = 15.6 TB`**; 25% working headroom gives **19.5 TB**, rounded to a conditional **20 TB** target. That target requires completed products to move to Acacia before another full pass occupies scratch. [Pawsey filesystem guidance](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925876/Pawsey+Filesystems+and+their+Use) describes scratch as temporary.

## 5. What the 22 GiB size cut does and does not establish

The **22 GiB split** is a two-group runtime-imputation choice suggested by the current timing/memory plot. It is **not a verified rule for Slurm partition assignment**. The live [Setonix creation script](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/27_creation.sh) lists highmem cube IDs explicitly. For example, NGC4330 is about 21.1 GiB yet is on that highmem list, while NGC4689 is about 22.5 GiB but is not. Queue selection also depends on actual memory and wall time. The size inventory's field named `size_GB` is calculated as bytes divided by `1024**3`, so this report consistently labels it **GiB**. The two assumed input sizes are 50 GiB directly.

The 8900-Angstrom diagnostic figure still has **32** clean timing points and **29** matched memory points because only those have direct evidence; the separate proposal figure and Appendix A contain **all 39** budgeted cube IDs. Three regular jobs finished after the May/July efficiency snapshots and lack a matched regular MaxRSS record. The maximum matched 8900-Angstrom MaxRSS is **799 GiB**, so memory requirements should be checked per cube before assigning a queue.

## 6. Measurements that would sharpen this estimate

1. Obtain the actual input FITS sizes for NGC4548 and NGC4579; their 50-GiB values and large-group runtimes are planning assumptions.
2. Collect `sacct` elapsed times, requested physical cores, and charged-SU records for completed, interrupted, and restarted 8900-Angstrom jobs. This would replace LOGFILE wall-time approximation and quantify retry cost.
3. Benchmark representative future wavelength and separate gas-binning passes before treating the 8900-Angstrom cost as an exact price for every configuration. The current longest clean 8900-Angstrom run is 90.25 hours, close to the 96-hour long/highmem queue limit.
4. Audit output-family sizes, retention periods, and peak scratch occupancy before converting the single-pass Acacia formula into a multi-pass storage request.

## 7. Proposal-ready methodology paragraph

We estimate Setonix CPU time from one complete 4800-8900 Angstrom nGIST v7.6.8 pass over 39 datacubes representing 40 MAUVE galaxies. Thirty-two uninterrupted full runs contribute 633.4 measured node-hours. Five known-size cubes without comparable full runs add 329.1 hours using the largest measured 8900-Angstrom runtime in their matching small or large input-size group, divided at 22 GiB; NGC4548 and NGC4579, each provisionally assigned 50 GiB, add 180.5 hours using the large-group maximum. The resulting pass is 1,143.0 node-hours. At 128 requested physical CPU cores and one SU per core-hour, it costs 146,306 SU; one 20% allowance raises the buffered reference to **175,568 SU per pass**. Once the future number N of comparable passes is chosen, its CPU estimate is N times this buffered reference. One retained output family plus all inputs needs 10.41 TB on Acacia before headroom, and the storage request scales with the number of output families retained together.

## Appendix A. Per-cube 8900-Angstrom reference ledger

This table is generated from the same unrounded notebook calculations. `Full LOGFILE` is a clean complete attempt; the small/large observed maxima are explicit planning proxies. Input size is uncompressed FITS **GiB**; the two assumed sizes are **50 GiB each**.

| Cube ID | Input size (GiB) | Size group | Runtime basis | Hours | SU at 128 cores |
|:--|--:|:--|:--|--:|--:|
| NGC4607 | 5.94 | small (<22 GiB) | Full LOGFILE | 3.54 | 453 |
| NGC4606 | 5.95 | small (<22 GiB) | Full LOGFILE | 3.22 | 413 |
| NGC4189 | 6.01 | small (<22 GiB) | Full LOGFILE | 3.46 | 443 |
| NGC4351 | 6.02 | small (<22 GiB) | Full LOGFILE | 5.03 | 644 |
| IC3392 | 6.03 | small (<22 GiB) | Full LOGFILE | 5.11 | 654 |
| NGC4424 | 6.04 | small (<22 GiB) | Full LOGFILE | 5.40 | 691 |
| NGC4694 | 6.10 | small (<22 GiB) | Full LOGFILE | 5.75 | 736 |
| NGC4388 | 8.86 | small (<22 GiB) | Full LOGFILE | 8.51 | 1,089 |
| NGC4405 | 10.91 | small (<22 GiB) | Full LOGFILE | 6.56 | 840 |
| NGC4402 | 11.55 | small (<22 GiB) | Full LOGFILE | 7.99 | 1,022 |
| NGC4294 | 11.85 | small (<22 GiB) | Full LOGFILE | 12.97 | 1,660 |
| NGC4064 | 11.87 | small (<22 GiB) | Full LOGFILE | 10.15 | 1,300 |
| NGC4522 | 12.26 | small (<22 GiB) | Full LOGFILE | 6.72 | 860 |
| NGC4580 | 12.63 | small (<22 GiB) | Full LOGFILE | 16.07 | 2,057 |
| NGC4394 | 12.73 | small (<22 GiB) | Full LOGFILE | 6.51 | 833 |
| NGC4419 | 12.75 | small (<22 GiB) | Full LOGFILE | 16.30 | 2,086 |
| NGC4383 | 14.07 | small (<22 GiB) | Full LOGFILE | 21.25 | 2,719 |
| NGC4302 | 14.41 | small (<22 GiB) | Full LOGFILE | 13.39 | 1,714 |
| NGC4293 | 15.48 | small (<22 GiB) | Full LOGFILE | 20.03 | 2,564 |
| NGC4450 | 16.45 | small (<22 GiB) | Small-cube maximum | 29.17 | 3,734 |
| NGC4457 | 16.57 | small (<22 GiB) | Full LOGFILE | 12.61 | 1,614 |
| NGC4535 | 17.27 | small (<22 GiB) | Full LOGFILE | 9.80 | 1,255 |
| NGC4222 | 20.99 | small (<22 GiB) | Small-cube maximum | 29.17 | 3,734 |
| NGC4396 | 21.02 | small (<22 GiB) | Full LOGFILE | 29.17 | 3,734 |
| NGC4330 | 21.11 | small (<22 GiB) | Full LOGFILE | 7.53 | 964 |
| NGC4689 | 22.54 | large (>=22 GiB) | Large-cube maximum | 90.25 | 11,552 |
| NGC4698 | 22.79 | large (>=22 GiB) | Full LOGFILE | 10.90 | 1,396 |
| NGC4380 | 31.83 | large (>=22 GiB) | Full LOGFILE | 20.78 | 2,659 |
| NGC4567_8 | 32.42 | large (>=22 GiB) | Full LOGFILE | 72.30 | 9,254 |
| NGC4298 | 32.68 | large (>=22 GiB) | Large-cube maximum | 90.25 | 11,552 |
| NGC4569 | 35.27 | large (>=22 GiB) | Full LOGFILE | 22.46 | 2,875 |
| NGC4254 | 37.96 | large (>=22 GiB) | Full LOGFILE | 34.03 | 4,355 |
| NGC4321 | 38.77 | large (>=22 GiB) | Full LOGFILE | 61.72 | 7,900 |
| NGC4654 | 40.81 | large (>=22 GiB) | Full LOGFILE | 90.25 | 11,552 |
| NGC4501 | 46.95 | large (>=22 GiB) | Full LOGFILE | 45.91 | 5,877 |
| NGC4548 | 50.00 | large (>=22 GiB) | Assumed 50 GiB; large maximum | 90.25 | 11,552 |
| NGC4579 | 50.00 | large (>=22 GiB) | Assumed 50 GiB; large maximum | 90.25 | 11,552 |
| NGC4192 | 56.08 | large (>=22 GiB) | Full LOGFILE | 37.98 | 4,862 |
| NGC4216 | 58.55 | large (>=22 GiB) | Large-cube maximum | 90.25 | 11,552 |

## Source files and policy

- [Executed budget notebook](/Users/Igniz/Desktop/ICRAR/further/20260929_Pawsey_proposal_budget_accounting.ipynb): editable assumptions, calculations, diagnostic and 39-cube plots, and this same report narrative.
- [Regular full-run products](/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8): per-cube `LOGFILE` and `CONFIG` files.
- [Input-size inventory](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/cube_sizes_v3tk.csv), [39-cube allowlist](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/27_galaxies.sh), and [queue-assignment script](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/27_creation.sh).
- [Latest product status](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/27_status_log_20260929_105522.txt), [latest finished-job report](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/finished_efficiency_20260929_130155.txt), and the [May](/Users/Igniz/Desktop/ICRAR/further/finished_efficiency_20260528_111829.txt) and [July](/Users/Igniz/Desktop/ICRAR/nGIST_v7.6.8_setup/config_setup/finished_efficiency_20260728_110558.txt) regular efficiency snapshots.
- [Pawsey 2027 Partner scheme](https://pawsey.atlassian.net/wiki/spaces/US/pages/51931328/The+Pawsey+Partner+Merit+Allocation+Scheme), [Setonix scheduling](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925964/Job+Scheduling), [filesystem use](https://pawsey.atlassian.net/wiki/spaces/US/pages/51925876/Pawsey+Filesystems+and+their+Use), and [Acacia guidance](https://pawsey.atlassian.net/wiki/spaces/US/pages/51924510/Acacia+-+Troubleshooting).
