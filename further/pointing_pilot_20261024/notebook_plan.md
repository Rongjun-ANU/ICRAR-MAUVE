# All-galaxy notebook plan

1. Complete: validate NGC4321 variance segmentation outside AO exclusion 5770-6050 A.
2. Implement reusable conservative geometry detector (no fixed 2x3 assumption),
   explicit ambiguous status, shared footprint masks and interior extraction masks.
3. Discover products via live CONFIG; compare 7000 and 8900 baseline runs, skip
   single pointings, flag uncertain geometry rather than force ratios.
4. Write self-contained explanatory notebook under /home/rhuang; cache small maps
   and spectra, save diagnostic geometry and normalized adjacent-region ratios.
5. Verify key geometry and ratio invariants, execute NGC4321/NGC4689 and controls,
   then execute notebook population pass with clear per-galaxy status reporting.

AO mask is observed-frame 5770-6050 A from live local specMask_KIN_narrow10_AO.
Avoid this interval in detection, extraction and normalization for AO datasets.
No science products overwritten. All new remote outputs under dedicated home folder.

Implementation complete; seven focused invariants pass locally and remotely.
All four real demonstrations complete, both-run ratio figures generated. Actual
NGC4321 output numerically verified blank across AO interval. Schema/all-cell
compile checks pass. Final remote population execution in progress as job 51487786.
First job stopped after unsupported CONT endpoint was diagnosed; corrected
common-channel preflight and versioned spectrum cache added. A gallery string
escape caught by whole-book compile was corrected before final population job.
