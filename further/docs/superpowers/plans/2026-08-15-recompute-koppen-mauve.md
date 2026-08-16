# Recompute the Köppen criterion for the MAUVE sample

**Goal:** Create and execute a dated notebook that implements the revised controlled-reproduction roadmap for all 40 physical MAUVE galaxies.

**Inputs:** Use `mauve_master_wiki_newclass.fits` for the local MAUVE quantities and stage labels, `MAUVE_galaxy_colors.dat` for stable galaxy identity colours, the 2026-08-15 HyperLeda snapshot embedded in the notebook for `logd25`, `T`, and `vrot`, and the local SIMBAD cross-match for sky positions.

**Method:** Keep the Köppen et al. (2018) H I profile, restoring-force approximation, Virgo ICM beta model, and dark-matter beta model fixed. Replace the historical expected H I calibration with the Jones et al. (2018) AMIGA morphology-dependent relation. Compute deterministic values first, then propagate catalogue errors, AMIGA intrinsic scatter, a stated H I-mass uncertainty floor, a rotation-speed systematic floor, and the factor-of-two local-pressure uncertainty with Monte Carlo draws.

**Interpretation:** Preserve the continuous ratio `p_def / p_loc` and `Delta_RP = log10(p_def / p_loc)`. Treat `p_def > 2 p_loc` as the historical past-stripping threshold, report its probability, and keep non-positive H I deficiencies outside the stripping-radius inversion.

**Validation:** Assert complete one-to-one sample joins, test analytic limiting behaviour, reproduce published Köppen `p_loc` values from a compact paper snapshot, execute every notebook cell, inspect outputs for exceptions/non-finite values, and confirm the dated notebook is the only scientific artifact added.
