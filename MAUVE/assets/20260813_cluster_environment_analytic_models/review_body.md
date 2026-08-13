# Analytical Models of Galaxy Evolution in Cluster Environments and Ram Pressure Stripping

**A systematic-style scoping review with a MAUVE implementation roadmap**

**Search and verification date:** 2026-08-13 (Australia/Perth)  
**Scope:** peer-reviewed astronomical models and the validation papers needed to judge them  
**Corpus:** 80 DOI-verified papers plus one foundational DOI-less ICM-profile paper  
**Primary application:** nearby cluster galaxies, especially spatially resolved Virgo systems observed by MAUVE

## Abstract

Analytical models in this field do three different jobs, and the literature is easiest to understand when those jobs are not conflated. First, force and transport models estimate whether gas at a given radius is susceptible to removal. Second, population and phase-space models estimate distributions of infall histories, gas loss, and quenching across satellite populations. Third, tailored forward models combine an orbit, a time-dependent intracluster-medium (ICM) wind, a galaxy potential, and resolved observations to constrain the history of an individual galaxy. No single equation uniquely determines the environmental history of an observed cluster galaxy because its true cluster-centric distance, transverse velocity, wind-disk angle, pre-infall gas distribution, ICM substructure, and prior interactions are normally latent.

The most useful fast model for a cold disk remains the Gunn-Gott restoring-force comparison. Its main modern extension is the finite-pulse treatment of Koppen et al., which shows that long pressure pulses are governed mainly by peak pressure while short pulses are governed by integrated momentum. For a hot circumgalactic halo, the McCarthy et al. spherical criterion plus a finite stripping time is the standard analytical prescription. Semi-analytic galaxy-formation models demonstrate why gradual hot-halo removal, shielding of the cold disk, local ICM estimates, and stripping outside a sharp virial boundary all matter for population predictions. Projected phase space is valuable for statistical infall-time inference, but its several-gigayear scatter and projection degeneracy rule out deterministic stage labels for isolated objects. The most informative individual-galaxy results come from multi-observable forward modelling, exemplified by the Virgo sequence of Vollmer and collaborators and the resolved star-formation-history work on NGC 4330.

For MAUVE, the recommended next step is a hierarchical Bayesian forward model that combines (i) a spatial restoring-pressure map from stellar and gas surface densities, (ii) an orbit-conditioned pressure pulse with explicit ICM and geometry priors, (iii) neutral and molecular gas removal, (iv) a non-equilibrium gas-regulator and chemical-evolution layer, and (v) likelihoods for H I, CO, Halpha, stellar populations, metallicity, and gas/stellar kinematics. Its outputs should be posterior distributions for stripping radius, stripped fraction, peak pressure, impulse, time to or since peak, infall time, and mechanism mixture - not a single hard environmental class.

## Executive answer

### What analytical models exist?

The review identifies eight practically distinct families:

1. **Cold-disk ram-pressure criteria** - compare ICM momentum flux with the vertical gravitational restoring force. These predict a stripping radius or susceptible gas fraction.
2. **Time-dependent pressure-pulse models** - distinguish quasi-static stripping from impulsive stripping and add a clock through the orbit-dependent pressure history.
3. **Hot-halo or circumgalactic-medium stripping models** - predict gradual loss of the replenishing reservoir and therefore starvation.
4. **Transport and tidal models** - estimate viscous/turbulent stripping, thermal evaporation, tidal truncation, and cumulative harassment.
5. **Semi-analytic galaxy-formation models** - embed prescriptions for gas removal, cooling, star formation, feedback, tides, and orbital evolution in cosmological merger trees.
6. **Projected-phase-space and quenching-clock models** - map observed projected radius and line-of-sight velocity to statistical infall histories and quenching times.
7. **Gas-regulator and chemical clocks** - infer suppressed inflow, gas consumption, outflow, and enrichment histories from SFR and metallicity.
8. **Tailored individual-galaxy inversions** - forward-model a particular orbit and pressure history, then compare resolved gas, radio, UV, optical, and spectroscopic observables.

### Can an analytical model determine the environment or history of one galaxy?

**It can constrain, but usually cannot uniquely determine, the history.** The answer depends on what "determine" means:

- **Cluster membership or gravitational binding:** dynamical models, caustics, and a cluster potential can provide a membership probability or bound-region estimate.
- **Current susceptibility to ram pressure:** a restoring-force map can identify which disk radii should be vulnerable for an assumed ICM density and relative velocity.
- **Time since infall:** projected phase-space models provide a broad probability distribution, not an exact clock. Oman et al. found a 68 percent error of about 2.6 Gyr using both observable projected coordinates, improving to about 1.9 Gyr only in favorable regions.
- **Pre-peak, near-peak, or post-peak RPS:** this becomes credible only when projected phase space is combined with resolved gas morphology, wind direction, truncation, extraplanar material, stellar-population ages, and a forward model.
- **Dominant physical mechanism:** no one observable is sufficient. Stellar asymmetry suggests tides; gas-only disturbance suggests hydrodynamics; metallicity and SFH constrain gas-supply evolution; all have degeneracies.

### Bottom line for MAUVE

The scientifically defensible product is a **posterior over histories and mechanisms**. A hard label can be reported as a summary of that posterior, but the model should retain uncertainty and competing explanations. A fast Gunn-Gott calculation is an excellent first layer, not a complete evolutionary model.

![Literature timeline](assets/20260813_cluster_environment_analytic_models/figure_01_literature_timeline.png)

*Figure 1. Growth of the included analytical and validation literature. Counts come from the DOI-verified source registry; the apparent recent dominance of diagnostic work reflects the review scope rather than the total publication rate of cluster-galaxy papers.*

## 1. Review questions and definitions

This review addresses four questions:

1. Which analytical or semi-analytical models describe galaxy evolution in groups and clusters?
2. What are their inputs, equations or prescriptions, outputs, calibration domains, and failure modes?
3. Which models can constrain an individual observed galaxy, and which are only reliable statistically?
4. Which combination is most useful for the resolved MAUVE Virgo sample?

### 1.1 What counts as analytical here?

The literature uses "analytical" in several ways. This review labels each use explicitly:

- **Analytical physics model:** a force, momentum, energy, continuity, or chemical-evolution equation with interpretable parameters.
- **Semi-analytical population model (SAM):** analytic sub-grid prescriptions evolved along dark-matter merger trees. The orbit catalogue is numerical, but gas removal and galaxy evolution are prescribed rather than hydrodynamically resolved.
- **Simulation-calibrated analytic fit:** a compact fitting function derived from hydrodynamic simulations, such as an ICM ram-pressure profile.
- **Parametric or empirical clock:** a delay, exponential decline, or onset-radius model fit to observed populations plus simulated orbital histories.
- **Tailored forward model:** numerical or semi-analytical dynamics for a single galaxy, constrained by multiple observed maps. These are included because they are the main route from general equations to individual histories.

Pure observational catalogues and pure hydrodynamic simulations were excluded unless they directly validated, falsified, calibrated, or operationalized an analytical model. Reviews were retained for context and citation chaining. This is a systematic-style scoping review, not a formal PRISMA review and not a claim that every astronomical paper containing the word "model" was captured.

## 2. Search, selection, and verification protocol

### 2.1 Search surfaces

The search used NASA ADS/SciX and ADS full-text records, arXiv, publisher pages (ApJ/AAS, MNRAS/Oxford, A&A, PASA/Cambridge, Nature, and Springer), DOI resolution, and backward/forward citation chaining from the two modern reviews by [Cortese et al. (2021)](https://doi.org/10.1017/pasa.2021.18) and [Boselli, Fossati & Sun (2022)](https://doi.org/10.1007/s00159-022-00140-3). A previous MAUVE environmental atlas was used only as a seed for Virgo-specific model papers; the model literature was searched and verified again for this report.

Representative concept blocks were:

```text
(analytic* OR semi-analytic OR toy model OR gas regulator OR phase space OR orbit)
AND (cluster galax* OR satellite OR Virgo)
AND (ram pressure OR stripping OR starvation OR strangulation OR quenching
     OR tidal OR harassment OR evaporation OR metallicity)
```

Additional searches targeted "stripping radius", "pressure pulse", "impulse", "hot gaseous halo", "local background environment", "quenching time-scale", "time since infall", and named Virgo galaxies with published dynamical models.

### 2.2 Inclusion logic

A work was included when it supplied at least one of the following:

- an equation or compact prescription for an environmental mechanism;
- a simulation-calibrated analytic environmental profile;
- a semi-analytic environmental implementation with documented assumptions;
- a probabilistic mapping from observables to infall or quenching history;
- a tailored inverse/forward model with individual-galaxy timing constraints;
- a validation result that materially changes how an analytical model should be used.

The final corpus contains **80 DOI-bearing works, all 80 resolved through Crossref on 2026-08-13**, plus [Cavaliere & Fusco-Femiano (1976)](https://ui.adsabs.harvard.edu/abs/1976A%26A....49..137C/abstract), whose foundational beta-model paper has no DOI in the checked bibliographic record. The exact registry, DOI audit, source roles, and reproducible scripts accompany this report.

### 2.3 Evidence rules

- Model-derived times are not treated as direct measurements.
- A fit to a population is not assumed to identify the mechanism in one object.
- A simulation-calibrated relation is not called first-principles analytical.
- Current projected position is not treated as true three-dimensional radius.
- Halpha is not automatically equated with in-situ star formation in shocked, diffuse-ionized, or AGN-contaminated gas.
- The same model fitting several correlated observables is weaker than genuinely independent constraints.
- Numerical values affected by an erratum are flagged. In particular, [Jaffe et al. (2018)](https://doi.org/10.1093/mnras/sty500) has a [published unit-related erratum](https://doi.org/10.1093/mnras/sty2774); the authors state that the main conclusions remain, but several model values and figures were corrected.

## 3. The model stack

An individual-galaxy analysis is not one equation. It is a stack linking observables, the cluster, gas-removal physics, subsequent evolution, and statistical inference.

![Model stack](assets/20260813_cluster_environment_analytic_models/figure_02_model_stack.png)

*Figure 2. Recommended separation of model layers. Each layer has different latent variables and should be validated separately.*

The main identifiability problem is visible in this stack. A measured projected radius and line-of-sight velocity do not specify the local ICM density or full relative velocity. Even if the instantaneous ram pressure were known, the disk response depends on inclination, its vertical potential, multiphase structure, and the duration of the pulse. Finally, gas loss does not map one-to-one onto SFR or metallicity because gas consumption, feedback, fallback, mixing, and halted dilution operate on different timescales.

## 4. Core analytical building blocks

### 4.1 Cluster potential and ICM density

Most practical RPS calculations begin with a spherically averaged ICM density. The classic [Cavaliere & Fusco-Femiano beta model (1976)](https://ui.adsabs.harvard.edu/abs/1976A%26A....49..137C/abstract) is

```text
rho_ICM(r) = rho_0 [1 + (r/r_c)^2]^(-3 beta / 2).
```

The ram pressure is then

```text
P_ram(t) = rho_ICM[r(t)] |v_gal(t) - v_ICM(t)|^2.
```

The cluster gravitational potential is commonly represented with the [Navarro-Frenk-White profile](https://doi.org/10.1086/304888). The analytic caustic construction of [Diaferio & Geller (1997)](https://doi.org/10.1086/304075) identifies the projected infall envelope, while [Diaferio (1999)](https://doi.org/10.1046/j.1365-8711.1999.02864.x) turns its amplitude into a mass-profile estimator without assuming dynamical equilibrium. These are environmental building blocks rather than galaxy-evolution models by themselves, and a caustic envelope does not supply an individual infall time.

The spherical approximation is often the dominant uncertainty. [Tonnesen & Bryan (2008)](https://doi.org/10.1086/592066) found that ICM density and velocity substructure can produce more than an order-of-magnitude spread in ram pressure at a fixed cluster radius. [Tecce et al. (2010)](https://doi.org/10.1111/j.1365-2966.2010.17262.x) found that a smooth analytic ICM profile overestimated ram pressure by more than 50 percent at redshift above 0.5 in their comparison. [Tecce, Cora & Tissera (2011)](https://doi.org/10.1111/j.1365-2966.2011.19267.x) supplied a compact mass-, radius-, and redshift-dependent pressure profile for use in dark-matter-only models, while [Vega-Martinez et al. (2022)](https://doi.org/10.1093/mnras/stab2908) recalibrated that idea and showed large radius/redshift biases in the older fit. Therefore, any such profile is a useful prior mean, not a precise local pressure measurement.

### 4.2 Cold-disk restoring-force criterion

The foundational [Gunn & Gott (1972)](https://doi.org/10.1086/151605) face-on condition is

```text
P_ram >= P_restore(R) approximately 2 pi G Sigma_star(R) Sigma_gas(R).
```

The approximation treats the stellar disk as the main source of vertical restoring gravity. A more general and preferable resolved form is

```text
P_restore(R) = Sigma_gas(R) max_z |dPhi(R,z)/dz|,
```

which can include the gas, bulge, dark halo, and finite disk thickness. [Abadi, Moore & Bower (1999)](https://doi.org/10.1046/j.1365-8711.1999.02715.x) showed that the analytical stripping radius agreed reasonably with hydrodynamic simulations, while also demonstrating that the bulge can protect central gas. [Roediger & Hensler (2005)](https://doi.org/10.1051/0004-6361:20042131) formalized the more general maximum-restoring-acceleration comparison and separated the response into displacement, rapid removal, and longer continuous stripping phases.

For exponential stellar and gas disks,

```text
Sigma_star = Sigma_star,0 exp(-R/R_star)
Sigma_gas  = Sigma_gas,0  exp(-R/R_gas),
```

the thin-disk equality yields a convenient closed-form stripping radius:

```text
R_strip = ln[(2 pi G Sigma_star,0 Sigma_gas,0) / P_ram]
          / (1/R_star + 1/R_gas).
```

For an exponential gas disk with x = R_strip/R_gas, the gas fraction initially outside that radius is

```text
f_out = exp(-x) (1 + x).
```

These equations are transparent and fast. Their output is **susceptibility under an assumed pressure**, not the observed mass already removed. Negative R_strip means the assumed pressure exceeds the central thin-disk restoring estimate; an R_strip beyond the real gas edge means negligible direct disk stripping. Both limits require careful treatment rather than blind extrapolation.

### 4.3 Time dependence: peak pressure versus impulse

An orbit through a cluster produces a pressure pulse, not a permanently uniform wind. [Takeda, Nulsen & Fabian (1984)](https://doi.org/10.1093/mnras/208.2.261) analyzed stripping in a changing environment and distinguished response regimes. [Roediger & Brueggen (2007)](https://doi.org/10.1111/j.1365-2966.2007.12241.x) found that the instantaneous criterion tracked the remaining disk radius fairly well along realistic orbits but predicted mass loss too rapidly.

[Koppen et al. (2018)](https://doi.org/10.1093/mnras/sty1610) made the distinction operational:

- For a **long pulse**, the gas has time to respond and the peak pressure largely controls stripping. The Gunn-Gott limit is appropriate.
- For a **short pulse**, momentum transfer controls the outcome. Define the impulse per area

```text
J = integral P_ram(t) dt.
```

Gas is unbound when J/Sigma_gas exceeds the vertical velocity increment required to escape the galaxy potential. Thus two orbits with the same peak pressure can remove different gas masses if their pulse durations differ.

Koppen et al. showed that the pair `(P_max, integral P dt)` compactly describes stripped fractions across pulse durations and applied the framework to 232 Virgo galaxies. This is the most directly useful analytical upgrade to Gunn-Gott for a nearby-cluster project.

![Ram-pressure regimes](assets/20260813_cluster_environment_analytic_models/figure_03_rps_regimes.png)

*Figure 3. Schematic long- and short-pulse limits. The diagonal is not a universal numerical boundary: the galaxy potential, gas column, radius, and pulse shape set the physical normalization.*

### 4.4 Inclination, annuli, and non-radial orbits

[Roediger & Brueggen (2006)](https://doi.org/10.1111/j.1365-2966.2006.10335.x) found that inclination has a modest effect on total mass loss except close to edge-on configurations, although morphology and tails remain geometry-sensitive. [Hester (2006)](https://doi.org/10.1086/505614) developed an analytical treatment for disk and halo gas across groups and clusters and emphasized host-to-satellite mass ratio as a major control on susceptibility. [Singh, Gulati & Bagla (2019)](https://doi.org/10.1093/mnras/stz2523) evolved gas annuli under the net acceleration

```text
a_z(R,t) = [P_ram(t) - P_restore(R)] / Sigma_gas(R),
```

for idealized disks and halo environments. [Singh et al. (2024)](https://doi.org/10.1093/mnras/stae730) extended this family with non-radial EAGLE orbits and an inclination-dependent perpendicular pressure component. Their results reinforce an important inference rule: the face-on, radial-infall calculation is a high-efficiency limiting case, not a neutral default.

### 4.5 Hot-halo stripping and starvation

Cold-disk stripping and removal of the hot circumgalactic reservoir are different problems. The foundational supply-cutoff model of [Larson, Tinsley & Caldwell (1980)](https://doi.org/10.1086/157917) showed how loss of external gas supply followed by consumption of the remaining disk gas can fade a spiral toward an S0. [Balogh, Navarro & Morris (2000)](https://doi.org/10.1086/309323) attached a simple post-accretion gas-consumption history to simulated cluster assembly and recovered smooth star-formation gradients plus an important role for group preprocessing. These models describe the downstream evolution after supply loss; they do not derive the hydrodynamic removal itself.

[McCarthy et al. (2008)](https://doi.org/10.1111/j.1365-2966.2007.12577.x) derived a spherical analogue of Gunn-Gott:

```text
P_ram > alpha G M_sat(<r) rho_hot(r) / r,
```

with a geometry factor alpha of order unity (their calibration favored about 2). Crucially, they did not assume instantaneous removal. Gas outside the limiting radius is removed over

```text
t_ram = beta t_sound,
```

with beta also of order unity. Their calibrated model matched their controlled simulations to about 10 percent for most tested orbits and structures, while a typical satellite could retain about 30 percent of its hot halo after 10 Gyr.

[Font et al. (2008)](https://doi.org/10.1111/j.1365-2966.2008.13698.x) inserted this gradual prescription into GALFORM. Satellites retained hot gas for several Gyr, continuing to replenish their cold gas and correcting the excessive population of faint red satellites produced by instantaneous strangulation. [Zinger et al. (2018)](https://doi.org/10.1093/mnras/stx3329) used a simple analytic halo model plus simulated cluster profiles to argue that much of the hot reservoir can be stripped around or beyond the virial radius, setting up 2-3 Gyr starvation before strong inner-cluster disk stripping.

The recent validation by [Ghosh, Dutta & Sharma (2024)](https://doi.org/10.1093/mnras/stae1345) adds nuance: their simulations support a restoring-pressure description for the ISM, but CGM evolution depends strongly on hydrodynamic drag and momentum transfer rather than gravity alone. Therefore, a one-radius hot-halo cut is useful in a SAM but incomplete for an individual CGM history.

### 4.6 Viscous stripping, turbulence, evaporation, and conduction

The classic RPS balance is not the only gas-removal channel. [Nulsen (1982)](https://doi.org/10.1093/mnras/198.4.1007) showed that viscosity, thermal conduction, and turbulence can strip gas at rates comparable to or exceeding direct ram pressure in some regimes. A common order-of-magnitude turbulent/viscous mass-loss scaling is

```text
dot M_transport approximately pi R_gas^2 rho_ICM v,
```

multiplied by an efficiency controlled by Reynolds number, magnetic suppression, compressibility, and multiphase structure.

[Cowie & McKee (1977)](https://doi.org/10.1086/154911) derived classical and saturated evaporation rates for clouds in a hot medium; [Cowie & Songaila (1977)](https://doi.org/10.1038/266501a0) applied thermal evaporation to gas in galaxies. These models are analytically elegant but difficult to use as precise clocks because magnetic topology and conduction suppression are poorly known. Their main modern value is to prevent a false inference: gas that is not removed by instantaneous Gunn-Gott stripping is not necessarily safe over an orbital time.

[Fillingham et al. (2016)](https://doi.org/10.1093/mnras/stw2131) combined observed H I surface-density profiles with analytic ram-pressure and turbulent-viscous stripping calculations for low-mass satellites. Stripping becomes much more effective below stellar masses of roughly 10^8-10^9 Msun, but a smooth canonical halo removed only about 40-60 percent of the cold gas in their calculation; clumpy high-density gas was needed to reproduce the strongest Local Group quenching. The exact mass scale should not be transplanted unmodified to Virgo, but the analysis shows why diffuse dwarfs and ICM/CGM clumping require special priors.

### 4.7 Tides and harassment

An approximate tidal radius for a satellite of enclosed mass m at distance R from a host of enclosed mass M is

```text
r_t approximately R [m(<r_t) / (3 M(<R))]^(1/3),
```

with the numerical factor modified by the orbit and host mass slope. [Merritt (1983)](https://doi.org/10.1086/160571) analyzed relaxation and tidal stripping in rich clusters. [Moore et al. (1996)](https://doi.org/10.1038/379613a0) introduced the galaxy-harassment picture of repeated high-speed encounters plus the cluster tide. [Gnedin (2003)](https://doi.org/10.1086/344636) found that peaks in tidal forcing need not coincide with closest cluster approach because local substructure and encounters matter; the integrated heating could make an S0 but was insufficient in those models to turn a spiral into an elliptical.

For individual diagnosis, tides are best constrained by the collisionless component: disturbed old stellar isophotes, bars, bridges, shells, double nuclei, or coherent stellar-kinematic asymmetry. Gas-only disturbance with an otherwise regular old stellar disk favors hydrodynamics, but it does not prove RPS because minor interactions can also redistribute gas.

## 5. Detailed review of the major model families

### 5.1 Gunn-Gott and its validation lineage

**Model type:** physically motivated analytic static force-balance approximation.

**Inputs:** ICM density, relative velocity, gas surface density, stellar surface density or full vertical potential, and wind-disk angle.

**Outputs:** susceptible radius, remaining gas radius, and - with an assumed gas profile - an initially susceptible mass fraction.

**Best-supported use:** estimating the current or peak stripping radius of a resolved disk.

**Validation:** [Abadi et al. (1999)](https://doi.org/10.1046/j.1365-8711.1999.02715.x), [Mori & Burkert (2000)](https://doi.org/10.1086/309140), [Roediger & Hensler (2005)](https://doi.org/10.1051/0004-6361:20042131), [Roediger & Brueggen (2007)](https://doi.org/10.1111/j.1365-2966.2007.12241.x), and [Kulier et al. (2023)](https://doi.org/10.3847/1538-4357/aceda3) broadly agree that the criterion captures susceptibility or a disk boundary better than it predicts a detailed mass-loss history. Mori & Burkert's dwarf calculation highlights the importance of the galaxy's central gas density and potential; Kulier et al. found that analytical models generally underpredicted neutral-gas stripping and overpredicted ionized-halo stripping in EAGLE, with feedback enhancing stripping and gas compaction opposing it.

**Main limitations:** unknown 3D velocity and ICM density; disk thickness and non-axisymmetry; multiphase gas; pressure duration; feedback; fallback; inclination; nonthermal ICM; and a distinction between gas being accelerated and actually becoming unbound.

**Individual-galaxy readiness:** high for a spatial susceptibility map, low for an unqualified evolutionary stage.

### 5.2 Finite-pulse and orbit-aware analytical models

**Model type:** time-dependent mechanics of gas parcels or annuli.

**Inputs:** a pressure history P(t), galaxy potential, gas column profile, initial vertical state, and geometry.

**Outputs:** stripped fraction, displacement, velocity, sensitivity to peak pressure and impulse, and sometimes a predicted fallback component.

The trajectory from [Takeda et al. (1984)](https://doi.org/10.1093/mnras/208.2.261) to [Koppen et al. (2018)](https://doi.org/10.1093/mnras/sty1610) is a conceptual advance because it explains why instantaneous pressure alone is not a universal predictor. Koppen et al. used a Lorentz-like pulse and gas parcels in a fixed galaxy potential to bridge long and short durations. [Vollmer et al. (2001)](https://doi.org/10.1086/323368) supplied an influential Virgo implementation with time-varying ram pressure, galaxy orbits, inclination, and fallback in a sticky-particle ISM.

**Strength:** these models turn an environmental force into a history-sensitive prediction while retaining interpretable parameters.

**Limitation:** P(t) is itself inferred from a cluster model and a latent orbit. A fitted pulse is not a direct measurement of the ICM history.

**Individual-galaxy readiness:** medium in isolation; high when multiwavelength maps constrain wind direction and time.

### 5.3 Analytic ICM-pressure profiles

The beta model is simple but its parameters and spherical symmetry may not represent a dynamically young cluster. [Vega-Martinez et al. (2022)](https://doi.org/10.1093/mnras/stab2908) calibrated a universal pressure fit from hydrodynamic cluster resimulations:

```text
P_ram(M,z) = P0(z) [(r/R200) / xi(z)]^[-(3/2) alpha(M200,z)],

alpha(M200,z) = alpha_M(z) log(M200 h^-1 Msun) - 5.5.
```

Their fit improves the mass, radius, and redshift dependence relative to an earlier beta-profile fit, especially above redshift 1. It is useful for a SAM or for constructing priors when direct X-ray constraints are unavailable.

**Important caveat:** the fit predicts a typical radial pressure, not the instantaneous pressure on a specific galaxy in a clumpy, moving ICM. For Virgo, X-ray substructure, the M49/M86 regions, shocks, and bulk flows should be represented by mixture or scatter terms rather than absorbed into a single beta normalization.

### 5.4 Semi-analytic population models

SAMs answer a population question: given merger trees and prescriptions for cooling, star formation, feedback, stripping, and tides, what distribution of satellites is produced? They are not inverse models for an arbitrary individual galaxy.

| Model | Environmental prescription | Main result for interpretation | Principal limitation |
|---|---|---|---|
| [Font et al. 2008](https://doi.org/10.1111/j.1365-2966.2008.13698.x) | McCarthy-style gradual hot-halo stripping | Retaining hot gas for several Gyr prevents satellites becoming too red too quickly | Sensitive to treatment of supernova-ejected gas |
| [Tecce et al. 2010](https://doi.org/10.1111/j.1365-2966.2010.17262.x) | Local ICM density and velocity from gas particles; Gunn-Gott cold disk | Smooth analytic ICM profiles can overestimate pressure; gas depletion depends on host mass and epoch | Non-radiative cluster gas and instantaneous snapshot stripping |
| [Luo et al. 2016](https://doi.org/10.1093/mnras/stw268) | Resolution-aware hot gas plus optional annular cold-gas RPS | Improves convergence but still disagrees with some observed quenched and gas fractions | Simplified ICM and cold-disk response |
| [Cora et al. 2018](https://doi.org/10.1093/mnras/sty1131) | Gradual hot-gas RPS/starvation, cold-disk RPS after shielding, tidal stripping, integrated orphan orbits | Provides a coherent environmental-plus-feedback population model | Calibrated prescriptions can compensate for one another |
| [Ayromlou et al. 2019](https://doi.org/10.1093/mnras/stz1549) | Local Background Environment estimator; hot-gas RPS for centrals and satellites | No physical discontinuity at R200; environment acts beyond the nominal halo edge | Requires simulation particle environment; not directly observable |
| [Ayromlou et al. 2021](https://doi.org/10.1093/mnras/stab1245) | L-GALAXIES with continuous Local Background Environment stripping and global recalibration | Hot-gas removal can begin while objects are still centrals and several R200 from clusters | Improvement combines a new environmental prescription with MCMC recalibration |
| [Xie et al. 2020](https://doi.org/10.1093/mnras/staa2370) | Gradual hot stripping; annular cold RPS after hot-halo shielding; finite 400 Myr removal | Improves quenched, H I, and H2 trends but retains low-mass discrepancies | Assumes H I stripped before H2 and no gas metallicity gradient |
| [Stevens & Brown 2017](https://doi.org/10.1093/mnras/stx1596) | Dark Sage annular disk evolution with environmental and secular processes | Reproduces the relative halo-mass dependence of satellite H I fractions and quiescence | Global calibration and compensating secular/feedback processes limit mechanism uniqueness |
| [Wang et al. 2024](https://doi.org/10.1093/mnras/stae162) | Comparison of GAEA/L-Galaxies with EAGLE/TNG and observations | Environment changes H I earlier and more strongly than sSFR; cold-RPS models perform better in clusters | Model agreement varies with halo mass and observable |

The consistent lessons are:

1. Instantaneous removal of all hot gas is too aggressive.
2. Cold gas should not always be stripped at the moment a galaxy becomes a satellite.
3. The hot halo can shield the cold disk and determines how long cooling continues.
4. Environmental forcing has no sharp onset at R200.
5. H I is a more sensitive early diagnostic than integrated SFR.
6. Population agreement does not validate a unique decomposition into RPS, starvation, feedback, and tides because parameters can compensate.

[Brown et al. (2017)](https://doi.org/10.1093/mnras/stw2991) provides an important observational constraint: starvation by itself does not reproduce the full environmental H I depletion from pairs to clusters. [Wang et al. (2024)](https://doi.org/10.1093/mnras/stae162) similarly shows that models with explicit cold-gas stripping better reproduce the strong central decline of H I in massive haloes.

### 5.5 Projected phase-space models

The observed coordinates are usually

```text
x = R_proj/R200,
y = |Delta v_los|/sigma_cluster.
```

These are a lossy projection of six-dimensional phase space. A caustic or velocity cut is useful for membership, but contamination remains important: [Mamon, Biviano & Murante (2010)](https://doi.org/10.1051/0004-6361/200913948) found that interlopers retain a substantial presence in projected virial cones even after local velocity clipping. [Mahajan, Mamon & Raychaudhury (2011)](https://doi.org/10.1111/j.1365-2966.2011.19236.x) quantified virialized, infall, and backsplash mixtures and showed why low line-of-sight velocity outside the virial radius is backsplash-enriched but not uniquely backsplash.

[Oman, Hudson & Behroozi (2013)](https://doi.org/10.1093/mnras/stt328) used about 570,000 simulated satellite orbits to build conditional distributions of time since infall. Infalling, backsplash, and virialized populations occupy different but overlapping regions. There is no clean universal cut in projected phase space. [Hernandez-Fernandez et al. (2014)](https://doi.org/10.1093/mnras/stt2354) demonstrated how observed star-formation activity can be compared with an analytic RPS-intensity map in cluster projected phase space, but stacking clusters and projecting the orbit still make the resulting zones statistical.

[Rhee et al. (2017)](https://doi.org/10.3847/1538-4357/aa6d6c) associated projected regions with both time since infall and tidal mass loss and found weak dependence on galaxy or host mass in their simulations. [Pasquali et al. (2019)](https://doi.org/10.1093/mnras/sty3530) used constant-mean-infall-time zones to connect SDSS satellite properties with orbital history. These are valuable population tools, but a zone mean is not the time since infall of every galaxy in that zone.

[Jaffe et al. (2015)](https://doi.org/10.1093/mnras/stv100) went one step further by embedding a simple Gunn-Gott RPS prescription into cosmological orbits and predicting a stripping region in phase space. The distribution of H I detections in the z=0.2 BUDHIES cluster agreed well with the model, supporting strong first-infall gas removal. [Jaffe et al. (2018)](https://doi.org/10.1093/mnras/sty500) applied analogous analytic anchoring-force and ICM models to GASP jellyfish galaxies, which preferentially occupy high-velocity, small-radius regions consistent with strong first-infall RPS. The numerical corrections in the [published erratum](https://doi.org/10.1093/mnras/sty2774) should be used for quantitative reuse.

[Yoon et al. (2017)](https://doi.org/10.3847/1538-4357/aa6579) is particularly relevant to MAUVE because it combines VIVA H I morphology with Virgo projected phase space. Early stripping is preferentially found in first-infall regions, active stripping near the core on infall or exit, and severely truncated symmetric disks either deep in the cluster or in a low-velocity outskirts backsplash region. The classification is an evidence synthesis, not an analytical orbit solution.

A newer continuous estimator by [Masson & Parker (2026)](https://doi.org/10.3847/1538-4357/ae4e21) uses symbolic regression on IllustrisTNG projected radius, line-of-sight velocity, stellar mass, and redshift. Its redshift-zero scatter is still about 2.5 Gyr, and its time origin is first crossing of 3 R200. This is not interchangeable with studies that define infall at R200, 2.5 Rvir, satellite transition, or first pericentre. The lesson is not that continuous estimators fail, but that their posterior scatter and clock definition must travel with every quoted age.

### 5.6 Parametric quenching clocks

[Wetzel et al. (2013)](https://doi.org/10.1093/mnras/stt469) popularized the delayed-then-rapid empirical satellite SFH: after infall, a galaxy remains near a central-like SFR for a delay, then quenches rapidly. It fits population statistics but does not uniquely identify whether the trigger is starvation, RPS, feedback, or their sequence.

[Oman & Hudson (2016)](https://doi.org/10.1093/mnras/stw2195) combined projected phase space with simulated orbits and a post-infall delay model to infer cluster quenching times. [Oman et al. (2021)](https://doi.org/10.1093/mnras/staa3845) then convolved each galaxy's orbit-time probability distribution with a simple processed-fraction model and found that H I loss and star-formation quenching are offset clocks; their mock tests also showed that systematic timing shifts can reach several Gyr. [Baxter et al. (2022)](https://doi.org/10.1093/mnras/stac2149) matched GOGREEN clusters to IllustrisTNG accretion histories at redshift about 1 and found mass-dependent quenching times of roughly 1.6 Gyr at 10^10 Msun and 0.6-1 Gyr at 10^11 Msun, consistent with cold-gas depletion and a starvation interpretation for that sample. [Reeves, Hudson & Oman (2023)](https://doi.org/10.1093/mnras/stad1069) jointly forward-modelled stellar ages and quiescent fractions and placed complete quenching several Gyr after first pericentre, while showing sensitivity to the assumed stochastic star-formation histories. [Baxter et al. (2023)](https://doi.org/10.1093/mnras/stad2995) fit both an onset radius and a quenching timescale with MCMC and found distinct acceptable regions corresponding to a starvation-like path and rapid core quenching; their comparison with transition galaxies subsequently favored the starvation branch. The competing acceptable regions still demonstrate that quenched fractions alone admit different physical histories.

These timing results cannot be compared safely without harmonizing zero-points. Published clocks begin at first crossing of R200, 2.5 Rvir, or 3 R200; at satellite transition; or at first pericentre. Differences of several Gyr may therefore be definitional rather than physical. MAUVE should store each timing posterior in both the source definition and a common orbit-model definition.

[Broderick, Roberts & Hudson (2025)](https://doi.org/10.3847/1538-4357/add739) provides a recent resolved application. A slow-then-rapid outside-in RPS toy model reproduced the difference between field and Coma sSFR profiles for 45 MaNGA galaxies. The model supports RPS as an important population mechanism, but the fit remains parametric and does not deliver unique orbits for the individual galaxies.

The possibility of a temporary SFR enhancement is model-dependent rather than a contradiction of outside-in quenching. [Fujita & Nagashima (1999)](https://doi.org/10.1086/307139) coupled a simple molecular-cloud model to radial cluster infall and predicted at most an approximately factor-two compression-driven SFR rise before rapid disk fading in a dense cluster. This idealized result provides a useful burst prior, not a generic expectation: the observable response depends on gas phase, disk structure, pressure rise time, and the selection window.

Preprocessing also has an explicit analytical lineage. [Fujita (2004)](https://doi.org/10.1093/pasj/56.1.29) placed RPS, strangulation, and evaporation in a hierarchical cluster/subcluster model and showed that substantial environmental evolution can precede arrival in the main cluster. For Virgo, this argues for subcluster-association and prior-host terms in the orbit prior rather than treating every galaxy as a pristine first infaller into a single M87 halo.

### 5.7 Gas regulators and chemical clocks

The gas-regulator framework of [Lilly et al. (2013)](https://doi.org/10.1088/0004-637X/772/2/119) and the time-dependent analytic solutions of [Peng & Maiolino (2014)](https://doi.org/10.1093/mnras/stu1288) describe a galaxy as a gas reservoir regulated by inflow, star formation, recycling, and outflow. A cluster-ready local surface-density version is

```text
d Sigma_g/dt = Sigma_in
               - (1 - R + lambda) Sigma_SFR
               - Sigma_strip,

Sigma_SFR = epsilon Sigma_g,

d(Z_g Sigma_g)/dt = Z_in Sigma_in
                     + y (1 - R) Sigma_SFR
                     - Z_g (1 - R + lambda) Sigma_SFR
                     - Z_strip Sigma_strip.
```

This last equation is an explicit synthesis for MAUVE, not a verbatim equation from one source. If stripped gas has the local ISM metallicity, Z_strip = Z_g and stripping does not instantly enrich the gas that remains at that radius. Metallicity changes later because inflow/dilution stops, star formation continues, gas mixes, and low-metallicity outer gas may be preferentially removed. An integrated metallicity can rise simply because the observed gas weighting changes, even without additional local enrichment.

[Carton et al. (2015)](https://doi.org/10.1093/mnras/stv967) showed that a simple local equilibrium model could approximate metallicity-profile shapes through the local gas-to-stellar ratio, but their H I-rich/normal disk calibration is not a cluster stripping model. Applying it to a disturbed Virgo disk requires the non-equilibrium stripping term above.

[Peng, Maiolino & Cochrane (2015)](https://doi.org/10.1038/nature14439) used the stellar-metallicity difference between star-forming and passive galaxies with a closed-box model to infer strangulation over about 4 Gyr. [Trussler et al. (2020)](https://doi.org/10.1093/mnras/stz3286) expanded this to closed- and leaky-box models: their pure-starvation solutions required an extended enrichment phase plus a final heating/ejection step, while mixed starvation and outflow solutions better matched both SFR and metallicity. [Trussler et al. (2021)](https://doi.org/10.1093/mnras/staa3545) then showed that much of the apparent environment dependence of stellar populations weakens after star-forming, green-valley, and passive populations are separated. This is an important selection warning for cluster chemical clocks.

### 5.8 Tailored Virgo inversions

Tailored modelling is the only reviewed family that can jointly constrain a particular galaxy's wind geometry, peak pressure, and time to or since peak. It is also the most assumption-heavy.

[Vollmer et al. (2001)](https://doi.org/10.1086/323368) linked Virgo orbits to time-dependent ram pressure and introduced a modelling lineage based on a sticky-particle multiphase ISM in an external wind. [Vollmer (2009)](https://doi.org/10.1051/0004-6361/200911892) synthesized five Virgo spirals into a model-based sequence: NGC 4501 and NGC 4330 before peak pressure, NGC 4522 shortly after, NGC 4388 later after peak, and NGC 4569 still later. These are model-family estimates rather than uniform posterior measurements, and the pressure normalizations implied by a smooth Virgo ICM can differ by about a factor of two from those required by the galaxy models.

| Galaxy/model | Principal constraints | Representative timing | Important caveat |
|---|---|---|---|
| [NGC 4501: Vollmer et al. 2008](https://doi.org/10.1051/0004-6361:20078139) | Compressed leading edge, diffuse opposite-side gas, multiwavelength asymmetry | about 190-320 Myr before peak; representative about 250 Myr | Wind angle and time are covariant |
| [NGC 4330: Vollmer et al. 2012](https://doi.org/10.1051/0004-6361/201117680) | H I, UV, Halpha, leading-edge truncation, tail | about 100 Myr before peak | SF prescription and multiphase response are simplified |
| [NGC 4330: Fossati et al. 2018](https://doi.org/10.1051/0004-6361/201732373) | Resolved UV-optical SFH plus VESTIGE Halpha | outer suppression began about 400-600 Myr ago; sharp truncation about 100 Myr ago | Different disk zones and SFH forms give different clocks |
| [NGC 4330: Vollmer et al. 2021](https://doi.org/10.1051/0004-6361/202038507) | H I, UV, Halpha, SEDs, spectra, and a new diffuse-gas stripping prescription | jointly tests the radial 500-to-100 Myr quenching gradient rather than supplying one global clock | Dense sticky-particle and diffuse warm/hot gas recipes remain model choices |
| [NGC 4522: Vollmer et al. 2006](https://doi.org/10.1051/0004-6361:20064954) | Gas dynamics and polarized-radio ridge | about 50 Myr after peak | Requires local pressure enhancement and geometry choice |
| [NGC 4388: Vollmer et al. 2018](https://doi.org/10.1051/0004-6361/201731910) | H I, CO, radio, UV, disk and far tail | about 240 Myr after peak in the favored model | No one realization fits every disk and tail observable |
| [NGC 4569: Vollmer et al. 2004](https://doi.org/10.1051/0004-6361:20034552) | Strongly truncated gas disk and dynamical response | about 300 Myr after peak | Cluster profile and adopted orbit control peak pressure |

The multizone chemo-spectrophotometric models of [Boselli et al. (2006)](https://doi.org/10.1086/507766) discriminate starvation from RPS using the radial SFH of NGC 4569: outside-in truncation favors rapid gas removal over a spatially uniform slow decline. [Boselli et al. (2008)](https://doi.org/10.1086/525513) extended the framework to Virgo dwarfs and predicted that repeated/strong RPS can transform star-forming dwarfs rapidly, whereas starvation evolves them more slowly.

These tailored studies establish the right philosophy for MAUVE: fit morphology, kinematics, gas phases, and stellar clocks simultaneously, and retain tensions rather than selecting a visually attractive snapshot as a unique solution.

### 5.9 Morphological candidate diagnostics

[Consolandi et al. (2024)](https://doi.org/10.1093/mnras/stad3881) used asymmetry, concentration, Sersic index, Gini-M20-derived measures, and diagnostic boundaries to find RPS candidates in optical imaging. The method can efficiently triage a large survey, but tidal interactions contaminate the same morphological transition zone. It identifies candidates, not a physical RPS probability unless the training labels, selection function, and contamination model are propagated.

## 6. What the models agree on

1. **RPS is fundamentally radial within a galaxy.** Low restoring pressure in the outer disk makes outside-in gas removal the generic expectation, even though the subsequent SFR response can be complex.
2. **Pressure history matters.** Peak pressure controls slow pulses; integrated impulse matters for short pulses; removal takes finite time.
3. **The hot halo and cold disk respond differently.** Hot-halo loss can start early and drive starvation; dense cold or molecular gas is more resilient.
4. **H I changes before integrated SFR.** Population comparisons consistently find neutral gas to be an early and sensitive environmental tracer.
5. **There is no physical switch at R200.** Accretion shocks, local gas, preprocessing, and backsplash extend environmental effects beyond the nominal virial radius.
6. **Projection makes individual orbital histories broad.** Phase space is probabilistic.
7. **ICM structure is not a small correction.** Local density and velocity fields can move a galaxy far from a spherically averaged prediction.
8. **RPS alone need not complete quenching.** Cold stripping, halted replenishment, consumption, feedback, and sometimes tides form a sequence rather than mutually exclusive explanations. This conclusion is emphasized by the reviews of [Cortese et al. (2021)](https://doi.org/10.1017/pasa.2021.18) and [Boselli et al. (2022)](https://doi.org/10.1007/s00159-022-00140-3).

## 7. Where the models disagree or fail

### 7.1 Instantaneous versus gradual removal

The simplest SAMs removed a satellite's hot gas immediately. Font and McCarthy showed this quenches satellites too fast and misses hot-gas survival. Conversely, slow hot-halo stripping without direct cold-gas removal can leave too much H I in inner-cluster satellites. The defensible prescription is phase-dependent and finite-time, not universally instantaneous or universally slow.

### 7.2 Smooth ICM versus local ICM

Analytic beta/NFW profiles enable transparent calculations. Tecce, Tonnesen & Bryan, Ayromlou, and Vega-Martinez show why local deviations matter. A practical model should treat the smooth profile as a mean plus intrinsic scatter and, where data allow, explicit subcluster/shock components.

### 7.3 Disk radius versus mass-loss history

Gunn-Gott often predicts a boundary well, but gas outside the boundary is not removed at infinite speed. Roediger & Brueggen and Kulier show that the time evolution and phase dependence can differ substantially even when the radius looks correct.

### 7.4 Inclination

Moderate inclination often has less effect on total stripped mass than naive projected-area arguments imply, while near-edge-on flows, tail morphology, and disk compression remain geometry-sensitive. A scalar cosine correction is inadequate for detailed morphology.

### 7.5 RPS versus starvation

These are not exclusive. Removing the hot halo produces starvation; later cold-disk stripping accelerates outside-in quenching. Chemical clocks can show halted dilution but cannot identify whether the cutoff was caused by halo stripping, virial shock heating, AGN heating, or cosmological supply decline without environmental evidence.

### 7.6 Model calibration degeneracy

A SAM can match passive fractions by retuning feedback, cooling, stripping, and orphan survival. A parametric quenching clock can match the same quenched fraction with different onset radii and timescales. Agreement with one population statistic is therefore not proof of the physical mechanism.

### 7.7 Multiphase and tracer mismatch

H I, H2, ionized gas, hot CGM, and diffuse ionized gas are dynamically and observationally different. A tail visible in Halpha may trace shocks or leaked photons, not local massive-star formation. CO survival does not imply normal star-formation efficiency. A model must predict the phase corresponding to the data being fitted.

## 8. Which model should be used for which question?

![Applicability matrix](assets/20260813_cluster_environment_analytic_models/figure_04_applicability_matrix.png)

*Figure 4. Qualitative scope of model families. "Strong" means the quantity is a direct intended output when adequate inputs are available; it does not mean unbiased or precisely calibrated.*

| Scientific question | Minimum useful model | Required evidence | Honest output |
|---|---|---|---|
| Is the outer gas disk currently vulnerable? | Resolved Gunn-Gott or full vertical restoring force | Sigma_star, Sigma_gas, ICM prior, relative-speed prior | R_strip distribution and vulnerable gas fraction |
| Was the encounter brief or quasi-static? | Koppen finite-pulse model | Orbit/pressure-pulse prior and disk vertical timescale | P_max-impulse posterior |
| Has replenishment been cut off? | McCarthy hot-halo model plus regulator | CGM/hot-halo prior, gas fraction, SFR, metallicity | hot-gas survival and starvation probability |
| Is the galaxy a recent infaller or backsplash object? | Simulation-calibrated projected phase space | R_proj, Delta v_los, cluster mass and membership | broad infall-time/backsplash probability |
| What is time to/since peak RPS? | Tailored time-dependent forward model | resolved H I/CO/Halpha, wind geometry, stellar/SFH clocks | posterior over time and peak pressure |
| What mechanism dominates? | Joint mixture model | gas and stellar morphology/kinematics plus SFH/metallicity | probability for RPS, tide, starvation, and interaction mixtures |
| What happens to the population? | SAM or empirical quenching model | well-defined selection function and population statistics | distributions and trends, not individual histories |

## 9. A recommended MAUVE model

### 9.1 Level 0: evidence tags

Before fitting, encode independent observed features:

- gas truncation relative to the old stellar disk;
- one-sided H I/CO/Halpha tail;
- leading-edge compression;
- extraplanar gas;
- stellar asymmetry or bridge;
- gas-stellar kinematic decoupling;
- radial SFH truncation age;
- H II-selected metallicity gradient;
- diffuse-ionized/shock/AGN contamination;
- projected phase-space coordinates and subcluster association.

These tags define which likelihood terms are valid. For example, a Seyfert cone or shock-dominated tail must not be used as a simple Halpha-SFR likelihood.

### 9.2 Level 1: resolved restoring-pressure map

For every spatial bin or radial annulus, compute

```text
P_restore(R,phi) = Sigma_gas(R,phi) max_z |dPhi/dz|.
```

Use the observed stellar mass map and, where possible, the gas contribution and a dark-halo prior. Sample over the CO-to-H2 conversion factor, H I opacity, stellar mass-to-light ratio, vertical scale heights, distance, and inclination. This produces a susceptibility map rather than a single radius.

### 9.3 Level 2: orbit and pressure pulse

Define latent variables

```text
theta_orbit = {z_los, v_tan, orbital phase, subcluster, ICM bulk velocity}
theta_wind  = {inclination to disk, azimuth, P_max, tau_pulse, pulse shape}.
```

The cluster potential and X-ray gas model supply priors. Propagate the orbit to P_ram(t). For a Virgo-like non-relaxed environment, use a mixture of ICM components plus lognormal or heavy-tailed local scatter. The likelihood should include the projected position and line-of-sight velocity but should not force a unique 3D orbit.

### 9.4 Level 3: gas removal and response

For each gas phase:

- use the long-pulse restoring condition and short-pulse impulse condition;
- include a finite removal time;
- allow fallback for displaced but bound gas;
- treat H I and H2 separately;
- allow compression to change molecular fraction and SFR before removal;
- include hot-halo removal and reduced inflow.

The local regulator equations then evolve gas, SFR, and metallicity. The stripping term should carry the metallicity of the removed zone rather than the galaxy-wide mean.

### 9.5 Level 4: joint likelihood

A schematic posterior is

```text
p(theta | data) proportional to
    L_HI morphology,kinematics
  x L_CO morphology,kinematics
  x L_Halpha HII-only morphology,kinematics
  x L_SFH radial ages
  x L_Z gas metallicity profile
  x L_stars morphology,kinematics
  x L_phase-space
  x priors(cluster, orbit, disk, feedback).
```

The products are schematic: correlated measurements require a covariance model and should not be double-counted. Posterior predictive maps are essential. If a parameter set fits the H I truncation but fails the CO ridge or stellar ages, the tension should remain visible.

### 9.6 Recommended outputs

For each galaxy, report medians and credible intervals for:

- current and peak P_ram;
- pressure impulse;
- stripping radius by gas phase;
- fraction removed, displaced-bound, and remaining;
- time to/since peak pressure;
- time since first infall and backsplash probability;
- hot-halo survival fraction and inflow suppression;
- radial quenching onset and e-folding time;
- posterior probabilities for RPS-dominated, tide-dominated, starvation-dominated, mixed, and unconstrained cases.

The stage label `pre-peak / near-peak / post-peak` should be derived from this posterior, for example by integrating the probability over time relative to peak, and accompanied by an `unconstrained` option.

### 9.7 Minimal and ambitious data configurations

**Minimal:** optical stellar mass map, Halpha map with excitation masking, global H I mass/deficiency, projected cluster coordinates, line-of-sight velocity, disk inclination/PA, and an X-ray ICM prior. This supports a probabilistic susceptibility and broad stage analysis.

**Good:** resolved H I, CO, Halpha, stellar population ages, gas and stellar velocity fields, resolved H II metallicity, and diffuse-ionized/shock masks. This supports a finite-pulse and response model.

**Best:** the above plus cluster distance constraints, radio continuum/polarization, UV, deep old-stellar imaging, and a substructured ICM velocity/density model. This is the regime where a tailored dynamical inversion can be strongly informative.

## 10. Practical ranking of reviewed models for MAUVE

### Tier A - implement first

1. **Resolved generalized Gunn-Gott map** - fastest falsifiable check and direct connection to MUSE plus H I/CO.
2. **Koppen peak-pressure/impulse treatment** - provides the missing duration axis and is already demonstrated on Virgo galaxies.
3. **Projected-phase-space priors from Oman/Rhee/Yoon** - useful only as priors, never as exact clocks.
4. **Non-equilibrium local gas regulator with stripping** - connects gas loss to SFR and metallicity without assuming equilibrium.

### Tier B - add for physical completeness

5. **McCarthy hot-halo finite-time stripping** - required for starvation and delayed quenching.
6. **Vega-Martinez or observation-based ICM pressure priors** - useful for missing or incomplete X-ray profiles.
7. **Tidal-radius and stellar-disturbance likelihood** - prevents gas-only models from absorbing interaction-driven systems.

### Tier C - calibration and population context

8. **SAM comparisons (SAG, GAEA, L-Galaxies)** - define plausible priors and population-level checks.
9. **Empirical delayed-then-rapid clocks** - compare inferred SFHs with known population timescales.
10. **Tailored sticky-particle/hydrodynamic libraries** - use for high-information benchmark galaxies, not as the first model for all targets.

## 11. Critical uncertainties for Virgo and MAUVE

- **Distance and line-of-sight placement:** Virgo has substantial depth and substructure; projected M87-centric radius is not true radius.
- **ICM dynamics:** a static beta model misses shocks, sloshing, subcluster gas, and bulk velocity.
- **Galaxy interactions:** close pairs such as NGC 4567/4568 require a joint tide-plus-RPS model.
- **Initial conditions:** field scaling relations do not reconstruct the exact pre-infall H I/H2 disk or metallicity gradient.
- **CO conversion factor:** X_CO may vary with metallicity, pressure, and excitation, affecting both restoring force and depletion time.
- **Diffuse ionized gas:** shocks, mixing layers, leaked photons, and AGN can bias Halpha SFR and gas-phase metallicity diagnostics.
- **Feedback:** stellar feedback can make gas easier to strip; compression can make it denser and temporarily harder to remove.
- **Fallback:** extraplanar gas is not synonymous with permanently lost gas.
- **Model discrepancy:** even the best tailored Virgo models do not fit all phases and spatial scales simultaneously.
- **Selection:** jellyfish and optically asymmetric samples emphasize visible, often high-pressure events and do not represent all environmental evolution.

## 12. Falsifiable tests enabled by MAUVE

1. **Restoring-force residual test:** after conditioning on P_ram/P_restore, do H I truncation, molecular fraction, and Halpha suppression follow the predicted radial ordering?
2. **Pulse-regime test:** at fixed peak-pressure estimate, do galaxies with larger orbit-conditioned impulse show larger stripped fractions?
3. **Hot-before-cold test:** do galaxies outside the core show suppressed inferred inflow/CGM proxies while retaining relatively normal cold disks?
4. **Gas-before-SFR test:** is H I deficiency enhanced before resolved sSFR suppression, as predicted by SAM comparisons?
5. **Metallicity-selection test:** does the integrated metallicity enhancement disappear after controlling for radial coverage and H II selection?
6. **Tail clock test:** do UV/optical stellar ages in stripped tails agree with the dynamical time since gas displacement?
7. **Tide-versus-RPS test:** can the old stellar velocity field predict which gas asymmetries cannot be explained by hydrodynamics alone?
8. **Posterior predictive stage test:** do pre-peak posterior systems show leading-edge compression, while post-peak systems show stronger truncation/fallback and older outer-disk quenching ages?

## 13. Research gaps

1. There is no widely adopted open Bayesian code that jointly fits a substructured cluster, time-dependent RPS, multiphase gas, radial SFH, and metallicity for individual nearby galaxies.
2. Analytical models usually predict gas mass or radius, not the emission-line observables used to infer them.
3. ICM magnetic fields, anisotropic conduction, viscosity, and cosmic rays are rarely propagated into individual-history uncertainty.
4. SAMs often simplify or omit spatial metallicity gradients in stripped gas.
5. Projected-phase-space mappings are simulation- and definition-dependent; infall may mean crossing 2.5 Rvir, R200, or another boundary.
6. Tailored model times often come from a discrete library and visual/multi-observable matching rather than a fully sampled posterior.
7. Joint tides plus RPS remain especially underconstrained for interacting cluster pairs.
8. Molecular-gas compression, conversion, and removal are less cleanly described analytically than diffuse H I stripping.

## 14. Conclusions

The field has a strong analytical backbone, but its components answer different questions. Gunn-Gott remains the correct first calculation for a cold disk. Koppen's finite-pulse formulation is the most useful analytical extension for orbital history. McCarthy provides the standard finite-time hot-halo prescription. SAMs show how these pieces affect populations and expose the importance of gradual removal, shielding, and environment beyond R200. Projected phase space supplies statistical history priors, while gas regulators and chemical clocks constrain the response after gas supply changes. Tailored multiwavelength forward models are necessary for individual-galaxy timing.

The central methodological conclusion is simple: **no analytical model can uniquely determine the cluster history of a galaxy from projected phase space or a single gas diagnostic.** A robust MAUVE analysis should combine complementary equations in a probabilistic forward model, explicitly marginalize latent orbit and ICM variables, and use resolved gas, stars, SFH, metallicity, and kinematics to break degeneracies. The result should be a distribution over histories and mechanisms, backed by posterior predictive maps and explicit model failures.

## 15. Reproducibility and disclosure

This report was produced with AI-assisted literature searching, metadata checking, synthesis, figure generation, and PDF layout. Scientific claims were restricted to source-supported statements, while proposed MAUVE equations and the implementation roadmap are explicitly identified as synthesis. Crossref metadata resolution verifies bibliographic identity, not scientific correctness. Publisher/arXiv/ADS records were spot-checked for the central model claims. The search is broad and source-traceable but cannot guarantee universal completeness across all historical proceedings, theses, non-English literature, or papers that use environmental prescriptions without searchable terminology.

The following appendices are generated directly from the verified source registry.
