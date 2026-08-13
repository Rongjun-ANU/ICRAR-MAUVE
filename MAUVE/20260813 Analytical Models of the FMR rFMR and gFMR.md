---
title: "Analytical Models of the FMR, rFMR, and gFMR"
subtitle: "A physics-first literature review and a two-driver route to positive SFR–metallicity correlations at fixed mass"
date: "13 August 2026"
lang: en-AU
toc: true
toc-depth: 3
---

# Executive summary

This review asks a narrower question than a general review of the mass–metallicity relation: **which astronomical papers provide an analytical or semi-analytical explanation of the three-variable relation between gas-phase metallicity, mass, and star formation or gas content?** It distinguishes three related objects:

1. the **global FMR**, $Z_{\mathrm{g}}(M_\star,{\mathrm{SFR}})$;
2. the **resolved FMR (rFMR)**, usually $Z_{\mathrm{g}}(\Sigma_\star,\Sigma_{\mathrm{SFR}})$ within galaxies; and
3. the **gaseous FMR (gFMR)** in the specific recent sense $Z_{\mathrm{g}}(M_\star,M_{\mathrm{g}})$, not merely “a gas-phase FMR.”

The main conclusions are:

* There is no single canonical analytical FMR model. Most physical models are variants of one gas-and-metal continuity system. They differ in which quantity is assumed to fluctuate, which quantities are held fixed, and whether equilibrium is imposed.
* The classical negative SFR–$Z_{\mathrm{g}}$ residual at fixed $M_\star$ follows naturally when **metal-poor inflow is the dominant short-timescale driver**. A rise in inflow first raises gas mass and SFR and simultaneously dilutes the gas. Dayal et al. (2013), Forbes et al. (2014), Lilly et al. (2013), and Wang et al. (2026) provide distinct analytical versions of this idea. ([Dayal et al. 2013](https://arxiv.org/abs/1202.4770))<!--ref:dayal2013--><!--anchor:section:2--> ([Forbes et al. 2014](https://arxiv.org/abs/1311.1509))<!--ref:forbes2014--><!--anchor:section:2--> ([Lilly et al. 2013](https://arxiv.org/abs/1303.5059))<!--ref:lilly2013--><!--anchor:abstract--> ([Wang et al. 2026](https://arxiv.org/abs/2608.04784))<!--ref:wang2026--><!--anchor:section:4-->
* A positive SFR–$Z_{\mathrm{g}}$ residual at fixed mass is not pathological. It arises when **star-formation-efficiency variations**, rather than inflow variations, dominate: SFR rises immediately at nearly fixed gas supply, while metal production drives $Z_{\mathrm{g}}$ upward after a finite response time. Wang & Lilly (2021) derive this sign analytically; Huang et al. (2026) express the resolved sign as competition between depletion-time and inflow-time fluctuations. ([Wang & Lilly 2021](https://arxiv.org/abs/2009.01935))<!--ref:wanglilly2021--><!--anchor:section:2.3--> ([Huang et al. 2026](https://arxiv.org/abs/2605.31412))<!--ref:huang2026--><!--anchor:section:5.2-->
* The requested modifications of **Wang et al. (2026)** and **Forbes et al. (2014)** are therefore possible, but the physically clean modification is a two-driver regulator: add a stochastic or coherent efficiency field $\epsilon(t)$ alongside the inflow $\Phi(t)$. This is an extension of those two particular models, not a newly discovered mechanism.
* A reproducible numerical experiment in this report demonstrates the sign change. With all other parameters fixed, an inflow-only stochastic regulator gives $r(\Delta\log {\mathrm{SFR}},\Delta\log Z_{\mathrm{g}})=-0.38$; an efficiency-only regulator gives $r=+0.58$; and a mixed case in which efficiency fluctuations dominate gives $r=+0.25$. Across 12 stochastic realisations per amplitude ratio, the median correlation first becomes positive near $\sigma_{\ln\epsilon}/\sigma_{\ln\Phi}=1.125$ for the chosen response and coherence times. This threshold is illustrative, not universal.
* For MAUVE, the discriminating observation is not the sign by itself. A positive branch caused by efficiency fluctuations should be associated with shorter depletion time or higher SFE at fixed $\Sigma_\star$, whereas an inflow-driven negative branch should be associated with larger gas fraction and dilution. CO/HI information, response-time proxies, consistent abundance calibrations, and explicit DIG/shock control are therefore essential.

# 1. Scope, terminology, and search method

## 1.1 What counts as an analytical model here

Included papers do at least one of the following:

* derive a closed-form or low-dimensional differential-equation relation among metallicity, stellar mass, gas mass, SFR, inflow, outflow, or efficiency;
* derive the sign or response amplitude of metallicity–SFR residuals;
* construct a resolved annular or local regulator with an explicit analytical metallicity solution; or
* provide an empirical analytical surface that is often called an FMR or rFMR and therefore must be distinguished from a causal model.

Hydrodynamic simulations and full semi-analytic galaxy-formation codes are used only as context. They can test whether an FMR emerges, but they are not “analytical models” in the narrow sense used here. The review also separates **foundational regulator papers** from papers that explicitly model the FMR.

## 1.2 Search strategy and limits

The online search was conducted on 13 August 2026. It combined exact-phrase and concept searches for “fundamental metallicity relation,” “analytical/analytic model,” “gas regulator,” “resolved FMR,” “local metallicity–SFR relation,” “gaseous FMR,” “gas fraction metallicity,” “positive SFR metallicity correlation,” and citation chaining from the principal model papers. Primary arXiv records, journal pages, and author-hosted manuscripts were preferred. The equations and model claims below were checked against the papers, not against search-result summaries.

This is a **deep, structured scoping review**, not a formally registered systematic review. “Any paper” cannot be proven literally exhaustive because terminology is inconsistent: some relevant papers never use FMR, rFMR, or gFMR in the title; “semi-analytic” may mean either algebra or a numerical galaxy-formation framework; and gFMR has multiple informal meanings. The evidence map therefore labels direct models, foundations, empirical surfaces, and contextual numerical work instead of inflating all of them into one category.

## 1.3 Definitions that prevent category errors

* **Global FMR:** a relation or surface in $(M_\star,{\mathrm{SFR}},Z_{\mathrm{g}})$. Mannucci et al. (2010) popularised this form; Lara-López et al. (2010) proposed a related plane. These are observational parameterisations, not continuity-equation explanations. ([Mannucci et al. 2010](https://arxiv.org/abs/1005.0006))<!--ref:mannucci2010--><!--anchor:abstract--> ([Lara-López et al. 2010](https://arxiv.org/abs/1005.0509))<!--ref:laralopez2010--><!--anchor:abstract-->
* **rFMR:** a relation among local surface-density quantities, most commonly $(\Sigma_\star,\Sigma_{\mathrm{SFR}},Z_{\mathrm{g}})$. The word “resolved” does not by itself imply a causal local model; a fitted surface and an annular conservation law are different evidential objects.
* **gFMR:** in Wang et al. (2026), the “gaseous FMR” is $(M_\star,M_{\mathrm{g}},Z_{\mathrm{g}})$. Older HI- or H$_2$-based observational relations are close conceptual ancestors, but they do not necessarily use that name. ([Wang et al. 2026](https://arxiv.org/abs/2608.04784))<!--ref:wang2026--><!--anchor:section:4-->
* **The sign at fixed mass:** throughout this report, “positive” means $\partial Z_{\mathrm{g}}/\partial\log {\mathrm{SFR}}>0$ at fixed $M_\star$ globally or fixed $\Sigma_\star$ locally. A projected mass trend, radial gradient, or correlation driven by mixing different masses is not sufficient.

# 2. The shared analytical skeleton

Most of the physical literature can be written as a one-zone or one-patch gas regulator. Let $M_{\mathrm{g}}$ be gas mass, $Z$ its metal mass fraction, $\Phi$ the inflow rate, $\epsilon\equiv {\mathrm{SFR}}/M_{\mathrm{g}}$ the star-formation efficiency, $R$ the returned fraction, $\eta$ the outflow mass-loading factor, $y$ the yield in the adopted convention, and $Z_{\mathrm{in}}$ the inflow metallicity. For a well-mixed outflow,

$$
\frac{{\mathrm{d}}M_{\mathrm{g}}}{{\mathrm{d}}t}=\Phi-(1-R+\eta)\epsilon M_{\mathrm{g}},
$$

$$
\frac{{\mathrm{d}}(M_{\mathrm{g}}Z)}{{\mathrm{d}}t}=y\epsilon M_{\mathrm{g}}+Z_{\mathrm{in}}\Phi
-Z(1-R+\eta)\epsilon M_{\mathrm{g}}.
$$

Combining them removes the explicit well-mixed outflow term from the instantaneous abundance equation:

$$
\frac{{\mathrm{d}}Z}{{\mathrm{d}}t}=y\epsilon+(Z_{\mathrm{in}}-Z)\frac{\Phi}{M_{\mathrm{g}}}
=\frac{Z^\dagger-Z}{\tau_\Phi},
$$

where

$$
\tau_\Phi\equiv\frac{M_{\mathrm{g}}}{\Phi},\qquad
Z^\dagger\equiv Z_{\mathrm{in}}+y\frac{\epsilon M_{\mathrm{g}}}{\Phi}
=Z_{\mathrm{in}}+y\frac{{\mathrm{SFR}}}{\Phi}.
$$

This compact form exposes the driver competition. A first-order perturbation of the target metallicity gives

$$
\delta Z^\dagger=(Z^\dagger-Z_{\mathrm{in}})
\left(\delta\ln {\mathrm{SFR}}-\delta\ln\Phi\right)+\delta Z_{\mathrm{in}}.
$$

Consequently:

* **inflow-driven excursion:** $\delta\ln\Phi$ is large; dilution is immediate, while gas consumption and enrichment lag, producing the usual negative residual;
* **efficiency-driven excursion:** $\Phi$ is nearly fixed while $\epsilon$ and SFR rise; the target metallicity rises, producing a positive residual after the metal-response filter;
* **enriched/recycled inflow:** a positive $\delta Z_{\mathrm{in}}$ correlated with SFR can also make the residual positive;
* **variable mass loading:** in equilibrium, $Z\sim Z_{\mathrm{in}}+y/(1-R+\eta)$, so an anti-correlation of $\eta$ with SFR can yield a positive residual, though this needs an independent physical justification.

The key lesson is that the sign diagnoses the **dominant perturbation and its timing**, not merely the existence of “self-regulation.” Wang & Lilly (2021) and Huang et al. (2026) derive this logic explicitly at global and resolved scales. ([Wang & Lilly 2021](https://arxiv.org/abs/2009.01935))<!--ref:wanglilly2021--><!--anchor:section:2--> ([Huang et al. 2026](https://arxiv.org/abs/2605.31412))<!--ref:huang2026--><!--anchor:section:5.2-->

# 3. Evidence map

| Model family | Representative papers | Scale / variables | Mathematical role | Typical fixed-mass SFR–$Z$ sign |
|---|---|---|---|---|
| Empirical FMR surfaces | Mannucci et al. 2010; Lara-López et al. 2010; Curti et al. 2020 | Global $M_\star$, SFR, $Z$ | Polynomial, plane, or saturating fit | Usually negative at low/intermediate mass; fitted rather than predicted |
| Equilibrium/bathtub foundations | Davé et al. 2012; Feldmann 2013; Lilly et al. 2013 | Global gas, stars, metals | Continuity equations and equilibrium/slow evolution | Negative if high SFR traces high gas supply; none in strict equilibrium without a residual driver |
| Inflow/outflow analytic FMR | Dayal et al. 2013; Hunt et al. 2016b | Global FMR | Closed chemical-evolution solution with rates tied to SFR | Negative below the high-mass saturation regime |
| Stochastic non-equilibrium regulator | Forbes et al. 2014 | Global residuals at fixed mass | Linear response to stochastic inflow | Negative in the published inflow-only model |
| Generalised regulator solutions | Pipino et al. 2014; Peng & Maiolino 2014 | Global | Time-dependent or quasi-steady analytic solutions | Depends on driver; commonly used to explain negative FMR |
| Non-equilibrium calibrated evolution | Lin & Zu 2023 | Global FMR | ODE model along mean SF histories with variable $\eta$ | Reproduces observed global surface; not a simple residual impulse model |
| Inflow-driven closed-form gFMR | Wang et al. 2026 | $M_\star$, $M_{\mathrm{g}}$, $Z$ and mapped FMR | Closed-form ideal model plus cosmological calibration | gFMR has $Z\uparrow$ with $M_\star/M_{\mathrm{g}}$; mapped FMR is negative for its calibrated efficiency law |
| Local/annular equilibrium | Carton et al. 2015 | Radial $M_{\mathrm{g}}/M_\star$, sSFR, $Z$ | Annular regulator metallicity profile | Usually negative when sSFR represents gas-rich inflow |
| Stellar-fraction metallicity relation | Ascasibar et al. 2015 | Global and kpc patches, $M_\star/M_{\mathrm{g}}$ | Closed form linking gas fraction to metallicity | Indirect; maps to negative FMR if SFR tracks gas mass |
| Local-to-global residual mapping | Sánchez Almeida & Sánchez-Menguiano 2019 | Local and integrated residuals | Analytical integration of an assumed local slope | Preserves the assumed local sign under stated conditions |
| Driver-response model | Wang & Lilly 2021 | Global or local fluctuations | Sinusoidal linear response and amplitudes | Inflow negative; efficiency positive |
| Empirical rFMR surface | Baker et al. 2023 | $\Sigma_\star$, $\Sigma_{\mathrm{SFR}}$, $Z$ | Saturating fitted surface and partial correlations | Negative on average in their sample |
| Spatial transport foundation | Sharda et al. 2021 | Radial metallicity profiles | Advection–diffusion–production equation | No unique local SFR-residual sign |
| Two-timescale resolved model | Huang et al. 2026 | $\Sigma_\star$, $\Sigma_{\mathrm{SFR}}$, $Z$ | Perturbation and temporal-filter solution | Positive if efficiency variations dominate; negative if inflow dominates |

The table intentionally does not list every simulation that reproduces an FMR. Fontanot et al. (2021), Yates et al. (2012, 2014), EAGLE, FIRE, and IllustrisTNG studies are important checks, but their galaxy-formation frameworks are numerically semi-analytic or hydrodynamic rather than the closed analytical models requested here.

# 4. Global FMR models in detail

## 4.1 Observational surfaces: Mannucci, Lara-López, and Curti

Mannucci et al. (2010) defined a low-scatter surface in stellar mass, SFR, and oxygen abundance and introduced $\mu_\alpha=\log M_\star-\alpha\log {\mathrm{SFR}}$, with $\alpha\simeq0.32$ for their calibration and sample. Lara-López et al. (2010) fitted a plane-like relation. Curti et al. (2020) later provided direct-$T_{\mathrm{e}}$-anchored empirical parameterisations. These formulae are analytically evaluable, but they are **descriptions of a selected and calibrated data set**, not derivations of inflow, enrichment, or outflow physics. ([Mannucci et al. 2010](https://arxiv.org/abs/1005.0006))<!--ref:mannucci2010--><!--anchor:section:4--> ([Lara-López et al. 2010](https://arxiv.org/abs/1005.0509))<!--ref:laralopez2010--><!--anchor:abstract--> ([Curti et al. 2020](https://arxiv.org/abs/1910.00597))<!--ref:curti2020--><!--anchor:section:4-->

Their continuing value is as a target that physical models must reproduce. Their limitation is equally important: the fitted coefficient of SFR changes with abundance indicator, aperture, S/N selection, mass range, and the assumed functional shape. A polynomial surface cannot decide whether the third-parameter dependence is dilution, SFE, mass loading, enriched recycling, or a measurement covariance.

## 4.2 Davé et al. (2012), Feldmann (2013), and the equilibrium foundation

The “equilibrium” or “bathtub” view balances gas supply with star formation and outflow. Davé et al. (2012) showed how ejective feedback, preventive feedback, and wind recycling set galaxy gas fractions, SFRs, and metallicities. Feldmann (2013) derived equilibrium scaling relations for gas and stellar fractions, sSFR, and metallicity. These papers supply the physical grammar later FMR models use. ([Davé et al. 2012](https://arxiv.org/abs/1108.0426))<!--ref:dave2012--><!--anchor:section:2--> ([Feldmann 2013](https://arxiv.org/abs/1212.2973))<!--ref:feldmann2013--><!--anchor:section:2-->

A strict steady state is nevertheless insufficient to explain a residual covariance by itself: if all parameters are fixed, every object sits on one equilibrium point. An FMR appears only after allowing mass-dependent parameters, non-zero gas growth, different evolutionary stages, or fluctuations around equilibrium. This distinction is sometimes obscured when “gas regulation explains the FMR” is stated without naming the variable that produces the scatter.

## 4.3 Dayal, Ferrara & Dunlop (2013): a direct closed FMR

Dayal et al. construct one of the cleanest direct analytical FMR models. They assume

$$
\dot M_\star=\psi=\epsilon_\star M_{\mathrm{g}},\qquad
\dot M_{\mathrm{g}}=-(1-R)\psi+(a-w)\psi,
$$

$$
M_{\mathrm{g}}\dot X=y(1-R)\psi-aX\psi,
$$

where $a\psi$ and $w\psi$ are pristine inflow and enriched outflow and $X$ is the oxygen mass fraction. The solution may be expressed as

$$
X(\psi)=\frac{y(1-R)}{a}\left[1-\mu^{-\alpha}\right],\qquad
\mu=\frac{M_{\mathrm{g}}}{M_{{\mathrm{g}},0}}=\frac{\psi}{\epsilon_\star M_{{\mathrm{g}},0}},
$$

with $\alpha=a/(R-1+a-w)$. Their mass-dependent fitted inflow and outflow efficiencies reproduce the local Mannucci surface and, without an additional fit, the gas-fraction–metallicity relation. Low-mass galaxies show decreasing metallicity with increasing SFR; massive galaxies saturate because enrichment and dilution approach balance while outflows weaken. ([Dayal et al. 2013](https://academic.oup.com/mnras/article/430/4/2891/1101320))<!--ref:dayal2013--><!--anchor:section:2-->

**Strengths.** It is explicit, compact, and directly connects the FMR to inflow, outflow, and yield. It also predicts a second observable relation.

**Limitations.** Inflow and outflow are imposed proportional to SFR; $a(M_\star)$ and $w(M_\star)$ are calibrated to the FMR; instantaneous recycling/mixing and a constant efficiency are assumed; and time-correlated residuals are not modelled. It is therefore more a family of deterministic evolutionary tracks than a dynamical explanation of short-timescale scatter.

## 4.4 Lilly et al. (2013): the gas regulator

Lilly et al. derive metallicity, gas content, and star formation in a regulator whose reservoir is fed by inflow and drained by star formation and outflow. Their key advance for FMR work is allowing the gas reservoir to grow rather than demanding exact gas-mass equilibrium, while adopting a slowly varying or equilibrium-like abundance solution. Metallicity becomes a function of efficiency, mass loading, gas-to-stellar ratio, and sSFR. If efficiency and loading do not explicitly evolve, changes in the main sequence and gas content can generate an approximately redshift-invariant $Z(M_\star,{\mathrm{SFR}})$ surface. ([Lilly et al. 2013](https://arxiv.org/abs/1303.5059))<!--ref:lilly2013--><!--anchor:abstract-->

The model is conceptually powerful because it links the FMR to independently measurable gas-regulator quantities. Its danger is approximation mixing: setting $\dot Z\simeq0$ while allowing $\dot M_{\mathrm{g}}\ne0$ is not the same as solving the fully time-dependent metal equation. Forbes et al. explicitly criticise that mixed equilibrium, and Pipino et al. quantify where the approximation is accurate. ([Forbes et al. 2014](https://arxiv.org/abs/1311.1509))<!--ref:forbes2014--><!--anchor:section:4.2--> ([Pipino et al. 2014](https://arxiv.org/abs/1403.6146))<!--ref:pipino2014--><!--anchor:section:3-->

## 4.5 Forbes et al. (2014): stochastic inflow and covariance

Forbes et al. focus on **scatter at fixed mass**, making the paper especially relevant here. Their gas equation is

$$
\dot M_{\mathrm{g}}=\dot M_{\mathrm{ext}}(t)-\frac{M_{\mathrm{g}}}{t_{\mathrm{loss}}},
\qquad
\dot M_{\mathrm{ext}}=\exp[\mu+\sigma x(t)],
$$

where a new Gaussian inflow deviate is drawn on a coherence time $t_{\mathrm{coherence}}$. Gas loss combines long-lived star formation and outflow,

$$
\frac{M_{\mathrm{g}}}{t_{\mathrm{loss}}}=(f_R+\eta){\mathrm{SFR}}.
$$

With $\Psi=(M_{\mathrm{g}}/t_{\mathrm{loss}})e^{-\mu}$ and $\tau=t/t_{\mathrm{loss}}$, the gas response is

$$
\frac{{\mathrm{d}}\Psi}{{\mathrm{d}}\tau}=-\Psi+e^{\sigma x}.
$$

They add a metal-continuity equation with inflow metallicity and an effective yield $q$. After an upward inflow step, metallicity drops immediately through dilution while gas mass and SFR rise more slowly; after a downward step, continued stellar processing raises metallicity while SFR falls. The anti-correlation is strongest when $t_{\mathrm{coherence}}/t_{\mathrm{loss}}$ is of order unity. The same two dimensionless parameters—accretion scatter $\sigma$ and the coherence-to-loss-time ratio—control the scatter and covariance of multiple scaling relations. ([Forbes et al. 2014](https://www.mso.anu.edu.au/~krumholz/publications/2014/forbes14a.pdf))<!--ref:forbes2014--><!--anchor:page:4-->

**Strengths.** The model predicts a covariance, not merely a mean surface; retains time dependence; and ties the FMR sign to an observable response-time hierarchy.

**Limitations.** Inflow is the only stochastic driver. Depletion time, mass loading, yield, and inflow metallicity are fixed within a mass bin. The published negative sign is therefore conditional on the model design, not a theorem that every regulator must be negative.

## 4.6 Pipino, Lilly & Carollo (2014) and Peng & Maiolino (2014): time dependence clarified

Pipino et al. derive exact and approximate metallicity–sSFR relations for regulator histories, test the Lilly approximation, and show that it remains accurate to roughly 0.1 dex over a wide sSFR range for their explored models. They also predict stellar metallicity to lag gas metallicity by roughly 0.1–0.2 dex. ([Pipino et al. 2014](https://arxiv.org/abs/1403.6146))<!--ref:pipino2014--><!--anchor:abstract-->

Peng & Maiolino develop a general time-dependent regulator without requiring instantaneous equilibrium and analyse relaxation times and trajectories. It is a foundation for interpreting FMR evolution rather than a fitted FMR surface in its own right. ([Peng & Maiolino 2014](https://arxiv.org/abs/1402.5964))<!--ref:pengmaiolino2014--><!--anchor:abstract-->

Together these papers show why “equilibrium” must be specified: gas-mass equilibrium, abundance equilibrium, statistical equilibrium, and slow evolution are different approximations.

## 4.7 Harwit & Brisbin (2015), Hunt et al. (2016), and other compact extensions

Harwit & Brisbin build an equilibrium accounting model in which infalling gas mixes with native gas and stellar yields constrain the enrichment. Applied to a large SDSS sample, it infers that external accretion supplies a substantial fraction of the gas required for star formation. Its value is a directly interpretable census; its limitations are equilibrium selection and sensitivity to yield, aperture, and abundance calibration. ([Harwit & Brisbin 2015](https://arxiv.org/abs/1412.2436))<!--ref:harwit2015--><!--anchor:abstract-->

Hunt et al. (2016a) provide an empirical “fundamental plane of metallicity,” approximately linear in $\log M_\star$ and $\log {\mathrm{SFR}}$ for their heterogeneous MEGA compilation. Hunt et al. (2016b) then update the Dayal physical model, fitting inflow and outflow scalings and interpreting high-redshift offsets mainly through greater accretion and gas content. The two papers should not be conflated: one is a data reduction surface; the other is a chemical-evolution model. ([Hunt et al. 2016a](https://arxiv.org/abs/1608.05417))<!--ref:hunt2016a--><!--anchor:abstract--> ([Hunt et al. 2016b](https://arxiv.org/abs/1608.05418))<!--ref:hunt2016b--><!--anchor:abstract-->

Magrini et al. (2012) compare multi-phase chemical-evolution modes for metal-poor starbursts and scaling relations. Zahid et al. (2012) link empirically inferred enrichment histories, SFR histories, and oxygen outflows. These broaden the chemical-evolution context, but neither is as minimal or directly diagnostic of FMR residual covariance as Forbes et al. ([Magrini et al. 2012](https://doi.org/10.1111/j.1365-2966.2012.22055.x))<!--ref:magrini2012--><!--anchor:abstract--> ([Zahid et al. 2012](https://arxiv.org/abs/1207.5509))<!--ref:zahid2012--><!--anchor:abstract-->

## 4.8 Lin & Zu (2023): a calibrated non-equilibrium chemical-evolution model

Lin & Zu construct a non-equilibrium chemical-evolution model (NE-CEM) that follows an average star-formation history and a time-dependent mass-loading factor. For pristine inflow, their abundance evolution includes

$$
\dot Z_{\mathrm{O}}=\frac{\dot M_{\mathrm{O}}}{M_{\mathrm{g}}}-rac{\dot M_{\mathrm{g}}}{M_{\mathrm{g}}}Z_{\mathrm{O}}.
$$

The functional dependence of $\eta(M_\star,{\mathrm{sSFR}})$ is informed by EAGLE and then constrained with SDSS. The model reproduces the local metallicity–mass–SFR relation and interprets its evolution as coherent non-equilibrium enrichment along galaxy growth tracks. ([Lin & Zu 2023](https://arxiv.org/abs/2212.01402))<!--ref:linzu2023--><!--anchor:section:2-->

This is richer than a one-parameter bathtub and directly useful for outflow inference. The trade-off is identifiability: mean SF histories, the EAGLE-motivated form, yield choices, and the inferred gas history all contribute. It does not by itself isolate whether short-timescale positive residuals come from SFE, recycled inflow, or a population-selection effect.

## 4.9 Wang et al. (2026): inflow-driven FMR and the gFMR

Wang et al. present the newest and most explicit analytical treatment in this search. It is currently an arXiv v1 preprint submitted 5 August 2026, so its claims should be treated as pre-peer-review. Their ideal model adopts pristine inflow and constant $\Phi$, $\epsilon$, and $\eta$:

An adjacent single-author paper, Wang (2026), introduced the same non-equilibrium inflow-driven regime while studying the **galaxy size–stellar-metallicity relation**. Its regulator shows why SFE can affect stellar metallicity before long-term equilibrium and motivates some of the later FMR algebra, but it is not itself an FMR/rFMR/gFMR paper. It is included here to prevent the two 2026 Wang works from being conflated. ([Wang 2026](https://arxiv.org/abs/2510.02573))<!--ref:wangsize2026--><!--anchor:section:5-->

$$
\dot M_{\mathrm{g}}=\Phi-(1-R+\eta)\epsilon M_{\mathrm{g}},
$$

$$
\frac{{\mathrm{d}}(M_{\mathrm{g}}Z_{\mathrm{g}})}{{\mathrm{d}}t}
=y\epsilon M_{\mathrm{g}}-Z_{\mathrm{g}}(1-R+\eta)\epsilon M_{\mathrm{g}}.
$$

For $\tau_{\mathrm{eq}}=[(1-R+\eta)\epsilon]^{-1}$, they derive

$$
M_{\mathrm{g}}=\Phi\tau_{\mathrm{eq}}(1-e^{-t/\tau_{\mathrm{eq}}}),
$$

$$
Z_{\mathrm{g}}=y\epsilon\left[\tau_{\mathrm{eq}}
-\frac{t}{e^{t/\tau_{\mathrm{eq}}}-1}\right],
$$

$$
M_\star=(1-R)\epsilon\tau_{\mathrm{eq}}\Phi
\left[t-\tau_{\mathrm{eq}}(1-e^{-t/\tau_{\mathrm{eq}}})\right].
$$

In the early inflow-driven limit, $M_{\mathrm{g}}\simeq\Phi t$, $Z_{\mathrm{g}}\simeq y\epsilon t/2$, and $M_\star\simeq(1-R)\Phi\epsilon t^2/2$, giving the particularly transparent gFMR

$$
Z_{\mathrm{g}}=\frac{y}{1-R}\frac{M_\star}{M_{\mathrm{g}}}.
$$

More generally they derive correction functions of $x=t/\tau_{\mathrm{eq}}$ that connect $Z_{\mathrm{g}}$, $M_{\mathrm{g}}/M_\star$, and $\eta$. Their cosmological model is calibrated to the mass–metallicity relation, main sequence, and stellar-to-halo mass relation, then predicts an FMR. To map gas fraction into SFR they allow an empirical efficiency law $\epsilon\propto M_\star^a{\mathrm{SFR}}^b$, yielding a projected coordinate $\log M_\star-\alpha\log {\mathrm{SFR}}$ with $\alpha=(1-b)/(1+a)$ and a calibrated value near 0.55. ([Wang et al. 2026](https://arxiv.org/pdf/2608.04784))<!--ref:wang2026--><!--anchor:page:12-->

**Strengths.** It cleanly separates the more fundamental gas-fraction relation from the SFR projection; gives closed forms beyond the earliest limit; and embeds them in a cosmological calibration.

**Limitations relevant to this report.** The ideal closed solution fixes efficiency and cannot generate a positive SFR–metallicity residual at fixed mass. The broader population model includes a mean efficiency mapping but deliberately does not model efficiency scatter at fixed mass. Recycling, enriched inflow, and explicitly correlated short-timescale drivers are omitted. Stochastic inflow can generate individual anti-correlated excursions, but a population-level, redshift-invariant FMR requires additional calibrated structure. These limitations are acknowledged in the preprint and define the most natural extension.

# 5. The gFMR: gas content as the more direct coordinate

The physical argument for a gas-based relation predates the label gFMR. Closed-box and leaky-box models have always related abundance to gas fraction; the newer question is whether a tight $Z(M_\star,M_{\mathrm{g}})$ surface is more fundamental than $Z(M_\star,{\mathrm{SFR}})$.

## 5.1 Empirical gas-FMR precursors

Bothwell et al. (2013) reported that HI mass is a strong third parameter of the mass–metallicity relation and argued that the apparent SFR dependence may be secondary to the gas reservoir. Bothwell et al. (2016) extended this logic to molecular gas. These are empirical conditional relations, not closed dynamical solutions; gas-selection limits and conversion factors matter. ([Bothwell et al. 2013](https://arxiv.org/abs/1304.4940))<!--ref:bothwell2013--><!--anchor:abstract--> ([Bothwell et al. 2016](https://arxiv.org/abs/1507.01004))<!--ref:bothwell2016--><!--anchor:abstract-->

Their causal interpretation is nevertheless compelling: gas mass is a state variable in the continuity equations, whereas SFR is a rate inferred from that state through an efficiency law. Replacing gas by SFR inevitably imports depletion-time scatter.

## 5.2 Ascasibar et al. (2015): stellar fraction and metallicity

Ascasibar et al. derive a general relation between metal abundance and the stellar-to-gas ratio, applicable globally and on kpc scales. In their notation the abundance-like variable $\chi$ has the rational form

$$
\chi=\frac{(1-\bar\epsilon_w)\bar\chi_{\mathrm{SN}}
\left[\frac{R}{1-R}\frac{M_\star}{M_{\mathrm{g}}}\right]}
{1+\Upsilon\left[\frac{R}{1-R}\frac{M_\star}{M_{\mathrm{g}}}\right]}.
$$

It is linear in $M_\star/M_{\mathrm{g}}$ in the gas-rich limit and saturates in the gas-poor limit. The relation is argued to be relatively insensitive to the detailed infall and star-formation history and is tested with integrated galaxies and resolved regions in NGC 628 and NGC 3184. ([Ascasibar et al. 2015](https://arxiv.org/abs/1406.6397))<!--ref:ascasibar2015--><!--anchor:section:2-->

This is a direct analytical precursor to the Wang gFMR. Its strengths are the transparent limiting behaviour and global/local applicability. Its parameters combine yields, stellar mass return, and metal loss, so calibration degeneracies remain. It contains no explicit SFR; an FMR sign appears only after imposing a star-formation law.

## 5.3 Wang et al. (2026): why the gFMR is cleaner

The Wang inflow-limit result,

$$
Z_{\mathrm{g}}=\frac{y}{1-R}\frac{M_\star}{M_{\mathrm{g}}},
$$

is the simplest analytical statement found in this search. It says that gas metallicity measures integrated star formation per unit current gas, before equilibrium corrections. The full solution adds two dimensionless functions of $t/\tau_{\mathrm{eq}}$ and the loading factor. Their model finds the gFMR tighter and more directly physical than the standard FMR, while the observed FMR is a projection through $\epsilon(M_\star,{\mathrm{SFR}})$. ([Wang et al. 2026](https://arxiv.org/abs/2608.04784))<!--ref:wang2026--><!--anchor:section:4-->

This hierarchy has an immediate test: at fixed stellar mass, metallicity should anti-correlate more fundamentally with gas mass than with SFR. If depletion-time variations drive a positive SFR–$Z$ branch, the gas-mass residual can remain negative even while the SFR residual turns positive. Wang & Lilly (2021) explicitly identify this difference between gas and SFR responses under efficiency forcing. ([Wang & Lilly 2021](https://arxiv.org/abs/2009.01935))<!--ref:wanglilly2021--><!--anchor:section:2.3-->

# 6. Resolved and local analytical models

## 6.1 Carton et al. (2015): annular equilibrium regulator

Carton et al. apply a Lilly-like regulator independently to radial annuli. With $r_{\mathrm{gas}}=M_{\mathrm{gas}}/M_\star=\epsilon^{-1}{\mathrm{sSFR}}$, their gas continuity gives an inferred inflow

$$
\dot M_{\mathrm{in}}=\left[(1-R)(1+r_{\mathrm{gas}})+\lambda
+\epsilon^{-1}\frac{{\mathrm{d}}\ln r_{\mathrm{gas}}}{{\mathrm{d}}t}\right]{\mathrm{SFR}},
$$

and the local equilibrium abundance is

$$
Z_{\mathrm{eq}}=Z_0+
\frac{y}{1+r_{\mathrm{gas}}+(1-R)^{-1}
\left(\lambda+\epsilon^{-1}{\mathrm{d}}\ln r_{\mathrm{gas}}/{\mathrm{d}}t\right)}.
$$

They fit metallicity profiles of 50 Bluedisk galaxies and show that the radial gas-to-stellar ratio is central to the gradient. ([Carton et al. 2015](https://arxiv.org/abs/1505.02797))<!--ref:carton2015--><!--anchor:section:2.1-->

The model is analytical and resolved, but it is better described as a **local metallicity-profile regulator** than a complete rFMR residual model. It assumes annular equilibrium, zero radial exchange between annuli in the basic form, and a prescribed efficiency; yield, inflow metallicity, and loading are highly degenerate. Because $r_{\mathrm{gas}}=\epsilon^{-1}{\mathrm{sSFR}}$, interpreting an sSFR residual as gas supply or as efficiency requires independent gas data.

## 6.2 Sánchez Almeida & Sánchez-Menguiano (2019): local-to-global mapping

This paper starts from a local residual law

$$
\Delta\log Z_{\mathrm{g}}=m\,\Delta\log\Sigma_{\mathrm{SFR}}
$$

at fixed local stellar surface density. Under small perturbations, an approximately linear star-formation law, and negligible covariance between local SFR residuals and the large-scale metallicity field, the spatial integral preserves the slope:

$$
\Delta\log\langle Z_{\mathrm{g}}\rangle
\simeq m\,\Delta\log\langle{\mathrm{SFR}}\rangle.
$$

They test the correspondence with 736 MaNGA galaxies. ([Sánchez Almeida & Sánchez-Menguiano 2019](https://arxiv.org/abs/1905.05826))<!--ref:sanchezalmeida2019--><!--anchor:section:2-->

This is an important analytical bridge, but it does not derive $m$ from physical conservation laws. Its assumptions can fail when radial gradients, DIG, selection masks, nonlinear star-formation laws, or coherent galaxy-wide events covary with the local residuals. The paper also notes that local slopes can become weakly positive at high mass, a warning against imposing one universal sign.

## 6.3 Wang & Lilly (2021): the decisive driver experiment

Wang & Lilly use the same gas and metal continuity system but force it sinusoidally in two different ways.

* With **time-variable inflow and fixed SFE**, SFR rises after added gas while metallicity is diluted. The metallicity–SFR residual is negative.
* With **time-variable SFE and fixed inflow**, SFR responds immediately to $\epsilon$ and enrichment raises metallicity. The residual is strongly positive over the relevant forcing range, although phase loops and step-like forcing can complicate a single fitted slope.

For a sinusoidal efficiency driver, their linear solution gives a metallicity-to-SFR response-amplitude ratio of the form

$$
\frac{\sigma(\log Z)}{\sigma(\log {\mathrm{SFR}})}
\simeq\frac{1}{\sqrt{1+\xi^2}}
\left(1+\frac{Z_0}{y_{\mathrm{eff}}}\right)^{-1},
$$

where $\xi=2\pi\tau_{\mathrm{dep,eff}}/T_p$ compares the regulator response time with the forcing period. Rapid forcing is filtered out; slow coherent efficiency variations generate a larger positive abundance response. ([Wang & Lilly 2021](https://arxiv.org/abs/2009.01935))<!--ref:wanglilly2021--><!--anchor:page:7-->

This is the clearest prior analytical answer to the requested positive-correlation problem. The interpretation depends on treating spatial scatter as an ensemble of regions at different phases or drivers. On very small scales, ionisation-source age, metal mixing delay, and SFR-tracer averaging times become part of the transfer function.

## 6.4 Baker et al. (2023): an empirical rFMR surface

Baker et al. fit a resolved saturating relation,

$$
12+\log({\mathrm{O/H}})=Z_0-\frac{\gamma}{\phi}
\log\left[1+\left(\frac{\Sigma_\star}{M_0(\Sigma_{\mathrm{SFR}})}\right)^{-\phi}\right],
$$

with $\log M_0=m_0+m_1\log\Sigma_{\mathrm{SFR}}$. They also use $\mu_\alpha=\log\Sigma_\star-\alpha\log\Sigma_{\mathrm{SFR}}$ and find $\alpha\simeq0.54$ and a combined scatter near 0.060 dex. Partial correlations and random forests identify $\Sigma_\star$ as the dominant local predictor and $\Sigma_{\mathrm{SFR}}$ as a secondary anti-correlated predictor; a global $M_\star$ dependence remains. ([Baker et al. 2023](https://arxiv.org/abs/2210.03755))<!--ref:baker2023--><!--anchor:section:3.3-->

This is a strong empirical benchmark, not a causal regulator. Its coefficients are sample- and calibrator-dependent, and a global-mass term means that spaxels are not exchangeable independent local systems. Extrapolating its negative average slope into the high-$\Sigma_\star$ MAUVE regime would therefore beg the question.

## 6.5 Sharda et al. (2021): radial transport as missing physics

Sharda et al. construct an analytical metallicity-distribution model for galactic discs that includes cosmological accretion, metal production, metal-enriched outflow, radial advection, and turbulent diffusion. It explains mass-dependent metallicity gradients and provides a more complete spatial conservation equation than independent annuli. ([Sharda et al. 2021](https://arxiv.org/abs/2102.09733))<!--ref:sharda2021--><!--anchor:section:2-->

It is not a direct rFMR formula, but it matters because radial transport can decorrelate local SFR and local enrichment. A one-patch positive relation can be weakened or reversed if enriched material is advected or mixed over a scale larger than the SFR aperture. Thus Carton-like local equilibrium and Wang–Lilly local forcing bracket two limits; Sharda supplies the spatial coupling omitted by both.

## 6.6 Huang et al. (2026): a mass-dependent resolved sign inversion

The MAUVE–MUSE analysis reports a positive $Z_{\mathrm{g}}$–$\Sigma_{\mathrm{SFR}}$ residual at high $\Sigma_\star$ and a negative one at low $\Sigma_\star$, with an inversion around $\log\Sigma_\star\sim7.5$–8.0 in the adopted units and an important dependence on abundance indicator. It derives the surface-density regulator

$$
\frac{{\mathrm{d}}Z}{{\mathrm{d}}t}
=\frac{1}{\Sigma_{\mathrm{g}}}
\left[y(1-R)\Sigma_{\mathrm{SFR}}+(Z_\Phi-Z)\Sigma_\Phi\right],
$$

or

$$
\frac{{\mathrm{d}}Z}{{\mathrm{d}}t}=\frac{Z^\dagger-Z}{\tau_\Phi},\qquad
Z^\dagger=Z_\Phi+y(1-R)\frac{\Sigma_{\mathrm{SFR}}}{\Sigma_\Phi}.
$$

The perturbation is

$$
\delta Z^\dagger=(Z^\dagger-Z_\Phi)
\left(\delta\log\Sigma_{\mathrm{SFR}}-\delta\log\Sigma_\Phi\right),
$$

and the actual abundance is the exponential time-filtered response over $\tau_\Phi$. Equivalently, the sign is governed by the competition

$$
\delta Z\ \propto\ \delta\log\tau_\Phi-\delta\log\tau_{\mathrm{dep}}.
$$

Efficiency/depletion-time variations dominate the positive high-density branch; supply/dilution variations dominate the negative low-density branch. ([Huang et al. 2026](https://arxiv.org/abs/2605.31412))<!--ref:huang2026--><!--anchor:section:5.2-->

This model is directly aligned with Wang & Lilly (2021), but adds a resolved, mass-dependent observational target. Its main uncertainties are metallicity-calibrator dependence, DIG/shock contamination, gas-surface-density availability, the conversion from spatial ensembles to temporal forcing, and the fact that environmental processes in Virgo can change both supply and efficiency.

# 7. Positive correlations already seen or predicted elsewhere

A positive branch should not be introduced as though no precedent exists.

* Yates, Kauffmann & Guo (2012) find opposite SFR–metallicity signs at low and high stellar mass in SDSS and L-GALAXIES. Their massive, low-SFR, low-metallicity systems have exhausted gas after a burst and are subsequently diluted by weak accretion. Yates & Kauffmann (2014) develop this “slow dilution after quenching” explanation and identify correlated morphology, black-hole mass, gas fraction, and $Z_{\mathrm{g}}-Z_\star$ signatures. This is a semi-analytic evolutionary population mechanism, not the same as the short-timescale SFE driver. ([Yates et al. 2012](https://arxiv.org/abs/1107.3145))<!--ref:yates2012--><!--anchor:abstract--> ([Yates & Kauffmann 2014](https://arxiv.org/abs/1310.5151))<!--ref:yates2014--><!--anchor:abstract-->
* Forbes et al. (2014) explicitly note that galaxy-to-galaxy scatter in mass loading can generate a **positive** equilibrium SFR–metallicity slope, because both quantities fall as $\eta$ rises. Their worked stochastic model holds $\eta$ fixed and therefore isolates the negative inflow route. ([Forbes et al. 2014](https://www.mso.anu.edu.au/~krumholz/publications/2014/forbes14a.pdf))<!--ref:forbes2014--><!--anchor:page:14-->
* Wang & Lilly (2021) analytically predict a positive response to SFE forcing, and Huang et al. (2026) use the same driver competition for a resolved sign inversion. These are the closest precedents to the proposed modifications.

The mechanisms make different predictions. Slow dilution predicts positive SFR–$Z$ because the low-SFR tail is abnormally metal poor and often morphologically distinct. Efficiency forcing predicts the high-SFR tail to be metal enhanced at otherwise comparable gas supply. Mass-loading scatter predicts both SFR and equilibrium $Z$ to track $1/(1-R+\eta)$. They should not be collapsed into one “positive FMR” label.

# 8. Modifying Wang et al. (2026) to obtain a positive fixed-mass slope

## 8.1 Why the published ideal mapping is negative

In the inflow-driven limit, Wang et al. obtain

$$
Z_{\mathrm{g}}=\frac{y\epsilon}{1-R}\frac{M_\star}{{\mathrm{SFR}}}.
$$

If $M_\star$ and $\epsilon$ are fixed, $Z_{\mathrm{g}}\propto {\mathrm{SFR}}^{-1}$: the slope is necessarily negative. The gFMR says the same thing in state-variable form: high gas fraction means low enrichment per unit gas, and fixed efficiency maps high gas mass into high SFR.

## 8.2 Minimal residual-efficiency extension

Do **not** alter the metal-yield sign or introduce an ad hoc positive term. Retain the continuity equations and allow efficiency to vary between systems or in time:

$$
\ln\epsilon=\ln\bar\epsilon(M_\star,z)+\delta_\epsilon(t),
\qquad
\ln\Phi=\ln\bar\Phi(M_\star,z)+\delta_\Phi(t).
$$

At fixed $M_\star$ in the early inflow-driven solution,

$$
t=\left[\frac{2M_\star}{(1-R)\Phi\epsilon}\right]^{1/2},
$$

$$
{\mathrm{SFR}}=\left[\frac{2\epsilon\Phi M_\star}{1-R}\right]^{1/2},
\qquad
Z_{\mathrm{g}}=\frac{y}{2}
\left[\frac{2\epsilon M_\star}{(1-R)\Phi}\right]^{1/2}.
$$

Thus, at fixed inflow,

$$
Z_{\mathrm{g}}=\frac{y}{2}\frac{{\mathrm{SFR}}}{\Phi}
$$

in the yield convention of the FMR preprint: increasing efficiency moves both SFR and $Z_{\mathrm{g}}$ upward. In logarithmic residuals,

$$
\delta\ln {\mathrm{SFR}}=\tfrac12(\delta\ln\epsilon+\delta\ln\Phi),
$$

$$
\delta\ln Z_{\mathrm{g}}=\tfrac12(\delta\ln\epsilon-\delta\ln\Phi).
$$

For uncorrelated drivers,

$$
{\mathrm{Cov}}(\delta\ln Z_{\mathrm{g}},\delta\ln {\mathrm{SFR}})
=\frac14\left[{\mathrm{Var}}(\delta\ln\epsilon)
-{\mathrm{Var}}(\delta\ln\Phi)\right].
$$

The fixed-mass sign is positive precisely when efficiency variance exceeds inflow variance in this limit. With correlated drivers, the covariance terms cancel from this particular product but affect the marginal variances and evolutionary interpretation. This analytic result is stronger than merely fitting a negative or positive $\alpha$.

The same conclusion follows from Wang's projected equation: if $\epsilon\propto M_\star^a{\mathrm{SFR}}^b$, then $\alpha=(1-b)/(1+a)$. For $1+a>0$, a positive fixed-mass slope requires $b>1$. The early fixed-$M_\star$, fixed-$\Phi$ ensemble above has ${\mathrm{SFR}}\propto\epsilon^{1/2}$, hence $b=2$. Treating this as a causal residual model is preferable to inserting a freely fitted $b>1$ with no driver.

## 8.3 Finite-time version beyond the inflow limit

Integrate the original Wang equations with $\epsilon(t)$ and $\Phi(t)$ rather than substituting a static SFR mapping. Use two stochastic processes with amplitudes $(\sigma_\epsilon,\sigma_\Phi)$ and coherence times $(t_\epsilon,t_\Phi)$. The abundance equation is a low-pass filter with response time $\tau_\Phi=M_{\mathrm{g}}/\Phi$. A positive relation requires not only sufficient efficiency variance but also coherence long enough for enrichment to respond. Very rapid $\epsilon$ flickering can raise an instantaneous SFR tracer without measurably changing $Z$.

This is the recommended Wang modification because it:

1. leaves the published conservation laws intact;
2. reduces to the published negative FMR when $\sigma_\epsilon=0$;
3. reduces to the Wang–Lilly positive case when $\sigma_\Phi=0$; and
4. predicts a continuous sign transition rather than imposing two unrelated formulae.

## 8.4 Other Wang modifications, ranked

1. **Variable efficiency — preferred.** Direct precedent, observable through depletion time, and mathematically sufficient.
2. **Correlated enriched/recycled inflow.** Add $Z_{\mathrm{in}}(t)$; positive if enriched return rises with SFR strongly enough to beat dilution. Appropriate in fountains or bar/spiral recycling, but requires metal-tagging or a recycling model.
3. **Variable mass loading.** Near equilibrium, lower $\eta$ raises both SFR and $Z$. Plausible across potential depth or feedback state, but difficult to vary independently at fixed mass and local surface density.
4. **Purely phenomenological $b>1$.** Algebraically works but is least explanatory and risks circularity because efficiency is defined using SFR.

# 9. Modifying Forbes et al. (2014)

## 9.1 Route already identified by Forbes: mass-loading scatter

At equilibrium, for well-mixed outflow and fixed inflow,

$$
{\mathrm{SFR}}_{\mathrm{eq}}=\frac{\Phi}{f_R+\eta},
$$

$$
Z_{\mathrm{eq}}-Z_{\mathrm{IGM}}=\frac{y f_R}{f_R+\eta}
=\frac{y f_R}{\Phi}{\mathrm{SFR}}_{\mathrm{eq}}.
$$

Scatter in $\eta$ therefore moves equilibrium systems along an exactly positive line in this ideal limit. Forbes et al. discuss this possibility but constrain their main stochastic calculation to fixed $\eta$ in order to infer accretion statistics. ([Forbes et al. 2014](https://arxiv.org/abs/1311.1509))<!--ref:forbes2014--><!--anchor:section:5.1-->

This route is legitimate but describes cross-system equilibrium variation. It is not the same temporal loop as a burst of efficiency. It predicts that independent outflow diagnostics should correlate with both residuals.

## 9.2 Add a second stochastic efficiency driver

The original dimensionless recurrence assumes a fixed loss time. Replace the constant depletion law by

$$
{\mathrm{SFR}}(t)=\epsilon(t)M_{\mathrm{g}}(t),
$$

$$
\dot M_{\mathrm{g}}=\Phi(t)-[f_R+\eta]\epsilon(t)M_{\mathrm{g}},
$$

and retain the full physical-time metal equation. Let

$$
\ln\Phi(t)=\mu_\Phi+\sigma_\Phi x_\Phi(t),\qquad
\ln\epsilon(t)=\mu_\epsilon+\sigma_\epsilon x_\epsilon(t),
$$

where the two stochastic processes can have different coherence times and a fitted cross-correlation. The original Forbes model is recovered exactly for $\sigma_\epsilon=0$.

A static galaxy-to-galaxy scatter in depletion time does not change equilibrium SFR or equilibrium metallicity at fixed $\Phi$ and $\eta$; it only changes gas mass and the response time. A **time-variable or phase-selected** efficiency is required for the Wang–Lilly positive transient. This is why simply drawing a new $t_{\mathrm{dep}}$ once per galaxy and reusing the Forbes equilibrium formula would not test the intended physics.

## 9.3 Expected sign and identifiability

With independent inflow and efficiency fluctuations, the observed covariance is a sum of filtered components:

$$
{\mathrm{Cov}}(\Delta\log Z,\Delta\log {\mathrm{SFR}})
=C_\epsilon(\sigma_\epsilon,t_\epsilon)
-C_\Phi(\sigma_\Phi,t_\Phi)
+C_\eta+C_{Z_{\mathrm{in}}}+C_{\mathrm{cross}}.
$$

$C_\epsilon$ and $C_\Phi$ are non-negative response amplitudes before their displayed signs. A positive slope only establishes that the sum of enrichment-linked terms exceeds dilution-linked terms. Without gas mass, depletion time, inflow/recycling, or outflow constraints, amplitudes and coherence times are degenerate.

# 10. Reproducible two-driver experiment

![Fixed-mass gas-regulator experiment. The first three panels show residual tracks for inflow-only, efficiency-only, and mixed forcing. The final panel gives the correlation as a function of the rms efficiency-to-inflow fluctuation ratio; shading is the 16th–84th percentile over 12 realisations.](20260813_FMR_driver_competition_toy_model.png){width=100%}

## 10.1 Setup

The accompanying Python calculation integrates the common gas and metal equations at fixed stellar mass. It uses $R=0.4$, $\eta=1$, $\epsilon_0=1$ in arbitrary inverse-time units, $\Phi_0=1$, $y=0.015$, and $Z_{\mathrm{in}}=0.002$. This gives $\tau_{\mathrm{eq}}=0.625$. Independent log-inflow and log-efficiency fluctuations are Ornstein–Uhlenbeck processes with coherence time $2\tau_{\mathrm{eq}}$. For the amplitude scan, $\sigma_{\ln\Phi}=0.28$ and $\sigma_{\ln\epsilon}/\sigma_{\ln\Phi}$ ranges from 0 to 2.5.

The experiment samples a regulator patch after burn-in and correlates residual $\log {\mathrm{SFR}}$ with residual $\log Z$. It does **not** evolve $M_\star$ or generate a cosmological galaxy population, so it tests the local fixed-mass response sign only.

## 10.2 Results

| Driver case | $\sigma_{\ln\Phi}$ | $\sigma_{\ln\epsilon}$ | Pearson $r$ | OLS $\Delta\log Z/\Delta\log {\mathrm{SFR}}$ |
|---|---:|---:|---:|---:|
| Inflow only (Forbes limit) | 0.28 | 0.00 | -0.385 | -0.122 |
| Efficiency only | 0.00 | 0.28 | +0.583 | +0.279 |
| Both, efficiency-dominated | 0.28 | 0.42 | +0.245 | +0.097 |

For these parameter choices, equal rms efficiency and inflow fluctuations produce a median correlation close to zero; the first sampled positive median occurs at a ratio of 1.125. The transition is broad across realisations. Its location changes with coherence time, inflow metallicity, equilibrium state, driver correlation, SFR averaging window, and response-time definition. The robust result is qualitative: **the unmodified continuity laws admit either sign once the dominant driver is allowed to vary.**

## 10.3 What this experiment does and does not establish

It establishes mathematical sufficiency and verifies the numerical implementation against the analytically expected signs. It does not show that MAUVE’s positive branch is caused by efficiency fluctuations. That requires testing the predicted gas, time-scale, environment, and ionisation signatures.

# 11. A testable MAUVE implementation

## 11.1 Minimum observational model

For each spatial element or resolution-matched region, model residuals after removing the chosen $\Sigma_\star$ trend:

$$
\Delta\log\Sigma_{\mathrm{SFR}}=A_\Phi u_\Phi+A_\epsilon u_\epsilon+e_{\mathrm{SFR}},
$$

$$
\Delta\log Z=-B_\Phi\,\mathcal{F}_\Phi[u_\Phi]
+B_\epsilon\,\mathcal{F}_\epsilon[u_\epsilon]+e_Z,
$$

where $\mathcal{F}$ are response-time filters. If gas data exist, add

$$
\Delta\log\Sigma_{\mathrm{g}}=G_\Phi u_\Phi-G_\epsilon\mathcal{F}_{\mathrm{g}}[u_\epsilon]+e_{\mathrm{g}}.
$$

The hierarchy can allow $(A,B,G)$ or the driver variances to change smoothly with $\Sigma_\star$, radius, arm/inter-arm status, or environmental class. The scientific parameter is the crossing where the efficiency-linked covariance exceeds the inflow-linked covariance—not an arbitrary broken-line sign fit.

## 11.2 Discriminating predictions

| Observation at fixed $\Sigma_\star$ | Inflow/dilution driver | Efficiency driver | Slow post-burst dilution | Mass-loading scatter |
|---|---|---|---|---|
| $\Delta Z$ versus $\Delta\Sigma_{\mathrm{SFR}}$ | Negative | Positive | Positive mainly through low-SFR, low-$Z$ tail | Positive near equilibrium |
| $\Delta Z$ versus $\Delta\Sigma_{\mathrm{g}}$ | Negative | Often remains negative or weak | Depends on depleted state | Indirect |
| Depletion time | High SFR mostly from more gas | High SFR at shorter $\tau_{\mathrm{dep}}$ | Long/quenched in diluted systems | Not the primary axis |
| Morphology / stellar age | No unique requirement | May track arms, compression, dense gas | Old/bulge-dominated or post-burst | Tracks feedback/wind state |
| Outflow diagnostic | Consequence, not sole driver | May respond after burst | Historical | Directly correlated with residuals |
| $Z_{\mathrm{g}}-Z_\star$ | No unique sign | Modest short-time response | Especially low after late dilution | Depends on history |

## 11.3 Practical sequence

1. **Re-measure the sign with matched masks.** Hold $\Sigma_\star$ fixed flexibly; avoid a single global linear correction if the rMZR is curved. Use identical spaxels for all metallicity indicators.
2. **Stress-test abundance systematics.** Repeat with O3N2, N2, and any available multi-line or Bayesian calibration; control N/O sensitivity, ionisation parameter, DIG, shocks, and Balmer extinction. The reported MAUVE inversion is indicator dependent, so no physical fit should hide this branch.
3. **Add gas state variables.** Convolve CO and HI products to the MUSE/analysis resolution and compare $Z$ with $\Sigma_{\mathrm{g}}$ and $\tau_{\mathrm{dep}}=\Sigma_{\mathrm{g}}/\Sigma_{\mathrm{SFR}}$. Even incomplete CO coverage is more informative than substituting a deterministic Kennicutt–Schmidt law, which would erase the very efficiency scatter being tested.
4. **Exploit spatial context.** Compare arms/inter-arm regions, leading/trailing sides, bar/ring regions, projected ICM-wind direction, and truncation fronts. Compression can increase dense-gas fraction and efficiency; fresh accretion or radial transport can increase supply.
5. **Fit negative-only and two-driver models.** Use held-out galaxies, not random spaxels, for cross-validation. Spatially adjacent spaxels are not independent. Compare predictive performance and posterior driver amplitudes rather than only an information criterion on all pixels.
6. **Test response-time sensitivity.** H$\alpha$ SFR averages a different time window from FUV or resolved stellar-population SFR. The positive efficiency signal should weaken if the SFR tracer varies faster than enrichment/mixing can respond.
7. **Report galaxies individually and pooled.** A pooled sign can be produced by between-galaxy offsets. Estimate within-galaxy slopes, a hierarchical population distribution, and the fraction of galaxies genuinely requiring a positive component.

## 11.4 Falsification criteria

The efficiency interpretation would be weakened if:

* the positive branch disappears under metallicity diagnostics insensitive to N/O or ionisation parameter;
* positive $Z$–SFR residuals are entirely explained by between-galaxy zero-point offsets;
* high-SFR, high-$Z$ regions have higher gas mass but no shorter depletion time;
* the sign tracks DIG/shock classification rather than star-forming gas;
* a slow-dilution population model explains the residuals and their stellar-age/morphological signatures better; or
* the inferred efficiency fluctuation time is much shorter than the metal-mixing/enrichment response required to create the observed $Z$ amplitude.

# 12. Critical synthesis

## 12.1 What is robust across models

1. Gas fraction or $M_\star/M_{\mathrm{g}}$ is the most direct state variable connecting integrated enrichment to current dilution.
2. The SFR coordinate is a projection through efficiency, so depletion-time scatter is not nuisance scatter—it can change the sign.
3. Equilibrium sets mean relations; non-equilibrium response and parameter scatter set covariance.
4. A negative FMR is expected when stochastic metal-poor inflow dominates on a timescale comparable to the gas-loss/response time.
5. A positive branch has at least three physically distinct routes: efficiency forcing, mass-loading scatter, and slow post-burst dilution. Enriched recycling is a fourth plausible route.

## 12.2 What is not robust

* The numerical value of the FMR projection coefficient $\alpha$ is not universal. It depends on functional form, sample, SFR scale, and abundance calibration.
* The existence of one non-evolving global surface across all redshifts is debated and can be obscured by calibration evolution and selection.
* A resolved correlation is not automatically a temporal evolutionary track. Spatial mixing, radial gradients, and region age violate naive ergodicity.
* Instantaneous recycling and homogeneous mixing are weakest exactly where sub-kpc measurements are most interesting.
* Strong-line metallicity indicators share emission lines with SFR and ionisation diagnostics, creating non-physical covariance if errors and selection are not propagated jointly.

## 12.3 Ranking the models for the present question

1. **Best sign physics:** Wang & Lilly (2021), because it changes one driver at a time and derives both signs.
2. **Best fixed-mass stochastic foundation:** Forbes et al. (2014), after adding an efficiency process or using its own mass-loading-scatter argument.
3. **Best current gas-state hierarchy:** Wang et al. (2026), because its gFMR makes the gas-to-SFR projection explicit; preprint status remains a caveat.
4. **Best resolved observationally anchored formulation:** Huang et al. (2026), because it expresses the sign as a competition of local response times and targets the MAUVE inversion.
5. **Best annular metallicity-profile model:** Carton et al. (2015), especially when gas profiles are available, though not sufficient alone for residual sign inversion.
6. **Best empirical rFMR benchmark:** Baker et al. (2023), useful as a flexible null surface but not as causal evidence.

# 13. Conclusions

The requested positive correlation can be produced without violating chemical conservation. The most defensible modification to both Wang et al. (2026) and Forbes et al. (2014) is to promote star-formation efficiency from a fixed mapping to a second stochastic or coherent physical driver. In the Wang inflow limit, a fixed-mass ensemble has a positive covariance whenever $\operatorname{Var}(\ln\epsilon)>\operatorname{Var}(\ln\Phi)$ under the stated independent-driver assumptions. In the time-dependent Forbes extension, the exact threshold is filtered by the efficiency and inflow coherence times relative to the gas/metal response. Forbes’s alternative equilibrium route—scatter in mass loading—also produces a positive slope and must be tested rather than overlooked.

For MAUVE, the high-value result would not be merely another positive fitted slope. It would be evidence that the sign inversion is accompanied by the predicted transition from gas-supply-driven to efficiency-driven residuals, with consistent gas fractions, depletion times, spatial context, and response-time behaviour. That is a falsifiable physical interpretation and a natural bridge between the gFMR, global FMR, and rFMR literatures.

# References

Ascasibar, Y. et al. 2015, *MNRAS*, 448, 2126, “Understanding chemical evolution in resolved galaxies—I. The local star fraction–metallicity relation,” [arXiv:1406.6397](https://arxiv.org/abs/1406.6397), [doi:10.1093/mnras/stv098](https://doi.org/10.1093/mnras/stv098).

Baker, W. M. et al. 2023, *MNRAS*, 519, 1149, “The metallicity's fundamental dependence on both local and global galactic quantities,” [arXiv:2210.03755](https://arxiv.org/abs/2210.03755), [doi:10.1093/mnras/stac3594](https://doi.org/10.1093/mnras/stac3594).

Bothwell, M. S. et al. 2013, *MNRAS*, 433, 1425, “A fundamental relation between the metallicity, gas content and stellar mass of local galaxies,” [arXiv:1304.4940](https://arxiv.org/abs/1304.4940).

Bothwell, M. S. et al. 2016, *MNRAS*, 455, 1156, “The fundamental relation between metallicity and molecular gas content,” [arXiv:1507.01004](https://arxiv.org/abs/1507.01004).

Carton, D. et al. 2015, *MNRAS*, 451, 210, “Gas-phase metallicity profiles of the Bluedisk galaxies: is metallicity in a local star-formation regulated equilibrium?,” [arXiv:1505.02797](https://arxiv.org/abs/1505.02797), [doi:10.1093/mnras/stv967](https://doi.org/10.1093/mnras/stv967).

Curti, M. et al. 2020, *MNRAS*, 491, 944, “The mass–metallicity and the fundamental metallicity relation revisited,” [arXiv:1910.00597](https://arxiv.org/abs/1910.00597).

Davé, R., Finlator, K. & Oppenheimer, B. D. 2012, *MNRAS*, 421, 98, “An analytic model for the evolution of the stellar, gas and metal content of galaxies,” [arXiv:1108.0426](https://arxiv.org/abs/1108.0426).

Dayal, P., Ferrara, A. & Dunlop, J. S. 2013, *MNRAS*, 430, 2891, “The physics of the fundamental metallicity relation,” [arXiv:1202.4770](https://arxiv.org/abs/1202.4770), [doi:10.1093/mnras/stt083](https://doi.org/10.1093/mnras/stt083).

Feldmann, R. 2013, *MNRAS*, 433, 1910, “The equilibrium view on galaxy evolution,” [arXiv:1212.2973](https://arxiv.org/abs/1212.2973).

Forbes, J. C., Krumholz, M. R., Burkert, A. & Dekel, A. 2014, *MNRAS*, 443, 168, “On the origin of the fundamental metallicity relation and the scatter in galaxy scaling relations,” [arXiv:1311.1509](https://arxiv.org/abs/1311.1509), [doi:10.1093/mnras/stu1142](https://doi.org/10.1093/mnras/stu1142).

Harwit, M. & Brisbin, D. 2015, *ApJ*, 800, L26, “Origin of galaxy mass–metallicity–star-formation relation,” [arXiv:1412.2436](https://arxiv.org/abs/1412.2436).

Huang, R. et al. 2026, *MNRAS*, 549, 1019, “MAUVE–MUSE: When metallicity follows or fights star formation—a mass-dependent inversion in Virgo galaxies,” [arXiv:2605.31412](https://arxiv.org/abs/2605.31412), [doi:10.1093/mnras/stag1019](https://doi.org/10.1093/mnras/stag1019).

Hunt, L. K. et al. 2016a, *MNRAS*, “Coevolution of metallicity and star formation in galaxies to $z=3.7$—I. A fundamental plane,” [arXiv:1608.05417](https://arxiv.org/abs/1608.05417), [doi:10.1093/mnras/stw1993](https://doi.org/10.1093/mnras/stw1993).

Hunt, L. K. et al. 2016b, *MNRAS*, 463, 2020, “Coevolution of metallicity and star formation in galaxies to $z=3.7$—II. A theoretical model,” [arXiv:1608.05418](https://arxiv.org/abs/1608.05418), [doi:10.1093/mnras/stw2091](https://doi.org/10.1093/mnras/stw2091).

Lara-López, M. A. et al. 2010, *A&A*, 521, L53, “A fundamental plane for field star-forming galaxies,” [arXiv:1005.0509](https://arxiv.org/abs/1005.0509).

Lilly, S. J. et al. 2013, *ApJ*, 772, 119, “Gas regulation of galaxies: the evolution of the cosmic specific star formation rate, the metallicity–mass–star-formation rate relation, and the stellar content of haloes,” [arXiv:1303.5059](https://arxiv.org/abs/1303.5059), [doi:10.1088/0004-637X/772/2/119](https://doi.org/10.1088/0004-637X/772/2/119).

Lin, Y.-T. & Zu, Y. 2023, *MNRAS*, 521, 411, “Constraints on galactic outflows from the metallicity–stellar mass–SFR relation of EAGLE simulation and SDSS galaxies,” [arXiv:2212.01402](https://arxiv.org/abs/2212.01402), [doi:10.1093/mnras/stad502](https://doi.org/10.1093/mnras/stad502).

Magrini, L. et al. 2012, *MNRAS*, “Scaling relations of metallicity, stellar mass and star formation rate in metal-poor starbursts—II. Theoretical models,” [doi:10.1111/j.1365-2966.2012.22055.x](https://doi.org/10.1111/j.1365-2966.2012.22055.x).

Mannucci, F. et al. 2010, *MNRAS*, 408, 2115, “A fundamental relation between mass, star formation rate and metallicity in local and high-redshift galaxies,” [arXiv:1005.0006](https://arxiv.org/abs/1005.0006), [doi:10.1111/j.1365-2966.2010.17291.x](https://doi.org/10.1111/j.1365-2966.2010.17291.x).

Peng, Y.-j. & Maiolino, R. 2014, *MNRAS*, 443, 3643, “The gas regulator model of galaxies: the role of inflows, outflows and star formation,” [arXiv:1402.5964](https://arxiv.org/abs/1402.5964).

Pipino, A., Lilly, S. J. & Carollo, C. M. 2014, *MNRAS*, 441, 1444, “On the metallicity of galaxies in the gas-regulator model,” [arXiv:1403.6146](https://arxiv.org/abs/1403.6146).

Sánchez Almeida, J. & Sánchez-Menguiano, L. 2019, *ApJL*, “The fundamental metallicity relation emerges from the local anti-correlation between star formation rate and gas-phase metallicity existing in disk galaxies,” [arXiv:1905.05826](https://arxiv.org/abs/1905.05826), [doi:10.3847/2041-8213/ab218d](https://doi.org/10.3847/2041-8213/ab218d).

Sharda, P. et al. 2021, *MNRAS*, 502, 5935, “On the origin of the mass–metallicity gradient relation in the local Universe,” [arXiv:2102.09733](https://arxiv.org/abs/2102.09733).

Wang, E. & Lilly, S. J. 2021, *ApJ*, 910, 137, “Gas-phase metallicity as a diagnostic of the drivers of star formation on different spatial scales,” [arXiv:2009.01935](https://arxiv.org/abs/2009.01935), [doi:10.3847/1538-4357/abe413](https://doi.org/10.3847/1538-4357/abe413).

Wang, K. 2026, *MNRAS*, 545, staf2113, “The origin of the galaxy size–stellar metallicity relation—I. A semi-analytical perspective,” [arXiv:2510.02573](https://arxiv.org/abs/2510.02573), [doi:10.1093/mnras/staf2113](https://doi.org/10.1093/mnras/staf2113).

Wang, K. et al. 2026, arXiv preprint v1, “Inflow-driven galaxy evolution—I. Revealing the physics of the fundamental metallicity relation,” [arXiv:2608.04784](https://arxiv.org/abs/2608.04784).

Yates, R. M., Kauffmann, G. & Guo, Q. 2012, *MNRAS*, 422, 215, “The relation between metallicity, stellar mass and star formation in galaxies,” [arXiv:1107.3145](https://arxiv.org/abs/1107.3145).

Yates, R. M. & Kauffmann, G. 2014, *MNRAS*, “Dilution in elliptical galaxies: implications for the relation between metallicity, stellar mass and star formation rate,” [arXiv:1310.5151](https://arxiv.org/abs/1310.5151).

Zahid, H. J. et al. 2012, *ApJ*, 757, 54, “A census of oxygen in star-forming galaxies: an empirical model linking metallicities, star formation rates, and outflows,” [arXiv:1207.5509](https://arxiv.org/abs/1207.5509).

# AI and reproducibility disclosure

OpenAI Codex assisted with the online literature search, equation extraction, synthesis, numerical experiment, figure generation, and document production. The review prioritised primary papers and marks Wang et al. (2026) as a new preprint. The numerical experiment is deterministic given the recorded random seeds and parameters; its source script and CSV scan are retained with the working materials. Scientific interpretation, paper selection, and any use in a manuscript should receive domain-expert review, particularly for yield conventions, metallicity calibration, and the equivalence—or non-equivalence—of spatial and temporal residuals.
