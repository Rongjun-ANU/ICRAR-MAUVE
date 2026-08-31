---
title: "Ram-pressure stripping and local star formation"
subtitle: "Theory, step-by-step analytical models, observable predictions, and a directional gas regulator"
author: "Research report prepared for Rongjun"
date: "30 August 2026"
lang: en
---

# 1. Main answer and recommended approach

**Yes: there are analytical models of pressure-triggered star formation, pressure-confined cloud collapse, and gravitational instability in compressed rotating discs.** There are also models that combine pressure-dependent star formation with gas removal. The proposed research is therefore physically well motivated. The most useful extension is to connect these established ingredients to **resolved, directional, time-dependent gas conservation**, with explicit predictions for quantities that can be observed.

The key distinction is between three changes:

1. **More gas in an aperture:** in-plane accumulation, inward transport, fallback, or conversion from atomic to molecular gas.
2. **Faster star formation per unit relevant gas:** a shorter collapse time, a larger bound/dense fraction, or a change in efficiency per free-fall time.
3. **Less gas available to form stars:** stripping, destruction/mixing of cold material, reduced supply, and consumption.

A leading-side SF excess can arise from the first mechanism without any increase in molecular-gas star formation efficiency. A trailing disc can remain active, be enhanced after transport or a delay, or be suppressed. A detached tail is a third physical environment; it is not interchangeable with the trailing half of a rotating disc.

Three particularly relevant starting points are the pressure-dependent cloud-population model of [Fujita & Nagashima (1999)][fujita], the jellyfish-galaxy model of [Safarzadeh & Loeb (2019)][safarzadeh], and the pressure-confined rotating-disc stability analysis of [Kim et al. (2012)][kim]. For a framework closer to resolved Virgo observations, see [Nehlig et al. (2016)][nehlig] and [Lizee et al. (2021)][lizee]. The latter includes both tidal forcing and RPS, so it is not a clean calibration of ram pressure alone.

**Recommended model:** start with a local gas regulator in which pressure changes the star formation rate per unit gas, and stripping and transport change the gas reservoir. Use Jeans/Bonnor-Ebert and Toomre calculations as physical checks on the assumed response. Do not use a single modified $Q$ as both a collapse criterion and an SFR law.

The central conditional relation derived below is

$$
\boxed{\mathcal E=g\,f\,e\,\frac{\sqrt{d\,p}}{s}.}
\tag{1.1}
$$

Here $\mathcal E$ is the local SFR surface-density ratio relative to a matched reference; $g$ is the cold-gas surface-density ratio; $p$ is the **internal midplane-pressure ratio**, not automatically the ram-pressure ratio; $s$ is the vertical support-speed ratio; and $f,e,d$ describe changes in star-forming gas fraction, efficiency per free-fall time, and density contrast. If those last three factors are fixed, enhancement requires $g\sqrt p>s$.

This relation is an explicitly stated closure, not a validated universal RPS prescription. Its value is that it separates competing effects, supplies a null prediction, and exposes the degeneracies in any pressure inference from SF.

The report provides six original figures, exact solutions for a simplified regulator, independent numerical checks, and a practical observing/inference strategy. **No MAUVE measurements have been fitted or analysed here.** All numerical galaxy examples are deliberately illustrative. A literature-supported mechanism, an algebraically correct model, and an observationally validated model are different levels of evidence.

# 2. Define the question before choosing a model

## 2.1 Spatial regions and the meaning of enhancement

Let the **windward/leading disc** be the part facing the incoming ICM flow, and the **leeward/trailing disc** the opposite region that remains in, or close to, the rotating galaxy. A **detached tail** consists of extraplanar/downstream material whose dynamics need not follow disc rotation. The projected tail direction is a useful directional indicator, but it can reflect past motion, ICM flows, bending, and projection rather than the instantaneous three-dimensional wind vector.

Use a fixed, specified physical aperture or common deprojected area. Define

$$
\mathcal E(R,\phi)=
\frac{\Sigma_{\rm SFR}(R,\phi)}{\Sigma_{{\rm SFR},0}(R,\phi)},
\qquad
t_{\rm dep,mol}=\frac{\Sigma_{\rm mol}}{\Sigma_{\rm SFR}}.
\tag{2.1}
$$

The reference should match relevant conditions such as stellar surface density, galactocentric radius, galaxy mass, resolution, and selection. A reference at fixed stellar surface density is appropriate for a resolved SF main-sequence residual. A reference at fixed molecular-gas surface density asks a different question: whether molecular-gas efficiency changed.

For the molecular reservoir,

$$
\Delta\log_{10}\Sigma_{\rm SFR}
=\Delta\log_{10}\Sigma_{\rm mol}
-\Delta\log_{10}t_{\rm dep,mol}.
\tag{2.2}
$$

This identity is a useful first diagnostic. A positive SF residual alone cannot separate the two terms on its right-hand side.

Likewise, a leading/trailing contrast,

$$
\Delta_{\rm LT}=\log_{10}
\frac{\Sigma_{{\rm SFR},L}}{\Sigma_{{\rm SFR},T}},
\tag{2.3}
$$

does not establish absolute enhancement. If the two regions have identical reference SFRs, $\mathcal E_L=0.5$ and $\mathcal E_T=0.2$ give $\Delta_{\rm LT}=0.398$ dex even though both are suppressed. With unequal references, compare their residuals as well as their raw SFRs.

Compression of a moving parcel and a rise in fixed-aperture surface density are also different. A parcel can become denser because its thickness decreases while its projected surface density stays constant. Conversely, an aperture can accumulate molecular material without changing individual cloud properties. The model below specifies where each effect enters.

## 2.2 Why severe stripping and local enhancement can coexist

The diffuse outer disc and dense inner clouds have different restoring forces, columns, cooling rates, and exposure times. Substantial global gas loss can coexist with a short-lived burst in a retained region. Conversely, a high instantaneous external pressure can coincide with low SF because the relevant cold gas has already been removed. The physically meaningful question is therefore not simply whether pressure is high, but **which gas remains, for how long, and with what density and internal support**.

Neither an instantaneous SF map nor a morphological stripping label uniquely supplies the pressure history. In particular, a present-day leading-side excess need not occur at the epoch of maximum ram pressure.

# 3. What the theoretical and observational literature actually establishes

## 3.1 Analytical and semi-analytical foundations

**Pressure-dependent cloud populations.** [Fujita & Nagashima (1999)][fujita] evolve molecular-cloud mass bins with pressure-dependent cloud lifetimes and efficiencies, coupled to the supply of gas. Their adopted cluster models show a transient, modest global increase followed by decline when stripping removes the replenishing reservoir. Their result is not a universal upper limit on local enhancement, and their selected calculations retain the molecular clouds. This is an early, directly relevant precedent for a pressure-modified gas regulator.

**An integrated jellyfish model.** [Safarzadeh & Loeb (2019)][safarzadeh] combine a pressure-dependent cloud efficiency with an assumed cold-gas loss history to model elevated galaxy-integrated SF. The cloud-efficiency ingredient follows [Elmegreen & Efremov (1997)][elmegreen]. This is useful motivation, but neither the gas-loss law nor an azimuthally resolved windward/leeward response is derived from first principles. An efficiency without an associated timescale is not by itself an SFR.

**Resolved multiphase-disc models.** [Vollmer & Leroy (2011)][vollmer] provide a turbulent, clumpy gas-disc framework with accretion, feedback, and cloud timescales. [Nehlig et al. (2016)][nehlig] apply environmental compression to the multiphase ISM and find that compressed regions need not have greater molecular SFE. [Lizee et al. (2021)][lizee] model NGC4654 using observations and a turbulent-disc description; their compressed region combines changed gas content, stability, and turbulent support. CO conversion factors and the galaxy's tidal interaction matter to its interpretation.

**Pressure-confined gravitational instability.** [Kim et al. (2012)][kim] analyse rotating, vertically stratified, pressure-bounded gas slabs and construct a generalized stability parameter. Free boundaries introduce surface-distortion modes, so the relevant mode speed is not simply the observed gas linewidth. Strong confinement can favour instability without establishing the dense collapse required for SF; [Lubow & Pringle (1993)][lubow] make the related caution particularly clear.

**Cloud collapse.** [Bonnor (1956)][bonnor] supplies the pressure-confined isothermal-sphere instability underlying the Bonnor-Ebert threshold. This is usually a more direct starting point for a cloud-scale pressure-triggering argument than inserting a pressure term into a galaxy-scale $Q$ by inspection.

**Gas supply and loss.** The regulator approach of [Lilly et al. (2013)][lilly] provides the conservation structure. The RPS threshold of [Gunn & Gott (1972)][gunn] and the finite-pulse treatment of [Koppen et al. (2018)][koppen] constrain removal, but are not prescriptions for the local star formation rate or the instantaneous stripping rate.

Two related papers need careful scope labels. [Whitworth (2016)][whitworth], despite the title *A ram-pressure threshold for star formation*, studies converging molecular flows and cooling/fragmentation; its threshold should not be transplanted as a universal ICM-stripping criterion. [Pasetto et al. (2015)][pasetto] analyse environmental fluid instabilities for idealized systems including dwarfs and clusters; that geometry is not a resolved rotating spiral-disc SF model.

## 3.2 What the resolved observations require a model to explain

[Roberts et al. (2022)][roberts] find leading-side SF enhancement in a LoTSS/MaNGA jellyfish sample. In their detailed IC3949 example, the leading region has more molecular gas and a depletion time consistent with the galaxy median. This is an especially relevant demonstration that enhanced local SF need not imply enhanced molecular SFE.

[Vulcani et al. (2020)][vulcani] find elevated resolved SF at fixed stellar surface density in a stripping-galaxy sample. That comparison establishes an SF excess relative to the stellar distribution, but does not alone distinguish gas accumulation from faster conversion of molecular gas into stars.

[Brown et al. (2023)][brown] explicitly separate molecular content and efficiency in Virgo: early stripping can include locally enhanced molecular content, while more environmentally affected systems can show deficiencies in both gas and SF efficiency. The compressed regions studied by [Nehlig et al. (2016)][nehlig] also warn against assuming pressure and SFE rise together.

There are population-level counterexamples. [Cakir et al. (2026)][cakir] report outer suppression without a central SF excess relative to controls in their selected SAMI cluster populations. These results need not contradict a temporary enhancement in selected jellyfish systems: sample selection, stripping stage, aperture, reference population, and tracer timescale differ. A model must allow enhancement, no response, and suppression.

## 3.3 Numerical evidence: useful tests, not analytical proofs

Simulations support several channels, but their SF prescriptions and apertures matter:

- [Troncoso-Iribarren et al. (2020)][eagle] find an asymmetric response in EAGLE. Its adopted SF recipe already depends on pressure, so an SF-pressure correlation is not an independent validation of pressure-triggered cloud collapse. A trailing hemisphere is not synonymous with a detached tail.
- [Zhu et al. (2024)][zhu] demonstrate central enhancement through gas transport and increased dense-gas content, without requiring a changed SF law at fixed gas surface density. Their time-averaged disc halves can have similar SFRs; this is a useful counterexample to a mandatory leading-only response.
- [Bekki (2014)][bekki] finds outcomes that depend on galaxy and wind properties, including enhancement and quenching. [Lee et al. (2020)][lee] likewise find different outcomes for moderate and strong stripping in their adopted models. Their amplitudes are not universal calibration constants.
- [Tonnesen & Bryan (2009)][tonnesen2009] study multiphase gas stripping **without implementing star formation**. Their dense-gas results should not be reported as directly simulated SFRs. [Tonnesen & Bryan (2012)][tonnesen2012] do model SF and distinguish tail gas delivery by ram pressure from confinement associated with thermal ICM pressure.
- [Kapferer et al. (2009)][kapferer] obtain strong SF in stripped wakes in some calculations. Differences from other tail models emphasize sensitivity to gas treatment, resolution, feedback, and the definition of system versus disc SF.
- [Goller et al. (2023)][goller] find no population-wide SF enhancement for TNG50 jellyfish galaxies relative to their controls, despite individual time-dependent events. An individual burst and a population excess are different claims.

These findings motivate a **conditional** model: pressure, transport, survival, support, and time history must be allowed to vary. They do not justify fixing the sign of the trailing-side response in advance.

# 4. The external forcing: pressure, geometry, and removal

## 4.1 Start with the relative velocity and momentum flux

The physical wind is relative to the local ICM:

$$
\boldsymbol v_w=\boldsymbol v_{\rm gal}-\boldsymbol v_{\rm ICM},
\qquad P_{\rm ram}=\rho_{\rm ICM}|\boldsymbol v_w|^2.
\tag{4.1}
$$

There is no factor of $1/2$ in this momentum-flux convention. For an ideal fluid, the incident momentum-flux tensor is

$$
T_{ij}=P_{\rm th}\delta_{ij}+\rho_{\rm ICM}v_{w,i}v_{w,j}.
\tag{4.2}
$$

The normal stress on a surface with unit normal $\boldsymbol n$ is therefore

$$
P_n=P_{\rm th}+\rho_{\rm ICM}(\boldsymbol v_w\cdot\boldsymbol n)^2.
\tag{4.3}
$$

Thermal pressure is isotropic; it does not acquire a single cosine factor. If $c_s$ is the adiabatic ICM sound speed, $P_{\rm th}=\rho c_s^2/\gamma$. The stress incident on a surface is not automatically the static pressure transmitted to a molecular cloud. Shocks, deflection, shielding, magnetic stresses, and porosity mediate that transmission. Consequently, the phenomenological pressure expression used in an integrated model such as [Safarzadeh & Loeb (2019)][safarzadeh] should not be treated as an exact projection law.

For scale, assuming mass per electron $\mu_e m_p$,

$$
\frac{P_{\rm ram}}{k_B}=1.417\times10^5
\left(\frac{\mu_e}{1.17}\right)
\left(\frac{n_e}{10^{-3}\,{\rm cm}^{-3}}\right)
\left(\frac{v_w}{1000\,{\rm km\,s}^{-1}}\right)^2
\,{\rm K\,cm}^{-3}.
\tag{4.4}
$$

Here $n_e$ is the **ICM electron density**, not a density inferred from an H II-region line ratio. The numerical coefficient follows from the stated composition assumption and physical constants.

## 4.2 Confinement and net acceleration are different

For a column subject to bottom and top boundary stresses $P_-$ and $P_+$, its centroid approximately obeys

$$
\Sigma_g\ddot z_c=P_--P_+
-\int\rho\frac{\partial\Phi}{\partial z}\,dz.
\tag{4.5}
$$

The pressure **difference** accelerates the column. Equal confining pressures can compress it without producing that net acceleration. The equilibrium relation $P_{\rm mid}=P_{\rm ext}+W$ derived in Section 5 assumes a supported layer; it is not an exact identity for a rapidly accelerating, one-sided wind interaction.

For a minimal directional closure, introduce an effective confining pressure. One possible illustrative form is

$$
\begin{aligned}
P_{\rm conf}(R,\phi,t)=P_{\rm bg}(R,t)+P_{\rm ram}(t)
\big[&\chi_\perp\cos^2 i\\
&+\chi_\parallel\sin^2 i\,[\max(\cos\phi,0)]^2\big].
\end{aligned}
\tag{4.6}
$$

$i$ is the angle between the wind and the disc normal; $\phi=0$ is upstream in the disc plane. $P_{\rm bg}$ is the effective background boundary pressure. The dimensionless $\chi$ values summarize coupling and geometry and must be constrained or marginalized over. They are not known drag coefficients. Equation (4.6) is a chosen angular interpolation, not a solution of the wind/ISM boundary problem. A leeward pressure component can be added when justified.

For an exactly face-on wind in an initially axisymmetric disc, this closure has no preferred azimuth. For an inclined wind, pressure on the rim and cloud surfaces cannot be represented solely by projecting onto the large-scale disc normal. The two-region numerical example later uses independent effective couplings rather than claiming that equation (4.6) determines them.

## 4.3 A stripping condition does not determine a stripping rate

For slow, approximately face-on forcing, a useful diagnostic is

$$
\Pi=\frac{P_{{\rm ram},n}}{\Sigma_g g_{z,\max}},
\qquad
\Sigma_g g_{z,\max}\simeq2\pi G\Sigma_*\Sigma_g
\quad\hbox{in a thin stellar-sheet approximation}.
\tag{4.7}
$$

The last expression is the familiar [Gunn & Gott (1972)][gunn] estimate. It has geometric and timescale restrictions. It should not be equated to the internal midplane weight, nor interpreted as a universal criterion for parallel ablation or every gas phase.

For a short pulse, the delivered momentum per unit area,

$$
J=\int P_{{\rm ram},n}(t)\,dt,
\qquad \Delta v\sim J/\Sigma_g,
\tag{4.8}
$$

becomes important. [Koppen et al. (2018)][koppen] distinguish these slow- and fast-forcing regimes. Whether a parcel escapes still depends on its potential and initial velocity; arbitrary wind directions can exert torque. A scalar pressure criterion assuming conserved azimuthal angular momentum is therefore not generally valid.

The regulator will use $\dot\Sigma_{\rm strip}=k(t)\Sigma_g$ as a transparent **loss closure**. Neither equation (4.7) nor equation (4.8) fixes $k(t)$. A law such as $k\propto[\Pi-1]_+$ would require an additional removal timescale and calibration; it also omits stripping mechanisms below or outside that idealized threshold.

# 5. From confinement to a local star formation law

## 5.1 Derive the vertical weight

Consider a symmetric, approximately hydrostatic gas layer. For $z>0$ let $g_z$ denote the magnitude of the downward gravitational acceleration. Then

$$
\frac{dP}{dz}=-\rho g_z,
\qquad
P_{\rm mid}-P_{\rm conf}=\int_0^{z_b}\rho g_z\,dz\equiv W.
\tag{5.1}
$$

For a purely self-gravitating plane-parallel layer,

$$
\frac{dg_z}{dz}=4\pi G\rho,
\qquad g_z(z_b)=2\pi G\Sigma_g.
\tag{5.2}
$$

Changing the integration variable from height to acceleration gives

$$
W_g=\frac{1}{4\pi G}\int_0^{2\pi G\Sigma_g}g_z\,dg_z
=\frac{\pi G}{2}\Sigma_g^2.
\tag{5.3}
$$

This self-gravity result is exact for the stated geometry, including a pressure-truncated layer. A commonly used approximation when stellar gravity also matters is

$$
W\simeq\frac{\pi G}{2}\Sigma_g
\left(\Sigma_g+\frac{c_z}{\sigma_{*,z}}\Sigma_*\right),
\tag{5.4}
$$

where $c_z$ is the effective vertical gas support speed and $\sigma_{*,z}$ the stellar vertical dispersion. This approximation assumes compatible vertical distributions; it is not a universal relation for a disturbed gas layer.

Alternatively, in a locally uniform external stellar-plus-dark-matter density $\rho_{\rm sd}$, the harmonic vertical acceleration is $g_z=4\pi G\rho_{\rm sd}z$. An isothermal gas layer is Gaussian, with $\Sigma_g=\sqrt{2\pi}\rho_0 c_z/\sqrt{4\pi G\rho_{\rm sd}}$. Hence its externally generated weight is approximately $\Sigma_g c_z\sqrt{2G\rho_{\rm sd}}$. A useful interpolation is

$$
W\simeq\frac{\pi G}{2}\Sigma_g^2+
\Sigma_g c_z\sqrt{2G\rho_{\rm sd}}.
\tag{5.5}
$$

Equations (5.4) and (5.5) are alternative approximations; their stellar terms must not be added together. Feedback-regulated vertical-equilibrium models and molecular-ISM observations provide context for using such weight estimates, but do not establish equilibrium during a rapid RPS event. See [Ostriker, McKee & Leroy (2010)][ostriker] and [Sun et al. (2020)][sun].

## 5.2 Relate pressure, density, and free-fall time

If thermal and turbulent support can be represented by $P=\rho c_z^2$, equations (5.1)-(5.4) imply

$$
P_{\rm mid}=W+P_{\rm conf},
\qquad
\rho_{\rm mid}=\frac{P_{\rm mid}}{c_z^2},
\qquad
H_{\rm eff}=\frac{\Sigma_g}{2\rho_{\rm mid}}.
\tag{5.6}
$$

Do not silently include magnetic or cosmic-ray pressure in $P_{\rm mid}$ and still interpret a measured turbulent linewidth as all of $c_z$. Those support components need their own closure. Also distinguish vertical random support from coherent streaming and wind-driven bulk acceleration.

The characteristic density of star-forming material is generally greater than the mean density of a resolution element. Define it explicitly by

$$
\rho_{\rm sf}=\xi\rho_{\rm mid},
\qquad
t_{\rm ff,sf}=\left(\frac{3\pi}{32G\rho_{\rm sf}}\right)^{1/2}.
\tag{5.7}
$$

$\xi$ is an effective density contrast, not a measured universal clumping factor. A broad density distribution can make a one-density free-fall time inadequate; then replace the following closure by a mass-weighted density integral.

Let $f_{\rm sf}$ be the fraction of total cold-gas mass assigned to the star-forming reservoir. Adopt

$$
\Sigma_{\rm SFR}=\epsilon_{\rm ff}\frac{f_{\rm sf}\Sigma_g}{t_{\rm ff,sf}}.
\tag{5.8}
$$

The choice of reservoir is part of the definition. If it is molecular gas, $f_{\rm sf}=f_{\rm mol}$ and $\epsilon_{\rm ff}$ must be interpreted consistently. If it is only bound dense gas, its density and efficiency cannot be borrowed unchanged from a total-molecular-gas calibration. The molecular-fraction/free-fall form has a well-established precedent in [Krumholz, Dekel & McKee (2012)][krumholz]. Their [2013 erratum][krumholzerratum] changes the fitted normalization to $\epsilon_{\rm ff}=0.015$; it does not remove the phase-fraction factor. No such fitted normalization is imposed on the relative toy curves here.

Substitution yields a local rate per unit cold gas, $\nu\equiv\Sigma_{\rm SFR}/\Sigma_g$:

$$
\boxed{
\nu=\frac{\epsilon_{\rm ff}f_{\rm sf}}{c_z}
\left(\frac{32G\xi P_{\rm mid}}{3\pi}\right)^{1/2}.}
\tag{5.9}
$$

Its dimensions are inverse time: $(GP)^{1/2}/c_z$ has units of ${\rm s}^{-1}$. Pressure has not supplied extra gas. It has entered the assumed density and collapse timescale of the gas that remains.

## 5.3 Derive the enhancement condition

Compare a perturbed aperture to a matched reference, with

$$
\begin{gathered}
g=\Sigma_g/\Sigma_{g,0},\quad p=P_{\rm mid}/P_{{\rm mid},0},
\quad s=c_z/c_{z,0},\\
f=f_{\rm sf}/f_{{\rm sf},0},\quad
e=\epsilon_{\rm ff}/\epsilon_{{\rm ff},0},\quad d=\xi/\xi_0.
\end{gathered}
\tag{5.10}
$$

Taking the ratio of equation (5.8) gives

$$
\mathcal E=gfe\frac{\sqrt{dp}}{s},
\qquad
\mathcal E>1\ \Longleftrightarrow\
p>\frac{s^2}{g^2f^2e^2d}.
\tag{5.11}
$$

This result separates gas abundance from efficiency. With fixed $f,e,d$, doubling midplane pressure only increases $\nu$ by $\sqrt2$ if $c_z$ stays fixed. Losing half the gas would then leave $\mathcal E=0.707$, despite a shorter free-fall time. Conversely, gas accumulation can raise $\mathcal E$ even with unchanged molecular depletion time.

![Original illustrative regime diagram. The bold contour marks unchanged SFR. Gas loss and stronger vertical support shift a region toward suppression even when its internal pressure rises. Here the phase fraction, efficiency per free-fall time, and density contrast are held fixed; pressure and gas ratio are treated as independent diagnostic coordinates, not as a calibrated RPS trajectory.](assets/20260830_rps_sf/01_enhancement_regimes.png){width=100%}

In a time-dependent application, $W$, $P_{\rm mid}$, and sometimes $c_z$ change when gas is removed. They must be recomputed together. Holding $p$ fixed while arbitrarily changing $g$ is useful for the diagnostic diagram, but is not a self-consistent dynamical evolution unless the boundary pressure compensates.

## 5.4 Molecular conversion and turbulence need separate closures

A pressure-dependent molecular fraction is one possible ingredient. For example,

$$
R_{\rm mol}=\left(\frac{P_h}{P_0}\right)^\alpha,
\qquad
f_{\rm mol}=\frac{R_{\rm mol}}{1+R_{\rm mol}}.
\tag{5.12}
$$

The empirical relation of [Blitz & Rosolowsky (2006)][blitz] uses hydrostatic-pressure estimates. Replacing $P_h$ by an anisotropic, transient ram-pressure term is a testable extrapolation, not part of the original calibration. Molecular formation also requires time and depends on shielding and chemical conditions. At high molecular fraction, this route saturates.

Using both $f_{\rm mol}(P)$ and a pressure-shortened free-fall time is permissible only if they describe distinct phase and density changes. Otherwise it can count the same empirical SF response twice. The two-phase regulator in Section 8 makes that distinction explicit.

Turbulence has competing roles: it supplies support and changes the density distribution, including compressive density fluctuations. An increase in linewidth therefore cannot be assigned a universal suppressive or enhancing sign. [Federrath & Klessen (2012)][federrath] analyse these dependencies, including forcing and magnetic effects. The factor $s$ in equation (5.11) isolates support at fixed $f,e,d$; it does not assert that those factors remain fixed in a real turbulent medium.

For an optional energy closure, let the turbulent energy per area be $U_t=(3/2)\Sigma_g\sigma_t^2$ under an isotropic convention. A schematic local budget is

$$
\frac{dU_t}{dt}=\epsilon_w\frac{\rho_{\rm ICM}v_w^3}{2}
+\mathcal I_{\rm fb}-\frac{U_t}{t_{\rm diss}}
+\mathcal W_{\rm comp}-\mathcal L_{\rm adv}.
\tag{5.13}
$$

$\epsilon_w$ includes intercepted area and coupling, $\mathcal I_{\rm fb}$ is feedback input, and the final terms represent compression work and advective losses. With $t_{\rm diss}\sim L/\sigma_t$, dissipation scales as $\Sigma_g\sigma_t^3/L$. Wind power can thus raise support independently of SF. One cannot infer an SFR simply by assigning every extra pressure contribution to stellar-feedback production.

# 6. A modified Jeans criterion: pressure-confined clouds

## 6.1 Recover the Jeans result and identify its limitation

For a homogeneous, non-rotating medium with effective isotropic sound speed $c_c$, the linearized continuity, momentum, and Poisson equations give

$$
\frac{\partial^2\delta\rho}{\partial t^2}
=c_c^2\nabla^2\delta\rho+4\pi G\rho\,\delta\rho.
\tag{6.1}
$$

For a perturbation proportional to $\exp[i(\boldsymbol k\cdot\boldsymbol x-\omega t)]$,

$$
\omega^2=c_c^2k^2-4\pi G\rho.
\tag{6.2}
$$

An unstable mode has $\omega^2<0$. The marginal wavelength and the mass in a sphere of radius half that wavelength are

$$
\lambda_J=c_c\left(\frac{\pi}{G\rho}\right)^{1/2},
\qquad
M_J=\frac{\pi^{5/2}}{6}\frac{c_c^3}{G^{3/2}\rho^{1/2}}.
\tag{6.3}
$$

If compression raises density while support remains fixed, $M_J$ decreases. But an infinite uniform medium is not a pressure-bounded molecular cloud. The boundary problem matters; simply substituting an ICM pressure into $\rho=P/c_c^2$ does not establish that the cold cloud is actually confined at that pressure.

## 6.2 Show how external pressure enters the cloud virial balance

For a uniform spherical cloud of mass $M$ and radius $R_c$, omitting magnetic stresses, mass fluxes, and tides, a schematic scalar virial equation is

$$
\frac12\ddot I=3Mc_c^2-\frac{3GM^2}{5R_c}
-4\pi P_cR_c^3.
\tag{6.4}
$$

$P_c$ is the boundary pressure acting on the cloud. The last term is compressive. Its sign explains why an isolated-cloud virial parameter alone can misclassify a pressure-confined object. Equation (6.4) uses a uniform-density approximation; it does not give the exact critical coefficient of an isothermal sphere.

For a hydrostatic isothermal sphere, write

$$
\rho=\rho_c e^{-\psi},\quad r=r_0 x,\quad
r_0=\frac{c_c}{\sqrt{4\pi G\rho_c}}.
\tag{6.5}
$$

Combining hydrostatic balance with mass conservation produces the isothermal Lane-Emden equation,

$$
\frac{1}{x^2}\frac{d}{dx}\left(x^2\frac{d\psi}{dx}\right)=e^{-\psi},
\qquad \psi(0)=\psi'(0)=0.
\tag{6.6}
$$

At a boundary $x_b$, $P_c=c_c^2\rho_c e^{-\psi_b}$. Eliminating $\rho_c$ from the mass gives

$$
M=\frac{c_c^4}{G^{3/2}P_c^{1/2}}
\left[\frac{x_b^2\psi'_b e^{-\psi_b/2}}{\sqrt{4\pi}}\right].
\tag{6.7}
$$

The term in brackets reaches a maximum. Independent integration for this report gives $x_b=6.45075$ and coefficient $1.18223$, recovering the standard [Bonnor (1956)][bonnor] critical mass:

$$
\boxed{M_{\rm BE}=1.182\frac{c_c^4}{G^{3/2}P_c^{1/2}}.}
\tag{6.8}
$$

At fixed $M,c_c$, no stable isothermal equilibrium on the usual stable branch remains once

$$
P_c>P_{c,\rm crit}
\simeq1.397\frac{c_c^8}{G^3M^2}.
\tag{6.9}
$$

This is a concrete pressure-triggering criterion, subject to the cloud's boundary, thermal, and dynamical assumptions. Treating a turbulent linewidth as $c_c$ is an effective approximation, not an exact turbulent-cloud theorem.

## 6.3 Pressure and support compete very strongly

Let $p_c=P'_c/P_c$ and $s_c=c'_c/c_c$. Then

$$
\frac{M'_{\rm BE}}{M_{\rm BE}}=\frac{s_c^4}{\sqrt{p_c}},
\qquad
M'_{\rm BE}<M_{\rm BE}\ \Longleftrightarrow\ p_c>s_c^8.
\tag{6.10}
$$

A 20 percent increase in effective cloud support requires a pressure increase exceeding $1.2^8=4.30$ merely to lower the critical mass. By contrast, at a specified pressure-equilibrium density, shortening the free-fall time requires only $p_c>s_c^2$. These are different questions: a shorter timescale does not necessarily turn a previously stable cloud into an unstable one.

For illustration, $c_c=2\,{\rm km\,s^{-1}}$ and $P_c/k_B=10^5\,{\rm K\,cm^{-3}}$ give $M_{\rm BE}\simeq1.485\times10^4\,M_\odot$. The support speed here includes the chosen effective cloud support; it is not the thermal sound speed of cold molecular gas alone. It must not be confused with a kiloparsec-scale gas dispersion or an ionized-gas linewidth.

## 6.4 Turn a collapse threshold into a simple cloud-population prediction

A threshold alone is not an SFR. An optional bridge is to specify a cloud mass function and an available collapse time. Suppose, only as a toy assumption, that $dN/dM\propto M^{-2}$ between $M_{\min}$ and $M_{\max}$ and that all clouds above $M_{\rm BE}$ are eligible to collapse. Equal logarithmic mass intervals then contain equal mass, so

$$
f_{\rm eligible}=\left[
\frac{\ln\{M_{\max}/\max[M_{\min},M_{\rm BE}]\}}
{\ln(M_{\max}/M_{\min})}\right]_0^1.
\tag{6.11}
$$

The brackets truncate the result to $[0,1]$. At fixed support and away from the endpoints, increasing boundary pressure by $p_c$ changes this fraction by $\ln p_c/[2\ln(M_{\max}/M_{\min})]$. The response saturates when all clouds are eligible. This simple construction demonstrates why a power-law boost in SF need not continue indefinitely with pressure.

Eligibility is not actual conversion into stars. A finite collapse time, cloud disruption, turbulent density structure, magnetic criticality, and feedback still intervene. The cloud-population approach in [Fujita & Nagashima (1999)][fujita] is a more developed precedent; equation (6.11) is an explicitly simplified construction for this report, not its calibrated prescription.

## 6.5 Require time to collapse

At minimum, test whether

$$
t_{\rm ff,sf}\lesssim t_{\rm residence},\quad
t_{\rm ff,sf}\lesssim t_{\rm destruction},
\quad t_{\rm cool}\lesssim t_{\rm compression}
\tag{6.12}
$$

in the proposed cold phase. These are necessary plausibility checks, not jointly sufficient conditions for SF. A static pressure criterion is inappropriate when the pressure pulse ends or the cloud is mixed away before collapse develops. Anisotropic magnetic support, tides, and open cloud boundaries require a more complete virial or dynamical treatment.

# 7. Can Toomre Q be modified for ram pressure?

## 7.1 Derive the classical gaseous-disc criterion

In a thin, locally steady rotating gas disc, consider an axisymmetric radial disturbance $\propto e^{i(kx-\omega t)}$. The linearized surface continuity and in-plane momentum equations are

$$
\begin{aligned}
-i\omega\delta\Sigma+ik\Sigma_g\delta u&=0,\\
-i\omega\delta u-2\Omega\delta v&=
-ik\left(c_R^2\frac{\delta\Sigma}{\Sigma_g}+\delta\Phi\right),\\
-i\omega\delta v+\frac{\kappa^2}{2\Omega}\delta u&=0.
\end{aligned}
\tag{7.1}
$$

Here $c_R$ is effective radial random support, and

$$
\kappa^2=4\Omega^2+R\frac{d\Omega^2}{dR}.
\tag{7.2}
$$

For an infinitesimally thin sheet, $\delta\Phi=-2\pi G\delta\Sigma/|k|$. Eliminating the perturbed velocities gives

$$
\omega^2=\kappa^2+c_R^2 k^2-2\pi G\Sigma_g|k|.
\tag{7.3}
$$

Minimizing over positive $k$ gives $k_{\min}=\pi G\Sigma_g/c_R^2$ and

$$
\omega_{\min}^2=\kappa^2\left(1-Q_g^{-2}\right),
\qquad
\boxed{Q_g=\frac{\kappa c_R}{\pi G\Sigma_g}.}
\tag{7.4}
$$

The razor-thin fluid is axisymmetrically unstable when $Q_g<1$. This is the gaseous-disc criterion associated with the classical disc-stability literature, including [Goldreich & Lynden-Bell (1965)][goldreich]. [Toomre (1964)][toomre] treats a stellar disc; its familiar collisionless normalization uses 3.36 rather than the gas formula's $\pi$. Gas and stars cannot be combined by swapping these constants without a multicomponent response calculation.

For $Q_g<1$, the most unstable thin-disc mode has

$$
t_{\rm grow}=\frac{1}{\kappa\sqrt{Q_g^{-2}-1}},
\qquad
\lambda_{\rm fastest}=\frac{2c_R^2}{G\Sigma_g}.
\tag{7.5}
$$

Compare this wavelength with the aperture, radial background scale, and gas thickness. Compare the growth time with the forcing, advection, and removal times. A formally unstable mode that does not fit or cannot grow during the exposure is not a sufficient prediction of an SF enhancement.

## 7.2 The first important result: pure vertical compression leaves classical Q unchanged

External pressure does not appear explicitly in equation (7.4). If it only decreases the thickness at fixed $\Sigma_g,c_R,\kappa$, then classical $Q_g$ does not change. This does **not** mean that all three-dimensional gravitational stability is unchanged; it means the razor-thin formula has no remaining thickness information.

If the wind also changes gas column or in-plane support,

$$
\frac{Q'_g}{Q_g}=
\frac{\kappa'}{\kappa}\frac{c'_R}{c_R}\frac{\Sigma_g}{\Sigma'_g}.
\tag{7.6}
$$

In-plane accumulation lowers $Q_g$ at fixed other quantities. Gas removal raises it. Increased random support raises it. The change in $c_R$ need not equal the change in vertical $c_z$ or in internal cloud $c_c$.

For the illustrative values $\Sigma_g=20\,M_\odot\,{\rm pc}^{-2}$, $c_R=8\,{\rm km\,s^{-1}}$, and $\kappa=40\,{\rm km\,s^{-1}\,kpc^{-1}}$, $Q_g=1.184$. A 50 percent column increase with a 20 percent support increase gives $Q'_g=0.947$. This crosses the razor-thin threshold only; it is not a demonstration that a real two-component, thick, magnetized disc becomes star forming.

## 7.3 A transparent finite-thickness approximation

One approximate way to retain thickness is to reduce the self-gravity term by $F(kH)=1/(1+|k|H)$:

$$
\omega^2\simeq\kappa^2+c_R^2k^2-
\frac{2\pi G\Sigma_g|k|}{1+|k|H}.
\tag{7.7}
$$

This is a fixed-profile approximation, not a solution for deformable pressure boundaries. Set $x=kc_R/\kappa$ and $\eta_H=\kappa H/c_R$. Then

$$
\frac{\omega^2}{\kappa^2}=1+x^2-
\frac{2x}{Q_g(1+\eta_H x)}.
\tag{7.8}
$$

Requiring nonnegative frequency squared at every wavelength yields

$$
Q_{\rm crit}(\eta_H)=\max_{x>0}
\frac{2x}{(1+x^2)(1+\eta_H x)}.
\tag{7.9}
$$

Differentiation reduces the maximizing condition to

$$
1-x^2-2\eta_Hx^3=0.
\tag{7.10}
$$

The limit $H\rightarrow0$ recovers $Q_{\rm crit}=1$; finite thickness gives a smaller threshold. Thus vertical compression can destabilize this finite-thickness model by making self-gravity more effective even when the numerical value of classical $Q_g$ is unchanged. The profile and pressure boundary must still be physically appropriate.

![Original threshold calculations. Left: increased effective cloud support can offset a large confining-pressure increase in the Bonnor-Ebert threshold. Right: the critical gas Q in the stated frozen-profile thickness approximation approaches unity as the layer becomes thinner. These are distinct stability problems, not interchangeable SF laws.](assets/20260830_rps_sf/04_cloud_and_disk_thresholds.png){width=100%}

## 7.4 A pressure-dependent Q already exists: understand the boundary mode

An independent limiting calculation clarifies the pressure-confined-slab framework. Take an isothermal, self-gravitating equilibrium slab with symmetric constant external pressure $P_{\rm ext}$ and free surfaces at $z=\pm z_b$. Its density is

$$
\rho(z)=\rho_0\operatorname{sech}^2(z/h),
\quad h=\frac{c_s}{\sqrt{2\pi G\rho_0}},
\quad A=\tanh(z_b/h).
\tag{7.11}
$$

Integration and pressure matching give

$$
\Sigma_g=2\rho_0hA,
\qquad
A^2=\left(1+\frac{2P_{\rm ext}}{\pi G\Sigma_g^2}\right)^{-1}.
\tag{7.12}
$$

To extract a long-wavelength pressure-restoring speed, consider the auxiliary acoustic/surface response before adding perturbation self-gravity. The fundamental enthalpy perturbation $\delta h=\delta P/\rho$ is approximately vertically uniform. Internal compression contributes $\delta\Sigma_{\rm int}=\Sigma_g\delta h/c_s^2$. A free boundary can also move: constant exterior pressure requires $\delta h=g_b\xi_b$, with $g_b=2\pi G\Sigma_g$. The two moving surfaces therefore contribute $2\rho_b\delta h/g_b$ to the column perturbation.

The horizontal pressure force, including the external-pressure traction at those surfaces, has an effective perturbation $\delta\Pi=\Sigma_g\delta h$. Consequently,

$$
\begin{aligned}
c_{\rm mode}^2
&=\frac{\delta\Pi}{\delta\Sigma_{\rm int}+\delta\Sigma_{\rm surf}}\\
&=\frac{c_s^2}{1+2P_{\rm ext}/(g_b\Sigma_g)}
=\frac{c_s^2}{1+P_{\rm ext}/(\pi G\Sigma_g^2)}.
\end{aligned}
\tag{7.13}
$$

Defining a long-wavelength modal parameter gives

$$
Q_P\equiv\frac{\kappa c_{\rm mode}}{\pi G\Sigma_g}
=\frac{Q_g}{\sqrt{1+P_{\rm ext}/(\pi G\Sigma_g^2)}}.
\tag{7.14}
$$

This limiting construction is consistent with the free-boundary treatment in [Kim et al. (2012)][kim]. Its full dispersion relation requires the appropriate gravity-reduction factor and wavelength-dependent response. Their mixed-mode critical values are approximately 0.68-0.76 in their model family, rather than a universal threshold of one. Strongly confined modes can largely distort surfaces instead of strongly compressing gas; instability alone does not establish SF.

![Original illustration of the long-wavelength free-boundary factor in equation (7.14). The reduction describes a mode involving movable boundaries. It does not say that external pressure reduces a measured turbulent linewidth, and it does not supply a universal wind-pressure threshold for star formation.](assets/20260830_rps_sf/05_pressure_confined_mode.png){width=95%}

Several restrictions are decisive. The slab is symmetric, hydrostatic, isothermal, and vertically self-gravitating; its exterior pressure is prescribed. A stellar potential, asymmetric wind, magnetic stress, cooling, mass exchange, and finite response time change the problem. In particular, replacing $P_{\rm ext}$ by $\rho_{\rm ICM}v_w^2$ in equation (7.14), retaining a measured turbulent $c_R$, and declaring SF wherever $Q_P<1$ would combine incompatible assumptions.

## 7.5 What role should Q play in the proposed project?

Use classical or generalized $Q$ to ask whether a plausible gas layer can develop the relevant gravitational mode. Use a cloud criterion to ask whether dense material can collapse. Use the regulator to predict how much gas participates and how long the resulting SF lasts. A single $Q$ map cannot replace those three steps.

For strongly disturbed gas, infer $\kappa$ from a defensible underlying rotation model rather than interpreting wind-driven velocity gradients as circular rotation. Include stellar response where it matters. Treat regions that violate local steadiness or geometry assumptions as outside the simple stability model's domain, rather than assigning them an apparently precise $Q$.

# 8. A pressure-modified bathtub model with gas removal

## 8.1 Begin with local mass conservation

For a fixed disc-plane area, the cold-gas surface density obeys

$$
\frac{\partial\Sigma_g}{\partial t}
+\nabla_\parallel\cdot(\Sigma_g\boldsymbol u_\parallel)
=\Phi-(1-\mathcal R+\eta)\Sigma_{\rm SFR}
-\dot\Sigma_{\rm strip}+\dot\Sigma_{\rm fallback}.
\tag{8.1}
$$

$\Phi$ is gas supplied from outside the tracked cold reservoir; $\mathcal R$ is a returned stellar mass fraction under an instantaneous-recycling approximation; $\eta$ is the feedback-outflow mass loading relative to SFR. The divergence term describes in-plane transport. It can raise gas content in one region while reducing it elsewhere. Fallback is distinct from new supply.

This is the regulator logic of [Lilly et al. (2013)][lilly] extended here to a resolved open region. **Pressure is not an extra mass-source term.** Its role enters the relation between SFR and available gas, phase conversion, transport, or loss. Counting both compression and a gas source for the same advective accumulation would violate the bookkeeping.

For a one-region toy model, absorb specified net transport and fallback into an effective supply $\Phi_{\rm eff}$. Define

$$
a=1-\mathcal R+\eta,
\quad \Sigma_{\rm SFR}=\nu(t)\Sigma_g,
\quad \dot\Sigma_{\rm strip}=k(t)\Sigma_g.
\tag{8.2}
$$

Then

$$
\boxed{\frac{d\Sigma_g}{dt}=\Phi_{\rm eff}(t)
-[a\nu(t)+k(t)]\Sigma_g.}
\tag{8.3}
$$

Equation (5.9) is one possible closure for $\nu$. A measured molecular depletion time or a two-phase prescription can replace it. Each choice should be tested separately before adding parameters.

On very short timescales, instantaneous recycling is only an approximation. A more detailed model replaces $\mathcal R\Sigma_{\rm SFR}$ by the convolution of past SF with a stellar return kernel. The exact solutions below assume constant $a$ and should not be used to conceal that approximation.

## 8.2 Solve the constant-coefficient case step by step

Assume that before the perturbation the aperture is in equilibrium with supply $\Phi_0$, rate $\nu_0$, and no stripping:

$$
\Sigma_{g,0}=\frac{\Phi_0}{a\nu_0},
\qquad \Sigma_{{\rm SFR},0}=\nu_0\Sigma_{g,0}.
\tag{8.4}
$$

At $t=0$, change the coefficients to constants $\nu=b\nu_0$, $\Phi_{\rm eff}=F\Phi_0$, and $k\geq0$. Define

$$
g=\frac{\Sigma_g}{\Sigma_{g,0}},
\qquad \tau_0=\frac{1}{a\nu_0},
\qquad K=k\tau_0.
\tag{8.5}
$$

Equation (8.3) becomes

$$
\frac{dg}{dt}=\frac{F-(b+K)g}{\tau_0}.
\tag{8.6}
$$

Multiplying by the integrating factor $\exp[(b+K)t/\tau_0]$ and imposing $g(0)=1$ gives

$$
\begin{aligned}
g(t)&=g_{\rm eq}+(1-g_{\rm eq})e^{-t/\tau_{\rm eq}},\\
g_{\rm eq}&=\frac{F}{b+K},\qquad
\tau_{\rm eq}=\frac{\tau_0}{b+K}.
\end{aligned}
\tag{8.7}
$$

The SF response follows directly:

$$
\boxed{
\mathcal E(t)=b\left[\frac{F}{b+K}
+\left(1-\frac{F}{b+K}\right)e^{-t/\tau_{\rm eq}}\right].}
\tag{8.8}
$$

Initially the gas mass is unchanged, so $\mathcal E(0)=b$. A step in efficiency can produce a prompt burst. A physical finite response time smooths this idealized discontinuity.

At late times,

$$
\mathcal E_\infty=\frac{bF}{b+K}.
\tag{8.9}
$$

If supply has not increased ($F\leq1$), stripping is nonnegative, and the other assumptions hold, $\mathcal E_\infty\leq1$. With unchanged supply and no stripping, the long-term SFR returns to its previous equilibrium value; faster conversion reduces the equilibrium reservoir instead of indefinitely raising the SFR. With stripping, the asymptotic SFR is lower. Extra inward supply, $F>1$, can change this conclusion for a local aperture.

If $b>1$ but $\mathcal E_\infty<1$, setting equation (8.8) equal to one gives the duration above the original reference:

$$
t_{\rm cross}=\tau_{\rm eq}
\ln\left(\frac{b-\mathcal E_\infty}{1-\mathcal E_\infty}\right).
\tag{8.10}
$$

This is an observable-timescale prediction **conditional on the adopted coefficients and initial equilibrium**. It is not a universal burst lifetime.

## 8.3 An exact illustrative comparison

Take a baseline total-cold-gas depletion time $1/\nu_0=2$ Gyr, $\mathcal R=0.4$, and $\eta=0$, giving $a=0.6$ and $\tau_0=3.333$ Gyr. Keep $F=1$. The three cases in the next figure give:

| Imposed response | Late SFR ratio | Time above the reference |
|:--|:--|:--|
| $b=2$, no stripping | 1.000 | Approaches one from above |
| $b=2$, $k^{-1}=500$ Myr | 0.231 | 320 Myr |
| $b=1.1$, $k^{-1}=250$ Myr | 0.076 | 23.7 Myr |

![Original exact regulator solutions. A fixed pressure-associated efficiency boost can produce a transient SF excess while the gas reservoir declines. The curves keep their imposed coefficients constant; they are not orbit models. The analytic gas histories agree with direct numerical integration to better than three parts in a trillion in absolute normalized gas mass.](assets/20260830_rps_sf/02_exact_regulator.png){width=95%}

The labels describe chosen loss times, not a claim that one part of every stripped galaxy has those values. A constant-coefficient solution is useful because its assumptions and limiting behaviour can be checked exactly. Once $\nu$ depends on the changing gas weight, equation (8.8) is no longer the solution of that nonlinear problem.

## 8.4 An integrated constraint that does not depend on the SF recipe

Write equation (8.3) more generally as $\dot\Sigma_g=\Phi_{\rm eff}-a\Sigma_{\rm SFR}-L$, where $L$ is a specified gas loss rate. Subtract a reference evolution and integrate from $t_i$ to $t_f$. If initial gas contents are equal and $a$ is the same constant in both,

$$
\begin{aligned}
a\int_{t_i}^{t_f}\left(\Sigma_{\rm SFR}-\Sigma_{{\rm SFR},0}\right)dt
={}&\int_{t_i}^{t_f}(\Phi_{\rm eff}-\Phi_{{\rm eff},0})dt\\
&-\int_{t_i}^{t_f}(L-L_0)dt\\
&-[\Sigma_g(t_f)-\Sigma_{g,0}(t_f)].
\end{aligned}
\tag{8.11}
$$

An excess of newly formed stars must be supplied by extra input or a smaller remaining reservoir, after accounting for extra losses. If initial and final gas masses and total supply are the same, extra stripping reduces the integrated stellar production in this model. A temporary burst can rearrange when stars form without increasing their eventual total mass.

This conservation test is more robust than the chosen pressure-SFR exponent. It fails in the stated simple form if $a$ differs or stellar return is delayed, but the corresponding full mass budget must still close.

## 8.5 Extend to leading disc, trailing disc, and a tail without losing mass

Use masses $M_L,M_T,M_{\rm tail}$ in explicitly defined regions and total SFRs $\Psi_i=\nu_iM_i$. Let $\mathcal T_{LT}$ and $\mathcal T_{TL}$ be mass-transfer rates between the disc regions, $B_i$ fallback from the tracked tail, and $I_i$ other input. A conservative extension is

$$
\begin{aligned}
\dot M_L={}&I_L+B_L-a_L\Psi_L-k_LM_L
-\mathcal T_{LT}+\mathcal T_{TL},\\
\dot M_T={}&I_T+B_T-a_T\Psi_T-k_TM_T
+\mathcal T_{LT}-\mathcal T_{TL},\\
\dot M_{\rm tail}={}&k_LM_L+k_TM_T+\mathcal C_{\rm ICM}
-B_L-B_T-a_{\rm tail}\Psi_{\rm tail}-k_{\rm mix}M_{\rm tail}.
\end{aligned}
\tag{8.12}
$$

All stripped cold gas enters the tracked tail in this particular bookkeeping choice; material that mixes immediately can instead be assigned explicitly to an untracked hot component. $\mathcal C_{\rm ICM}$ allows new cold gas from cooling ambient/mixed material. Summing the equations cancels internal transfers, stripping-to-tail transfer, and fallback. Only external input, stellar production/outflows, and mixing out of the cold system remain.

This makes two points clear. First, a high trailing-disc SFR could come from transport or fallback rather than direct pressure triggering. Second, a bright tail need not imply that the disc itself is enhanced. The SF law and survival conditions for the tail must be specified separately.

## 8.6 An optional two-phase bathtub: atomic gas versus molecular gas

If the primary observational question is whether pressure creates more molecular gas while leaving its depletion time unchanged, split the cold reservoir into $\Sigma_a$ and $\Sigma_m$. Let $t_{\rm form}$ and $t_{\rm dis}$ describe formation and dispersal/dissociation of the molecular phase. A simple conservative system is

$$
\begin{aligned}
\dot\Sigma_a={}&\Phi-\frac{\Sigma_a}{t_{\rm form}}
+\frac{\Sigma_m}{t_{\rm dis}}
+\mathcal R\Sigma_{\rm SFR}-\eta_a\Sigma_{\rm SFR}-k_a\Sigma_a,\\
\dot\Sigma_m={}&\frac{\Sigma_a}{t_{\rm form}}
-\frac{\Sigma_m}{t_{\rm dis}}
-(1+\eta_m)\Sigma_{\rm SFR}-k_m\Sigma_m,\\
\Sigma_{\rm SFR}={}&\epsilon_{\rm ff}\Sigma_m/t_{\rm ff,m}.
\end{aligned}
\tag{8.13}
$$

Stellar return has been assigned to the atomic reservoir for this illustrative closure. Summing recovers the total-cold-gas regulator with $\eta=\eta_a+\eta_m$. Phase exchange cancels. Any phase-specific outflow prescription must preserve nonnegative reservoirs.

Pressure could reduce $t_{\rm form}$ while leaving $t_{\rm ff,m}/\epsilon_{\rm ff}$ unchanged. Then molecular content and SFR rise together at constant molecular SFE, as one possible explanation of a leading molecular excess. Alternatively, pressure can change $t_{\rm ff,m}$ or the bound fraction, producing a shorter molecular depletion time. The two mechanisms predict different CO/SFR behaviour even if their H-alpha maps look similar.

For the first observational test, this two-phase distinction is likely more informative than fitting several poorly constrained pressure-dependent stability parameters at once.

# 9. Time-dependent forcing, rotation, and the downstream response

## 9.1 A fully specified pressure-pulse toy example

The next example integrates equation (8.3), recalculating equation (5.4) as gas and support change. It has no between-region transfer and no fitted data. Let

$$
q(t)=\exp\left[-\frac{t^2}{2(80\,{\rm Myr})^2}\right],
\qquad
P_{\rm ram}(t)/k_B=1.2\times10^5q(t)\,{\rm K\,cm^{-3}}.
\tag{9.1}
$$

In each region set $P_{\rm conf}=P_{\rm bg}+\chi P_{\rm ram}$ and

$$
s(t)=\sqrt{1+u q(t)},\quad
k(t)=k_{\rm peak}q(t),\quad
b_{\rm target}=\frac{\sqrt{p(t)}}{s(t)},\quad
\dot b=\frac{b_{\rm target}-b}{15\,{\rm Myr}}.
\tag{9.2}
$$

The gas equation is $\dot g=a\nu_0[1-bg]-kg$, and $\mathcal E=bg$. The phase, intrinsic efficiency, and clumping ratios are fixed to one. Baseline values are $\Sigma_{g,0}=20\,M_\odot\,{\rm pc^{-2}}$, $\Sigma_*=100\,M_\odot\,{\rm pc^{-2}}$, $c_{z,0}=8\,{\rm km\,s^{-1}}$, $\sigma_{*,z}=30\,{\rm km\,s^{-1}}$, $P_{\rm bg}/k_B=10^4\,{\rm K\,cm^{-3}}$, $\nu_0^{-1}=2$ Gyr, and $a=0.6$. Supply stays at its baseline value.

| Region label | Confinement and support | Peak removal time |
|:--|:--|:--|
| Windward disc | $\chi=0.60$, $u=0.25$ | $k_{\rm peak}^{-1}=350$ Myr |
| Leeward disc | $\chi=0.12$, $u=1.00$ | $k_{\rm peak}^{-1}=160$ Myr |

**These choices intentionally construct a case with relatively retained, compressed windward gas and greater leeward loss/support.** They do not derive that ordering from the wind direction. Different couplings, transport, or loss rates can erase or reverse the asymmetry. A physical predictive application must constrain these inputs rather than choose them to match the desired answer.

![Original time-dependent example. Pressure peaks near time zero, while stripping changes the reservoir and therefore its weight. The windward SF maximum occurs before the ram-pressure maximum because subsequent gas loss overtakes the efficiency increase. The dotted curves apply an illustrative 100 Myr exponential response kernel; they are not calibrated FUV predictions.](assets/20260830_rps_sf/03_pressure_pulse.png){width=88%}

The computed windward SF maximum is $\mathcal E_L=1.164$ at $t=-65.3$ Myr. At maximum ram pressure, $\mathcal E_L=1.062$ while $\mathcal E_T=0.401$. At $t=200$ Myr, the values are 0.490 and 0.216: their contrast remains about 2.27, but both are suppressed. This demonstrates that an azimuthal contrast and absolute enhancement are separable, and that SF need not peak at the maximum pressure even with a positive response delay.

The leeward region never exceeds its baseline in this particular calculation. That is a consequence of the selected support, coupling, and loss histories, not a theorem of RPS. Conservation of gas plus accumulated sinks and supply closes to approximately $2.1\times10^{-15}$ in normalized mass in the numerical implementation. This verifies bookkeeping, not the physical validity of the chosen coefficients.

## 9.2 A simple prediction for the azimuthal offset of young SF

Suppose a small fractional perturbation $\delta=\delta\nu/\nu_0$ relaxes toward the first angular harmonic of a stationary pressure pattern. It also describes the fractional SFR response when gas column is fixed. At fixed radius,

$$
t_r\left(\frac{\partial\delta}{\partial t}
+\Omega_{\rm rel}\frac{\partial\delta}{\partial\phi}\right)
=A_1\cos\phi-\delta,
\tag{9.3}
$$

where $\Omega_{\rm rel}$ is the gas angular speed relative to the wind pattern. $A_1$ is a small forcing amplitude, not the full positive pressure field. This linearization treats the background reservoir and coefficients as slowly varying.

For a steady pattern, set the time derivative to zero and try $\delta=C\cos\phi+D\sin\phi$. Matching sine and cosine coefficients gives

$$
C=\frac{A_1}{1+(\Omega_{\rm rel}t_r)^2},\qquad
D=\frac{A_1\Omega_{\rm rel}t_r}{1+(\Omega_{\rm rel}t_r)^2}.
\tag{9.4}
$$

Equivalently,

$$
\boxed{
\delta=\frac{A_1}{\sqrt{1+(\Omega_{\rm rel}t_r)^2}}
\cos(\phi-\phi_{\rm lag}),\quad
\phi_{\rm lag}=\arctan(\Omega_{\rm rel}t_r).}
\tag{9.5}
$$

The response is attenuated and shifted downstream in the sense of rotation. For $v_c=200\,{\rm km\,s^{-1}}$ at $R=5$ kpc and a stationary wind pattern, $\Omega_{\rm rel}=0.04091\,{\rm Myr}^{-1}$. A 15 Myr response gives a $31.5$ degree shift and amplitude factor 0.852.

![Original first-harmonic response calculation. A longer response time moves the response peak away from the upstream direction and reduces the asymmetry amplitude. The offset is a conditional clock only when rotation, wind geometry, forcing evolution, and tracer response are independently constrained.](assets/20260830_rps_sf/06_azimuthal_lag.png){width=100%}

For harmonic number $m$, replace $\Omega_{\rm rel}t_r$ by $m\Omega_{\rm rel}t_r$ and divide the phase angle by $m$. A sharp pressure ridge is therefore smoothed as well as shifted. The half-wave pressure profile in equation (4.6) contains several harmonics; the first-harmonic result is not its exact full response.

An observed young-SF ridge need not lie exactly at the present pressure ridge. Gas rotation, collapse/chemistry delays, stellar motion, and a changing wind pattern can all contribute. A downstream ridge is consequently not sufficient evidence that RPS has failed to trigger SF. If gas column varies strongly with azimuth, fit the efficiency response or model the column separately; its SF phase need not equal the phase of $\nu$.

## 9.3 Observational tracers convolve the history

A measured SF indicator responds to a time-weighted history,

$$
\Sigma_{{\rm SFR},X}(t)=\int_0^\infty
K_X(\tau)\Sigma_{\rm SFR}(t-\tau)\,d\tau,
\qquad \int_0^\infty K_X(\tau)d\tau=1.
\tag{9.6}
$$

Spatially resolved measurements also require a kernel for drift and projection. The illustrative curves use an exponential kernel with a 100 Myr time constant; the supplied numerical table also includes a 5 Myr kernel. These are mathematical demonstrations, not detailed stellar-population response functions for FUV or H-alpha.

Comparing a short-lived massive-star tracer, a longer-lived UV tracer, young-cluster ages, and gas structure can help constrain whether an aperture is rising or fading. Dust, IMF sampling, mixed stellar populations, shocks, and gas/star displacement must be treated. [Kruijssen & Longmore (2014)][kruijssen] explain why small-aperture gas/SF ratios can be strongly affected by region lifecycles and sampling rather than an intrinsic efficiency change.

## 9.4 Why a tail need not be suppressed

For a cloud moving through hot gas, a useful hydrodynamical interaction time is

$$
t_{\rm cc}\sim\left(\frac{\rho_c}{\rho_{\rm hot}}\right)^{1/2}
\frac{R_c}{v_{\rm rel}}.
\tag{9.7}
$$

The fate of cold material depends on cooling and mixing, not only on this cloud-crushing time. [Gronke & Oh (2018)][gronke] show that cooling of mixed material can allow cold-gas growth under suitable conditions. Cold-gas survival or growth is not itself a criterion for stellar collapse.

Tail SF additionally needs sufficiently dense, cold, gravitationally unstable material that remains available long enough. Thermal ICM confinement can help maintain density; bulk ram pressure can supply or disturb the material. The tail calculations of [Tonnesen & Bryan (2012)][tonnesen2012] illustrate this distinction. A tail can contain young stellar knots, diffuse ionized gas with little SF, or both. Bright H-alpha outside the disc is not automatically evidence of SF because shocks and other ionization sources can contribute.

The model should therefore predict disc and tail SFRs separately. It should not impose a negative tail response as the complement of a positive leading response.

# 10. Observable predictions and what can be inferred

## 10.1 Predictions that distinguish mechanisms

The following are conditional, falsifiable predictions of specific closures, not universal signatures of RPS.

**Gas abundance or phase conversion dominates.** At matched radius and stellar background, $\Delta\log\Sigma_{\rm SFR}\simeq\Delta\log\Sigma_{\rm mol}$ and $\Delta\log t_{\rm dep,mol}\simeq0$. A molecular excess with no total-cold-gas excess favours phase conversion over net accumulation, provided both H I and molecular gas are measured consistently. The IC3949 analysis in [Roberts et al. (2022)][roberts] motivates this distinction.

**A pressure-shortened collapse time dominates.** With the other factors in equation (5.11) fixed,

$$
\Delta\log_{10}\nu=\frac12\Delta\log_{10}P_{\rm mid}
-\Delta\log_{10}c_z.
\tag{10.1}
$$

At fixed gas content, independently measured higher pressure with unchanged support predicts a shorter depletion time. A slope of one half is the hypothesis of this closure; it is not an empirical law already established for severely stripped discs. Changing the dense fraction, magnetic support, or $\epsilon_{\rm ff}$ changes the prediction.

**Support or disruption offsets compression.** The density or pressure can rise while depletion time does not shorten, or even lengthens. Cloud-scale linewidths, virial estimates, and dense-gas tracers help determine whether this is support, phase structure, or destruction. A larger ionized-gas linewidth alone cannot do so.

**Gas accumulation changes gravitational stability.** At fixed $\kappa$ and $c_R$, classical $Q_g$ decreases inversely with total gas column. Pure vertical compression leaves its value unchanged, while changing thickness or boundary modes. Atomic-to-molecular conversion can alter phase-specific support even if total gas mass is fixed; a single-fluid calculation then needs an appropriate effective response.

**The response has a finite clock.** The pressure/compression ridge and the youngest stellar ridge can have a rotation-dependent offset, with amplitude reduction as in equation (9.5). The short- and long-timescale SF tracers can disagree during a rise or decline. A transient model predicts a loop, rather than necessarily a single-valued relation, in SFR versus instantaneous pressure.

**Radial response is non-monotonic when survival competes with pressure.** At fixed coupling, $P_{\rm ram}/W$ is often more influential where internal weight is small, but those regions can also lose gas most readily. A retained annulus inside a depleted outer disc can therefore show the strongest SF response. Its position is not set by external pressure alone.

**A directional contrast persists after enhancement ends.** The example in Section 9 predicts $\mathcal E_L/\mathcal E_T>1$ even when both are below one. A sample analysis should therefore report absolute residuals and relative asymmetry separately.

## 10.2 Required measurements and their limits

| Measurement | Main constraint | Main caution |
|:--|:--|:--|
| Dust-corrected H-alpha with excitation classification | Recent SF distribution | Shocks, AGN, diffuse emission, and nondetections |
| Stellar continuum and population modelling | Stellar mass, reference structure, recent history | Age/dust degeneracy and mixed populations |
| Resolved CO plus a justified conversion | Molecular column and depletion time | Conversion factor, excitation, beam filling |
| H I at matched resolution | Atomic reservoir and total cold-gas retention | Beam size, sensitivity, projection |
| Resolved cold-gas line profiles | Random support and bulk flows | Beam smearing, unresolved components, anisotropy |
| Rotation model and stellar vertical structure | $\kappa$ and approximate weight | Disturbed streaming is not circular rotation |
| X-ray ICM density and temperature | Ambient density and thermal pressure | Three-dimensional location and substructure |
| Tail/radio geometry and young-region ages | Wind-direction priors and response timing | Historical winds, drift, projection |

Dense-gas lines such as HCN or HCO+ can add information, but their intensities depend on excitation, optical depth, abundances, and conversion factors. They are not automatic direct measurements of the bound mass fraction $f_{\rm sf}$.

For a MAUVE/MUSE-centred project, optical data can establish spatial SF residuals, excitation structure, stellar context, and some kinematics. **Optical data alone generally do not determine total cold-gas surface density, cloud support, gas Toomre Q, and ICM ram pressure separately.** The strongest physical discrimination comes from adding matched-resolution cold-gas measurements and independent environmental constraints.

## 10.3 Why H II-region pressure is not an ICM barometer

For a simple ionized-hydrogen plasma, an approximate thermal pressure is $P_{\rm HII}\simeq2n_e k_BT$. But an [S II] line ratio samples emitting warm gas; it can be density-weighted, affected by low-density sensitivity limits, and influenced by internal stellar feedback, shocks, and unresolved structure.

Thus

$$
P_{\rm HII}\ne P_{\rm mid}\ne P_{\rm conf}\ne P_{\rm ram}
\quad\hbox{in general}.
\tag{10.2}
$$

Likewise an observed ionized-gas linewidth contains instrumental, thermal, turbulent, beam-smearing, and coherent-flow contributions. It cannot be inserted directly as $c_R$, $c_z$, or $c_c$ in the three different models. An optical pressure or linewidth asymmetry is valuable evidence about the ionized phase; the conversion to a cold-gas or external pressure requires an additional physical model.

## 10.4 Conditional inference of an effective internal pressure

Invert equation (5.9):

$$
\boxed{
P_{\rm mid}=\frac{3\pi}{32G\xi}
\left(\frac{c_z\Sigma_{\rm SFR}}
{\epsilon_{\rm ff}f_{\rm sf}\Sigma_g}\right)^2.}
\tag{10.3}
$$

For ratios, equation (5.11) gives

$$
p=\frac{1}{d}\left(\frac{\mathcal E s}{gfe}\right)^2.
\tag{10.4}
$$

This estimates the pressure required by the assumed density/SF closure. It is not an independent measurement of that closure's validity. If all factors other than $p$ are fixed, an efficiency increase of 1.5 requires a pressure ratio of 2.25. That conclusion changes quadratically with an error in the assumed efficiency or support ratio.

In differential form,

$$
\delta\ln P_{\rm mid}=2\delta\ln\nu+2\delta\ln c_z
-2\delta\ln\epsilon_{\rm ff}-2\delta\ln f_{\rm sf}-\delta\ln\xi.
\tag{10.5}
$$

Uncertainties are often correlated. For example, $\Sigma_g$ and $f_{\rm mol}$ can share the same CO/H I measurements; they must not be treated as independent errors. Even independent 0.1 dex uncertainties in $\nu$ and $c_z$ already give about 0.28 dex uncertainty in inferred pressure before uncertainty in cloud structure, efficiency, geometry, and equilibrium is included.

To infer an effective wind contribution, one would further need

$$
P_{\rm dyn,eff}=P_{\rm mid}-W-P_{\rm bg}
\simeq\chi_{\rm eff}\rho_{\rm ICM}v_w^2.
\tag{10.6}
$$

Even a well-measured left-hand side constrains the product $\chi_{\rm eff}\rho v_w^2$, not its three factors individually. A negative inferred excess indicates inconsistency, uncertainty, or a breakdown of assumptions; it should not be silently clipped and presented as a measured zero pressure.

If ICM density and coupling are independently known, one can conditionally estimate $v_w=[P_{\rm dyn,eff}/(\chi_{\rm eff}\rho_{\rm ICM})]^{1/2}$. Usually coupling and three-dimensional location are not known well enough for a unique velocity. An SF map alone cannot recover a unique orbit or ram-pressure history.

## 10.5 Other potentially useful inferred quantities

**Response time from an azimuthal offset.** Under the first-harmonic, steady-pattern assumptions,

$$
t_r=\frac{\tan\phi_{\rm lag}}{\Omega_{\rm rel}}.
\tag{10.7}
$$

This offers a route to a response-time constraint from geometry plus rotation, rather than an unobservable instantaneous compression clock. It becomes poorly constrained if the wind direction, pattern evolution, deprojection, drift, or tracer kernel is uncertain. A displaced spiral arm can also mimic such a phase offset.

**Loss rate from a time-resolved reservoir estimate.** Rearranging equation (8.3) gives

$$
k=\frac{\Phi_{\rm eff}}{\Sigma_g}-a\nu-
\frac{d\ln\Sigma_g}{dt}.
\tag{10.8}
$$

A single image supplies neither $\Phi_{\rm eff}$ nor the gas-mass derivative. Inferring $k$ therefore requires additional information about transport, fallback, history, or a calibrated orbit model. Present SFR and present gas content do not uniquely give a stripping timescale.

**Cloud pressure benchmark.** Measured cloud masses, sizes, and internal support can be compared with equation (6.9). But $P_{c,\rm crit}\propto c_c^8/M^2$ is extremely sensitive to support: a 0.1 dex error in $c_c$ corresponds to 0.8 dex in the pressure benchmark at fixed mass. The presence of young stars does not prove that pressure crossed this threshold recently; the cloud might have been unstable before exposure.

**Dense-phase density under an SF closure.** Equation (5.8) gives $\rho_{\rm sf}=[3\pi/(32G)]\,[\nu/(\epsilon_{\rm ff}f_{\rm sf})]^2$. This can characterize the density required to explain an observed rate under an assumed efficiency. It is an effective model parameter, not a replacement for a density diagnostic.

## 10.6 Avoid circular tests

Do not infer $P_{\rm mid}$ from SFR using equation (10.3), then correlate that same SFR with the inferred pressure and call the result confirmation of equation (5.9). The relation was imposed in the inference.

To test the pressure hypothesis, estimate pressure or confinement from independent gas dynamics, cloud properties, vertical structure, or environmental measurements, with their covariances carried through. Similarly, an SF prescription that already contains pressure in a simulation cannot independently establish the same exponent from its resulting SF-pressure correlation; the EAGLE example in Section 3 illustrates the issue.

# 11. A dimensional worked example

This example is independent of the pressure-pulse parameters and is not a fit to a named galaxy. Its purpose is to show the scale of the quantities and how easily a modest gas loss can offset a pressure response.

Assume $n_{e,\rm ICM}=3\times10^{-4}\,{\rm cm^{-3}}$, mass per electron $1.17m_p$, $v_w=1000\,{\rm km\,s^{-1}}$, and $T_{\rm ICM}=2\times10^7$ K. Take the total particle density to be $1.95n_e$ for the thermal-pressure estimate. Equations (4.1) and (4.4) give

$$
P_{\rm ram}/k_B=4.252\times10^4\,{\rm K\,cm^{-3}},
\qquad
P_{\rm th}/k_B=1.170\times10^4\,{\rm K\,cm^{-3}}.
\tag{11.1}
$$

Assume, only for this calculation, that the thermal pressure is transmitted as the background confining pressure. For $\Sigma_g=20\,M_\odot\,{\rm pc^{-2}}$, $\Sigma_*=100\,M_\odot\,{\rm pc^{-2}}$, $c_z=8\,{\rm km\,s^{-1}}$, and $\sigma_{*,z}=30\,{\rm km\,s^{-1}}$, equation (5.4) gives

$$
W_0/k_B=3.091\times10^4\,{\rm K\,cm^{-3}},
\qquad
P_{{\rm mid},0}/k_B=4.261\times10^4\,{\rm K\,cm^{-3}}.
\tag{11.2}
$$

Now transmit an additional $0.5P_{\rm ram}$ while increasing vertical support by 10 percent, keeping the gas column fixed. The stellar-weight term must be updated because it contains $c_z$. It becomes $W'/k_B=3.268\times10^4\,{\rm K\,cm^{-3}}$, so

$$
P'_{\rm mid}/k_B=6.564\times10^4\,{\rm K\,cm^{-3}},
\quad p=1.540,
\quad b=\sqrt p/1.1=1.128.
\tag{11.3}
$$

With fixed phase fraction, cloud contrast, and efficiency per free-fall time, the fixed-column SF increase is only about 13 percent. If one compared a diagnostic point with $g=0.8$ at that same $p,s$, equation (5.11) would give $\mathcal E=0.903$. A true evolution to $g=0.8$ would also require recomputing its weight and therefore its pressure; it cannot simply reuse equation (11.3).

The baseline effective half-thickness is 73.6 pc and vertical crossing time is about 9.0 Myr. The mean midplane free-fall time is 22.0 Myr. This is not the free-fall time of dense clouds: $\rho_{\rm sf}=\xi\rho_{\rm mid}$ gives $t_{\rm ff,sf}=22.0/\sqrt\xi$ Myr. Those distinctions determine whether quasi-equilibrium and a 15 Myr response are plausible for a chosen phase.

This exercise demonstrates three safeguards: compare the extra forcing with existing internal weight; update support and weight consistently; and do not equate a pressure ratio with an SFR ratio.

# 12. A practical research programme for resolved data

## 12.1 A focused first question

A tractable initial question is:

> At matched stellar surface density and radius, is the windward SF excess explained by more molecular gas, or is an additional change in molecular depletion time required? Does either component vary with an independent measure of confinement relative to internal weight?

This question is more identifiable than trying to infer every galaxy's three-dimensional velocity from an H-alpha asymmetry. It also admits a meaningful result if the molecular depletion time does not change: gas accumulation or phase conversion is still a physically informative environmental response.

## 12.2 Analysis sequence

1. **Fix geometry independently of the SF result.** Define disc centre, projected position angle, inclination, upstream-direction proxy, disc/tail boundary, radial bins, and common resolution before searching for maximum SF contrast. Carry uncertainty in the wind direction. A divider chosen to maximize the measured SF difference needs a selection-aware null test.

2. **Preserve the project's actual SF selection.** Use the approved excitation and quality criteria consistently, document dust correction and sensitivity, and inspect whether a linewidth or equivalent-width cut preferentially removes compressed or shocked regions. Do not replace the existing MAUVE selection with guessed thresholds. No current MAUVE masks or FITS products were inspected for this report.

3. **Measure three distinct outcomes.** Report total SF per usable disc area, the area/fraction with detected SF under the stated selection, and the distribution among surviving SF regions. The usable-disc denominator must not be defined only by a detected SF tracer. Otherwise disappearance of SF regions is hidden by construction. Treat nondetections, non-SF excitation, and missing coverage distinctly.

4. **Construct both absolute and within-galaxy comparisons.** Match control conditions using stellar density and radius, and compare leading/trailing annuli within a galaxy. Internal pairing controls some galaxy-to-galaxy variation, but the opposite half can also be environmentally affected; it is not automatically a pristine control. Bars, spiral arms, interactions, and pre-existing asymmetry need consideration.

5. **Add gas before interpreting efficiency.** Convolve SF, CO, and H I products to compatible resolution and aperture definitions. Distinguish molecular and total-cold-gas depletion times. Carry CO conversion, nondetection, and geometry uncertainties. Optical-only results should remain SF-residual or excitation findings, not claims of measured cold-gas efficiency.

6. **Test models in increasing complexity.** Start with gas/phase abundance at unchanged molecular depletion time. Next test a pressure/support-dependent efficiency using independent pressure information. Add a finite response time or transport only if the data constrain it and improve predictive performance. Use cloud and disc stability calculations to reject physically inconsistent parameter combinations.

7. **Use independent statistical units.** Galaxies, not their many correlated spaxels, are the independent environmental realizations. Aggregate or fit a hierarchical model with galaxy-level effects; bootstrap whole galaxies for population uncertainty. Report galaxy scatter separately from uncertainty in a population median. Within-galaxy nulls can rotate a prescribed direction or permute labels while preserving radial structure and spatial correlation, with the limitations of such a null stated.

8. **Validate on information not used to set the parameters.** Predict the other disc half, another radial interval, young-region phase offsets, a gas map, or a withheld galaxy. Fitting arbitrary $\chi$, $k$, $f$, and $\epsilon_{\rm ff}$ independently in every aperture can reproduce almost any SF map and is not a discriminating test.

## 12.3 A minimal model-comparison hierarchy

**Model A: gas abundance.** Predict SF from the measured molecular gas and a matched-reference molecular depletion time. Its residual directly asks whether changed gas amount is sufficient.

**Model B: abundance plus a pressure/support response.** Test equation (10.1), or a fitted exponent around it, using independent pressure/support estimates. Shared measurement errors matter because estimated weight contains gas column. This model should be compared with Model A rather than assumed in advance.

**Model C: a time-dependent regulator.** Use a small set of galaxy-level forcing, response, and loss parameters, with radial or azimuthal dependencies tied to explicit geometry. Fit or constrain a common history using SF tracers of different response times, gas content, and morphology. Marginalize over unknown coupling rather than interpreting its best-fit value as a direct measurement.

The most useful initial output is likely a decomposition into gas-content, efficiency, and SF-area contributions, followed by conditional pressure or timescale constraints. A failure to identify pressure uniquely is an inference result, not a reason to hide the degeneracy.

## 12.4 What would constitute a theoretical advance?

Pressure-dependent SF and pressure-confined disc instability are already established ideas. Merely renaming equation (7.14) as an RPS-modified Toomre parameter would not establish novelty. Nor would multiplying a standard SF law by an unconstrained azimuthal factor.

A stronger contribution would connect independently constrained gas retention and pressure/support to **specific resolved predictions**, such as the sign and amplitude of residual molecular efficiency, the radial location of retained enhancement, or a response-time-dependent azimuthal offset. Demonstrating which mechanism explains the data, where a simple model fails, and which parameters cannot be inferred is scientifically valuable. Priority or novelty would still require a project-specific literature check before making a publication claim.

# 13. Conclusions and limits

The proposed project has a sound theoretical basis. The closest analytical precedents are pressure-dependent cloud/regulator models and rotating pressure-confined slabs, while resolved multiphase-disc work already demonstrates that compression, turbulence, and molecular content must be considered together.

For leading-side SF enhancement during severe RPS, the most useful starting point is **a conserving local regulator with separate pressure, support, gas-loss, and transport terms**. It can include a two-phase atomic/molecular reservoir. Bonnor-Ebert and disc-stability calculations constrain whether its assumed gas response is plausible; they do not provide the complete SFR prediction themselves.

The calculations here establish that pressure can lower a cloud-collapse threshold, that support can strongly offset this effect, that a gas reservoir can undergo a transient burst before suppression, and that rotation or gas loss can displace the SF maximum from the pressure maximum. None of these mechanisms forces the trailing disc or detached tail to be suppressed in every galaxy.

The most direct measurable distinction is between more molecular gas and a shorter molecular depletion time. A model can also provide a conditional effective pressure or response-time estimate. Recovering the actual ICM ram pressure, orbital speed, or stripping rate requires additional information about geometry, coupling, cloud structure, and history.

The numerical examples validate algebra and implementation under stated assumptions. They do not validate those assumptions for MAUVE. The recommended next scientific step is an observable decomposition and controlled model comparison, not an unqualified inversion of H-alpha into ram pressure.

# Appendix A. Notation and conventions

| Symbol | Meaning |
|:--|:--|
| $\Sigma_g$ | Total tracked cold-gas mass per fixed disc-plane area |
| $\Sigma_a,\Sigma_m$ | Atomic and molecular components of that reservoir |
| $\Sigma_{\rm SFR}$ | Mass converted into stars per time per area, before applying a returned fraction |
| $\nu$ | $\Sigma_{\rm SFR}/\Sigma_g$, inverse depletion time of the specified reservoir |
| $\Phi_{\rm eff}$ | Supply plus specified net transport/fallback into the aperture |
| $a$ | $1-\mathcal R+\eta$, under the stated recycling/outflow approximation |
| $k$ | Cold-gas stripping/removal coefficient, in inverse time |
| $P_{\rm ram}$ | Incident ICM momentum flux $\rho_{\rm ICM}v_w^2$ |
| $P_{\rm conf},P_{\rm bg}$ | Effective boundary confinement and its background component |
| $P_{\rm mid},W$ | Internal midplane pressure and vertical gravitational weight per area |
| $P_c$ | Boundary pressure of the cloud used in the Bonnor-Ebert problem |
| $c_z,c_R,c_c$ | Vertical gas support, radial gas support, internal cloud support; distinct quantities |
| $c_{\rm mode}$ | Speed entering an acoustic/surface mode, not an observed turbulent linewidth |
| $f_{\rm sf},\epsilon_{\rm ff},\xi$ | Star-forming mass fraction, conversion per free-fall time, density contrast |
| $g,p,s,f,e,d$ | Perturbed/reference ratios defined in equation (5.10) |
| $b$ | $\nu/\nu_0$; it equals $fe\sqrt{dp}/s$ only under the adopted closure |
| $F,K,\tau_0$ | Supply ratio, dimensionless loss coefficient, baseline regulator time |
| $\Omega,\kappa$ | Angular and epicyclic frequencies of the underlying rotation model |
| $Q_g,Q_P$ | Classical gaseous-disc parameter and the stated long-wavelength modal parameter |

Use consistent helium conventions when combining H I and CO-derived masses. The numerical examples use mass surface density directly and do not apply an additional hidden helium multiplier. Pressure is in ${\rm dyn\,cm^{-2}}$; quoted $P/k_B$ is in ${\rm K\,cm^{-3}}$. The support speeds are one-dimensional quantities under their respective effective assumptions. $1\,{\rm km\,s^{-1}\,kpc^{-1}}=1.02271\times10^{-3}\,{\rm Myr^{-1}}$.

Natural logarithms appear in time and threshold derivations; observational residuals marked $\log_{10}$ are in dex. The orbital wind inclination and the observer's line-of-sight disc inclination are different angles. All ratios require the same aperture/area convention and a stated reference.

# Appendix B. Reproducibility and numerical verification

The accompanying [reproduction script](assets/20260830_rps_sf/reproduce_models.py) generates the six figures, figure data tables, and [machine-readable checks](assets/20260830_rps_sf/numerical_checks.json). It requires Python, NumPy, SciPy, and Matplotlib. Running it writes only its own outputs beside the script. No observational files or external services are used by these calculations.

Verification was performed with Python 3.13.3, NumPy 2.2.6, SciPy 1.15.2, and Matplotlib 3.10.3. The following checks have distinct purposes:

| Check | Result and scope |
|:--|:--|
| Isothermal Lane-Emden integration | Critical radius 6.450751; mass coefficient 1.182227; density contrast 14.0420 |
| Thin-disc dispersion minimization | Maximum absolute residual relative to the analytic minimum below $9\times10^{-16}$ |
| Finite-thickness critical curve | Minimum-frequency residual below $3.5\times10^{-12}$ at sampled thresholds |
| Exact regulator versus ODE | Maximum absolute normalized gas difference below $2.4\times10^{-12}$ in the three cases |
| Pressure-pulse mass conservation | Maximum normalized budget residual below $2.1\times10^{-15}$ |
| Azimuthal response substitution | Maximum equation residual below $2.8\times10^{-16}$ |
| Pressure inversion over varied factors | Maximum relative recovery error below $8\times10^{-16}$ |

These checks confirm formula implementation, limiting behaviour, and mass bookkeeping. They do not calibrate coupling, stripping rates, efficiency, molecular conversion, or cloud structure. The script does not solve a full three-dimensional wind/ISM flow, magnetohydrodynamics, radiative chemistry, or the complete eigenvalue problem of a pressure-confined rotating slab. The plot of $Q_P/Q_g$ evaluates only the stated long-wavelength modal limit.

The Markdown and PDF are generated from the same report text. Mathematical expressions are typeset in the PDF, and all figures are original calculations for this report rather than reproduced journal figures.

# Appendix C. Research scope and source access

The literature search and source checks were completed on **30 August 2026**. This is a targeted theoretical review and model construction, not a formal exhaustive systematic review. Searches covered pressure-dependent SF, RPS and gas-regulator models, Jeans/Bonnor-Ebert collapse, pressure-confined rotating slabs, resolved leading-side observations, and simulations including contrary outcomes.

Primary journal pages, author-hosted manuscripts, original scans, and arXiv records were preferred. The source notes below distinguish inspected full-text sections from abstract/indexed-excerpt or bibliographic checks. Where the full paper could not be accessed, the discussion was restricted to supported claims; unavailable equations or numerical calibrations were not reconstructed from a secondary summary.

The supplied ChatGPT discussion's title, user questions, and the beginning of an embedded research response were recovered. The complete embedded English report could not be retrieved. It was used only to recover the intended topic and the request for English text with properly rendered mathematics; the scientific claims in this report were independently checked against primary sources.

A prior local [arbitrary-orientation stripping derivation](</Users/Igniz/Desktop/ICRAR/MAUVE/20260823 Finite-Thickness Arbitrary-Orientation Ram Pressure Stripping - Full Derivation.md>) was consulted for scope and the caveat about wind-driven torque. Its assumptions were not imported as universal constraints on the present SF models. No existing research notebook, FITS product, or repository was changed.

One bounded read-only assistant audit checked simulation evidence and selected bibliographic/erratum details. The main analysis owns the derivations, interpretation, numerical examples, and acceptance of the report. Full-text access does not imply that every page of every source was reviewed.

# References

All access dates are 30 August 2026. **FT** denotes primary full-text sections inspected; **A** denotes a primary abstract or indexed excerpt; **M** denotes bibliographic metadata checked. A source marked A was not treated as if its complete methods had been inspected. Publication year is used when it differs from the preprint year.

## Analytical models and physical foundations

1. **Fujita, Y., & Nagashima, M. (1999).** Effects of Ram Pressure from the Intracluster Medium on the Star Formation Rate of Disk Galaxies in Clusters of Galaxies. *ApJ*, 516, 619-625. DOI: 10.1086/307139. [Primary manuscript][fujita]. **FT:** cloud-population model and its gas-removal assumptions.

2. **Safarzadeh, M., & Loeb, A. (2019).** Explaining the enhanced star formation rate of Jellyfish galaxies in galaxy clusters. *MNRAS Letters*, 486, L26-L30. DOI: 10.1093/mnrasl/slz053. [Publisher article][safarzadeh]. **FT:** integrated pressure/efficiency and reservoir assumptions.

3. **Elmegreen, B. G., & Efremov, Y. N. (1997).** A Universal Formation Mechanism for Open and Globular Clusters in Turbulent Gas. *ApJ*, 480, 235-245. DOI: 10.1086/303966. [Publisher record][elmegreen]. **A/M:** author-institution abstract and publisher-deposited metadata; cloud-efficiency equations not independently rechecked.

4. **Vollmer, B., & Leroy, A. K. (2011).** Sustaining Star Formation Rates in Spiral Galaxies: Supernova-driven Turbulent Accretion Disk Models Applied to THINGS Galaxies. *AJ*, 141, 24. DOI: 10.1088/0004-6256/141/1/24. [Primary preprint][vollmer]. **A:** framework and physical ingredients.

5. **Nehlig, F., Vollmer, B., & Braine, J. (2016).** Effects of environmental gas compression on the multiphase ISM and star formation. *A&A*, 587, A108. DOI: 10.1051/0004-6361/201527021. [Primary preprint record][nehlig]. **A/M:** environmental findings and verified journal metadata; full PDF retrieval unavailable in this review.

6. **Lizee, T., Vollmer, B., Braine, J., & Nehlig, F. (2021).** Gas compression and stellar feedback in the tidally interacting and ram-pressure stripped Virgo spiral galaxy NGC4654. *A&A*, 645, A111. DOI: 10.1051/0004-6361/202038910. [Primary preprint][lizee]. **A:** abstract and indexed analytical-model section; full methods not comprehensively reviewed.

7. **Bonnor, W. B. (1956).** Boyle's Law and Gravitational Instability. *MNRAS*, 116, 351-359. DOI: 10.1093/mnras/116.3.351. [Original article scan][bonnor]. Historical foundation; the critical coefficient is independently reproduced numerically in this report.

8. **Goldreich, P., & Lynden-Bell, D. (1965).** I. Gravitational Stability of Uniformly Rotating Disks. *MNRAS*, 130, 97-124. DOI: 10.1093/mnras/130.2.97. [Publisher record][goldreich]. **A/M:** classical fluid-disc foundation; the thin-disc derivation is given explicitly here.

9. **Toomre, A. (1964).** On the Gravitational Stability of a Disk of Stars. *ApJ*, 139, 1217. DOI: 10.1086/147861. [Original article scan][toomre]. Historical stellar-disc reference; distinguished from the gaseous normalization in the text.

10. **Kim, J.-G., Kim, W.-T., Seo, Y. M., & Hong, S. S. (2012).** Gravitational Instability of Rotating, Pressure-Confined, Polytropic Gas Disks With Vertical Stratification. *ApJ*, 761, 131. DOI: 10.1088/0004-637X/761/2/131. [Primary manuscript][kim]. **FT:** equilibrium, boundary modes, effective speed, and stability definitions.

11. **Lubow, S. H., & Pringle, J. E. (1993).** The gravitational stability of a compressed slab of gas. *MNRAS*, 263, 701-706. DOI: 10.1093/mnras/263.3.701. [Publisher article][lubow]. **A:** interpretation of strongly confined distortion modes.

12. **Lilly, S. J., Carollo, C. M., Pipino, A., Renzini, A., & Peng, Y. (2013).** Gas Regulation of Galaxies: The Evolution of the Cosmic Specific Star Formation Rate, the Metallicity-Mass-Star-formation Rate Relation, and the Stellar Content of Halos. *ApJ*, 772, 119. DOI: 10.1088/0004-637X/772/2/119. [Primary preprint][lilly]. **A:** regulator framework; resolved directional extensions here are explicit constructions.

13. **Blitz, L., & Rosolowsky, E. (2006).** The Role of Pressure in GMC Formation. II. The H2-Pressure Relation. [Primary preprint, astro-ph/0605035][blitz]. **A:** empirical hydrostatic-pressure/molecular-fraction relation; no transient-wind calibration assumed.

14. **Ostriker, E. C., McKee, C. F., & Leroy, A. K. (2010).** Regulation of Star Formation Rates in Multiphase Galactic Disks: A Thermal/Dynamical Equilibrium Model. *ApJ*, 721, 975. DOI: 10.1088/0004-637X/721/2/975. [Primary preprint][ostriker]. **A:** equilibrium and feedback context.

15. **Krumholz, M. R., Dekel, A., & McKee, C. F. (2012).** A Universal, Local Star Formation Law in Galactic Clouds, Nearby Galaxies, High-Redshift Disks, and Starbursts. *ApJ*, 745, 69. DOI: 10.1088/0004-637X/745/1/69. [Author-hosted manuscript][krumholz]. **FT:** local law including the molecular fraction.

16. **Krumholz, M. R., Dekel, A., & McKee, C. F. (2013).** Erratum to the preceding paper. *ApJ*, 779, 89. DOI: 10.1088/0004-637X/779/1/89. [Author-hosted erratum][krumholzerratum]. **A:** indexed primary erratum text; corrected fitted efficiency 0.015, with the generic functional form retained.

17. **Federrath, C., & Klessen, R. S. (2012).** The Star Formation Rate of Turbulent Magnetized Clouds: Comparing Theory, Simulations, and Observations. *ApJ*, 761, 156. DOI: 10.1088/0004-637X/761/2/156. [Primary preprint][federrath]. **A:** turbulent forcing/support and magnetic dependencies.

18. **Gunn, J. E., & Gott, J. R. III (1972).** On the Infall of Matter Into Clusters of Galaxies and Some Effects on Their Evolution. *ApJ*, 176, 1-19. DOI: 10.1086/151605. [Original article scan][gunn]. **FT:** momentum flux and the stellar-sheet restoring estimate.

19. **Koppen, J., Jachym, P., Taylor, R., & Palous, J. (2018).** Ram pressure stripping made easy: an analytical approach. *MNRAS*, 479, 4367-4390. DOI: 10.1093/mnras/sty1610. [Primary preprint][koppen]. **A:** long-pulse versus impulse regimes.

20. **Whitworth, A. P. (2016).** A ram-pressure threshold for star formation. *MNRAS*, 458, 1815-1832. DOI: 10.1093/mnras/stw442. [Publisher article][whitworth]. **A:** converging molecular-flow context; not a universal ICM threshold.

21. **Pasetto, S., Cropper, M., Fujita, Y., Chiosi, C., & Grebel, E. K. (2015).** Environmental effects on star formation in dwarf galaxies and star clusters. *A&A*, 573, A48. DOI: 10.1051/0004-6361/201424502. [Primary preprint][pasetto]. **A:** idealized environmental-instability framework.

## Resolved observational constraints and interpretation

22. **Roberts, I. D., et al. (2022).** LoTSS Jellyfish Galaxies. IV. Enhanced Star Formation on the Leading Half of Cluster Galaxies and Gas Compression in IC3949. *ApJ*, 941, 77. DOI: 10.3847/1538-4357/ac9e9f. [Primary manuscript][roberts]. **FT:** leading-side evidence, aperture context, and molecular-gas interpretation.

23. **Vulcani, B., et al. (2020).** GASP. XXX. The Spatially Resolved SFR-Mass Relation in Stripping Galaxies in the Local Universe. *ApJ*, 899, 98. DOI: 10.3847/1538-4357/aba4ae. [Primary preprint][vulcani]. **A:** resolved SF excess at fixed stellar surface density.

24. **Brown, T., et al. (2023).** VERTICO. VII. Environmental Quenching Caused by Suppression of Molecular Gas Content and Star Formation Efficiency in Virgo Cluster Galaxies. [Primary preprint, arXiv:2308.10943][brown]. **A:** separation of gas content and efficiency; first author is Brown.

25. **Cakir, O., et al. (2026).** The SAMI Galaxy Survey: Quenching of Star Formation in Clusters III. Ram-Pressure-Affected Galaxy Populations. *PASA*, 43, e029. DOI: 10.1017/pasa.2026.10157. [Primary preprint][cakir]. **A/M:** outer suppression and absence of central excess in the selected comparison.

26. **Sun, J., et al. (2020).** Dynamical Equilibrium in the Molecular ISM in 28 Nearby Star-forming Galaxies. *ApJ*, 892, 148. DOI: 10.3847/1538-4357/ab781c. [Primary preprint][sun]. **A:** molecular pressure/weight context, not an ICM-pressure calibration.

27. **Kruijssen, J. M. D., & Longmore, S. N. (2014).** An uncertainty principle for star formation - I. Why galactic star formation relations break down below a certain spatial scale. *MNRAS*, 439, 3239. DOI: 10.1093/mnras/stu098. [Primary preprint][kruijssen]. **A:** lifecycle, sampling, and spatial-scale interpretation.

## Numerical evidence and cold-gas survival

28. **Troncoso-Iribarren, P., et al. (2020).** The better half - asymmetric star formation due to ram pressure in the EAGLE simulations. *MNRAS*, 497, 4145-4161. DOI: 10.1093/mnras/staa274. [Publisher article][eagle]. **FT:** pressure-dependent SF prescription, directional definitions, and aperture scope.

29. **Zhu, J., Tonnesen, S., & Bryan, G. L. (2024).** When and How Ram Pressure Stripping in Low-mass Satellite Galaxies Enhances Star Formation. *ApJ*, 960, 54. DOI: 10.3847/1538-4357/acfe6f. [Primary preprint text][zhu]. **FT:** version 1 of arXiv:2309.07037 inspected for methods and directional/time behaviour; final journal title/year used here.

30. **Bekki, K. (2014).** Galactic star formation enhanced and quenched by ram pressure in groups and clusters. *MNRAS*, 438, 444-462. DOI: 10.1093/mnras/stt2216. [Publisher article][bekki]. **FT:** dependence on wind/galaxy configuration and SF prescription.

31. **Tonnesen, S., & Bryan, G. L. (2009).** Gas Stripping in Simulated Galaxies with a Multiphase Interstellar Medium. *ApJ*, 694, 789-804. DOI: 10.1088/0004-637X/694/2/789. [Primary manuscript][tonnesen2009]. **FT:** gas-only treatment; no implemented SF.

32. **Tonnesen, S., & Bryan, G. L. (2012).** Star formation in ram pressure stripped galactic tails. *MNRAS*, 422, 1609-1624. DOI: 10.1111/j.1365-2966.2012.20737.x. [Publisher article][tonnesen2012]. **FT:** tail SF and thermal versus ram-pressure roles.

33. **Kapferer, W., et al. (2009).** The effect of ram pressure on the star formation, mass distribution and morphology of galaxies. *A&A*, 499, 87-102. [Publisher article][kapferer]. **FT:** disc versus wake SF and simulation assumptions, using the author-hosted manuscript.

34. **Lee, J., et al. (2020).** Dual Effects of Ram Pressure on Star Formation in Multiphase Disk Galaxies with Strong Stellar Feedback. *ApJ*, 905, 31. DOI: 10.3847/1538-4357/abc3b8. [Primary preprint][lee]. **A:** abstract-level evidence only; no detailed calibration adopted.

35. **Goller, J., Joshi, G. D., Rohr, E., Zinger, E., & Pillepich, A. (2023).** Jellyfish galaxies with the IllustrisTNG simulations - No enhanced population-wide star formation according to TNG50. *MNRAS*, 525, 3551-3570. DOI: 10.1093/mnras/stad2551. [Primary preprint][goller]. **FT:** primary journal manuscript checked for population-versus-history interpretation.

36. **Gronke, M., & Oh, S. P. (2018).** The growth and entrainment of cold gas in a hot wind. *MNRAS Letters*, 480, L111-L115. DOI: 10.1093/mnrasl/sly131. [Primary preprint][gronke]. **A:** cooling/mixing condition for cold-gas survival and growth, not an SF criterion.

[fujita]: https://arxiv.org/pdf/astro-ph/9812378
[safarzadeh]: https://academic.oup.com/mnrasl/article/486/1/L26/5454771
[elmegreen]: https://doi.org/10.1086/303966
[vollmer]: https://arxiv.org/abs/1009.3722
[nehlig]: https://arxiv.org/abs/1601.04883
[lizee]: https://arxiv.org/abs/2011.10531
[bonnor]: https://articles.adsabs.harvard.edu/pdf/1956MNRAS.116..351B
[goldreich]: https://academic.oup.com/mnras/article/130/2/97/2600978
[toomre]: https://adsabs.harvard.edu/pdf/1964ApJ...139.1217T
[kim]: https://arxiv.org/pdf/1210.6207
[lubow]: https://academic.oup.com/mnras/article/263/3/701/1206048
[lilly]: https://arxiv.org/abs/1303.5059
[blitz]: https://arxiv.org/abs/astro-ph/0605035
[ostriker]: https://arxiv.org/abs/1008.0410
[krumholz]: https://www.mso.anu.edu.au/~krumholz/publications/2012/krumholz12a.pdf
[krumholzerratum]: https://www.mso.anu.edu.au/~krumholz/publications/2012/krumholz12a_erratum.pdf
[federrath]: https://arxiv.org/abs/1209.2856
[gunn]: https://adsabs.harvard.edu/pdf/1972ApJ...176....1G
[koppen]: https://arxiv.org/abs/1806.05887
[whitworth]: https://academic.oup.com/mnras/article/458/2/1815/2589125
[pasetto]: https://arxiv.org/abs/1407.0024
[roberts]: https://arxiv.org/pdf/2210.16013
[vulcani]: https://arxiv.org/abs/2007.04996
[brown]: https://arxiv.org/abs/2308.10943
[cakir]: https://arxiv.org/abs/2602.03088
[sun]: https://arxiv.org/abs/2002.08964
[kruijssen]: https://arxiv.org/abs/1401.4459
[eagle]: https://academic.oup.com/mnras/article/497/4/4145/5729048
[zhu]: https://arxiv.org/html/2309.07037v1
[bekki]: https://academic.oup.com/mnras/article/438/1/444/1034466
[tonnesen2009]: https://arxiv.org/pdf/0901.2115
[tonnesen2012]: https://academic.oup.com/mnras/article/422/2/1609/1039842
[kapferer]: https://www.aanda.org/articles/aa/abs/2009/19/aa11551-08/aa11551-08.html
[lee]: https://arxiv.org/abs/2010.11028
[goller]: https://arxiv.org/abs/2304.09199
[gronke]: https://arxiv.org/abs/1806.02728
