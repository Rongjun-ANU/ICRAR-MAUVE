---
title: "Connecting H I stripping, star formation, and HOLMES line-ratio evolution"
subtitle: "Explicit derivation and photon-budget tests for MAUVE"
author: "Research report prepared for Rongjun Huang"
date: "1 October 2026"
lang: en
---

# 1. What needs to change in the interpretation

The proposed connection is physically meaningful under stated conditions. Removing atomic gas reduces the supply of molecular gas; the molecular reservoir subsequently declines through net consumption; the young stellar ionizing luminosity follows the declining SFR; and a retained old stellar population can provide a more slowly changing ionizing contribution. The old-star fraction of the Balmer luminosity then increases. If its emitting component has a higher forbidden-line-to-Balmer ratio, the observed ratio increases even while both emission lines become fainter.

The previous objection concerned two components powered by the **same young stellar population**, with a fixed division of its photons. Their luminosities fade proportionally, so their relative weights do not change. That objection does not apply when we add an independently supplied HOLMES component. Here HOLMES means hot low-mass evolved stars, including the post-asymptotic-giant-branch populations relevant to old-star photoionization.

However, four corrections are necessary before turning the proposal into a model.

1. Equality of the compact-HII and leaked-OB line ratios is an optional limiting approximation. Even an unchanged incident stellar spectrum can produce different ratios at different ionization parameters and gas conditions. Filtering can also change the leaked spectrum. We retain separate ratios; their equality is not needed for the derivation.
2. A harder ionizing spectrum does not by itself guarantee that every forbidden-line-to-Balmer ratio is larger. The ordering must be specified and tested separately for [N II], [S II], and [O III]. Metallicity, N/O, ionization parameter, and temperature matter; post-AGB calculations demonstrate these dependencies ([Byler et al. 2019](#ref-byler)).
3. An approximately constant old-star ionizing photon supply does not guarantee constant old-star Halpha emission. The gas must remain present, absorb the photons, and maintain the appropriate recombination conditions.
4. SF and NSF are observational selections. They are not switches that turn a physical source on or off. All three terms remain in the budget; neglecting one requires a quantitative small-contribution test.

There is a further distinction between **adding luminosities from separate emitting regions** and **illuminating the same gas with a mixed stellar spectrum**. The former permits a fixed-template linear mixture. The latter requires a joint ionization and thermal calculation. Belfiore et al. (2022), sections 5.4--5.5, use the latter approach for their spectral models. Their small galaxy-integrated HOLMES Halpha fraction does not imply a fixed small contribution to every forbidden line, or a fixed fraction in every local region ([Belfiore et al. 2022](#ref-belfiore)).

This report therefore develops a restricted, explicit three-component **emitting-region model**, connects it to the gas solution, and tests its normalization. It does not assume that obtaining the correct trend establishes a successful fit. In the MAUVE example below, the mechanism gives the expected direction for [N II]/Halpha, but its fiducial photon budget is insufficient to reproduce the full difference between the selected pre-peak and post-peak NSF means with the illustrative fixed spectra.

# 2. Definition of the local two-reservoir model

## 2.1 Quantities, mass conventions, and scope

We consider a projected position $\boldsymbol{x}$ in a galactic disc. Every surface density and parameter may depend on $\boldsymbol{x}$, but we suppress this argument while deriving the evolution of one local region. Time $t=0$ denotes the onset of the imposed environmental and supply conditions. It is not an observed time assigned to an infall-stage category.

All surface quantities must use the same area convention. The numerical MAUVE export inherits the inclination-corrected surface densities in the pipeline: both `SFR+Z.py` and `Mass.py` multiply the sky-projected density by their adopted axial-ratio factor. We therefore interpret the numerical quantities per that adopted disc area, not per uncorrected sky-projected area. The local photon-budget relation is unchanged when luminosity and stellar mass receive the same area correction. This source-code check does not independently revalidate every previously generated FITS header or stellar mass-to-light convention.

We denote the atomic and molecular phase mass surface densities by $\Sigma_{\mathrm{HI}}(t)$ and $\Sigma_{\mathrm{H_2}}(t)$. For the numerical gas calculation, both include associated helium consistently; thus the labels identify phases rather than hydrogen-only masses. A conversion to hydrogen-only columns would require the same conversion in the depletion time and all mass fluxes. Gas columns are quoted in $M_\odot\,\mathrm{pc}^{-2}$, times in Gyr, and observed SFR surface densities in $M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}$.

Following the resolved regulator notation in [Huang et al. (2026), section 4](#ref-huang), we define

$$
\Sigma_{\mathrm{SFR}}(t)\equiv
\frac{\Sigma_{\mathrm{H_2}}(t)}{\tau_{\mathrm{dep}}},
\qquad
\Sigma_\Phi(t)\equiv\frac{\Sigma_{\mathrm{HI}}(t)}{\tau_{\mathrm{conv}}}.
\tag{1}
$$

Here $\tau_{\mathrm{dep}}$ is the molecular depletion time defined relative to the total rate of star formation. The gas supply-rate surface density $\Sigma_\Phi$ feeds the molecular reservoir. The conversion time $\tau_{\mathrm{conv}}$ describes the assumed net transfer from the local atomic phase. The second equality is **our linear closure**, not a measured conversion law or an equation established by Huang et al. Their replenishment timescale is instead

$$
\tau_\Phi(t)\equiv
\frac{\Sigma_{\mathrm{H_2}}(t)}{\Sigma_\Phi(t)}
=\tau_{\mathrm{conv}}
\frac{\Sigma_{\mathrm{H_2}}(t)}{\Sigma_{\mathrm{HI}}(t)}.
\tag{2}
$$

Consequently, $\Sigma_\Phi$ is a rate, whereas $\tau_\Phi$ and $\tau_{\mathrm{conv}}$ are times. Even when $\tau_{\mathrm{conv}}$ is constant, $\tau_\Phi$ generally evolves. We do not impose $\tau_\Phi\leq\tau_{\mathrm{dep}}$ on a supply-starved system; that inequality is not a consequence of the definition.

Let $R$ be the prompt stellar mass return fraction and $\lambda$ the feedback mass-loading factor, so that the feedback mass loss is $\lambda\Sigma_{\mathrm{SFR}}$. The net removal associated with star formation and feedback is $(1-R+\lambda)\Sigma_{\mathrm{SFR}}$. This is the usual regulator bookkeeping, here assigned effectively to the molecular reservoir ([Lilly et al. 2013](#ref-lilly); [Huang et al. 2026](#ref-huang)). Treating recycled material as promptly available to this reservoir is an approximation; explicit phase-dependent recycling would require additional terms.

We assume positive initial gas columns and conversion/depletion times, $0\leq R<1$, $\lambda\geq0$, and $\gamma_{\mathrm{strip}}\geq0$. The resulting reservoir response rates are positive, as required for the decline and peak statements below.

We use only $\gamma_{\mathrm{strip}}$ for the direct RPS loss coefficient. It has units of inverse time and acts only on HI. The baseline has no direct molecular stripping, no external supply to HI after $t=0$, and no explicit lateral transport. All coefficients are initially constant. A linear molecular law is an empirical first approximation in nearby discs, with substantial environmental and scale-dependent limitations ([Leroy et al. 2013](#ref-leroy)). The HI-only stripping choice is a hypothesis for this calculation, not a claim that molecular stripping never occurs; observations provide counterexamples ([Boselli et al. 2014](#ref-boselli)).

## 2.2 Continuity equations and the meaning of the rates

Removing each mass flux from the reservoir that supplies it gives

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{HI}}}{dt}
&=-\Sigma_\Phi-\gamma_{\mathrm{strip}}\Sigma_{\mathrm{HI}},\\
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
&=\Sigma_\Phi-(1-R+\lambda)\Sigma_{\mathrm{SFR}}.
\end{aligned}
\tag{3}
$$

These are the model's local mass balances. Their sum explicitly cancels the internal transfer:

$$
\frac{d}{dt}(\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H_2}})
=-\gamma_{\mathrm{strip}}\Sigma_{\mathrm{HI}}
-(1-R+\lambda)\Sigma_{\mathrm{SFR}}.
\tag{4}
$$

Thus conversion does not destroy gas. The total cold reservoir declines through stripping and net stellar/feedback consumption. We define two physically labelled response rates,

$$
\gamma_{\mathrm{HI}}\equiv\frac{1}{\tau_{\mathrm{conv}}}+\gamma_{\mathrm{strip}},
\qquad
\gamma_{\mathrm{H_2}}\equiv\frac{1-R+\lambda}{\tau_{\mathrm{dep}}}.
\tag{5}
$$

The first is the total fractional removal rate from the HI reservoir, including conversion. The second is the net molecular consumption rate. It is not an H2 stripping rate. Both use the same rate notation; $\gamma_{\mathrm{strip}}$ remains the only direct stripping coefficient.

Equations (1)--(5) use consistent physical units. The explicit conversion needed in the numerical tables is

$$
\frac{\Sigma_{\mathrm{SFR}}}{M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}}
=10^{-3}
\frac{\Sigma_{\mathrm{H_2}}/(M_\odot\,\mathrm{pc}^{-2})}
{\tau_{\mathrm{dep}}/\mathrm{Gyr}}.
\tag{6}
$$

The factor is $10^6/10^9$: square parsecs to square kiloparsecs, then Gyr to years.

# 3. Explicit solution and the conditions for SFR decline or enhancement

## 3.1 Solve the atomic reservoir

For positive initial column $\Sigma_{\mathrm{HI},0}$, substituting equation (1) into equation (3), separating the variables, and integrating from the initial state gives

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{HI}}}{dt}&=-\gamma_{\mathrm{HI}}\Sigma_{\mathrm{HI}},\\
\int_{\Sigma_{\mathrm{HI},0}}^{\Sigma_{\mathrm{HI}}(t)}
\frac{d\widetilde\Sigma_{\mathrm{HI}}}{\widetilde\Sigma_{\mathrm{HI}}}
&=-\gamma_{\mathrm{HI}}\int_0^t du,\\
\ln\frac{\Sigma_{\mathrm{HI}}(t)}{\Sigma_{\mathrm{HI},0}}
&=-\gamma_{\mathrm{HI}}t.
\end{aligned}
\tag{7}
$$

Here $u$ and $\widetilde\Sigma_{\mathrm{HI}}$ are integration variables. Exponentiating and applying the conversion law yields

$$
\boxed{\Sigma_{\mathrm{HI}}(t)=\Sigma_{\mathrm{HI},0}e^{-\gamma_{\mathrm{HI}}t}},
\qquad
\Sigma_\Phi(t)=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\tag{8}
$$

The decline of $\Sigma_\Phi$ is the indirect route by which atomic stripping affects molecular gas. The physical possibility of molecular depletion following atomic deficiency has observational and theoretical antecedents ([Fumagalli et al. 2009](#ref-fumagalli)); the exact exponential form here follows from our stated closure.

## 3.2 Solve the molecular reservoir without omitting the integrating-factor algebra

First, insert equation (8) into the second balance:

$$
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
+\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2}}
=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\tag{9}
$$

Multiplication by $e^{\gamma_{\mathrm{H_2}}t}$ makes the left-hand side a product derivative. To see this explicitly,

$$
\begin{aligned}
\frac{d}{dt}\left[e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2}}(t)\right]
&=e^{\gamma_{\mathrm{H_2}}t}\frac{d\Sigma_{\mathrm{H_2}}}{dt}
+\gamma_{\mathrm{H_2}}e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2}},\\
&=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})t}.
\end{aligned}
\tag{10}
$$

Integrating the product derivative from $0$ to $t$ gives

$$
e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2}}(t)-\Sigma_{\mathrm{H_2},0}
=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\int_0^t e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})u}\,du.
\tag{11}
$$

For unequal rates, the integral is

$$
\begin{aligned}
\int_0^t e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})u}\,du
&=\left[
\frac{e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})u}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}\right]_0^t\\
&=\frac{e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})t}-1}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}
\tag{12}
$$

Now add the initial column, multiply by $e^{-\gamma_{\mathrm{H_2}}t}$, and distribute the exponential:

$$
\begin{aligned}
\Sigma_{\mathrm{H_2}}(t)
&=e^{-\gamma_{\mathrm{H_2}}t}
\left[\Sigma_{\mathrm{H_2},0}
+\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\frac{e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})t}-1}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}\right]\\
&=\Sigma_{\mathrm{H_2},0}e^{-\gamma_{\mathrm{H_2}}t}
+\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}
\tag{13}
$$

The first term is the surviving initial molecular gas. The second is molecular gas supplied after $t=0$, with its subsequent consumption included. Although individual coefficients in an exponential expansion can be negative, this supplied-gas term is nonnegative: its numerator and denominator always have the same sign.

Finally, dividing every term by the constant depletion time gives the requested explicit SFR solution:

$$
\boxed{
\begin{aligned}
\Sigma_{\mathrm{SFR}}(t)
&=\frac{\Sigma_{\mathrm{H_2},0}}{\tau_{\mathrm{dep}}}
e^{-\gamma_{\mathrm{H_2}}t}\\
&\quad+\frac{\Sigma_{\mathrm{HI},0}}
{\tau_{\mathrm{conv}}\tau_{\mathrm{dep}}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}}
\tag{14}
$$

All results in equations (7)--(14) are direct integrations of this report's balances. They are not additional empirical laws.

If the two rates are equal, the integrand in equation (11) is unity. There is no physical divergence:

$$
\begin{aligned}
\Sigma_{\mathrm{H_2}}(t)
&=\left(\Sigma_{\mathrm{H_2},0}
+\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}t\right)
e^{-\gamma_{\mathrm{H_2}}t},\\
\Sigma_{\mathrm{SFR}}(t)
&=\frac{1}{\tau_{\mathrm{dep}}}
\left(\Sigma_{\mathrm{H_2},0}
+\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}t\right)
e^{-\gamma_{\mathrm{H_2}}t}.
\end{aligned}
\tag{15}
$$

## 3.3 What determines whether the SFR initially rises or falls?

Define the initial replenishment time $\tau_{\Phi,0}=\Sigma_{\mathrm{H_2},0}/\Sigma_{\Phi,0}$ and the dimensionless remaining-SFR fraction $F_{\mathrm{SFR}}(t)=\Sigma_{\mathrm{SFR}}(t)/\Sigma_{\mathrm{SFR},0}$. Equation (14) becomes

$$
F_{\mathrm{SFR}}(t)=e^{-\gamma_{\mathrm{H_2}}t}
+\frac{1}{\tau_{\Phi,0}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\tag{16}
$$

More directly, divide equation (9) by the molecular column and use constant $\tau_{\mathrm{dep}}$:

$$
\begin{aligned}
\frac{1}{\Sigma_{\mathrm{SFR}}}\frac{d\Sigma_{\mathrm{SFR}}}{dt}
&=\frac{1}{\Sigma_{\mathrm{H_2}}}\frac{d\Sigma_{\mathrm{H_2}}}{dt}\\
&=\frac{\Sigma_\Phi}{\Sigma_{\mathrm{H_2}}}-\gamma_{\mathrm{H_2}}\\
&=\frac{1}{\tau_\Phi(t)}-\gamma_{\mathrm{H_2}}.
\end{aligned}
\tag{17}
$$

Thus the sign follows from supply relative to consumption:

$$
\begin{aligned}
\left.\frac{d\Sigma_{\mathrm{SFR}}}{dt}\right|_0>0
&\ \Longleftrightarrow\ 
\Sigma_{\Phi,0}>\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0}
\ \Longleftrightarrow\ 
\tau_{\Phi,0}<\frac{\tau_{\mathrm{dep}}}{1-R+\lambda},\\
\left.\frac{d\Sigma_{\mathrm{SFR}}}{dt}\right|_0<0
&\ \Longleftrightarrow\ 
\tau_{\Phi,0}>\frac{\tau_{\mathrm{dep}}}{1-R+\lambda}.
\end{aligned}
\tag{18}
$$

Equality gives zero initial slope, not permanent equilibrium. If the initial molecular reservoir is balanced, substitution of $1/\tau_{\Phi,0}=\gamma_{\mathrm{H_2}}$ gives

$$
\begin{aligned}
F_{\mathrm{SFR}}(t)
&=\frac{\gamma_{\mathrm{HI}}e^{-\gamma_{\mathrm{H_2}}t}
-\gamma_{\mathrm{H_2}}e^{-\gamma_{\mathrm{HI}}t}}
{\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}}},\\
\frac{dF_{\mathrm{SFR}}}{dt}
&=\frac{\gamma_{\mathrm{HI}}\gamma_{\mathrm{H_2}}}
{\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}}}
\left(e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}\right)<0
\quad(t>0).
\end{aligned}
\tag{19}
$$

The derivative is negative for either ordering of the positive rates. At the onset, the first derivative vanishes but the curvature is negative:

$$
F_{\mathrm{SFR}}(t)
=1-\frac{\gamma_{\mathrm{HI}}\gamma_{\mathrm{H_2}}}{2}t^2
+O(t^3).
\tag{20}
$$

This explains the delayed decrease. Atomic removal immediately reduces future supply, but it does not instantly destroy the molecular reservoir or its ongoing star formation.

At late times the slower exponential usually dominates. In the common regime $\gamma_{\mathrm{HI}}>\gamma_{\mathrm{H_2}}$, the molecular response approaches an exponential with timescale $1/\gamma_{\mathrm{H_2}}=\tau_{\mathrm{dep}}/(1-R+\lambda)$. It is therefore exponential-like, but generally not one exponential from the onset. For any nonnegative supply, equation (13) also implies

$$
F_{\mathrm{SFR}}(t)\geq e^{-\gamma_{\mathrm{H_2}}t}.
\tag{21}
$$

Increasing HI stripping cannot make retained H2 disappear faster than the no-supply consumption solution under these assumptions. This bound is an important check on attempts to fit rapid quenching with HI-only loss.

## 3.4 Derive the enhancement peak explicitly

Suppose the initial supply exceeds molecular consumption. Differentiating equation (16) and collecting terms gives

$$
\frac{dF_{\mathrm{SFR}}}{dt}
=\frac{
\gamma_{\mathrm{H_2}}[\tau_{\Phi,0}^{-1}+\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}}]
e^{-\gamma_{\mathrm{H_2}}t}
-\gamma_{\mathrm{HI}}\tau_{\Phi,0}^{-1}e^{-\gamma_{\mathrm{HI}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\tag{22}
$$

At a stationary point $t_{\mathrm{peak}}$, the numerator vanishes. Moving one term to the other side and dividing by the positive exponential gives

$$
\begin{aligned}
\gamma_{\mathrm{H_2}}
[1+(\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}})\tau_{\Phi,0}]
e^{-\gamma_{\mathrm{H_2}}t_{\mathrm{peak}}}
&=\gamma_{\mathrm{HI}}e^{-\gamma_{\mathrm{HI}}t_{\mathrm{peak}}},\\
e^{(\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}})t_{\mathrm{peak}}}
&=\frac{\gamma_{\mathrm{HI}}}
{\gamma_{\mathrm{H_2}}[1+(\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}})\tau_{\Phi,0}]}.
\end{aligned}
\tag{23}
$$

Taking the natural logarithm yields

$$
\boxed{
t_{\mathrm{peak}}=
\frac{1}{\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}}}
\ln\left[
\frac{\gamma_{\mathrm{HI}}}
{\gamma_{\mathrm{H_2}}[1+(\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}})\tau_{\Phi,0}]}
\right]}.
\tag{24}
$$

This is a positive-time peak when $\tau_{\Phi,0}^{-1}>\gamma_{\mathrm{H_2}}$. For equal rates, differentiating $(1+t/\tau_{\Phi,0})e^{-\gamma_{\mathrm{H_2}}t}$ gives the finite limit $t_{\mathrm{peak}}=1/\gamma_{\mathrm{H_2}}-\tau_{\Phi,0}$.

To verify that the stationary point is a maximum, differentiate the molecular balance once more:

$$
\begin{aligned}
\frac{d^2\Sigma_{\mathrm{H_2}}}{dt^2}
&=-\gamma_{\mathrm{HI}}\Sigma_\Phi
-\gamma_{\mathrm{H_2}}\frac{d\Sigma_{\mathrm{H_2}}}{dt},\\
\left.\frac{d^2\Sigma_{\mathrm{H_2}}}{dt^2}\right|_{t_{\mathrm{peak}}}
&=-\gamma_{\mathrm{HI}}\Sigma_\Phi(t_{\mathrm{peak}})<0.
\end{aligned}
\tag{25}
$$

Every stationary point has negative curvature. Therefore a solution starting with positive slope has one maximum followed by decline; a solution starting with nonpositive slope cannot develop a later minimum and rise within this constant-coefficient closed model.

## 3.5 What could produce the initial excess supply?

Turning on $\gamma_{\mathrm{strip}}$ alone does not change the initial $\Sigma_\Phi$ at fixed atomic column and conversion time. It therefore cannot produce an initial enhancement from an exactly balanced molecular state. A separate response to the environment is needed.

One possibility is a shorter conversion time in compressed gas. Let the pre-onset state satisfy $\Sigma_{\mathrm{HI},0}/\tau_{\mathrm{conv,pre}}=\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0}$. If conversion becomes faster while gas columns and $\tau_{\mathrm{dep}}$ remain initially continuous, then

$$
\begin{aligned}
\frac{1}{\Sigma_{\mathrm{SFR},0}}
\left.\frac{d\Sigma_{\mathrm{SFR}}}{dt}\right|_{0^+}
&=\frac{\Sigma_{\mathrm{HI},0}}
{\tau_{\mathrm{conv,post}}\Sigma_{\mathrm{H_2},0}}-\gamma_{\mathrm{H_2}}\\
&=\gamma_{\mathrm{H_2}}
\left(\frac{\tau_{\mathrm{conv,pre}}}{\tau_{\mathrm{conv,post}}}-1\right).
\end{aligned}
\tag{26}
$$

The slope becomes positive if $\tau_{\mathrm{conv,post}}<\tau_{\mathrm{conv,pre}}$. This is an explicit phenomenological compression prescription. The gas equations do not derive the response of $\tau_{\mathrm{conv}}$ to ram pressure.

A second possibility is a decrease in $\tau_{\mathrm{dep}}$. For a time-dependent depletion time, differentiating equation (1) gives

$$
\frac{d\ln\Sigma_{\mathrm{SFR}}}{dt}
=\frac{1}{\tau_\Phi}
-\frac{1-R+\lambda}{\tau_{\mathrm{dep}}}
-\frac{d\ln\tau_{\mathrm{dep}}}{dt}.
\tag{27}
$$

A sufficiently rapid decrease of $\tau_{\mathrm{dep}}$ can raise SFR even before the molecular column increases. An idealized instantaneous change from 2 to 1.5 Gyr at fixed column raises SFR by $2/1.5$, or 0.125 dex; a real response should be smoothed over a physical compression/cloud-evolution time.

There is also a mass-budget ceiling for conversion-driven enhancement at fixed depletion time. Integrating equation (4), with all sinks nonnegative, gives

$$
\Sigma_{\mathrm{H_2}}(t)\leq
\Sigma_{\mathrm{H_2},0}+\Sigma_{\mathrm{HI},0}
\quad\Longrightarrow\quad
F_{\mathrm{SFR}}(t)\leq
1+\frac{\Sigma_{\mathrm{HI},0}}{\Sigma_{\mathrm{H_2},0}}.
\tag{28}
$$

This ceiling is generally not attained because gas is being stripped and consumed. A molecular-dominated patch has little atomic material available for conversion, so faster conversion alone may yield only a weak enhancement.

These conditions can accommodate a potential local enhancement such as that discussed for NGC4654. They do not establish its cause or amplitude. A facing-to-opposite contrast can also arise from stronger suppression on the opposite side. The earlier gradient diagnostic is motivation; it is not used as a fitted constraint in this revision.

# 4. From the SFR solution to young-star Halpha emission

## 4.1 The young stellar population has a finite response time

Halpha traces the ionizing photons supplied by short-lived massive stars, not an infinitely instantaneous SFR. A stellar population produces an age-dependent ionizing output; summing stellar cohorts is therefore a convolution. This is the population-synthesis basis of recombination-line SFR indicators ([Kennicutt & Evans 2012](#ref-ke)). We approximate the response with a normalized kernel $K_\alpha(a)$, where $a\geq0$ is stellar age:

$$
\overline\Sigma_{\mathrm{SFR},\alpha}(t)
=\int_0^\infty K_\alpha(a)\Sigma_{\mathrm{SFR}}(t-a)\,da,
\qquad \int_0^\infty K_\alpha(a)\,da=1.
\tag{29}
$$

The overbar denotes this Halpha response average. For transparent analytic calculations we assume

$$
K_\alpha(a)=\frac{1}{\tau_{\mathrm{ion}}}e^{-a/\tau_{\mathrm{ion}}},
\qquad
\tau_{\mathrm{ion}}\frac{d\overline\Sigma_{\mathrm{SFR},\alpha}}{dt}
+\overline\Sigma_{\mathrm{SFR},\alpha}=\Sigma_{\mathrm{SFR}}.
\tag{30}
$$

The exponential kernel and numerical choice $\tau_{\mathrm{ion}}=3$ Myr are our approximations, not a fitted stellar-population model. For the continuous-onset gas solutions we assume a constant pre-onset SFR, so $\overline\Sigma_{\mathrm{SFR},\alpha}(0)=\Sigma_{\mathrm{SFR},0}$. An instantaneous depletion-time change instead requires the actual pre-change stellar history in equation (29); its filtered initial value need not equal the new instantaneous SFR.

To show the algebra, consider a unit-amplitude input $e^{-\gamma t}$ with initial filtered value one. Multiplication of equation (30) by $e^{t/\tau_{\mathrm{ion}}}$ and integration gives

$$
\begin{aligned}
e^{t/\tau_{\mathrm{ion}}}\mathcal F_\alpha(t;\gamma)-1
&=\frac{1}{\tau_{\mathrm{ion}}}
\int_0^t e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)u}\,du\\
&=\frac{e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)t}-1}
{1-\gamma\tau_{\mathrm{ion}}},\\
\mathcal F_\alpha(t;\gamma)
&=\frac{e^{-\gamma t}-\gamma\tau_{\mathrm{ion}}e^{-t/\tau_{\mathrm{ion}}}}
{1-\gamma\tau_{\mathrm{ion}}}.
\end{aligned}
\tag{31}
$$

Here $\mathcal F_\alpha$ is the dimensionless filtered response to one exponential mode, not a new emitting component. At $\gamma\tau_{\mathrm{ion}}=1$, direct integration gives $(1+t/\tau_{\mathrm{ion}})e^{-t/\tau_{\mathrm{ion}}}$.

Linearity then supplies the explicit filtered version of equation (16):

$$
\begin{aligned}
F_\alpha(t)&\equiv
\frac{\overline\Sigma_{\mathrm{SFR},\alpha}(t)}{\Sigma_{\mathrm{SFR},0}}\\
&=\mathcal F_\alpha(t;\gamma_{\mathrm{H_2}})
+\frac{\mathcal F_\alpha(t;\gamma_{\mathrm{HI}})
-\mathcal F_\alpha(t;\gamma_{\mathrm{H_2}})}
{\tau_{\Phi,0}(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})}.
\end{aligned}
\tag{32}
$$

The coefficients add to unity, satisfying the specified initial value. The equal-gas-rate case can be evaluated from equation (29), avoiding a numerically singular difference. On Gyr-scale gas evolution, the 3-Myr smoothing is small, but it prevents an unphysical instantaneous Halpha response to a sharp SFR change.

## 4.2 Allocate the young photons without counting leakage twice

Let $\mathcal Q_{\mathrm{OB}}$ be the young-star hydrogen-ionizing photon production rate per adopted area. Let $f_{\mathrm{HII}}$ be the fraction absorbed by hydrogen in compact HII regions, and $f_{\mathrm{leak}}$ the fraction that escapes those regions and is subsequently absorbed by hydrogen in the diffuse gas assigned to the modeled area. Escaping photons that never ionize this gas and ionizing photons absorbed by dust occupy the remaining budget. In the local approximation,

$$
0\leq f_{\mathrm{HII}}+f_{\mathrm{leak}}\equiv f_{\mathrm{young}}\leq1,
\qquad
\begin{aligned}
\mathcal L_\alpha^{\mathrm{HII}}
&=\epsilon_\alpha f_{\mathrm{HII}}\mathcal Q_{\mathrm{OB}},\\
\mathcal L_\alpha^{\mathrm{leak}}
&=\epsilon_\alpha f_{\mathrm{leak}}\mathcal Q_{\mathrm{OB}}.
\end{aligned}
\tag{33}
$$

$\mathcal L$ denotes luminosity per adopted area, in $\mathrm{erg\,s^{-1}\,kpc^{-2}}$. The photon-to-Halpha energy factor $\epsilon_\alpha$ is derived in the next section. Leaked photons are a redistribution of the young photon budget, not a second independent supply of OB photons. Radiation transport can make leakage nonlocal; Appendix A states the corresponding limitation.

The SFR calibration coefficient $C_\alpha$ is defined by $\Sigma_{\mathrm{SFR}}=C_\alpha\mathcal L_\alpha$ for the adopted fully absorbed young-star calibration. Consequently,

$$
\mathcal L_\alpha^{\mathrm{young}}(t)
\equiv\mathcal L_\alpha^{\mathrm{HII}}+\mathcal L_\alpha^{\mathrm{leak}}
=\frac{f_{\mathrm{young}}}{C_\alpha}
\overline\Sigma_{\mathrm{SFR},\alpha}(t).
\tag{34}
$$

Our numerical value, $C_\alpha=4.9835821\times10^{-42}\ M_\odot\,\mathrm{yr^{-1}}/(\mathrm{erg\,s^{-1}})$, follows the existing MAUVE pipeline convention checked by its saved source fingerprint. It is not assumed universal across IMFs or stellar populations. If $f_{\mathrm{young}}$ is constant, the remaining young-Halpha fraction is exactly $F_\alpha(t)$.

# 5. Deriving and checking the HOLMES Halpha term

## 5.1 The old-star photon supply

Let $\Sigma_*^{\mathrm{old}}$ be the current mass surface density of the old stellar population, including its associated remnants. Let $q_{\mathrm{H,HOLMES}}$ be its hydrogen-ionizing photon production rate per unit **current** mass. Then

$$
\mathcal Q_{\mathrm{HOLMES}}
=q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}}.
\tag{35}
$$

This proportionality is a population-normalization identity for a specified population. It is not a universal assertion that every stellar population produces the same number of ionizing photons per unit mass. A mixture of ages and metallicities requires a correspondingly weighted $q$.

For the numerical benchmark we use $q_{\mathrm{H,HOLMES}}=7\times10^{40}\ \mathrm{s^{-1}}\,M_\odot^{-1}$, the PEGASE value adopted by [Belfiore et al. (2022), section 3.2](#ref-belfiore). Their footnote 5 specifies current stars-plus-remnants mass. Their 10-Gyr solar-metallicity FSPS comparison gives $5\times10^{40}$ in that convention. These are model choices, not measured MAUVE photon rates. We do not silently identify the PEGASE normalization with the IMF of their separate FSPS spectral grid.

## 5.2 Convert absorbed photons into Halpha luminosity

In photoionization equilibrium, the absorbed hydrogen-ionizing photon rate equals the Case-B recombination rate. Let $f_{\mathrm{abs,HOLMES}}$ denote the fraction of the HOLMES photons assigned to the region that actually ionize its hydrogen. For a region with adopted physical area $A_{\mathrm{reg}}$,

$$
f_{\mathrm{abs,HOLMES}}\mathcal Q_{\mathrm{HOLMES}}A_{\mathrm{reg}}
=\alpha_B\int_{V_{\mathrm{reg}}}n_e n_p\,dV.
\tag{36}
$$

Here $n_e$ and $n_p$ are electron and proton densities, $\alpha_B$ is the recombination coefficient excluding direct recombinations to the ground state, and $V_{\mathrm{reg}}$ is the emitting volume. The area and volume must be converted to a consistent unit system. The emitted Halpha luminosity is

$$
L_\alpha^{\mathrm{HOLMES}}
=h_{\mathrm P}\nu_\alpha\alpha_\alpha^{\mathrm{eff}}
\int_{V_{\mathrm{reg}}}n_e n_p\,dV,
\tag{37}
$$

where $h_{\mathrm P}$ is Planck's constant, $\nu_\alpha=c/\lambda_\alpha$ is the line frequency, and $\alpha_\alpha^{\mathrm{eff}}$ counts recombinations that generate Halpha photons. Case-B coefficients depend on temperature and density ([Hummer & Storey 1987](#ref-hs)). Solving equation (36) for the volume integral, inserting it into equation (37), and dividing by $A_{\mathrm{reg}}$ gives

$$
\begin{aligned}
\mathcal L_\alpha^{\mathrm{HOLMES}}
&=h_{\mathrm P}\nu_\alpha
\frac{\alpha_\alpha^{\mathrm{eff}}}{\alpha_B}
f_{\mathrm{abs,HOLMES}}\mathcal Q_{\mathrm{HOLMES}}\\
&=\epsilon_\alpha f_{\mathrm{abs,HOLMES}}
q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}},\\
\epsilon_\alpha&\equiv h_{\mathrm P}\nu_\alpha p_\alpha,
\qquad p_\alpha\equiv\frac{\alpha_\alpha^{\mathrm{eff}}}{\alpha_B}.
\end{aligned}
\tag{38}
$$

This is the verified physical content of equation (34) in the 14 September report. The normalization requires specifying the stellar population, the mass convention, and the absorbed fraction. [Cid Fernandes et al. (2011), equation 2](#ref-cid), explicitly use $L_{\mathrm{Halpha}}=h\nu_\alpha Q_{\mathrm H}/2.206$ when escape and extinction are neglected. We adopt that approximate low-density nebular conversion, $p_\alpha=1/2.206$, giving $\epsilon_\alpha=1.3721\times10^{-12}$ erg per absorbed hydrogen-ionizing photon. Their population mass normalization elsewhere in that paper is formed stellar mass; we use the **Belfiore current-mass convention** for $q$ and $\Sigma_*^{\mathrm{old}}$ here.

## 5.3 When may this term be treated as constant?

Taking the logarithmic derivative of equation (38) exposes all conditions:

$$
\frac{d\ln\mathcal L_\alpha^{\mathrm{HOLMES}}}{dt}
=\frac{d\ln\epsilon_\alpha}{dt}
+\frac{d\ln f_{\mathrm{abs,HOLMES}}}{dt}
+\frac{d\ln q_{\mathrm{H,HOLMES}}}{dt}
+\frac{d\ln\Sigma_*^{\mathrm{old}}}{dt}.
\tag{39}
$$

For a short interval relative to the age of an old population, its mass and specific photon output may vary slowly. We additionally assume nearly fixed absorption and recombination conditions in the retained inner gas. Only with these assumptions do we set $\mathcal L_\alpha^{\mathrm{HOLMES}}(t)\simeq\mathcal L_{\alpha,0}^{\mathrm{HOLMES}}$.

The gas must also be able to absorb the available photons. Equation (36) implies a necessary emission measure. For example, if a maximum available ionized column is specified, its recombination capacity cannot be exceeded merely by increasing $q\Sigma_*^{\mathrm{old}}$. Gas removal can reduce the covering fraction and capacity, causing the HOLMES-powered emission to fall. The recombination time is approximately $\tau_{\mathrm{rec}}=(\alpha_B n_e)^{-1}$ at fixed density; sustained emission over longer intervals requires sustained ionization, not a remnant afterglow.

The constant term is therefore particularly relevant to **retained inner gas**. It is not a reason to predict a permanent Halpha floor in a completely stripped outer region.

# 6. The three-component line budget and its changing weights

## 6.1 Add luminosities first, then construct the ratio

Let $\ell$ denote a particular forbidden line and $B$ its Balmer denominator. We use [N II]6583/Halpha, the [S II]6716+6731 doublet sum/Halpha, and [O III]5007/Hbeta. The proposed budget is

$$
\mathcal L_\ell
=\mathcal L_\ell^{\mathrm{HII}}
+\mathcal L_\ell^{\mathrm{leak}}
+\mathcal L_\ell^{\mathrm{HOLMES}}.
\tag{40}
$$

Within the emitting-region approximation these are positive contributions from distinct gas parcels or effective templates. Shock and AGN contributions are set to zero as a baseline hypothesis requiring observational checks, rather than proved absent. DIG is a description of diffuse gas, which can receive both leaked OB and HOLMES photons; it is not synonymous with one of the two ionizing sources.

Define the component ratio $R_{\ell/B}^j=\mathcal L_\ell^j/\mathcal L_B^j$ and the Balmer weight $w_{j,B}=\mathcal L_B^j/\sum_i\mathcal L_B^i$, with $j\in\{\mathrm{HII,leak,HOLMES}\}$. Substitute $\mathcal L_\ell^j=R_{\ell/B}^j\mathcal L_B^j$:

$$
\begin{aligned}
R_{\ell/B}
&=\frac{\sum_j\mathcal L_\ell^j}{\sum_i\mathcal L_B^i}
=\frac{\sum_j R_{\ell/B}^j\mathcal L_B^j}{\sum_i\mathcal L_B^i}\\
&=\sum_j\left(\frac{\mathcal L_B^j}{\sum_i\mathcal L_B^i}\right)R_{\ell/B}^j\\
&=w_{\mathrm{HII},B}R_{\ell/B}^{\mathrm{HII}}
+w_{\mathrm{leak},B}R_{\ell/B}^{\mathrm{leak}}
+w_{\mathrm{HOLMES},B}R_{\ell/B}^{\mathrm{HOLMES}},\\
\sum_jw_{j,B}&=1.
\end{aligned}
\tag{41}
$$

This is an exact algebraic identity for additive component luminosities. HII/DIG Halpha-flux mixing has a direct antecedent in [Blanc et al. (2009), equations 7--8](#ref-blanc). Their [S II] template is for the single 6717 line; we do not import that numerical template for our doublet sum. Mixture addition occurs in linear ratios, followed by the logarithm used for a BPT plot. Logarithmic BPT coordinates must not be averaged with these weights.

## 6.2 Hbeta ratios require Hbeta weights

Let $\mathcal B_j=\mathcal L_\alpha^j/\mathcal L_\beta^j$ be the component Balmer decrement. From $\mathcal L_\beta^j=\mathcal L_\alpha^j/\mathcal B_j$,

$$
w_{j,\beta}
=\frac{w_{j,\alpha}/\mathcal B_j}
{\sum_iw_{i,\alpha}/\mathcal B_i}.
\tag{42}
$$

Only equal decrements make $w_{j,\beta}=w_{j,\alpha}$. We use intrinsic/de-reddened luminosities and equal $\mathcal B_j=2.86$ in the illustrative curves, consistent with a standard low-density, approximately $10^4$ K Case-B approximation. The empirical [O III]/Hbeta test uses the actual exported Hbeta denominator. A single dust correction applied to a mixture does not necessarily recover the individually corrected source components.

## 6.3 Different HII and leakage ratios do not prevent an analytic solution

For a fixed allocation of young photons, define the compact-HII fraction of **young** Balmer emission,

$$
\eta_B\equiv
\frac{\mathcal L_B^{\mathrm{HII}}}
{\mathcal L_B^{\mathrm{HII}}+\mathcal L_B^{\mathrm{leak}}},
\qquad
R_{\ell/B}^{\mathrm{young}}
\equiv\eta_B R_{\ell/B}^{\mathrm{HII}}
+(1-\eta_B)R_{\ell/B}^{\mathrm{leak}}.
\tag{43}
$$

For constant component decrements and photon allocation, $\eta_B$ is constant. It is $f_{\mathrm{HII}}/f_{\mathrm{young}}$ for Halpha under equation (33). Both young components then fade in the same proportion, leaving their combined ratio constant. It is unnecessary to assume that their individual ratios are equal.

Inserting $w_{\mathrm{HII},B}=(1-w_{\mathrm{HOLMES},B})\eta_B$ and $w_{\mathrm{leak},B}=(1-w_{\mathrm{HOLMES},B})(1-\eta_B)$ into equation (41) gives

$$
R_{\ell/B}(t)
=R_{\ell/B}^{\mathrm{young}}
+\left(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}}\right)
w_{\mathrm{HOLMES},B}(t).
\tag{44}
$$

This is a two-term mathematical reduction of a three-component physical budget. Neither compact HII nor leakage has been discarded.

## 6.4 The source of the changing weight

With constant HOLMES luminosity and young luminosity proportional to $F_\alpha(t)$,

$$
w_{\mathrm{HOLMES},B}(t)
=\frac{\mathcal L_{B,0}^{\mathrm{HOLMES}}}
{\mathcal L_{B,0}^{\mathrm{young}}F_\alpha(t)+\mathcal L_{B,0}^{\mathrm{HOLMES}}}
=\frac{w_{\mathrm{HOLMES},B}(0)}
{w_{\mathrm{HOLMES},B}(0)+[1-w_{\mathrm{HOLMES},B}(0)]F_\alpha(t)}.
\tag{45}
$$

The second expression follows by dividing numerator and denominator by the initial total Balmer luminosity. Differentiating the first expression, using a constant numerator, gives

$$
\begin{aligned}
\frac{dw_{\mathrm{HOLMES},B}}{dt}
&=-\frac{\mathcal L_{B,0}^{\mathrm{HOLMES}}
\mathcal L_{B,0}^{\mathrm{young}}}
{[\mathcal L_{B,0}^{\mathrm{young}}F_\alpha+
\mathcal L_{B,0}^{\mathrm{HOLMES}}]^2}\frac{dF_\alpha}{dt}\\
&=-w_{\mathrm{HOLMES},B}(1-w_{\mathrm{HOLMES},B})
\frac{d\ln F_\alpha}{dt}.
\end{aligned}
\tag{46}
$$

Therefore the HOLMES Balmer weight rises when the young emission declines. Differentiating equation (44) gives the conditional prediction

$$
\frac{dR_{\ell/B}}{dt}
=\left(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}}\right)
\frac{dw_{\mathrm{HOLMES},B}}{dt}>0
\tag{47}
$$

when $dF_\alpha/dt<0$ and $R_{\ell/B}^{\mathrm{HOLMES}}>R_{\ell/B}^{\mathrm{young}}$. The ratio enhancement is thus a consequence of unequal source evolution, not fading alone.

## 6.5 Why the forbidden line may still fade while its ratio increases

The total forbidden-line luminosity is

$$
\mathcal L_\ell(t)
=R_{\ell/B}^{\mathrm{young}}\mathcal L_{B,0}^{\mathrm{young}}F_\alpha(t)
+R_{\ell/B}^{\mathrm{HOLMES}}\mathcal L_{B,0}^{\mathrm{HOLMES}}.
\tag{48}
$$

For positive ratios, its derivative is $R_{\ell/B}^{\mathrm{young}}\mathcal L_{B,0}^{\mathrm{young}}\,dF_\alpha/dt<0$. Its fractional fading is smaller in magnitude than the Balmer fractional fading if its HOLMES-to-young ratio contrast is larger. Eliminating $F_\alpha$ between the total Balmer luminosity and equation (48) gives

$$
\begin{aligned}
\mathcal L_\ell
&=R_{\ell/B}^{\mathrm{young}}\mathcal L_B
+(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\mathcal L_{B,0}^{\mathrm{HOLMES}},\\
R_{\ell/B}
&=R_{\ell/B}^{\mathrm{young}}
+\frac{(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\mathcal L_{B,0}^{\mathrm{HOLMES}}}{\mathcal L_B},\\
\frac{d\ln R_{\ell/B}}{d\ln\mathcal L_B}
&=-\frac{R_{\ell/B}-R_{\ell/B}^{\mathrm{young}}}{R_{\ell/B}}.
\end{aligned}
\tag{49}
$$

For the higher-ratio HOLMES branch, the last slope lies between $-1$ and $0$. This supplies a directly testable ratio--surface-brightness relation at fixed old stellar density and fixed templates. It also clarifies what differential fading means here: the distinct components fade differently; each component's internal line ratios are held constant.

The HOLMES fraction of a forbidden line can exceed its Balmer fraction. Explicitly,

$$
\frac{\mathcal L_\ell^{\mathrm{HOLMES}}}{\mathcal L_\ell}
=\frac{w_{\mathrm{HOLMES},B}R_{\ell/B}^{\mathrm{HOLMES}}}
{(1-w_{\mathrm{HOLMES},B})R_{\ell/B}^{\mathrm{young}}
+w_{\mathrm{HOLMES},B}R_{\ell/B}^{\mathrm{HOLMES}}}.
\tag{50}
$$

At a Balmer fraction of 0.02, a forbidden-line ratio contrast greater than 49 would make the HOLMES component dominate that line. This is our algebraic illustration, not an extraction of a contrast from Belfiore. In a shared gas volume, even this source-by-source attribution of forbidden-line cooling is not unique; the mixed incident spectrum changes the gas state jointly.

# 7. The complete connection and the meaning of SF/NSF

Combining equations (32), (34), (38), and (44) provides the closed prediction for a Balmer-alpha ratio:

$$
\boxed{
\begin{aligned}
\mathcal L_\alpha(t)
&=\frac{f_{\mathrm{young}}\Sigma_{\mathrm{SFR},0}}{C_\alpha}F_\alpha(t)
+\epsilon_\alpha f_{\mathrm{abs,HOLMES}}
q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}},\\
R_{\ell/\alpha}(t)
&=R_{\ell/\alpha}^{\mathrm{young}}
+\left(R_{\ell/\alpha}^{\mathrm{HOLMES}}-R_{\ell/\alpha}^{\mathrm{young}}\right)
\frac{\epsilon_\alpha f_{\mathrm{abs,HOLMES}}
q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}}}
{\mathcal L_\alpha(t)}.
\end{aligned}}
\tag{51}
$$

The gas rates enter through $F_\alpha(t)$; the young luminosity normalization enters through the SFR calibration; and the constant component is constrained by an old-star photon budget. The line-ratio model and gas model are therefore connected explicitly. Halpha is the observable link, but the model additionally requires gas retention and old-star ionization.

If there is an initial SFR enhancement, the young contribution can initially increase, reducing the HOLMES fraction and moving ratios toward the young-component value. After the peak and the short ionizing response delay, the sense reverses. This is a conditional prediction of the same equations, without changing the spectral templates.

We should apply the three-component budget in both SF and NSF regions. Ignoring HOLMES in SF is acceptable only when it is negligible for the quantity in question: a small Balmer contribution alone does not guarantee a small [O III] contribution. Similarly, NSF selection does not prove $\mathcal L_\alpha^{\mathrm{HII}}=0$. A faint compact HII region can be unresolved or can fail an EW, dispersion, or quality threshold.

In the source analysis, SF requires the finite HII-selected SFR product, Halpha EW greater than 6 Angstrom, and intrinsic Halpha dispersion below $45\ \mathrm{km\,s^{-1}}$. ND is the joint Balmer non-detection category; NSF comprises the remaining Balmer-detected regions outside the SF selection. Therefore NSF is not identical to LIER, DIG, or HOLMES domination. The fraction of spatial area classified as NSF is also not a Balmer luminosity weight.

A model of NSF occupancy must forward-model the continuum, line fluxes, line widths, noise, and masks before applying the observational rules. That classification calculation is outside this revision. We predict the evolution of line emission within retained-gas regions, not an already validated SF-to-NSF transition rate.

# 8. Numerical predictions and comparison with MAUVE scales

## 8.1 Which observations are used?

We reuse the five-line, common-support export from the 14 September analysis. The example bin is centered on $\log_{10}[\Sigma_*/(M_\odot\,\mathrm{kpc}^{-2})]=8.625$. Within each category and bin, line surface luminosities are averaged within each contributing galaxy and then combined with equal galaxy weighting. Ratios below are ratios of those linear line means, not means of spaxel ratios. The source imposes the usable-disc selection and post-fit continuum S/N greater than 25. The present calculation regenerates these anchors from the export after checking its recorded notebook/pipeline fingerprints; it does not re-execute the maps or alter the selection.

**Table 1. NSF observational anchors.** Surface luminosities are in $\mathrm{erg\,s^{-1}\,kpc^{-2}}$. These cross-sectional galaxy samples are not a temporal sequence of identical patches.

| Stage | Galaxies on common support | $\mathcal L_\alpha$ | $\mathcal L_\beta$ | N2 | S2 | O3 |
|:--|--:|--:|--:|--:|--:|--:|
| Pre-peak | 6 | $3.558\times10^{39}$ | $1.244\times10^{39}$ | 0.3731 | 0.3773 | 0.8978 |
| Close-to-peak | 5 | $2.077\times10^{39}$ | $7.266\times10^{38}$ | 0.4367 | 0.3857 | 0.4678 |
| Post-peak | 13 | $1.111\times10^{39}$ | $3.892\times10^{38}$ | 0.5663 | 0.4168 | 0.7528 |

Here N2=[N II]6583/Halpha, S2=[S II] sum/Halpha, and O3=[O III]5007/Hbeta. The table shows increasing N2 and S2 but does **not** show a monotonic O3 increase. A model claiming that all BPT ratios rise must confront this distinction. These are descriptive mean-value checks; no significance is assigned without galaxy-level uncertainty propagation and matched support.

Pre-peak Virgo galaxies can serve as an operational, relatively less-processed reference. Calling them field-like is a modeling approximation, not an established environmental equivalence to an external field sample.

## 8.2 Normalize the HOLMES term before selecting line ratios

Take the bin center as a representative $\Sigma_*=10^{8.625}=4.217\times10^8\ M_\odot\,\mathrm{kpc}^{-2}$. For a deliberately generous **local benchmark**, assume all this mass is old, its mass convention matches current stars plus remnants, and all available HOLMES photons ionize the retained hydrogen. Equation (38) gives

$$
\begin{aligned}
\mathcal L_{\alpha}^{\mathrm{HOLMES}}
&=(1.3721\times10^{-12})(7\times10^{40})(4.217\times10^8)\\
&=4.050\times10^{37}\ \mathrm{erg\,s^{-1}\,kpc^{-2}}.
\end{aligned}
\tag{52}
$$

The corresponding pre-peak and post-peak NSF Halpha weights are 0.01138 and 0.03647. With the FSPS comparison normalization they would be smaller by $5/7$. A younger mass fraction or incomplete absorption also reduces them. The bin-center approximation, population uncertainty, mass-convention compatibility, and nonlocal photon transport prevent treating this as an absolute universal ceiling. It is the maximum **within the specified local fiducial population model**.

This numerical result is essential: a spatially broad HOLMES component may be real yet contribute little Halpha in a bright NSF sample. The relevant variable is $\mathcal L_\alpha/\Sigma_*^{\mathrm{old}}$, not Halpha brightness or stellar density separately.

## 8.3 Gas evolution, delayed decline, and a modest enhancement

Subtract equation (52) from the pre-peak NSF Halpha anchor, assume $f_{\mathrm{young}}=1$, and apply $C_\alpha$. This gives an illustrative initial young SFR of $0.01753\ M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}$. This is a model decomposition of a luminosity scale; it is **not an independent SFR measurement in NSF**.

**Table 2. Gas-model parameters.** The gas columns and times are illustrative, not CO/HI fits.

| Quantity | Value | Reason/status |
|:--|--:|:--|
| $\tau_{\mathrm{dep}}$ | 2 Gyr | Normal-disc scale; restricted constant-efficiency assumption |
| $R$, $\lambda$ | 0.4, 0 | Illustrative recycling and zero feedback outflow |
| $\gamma_{\mathrm{H_2}}$ | $0.3\ \mathrm{Gyr^{-1}}$ | Derived from equation (5) |
| $\Sigma_{\mathrm{H_2},0}$ | $35.055\ M_\odot\,\mathrm{pc^{-2}}$ | Derived from adopted SFR and depletion time |
| $\Sigma_{\mathrm{HI},0}$ | $10\ M_\odot\,\mathrm{pc^{-2}}$ | Assumed local atomic phase column, including helium |
| $\tau_{\mathrm{conv,pre}}$ | 0.95088 Gyr | Chosen for initial molecular balance |
| $\gamma_{\mathrm{strip}}$ | $3\ \mathrm{Gyr^{-1}}$ | Assumed atomic loss rate; not measured from an orbit |
| $\tau_{\mathrm{ion}}$ | 0.003 Gyr | Approximate young-photon response |

For the balanced branch, $\gamma_{\mathrm{HI}}=4.05166\ \mathrm{Gyr^{-1}}$. Equations (19) and (32) give $F_{\mathrm{SFR}}(1\ \mathrm{Gyr})=0.79867$ and $F_\alpha=0.79937$. The SFR falls by about 0.098 dex in one Gyr, despite atomic gas declining much faster. The molecular reservoir buffers the loss of its supply.

The closed comparison with $\gamma_{\mathrm{strip}}=0$ gives $F_{\mathrm{SFR}}=0.89706$ at the same time. This is important: shutting off external supply and converting a finite atomic reservoir causes some decline even without RPS. The additional change attributable to stripping in this particular comparison is a factor $0.79867/0.89706=0.8903$, or about $-0.050$ dex. A maintained-supply control is a different boundary condition, discussed in Appendix A.

For the enhancement example, reduce $\tau_{\mathrm{conv}}$ by a factor of four to 0.23772 Gyr while keeping the initial columns, depletion time, and stripping rate unchanged. Then $\tau_{\Phi,0}=0.83333$ Gyr and the initial supply is four times the net consumption rate. Equation (24) gives

$$
t_{\mathrm{peak}}=0.18368\ \mathrm{Gyr},
\qquad
F_{\mathrm{SFR}}(t_{\mathrm{peak}})=1.06458.
\tag{53}
$$

The enhancement is only 0.0272 dex. By one Gyr, SFR has declined to 0.86940 of its initial value. This explicitly demonstrates enhancement followed by decline, but also shows why an arbitrarily large burst cannot be obtained just by speeding conversion. The mass-budget bound is 1.2853, or 0.1090 dex, and ongoing stripping/consumption keeps the actual peak well below that bound. A larger initial atomic reservoir, changed depletion time, or external/compressive transport would be required for a larger response.

![Figure 1. The analytical SFR response and atomic reservoir for an initially balanced region, its closed no-RPS comparison, and a fourfold increase in the conversion rate. All curves use the same initial gas columns. The compression branch rises briefly, then declines; the purely stripped balanced branch declines after a zero initial slope.](assets/20261001_HI_stripping_HOLMES/figure01_SFR_response.png)

## 8.4 How much fading is required for the HOLMES fraction to matter?

Solving equation (45) for the remaining young fraction at a target HOLMES weight $w_{\mathrm{target},B}$ gives

$$
F_\alpha
=\frac{w_{\mathrm{HOLMES},B}(0)[1-w_{\mathrm{target},B}]}
{w_{\mathrm{target},B}[1-w_{\mathrm{HOLMES},B}(0)]}.
\tag{54}
$$

For the bright pre-peak NSF anchor, reaching a 10% HOLMES Halpha weight requires $F_\alpha\simeq0.104$, and reaching 50% requires $F_\alpha\simeq0.0115$. A large change in the mixture therefore requires substantial fading when the initial old-star contribution is only about 1%.

For comparison, an explicitly hypothetical faint patch with the same old stellar density but initial $\mathcal L_\alpha=2\times10^{38}$ has a HOLMES weight of 0.203. Reducing its young emission by a factor of ten raises that weight to 0.717. With an effective young N2 ratio of 0.35 and a HOLMES N2 ratio of 1.5, its total N2 rises from approximately 0.583 to 1.175 while [N II] itself becomes fainter. The spectra here are illustrative endpoints, not measurements of individual MAUVE components.

For Figure 2, the effective young N2 value 0.35 is realized by a 70:30 young Halpha allocation with $R_{\mathrm{N2}}^{\mathrm{HII}}=0.30$ and $R_{\mathrm{N2}}^{\mathrm{leak}}=0.4667$. Thus equality of the HII and leakage ratios is demonstrably unnecessary. The figure varies the remaining young emission directly; it does not claim that a hundredfold fading is reached within the plotted gas model's two-Gyr interval or while every constant-population assumption remains valid.

![Figure 2. Three-component fading at the bright NSF scale and at a hypothetical faint scale, both with the same fiducial HOLMES luminosity. Left: the HOLMES Halpha weight increases as young emission fades. Middle: N2 increases as total Halpha decreases. Right: [N II] still fades. The independent variable decreases from left to right. These are illustrative spectral mixtures, not fitted component measurements.](assets/20261001_HI_stripping_HOLMES/figure02_fading_and_ratio.png)

## 8.5 Test the amplitude against the NSF anchors

For each diagnostic, choose an illustrative fixed HOLMES ratio, normalize the effective young ratio to the pre-peak NSF mean, and then predict the post-peak ratio using its **observed Balmer luminosity** and the fixed photon-budget term. The initial normalization follows by rearranging equation (44):

$$
R_{\ell/B}^{\mathrm{young}}
=\frac{R_{\ell/B,0}^{\mathrm{obs}}
-w_{\mathrm{HOLMES},B}(0)R_{\ell/B}^{\mathrm{HOLMES}}}
{1-w_{\mathrm{HOLMES},B}(0)}.
\tag{55}
$$

This is one-point calibration, not an independent prediction of the initial spectrum. The endmember ratios 1.5, 1.0, and 3.0 below are sensitivity choices, not values extracted from Belfiore, a CLOUDY grid, or the MAUVE spectra. We subsequently invert the equations so the conclusion does not rest only on those choices.

**Table 3. Conditional post-peak predictions.** Residuals are $\log_{10}(R^{\mathrm{pred}}/R^{\mathrm{obs}})$; no uncertainty or formal fit significance is assigned.

| Ratio | Assumed HOLMES | Calibrated young | Predicted post-peak | Observed post-peak | Residual (dex) |
|:--|--:|--:|--:|--:|--:|
| N2 | 1.5 | 0.3601 | 0.4017 | 0.5663 | -0.1492 |
| S2 | 1.0 | 0.3701 | 0.3931 | 0.4168 | -0.0255 |
| O3 | 3.0 | 0.8736 | 0.9510 | 0.7528 | +0.1015 |

The N2 prediction has the desired sign but insufficient amplitude. The O3 prediction has the opposite sign to the pre-to-post change of these exported means. This is not a proof against HOLMES ionization. It is a failure of this particular common-population, fixed-spectrum, additive interpretation of the selected mean points.

![Figure 3. Conditional fixed-HOLMES predictions compared with the exported NSF mean ratios. The initial young ratio is calibrated to the pre-peak point. Post-peak Balmer luminosities are supplied from the data. The differences therefore test the constant-spectrum/photon-budget assumptions rather than the gas clock. Lines connect different galaxy samples for comparison and do not imply observed temporal tracks.](assets/20261001_HI_stripping_HOLMES/figure03_MAUVE_budget_test.png)

## 8.6 Invert the required spectrum or photon budget

Write equation (49) at two observed points indexed by 0 and 1. Subtracting the ratios removes the effective young ratio:

$$
\begin{aligned}
R_{\ell/B,1}-R_{\ell/B,0}
&=(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\mathcal L_B^{\mathrm{HOLMES}}
\left(\frac{1}{\mathcal L_{B,1}}-\frac{1}{\mathcal L_{B,0}}\right),\\
R_{\ell/B}^{\mathrm{young}}
&=\frac{R_{\ell/B,0}\mathcal L_{B,0}-R_{\ell/B,1}\mathcal L_{B,1}}
{\mathcal L_{B,0}-\mathcal L_{B,1}}.
\end{aligned}
\tag{56}
$$

The second line follows by subtracting the **line luminosities**, whose constant intercept cancels. Substituting it into either point gives the required HOLMES ratio at a specified old-star Balmer luminosity:

$$
R_{\ell/B}^{\mathrm{HOLMES,required}}
=R_{\ell/B}^{\mathrm{young}}
+\frac{R_{\ell/B,1}-R_{\ell/B,0}}
{\mathcal L_B^{\mathrm{HOLMES}}(\mathcal L_{B,1}^{-1}-\mathcal L_{B,0}^{-1})}.
\tag{57}
$$

At the fiducial HOLMES budget, the two NSF endpoints require N2$_{\mathrm{HOLMES}}=7.99$, S2$_{\mathrm{HOLMES}}=1.94$, and O3$_{\mathrm{HOLMES}}=-4.84$. The negative O3 value is unphysical for a positive emitting component. The extreme required N2 value is not an adopted spectrum; it quantifies the burden placed on the simple mixture and should be checked against a self-consistent photoionization grid rather than accepted as a free fitting coefficient.

Alternatively, retaining an assumed HOLMES ratio gives

$$
\mathcal L_B^{\mathrm{HOLMES,required}}
=\frac{R_{\ell/B,1}-R_{\ell/B,0}}
{(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
(\mathcal L_{B,1}^{-1}-\mathcal L_{B,0}^{-1})}.
\tag{58}
$$

For N2$_{\mathrm{HOLMES}}=1.5$, the required Halpha luminosity is $2.57\times10^{38}$, about 6.34 times the fiducial local budget. For S2$_{\mathrm{HOLMES}}=1.0$, it is $9.97\times10^{37}$, about 2.46 times that budget. These two required amplitudes also disagree with one another. Enlarging the old-star component arbitrarily is therefore not a self-consistent simultaneous solution.

This inversion is a conditional mean-value consistency test. The samples have different galaxies, imperfectly matched gas conditions, and uncertainties not propagated here. It identifies which assumptions need testing; it is not a statistical rejection of the mechanism across all NSF spaxels.

## 8.7 Does the gas model produce the observed amount of fading fast enough?

With the fixed HOLMES term removed, the young-Halpha fraction needed to connect the pre-peak and post-peak NSF anchors is

$$
F_{\alpha,\mathrm{required}}
=\frac{\mathcal L_{\alpha,1}-\mathcal L_\alpha^{\mathrm{HOLMES}}}
{\mathcal L_{\alpha,0}-\mathcal L_\alpha^{\mathrm{HOLMES}}}
=0.30426.
\tag{59}
$$

For a monotonically declining history with constant prehistory, the young-star response averages earlier, larger SFR values, so $F_\alpha(t)\geq F_{\mathrm{SFR}}(t)\geq e^{-\gamma_{\mathrm{H_2}}t}$. With $\gamma_{\mathrm{H_2}}=0.3\ \mathrm{Gyr^{-1}}$, reaching equation (59) therefore requires at least 3.97 Gyr; solving the chosen balanced filtered model gives 4.23 Gyr. This long extrapolation lies outside the intended short-interval, fixed-population interpretation and is not an inferred infall time.

If an independent orbital constraint required the change to occur within one Gyr, the no-supply bound with $R=0.4$ and $\lambda=0$ would require $\tau_{\mathrm{dep}}\lesssim0.504$ Gyr, with positive supply making the requirement stricter. A changing depletion time, changed photon absorption, transport, direct molecular loss, or different initial conditions could alter the conclusion, but each is an additional hypothesis. HI-only stripping with a fixed normal-disc depletion time is not automatically a rapid-quenching model.

# 9. What the model explains, and what remains unresolved

The revised analytical model provides a consistent route from atomic loss to a slowly declining molecular reservoir and young ionizing emission. An initial enhancement is possible when molecular supply temporarily exceeds consumption or the depletion time decreases; direct HI stripping alone does not create that enhancement from a balanced initial state. The exact peak time and gas-mass ceiling make this statement quantitative.

Adding HOLMES supplies the independent slowly varying source needed for changing Balmer weights. The ratio rise follows rigorously when the relevant HOLMES spectrum has a larger ratio and its luminosity fades more slowly. It does not require equal compact-HII and leaked-OB ratios, nor does it require declaring that compact HII vanishes in NSF. The same equations allow forbidden lines to fade less strongly than Balmer lines.

The numerical tests also set two substantive limits. First, the HI-only gas model can decline too slowly at a fixed two-Gyr depletion time. Second, the fiducial local HOLMES photon budget is too small to reproduce the selected bright NSF N2 change with the illustrated additive fixed spectra. The selected O3 means do not follow the assumed increasing-ratio branch. These limitations should remain visible in any interpretation; the model is a testable baseline, not a completed explanation of every observed stage trend.

The most direct next empirical test is to fit line luminosities against Balmer luminosity at matched old stellar density, rather than fitting ratios alone. Equation (49) predicts a slope equal to the young ratio and an intercept tied to $q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}}$. All five lines should share a physically consistent source normalization and component decrements. Stellar-population constraints, spatial transport, component spectra, and observational covariance should be propagated jointly, with galaxies treated as independent systems. No such full fit is claimed here.

For a next physical model, shared-gas photoionization with a declining OB field and retained HOLMES field is particularly relevant. It allows the evolving radiation hardness and ionization parameter to change the intrinsic gas spectrum. That is a more faithful implementation of mixed illumination, but it requires a photoionization grid; arbitrary free functions for each ratio would remove the model's explanatory power.

# Appendix A. What changes if the simplifying quantities evolve?

## A.1 Time-dependent gas coefficients and external atomic supply

Let $\Sigma_{\mathrm{in}}(t)$ be external supply to the atomic reservoir, distinct from the internal molecular supply $\Sigma_\Phi$. Define the accumulated response $G_i(t)=\int_0^t\gamma_i(u)\,du$ for $i=\mathrm{HI,H_2}$. Multiplying each linear balance by its integrating factor gives

$$
\begin{aligned}
\Sigma_{\mathrm{HI}}(t)
&=e^{-G_{\mathrm{HI}}(t)}
\left[\Sigma_{\mathrm{HI},0}+\int_0^t e^{G_{\mathrm{HI}}(u)}\Sigma_{\mathrm{in}}(u)\,du\right],\\
\Sigma_{\mathrm{H_2}}(t)
&=e^{-G_{\mathrm{H_2}}(t)}
\left[\Sigma_{\mathrm{H_2},0}
+\int_0^t e^{G_{\mathrm{H_2}}(u)}
\frac{\Sigma_{\mathrm{HI}}(u)}{\tau_{\mathrm{conv}}(u)}\,du\right],\\
\Sigma_{\mathrm{SFR}}(t)&=\frac{\Sigma_{\mathrm{H_2}}(t)}{\tau_{\mathrm{dep}}(t)}.
\end{aligned}
\tag{60}
$$

The product-rule steps are identical to equations (10)--(11), but the exponent contains an integral of the rate. This quadrature solution remains exact for specified time-dependent coefficients, provided the system stays linear. The simple two-exponential form does not.

With constant positive external supply, the steady columns satisfy $\Sigma_{\mathrm{HI},\infty}=\Sigma_{\mathrm{in}}/\gamma_{\mathrm{HI}}$ and $\Sigma_{\mathrm{H_2},\infty}=\Sigma_{\mathrm{in}}/(\tau_{\mathrm{conv}}\gamma_{\mathrm{HI}}\gamma_{\mathrm{H_2}})$. At fixed external supply and conversion time, the asymptotic SFR relative to the otherwise identical unstripped equilibrium is $[1+\gamma_{\mathrm{strip}}\tau_{\mathrm{conv}}]^{-1}$. The model then approaches a nonzero level rather than necessarily quenching completely.

## A.2 Evolving HOLMES luminosity and component spectra

Without setting the old-star luminosity constant, differentiation of its weight gives

$$
\frac{dw_{\mathrm{HOLMES},B}}{dt}
=w_{\mathrm{HOLMES},B}(1-w_{\mathrm{HOLMES},B})
\left[
\frac{d\ln\mathcal L_B^{\mathrm{HOLMES}}}{dt}
-\frac{d\ln\mathcal L_B^{\mathrm{young}}}{dt}
\right].
\tag{61}
$$

Thus the relevant condition is **slower fractional fading** of HOLMES-powered emission. Strict constancy is sufficient but unnecessary. If absorbing gas disappears and the HOLMES term fades faster, the proposed weight increase may fail.

If the spectra also evolve, the full derivative is

$$
\begin{aligned}
\frac{dR_{\ell/B}}{dt}
&=(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\frac{dw_{\mathrm{HOLMES},B}}{dt}\\
&\quad +(1-w_{\mathrm{HOLMES},B})\frac{dR_{\ell/B}^{\mathrm{young}}}{dt}
+w_{\mathrm{HOLMES},B}\frac{dR_{\ell/B}^{\mathrm{HOLMES}}}{dt}.
\end{aligned}
\tag{62}
$$

Furthermore,

$$
\begin{aligned}
\frac{dR_{\ell/B}^{\mathrm{young}}}{dt}
&=\eta_B\frac{dR_{\ell/B}^{\mathrm{HII}}}{dt}
+(1-\eta_B)\frac{dR_{\ell/B}^{\mathrm{leak}}}{dt}\\
&\quad +(R_{\ell/B}^{\mathrm{HII}}-R_{\ell/B}^{\mathrm{leak}})
\frac{d\eta_B}{dt}.
\end{aligned}
\tag{63}
$$

These equations separate three effects: changing old-versus-young luminosity weights, changing intrinsic spectra, and changing compact-versus-diffuse allocation of young photons. The fixed-spectrum baseline is useful because it isolates the first. If it fails, equations (62)--(63) show which additional physical terms might matter, but do not determine them from one ratio.

## A.3 Nonlocal leakage and mixed illumination

A local SFR cannot necessarily predict the leakage-powered emission at that same position. A transport kernel $\mathcal T_{\mathrm{leak}}(\boldsymbol{x},\boldsymbol{x}')$, with units of inverse adopted area, can map emitted photons at $\boldsymbol{x}'$ to absorbed photons per unit area at $\boldsymbol{x}$:

$$
\mathcal Q_{\mathrm{abs,leak}}(\boldsymbol{x},t)
=\int \mathcal T_{\mathrm{leak}}(\boldsymbol{x},\boldsymbol{x}')
\mathcal Q_{\mathrm{OB}}(\boldsymbol{x}',t)\,dA'.
\tag{64}
$$

Its area integral must respect the available escaped-photon fraction. Neighbouring young populations can maintain diffuse emission while a local patch fades, breaking the assumption of a common local $F_\alpha$. Similarly, old-star photons can propagate away from their birth positions. The local photon benchmark in section 8 excludes these transfers.

In a gas parcel illuminated simultaneously by OB stars and HOLMES, absorbed photon rates may be budgeted by source, but forbidden-line emissivities depend on the jointly determined ionic fractions and temperature. There is no unique decomposition of each collisionally excited line into an OB part and a HOLMES part independent of that solution. The fundamental prediction then takes the form

$$
R_{\ell/B}=\mathscr R_{\ell/B}
\left(\mathrm{SED}_{\mathrm{OB}}+\mathrm{SED}_{\mathrm{HOLMES}},
U,Z,\mathrm{N/O},n_{\mathrm H},\ldots\right),
\qquad U\equiv\frac{\Phi_{\mathrm H}}{n_{\mathrm H}c}.
\tag{65}
$$

$\mathscr R$ denotes the result of a photoionization calculation, not a fitted analytic function in this report; the SEDs include their radiation normalizations. $\Phi_{\mathrm H}$ is incident ionizing photon flux, $n_{\mathrm H}$ hydrogen density, and $U$ the dimensionless ionization parameter. Metallicity $Z$, abundance ratio N/O, gas density, geometry, and dust also affect the result. Holding every intrinsic line ratio constant is a controlled approximation to be tested, especially as the radiation field fades. It should not be defended as an exact consequence of a fixed stellar spectrum.

# Appendix B. Complete notation and parameter ledger

**Table B1. Gas and temporal quantities.** Initial values carry subscript 0; pre/post on conversion times refer to the imposed change at model onset, not directly to observed infall-stage membership.

| Symbol | Definition and units |
|:--|:--|
| $\boldsymbol{x}$, $t$, $u$, $a$ | Projected location; elapsed model time; dummy time variable; stellar age. Times are consistently converted between Gyr and yr. |
| $\Sigma_{\mathrm{HI}}$, $\Sigma_{\mathrm{H_2}}$ | Atomic/molecular phase mass per adopted area; include associated helium in the numerical model; quoted in $M_\odot\,\mathrm{pc^{-2}}$. |
| $\Sigma_{\mathrm{SFR}}$ | Total star formation rate per adopted area, usually quoted in $M_\odot\,\mathrm{yr^{-1}\,kpc^{-2}}$. |
| $\Sigma_\Phi$, $\Sigma_{\mathrm{in}}$ | Internal atomic-to-molecular supply rate and external supply to HI, respectively; mass per area per time. |
| $\tau_{\mathrm{dep}}$ | Molecular mass divided by SFR; 2 Gyr in the baseline. Net consumption takes longer when $R>0$ and $\lambda=0$. |
| $\tau_{\mathrm{conv}}$, $\tau_\Phi$ | Atomic-to-molecular transfer time in the linear closure; molecular replenishment time $\Sigma_{\mathrm{H_2}}/\Sigma_\Phi$. They are different quantities. |
| $R$, $\lambda$ | Prompt stellar mass return fraction and feedback mass-loading factor; dimensionless. They do not denote line ratios or wavelength here. |
| $\gamma_{\mathrm{strip}}$ | Direct atomic stripping coefficient, $\mathrm{Gyr^{-1}}$; no direct H2 stripping term. |
| $\gamma_{\mathrm{HI}}$, $\gamma_{\mathrm{H_2}}$ | Total HI fractional loss including conversion; net H2 consumption rate, both $\mathrm{Gyr^{-1}}$. |
| $F_{\mathrm{SFR}}$, $t_{\mathrm{peak}}$ | Remaining SFR relative to its initial value; time of the positive-time enhancement maximum. |
| $G_i(t)$ | Dimensionless accumulated rate $\int_0^t\gamma_i(u)du$, used only in the evolving-coefficient appendix. |

**Table B2. Ionizing populations and line emission.** A calligraphic luminosity is per adopted area; ordinary $L$ is an integrated luminosity.

| Symbol | Definition and units |
|:--|:--|
| $K_\alpha(a)$, $\tau_{\mathrm{ion}}$ | Normalized young-ionizing response kernel, inverse time; its exponential response time, 3 Myr in the examples. |
| $\overline\Sigma_{\mathrm{SFR},\alpha}$, $F_\alpha$ | SFR convolved with that kernel; its ratio to initial SFR. |
| $\mathcal F_\alpha(t;\gamma)$ | Dimensionless filtered response to one exponential input with rate $\gamma$. |
| $\mathcal Q_{\mathrm{OB}}$, $\mathcal Q_{\mathrm{HOLMES}}$ | Ionizing photon production per adopted area, $\mathrm{s^{-1}\,kpc^{-2}}$. |
| $f_{\mathrm{HII}}$, $f_{\mathrm{leak}}$, $f_{\mathrm{young}}$ | Young photon fractions absorbed by hydrogen in compact HII, then diffuse gas, and their sum. All are dimensionless. |
| $\Sigma_*^{\mathrm{old}}$, $\Sigma_*$ | Current old-population stellar-plus-remnant mass per area; total stellar mass surface-density control variable. Their equality is only a numerical benchmark. |
| $q_{\mathrm{H,HOLMES}}$ | Old-population ionizing photons per second per current stellar-plus-remnant solar mass; fiducial $7\times10^{40}$. |
| $f_{\mathrm{abs,HOLMES}}$ | Fraction of assigned HOLMES photons that ionize hydrogen in the modeled region; set to 1 for the generous local benchmark. |
| $\mathcal L_\ell^j$, $L_\ell^j$ | Component luminosity surface density, $\mathrm{erg\,s^{-1}\,kpc^{-2}}$, or total luminosity, $\mathrm{erg\,s^{-1}}$, of line $\ell$. |
| $j$, HII, leak, HOLMES, young | Component index; compact HII emission; diffuse gas powered by leaked OB photons; old-star powered emission; sum of HII and leak. |
| $C_\alpha$ | Fully absorbed young Halpha-to-SFR conversion coefficient, $4.9835821\times10^{-42}\ M_\odot\,\mathrm{yr^{-1}}/(\mathrm{erg\,s^{-1}})$. |
| $\alpha$, $\beta$, $B$, $\ell$ | Halpha, Hbeta, the relevant Balmer denominator, and a chosen forbidden emission line. |

**Table B3. Recombination and spectral quantities.**

| Symbol | Definition and units |
|:--|:--|
| $h_{\mathrm P}$, $c$, $\nu_\alpha$, $\lambda_\alpha$ | Planck's constant; speed of light; Halpha frequency; Halpha wavelength, 6562.8 Angstrom in the calculation. |
| $\alpha_B$, $\alpha_\alpha^{\mathrm{eff}}$ | Total Case-B and effective Halpha recombination coefficients, $\mathrm{cm^3\,s^{-1}}$. |
| $p_\alpha$, $\epsilon_\alpha$ | Halpha photons per absorbed hydrogen-ionizing photon; Halpha energy per such absorbed photon. Numerical values 0.45331 and $1.3721\times10^{-12}$ erg. |
| $n_e$, $n_p$, $n_{\mathrm H}$ | Electron, proton, and total hydrogen number densities, $\mathrm{cm^{-3}}$. |
| $A_{\mathrm{reg}}$, $V_{\mathrm{reg}}$, $\tau_{\mathrm{rec}}$ | Adopted physical area, emitting volume, and approximate recombination time. |
| $R_{\ell/B}^j$, $R_{\ell/B}$ | Component and total linear forbidden-line-to-Balmer ratios; dimensionless. |
| $w_{j,B}$, $w_{\mathrm{target},B}$ | Component fraction of Balmer luminosity; chosen target weight. These are not area fractions. |
| $\mathcal B_j$, $\eta_B$ | Component Halpha/Hbeta decrement; compact-HII share of young Balmer emission. |
| N2, S2, O3 | [N II]6583/Halpha; [S II]6716+6731/Halpha; [O III]5007/Hbeta. |
| $\mathcal T_{\mathrm{leak}}$, $\mathcal Q_{\mathrm{abs,leak}}$ | Photon transport kernel, inverse area; leakage photon absorption rate per area. |
| $\mathscr R$, SED, $U$, $\Phi_{\mathrm H}$ | Photoionization prediction; spectral energy distribution; ionization parameter; incident photon flux in $\mathrm{s^{-1}\,cm^{-2}}$. |
| $Z$, N/O | Gas metallicity and nitrogen-to-oxygen abundance ratio; the calculation does not evolve them. |

**Table B4. Observational and mathematical conventions.**

| Term | Meaning |
|:--|:--|
| SF, NSF, ND | The source notebook's star-forming, other Balmer-detected, and joint Balmer non-detection categories. |
| RPS, DIG, OB, AGN | Ram-pressure stripping; diffuse ionized gas; massive O/B stellar population; active galactic nucleus. |
| BPT, EW, S/N, IMF | Standard optical excitation diagrams; equivalent width; signal-to-noise; initial mass function. |
| $\ln$, $\log_{10}$, dex | Natural logarithm, base-ten logarithm, and a base-ten logarithmic difference. |
| $O(t^3)$ | Terms cubic or higher in the short-time expansion, not a new physical parameter. |
| Superscripts obs, pred, required | Observed mean, conditional prediction, and an algebraic value required under the stated assumptions. |

# Appendix C. Provenance, checks, and reproducibility

## C.1 Source-to-model ledger

| Ingredient | Source and status |
|:--|:--|
| Local regulator notation and molecular supply | Huang et al. (2026), section 4; $\Sigma_\Phi$ is a rate, $\tau_\Phi$ a timescale. |
| Return/feedback mass balance | Regulator convention in Lilly et al. (2013) and Huang et al. (2026); effective phase assignment specified here. |
| Constant molecular depletion approximation | Leroy et al. (2013); our numerical 2 Gyr is illustrative, with helium included consistently. |
| HI deficiency and molecular/SFR response | Fumagalli et al. (2009); our exact linear transfer law is additional. |
| Excluding direct H2 stripping | A deliberate restricted hypothesis; Boselli et al. (2014) motivates the limitation. |
| Halpha stellar response | Kennicutt & Evans (2012); exponential 3-Myr kernel is our chosen approximation. |
| Case-B conversion | Hummer & Storey (1987); adopted numerical factor explicitly given by Cid Fernandes et al. (2011), equation 2. |
| HOLMES normalization | Belfiore et al. (2022), section 3.2 and footnote 5; current-mass convention retained. |
| Linear emission-ratio mixing | Blanc et al. (2009), equations 7--8; generalized here to three components and the correct Balmer denominator. |
| Spectral non-universality | Byler et al. (2019); Belfiore et al. (2022). Fixed ratios remain a testable assumption. |
| Analytical solutions, peaks, bounds, and inversion | Derived explicitly in this report from its balances and luminosity definitions. |

The Belfiore HTML sections, Huang section 4, Cid Fernandes equation 2 and mass convention, Blanc equations 7--8, Leroy mass convention, and Byler spectral dependencies were checked against primary sources. The Hummer--Storey abstract confirms the Case-B framework; the adopted numerical factor was checked through Cid Fernandes rather than independently interpolating the original recombination tables. Lilly metadata and regulator context were checked, but its PDF did not load reliably; no exact Lilly equation number is claimed.

## C.2 Numerical work actually executed

The reproducible script is `assets/20261001_HI_stripping_HOLMES/model_predictions.py`. It rebuilds `mauve_anchor_table.csv`, saves the gas/Halpha and three-component fading curves, performs the line-budget inversion, generates the three figures, and writes `numerical_audit.json`. The command is:

```bash
MPLCONFIGDIR=/private/tmp/mauve_20261001/mpl \
/opt/miniconda3/envs/ICRAR/bin/python \
  /Users/Igniz/Desktop/ICRAR/MAUVE/assets/20261001_HI_stripping_HOLMES/model_predictions.py
```

**Table C1. Executed mathematical checks.**

| Check | Fresh result |
|:--|:--|
| Analytical gas columns versus independent ODE, three branches | Maximum relative difference $7.69\times10^{-14}$ |
| Analytical young-Halpha response versus independent ODE | Maximum absolute normalized difference $3.33\times10^{-11}$ |
| Non-balanced equal-rate solution versus ODE | Maximum absolute difference $2.45\times10^{-15}$ |
| Peak condition, positive supply bound, total initial mass bound | Passed for the implemented examples |
| Unequal-decrement Balmer-weight identity | Zero difference at printed floating-point precision |
| Rising N2 with falling [N II] in three-component demonstration | Verified along both illustrative trajectories |
| Local source fingerprints | Six recorded notebook/pipeline/catalogue files match the extraction record |

The observational input is `assets/20260914_resolved_RPS_academic_model/stage_bpt_line_profiles_with_hbeta.csv`; its SHA-256 is recorded in the numerical audit. Fingerprints validate the stated extraction provenance, not the present contents of every large emission-map FITS file. No full science pipeline, bootstrap, stellar-population fit, radiative-transfer calculation, or photoionization grid was run. No new formal fit to the complete MAUVE sample is claimed.

The report follows the equation-centered style of the user's `20251109 The origin of correlation between metallicity and SFR.md`, with the intermediate product-rule and integration steps made explicit. Prior dated reports remain unchanged. PDF acceptance includes mathematical rendering, internal links, text boundaries, and visual inspection of rendered pages; the final checks are recorded separately in `assets/20261001_HI_stripping_HOLMES/verification_summary.md`.

# References

<span id="ref-belfiore"></span>
**Belfiore, F., et al. (2022).** *A tale of two DIGs: The relative role of H II regions and low-mass hot evolved stars in powering the diffuse ionised gas in PHANGS-MUSE galaxies.* A&A, 659, A26. [DOI](https://doi.org/10.1051/0004-6361/202141859); [primary full text](https://arxiv.org/html/2111.14876v3).

<span id="ref-blanc"></span>
**Blanc, G. A., Heiderman, A., Gebhardt, K., Evans, N. J. II, & Adams, J. (2009).** *The Spatially Resolved Star Formation Law from Integral Field Spectroscopy: VIRUS-P Observations of NGC 5194.* ApJ, 704, 842--862. [DOI](https://doi.org/10.1088/0004-637X/704/1/842); [primary manuscript](https://arxiv.org/pdf/0908.2810).

<span id="ref-boselli"></span>
**Boselli, A., et al. (2014).** *Cold gas properties of the Herschel Reference Survey. III. Molecular gas stripping in cluster galaxies.* A&A, 564, A67. [DOI](https://doi.org/10.1051/0004-6361/201322313); [primary manuscript](https://arxiv.org/abs/1402.0326).

<span id="ref-byler"></span>
**Byler, N., et al. (2019).** *Self-consistent predictions for LIER-like emission lines from post-AGB stars.* AJ, 158, 2. [DOI](https://doi.org/10.3847/1538-3881/ab1b70); [primary manuscript](https://arxiv.org/pdf/1904.10978).

<span id="ref-cid"></span>
**Cid Fernandes, R., Stasinska, G., Mateus, A., & Vale Asari, N. (2011).** *A comprehensive classification of galaxies in the Sloan Digital Sky Survey: how to tell true from fake AGN?* MNRAS, 413, 1687--1699. [DOI](https://doi.org/10.1111/j.1365-2966.2011.18244.x); [primary paper](https://minerva.ufsc.br/starlight/files/papers/j.1365-2966.2011.18244.x.pdf).

<span id="ref-fumagalli"></span>
**Fumagalli, M., Krumholz, M. R., Prochaska, J. X., Gavazzi, G., & Boselli, A. (2009).** *Molecular hydrogen deficiency in HI-poor galaxies and its implications for star formation.* ApJ, 697, 1811--1821. [DOI](https://doi.org/10.1088/0004-637X/697/2/1811); [primary manuscript](https://arxiv.org/abs/0903.3950).

<span id="ref-huang"></span>
**Huang, R., et al. (2026).** *MAUVE-MUSE: When Metallicity Follows or Fights Star Formation--A Mass-Dependent Inversion in Virgo Galaxies.* MNRAS, 549, stag1019. [DOI](https://doi.org/10.1093/mnras/stag1019); [primary manuscript](https://arxiv.org/html/2605.31412v1).

<span id="ref-hs"></span>
**Hummer, D. G., & Storey, P. J. (1987).** *Recombination-line intensities for hydrogenic ions--I. Case B calculations for HI and HeII.* MNRAS, 224, 801--820. [DOI](https://doi.org/10.1093/mnras/224.3.801).

<span id="ref-ke"></span>
**Kennicutt, R. C., Jr., & Evans, N. J. II (2012).** *Star Formation in the Milky Way and Nearby Galaxies.* ARA&A, 50, 531--608. [DOI](https://doi.org/10.1146/annurev-astro-081811-125610); [primary manuscript](https://arxiv.org/abs/1204.3552).

<span id="ref-leroy"></span>
**Leroy, A. K., et al. (2013).** *Molecular Gas and Star Formation in Nearby Disk Galaxies.* AJ, 146, 19. [DOI](https://doi.org/10.1088/0004-6256/146/2/19); [primary manuscript](https://arxiv.org/pdf/1301.2328).

<span id="ref-lilly"></span>
**Lilly, S. J., Carollo, C. M., Pipino, A., Renzini, A., & Peng, Y. (2013).** *Gas Regulation of Galaxies: The Evolution of the Cosmic Specific Star Formation Rate, the Metallicity-Mass-Star-formation Rate Relation, and the Stellar Content of Halos.* ApJ, 772, 119. [DOI](https://doi.org/10.1088/0004-637X/772/2/119); [primary manuscript](https://arxiv.org/abs/1303.5059).
