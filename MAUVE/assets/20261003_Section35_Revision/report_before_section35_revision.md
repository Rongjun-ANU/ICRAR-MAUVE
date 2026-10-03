---
title: "From HI stripping to fading young-star emission and HOLMES line ratios"
subtitle: "Model 0: explicit derivation, spatial enhancement, and numerical limits"
author: "Research report prepared for Rongjun Huang"
date: "2 October 2026"
lang: en
---

# 1. Physical picture and scope

We consider a local region in which ram-pressure stripping removes atomic gas but does not directly remove molecular gas. The loss of HI reduces the subsequent molecular supply. The retained H2 reservoir then supports star formation while it is gradually consumed. Young-star Halpha emission follows this declining SFR. If old stars and absorbing gas remain, hot low-mass evolved stars (HOLMES) can provide a more slowly varying ionizing contribution, whose fraction of the total Balmer emission increases. A forbidden-line-to-Balmer ratio rises if the HOLMES-powered emitting component has the larger intrinsic ratio.

We call the restricted analytical calculation **Model 0**. Its conversion time, molecular depletion time, recycling fraction, feedback loading, and HI stripping coefficient are constant in time at a specified location. The main spectral calculation combines compact HII emission and leaked-OB-powered emission into one effective young-star component. It adds a constant HOLMES component and adopts the same intrinsic Halpha/Hbeta ratio, 2.86, for both. These assumptions isolate one mechanism and permit an explicit connection from the gas columns to the line ratios; they do not establish a complete explanation of the MAUVE trends.

The main derivation distinguishes spatial enhancement from temporal growth and uses a direct SFR-to-Halpha approximation, a population-normalized HOLMES luminosity, and common Balmer weights. The finite stellar response, the HII/leakage partition, variable coefficients, unequal decrements, and spatial transport are developed separately in the appendices.

The numerical tests are deliberately retained even where the model fails. For the selected MAUVE means, normal-disc molecular consumption is too slow to produce the required fading within one Gyr, and the fiducial local HOLMES budget is too small to reproduce the N2 change with the illustrative spectra. The O3 endpoints require an unphysical negative HOLMES ratio in the fixed-spectrum inversion. These are conditional consistency results, not a statistical rejection of HOLMES or an inferred orbital clock.


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

Here $\tau_{\mathrm{dep}}$ is the molecular depletion time defined relative to the total rate of star formation. The gas supply-rate surface density $\Sigma_\Phi(t)$ feeds the molecular reservoir. The conversion time $\tau_{\mathrm{conv}}$ describes the assumed net transfer from the local atomic phase. The second equality is **our linear closure**, not a measured conversion law or an equation established by Huang et al. Their replenishment timescale is instead

$$
\tau_\Phi(t)\equiv
\frac{\Sigma_{\mathrm{H_2}}(t)}{\Sigma_\Phi(t)}
=\tau_{\mathrm{conv}}
\frac{\Sigma_{\mathrm{H_2}}(t)}{\Sigma_{\mathrm{HI}}(t)}.
\tag{2}
$$

Consequently, $\Sigma_\Phi(t)$ is a rate, whereas $\tau_\Phi$ and $\tau_{\mathrm{conv}}$ are times. Even when $\tau_{\mathrm{conv}}$ is constant, $\tau_\Phi$ generally evolves. We do not impose $\tau_\Phi\leq\tau_{\mathrm{dep}}$ on a supply-starved system; that inequality is not a consequence of the definition.

Let $R$ be the prompt stellar mass return fraction and $\lambda$ the feedback mass-loading factor, so that the feedback mass loss is $\lambda\Sigma_{\mathrm{SFR}}(t)$. The net removal associated with star formation and feedback is $(1-R+\lambda)\Sigma_{\mathrm{SFR}}(t)$. This is the usual regulator bookkeeping, here assigned effectively to the molecular reservoir ([Lilly et al. 2013](#ref-lilly); [Huang et al. 2026](#ref-huang)). Treating recycled material as promptly available to this reservoir is an approximation; explicit phase-dependent recycling would require additional terms.

We assume positive initial gas columns and conversion/depletion times, $0\leq R<1$, $\lambda\geq0$, and $\gamma_{\mathrm{strip}}\geq0$. The resulting reservoir response rates are positive, as required for the decline and peak statements below.

We use only $\gamma_{\mathrm{strip}}$ for the direct RPS loss coefficient. It has units of inverse time and acts only on HI. The baseline has no direct molecular stripping, no external supply to HI after $t=0$, and no explicit lateral transport. In Model 0, $\tau_{\mathrm{conv}}$, $\tau_{\mathrm{dep}}$, $R$, $\lambda$, and $\gamma_{\mathrm{strip}}$ are constant in time at fixed $\boldsymbol{x}$. Their possible spatial dependence is retained conceptually. Constant $\tau_{\mathrm{dep}}$ means constant molecular efficiency, not constant SFR. A linear molecular law is an empirical first approximation in nearby discs, with substantial environmental and scale-dependent limitations ([Leroy et al. 2013](#ref-leroy)). The HI-only stripping choice is a hypothesis for this calculation, not a claim that molecular stripping never occurs; observations provide counterexamples ([Boselli et al. 2014](#ref-boselli)).

## 2.2 Continuity equations and the meaning of the rates

Removing each mass flux from the reservoir that supplies it gives

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{HI}}(t)}{dt}
&=-\Sigma_\Phi(t)-\gamma_{\mathrm{strip}}\Sigma_{\mathrm{HI}}(t),\\
\frac{d\Sigma_{\mathrm{H_2}}(t)}{dt}
&=\Sigma_\Phi(t)-(1-R+\lambda)\Sigma_{\mathrm{SFR}}(t).
\end{aligned}
\tag{3}
$$

These are the model's local mass balances. Their sum explicitly cancels the internal transfer:

$$
\frac{d}{dt}(\Sigma_{\mathrm{HI}}(t)+\Sigma_{\mathrm{H_2}}(t))
=-\gamma_{\mathrm{strip}}\Sigma_{\mathrm{HI}}(t)
-(1-R+\lambda)\Sigma_{\mathrm{SFR}}(t).
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
\frac{\Sigma_{\mathrm{SFR}}(t)}{M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}}
=10^{-3}
\frac{\Sigma_{\mathrm{H_2}}(t)/(M_\odot\,\mathrm{pc}^{-2})}
{\tau_{\mathrm{dep}}/\mathrm{Gyr}}.
\tag{6}
$$

The factor is $10^6/10^9$: square parsecs to square kiloparsecs, then Gyr to years.

# 3. Explicit solution and the conditions for SFR decline or enhancement

## 3.1 Solve the atomic reservoir

For positive initial column $\Sigma_{\mathrm{HI},0}$, substituting equation (1) into equation (3), separating the variables, and integrating from the initial state gives

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{HI}}(t)}{dt}&=-\gamma_{\mathrm{HI}}\Sigma_{\mathrm{HI}}(t),\\
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

The decline of $\Sigma_\Phi(t)$ is the indirect route by which atomic stripping affects molecular gas. The physical possibility of molecular depletion following atomic deficiency has observational and theoretical antecedents ([Fumagalli et al. 2009](#ref-fumagalli)); the exact exponential form here follows from our stated closure.

## 3.2 Molecular reservoir and explicit SFR solution

First, insert equation (8) into the second balance:

$$
\frac{d\Sigma_{\mathrm{H_2}}(t)}{dt}
+\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2}}(t)
=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\tag{9}
$$

Multiplication by $e^{\gamma_{\mathrm{H_2}}t}$ makes the left-hand side a product derivative. To see this explicitly,

$$
\begin{aligned}
\frac{d}{dt}\left[e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2}}(t)\right]
&=e^{\gamma_{\mathrm{H_2}}t}\frac{d\Sigma_{\mathrm{H_2}}(t)}{dt}
+\gamma_{\mathrm{H_2}}e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2}}(t),\\
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

For $\gamma_{\mathrm{HI}}\ne\gamma_{\mathrm{H_2}}$, the integral is

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

If the HI-reservoir response rate and molecular-consumption rate are equal, $\gamma_{\mathrm{HI}}=\gamma_{\mathrm{H_2}}$, then $1/\tau_{\mathrm{conv}}+\gamma_{\mathrm{strip}}=(1-R+\lambda)/\tau_{\mathrm{dep}}$. The integrand in equation (11) is unity. There is no physical divergence:

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
\frac{1}{\Sigma_{\mathrm{SFR}}(t)}\frac{d\Sigma_{\mathrm{SFR}}(t)}{dt}
&=\frac{1}{\Sigma_{\mathrm{H_2}}(t)}\frac{d\Sigma_{\mathrm{H_2}}(t)}{dt}\\
&=\frac{\Sigma_\Phi(t)}{\Sigma_{\mathrm{H_2}}(t)}-\gamma_{\mathrm{H_2}}\\
&=\frac{1}{\tau_\Phi(t)}-\gamma_{\mathrm{H_2}}.
\end{aligned}
\tag{17}
$$

Thus the sign follows from supply relative to consumption:

$$
\begin{aligned}
\left.\frac{d\Sigma_{\mathrm{SFR}}(t)}{dt}\right|_0>0
&\ \Longleftrightarrow\ 
\Sigma_{\Phi,0}>\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0}
\ \Longleftrightarrow\ 
\tau_{\Phi,0}<\frac{\tau_{\mathrm{dep}}}{1-R+\lambda},\\
\left.\frac{d\Sigma_{\mathrm{SFR}}(t)}{dt}\right|_0<0
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

## 3.4 A temporal maximum for an initially over-supplied region

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

This is a positive-time peak of the temporal SFR history, not a measurement of spatial enhancement. It occurs when $\tau_{\Phi,0}^{-1}>\gamma_{\mathrm{H_2}}$. For equal rates, differentiating $(1+t/\tau_{\Phi,0})e^{-\gamma_{\mathrm{H_2}}t}$ gives the finite limit $t_{\mathrm{peak}}=1/\gamma_{\mathrm{H_2}}-\tau_{\Phi,0}$.

To verify that the stationary point is a maximum, differentiate the molecular balance once more:

$$
\begin{aligned}
\frac{d^2\Sigma_{\mathrm{H_2}}(t)}{dt^2}
&=-\gamma_{\mathrm{HI}}\Sigma_\Phi(t)
-\gamma_{\mathrm{H_2}}\frac{d\Sigma_{\mathrm{H_2}}(t)}{dt},\\
\left.\frac{d^2\Sigma_{\mathrm{H_2}}(t)}{dt^2}\right|_{t_{\mathrm{peak}}}
&=-\gamma_{\mathrm{HI}}\Sigma_\Phi(t_{\mathrm{peak}})<0.
\end{aligned}
\tag{25}
$$

Every stationary point has negative curvature. Therefore a solution starting with positive slope has one maximum followed by decline; a solution starting with nonpositive slope cannot develop a later minimum and rise within this constant-coefficient closed model.



## 3.5 A spatial SFR excess is not a positive temporal derivative

The observable motivated by the NGC4654 gradient discussion is a spatial excess relative to another region or a reference relation. We define that comparison explicitly:

$$
\Delta_{\mathrm{spatial}}\log_{10}\Sigma_{\mathrm{SFR}}(t)
\equiv\log_{10}\frac{\Sigma_{\mathrm{SFR}}^{\mathrm{leading}}(t)}
{\Sigma_{\mathrm{SFR}}^{\mathrm{reference}}(t)}.
\tag{26}
$$

A positive value does not determine the sign of $d\Sigma_{\mathrm{SFR}}^{\mathrm{leading}}/dt$. A region can already be declining and still lie above its reference. Moreover, a facing-to-opposite contrast can arise from suppression on the opposite side. The potential NGC4654 signal is motivation, not a fitted constraint or a freshly established detection in this report.

Our working interpretation is that compression and/or gas transport can first establish elevated local columns. Model 0 begins after that unresolved phase and predicts the subsequent evolution. For a short idealized accumulation episode, let $c_{\mathrm{HI}}$ and $c_{\mathrm{H_2}}$ be the ratios of the post-episode initial columns to reference initial columns. They are dimensionless initial-condition factors, not changes in conversion efficiency:

$$
\begin{aligned}
\Sigma_{\mathrm{HI},0}^{\mathrm{leading}}&=c_{\mathrm{HI}}\Sigma_{\mathrm{HI},0}^{\mathrm{reference}},\\
\Sigma_{\mathrm{H_2},0}^{\mathrm{leading}}&=c_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0}^{\mathrm{reference}},\\
\frac{\Sigma_{\mathrm{SFR},0}^{\mathrm{leading}}}{\Sigma_{\mathrm{SFR},0}^{\mathrm{reference}}}
&=c_{\mathrm{H_2}},\qquad
\tau_{\Phi,0}^{\mathrm{leading}}
=\frac{c_{\mathrm{H_2}}}{c_{\mathrm{HI}}}\tau_{\Phi,0}^{\mathrm{reference}}.
\end{aligned}
\tag{27}
$$

The last two relations assume identical $\tau_{\mathrm{dep}}$ and $\tau_{\mathrm{conv}}$ in the comparison. Thus $c_{\mathrm{H_2}}>1$ produces an elevated initial SFR. If both phases increase by the same factor, the replenishment time and the initial fractional SFR slope are unchanged. If the reference was initially balanced, continuing temporal growth requires $c_{\mathrm{HI}}>c_{\mathrm{H_2}}$; an elevated SFR level does not. Adding these columns to a fixed patch requires transport or a change in its physical area; the initial-condition prescription does not create mass within the closed evolution equations.

This separation is consistent with the early-stage interpretation in [Brown et al. (2023), section 3.3 and Figure 5](#ref-brown): enhanced outer-disc SFR is associated with greater molecular gas surface density at fixed stellar density, while molecular SFE is consistent with the field. Their early-RPS subset contains four galaxies, and their later-stage results also show lower SFE. These observations motivate a fixed-efficiency baseline; they do not establish constant efficiency for every MAUVE region. Nor do they determine the compression history or prove that $\tau_{\mathrm{conv}}$ decreases. Turbulence and compression can affect several processes, so we impose neither sign of its environmental response here.

The omitted transport contribution has the sign $-\boldsymbol{\nabla}\cdot(\Sigma_i\boldsymbol v_i)$ on the right-hand side of the continuity equation. A negative mass-flux divergence contributes positively to the local gas column. Appendix E writes the complete balances. The main model does not solve for a velocity field or the accumulation episode.

# 4. From instantaneous SFR to young-star Halpha emission

Let $\mathcal L_\alpha^{\mathrm{young}}$ be the Halpha luminosity per adopted area powered by the young stellar population. It includes both compact HII emission and emission powered by OB photons absorbed outside compact HII regions. For the main calculation these are one effective component. Their separate luminosities and spectra are unnecessary until Appendix C.

The conversion from SFR to an ionizing population normally averages the recent star formation history over the lifetimes of massive stars ([Kennicutt & Evans 2012](#ref-ke)). For gas evolution much slower than that stellar response, we use

$$
\boxed{\mathcal L_\alpha^{\mathrm{young}}(t)
\simeq\frac{f_{\mathrm{young}}}{C_\alpha}\Sigma_{\mathrm{SFR}}(t)
=\mathcal L_{\alpha,0}^{\mathrm{young}}F_{\mathrm{SFR}}(t).}
\tag{28}
$$

Here $f_{\mathrm{young}}$ is the constant fraction of the young ionizing photon budget absorbed by hydrogen in the modeled region, with $0<f_{\mathrm{young}}\leq1$ in this local approximation. $C_\alpha$ is the fully absorbed Halpha-to-SFR calibration, not its inverse. We use the existing MAUVE value $C_\alpha=4.9835821\times10^{-42}\ M_\odot\,\mathrm{yr^{-1}}/(\mathrm{erg\,s^{-1}})$; division converts $M_\odot\,\mathrm{yr^{-1}\,kpc^{-2}}$ to $\mathrm{erg\,s^{-1}\,kpc^{-2}}$. The numerical benchmark sets $f_{\mathrm{young}}=1$. The calibration depends on the stellar population and IMF and is not universal.

Appendix A derives the normalized response kernel and its exact exponential solution. For the illustrative 3-Myr response and the present gas coefficients, the normalized one-Gyr response is 0.79937, compared with instantaneous $F_{\mathrm{SFR}}=0.79867$: a relative correction of 0.0881%. This supports equation (28) for the smooth Model 0 history. It does not justify an instantaneous Halpha jump at a sudden SFR discontinuity. If the initial gas accumulation was recent on a few-Myr timescale, the actual prehistory must enter the response calculation.

There is also a spatial condition: young emission in an NSF patch may be powered by photons from neighbouring star-forming regions. In that case a local Halpha luminosity need not trace the local SFR. Equation (28) is then a coarse-grained or local-absorption approximation. This limitation is especially relevant to the illustrative gas normalization in section 8; Appendix E gives the nonlocal form.

# 5. A compact, physically normalized HOLMES contribution

Let $\Sigma_*^{\mathrm{old}}$ be the current mass surface density in the old population, including its associated remnants, and let $q_{\mathrm{H,HOLMES}}$ be the production rate of hydrogen-ionizing photons per unit of that current mass. Their product is the emitted photon production per area. Multiplying by $f_{\mathrm{abs,HOLMES}}$ gives the photon rate actually absorbed by hydrogen in the gas assigned to the region:

$$
\begin{aligned}
\mathcal Q_{\mathrm{HOLMES}}&=q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}},\\
\mathcal Q_{\mathrm{abs,HOLMES}}&=f_{\mathrm{abs,HOLMES}}\mathcal Q_{\mathrm{HOLMES}}.
\end{aligned}
\tag{29}
$$

The units are $(\mathrm{s^{-1}}M_\odot^{-1})(M_\odot\,\mathrm{kpc^{-2}})=\mathrm{s^{-1}\,kpc^{-2}}$. In ionization equilibrium, one absorbed ionizing photon balances a Case-B recombination. The probability that such a recombination produces Halpha is $p_\alpha=\alpha_\alpha^{\mathrm{eff}}/\alpha_B$. Each emitted Halpha photon carries energy $h_{\mathrm P}\nu_\alpha$. Therefore

$$
\begin{aligned}
\mathcal L_\alpha^{\mathrm{HOLMES}}
&=h_{\mathrm P}\nu_\alpha p_\alpha\mathcal Q_{\mathrm{abs,HOLMES}}\\
&=\boxed{\epsilon_\alpha f_{\mathrm{abs,HOLMES}}
q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}}},\\
\epsilon_\alpha&\equiv h_{\mathrm P}\nu_\alpha p_\alpha,\qquad
p_\alpha\equiv\frac{\alpha_\alpha^{\mathrm{eff}}}{\alpha_B}\simeq\frac{1}{2.206}.
\end{aligned}
\tag{30}
$$

$\alpha_B$ excludes recombinations directly to the hydrogen ground state; $\alpha_\alpha^{\mathrm{eff}}$ counts recombinations yielding Halpha. Both have units $\mathrm{cm^3\,s^{-1}}$, so $p_\alpha$ is dimensionless. With $\lambda_\alpha=6562.8$ Angstrom, $\epsilon_\alpha=1.3721\times10^{-12}$ erg per absorbed photon. The recombination framework comes from [Hummer & Storey (1987)](#ref-hs); the adopted numerical conversion is explicitly given in [Cid Fernandes et al. (2011), equation 2](#ref-cid). Appendix F derives the same relation by eliminating the volume emission measure. That paper's population normalization uses formed mass, whereas our $q$ and $\Sigma_*^{\mathrm{old}}$ consistently use current mass.

For the fiducial population we adopt $q_{\mathrm{H,HOLMES}}=7\times10^{40}\ \mathrm{s^{-1}}M_\odot^{-1}$ from the PEGASE normalization used by [Belfiore et al. (2022), section 3.2 and footnote 5](#ref-belfiore). Their 10-Gyr solar-metallicity FSPS comparison gives $5\times10^{40}$ in the current-stars-plus-remnants convention. The proportionality is valid for a specified population, not a universal photon yield per unit stellar mass.

The proposed slowly varying contribution requires four separate assumptions:

| Factor | Meaning of the Model 0 approximation |
|:--|:--|
| $\epsilon_\alpha$ | Recombination conditions remain near the adopted low-density Case-B conditions. |
| $q_{\mathrm{H,HOLMES}}$ | The old population's specific ionizing output evolves slowly over the modeled interval. A 1--2 Gyr interval is not negligible for every age mixture. |
| $\Sigma_*^{\mathrm{old}}$ | The old stellar mass assigned to the region changes little; newly formed stars are not automatically added to this old component. |
| $f_{\mathrm{abs,HOLMES}}$ | Sufficient gas, covering fraction, and recombination capacity persist to absorb approximately the same fraction of the assigned old-star photons. |

Under these assumptions $\mathcal L_\alpha^{\mathrm{HOLMES}}(t)\simeq\mathcal L_{\alpha,0}^{\mathrm{HOLMES}}$. This is a controlled retained-gas approximation. Balmer detection demonstrates emitting gas; it proves neither HOLMES domination nor a nonzero or constant absorption fraction for HOLMES specifically. A nearly constant source cannot maintain a fixed Halpha floor after the absorbing gas has been removed.

No universal measured value of $f_{\mathrm{abs,HOLMES}}$ is available for the MAUVE NSF regions here. It stays symbolic in the derivation. Setting it to one, together with $\Sigma_*^{\mathrm{old}}=\Sigma_*$, is only a maximal local benchmark for the specified population in section 8. Lower values weaken the contribution. Appendix F shows that strict constancy is unnecessary: it suffices that the old contribution fades more slowly in fractional terms than the young contribution.

# 6. How the luminosity weights change

Add the two positive Halpha contributions before defining their weights:

$$
\begin{aligned}
\mathcal L_\alpha(t)&=\mathcal L_\alpha^{\mathrm{young}}(t)+\mathcal L_\alpha^{\mathrm{HOLMES}},\\
w_{\mathrm{HOLMES}}(t)&\equiv
\frac{\mathcal L_\alpha^{\mathrm{HOLMES}}}
{\mathcal L_\alpha^{\mathrm{young}}(t)+\mathcal L_\alpha^{\mathrm{HOLMES}}},\qquad
w_{\mathrm{young}}(t)=1-w_{\mathrm{HOLMES}}(t).
\end{aligned}
\tag{31}
$$

These are light fractions, not fractions of area or gas mass. Write $w_{\mathrm{HOLMES},0}=w_{\mathrm{HOLMES}}(0)$. Dividing numerator and denominator by the initial total Halpha luminosity and using equation (28) gives

$$
\boxed{w_{\mathrm{HOLMES}}(t)=
\frac{w_{\mathrm{HOLMES},0}}
{w_{\mathrm{HOLMES},0}+(1-w_{\mathrm{HOLMES},0})F_{\mathrm{SFR}}(t)}.}
\tag{32}
$$

For fixed positive HOLMES luminosity, the quotient rule gives

$$
\begin{aligned}
\frac{dw_{\mathrm{HOLMES}}(t)}{dt}
&=-\frac{\mathcal L_\alpha^{\mathrm{HOLMES}}}
{[\mathcal L_\alpha(t)]^2}\frac{d\mathcal L_\alpha^{\mathrm{young}}(t)}{dt}\\
&=-w_{\mathrm{HOLMES}}(t)[1-w_{\mathrm{HOLMES}}(t)]
\frac{d\ln\Sigma_{\mathrm{SFR}}(t)}{dt}\\
&=w_{\mathrm{HOLMES}}(t)[1-w_{\mathrm{HOLMES}}(t)]
\left[\gamma_{\mathrm{H_2}}-\frac{1}{\tau_\Phi(t)}\right].
\end{aligned}
\tag{33}
$$

The second line uses constant $f_{\mathrm{young}}$ and $C_\alpha$; the third uses equation (17). This connects the gas balance directly to the changing source weight. Molecular consumption exceeding replenishment makes SFR decline and the HOLMES fraction rise. An initially over-supplied region has the opposite response until its SFR maximum. Proportional fading of two young components alone would leave their mutual weight unchanged; the independently supplied HOLMES term is what changes this conclusion.

# 7. The forbidden-line ratios and the complete connection

For a forbidden line $\ell$, let $B$ denote its Balmer denominator and define the linear component ratio $R_{\ell/B}^j=\mathcal L_\ell^j/\mathcal L_B^j$, with $j$ equal to young or HOLMES. N2 is [N II]6583/Halpha, S2 is ([S II]6716+[S II]6731)/Halpha, and O3 is [O III]5007/Hbeta. Shock and AGN emission are neglected as a baseline hypothesis, not established absent by these equations.

In Model 0 we adopt $\mathcal L_\alpha^j/\mathcal L_\beta^j=2.86$ for each component. Thus

$$
\frac{\mathcal L_\beta^{\mathrm{HOLMES}}}{\mathcal L_\beta}
=\frac{\mathcal L_\alpha^{\mathrm{HOLMES}}/2.86}{\mathcal L_\alpha/2.86}
=w_{\mathrm{HOLMES}}.
\tag{34}
$$

The common intrinsic decrement is an approximation consistent with the Case-B convention used for the observational dust correction ([Hummer & Storey 1987](#ref-hs)). It is not an assertion that every exported corrected Halpha/Hbeta measurement is exactly 2.86; Appendix D preserves the actual Hbeta values and quantifies the small numerical difference.

Substitute $\mathcal L_\ell^j=R_{\ell/B}^j\mathcal L_B^j$ into the total ratio and separate the terms:

$$
\begin{aligned}
R_{\ell/B}(t)
&=\frac{R_{\ell/B}^{\mathrm{young}}\mathcal L_B^{\mathrm{young}}(t)
+R_{\ell/B}^{\mathrm{HOLMES}}\mathcal L_B^{\mathrm{HOLMES}}}
{\mathcal L_B^{\mathrm{young}}(t)+\mathcal L_B^{\mathrm{HOLMES}}}\\
&=[1-w_{\mathrm{HOLMES}}(t)]R_{\ell/B}^{\mathrm{young}}
+w_{\mathrm{HOLMES}}(t)R_{\ell/B}^{\mathrm{HOLMES}}\\
&=\boxed{R_{\ell/B}^{\mathrm{young}}+
(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})w_{\mathrm{HOLMES}}(t)}.
\end{aligned}
\tag{35}
$$

This luminosity-weighted identity has an HII/DIG antecedent in [Blanc et al. (2009), equations 7--8](#ref-blanc). Their numerical [S II] template is for a single line and is not imported for our doublet sum. Ratios mix linearly; logarithmic BPT coordinates are taken only after addition.

For constant component spectra, differentiate equation (35):

$$
\frac{dR_{\ell/B}(t)}{dt}
=(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\frac{dw_{\mathrm{HOLMES}}(t)}{dt}.
\tag{36}
$$

When young emission fades, the sign of the ratio change is the sign of the spectral contrast. Harder ionization alone does not require every ratio to increase; temperature, ionic fractions, metallicity, N/O, and ionization parameter also matter ([Byler et al. 2019](#ref-byler)). N2, S2, and O3 must be tested separately. Both Balmer and forbidden-line luminosities can decline while their ratio rises because the Balmer line fades faster; Appendix F gives the explicit luminosity algebra.

Finally, inserting the gas solution, the young-star conversion, and the old-star normalization gives the complete Model 0 prediction:

$$
\boxed{\begin{aligned}
\mathcal L_\alpha(t)
&=\frac{f_{\mathrm{young}}\Sigma_{\mathrm{SFR},0}}{C_\alpha}F_{\mathrm{SFR}}(t)
+\epsilon_\alpha f_{\mathrm{abs,HOLMES}}q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}},\\
R_{\ell/B}(t)
&=R_{\ell/B}^{\mathrm{young}}
+(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\frac{\epsilon_\alpha f_{\mathrm{abs,HOLMES}}q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}}}
{\mathcal L_\alpha(t)}.
\end{aligned}}
\tag{37}
$$

Here $F_{\mathrm{SFR}}(t)$ is explicitly given by equation (16), or by its equal-rate limit from equation (15). The second line applies to Hbeta ratios as well because of equation (34). All luminosity factors and rate coefficients are defined; no arbitrary residual luminosity function is introduced. In the restricted model, Halpha is the observable connection between declining molecular supply and evolving line-ratio weights.

Only after deriving that relation should we identify its observational scope. SF and NSF are selections, not ionizing sources. SF in the source analysis requires a finite HII-selected SFR product, Halpha EW greater than 6 Angstrom, and intrinsic Halpha dispersion below $45\ \mathrm{km\,s^{-1}}$. ND is the joint Balmer non-detection category; NSF is the remaining Balmer-detected category outside SF. Neither component is switched off merely because a region is called SF or NSF. NSF is not synonymous with DIG, LIER, or HOLMES domination, and NSF occupancy is not $w_{\mathrm{HOLMES}}$.

The present model concerns emission within retained gas. Predicting changes in SF/NSF/ND area fractions requires the stellar continuum, EW, noise, widths, masks, and line-detection rules to be modeled and reapplied. Completely stripped outer regions also require an evolving gas absorption capacity. Those calculations are outside Model 0.


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

Take the bin center as a representative $\Sigma_*=10^{8.625}=4.217\times10^8\ M_\odot\,\mathrm{kpc}^{-2}$. For a deliberately generous **local benchmark**, assume all this mass is old, its mass convention matches current stars plus remnants, and all available HOLMES photons ionize the retained hydrogen. Equation (30) gives

$$
\begin{aligned}
\mathcal L_{\alpha}^{\mathrm{HOLMES}}
&=(1.3721\times10^{-12})(7\times10^{40})(4.217\times10^8)\\
&=4.050\times10^{37}\ \mathrm{erg\,s^{-1}\,kpc^{-2}}.
\end{aligned}
\tag{38}
$$

The corresponding pre-peak and post-peak NSF Halpha weights are 0.01138 and 0.03647. With the FSPS comparison normalization they would be smaller by $5/7$. A younger mass fraction or incomplete absorption also reduces them. The bin-center approximation, population uncertainty, mass-convention compatibility, and nonlocal photon transport prevent treating this as an absolute universal ceiling. It is the maximum **within the specified local fiducial population model**.

This numerical result is essential: a spatially broad HOLMES component may be real yet contribute little Halpha in a bright NSF sample. The relevant variable is $\mathcal L_\alpha/\Sigma_*^{\mathrm{old}}$, not Halpha brightness or stellar density separately.



## 8.3 The gas response and an elevated initial spatial state

For comparability with the previous report, subtract the fiducial HOLMES luminosity from the pre-peak NSF Halpha scale and apply $C_\alpha$ with $f_{\mathrm{young}}=1$. The resulting effective young SFR is $0.01753\ M_\odot\,\mathrm{yr^{-1}\,kpc^{-2}}$. With $\tau_{\mathrm{dep}}=2$ Gyr, equation (6) assigns $\Sigma_{\mathrm{H_2},0}=35.055\ M_\odot\,\mathrm{pc^{-2}}$.

This is a **luminosity-scaled numerical example**, not a measurement of local SFR or molecular mass in NSF. Nonlocal leaked photons would invalidate that local interpretation. A future gas fit should use independently measured CO/HI columns with compatible mass and area conventions. No such gas fit is performed here.

**Table 2. Constant Model 0 parameters.**

| Quantity | Value | Status |
|:--|--:|:--|
| $\tau_{\mathrm{dep}}$ | 2 Gyr | Illustrative normal-disc molecular-consumption scale |
| $R$, $\lambda$ | 0.4, 0 | Assumed prompt recycling and no feedback loss |
| $\Sigma_{\mathrm{H_2},0}$ | $35.055\ M_\odot\,\mathrm{pc^{-2}}$ | Assigned from the luminosity scale, not measured CO |
| $\Sigma_{\mathrm{HI},0}$ | $10\ M_\odot\,\mathrm{pc^{-2}}$ | Assumed; both phases include helium |
| $\tau_{\mathrm{conv}}$ | 0.95088 Gyr | Chosen so initial supply equals net molecular consumption |
| $\gamma_{\mathrm{strip}}$ | $3\ \mathrm{Gyr^{-1}}$ | Assumed HI loss coefficient |
| $\gamma_{\mathrm{HI}}$, $\gamma_{\mathrm{H_2}}$ | $4.05166$, $0.30000\ \mathrm{Gyr^{-1}}$ | Derived response rates |

Equation (19) yields $F_{\mathrm{SFR}}(1\ \mathrm{Gyr})=0.79867$, or a 0.0976-dex decline. The initial slope is zero and subsequent slopes are negative. The HI reservoir declines rapidly, but H2 buffers the SFR. The closed comparison with $\gamma_{\mathrm{strip}}=0$, all other quantities unchanged, gives 0.89706. Relative to that declining comparison the extra suppression is 0.0505 dex. The comparison is closed and has no fresh external supply; it is not a permanently maintained field equilibrium.

Now raise both initial gas columns by 1.5 without changing $\tau_{\mathrm{conv}}$, $\tau_{\mathrm{dep}}$, or the stripping coefficient. Linearity makes the entire stripped SFR curve 1.5 times the baseline stripped curve. Its initial SFR is 0.1761 dex above the uncompressed initial reference, but it has the same zero initial fractional slope and subsequent decline. At one Gyr its SFR is $1.5\times0.79867=1.1980$ in units of the original reference SFR, and it remains 0.1256 dex above the contemporaneous closed no-RPS reference, 0.89706.

This example makes the distinction explicit: an elevated spatial SFR and a negative temporal derivative coexist. The factor 1.5 is a chosen initial-condition illustration, not a fitted compression amplitude for NGC4654. The model evolves the accumulated gas but does not calculate the accumulation process.

![Figure 1. Model 0 at fixed conversion and depletion times. Left: SFR in units of the uncompressed initial reference; the orange curve begins with both gas columns multiplied by 1.5 and subsequently declines. Right: normalized atomic columns for the stripped and closed no-stripping cases. The orange atomic fraction coincides with the blue stripped fraction and is omitted. The elevated orange SFR is a spatial-level example, not a burst generated by switching on the HI sink.](assets/20261002_Model0_Derivation/figure01_SFR_response.png)


## 8.4 How much fading is required for the HOLMES fraction to matter?

Solving equation (32) for the remaining young fraction at a target HOLMES weight $w_{\mathrm{target}}$ gives

$$
F_{\mathrm{SFR}}
=\frac{w_{\mathrm{HOLMES},0}[1-w_{\mathrm{target}}]}
{w_{\mathrm{target}}[1-w_{\mathrm{HOLMES},0}]}.
\tag{39}
$$

For the bright pre-peak NSF anchor, reaching a 10% HOLMES Halpha weight requires $F_{\mathrm{SFR}}\simeq0.104$, and reaching 50% requires $F_{\mathrm{SFR}}\simeq0.0115$. A large change in the mixture therefore requires substantial fading when the initial old-star contribution is only about 1%.

For comparison, an explicitly hypothetical faint patch with the same old stellar density but initial $\mathcal L_\alpha=2\times10^{38}$ has a HOLMES weight of 0.203. Reducing its young emission by a factor of ten raises that weight to 0.717. With an effective young N2 ratio of 0.35 and a HOLMES N2 ratio of 1.5, its total N2 rises from approximately 0.583 to 1.175 while [N II] itself becomes fainter. The spectra here are illustrative endpoints, not measurements of individual MAUVE components.

Figure 2 uses only the effective young N2 ratio 0.35 and HOLMES N2 ratio 1.5. These positive constant spectra illustrate the mechanism without splitting HII and leakage. The remaining young fraction is varied directly; a hundredfold fading is not claimed to occur within the two-Gyr gas interval or within the validity of every fixed-population approximation.

![Figure 2. Effective young-star plus HOLMES fading at the bright NSF scale and a hypothetical faint scale. Left: increasing HOLMES Halpha weight. Middle: increasing N2 as total Halpha decreases. Right: [N II] still declines. The remaining young fraction decreases from left to right. These are illustrative component spectra, not a fit.](assets/20261002_Model0_Derivation/figure02_fading_and_ratio.png)

## 8.5 Test the amplitude against the NSF anchors

For each diagnostic, choose an illustrative fixed HOLMES ratio, normalize the effective young ratio to the pre-peak NSF mean, and then predict the post-peak ratio using its **observed Halpha luminosity** to calculate the common Model 0 weight and the fixed photon-budget term. The initial normalization follows by rearranging equation (35):

$$
R_{\ell/B}^{\mathrm{young}}
=\frac{R_{\ell/B,0}^{\mathrm{obs}}
-w_{\mathrm{HOLMES},0}R_{\ell/B}^{\mathrm{HOLMES}}}
{1-w_{\mathrm{HOLMES},0}}.
\tag{40}
$$

This is one-point calibration, not an independent prediction of the initial spectrum. The endmember ratios 1.5, 1.0, and 3.0 below are sensitivity choices, not values extracted from Belfiore, a CLOUDY grid, or the MAUVE spectra. We subsequently invert the equations so the conclusion does not rest only on those choices.

**Table 3. Conditional post-peak predictions.** Residuals are $\log_{10}(R^{\mathrm{pred}}/R^{\mathrm{obs}})$; no uncertainty or formal fit significance is assigned.

| Ratio | Assumed HOLMES | Calibrated young | Predicted post-peak | Observed post-peak | Residual (dex) |
|:--|--:|--:|--:|--:|--:|
| N2 | 1.5 | 0.3601 | 0.4017 | 0.5663 | -0.1492 |
| S2 | 1.0 | 0.3701 | 0.3931 | 0.4168 | -0.0255 |
| O3 | 3.0 | 0.8736 | 0.9511 | 0.7528 | +0.1016 |

The N2 prediction has the desired sign but insufficient amplitude. The O3 prediction has the opposite sign to the pre-to-post change of these exported means. This is not a proof against HOLMES ionization. It is a failure of this particular common-population, fixed-spectrum, additive interpretation of the selected mean points.

![Figure 3. Conditional fixed-HOLMES predictions compared with the exported NSF mean ratios. The initial young ratio is calibrated to the pre-peak point. Post-peak Halpha luminosities are supplied from the data. The differences therefore test the constant-spectrum/photon-budget assumptions rather than the gas clock. Model 0 uses one common Halpha weight; the actual-Hbeta variant is in Appendix D. Lines connect different galaxy samples for comparison and do not imply observed temporal tracks.](assets/20261002_Model0_Derivation/figure03_MAUVE_budget_test.png)

## 8.6 Invert the required spectrum or photon budget

Using the common weight, write equation (35) as a function of total Halpha luminosity:

$$
R_{\ell/B}=R_{\ell/B}^{\mathrm{young}}+\frac{(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})\mathcal L_\alpha^{\mathrm{HOLMES}}}{\mathcal L_\alpha}.
\tag{41}
$$

Evaluate it at two observed points indexed by 0 and 1. Subtracting the ratios removes the effective young ratio:

$$
\begin{aligned}
R_{\ell/B,1}-R_{\ell/B,0}
&=(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\mathcal L_\alpha^{\mathrm{HOLMES}}
\left(\frac{1}{\mathcal L_{\alpha,1}}-\frac{1}{\mathcal L_{\alpha,0}}\right),\\
R_{\ell/B}^{\mathrm{young}}
&=\frac{R_{\ell/B,0}\mathcal L_{\alpha,0}-R_{\ell/B,1}\mathcal L_{\alpha,1}}
{\mathcal L_{\alpha,0}-\mathcal L_{\alpha,1}}.
\end{aligned}
\tag{42}
$$

The second line follows by subtracting the ratio-times-Halpha quantities (the actual forbidden luminosity for N2 and S2, and 2.86 times that luminosity for the idealized O3 model), whose constant intercept cancels. Substituting it into either point gives the required HOLMES ratio at a specified old-star Balmer luminosity:

$$
R_{\ell/B}^{\mathrm{HOLMES,required}}
=R_{\ell/B}^{\mathrm{young}}
+\frac{R_{\ell/B,1}-R_{\ell/B,0}}
{\mathcal L_\alpha^{\mathrm{HOLMES}}(\mathcal L_{\alpha,1}^{-1}-\mathcal L_{\alpha,0}^{-1})}.
\tag{43}
$$

At the fiducial HOLMES budget, the two NSF endpoints require N2$_{\mathrm{HOLMES}}=7.99$, S2$_{\mathrm{HOLMES}}=1.94$, and O3$_{\mathrm{HOLMES}}=-4.82$. The negative O3 value is unphysical for a positive emitting component. The extreme required N2 value is not an adopted spectrum; it quantifies the burden placed on the simple mixture and should be checked against a self-consistent photoionization grid rather than accepted as a free fitting coefficient.

Alternatively, retaining an assumed HOLMES ratio gives

$$
\mathcal L_\alpha^{\mathrm{HOLMES,required}}
=\frac{R_{\ell/B,1}-R_{\ell/B,0}}
{(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
(\mathcal L_{\alpha,1}^{-1}-\mathcal L_{\alpha,0}^{-1})}.
\tag{44}
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
\tag{45}
$$

Under the instantaneous Model 0 approximation, define $F_{\mathrm{required}}\equiv F_{\alpha,\mathrm{required}}$ as the target value of $F_{\mathrm{SFR}}(t)$. Because the supplied molecular contribution is nonnegative, equation (21) bounds the instantaneous decline. Taking logarithms at the time when that target is reached gives the minimum time step by step:

$$
F_{\mathrm{SFR}}(t)\geq e^{-\gamma_{\mathrm{H_2}}t},\qquad
\ln F_{\mathrm{required}}\geq-\gamma_{\mathrm{H_2}}t,
\qquad t\geq\frac{-\ln F_{\mathrm{required}}}{\gamma_{\mathrm{H_2}}}.
\tag{46}
$$

For $F_{\mathrm{required}}=0.30426$ and $\gamma_{\mathrm{H_2}}=0.3\ \mathrm{Gyr^{-1}}$, the no-supply minimum is 3.966 Gyr. Solving the full balanced **instantaneous** Model 0 gives 4.223 Gyr; retaining the Appendix A stellar response gives 4.226 Gyr. The approximation in section 4 therefore leaves the physical tension unchanged. These long extrapolations are diagnostics of the assumed consumption rate, not inferred infall times, and extend beyond the intended 1--2 Gyr fixed-population interval.

For a hypothetical one-Gyr constraint, the same inequality requires

$$
\gamma_{\mathrm{H_2}}\geq-\frac{\ln F_{\mathrm{required}}}{1\ \mathrm{Gyr}}
\quad\Longrightarrow\quad
\tau_{\mathrm{dep}}\leq\frac{(1-R+\lambda)(1\ \mathrm{Gyr})}{-\ln F_{\mathrm{required}}}.
\tag{47}
$$

With $R=0.4$ and $\lambda=0$, this is $\tau_{\mathrm{dep}}\lesssim0.504$ Gyr; positive replenishment makes the requirement stricter. This is a conditional bound, not evidence that the observed stage separation is one Gyr or that the actual depletion time has this value.

## 8.8 Check absolute line fading, not only line ratios

Multiplying the exported ratios by their actual denominators gives the following pre-to-post luminosity fractions. No component decomposition is needed for this arithmetic.

| Line | Post-peak / pre-peak luminosity |
|:--|--:|
| Halpha | 0.3122 |
| Hbeta | 0.3128 |
| [N II]6583 | 0.4739 |
| [S II] doublet | 0.3449 |
| [O III]5007 | 0.2623 |

[N II] fades more slowly than Halpha, [S II] only slightly more slowly, and [O III] faster than Hbeta between these means. The observation is not that forbidden emission remains constant. A successful extension must reproduce the absolute luminosities and all three ratios together, including the nonmonotonic close-to-peak O3 point. These factors remain descriptive cross-sectional comparisons without propagated galaxy-level uncertainty.


# 9. Interpretation and limitations

Model 0 establishes a self-consistent conditional sequence. Atomic stripping reduces the future molecular supply; existing H2 buffers the SFR decline; young-star Halpha subsequently fades; and a retained old-star photon budget can become a larger fraction of the remaining Balmer emission. Fixed positive spectra can then produce a rising forbidden-to-Balmer ratio while both lines fade. No changing HII-to-leakage partition is needed for that mechanism.

The model also separates two statements that should not be conflated. A spatially enhanced molecular column produces an elevated SFR at fixed efficiency. Its subsequent derivative can already be negative. The closed reservoir solution predicts that subsequent evolution but does not generate the compression or transport that supplied the initial column. Switching on the HI sink alone does not generate an enhancement from molecular balance.

Its numerical limitations are substantive. The chosen 2-Gyr depletion time cannot yield the required fading within one Gyr even after replenishment stops. The fiducial local HOLMES budget gives only about 1.14% and 3.65% of the pre/post NSF Halpha means. With the assumed spectra the N2 increase is too small and O3 changes in the wrong direction. Allowing freely chosen but fixed spectra still demands a negative O3 HOLMES contribution at this photon budget. Lowering the absorbed fraction does not repair these particular amplitude tests.

These results do not prove that HOLMES are absent or that all NSF regions obey the same spectrum. The points combine different galaxies, the bin-center stellar mass is only representative, and the photon yield and stellar mass conventions carry uncertainties. The model uses a local photon budget and independent emitting templates. If OB and HOLMES photons illuminate the same gas, the ionic fractions and temperature respond jointly, so the forbidden spectrum need not equal the sum of two separately fixed source spectra ([Belfiore et al. 2022, sections 5.4--5.5](#ref-belfiore)). An increasingly important HOLMES spectrum and a declining ionization parameter can act together; their quantitative effects require a physical photoionization calculation, not independent adjustment of three ratios.

For application to MAUVE, the next constraints are concrete: use CO/HI rather than NSF Halpha alone to normalize the gas reservoirs; estimate the old-population photon budget from the fitted stellar populations; test ratios and absolute luminosities against $\mathcal L_\alpha/\Sigma_*^{\mathrm{old}}$ within matched stellar-density/radius support; and assess leakage from nearby HII regions. Galaxy-level uncertainty and spatial support must accompany any comparison among stages. A successful spectral extension should predict N2, S2, O3 and Balmer emission with one consistent gas state. The present report supplies the analytical reference against which such extensions can be tested; it does not claim a fit to MAUVE occupancy fractions or the full data set.

# Appendix A. The finite Halpha response and its short-timescale limit


Halpha traces the ionizing photons supplied by short-lived massive stars, not an infinitely instantaneous SFR. A stellar population produces an age-dependent ionizing output; summing stellar cohorts is therefore a convolution. This is the population-synthesis basis of recombination-line SFR indicators ([Kennicutt & Evans 2012](#ref-ke)). We approximate the response with a normalized kernel $K_\alpha(a)$, where $a\geq0$ is stellar age:

$$
\overline\Sigma_{\mathrm{SFR},\alpha}(t)
=\int_0^\infty K_\alpha(a)\Sigma_{\mathrm{SFR}}(t-a)\,da,
\qquad \int_0^\infty K_\alpha(a)\,da=1.
\tag{48}
$$

The overbar denotes this Halpha response average. For transparent analytic calculations we assume

$$
K_\alpha(a)=\frac{1}{\tau_{\mathrm{ion}}}e^{-a/\tau_{\mathrm{ion}}},
\qquad
\tau_{\mathrm{ion}}\frac{d\overline\Sigma_{\mathrm{SFR},\alpha}}{dt}
+\overline\Sigma_{\mathrm{SFR},\alpha}=\Sigma_{\mathrm{SFR}}.
\tag{49}
$$

The exponential kernel and numerical choice $\tau_{\mathrm{ion}}=3$ Myr are our approximations, not a fitted stellar-population model. For the continuous-onset gas solutions we assume a constant pre-onset SFR, so $\overline\Sigma_{\mathrm{SFR},\alpha}(0)=\Sigma_{\mathrm{SFR},0}$. An instantaneous depletion-time change instead requires the actual pre-change stellar history in equation (48); its filtered initial value need not equal the new instantaneous SFR.

To show the algebra, consider a unit-amplitude input $e^{-\gamma t}$ with initial filtered value one. Multiplication of equation (49) by $e^{t/\tau_{\mathrm{ion}}}$ and integration gives

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
\tag{50}
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
\tag{51}
$$

The coefficients add to unity, satisfying the specified initial value. The equal-gas-rate case can be evaluated from equation (48), avoiding a numerically singular difference. On Gyr-scale gas evolution, the 3-Myr smoothing is small, but it prevents an unphysical instantaneous Halpha response to a sharp SFR change.



For completeness, the equivalent ODE follows by changing the integration variable to the stellar formation time $u=t-a$:

$$
\begin{aligned}
\overline\Sigma_{\mathrm{SFR},\alpha}(t)
&=\frac{1}{\tau_{\mathrm{ion}}}\int_{-\infty}^{t}
e^{-(t-u)/\tau_{\mathrm{ion}}}\Sigma_{\mathrm{SFR}}(u)\,du,\\
\frac{d\overline\Sigma_{\mathrm{SFR},\alpha}(t)}{dt}
&=\frac{\Sigma_{\mathrm{SFR}}(t)}{\tau_{\mathrm{ion}}}
-\frac{\overline\Sigma_{\mathrm{SFR},\alpha}(t)}{\tau_{\mathrm{ion}}}.
\end{aligned}
\tag{52}
$$

The boundary term comes from the upper limit; differentiating the exponential gives the second term. Normalization ensures that constant SFR is preserved. The exponential kernel is not a 10-Myr top-hat: for $\tau_{\mathrm{ion}}=3$ Myr it assigns $1-e^{-A/\tau_{\mathrm{ion}}}$ of the weight to ages below $A$, giving 90% below 6.91 Myr and 95% below 8.99 Myr.

For equal gas rates, $\gamma_{\mathrm{HI}}=\gamma_{\mathrm{H_2}}=\gamma$, the input is $(1+t/\tau_{\Phi,0})e^{-\gamma t}$. Separating the contribution from its constant prehistory yields a finite integral without a difference of nearly equal gas rates:

$$
F_\alpha(t)=e^{-t/\tau_{\mathrm{ion}}}
\left[1+\frac{1}{\tau_{\mathrm{ion}}}
\int_0^t e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)u}
\left(1+\frac{u}{\tau_{\Phi,0}}\right)du\right].
\tag{53}
$$

The integral of $u e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)u}$ follows by integration by parts. Writing the result without an additional physical parameter,

$$
\begin{aligned}
\int_0^t u e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)u}du
&=\frac{t e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)t}}{\tau_{\mathrm{ion}}^{-1}-\gamma}
-\frac{e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)t}-1}{(\tau_{\mathrm{ion}}^{-1}-\gamma)^2}.
\end{aligned}
\tag{54}
$$

If also $\gamma=\tau_{\mathrm{ion}}^{-1}$, the integrand in equation (53) is simply $1+u/\tau_{\Phi,0}$, giving

$$
F_\alpha(t)=e^{-t/\tau_{\mathrm{ion}}}
\left[1+\frac{t}{\tau_{\mathrm{ion}}}
+\frac{t^2}{2\tau_{\mathrm{ion}}\tau_{\Phi,0}}\right].
\tag{55}
$$

For smooth SFR evolution, rearranging the ODE and substituting its leading approximation on the derivative side gives

$$
\overline\Sigma_{\mathrm{SFR},\alpha}(t)
\simeq\Sigma_{\mathrm{SFR}}(t)
-\tau_{\mathrm{ion}}\frac{d\Sigma_{\mathrm{SFR}}(t)}{dt}.
\tag{56}
$$

The leading fractional correction is $-\tau_{\mathrm{ion}}\,d\ln\Sigma_{\mathrm{SFR}}/dt$. It is small when all relevant variation timescales exceed the response time and the prehistory has relaxed. A small first derivative at a single instant is insufficient if higher derivatives or an unresolved jump are large. For the present coefficients, $\gamma_{\mathrm{HI}}\tau_{\mathrm{ion}}=0.0122$ and $\gamma_{\mathrm{H_2}}\tau_{\mathrm{ion}}=0.0009$, and the exact calculation verifies the small one-Gyr correction.

# Appendix B. Variable gas coefficients and sensitivity tests

## B.1 Exact quadrature when the coefficients evolve


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
\tag{57}
$$

The product-rule steps are identical to equations (10)--(11), but the exponent contains an integral of the rate. This quadrature solution remains exact for specified time-dependent coefficients, provided the system stays linear. The simple two-exponential form does not.

With constant positive external supply, the steady columns satisfy $\Sigma_{\mathrm{HI},\infty}=\Sigma_{\mathrm{in}}/\gamma_{\mathrm{HI}}$ and $\Sigma_{\mathrm{H_2},\infty}=\Sigma_{\mathrm{in}}/(\tau_{\mathrm{conv}}\gamma_{\mathrm{HI}}\gamma_{\mathrm{H_2}})$. At fixed external supply and conversion time, the asymptotic SFR relative to the otherwise identical unstripped equilibrium is $[1+\gamma_{\mathrm{strip}}\tau_{\mathrm{conv}}]^{-1}$. The model then approaches a nonzero level rather than necessarily quenching completely.



Even with time-dependent depletion time, the SFR remains the molecular column divided by that time. Applying the quotient rule explicitly gives

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{SFR}}(t)}{dt}
&=\frac{1}{\tau_{\mathrm{dep}}(t)}\frac{d\Sigma_{\mathrm{H_2}}(t)}{dt}
-\frac{\Sigma_{\mathrm{H_2}}(t)}{\tau_{\mathrm{dep}}^2(t)}\frac{d\tau_{\mathrm{dep}}(t)}{dt},\\
\frac{d\ln\Sigma_{\mathrm{SFR}}(t)}{dt}
&=\frac{1}{\tau_\Phi(t)}-\frac{1-R+\lambda}{\tau_{\mathrm{dep}}(t)}
-\frac{d\ln\tau_{\mathrm{dep}}(t)}{dt}.
\end{aligned}
\tag{58}
$$

For piecewise-constant coefficients, solve each interval with its own constants and carry the terminal gas columns into the next interval as initial values. Gas masses remain continuous unless an explicit transport or removal impulse is imposed. If $\tau_{\mathrm{dep}}$ jumps, instantaneous SFR can jump even at continuous gas mass; Halpha must still use the stellar response in Appendix A.

## B.2 Faster conversion as a mathematical test only

Shortening $\tau_{\mathrm{conv}}$ is not the adopted explanation of RPS compression. Its sign of change is not established by Model 0 or by the cited observations. To retain the previous algebra as a sensitivity test, if initially balanced columns are held fixed while $\tau_{\mathrm{conv}}$ is reduced, the new fractional initial slope is

$$
\frac{1}{\Sigma_{\mathrm{SFR},0}}\left.\frac{d\Sigma_{\mathrm{SFR}}(t)}{dt}\right|_{0^+}
=\gamma_{\mathrm{H_2}}
\left(\frac{\tau_{\mathrm{conv,pre}}}{\tau_{\mathrm{conv,post}}}-1\right).
\tag{59}
$$

This follows by substituting the pre-change balance into equation (17). The labels pre/post in this equation refer to an imposed coefficient change, not the observational infall categories. Reducing the example conversion time by four gives a peak at 0.18368 Gyr with $F_{\mathrm{SFR}}=1.06458$, only 0.0272 dex. This is an alternative mathematical experiment, not one of the main Figure 1 curves.

Adding the two gas balances and integrating also gives a bound at fixed depletion time:

$$
\Sigma_{\mathrm{H_2}}(t)\leq\Sigma_{\mathrm{H_2},0}+\Sigma_{\mathrm{HI},0}
\quad\Longrightarrow\quad
F_{\mathrm{SFR}}(t)\leq1+\frac{\Sigma_{\mathrm{HI},0}}{\Sigma_{\mathrm{H_2},0}}.
\tag{60}
$$

For the illustrative columns this bound is 1.2853, or 0.1090 dex. It limits conversion of the existing local gas, not compression supplied by mass convergence, and not models with changing efficiency. It must not be used to rule out all RPS-induced spatial enhancement.

# Appendix C. Separating compact HII and leaked-OB emission


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
\tag{61}
$$

$\mathcal L$ denotes luminosity per adopted area, in $\mathrm{erg\,s^{-1}\,kpc^{-2}}$. The photon-to-Halpha energy factor $\epsilon_\alpha$ is defined in section 5. Leaked photons are a redistribution of the young photon budget, not a second independent supply of OB photons. Radiation transport can make leakage nonlocal; Appendix E states the corresponding limitation.

The SFR calibration coefficient $C_\alpha$ is defined by $\Sigma_{\mathrm{SFR}}=C_\alpha\mathcal L_\alpha$ for the adopted fully absorbed young-star calibration. Consequently,

$$
\mathcal L_\alpha^{\mathrm{young}}(t)
\equiv\mathcal L_\alpha^{\mathrm{HII}}+\mathcal L_\alpha^{\mathrm{leak}}
=\frac{f_{\mathrm{young}}}{C_\alpha}
\overline\Sigma_{\mathrm{SFR},\alpha}(t).
\tag{62}
$$

Our numerical value, $C_\alpha=4.9835821\times10^{-42}\ M_\odot\,\mathrm{yr^{-1}}/(\mathrm{erg\,s^{-1}})$, follows the existing MAUVE pipeline convention checked by its saved source fingerprint. It is not assumed universal across IMFs or stellar populations. If $f_{\mathrm{young}}$ is constant, the remaining young-Halpha fraction is exactly $F_\alpha(t)$.



For a fixed allocation of young photons, define the compact-HII fraction of **young** Balmer emission,

$$
\eta_B\equiv
\frac{\mathcal L_B^{\mathrm{HII}}}
{\mathcal L_B^{\mathrm{HII}}+\mathcal L_B^{\mathrm{leak}}},
\qquad
R_{\ell/B}^{\mathrm{young}}
\equiv\eta_B R_{\ell/B}^{\mathrm{HII}}
+(1-\eta_B)R_{\ell/B}^{\mathrm{leak}}.
\tag{63}
$$

For constant component decrements and photon allocation, $\eta_B$ is constant. It is $f_{\mathrm{HII}}/f_{\mathrm{young}}$ for Halpha under equation (61). Both young components then fade in the same proportion, leaving their combined ratio constant. It is unnecessary to assume that their individual ratios are equal.



A different young effective spectrum can arise if the partition changes or the emitting gas changes. These are additional terms, not a consequence of multiplying both young luminosities by the same fading factor. Compact and leaked emission need not have equal line ratios because radiation filtering, ionization parameter, and gas conditions can differ ([Belfiore et al. 2022](#ref-belfiore)). The main derivation needs only their fixed effective young spectrum.

# Appendix D. Unequal Balmer decrements and the actual Hbeta audit


Let $\mathcal B_j=\mathcal L_\alpha^j/\mathcal L_\beta^j$ be the component Balmer decrement. From $\mathcal L_\beta^j=\mathcal L_\alpha^j/\mathcal B_j$,

$$
w_{j,\beta}
=\frac{w_{j,\alpha}/\mathcal B_j}
{\sum_iw_{i,\alpha}/\mathcal B_i}.
\tag{64}
$$

Only equal decrements make $w_{j,\beta}=w_{j,\alpha}$. We use intrinsic/de-reddened luminosities and equal $\mathcal B_j=2.86$ in the illustrative curves, consistent with a standard low-density, approximately $10^4$ K Case-B approximation. The secondary numerical test in this appendix uses the actual exported Hbeta denominator.



For common component decrements, the model weights are exactly equal even when the total luminosity changes. Real corrected data need not satisfy the equality exactly. The pipeline's dust correction adopts 2.86 but does not force every low-decrement measurement onto that value; its nonnegative-extinction handling can leave smaller corrected values. The exported means give approximately 2.8595 pre-peak and 2.8535 post-peak. We preserve them rather than replacing the Hbeta column.

Using those Hbeta luminosities directly, with $\mathcal L_\beta^{\mathrm{HOLMES}}=\mathcal L_\alpha^{\mathrm{HOLMES}}/2.86$, gives O3$_{\mathrm{post,pred}}=0.95095$ and O3$_{\mathrm{HOLMES,required}}=-4.8361$. The main common-weight calculation gives slightly different numbers but the same sign failure. The CSV `line_budget_actual_denominator.csv` stores this secondary check; `line_budget_constraints.csv` stores Model 0. Thus the simplification is explicit and its consequence is quantified.

A single dust correction to a mixed spectrum also need not recover each individually corrected component. Differential extinction and temperature-dependent decrements belong in a more general model; they are not silently introduced into the main analytical chain.

# Appendix E. Spatial transport and mixed illumination

## E.1 Full gas continuity and the sign of compression

For a thin-disc surface density, local mass conservation takes the standard Eulerian form of the continuity equation. [Armitage (2022), equations 91 and 93--97](#ref-armitage), gives the volume equation and its vertical integration; the non-axisymmetric phase-specific sources and sinks below are our explicit extension:

$$
\frac{\partial\Sigma_i(\boldsymbol{x},t)}{\partial t}
+\boldsymbol\nabla\cdot[\Sigma_i(\boldsymbol{x},t)\boldsymbol v_i(\boldsymbol{x},t)]
=\mathcal S_i(\boldsymbol{x},t)-\mathcal D_i(\boldsymbol{x},t).
\tag{65}
$$

$\boldsymbol v_i$ is the in-plane phase velocity, and $\mathcal S_i$ and $\mathcal D_i$ are local surface source and sink rates. Inserting the phase transfers and losses used in section 2 yields

$$
\begin{aligned}
\frac{\partial\Sigma_{\mathrm{HI}}}{\partial t}
&=-\boldsymbol\nabla\cdot(\Sigma_{\mathrm{HI}}\boldsymbol v_{\mathrm{HI}})
-\Sigma_\Phi-\gamma_{\mathrm{strip}}\Sigma_{\mathrm{HI}},\\
\frac{\partial\Sigma_{\mathrm{H_2}}}{\partial t}
&=-\boldsymbol\nabla\cdot(\Sigma_{\mathrm{H_2}}\boldsymbol v_{\mathrm{H_2}})
+\Sigma_\Phi-(1-R+\lambda)\Sigma_{\mathrm{SFR}}.
\end{aligned}
\tag{66}
$$

Every field in this display depends on $(\boldsymbol{x},t)$; only here the arguments are suppressed to keep the spatial conservation equations readable. Any vertical removal represented by $\gamma_{\mathrm{strip}}$ is already included in that sink and must not be counted again as a boundary loss. The transport term has a minus sign on the right. It adds column where the divergence of the mass flux is negative. Velocity convergence alone is not identical to this condition, since $\boldsymbol\nabla\cdot(\Sigma_i\boldsymbol v_i)=\boldsymbol v_i\cdot\boldsymbol\nabla\Sigma_i+\Sigma_i\boldsymbol\nabla\cdot\boldsymbol v_i$.

Integrating the transport term over a patch converts it into minus the outward mass flux through the boundary. A net inward flux can establish the initial columns in section 3.5 without shortening $\tau_{\mathrm{conv}}$ or $\tau_{\mathrm{dep}}$. Continuing transport requires a specified velocity field or another closure, so the simple closed local ODE is no longer sufficient. Model 0 neglects that later transport; it does not claim that RPS has no compression phase.

## E.2 Nonlocal photons and shared gas


A local SFR cannot necessarily predict the leakage-powered emission at that same position. A transport kernel $\mathcal T_{\mathrm{leak}}(\boldsymbol{x},\boldsymbol{x}')$, with units of inverse adopted area, can map emitted photons at $\boldsymbol{x}'$ to absorbed photons per unit area at $\boldsymbol{x}$:

$$
\mathcal Q_{\mathrm{abs,leak}}(\boldsymbol{x},t)
=\int \mathcal T_{\mathrm{leak}}(\boldsymbol{x},\boldsymbol{x}')
\mathcal Q_{\mathrm{OB}}(\boldsymbol{x}',t)\,dA'.
\tag{67}
$$

Its area integral must respect the available escaped-photon fraction. Neighbouring young populations can maintain diffuse emission while a local patch fades, breaking the assumption of a common local $F_\alpha$. Similarly, old-star photons can propagate away from their birth positions. The local photon benchmark in section 8 excludes these transfers.

In a gas parcel illuminated simultaneously by OB stars and HOLMES, absorbed photon rates may be budgeted by source, but forbidden-line emissivities depend on the jointly determined ionic fractions and temperature. There is no unique decomposition of each collisionally excited line into an OB part and a HOLMES part independent of that solution. The fundamental prediction then takes the form

$$
R_{\ell/B}=\mathscr R_{\ell/B}
\left(\mathrm{SED}_{\mathrm{OB}}+\mathrm{SED}_{\mathrm{HOLMES}},
U,Z,\mathrm{N/O},n_{\mathrm H},\ldots\right),
\qquad U\equiv\frac{\Phi_{\mathrm H}}{n_{\mathrm H}c}.
\tag{68}
$$

$\mathscr R$ denotes the result of a photoionization calculation, not a fitted analytic function in this report; the SEDs include their radiation normalizations. $\Phi_{\mathrm H}$ is incident ionizing photon flux, $n_{\mathrm H}$ hydrogen density, and $U$ the dimensionless ionization parameter. Metallicity $Z$, abundance ratio N/O, gas density, geometry, and dust also affect the result. Holding every intrinsic line ratio constant is a controlled approximation to be tested, especially as the radiation field fades. It should not be defended as an exact consequence of a fixed stellar spectrum.



# Appendix F. Recombination details and differential fading

## F.1 Recovering the HOLMES normalization from the emission measure


In photoionization equilibrium, the absorbed hydrogen-ionizing photon rate equals the Case-B recombination rate. Let $f_{\mathrm{abs,HOLMES}}$ denote the fraction of the HOLMES photons assigned to the region that actually ionize its hydrogen. For a region with adopted physical area $A_{\mathrm{reg}}$,

$$
f_{\mathrm{abs,HOLMES}}\mathcal Q_{\mathrm{HOLMES}}A_{\mathrm{reg}}
=\alpha_B\int_{V_{\mathrm{reg}}}n_e n_p\,dV.
\tag{69}
$$

Here $n_e$ and $n_p$ are electron and proton densities, $\alpha_B$ is the recombination coefficient excluding direct recombinations to the ground state, and $V_{\mathrm{reg}}$ is the emitting volume. The area and volume must be converted to a consistent unit system. The emitted Halpha luminosity is

$$
L_\alpha^{\mathrm{HOLMES}}
=h_{\mathrm P}\nu_\alpha\alpha_\alpha^{\mathrm{eff}}
\int_{V_{\mathrm{reg}}}n_e n_p\,dV,
\tag{70}
$$

where $h_{\mathrm P}$ is Planck's constant, $\nu_\alpha=c/\lambda_\alpha$ is the line frequency, and $\alpha_\alpha^{\mathrm{eff}}$ counts recombinations that generate Halpha photons. Case-B coefficients depend on temperature and density ([Hummer & Storey 1987](#ref-hs)). 

Equations (69) and (70) contain the same volume integral. Solving the former for that integral, substituting into the latter, and dividing by $A_{\mathrm{reg}}$ gives equation (30). The area must have the same units on both sides before converting to kpc squared.


Taking the logarithmic derivative of equation (30) exposes all conditions:

$$
\frac{d\ln\mathcal L_\alpha^{\mathrm{HOLMES}}}{dt}
=\frac{d\ln\epsilon_\alpha}{dt}
+\frac{d\ln f_{\mathrm{abs,HOLMES}}}{dt}
+\frac{d\ln q_{\mathrm{H,HOLMES}}}{dt}
+\frac{d\ln\Sigma_*^{\mathrm{old}}}{dt}.
\tag{71}
$$

For a short interval relative to the age of an old population, its mass and specific photon output may vary slowly. We additionally assume nearly fixed absorption and recombination conditions in the retained inner gas. Only with these assumptions do we set $\mathcal L_\alpha^{\mathrm{HOLMES}}(t)\simeq\mathcal L_{\alpha,0}^{\mathrm{HOLMES}}$.

The gas must also be able to absorb the available photons. Equation (69) implies a necessary emission measure. For example, if a maximum available ionized column is specified, its recombination capacity cannot be exceeded merely by increasing $q\Sigma_*^{\mathrm{old}}$. Gas removal can reduce the covering fraction and capacity, causing the HOLMES-powered emission to fall. The recombination time is approximately $\tau_{\mathrm{rec}}=(\alpha_B n_e)^{-1}$ at fixed density; sustained emission over longer intervals requires sustained ionization, not a remnant afterglow.

The constant term is therefore particularly relevant to **retained inner gas**. It is not a reason to predict a permanent Halpha floor in a completely stripped outer region.


## F.2 Fixed spectra: luminosity can fall while the ratio rises


The total forbidden-line luminosity is

$$
\mathcal L_\ell(t)
=R_{\ell/B}^{\mathrm{young}}\mathcal L_{B,0}^{\mathrm{young}}F_{\mathrm{SFR}}(t)
+R_{\ell/B}^{\mathrm{HOLMES}}\mathcal L_{B,0}^{\mathrm{HOLMES}}.
\tag{72}
$$

For positive ratios, its derivative is $R_{\ell/B}^{\mathrm{young}}\mathcal L_{B,0}^{\mathrm{young}}\,dF_{\mathrm{SFR}}/dt<0$. Its fractional fading is smaller in magnitude than the Balmer fractional fading if its HOLMES-to-young ratio contrast is larger. Eliminating $F_{\mathrm{SFR}}$ between the total Balmer luminosity and equation (72) gives

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
\tag{73}
$$

For the higher-ratio HOLMES branch, the last slope lies between $-1$ and $0$. This supplies a directly testable ratio--surface-brightness relation at fixed old stellar density and fixed templates. It also clarifies what differential fading means here: the distinct components fade differently; each component's internal line ratios are held constant.

The HOLMES fraction of a forbidden line can exceed its Balmer fraction. Explicitly,

$$
\frac{\mathcal L_\ell^{\mathrm{HOLMES}}}{\mathcal L_\ell}
=\frac{w_{\mathrm{HOLMES},B}R_{\ell/B}^{\mathrm{HOLMES}}}
{(1-w_{\mathrm{HOLMES},B})R_{\ell/B}^{\mathrm{young}}
+w_{\mathrm{HOLMES},B}R_{\ell/B}^{\mathrm{HOLMES}}}.
\tag{74}
$$

At a Balmer fraction of 0.02, a forbidden-line ratio contrast greater than 49 would make the HOLMES component dominate that line. This is our algebraic illustration, not an extraction of a contrast from Belfiore. In a shared gas volume, even this source-by-source attribution of forbidden-line cooling is not unique; the mixed incident spectrum changes the gas state jointly.


## F.3 Evolving populations, absorption, and intrinsic spectra


Without setting the old-star luminosity constant, differentiation of its weight gives

$$
\frac{dw_{\mathrm{HOLMES},B}}{dt}
=w_{\mathrm{HOLMES},B}(1-w_{\mathrm{HOLMES},B})
\left[
\frac{d\ln\mathcal L_B^{\mathrm{HOLMES}}}{dt}
-\frac{d\ln\mathcal L_B^{\mathrm{young}}}{dt}
\right].
\tag{75}
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
\tag{76}
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
\tag{77}
$$

These equations separate three effects: changing old-versus-young luminosity weights, changing intrinsic spectra, and changing compact-versus-diffuse allocation of young photons. The fixed-spectrum baseline is useful because it isolates the first. If it fails, equations (76)--(77) show which additional physical terms might matter, but do not determine them from one ratio.


# Appendix G. Complete notation and parameter ledger

**Table G1. Gas and temporal quantities.** Initial values carry subscript 0; pre/post on conversion times refer to the imposed change at model onset, not directly to observed infall-stage membership.

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

**Table G2. Ionizing populations and line emission.** A calligraphic luminosity is per adopted area; ordinary $L$ is an integrated luminosity.

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

**Table G3. Recombination and spectral quantities.**

| Symbol | Definition and units |
|:--|:--|
| $h_{\mathrm P}$, $c$, $\nu_\alpha$, $\lambda_\alpha$ | Planck's constant; speed of light; Halpha frequency; Halpha wavelength, 6562.8 Angstrom in the calculation. |
| $\alpha_B$, $\alpha_\alpha^{\mathrm{eff}}$ | Total Case-B and effective Halpha recombination coefficients, $\mathrm{cm^3\,s^{-1}}$. |
| $p_\alpha$, $\epsilon_\alpha$ | Halpha photons per absorbed hydrogen-ionizing photon; Halpha energy per such absorbed photon. Numerical values 0.45331 and $1.3721\times10^{-12}$ erg. |
| $n_e$, $n_p$, $n_{\mathrm H}$ | Electron, proton, and total hydrogen number densities, $\mathrm{cm^{-3}}$. |
| $A_{\mathrm{reg}}$, $V_{\mathrm{reg}}$, $\tau_{\mathrm{rec}}$ | Adopted physical area, emitting volume, and approximate recombination time. |
| $R_{\ell/B}^j$, $R_{\ell/B}$ | Component and total linear forbidden-line-to-Balmer ratios; dimensionless. |
| $w_{j,B}$, $w_{\mathrm{target}}$ | Component fraction of Balmer luminosity; chosen target weight. These are not area fractions. |
| $\mathcal B_j$, $\eta_B$ | Component Halpha/Hbeta decrement; compact-HII share of young Balmer emission. |
| N2, S2, O3 | [N II]6583/Halpha; [S II]6716+6731/Halpha; [O III]5007/Hbeta. |
| $\mathcal T_{\mathrm{leak}}$, $\mathcal Q_{\mathrm{abs,leak}}$ | Photon transport kernel, inverse area; leakage photon absorption rate per area. |
| $\mathscr R$, SED, $U$, $\Phi_{\mathrm H}$ | Photoionization prediction; spectral energy distribution; ionization parameter; incident photon flux in $\mathrm{s^{-1}\,cm^{-2}}$. |
| $Z$, N/O | Gas metallicity and nitrogen-to-oxygen abundance ratio; the calculation does not evolve them. |

**Table G4. Observational and mathematical conventions.**

| Term | Meaning |
|:--|:--|
| SF, NSF, ND | The source notebook's star-forming, other Balmer-detected, and joint Balmer non-detection categories. |
| RPS, DIG, OB, AGN | Ram-pressure stripping; diffuse ionized gas; massive O/B stellar population; active galactic nucleus. |
| BPT, EW, S/N, IMF | Standard optical excitation diagrams; equivalent width; signal-to-noise; initial mass function. |
| $\ln$, $\log_{10}$, dex | Natural logarithm, base-ten logarithm, and a base-ten logarithmic difference. |
| $O(t^3)$ | Terms cubic or higher in the short-time expansion, not a new physical parameter. |
| Superscripts obs, pred, required | Observed mean, conditional prediction, and an algebraic value required under the stated assumptions. |


**Table G5. Additional quantities made explicit in this revision.**

| Symbol | Definition and units |
|:--|:--|
| $c_{\mathrm{HI}}$, $c_{\mathrm{H_2}}$ | Initial gas-column factors relative to a reference, dimensionless; the spatial example uses 1.5 for each. |
| $\Delta_{\mathrm{spatial}}\log_{10}\Sigma_{\mathrm{SFR}}$ | Logarithmic spatial SFR excess, in dex, not a temporal derivative. |
| $w_{\mathrm{HOLMES}}$, $w_{\mathrm{young}}$ | Common Halpha/Hbeta weights in Model 0; the former is $w_{\mathrm{HOLMES},B}$ when both component decrements equal 2.86. |
| $\mathcal Q_{\mathrm{abs,HOLMES}}$ | HOLMES photons absorbed by hydrogen per second per adopted area. |
| $\boldsymbol v_i$, $\mathcal S_i$, $\mathcal D_i$ | Phase velocity; surface source and sink rates. Velocity units must match the chosen length/time units. |
| $\boldsymbol\nabla$, $\partial/\partial t$ | In-plane spatial derivative; Eulerian time derivative at fixed position. |
| $\gamma$ in Appendix A | Rate of one exponential input; specializes to the common gas rate when $\gamma_{\mathrm{HI}}=\gamma_{\mathrm{H_2}}$. |
| $F_{\mathrm{required}}$ | Remaining young luminosity/SFR fraction needed for the specified endpoint comparison under Model 0. |


# Appendix H. Provenance, revision coverage, and verification

The source report is `20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.md`. The revision uses its existing analytical derivations and the saved observational export, not a fresh extraction from line maps. The [linked discussion](https://chatgpt.com/c/6abef07c-81f4-83ec-852a-a4e180a312ba) was read from its opening PDF-review request through the final 24-point revision summary: ten user/assistant exchanges. The first retrieval contained only the latest five exchanges; the earlier turns and the full final summary were subsequently read in the browser before this revision was written.

**Table H1. Where the requested changes are implemented.**

| Discussion request | Implementation |
|:--|:--|
| Explicit time dependence and constant Model 0 parameters | Section 2 and the gas derivation |
| Detailed HI-to-H2-to-SFR algebra, including equal rates | Sections 3.1--3.4 |
| Spatial enhancement from initial gas accumulation | Section 3.5; new Figure 1 and section 8.3 |
| Direct Halpha conversion in the main text | Section 4; response derivation in Appendix A |
| Compact, normalized HOLMES term and absorption conditions | Section 5; emission-measure details in Appendix F |
| Young + HOLMES and common Balmer weights | Sections 6--7; HII/leak and unequal decrements in Appendices C--D |
| Preserve the quantitative failure tests | Sections 8.5--8.8 and 9 |
| Keep extensions separate | Appendices B, E, and F |

Literature provenance is explicit at the point of use. The regulator bookkeeping and $\Sigma_\Phi$ notation follow Lilly et al. and Huang et al.; the linear atomic transfer law and the closed two-reservoir solution are assumptions and derivations of this report. Halpha calibration physics follows Kennicutt & Evans; the adopted numerical coefficient follows the MAUVE pipeline. Case B, the $1/2.206$ conversion, and the fiducial current-mass HOLMES photon yield have separate references. Linear luminosity mixing is algebra, with Blanc et al. as a related empirical construction. The new Brown citation was checked against the primary paper, including the size and definition of its early-stage subset. Belfiore section 3.2/footnote 5 and Cid Fernandes equation 2 were freshly checked for this revision; the remaining reference framework is retained from the source report, without claiming a new full literature audit.

The numerical script is `assets/20261002_Model0_Derivation/model_predictions.py`. It rebuilds the same equal-galaxy observational anchors, checks the saved source fingerprints, independently integrates the gas and response equations, evaluates the new spatial initial-state example, and calculates both common-weight and actual-Hbeta inversions. It writes the numerical audit, CSVs, and three figures. Run it with:

```bash
MPLCONFIGDIR=/private/tmp/mauve_20261002/mpl \
/opt/miniconda3/envs/ICRAR/bin/python \
  /Users/Igniz/Desktop/ICRAR/MAUVE/assets/20261002_Model0_Derivation/model_predictions.py
```


The executed gas-versus-ODE checks have maximum relative difference 7.68e-14; the filtered-response check has maximum absolute normalized difference 3.33e-11. The equal-gas-rate check differs by 2.44e-15; the scaled-initial-state ODE check differs by 2e-15. Peak, positive-supply, total-mass, and unequal-Balmer-weight identities pass for the implemented examples. These numerical tolerances test the equations, not the physical assumptions.


The observational input is `assets/20260914_resolved_RPS_academic_model/stage_bpt_line_profiles_with_hbeta.csv`. Its fingerprint and the six recorded source fingerprints are checked and saved in the new asset directory. This validates reuse of that extraction record, not the present contents of every large FITS map. No full map pipeline, bootstrap, fitted stellar population, new CO/HI analysis, radiation transport calculation, or photoionization grid was executed. No significance is assigned to the cross-sectional amplitude tests. Figure colors identify model cases, not galaxies.

The PDF is generated from the same Markdown and checked for equation numbering, internal links, text boundaries, and page rendering. Detailed acceptance evidence and the revision checklist are saved with the assets. The source report and other user files are preserved.


# References

<span id="ref-armitage"></span>
**Armitage, P. J. (2022).** *Lecture notes on accretion disk physics.* arXiv:2201.07262, sections II.A and III.A.1. [Primary manuscript](https://arxiv.org/pdf/2201.07262). Equations 9 and 91--97 define surface density and derive its continuity equation; the phase-specific sources and sinks are added explicitly in this report.

<span id="ref-belfiore"></span>
**Belfiore, F., et al. (2022).** *A tale of two DIGs: The relative role of H II regions and low-mass hot evolved stars in powering the diffuse ionised gas in PHANGS-MUSE galaxies.* A&A, 659, A26. [DOI](https://doi.org/10.1051/0004-6361/202141859); [primary full text](https://arxiv.org/html/2111.14876v3).


<span id="ref-blanc"></span>
**Blanc, G. A., Heiderman, A., Gebhardt, K., Evans, N. J. II, & Adams, J. (2009).** *The Spatially Resolved Star Formation Law from Integral Field Spectroscopy: VIRUS-P Observations of NGC 5194.* ApJ, 704, 842--862. [DOI](https://doi.org/10.1088/0004-637X/704/1/842); [primary manuscript](https://arxiv.org/pdf/0908.2810).


<span id="ref-boselli"></span>
**Boselli, A., et al. (2014).** *Cold gas properties of the Herschel Reference Survey. III. Molecular gas stripping in cluster galaxies.* A&A, 564, A67. [DOI](https://doi.org/10.1051/0004-6361/201322313); [primary manuscript](https://arxiv.org/abs/1402.0326).


<span id="ref-brown"></span>
**Brown, T., et al. (2023).** *VERTICO VII: Environmental Quenching Caused by Suppression of Molecular Gas Content and Star Formation Efficiency in Virgo Cluster Galaxies.* ApJ, 956, 37. [DOI](https://doi.org/10.3847/1538-4357/acf195); [primary manuscript](https://arxiv.org/pdf/2308.10943).


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

