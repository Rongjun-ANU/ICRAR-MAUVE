---
title: "The local response of gas, star formation, and emission lines during Virgo infall"
subtitle: "A spatially resolved derivation and worked MAUVE-scale example"
date: "29 September 2026"
author: "MAUVE-MUSE / ICRAR research note"
---

# Setup of our spatially resolved model

Here we want to understand whether one *local* gas-regulator calculation can accommodate three features seen in MAUVE: reduced star formation in the surviving SF regions, a possible small local enhancement, and fainter lines but larger forbidden/Balmer ratios in some NSF regions. We therefore proceed in the same order as the physical chain: write the H I and H2 continuity equations, solve for both reservoirs and the implied SFR, convert the SFR history into an ionizing-photon rate, divide those photons among emitting gas components, and finally calculate each line luminosity before taking ratios. We show the intermediate algebra whenever it changes the interpretation of a result.

All gas and stellar quantities are written per projected area of a fixed patch at disc position \(\boldsymbol{x}=(r,\varphi)\). We suppress \(\boldsymbol{x}\) in the intermediate equations only because the following ODEs contain no spatial derivative. The variable \(t\) means elapsed time after an assumed local perturbation; it is **not** an age assigned to a MAUVE infall-stage bin. Thus a time series from this model can suggest a mechanism, but three cross-sectional stage samples do not trace the same patch through time.

We address the questions in the [subsequent discussion](https://chatgpt.com/c/6aa7eb87-6498-83ec-93f7-fa5189f7274f): which assumption yields exponentials, how to distinguish instantaneous from cumulative suppression, whether an H-alpha age-response kernel matters, why *individual* forbidden lines can fade while their ratios rise, whether setting compact-gas photon absorption to unity is consistent with leakage, and when the intrinsic Balmer ratio 2.86 is appropriate. The linked discussion supplies questions, not independent scientific evidence. At each stage we distinguish a conservation identity, a literature-based relation, and a closure chosen only for the worked example.

# 1. Empirical scales and relevant literature

## 1.1 What the numerical example is anchored to

The scalar line products were re-extracted on 14 September from the live MAUVE maps under the strict post-fit continuum cut \(S/N_{\mathrm{POSTFIT}}>25\). We rechecked the recorded fingerprints of the three stage-analysis notebooks and the SFR pipeline on 29 September; their large FITS input maps were not individually rehashed. Table 1 uses a representative common-support stellar-density bin centred on \(\log_{10}(\Sigma_*/M_\odot\,{\mathrm{kpc}}^{-2})=8.625\). Each line surface luminosity is corrected by the pipeline and averaged according to its galaxy-level estimator on the **same five-line detected support**, then ratios are taken from those mean line luminosities. Values are rounded; they are not medians of individual-pixel ratios.

**Table 1.** MAUVE common-support scales used to anchor and challenge the example.

| Stage and observed class | Supported galaxy products | Mean H-alpha surface luminosity (erg s\(^{-1}\) kpc\(^{-2}\)) | [N II]/H-alpha | [S II]/H-alpha | [O III]/H-beta |
|:--|--:|--:|--:|--:|--:|
| Pre-peak SF | 7 | \(1.56\times10^{40}\) | 0.254 | 0.222 | 0.596 |
| Pre-peak NSF | 6 | \(3.56\times10^{39}\) | 0.373 | 0.377 | 0.898 |
| Post-peak SF | 12 | \(7.58\times10^{39}\) | 0.319 | 0.228 | 0.210 |
| Post-peak NSF | 13 | \(1.11\times10^{39}\) | 0.566 | 0.417 | 0.753 |

Here **SF** and **NSF** are the *observational classes* defined by the pipeline, not pure physical H II and non-H II emission components. SF additionally requires the pipeline's H II selection, H-alpha equivalent width above 6 Angstrom, intrinsic width below 45 km s\(^{-1}\), and finite H II SFR; NSF is Balmer-detected but fails the joint SF definition. Nondetected (ND) regions have a different selection. The representative bin is not a universal line-ratio template; the 14 September report contains the full mass/radius profiles, mask definitions, and whole-galaxy bootstrap comparisons. The numerical experiment below uses the first and fourth rows as *scale anchors*. They are different sets of galaxies and regions, so their difference is not an observed time derivative. In particular, a similar model value at the fourth row cannot establish that the first row physically evolved into it.

The pre-peak systems are retained as a **field-like control working assumption** even though they belong to Virgo. The use of this internal control does not prove that their gas is environmentally pristine. The measured post-peak surviving-SF intensity factor of approximately \(0.407\), with a \(0.341\)–\(0.504\) whole-galaxy interval, remains a separate descriptive fit from the earlier report. The present transformed-region example should not be equated with that surviving-SF estimator.

The updated 28 September execution of the specifically named NGC4654 gradient notebook still gives a facing-side median environmental residual of \(+0.0354\) dex against the pre-peak control and \(-0.0957\) dex on the opposite side. The side difference is \(0.1311\) dex. The fitted plane has \(R^2=0.0444\); neither its position-angle uncertainty nor a galaxy-level detection significance has been established. Thus NGC4654 motivates a possible **small local positive response** but does not calibrate the efficiency pulse introduced below.

## 1.2 What previous work justifies, and what it does not

The local molecular regulator and the symbols \(\Sigma_\Phi\), \(\tau_\Phi\), and \(\tau_{\mathrm{dep}}\) follow [Huang et al. (2026)](#ref-huang), with wider gas-regulator context from [Lilly et al. (2013)](#ref-lilly). [Köppen et al. (2018)](#ref-koppen) derive long- and short-pulse ram-pressure regimes; [Singh et al. (2019)](#ref-singh) calculate radial stripping in an idealized orbiting disc. These works motivate position- and time-dependent gas loss, but provide no universal local conversion from \(P_{\mathrm{ram}}\) to the fractional H I or H2 removal coefficient used here. The resolved Virgo comparison of [Brown et al. (2023)](#ref-brown) finds both molecular-content and efficiency changes, so an SFR change cannot automatically be assigned to only one of them.

[Lizee et al. (2021)](#ref-lizee) report a compressed NGC4654 region with a model-dependent molecular-efficiency increase of roughly 1.5–2 at approximately kpc resolution; the present optical element and its modest facing-side control offset are different measurements. [Vollmer et al. (2012)](#ref-vollmer12) and [Nehlig et al. (2016)](#ref-nehlig) show why environmental compression need not increase molecular efficiency in every region. The efficiency pulse below is therefore a conditional demonstration, not a direct RPS law.

For ionized gas, [Kennicutt & Evans (2012)](#ref-ke12) give the H-alpha tracer calibration and its short young-stellar age response. [Hummer & Storey (1987)](#ref-hs87) and [Storey & Hummer (1995)](#ref-sh95) give the Case-B recombination framework and effective hydrogen emissivities. [Belfiore et al. (2022)](#ref-belfiore22) find diffuse emission powered by both leaked young-star photons and hot evolved stars, with different importance for H-alpha and [O III]; [Belfiore et al. (2017)](#ref-belfiore16) emphasize that emitted ionizing photons make Balmer lines only if gas absorbs them. [Zhang et al. (2017)](#ref-zhang) show that diffuse emission changes diagnostic ratios. None of these sources establishes a universal ordering in which **all** non-H II emission fades more slowly than H II emission. That ordering must emerge from an identified source and gas-absorption history, or be stated as a conditional model result.

# 2. Two neutral reservoirs: definitions, solution, and quenching

## 2.1 The physical aperture and a restricted closure

Let \(\Sigma_{\mathrm{HI}}\) and \(\Sigma_{\mathrm{H_2}}\) be local atomic-associated and molecular-associated gas surface densities in \(M_\odot\,{\mathrm{pc}}^{-2}\), under one consistent helium convention. Let \(\Sigma_{\mathrm{SFR}}\) be the true formed-stellar-mass rate per area in \(M_\odot\,{\mathrm{yr}}^{-1}\,{\mathrm{kpc}}^{-2}\). All gas-rate equations below use one compatible area/time convention internally; the conversion \(1\,M_\odot\,{\mathrm{pc}}^{-2}\,{\mathrm{Gyr}}^{-1}=10^{-3}\,M_\odot\,{\mathrm{yr}}^{-1}\,{\mathrm{kpc}}^{-2}\) is applied in the numerical calculation.

For a fixed aperture, local phase continuity first states that the rate of change of each reservoir equals inputs minus outputs plus any net flux across the boundary. The numerical closure neglects **in-plane exchange with adjacent resolved patches**, external atomic supply *after* \(t=0\), reverse molecular dissociation, and reaccretion; environmental loss across the aperture boundary is instead represented by an effective stripping sink. Those omitted processes can be restored as explicit terms; their absence is not implied by spatial resolution. With \(\tau_{\mathrm{conv}}\) the **effective transfer time from the modelled H I-associated reservoir to H2**, and \(\gamma_{{\mathrm{strip}},i}\) the environmental fractional removal rate of reservoir \(i\in\{{\mathrm{HI}},{\mathrm{H_2}}\}\), define

$$
\Sigma_\Phi=\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}},
\qquad
\tau_\Phi=\frac{\Sigma_{\mathrm{H_2}}}{\Sigma_\Phi}
=\tau_{\mathrm{conv}}\frac{\Sigma_{\mathrm{H_2}}}{\Sigma_{\mathrm{HI}}}
\quad(\Sigma_\Phi>0).
\tag{1}
$$

\(\Sigma_\Phi\) is a **rate surface density**, with the units of \(\Sigma_{\mathrm{SFR}}\), and \(\tau_\Phi\) is a **time**. This preserves [Huang et al. (2026), equations 16–17 and 21](#ref-huang); \(\tau_{\mathrm{conv}}\) is an extra closure parameter introduced here and is neither \(\Sigma_\Phi\) nor a laboratory H2 formation timescale. When the atomic supply tends to zero, \(\tau_\Phi\) tends to infinity. The supply-time restriction appropriate to the prior paper's fitted regime is not imposed on this declining-reservoir experiment.

Define molecular depletion time by the constitutive law

$$
\Sigma_{\mathrm{SFR}}
=\frac{\Sigma_{\mathrm{H_2}}}{\tau_{\mathrm{dep}}},
\qquad
\tau_{\mathrm{dep}}=\frac{\Sigma_{\mathrm{H_2}}}{\Sigma_{\mathrm{SFR}}}
\quad(\Sigma_{\mathrm{SFR}}>0),
\tag{2}
$$

where the area and time units must first be made consistent. This is the standard molecular-regulator definition used by [Huang et al. (2026)](#ref-huang) and the resolved gas comparison of [Brown et al. (2023)](#ref-brown). Here \(R\) is the *effective promptly recycled fraction assigned back to the modelled cold reservoir*, and \(\eta\) is the dimensionless feedback mass-loading factor that removes \(\eta\Sigma_{\mathrm{SFR}}\) from that reservoir. Real stellar ejecta need not enter H2 immediately; that approximation is declared rather than hidden.

An incident ram-pressure scale is

$$
P_{\mathrm{ram}}(t)=\rho_{\mathrm{ICM}}(t)v_{\mathrm{rel}}^{\,2}(t).
\tag{3}
$$

The density \(\rho_{\mathrm{ICM}}\), relative speed \(v_{\mathrm{rel}}\), wind orientation, and gas binding vary along the orbit. [Köppen et al. (2018)](#ref-koppen) and [Singh et al. (2019)](#ref-singh) motivate a time-dependent forcing. Equation (3) alone does **not** determine either \(\gamma_{{\mathrm{strip}},{\mathrm{HI}}}\) or \(\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}\); those have units inverse time and summarize net removal conditional on the local geometry. No pressure-to-loss-rate mapping is asserted in this report.

**Step 1: write the balances.** With the above closure,

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{HI}}}{dt}
&=-\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}}
  -\gamma_{{\mathrm{strip}},{\mathrm{HI}}}(t)\Sigma_{\mathrm{HI}},\\
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
&=+\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}}
  -(1-R+\eta)\Sigma_{\mathrm{SFR}}
  -\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}(t)\Sigma_{\mathrm{H_2}}.
\end{aligned}
\tag{4}
$$

The transfer term enters once as an H I loss and once as an H2 gain. The stellar and feedback term applies once to H2. Environmental stripping is separate from feedback and can differ by phase. These are *new restricted source/loss closures* around the established regulator budget, not literal rate laws taken from the stripping papers. Both phases use the **same** scientifically named notation \(\gamma_{{\mathrm{strip}},i}\). We do not introduce an unrelated \(k_1/k_2\) pair or separate Greek alphabets for their total rates.

## 2.2 The instantaneous suppression and enhancement test

**Step 2: differentiate the star-formation law.** For positive \(\Sigma_{\mathrm{H_2}}\) and \(\tau_{\mathrm{dep}}\), equation (2) gives

$$
\frac{d\ln\Sigma_{\mathrm{SFR}}}{dt}
=\frac{1}{\Sigma_{\mathrm{H_2}}}\frac{d\Sigma_{\mathrm{H_2}}}{dt}
 -\frac{d\ln\tau_{\mathrm{dep}}}{dt}.
\tag{5}
$$

Now we substitute the molecular continuity equation into the *first* term of equation (5). To show which assumption enters at which line, we also substitute the depletion definition only after writing the gas terms:

$$
\begin{aligned}
\frac{d\ln\Sigma_{\mathrm{SFR}}}{dt}
&=\frac{1}{\Sigma_{\mathrm{H_2}}}
\left[
\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}}
-(1-R+\eta)\Sigma_{\mathrm{SFR}}
-\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}\Sigma_{\mathrm{H_2}}
\right]
-\frac{d\ln\tau_{\mathrm{dep}}}{dt}\\
&=\frac{\Sigma_\Phi}{\Sigma_{\mathrm{H_2}}}
-\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
-\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}
-\frac{d\ln\tau_{\mathrm{dep}}}{dt}\\
&=\frac{1}{\tau_\Phi}
-\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
-\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}
-\frac{d\ln\tau_{\mathrm{dep}}}{dt}.
\end{aligned}
\tag{6}
$$

Every term has units inverse time. The four terms represent supply per molecular mass, net stellar locking plus feedback, direct H2 removal, and change in molecular efficiency. Equation (6) is an algebraic consequence of the [Huang et al. molecular regulator](#ref-huang) under equation (4)'s declared losses. It applies locally even if its coefficients depend on position and time.

Define the **instantaneous logarithmic SFR decline rate** and its **integrated attenuation**:

$$
\begin{aligned}
\Gamma_{\mathrm{SFR}}(t)
&\equiv-\frac{d\ln\Sigma_{\mathrm{SFR}}}{dt}
=\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
 +\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}
 +\frac{d\ln\tau_{\mathrm{dep}}}{dt}
 -\frac{1}{\tau_\Phi},\\
\mathcal E(t)
&\equiv\int_0^t\Gamma_{\mathrm{SFR}}(t')\,dt'
=\ln\!\left[\frac{\Sigma_{\mathrm{SFR}}(0)}{\Sigma_{\mathrm{SFR}}(t)}\right].
\end{aligned}
\tag{7}
$$

\(\Gamma_{\mathrm{SFR}}\) is an inverse time; \(\mathcal E\) is dimensionless. Both are **defined from the SFR**, not separately observed physical loss coefficients. Integration gives

$$
\frac{\Sigma_{\mathrm{SFR}}(t)}{\Sigma_{\mathrm{SFR}}(0)}
=e^{-\mathcal E(t)}.
\tag{8}
$$

Thus \(\Gamma_{\mathrm{SFR}}>0\) means SFR is falling **now**, \(\Gamma_{\mathrm{SFR}}<0\) means it is rising **now**, while \(\mathcal E>0\) means it is below its value at \(t=0\) and \(\mathcal E<0\) means it remains above that value. A region can have \(\Gamma_{\mathrm{SFR}}>0\) but \(\mathcal E<0\): it is already declining after a burst, yet remains enhanced relative to its starting level. This distinction directly resolves the instantaneous-versus-accumulated question in the linked discussion.

The exact local decline condition is

$$
\boxed{\;
\frac{1}{\tau_\Phi}
<
\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
+\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}
+\frac{d\ln\tau_{\mathrm{dep}}}{dt}
\;}
\quad\Longleftrightarrow\quad
\Gamma_{\mathrm{SFR}}>0.
\tag{9}
$$

Reversing the inequality gives instantaneous enhancement; equality is instantaneous stationarity. For constant \(\tau_{\mathrm{dep}}\), the final term vanishes. An increasing depletion time suppresses SFR even if molecular mass stays nonzero. A decreasing depletion time can briefly make the right-hand side small enough for a positive response. Here “quenching” describes a *local falling SFR*, not an observational class boundary, a claim of zero SF, or a uniquely identified physical cause.

## 2.3 The exact constant-coefficient solution, line by line

Constant \(\tau_{\mathrm{conv}},\tau_{\mathrm{dep}},R,\eta,\gamma_{{\mathrm{strip}},i}\) over a *chosen interval* permit a closed solution. They do **not** claim that Virgo ram pressure is constant throughout an orbit. To keep both phase equations in one notation, define the two **composite loss rates**

$$
\gamma_{{\mathrm{loss}},{\mathrm{HI}}}
=\frac{1}{\tau_{\mathrm{conv}}}+\gamma_{{\mathrm{strip}},{\mathrm{HI}}},
\qquad
\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}
=\frac{1-R+\eta}{\tau_{\mathrm{dep}}}+\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}.
\tag{10}
$$

The index \(i\) specifies the reservoir. Both \(\gamma_{{\mathrm{loss}},i}\) are just **sums of already defined rates**; they add no free parameter. The H I sum includes transfer into H2; the H2 sum includes consumption/feedback. They must therefore not both be interpreted as environmental stripping.

**Step 3: integrate H I.** From \(d\Sigma_{\mathrm{HI}}/dt=-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}\Sigma_{\mathrm{HI}}\), divide by the positive reservoir and integrate:

$$
\int_{\Sigma_{{\mathrm{HI}},0}}^{\Sigma_{\mathrm{HI}}(t)}
\frac{d\Sigma_{\mathrm{HI}}}{\Sigma_{\mathrm{HI}}}
=-\int_0^t\gamma_{{\mathrm{loss}},{\mathrm{HI}}}\,dt'
\quad\Longrightarrow\quad
\Sigma_{\mathrm{HI}}(t)
=\Sigma_{{\mathrm{HI}},0}e^{-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}t}.
\tag{11}
$$

The simple \(e^{-\gamma t}\) form requires \(\gamma\) to be constant over the solved interval. The exponential itself does not imply that the real pressure or stripping rate is constant.

**Step 4: solve H2 by an integrating factor.** Let us first return to the *second* line of equation (4). Equation (2) converts its star-formation term into a coefficient multiplying \(\Sigma_{\mathrm{H_2}}\). Equation (10) gives that whole coefficient a physical name, and only then do we insert the H I solution:

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
&=\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}}
-\left[
\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
+\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}
\right]\Sigma_{\mathrm{H_2}}\\
&=\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}}
-\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}\Sigma_{\mathrm{H_2}}\\
\Longleftrightarrow\quad
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
+\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}\Sigma_{\mathrm{H_2}}
&=\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}t}.
\end{aligned}
\tag{12}
$$

The third line is a linear inhomogeneous ODE: the surviving H I supplies H2, while the existing H2 disappears at its own composite loss rate. We seek a factor \(I(t)\) for which \(I\,d\Sigma_{\mathrm{H_2}}/dt+I\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}\Sigma_{\mathrm{H_2}}\) is the derivative of one product. Under the **constant-coefficient assumption of this subsection**,

$$
\begin{aligned}
I(t)
&\equiv\exp\left[\int_0^t
\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}\,du\right]
=e^{\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t},\\
\frac{dI(t)}{dt}
&=\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}I(t),\\
\frac{d}{dt}\left[I(t)\Sigma_{\mathrm{H_2}}(t)\right]
&=I(t)\frac{d\Sigma_{\mathrm{H_2}}}{dt}
+\frac{dI(t)}{dt}\Sigma_{\mathrm{H_2}}\\
&=I(t)\left[
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
+\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}\Sigma_{\mathrm{H_2}}
\right].
\end{aligned}
$$

We now multiply the *entire* equation (12) by \(I(t)\). The product rule above turns its left-hand side into one derivative. On its right-hand side, the product of two exponentials adds their exponents:

$$
\begin{aligned}
\frac{d}{dt}\left[
e^{\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}
\Sigma_{\mathrm{H_2}}(t)\right]
&=e^{\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}
\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}t}\\
&=\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}
e^{(\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}-\gamma_{{\mathrm{loss}},{\mathrm{HI}}})t}.
\end{aligned}
\tag{13}
$$

Next, we integrate both sides between the known initial state and the desired time. The lower boundary on the left is \(I(0)\Sigma_{\mathrm{H_2}}(0)=\Sigma_{{\mathrm{H_2}},0}\). When the two composite rates differ, the remaining exponential integral has an elementary antiderivative:

$$
\begin{aligned}
e^{\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}
\Sigma_{\mathrm{H_2}}(t)
-\Sigma_{{\mathrm{H_2}},0}
&=\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}
\int_0^t
e^{(\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}
-\gamma_{{\mathrm{loss}},{\mathrm{HI}}})u}\,du\\
&=\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}
\frac{
e^{(\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}
-\gamma_{{\mathrm{loss}},{\mathrm{HI}}})t}-1
}{
\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}
-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}
}.
\end{aligned}
$$

Finally, we add \(\Sigma_{{\mathrm{H_2}},0}\) and multiply by \(e^{-\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}\). In the supplied-gas term, \(e^{-\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}e^{(\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}-\gamma_{{\mathrm{loss}},{\mathrm{HI}}})t}=e^{-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}t}\). We obtain

$$
\boxed{\;
\Sigma_{\mathrm{H_2}}(t)
=\Sigma_{{\mathrm{H_2}},0}e^{-\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}
+\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}
\frac{e^{-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}t}
      -e^{-\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}}
{\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}}
\;}
\tag{14}
$$

for unequal composite rates. The first term is surviving initial molecular gas. The second is transferred H I that remains after subsequent molecular loss; it is nonnegative because numerator and denominator always have the same sign. The result is easy to check at \(t=0\): both exponentials equal one, so the transfer term is zero and \(\Sigma_{\mathrm{H_2}}(0)=\Sigma_{{\mathrm{H_2}},0}\).

The user-facing quantity of interest is often \(\Sigma_{\mathrm{SFR}}\), not \(\Sigma_{\mathrm{H_2}}\) itself. Since \(\tau_{\mathrm{dep}}\) is **constant only in this special solution**, equation (2) lets us divide *each* term of equation (14) by the same depletion time. Thus we have the explicit SFR history

$$
\boxed{
\begin{aligned}
\Sigma_{\mathrm{SFR}}(t)
&=\frac{\Sigma_{\mathrm{H_2}}(t)}{\tau_{\mathrm{dep}}}\\
&=\frac{\Sigma_{{\mathrm{H_2}},0}}{\tau_{\mathrm{dep}}}
e^{-\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}\\
&\quad+
\frac{\Sigma_{{\mathrm{HI}},0}}
{\tau_{\mathrm{conv}}\tau_{\mathrm{dep}}}
\frac{
e^{-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}t}
-e^{-\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}
}{
\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}
-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}
}.
\end{aligned}
}
\tag{14a}
$$

Equation (14a) is written in **one internally consistent area/time unit system**. If the gas columns are printed in \(M_\odot\,{\mathrm{pc}}^{-2}\) and time in Gyr, then its right-hand side is initially in \(M_\odot\,{\mathrm{pc}}^{-2}\,{\mathrm{Gyr}}^{-1}\); multiply the resulting number by \(10^{-3}\) to print \(M_\odot\,{\mathrm{yr}}^{-1}\,{\mathrm{kpc}}^{-2}\). We can also divide the complete solution by \(\Sigma_{{\mathrm{SFR}},0}=\Sigma_{{\mathrm{H_2}},0}/\tau_{\mathrm{dep}}\) to see which combination of initial conditions controls the *shape*:

$$
\frac{\Sigma_{\mathrm{SFR}}(t)}
{\Sigma_{{\mathrm{SFR}},0}}
=e^{-\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}
+\frac{\Sigma_{{\mathrm{HI}},0}}
{\tau_{\mathrm{conv}}\Sigma_{{\mathrm{H_2}},0}}
\frac{
e^{-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}t}
-e^{-\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}t}
}{
\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}
-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}
}.
$$

The first coefficient is one because it carries the initial H2. The second coefficient is the *initial replenishment rate per unit initial H2*, \(1/\tau_\Phi(0)\), not another fitted stripping rate. Hence a declining H I supply changes the later SFR even if H2 was not directly stripped.

What if the two composite rates are equal? The apparent zero denominator is only an algebraic feature of writing the integral as a quotient. From the first line of the integrated equation above, its integrand becomes one when \(\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}=\gamma_{{\mathrm{loss}},{\mathrm{HI}}}\), so the integral is simply \(t\). Equivalently, the quotient in equation (14) tends to \(t e^{-\gamma_{\mathrm{loss}}t}\). Thus direct integration of equation (13) yields the continuous limit

$$
\Sigma_{\mathrm{H_2}}(t)
=e^{-\gamma_{\mathrm{loss}}t}
\left(\Sigma_{{\mathrm{H_2}},0}
+\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}\,t\right),
\qquad
\gamma_{\mathrm{loss}}
=\gamma_{{\mathrm{loss}},{\mathrm{HI}}}
=\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}.
\tag{15}
$$

The corresponding SFR is again obtained by dividing the whole expression by constant \(\tau_{\mathrm{dep}}\):

$$
\Sigma_{\mathrm{SFR}}(t)
=\frac{e^{-\gamma_{\mathrm{loss}}t}}{\tau_{\mathrm{dep}}}
\left(
\Sigma_{{\mathrm{H_2}},0}
+\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}\,t
\right)
\quad
\left(
\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}
=\gamma_{{\mathrm{loss}},{\mathrm{HI}}}
=\gamma_{\mathrm{loss}}
\right).
\tag{15a}
$$

This is the repeated-rate limit, not a physical divergence. Near that limit, evaluate the integral directly or use a stable exponential-difference function rather than subtracting two nearly equal floating-point numbers.

**Step 5: check mass conservation.** This is an important independent check of the signs. We add the two *original* equations (4), before replacing their common conversion term by a composite rate:

$$
\begin{aligned}
\frac{d}{dt}(\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H_2}})
&=\left[
-\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}}
-\gamma_{{\mathrm{strip}},{\mathrm{HI}}}\Sigma_{\mathrm{HI}}
\right]\\
&\quad+
\left[
+\frac{\Sigma_{\mathrm{HI}}}{\tau_{\mathrm{conv}}}
-(1-R+\eta)\Sigma_{\mathrm{SFR}}
-\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}\Sigma_{\mathrm{H_2}}
\right]\\
&=-(1-R+\eta)\Sigma_{\mathrm{SFR}}
-\gamma_{{\mathrm{strip}},{\mathrm{HI}}}\Sigma_{\mathrm{HI}}
-\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}\Sigma_{\mathrm{H_2}}.
\end{aligned}
\tag{16}
$$

The conversion term cancels exactly. This algebra is a meaningful double-counting check: neither molecular formation nor H I loss from that transfer can be entered again as an additional external sink. Equations (11)–(16) follow from the restricted balance (4); the cited gas-regulator and stripping papers motivate ingredients, not these particular fitted coefficients.

## 2.4 When loss rates vary through the orbit

We should not turn the special solution above into a claim that ram pressure stays constant. For time-dependent coefficients, define the cumulative exposure of each *composite* rate by \(\mathcal G_i(t)=\int_0^t\gamma_{{\mathrm{loss}},i}(u)\,du\). The H I equation then follows by integrating \(d\ln\Sigma_{\mathrm{HI}}/dt=-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}(t)\). For H2 we repeat the integrating-factor calculation, this time using \(I_{\mathrm{H_2}}(t)=e^{\mathcal G_{\mathrm{H_2}}(t)}\):

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
+\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}(t)
\Sigma_{\mathrm{H_2}}
&=\frac{\Sigma_{\mathrm{HI}}(t)}{\tau_{\mathrm{conv}}(t)},\\
\frac{d}{dt}
\left[e^{\mathcal G_{\mathrm{H_2}}(t)}
\Sigma_{\mathrm{H_2}}(t)\right]
&=e^{\mathcal G_{\mathrm{H_2}}(t)}
\frac{\Sigma_{\mathrm{HI}}(t)}{\tau_{\mathrm{conv}}(t)},\\
e^{\mathcal G_{\mathrm{H_2}}(t)}
\Sigma_{\mathrm{H_2}}(t)
-\Sigma_{{\mathrm{H_2}},0}
&=\int_0^t e^{\mathcal G_{\mathrm{H_2}}(u)}
\frac{\Sigma_{\mathrm{HI}}(u)}{\tau_{\mathrm{conv}}(u)}
\,du.
\end{aligned}
$$

Dividing the last line by \(e^{\mathcal G_{\mathrm{H_2}}(t)}\) makes the contribution of gas supplied at earlier time \(u\) explicit:

$$
\Sigma_{\mathrm{HI}}(t)
=\Sigma_{{\mathrm{HI}},0}e^{-\mathcal G_{\mathrm{HI}}(t)},
\qquad
\Sigma_{\mathrm{H_2}}(t)
=\Sigma_{{\mathrm{H_2}},0}e^{-\mathcal G_{\mathrm{H_2}}(t)}
+\int_0^t
\frac{\Sigma_{\mathrm{HI}}(u)}{\tau_{\mathrm{conv}}(u)}
e^{-[\mathcal G_{\mathrm{H_2}}(t)-\mathcal G_{\mathrm{H_2}}(u)]}\,du.
\tag{17}
$$

Here \(\gamma_{{\mathrm{loss}},{\mathrm{HI}}}(t)=\tau_{\mathrm{conv}}^{-1}(t)+\gamma_{{\mathrm{strip}},{\mathrm{HI}}}(t)\) and \(\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}(t)=(1-R+\eta)/\tau_{\mathrm{dep}}(t)+\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}(t)\). The exponential factor inside the integral is the fraction of molecular gas supplied at \(u\) that survives until \(t\) under this particular linear-loss closure. A variable rate thus produces an exponential of its **integrated exposure**, not generally a single fixed e-folding time.

Even with time-dependent molecular efficiency, we should display the SFR rather than leave it implicit. We apply equation (2) *after* solving the gas continuity equation:

$$
\boxed{
\Sigma_{\mathrm{SFR}}(t)
=\frac{1}{\tau_{\mathrm{dep}}(t)}
\left[
\Sigma_{{\mathrm{H_2}},0}e^{-\mathcal G_{\mathrm{H_2}}(t)}
+\int_0^t
\frac{\Sigma_{\mathrm{HI}}(u)}{\tau_{\mathrm{conv}}(u)}
e^{-[\mathcal G_{\mathrm{H_2}}(t)-\mathcal G_{\mathrm{H_2}}(u)]}
\,du
\right].
}
\tag{17a}
$$

Equation (17a) separates two effects: gas available at \(t\) inside the bracket, and the *current* conversion of that gas into stars through \(1/\tau_{\mathrm{dep}}(t)\). A cluster pressure pulse may change \(\gamma_{{\mathrm{strip}},i}\), but \(P_{\mathrm{ram}}(t)\) is not identified with \(\gamma_{{\mathrm{strip}},i}(t)\). [Köppen et al. (2018)](#ref-koppen) distinguish impulse and sustained-force stripping; their pulse widths cannot be assigned as universal values of either reservoir's \(\gamma\).

An initial molecular equilibrium requires

$$
\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}
=\left[
\frac{1-R+\eta}{\tau_{\mathrm{dep}}(0)}
+\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}(0^-)
\right]\Sigma_{{\mathrm{H_2}},0}.
\tag{18}
$$

The numerical example sets the pre-perturbation stripping contribution in this equation to zero. Keeping \(\Sigma_{{\mathrm{HI}},0}\) stationary *before* \(t=0\) would require an external atomic source balancing its conversion; that source is removed from equation (4) afterward. This initial condition is a model assumption, not a gas measurement in MAUVE.

We can use that initial balance to check the *first instant* after environmental removal begins. With constant \(\tau_{\mathrm{dep}}\), equation (6) evaluated at \(0^+\) becomes

$$
\begin{aligned}
\left.\frac{d\ln\Sigma_{\mathrm{SFR}}}{dt}\right|_{0^+}
&=\frac{\Sigma_{{\mathrm{HI}},0}}
{\tau_{\mathrm{conv}}\Sigma_{{\mathrm{H_2}},0}}
-\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
-\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}\\
&=\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
-\frac{1-R+\eta}{\tau_{\mathrm{dep}}}
-\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}\\
&=-\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}.
\end{aligned}
$$

The second equality uses only equation (18) with zero pre-perturbation stripping. Thus in this special equilibrium initialization, the initial **instantaneous** SFR decline directly equals the chosen direct H2 stripping rate. At later times the replenishment term evolves, so the measured SFR decline is no longer a direct measurement of that coefficient. This is a restricted diagnostic, not a general inversion of an observed stage difference.

# 3. A local positive response within the same gas model

Compression or changing cloud dynamics can increase H2 supply, increase the molecular fraction, or reduce \(\tau_{\mathrm{dep}}\). Those mechanisms have different observational consequences. [Brown et al. (2023)](#ref-brown) report an early Virgo-outskirts enhancement associated with molecular content, without the same increase in molecular efficiency; [Lizee et al. (2021)](#ref-lizee) discuss a region of NGC4654 with inferred efficiency enhancement under their gas model. The regulator identity

$$
\frac{d\ln\Sigma_{\mathrm{SFR}}}{dt}
=\frac{d\ln\Sigma_{\mathrm{H_2}}}{dt}
 -\frac{d\ln\tau_{\mathrm{dep}}}{dt}
\tag{19}
$$

shows directly why higher SFR alone cannot distinguish the two. Equation (9) states exactly how an increase in \(\Sigma_\Phi/\Sigma_{\mathrm{H_2}}\), a reduction in direct removal, or a decrease in \(\tau_{\mathrm{dep}}\) could reverse the local sign of \(d\Sigma_{\mathrm{SFR}}/dt\).

For one *smooth, illustrative* efficiency event, prescribe

$$
\tau_{\mathrm{dep}}(t)
=\frac{\tau_{\mathrm{dep,0}}}{1+g(t)},
\qquad
g(t)=A_\epsilon
\left(\frac{t}{t_{\mathrm{pk}}}\right)^m
\exp\!\left[m\!\left(1-\frac{t}{t_{\mathrm{pk}}}\right)\right]
\quad(t\geq0).
\tag{20}
$$

\(A_\epsilon\geq0\) is the **fractional peak increase in molecular efficiency**, \(t_{\mathrm{pk}}>0\) is its peak time after the chosen onset, and the dimensionless shape index \(m>1\) sets how concentrated the event is while giving a finite zero initial slope. This function satisfies \(g(0)=0\), \(g(t_{\mathrm{pk}})=A_\epsilon\), and later decays. It is a deliberately labelled phenomenological forcing, *not* a measured conversion of \(P_{\mathrm{ram}}\) into \(\tau_{\mathrm{dep}}\), and it is not required by the constant-coefficient gas solution.

Differentiation, which is needed for the exact sign test, gives

$$
\frac{dg}{dt}
=m g(t)\left(\frac{1}{t}-\frac{1}{t_{\mathrm{pk}}}\right),
\qquad
\frac{d\ln\tau_{\mathrm{dep}}}{dt}
=-\frac{1}{1+g(t)}\frac{dg}{dt}.
\tag{21}
$$

Equation (21) is evaluated for \(t>0\); its right-hand derivative has a continuous zero limit at \(t=0\) because we require \(m>1\). On the rising side, \(dg/dt>0\), so \(d\ln\tau_{\mathrm{dep}}/dt<0\). To decide whether this is enough to raise the **actual SFR**, we must still solve for the molecular gas that remains. For the pulse branch we keep \(\tau_{\mathrm{conv}}\) and both stripping coefficients constant but allow equation (20)'s \(\tau_{\mathrm{dep}}(t)\) to vary. The molecular composite exposure is then

$$
\begin{aligned}
\mathcal G_{\mathrm{H_2}}^{\mathrm{pulse}}(t)
&=\int_0^t\left[
\frac{1-R+\eta}{\tau_{\mathrm{dep}}(u)}
+\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}
\right]du\\
&=\left[
\frac{1-R+\eta}{\tau_{\mathrm{dep,0}}}
+\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}
\right]t
+\frac{1-R+\eta}{\tau_{\mathrm{dep,0}}}
\int_0^t g(u)\,du.
\end{aligned}
$$

The first part is the no-pulse molecular exposure; the integral is the *additional molecular consumption* associated with the temporary efficiency rise. Substitution into the already derived equation (17), followed by equation (2), gives the pulse branch without inventing a second gas model:

$$
\begin{aligned}
\Sigma_{\mathrm{H_2}}^{\mathrm{pulse}}(t)
&=\Sigma_{{\mathrm{H_2}},0}
e^{-\mathcal G_{\mathrm{H_2}}^{\mathrm{pulse}}(t)}
+\frac{\Sigma_{{\mathrm{HI}},0}}{\tau_{\mathrm{conv}}}
\int_0^t
e^{-\gamma_{{\mathrm{loss}},{\mathrm{HI}}}u}
e^{-[\mathcal G_{\mathrm{H_2}}^{\mathrm{pulse}}(t)
-\mathcal G_{\mathrm{H_2}}^{\mathrm{pulse}}(u)]}\,du,\\
\Sigma_{\mathrm{SFR}}^{\mathrm{pulse}}(t)
&=\frac{1+g(t)}{\tau_{\mathrm{dep,0}}}
\Sigma_{\mathrm{H_2}}^{\mathrm{pulse}}(t).
\end{aligned}
\tag{21a}
$$

In the numerical calculation we integrate the pulse ODE directly; equation (21a) states the analytical solution as a one-dimensional causal integral and makes its assumptions visible. We can now differentiate its *last* line and use the molecular balance once more:

$$
\begin{aligned}
\frac{d\ln\Sigma_{\mathrm{SFR}}^{\mathrm{pulse}}}{dt}
&=\frac{1}{\Sigma_{\mathrm{H_2}}^{\mathrm{pulse}}}
\frac{d\Sigma_{\mathrm{H_2}}^{\mathrm{pulse}}}{dt}
+\frac{1}{1+g}\frac{dg}{dt}\\
&=\frac{\Sigma_{\mathrm{HI}}}
{\tau_{\mathrm{conv}}\Sigma_{\mathrm{H_2}}^{\mathrm{pulse}}}
-\frac{1-R+\eta}{\tau_{\mathrm{dep}}(t)}
-\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}
+\frac{1}{1+g}\frac{dg}{dt}.
\end{aligned}
$$

Thus a *local SFR maximum* satisfies

$$
\frac{1}{1+g}\frac{dg}{dt}
=\frac{1-R+\eta}{\tau_{\mathrm{dep}}(t)}
+\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}
-\frac{\Sigma_{\mathrm{HI}}}
{\tau_{\mathrm{conv}}\Sigma_{\mathrm{H_2}}^{\mathrm{pulse}}}.
$$

At the maximum of the imposed efficiency profile, \(dg/dt=0\), but the right-hand side need not be zero: the gas has already changed. Therefore the SFR maximum need not occur at \(t_{\mathrm{pk}}\). At the maximum of the *SFR*, \(\Gamma_{\mathrm{SFR}}=0\); immediately afterward \(\Gamma_{\mathrm{SFR}}>0\) can coexist with \(\mathcal E<0\), because SFR may already be falling while still above its initial value. The worked calculation tests both statements.

No optical residual alone tells us whether the NGC4654 facing side acquired more H2, changed \(\tau_{\mathrm{dep}}\), changed its ionizing-photon capture, or reflects a different mixture of galaxies/regions. The 0.035-dex environmental offset is consequently used only as a scale for a *possible* positive branch, not as a fitted pulse amplitude.

# 4. From the gas solution to H-alpha

## 4.1 Why a stellar-age response appears, and when it is negligible

The gas model supplies the **instantaneous true** \(\Sigma_{\mathrm{SFR}}(t)\). H-alpha responds to ionizing photons from young stellar populations that formed at several preceding ages; it is not mathematically identical to the instantaneous gas-based SFR during an abrupt change. Define \(q_{\mathrm{H,y}}(a)\) as the hydrogen-ionizing photon production rate, in photons s\(^{-1}\,M_\odot^{-1}\), per unit formed stellar mass of age \(a\). A cohort formed in the age interval \(da\) contributed \(\Sigma_{\mathrm{SFR}}(t-a)da\) in mass per area, so adding the cohorts gives

$$
Q_{\mathrm{H,y}}(t)
=\int_0^\infty q_{\mathrm{H,y}}(a)
\Sigma_{\mathrm{SFR}}(t-a)\,da,
\qquad
K_\alpha(a)
=\frac{q_{\mathrm{H,y}}(a)}
{\int_0^\infty q_{\mathrm{H,y}}(u)\,du}.
\tag{22}
$$

\(Q_{\mathrm{H,y}}\) is the young-star ionizing-photon production rate per area, not an H-alpha luminosity; the age and SFR time units must match. \(K_\alpha\) is a nonnegative, unit-normalized response with units inverse time. Equation (22) is stellar-cohort summation; stellar synthesis determines \(q_{\mathrm{H,y}}\). [Kennicutt & Evans (2012), Table 1](#ref-ke12) give a mean H-alpha contributing age near 3 Myr and 90% of their adopted signal from ages below approximately 10 Myr, under their population-model assumptions. That age information motivates a short response; it is not evidence for a specific exponential kernel.

Define the recent-SFR equivalent measured by that response:

$$
\overline{\Sigma}_{\mathrm{SFR,\alpha}}(t)
=\int_0^\infty K_\alpha(a)
\Sigma_{\mathrm{SFR}}(t-a)\,da.
\tag{23}
$$

For a slowly varying SFR, the adopted MAUVE H-alpha conversion is approximately \(\Sigma_{\mathrm{SFR}}=C_\alpha\mathcal L_{\alpha}\) when the calibration's young photons are captured in the modelled hydrogen gas. The live pipeline uses \(C_\alpha=4.983582\times10^{-42}\ M_\odot\,{\mathrm{yr}}^{-1}/({\mathrm{erg\,s}}^{-1})\): it starts from its Kennicutt–Evans Kroupa coefficient and applies the stated Chabrier conversion in **SFR+Z.py**. This differs slightly from the directly tabulated Kroupa coefficient in [Kennicutt & Evans (2012)](#ref-ke12); the report uses the *pipeline* value for comparison, with consistent area factors.

To see the mathematics without claiming a population-synthesis prediction, choose \(K_\alpha(a)=\tau_{\mathrm{ion}}^{-1}e^{-a/\tau_{\mathrm{ion}}}\) for \(a\geq0\). We first change the integration variable in equation (23) from stellar age \(a\) to formation time \(u=t-a\). The upper age limit becomes the lower formation-time limit, so, assuming the prehistory is defined,

$$
\begin{aligned}
\overline{\Sigma}_{\mathrm{SFR,\alpha}}(t)
&=\frac{1}{\tau_{\mathrm{ion}}}
\int_{-\infty}^{t}e^{-(t-u)/\tau_{\mathrm{ion}}}
\Sigma_{\mathrm{SFR}}(u)\,du\\
&=\frac{e^{-t/\tau_{\mathrm{ion}}}}{\tau_{\mathrm{ion}}}
\int_{-\infty}^{t}e^{u/\tau_{\mathrm{ion}}}
\Sigma_{\mathrm{SFR}}(u)\,du.
\end{aligned}
$$

We now differentiate **both** the prefactor and the upper limit. The derivative of the integral is its integrand evaluated at \(u=t\), by the fundamental theorem of calculus:

$$
\begin{aligned}
\frac{d\overline{\Sigma}_{\mathrm{SFR,\alpha}}}{dt}
&=-\frac{1}{\tau_{\mathrm{ion}}}
\overline{\Sigma}_{\mathrm{SFR,\alpha}}(t)
+\frac{e^{-t/\tau_{\mathrm{ion}}}}{\tau_{\mathrm{ion}}}
e^{t/\tau_{\mathrm{ion}}}\Sigma_{\mathrm{SFR}}(t)\\
&=\frac{\Sigma_{\mathrm{SFR}}(t)
-\overline{\Sigma}_{\mathrm{SFR,\alpha}}(t)}{\tau_{\mathrm{ion}}}.
\end{aligned}
$$

Multiplication by \(\tau_{\mathrm{ion}}\) and rearrangement yield

$$
\tau_{\mathrm{ion}}\frac{d\overline{\Sigma}_{\mathrm{SFR,\alpha}}}{dt}
+\overline{\Sigma}_{\mathrm{SFR,\alpha}}
=\Sigma_{\mathrm{SFR}}(t),
\qquad
\overline{\Sigma}_{\mathrm{SFR,\alpha}}(0)
=\Sigma_{\mathrm{SFR}}(0)
\quad\text{if the earlier SFR was steady}.
\tag{24}
$$

The example takes \(\tau_{\mathrm{ion}}=3\) Myr. Its exponential kernel has mean 3 Myr and 90th percentile 6.9 Myr; it is deliberately approximate and does not reproduce the full stellar-population age distribution. To determine when we may omit this extra response equation, Taylor-expand the source about the observation time, \(\Sigma_{\mathrm{SFR}}(t-a)=\Sigma_{\mathrm{SFR}}(t)-a\,d\Sigma_{\mathrm{SFR}}/dt+(a^2/2)\,d^2\Sigma_{\mathrm{SFR}}/dt^2+\cdots\). The chosen normalized kernel has \(\int_0^\infty aK_\alpha(a)da=\tau_{\mathrm{ion}}\) and \(\int_0^\infty a^2K_\alpha(a)da=2\tau_{\mathrm{ion}}^2\). Inserting the series term by term into equation (23) gives

$$
\overline{\Sigma}_{\mathrm{SFR,\alpha}}(t)
=\Sigma_{\mathrm{SFR}}(t)
-\tau_{\mathrm{ion}}\frac{d\Sigma_{\mathrm{SFR}}}{dt}
+\tau_{\mathrm{ion}}^{\,2}
\frac{d^2\Sigma_{\mathrm{SFR}}}{dt^2}
+O\!\left(\tau_{\mathrm{ion}}^{\,3}
\frac{d^3\Sigma_{\mathrm{SFR}}}{dt^3}\right).
\tag{25}
$$

The kernel is therefore **not necessary for a useful slow-evolution approximation**: setting \(\overline{\Sigma}_{\mathrm{SFR,\alpha}}\simeq\Sigma_{\mathrm{SFR}}\) is justified when the change time is long compared with the few-Myr response. It matters near an abrupt burst or shutdown, and at approximately 100-pc resolution where finite stellar sampling and cloud lifetimes can matter; see the spatial-scale caution of [Kruijssen & Longmore (2014)](#ref-kl14). The numerical example explicitly computes both versions and finds a 0.76% difference at its illustrative 1-Gyr endpoint. That small result applies to *its* smooth trajectory only.

## 4.2 Case B, absorption, and a photon budget with no double counting

For recombination-dominated hydrogen at specified electron temperature and density, define \(n_e,n_p\) as electron and proton number densities, \(\alpha_B\) as the total Case-B recombination coefficient, and \(\alpha_\alpha^{\mathrm{eff}}\) as the effective coefficient producing H-alpha photons. Let \(A_{\mathrm{kpc}}\equiv9.521\times10^{42}\ {\mathrm{cm^2\,kpc^{-2}}}\) be the **surface-unit conversion factor**, the number of square centimetres in one square kiloparsec. Multiplying a flux per cm² by \(A_{\mathrm{kpc}}\) gives the same flux expressed per kpc². The volume H-alpha emissivity and absorbed-photon balance in the latter units are

$$
j_\alpha
=h\nu_\alpha\alpha_\alpha^{\mathrm{eff}}n_en_p,
\qquad
Q_{\mathrm{H,abs}}
=A_{\mathrm{kpc}}\alpha_B\int n_en_p\,dz,
\qquad
p_\alpha\equiv
\frac{\alpha_\alpha^{\mathrm{eff}}}{\alpha_B}.
\tag{26}
$$

Here \(h\nu_\alpha\) is the energy per H-alpha photon, \(Q_{\mathrm{H,abs}}\) is the rate of ionizing photons absorbed by the gas per projected kpc², and \(p_\alpha\) is the mean H-alpha-photon yield per absorbed ionization under steady Case-B balance. The symbol \(z\) is distance through the emitting gas along the line of sight, in cm inside the integral. We integrate the volume emissivity along that path and then convert from per cm² to per kpc² with \(A_{\mathrm{kpc}}\). Substituting the absorption balance in equation (26), **including the same line-of-sight integral**, gives the two equivalent forms:

$$
\begin{aligned}
\mathcal L_\alpha
&=A_{\mathrm{kpc}}\int j_\alpha\,dz\\
&=A_{\mathrm{kpc}}h\nu_\alpha\alpha_\alpha^{\mathrm{eff}}
\int n_en_p\,dz\\
&=h\nu_\alpha
\frac{\alpha_\alpha^{\mathrm{eff}}}{\alpha_B}
\left(A_{\mathrm{kpc}}\alpha_B\int n_en_p\,dz\right)\\
&=h\nu_\alpha p_\alpha Q_{\mathrm{H,abs}}.
\end{aligned}
\tag{27}
$$

The two forms are **photon-budget** and **emission-measure** forms of the same recombination relation. Their atomic coefficients and validity limits come from [Hummer & Storey (1987)](#ref-hs87) and [Storey & Hummer (1995)](#ref-sh95). This is direct recombination bookkeeping, rather than a distinct empirical “Draine equation.” The expressions require ionization balance and an appropriate Case-B regime. Shock excitation can supply or modify lines through further processes not represented by one stellar photon rate.

Now distinguish **where young-star photons go**. Let \(f_{\mathrm{HII}}\) be the fraction absorbed by compact H II gas; \(f_{\mathrm{leak}}\) the fraction of the *same young photon budget* that leaves those regions and is absorbed by more diffuse, non-H II gas; and \(f_{\mathrm{unabs}}\) the fraction not absorbed by the modelled hydrogen gas because of escape, dust absorption, or an unmodelled destination. To close this **local** photon budget, the numerical experiment assumes that photons produced inside the chosen aperture either are absorbed there or enter \(f_{\mathrm{unabs}}\); it does not add photons produced in neighbouring apertures. A resolved diffuse-gas model would need to transport leaked photons between positions, as discussed by [Belfiore et al. (2022)](#ref-belfiore22). Under this declared local simplification, the fractions are dimensionless, disjoint, nonnegative, and satisfy

$$
f_{\mathrm{HII}}+f_{\mathrm{leak}}+f_{\mathrm{unabs}}=1,
\qquad
Q_{\mathrm{H,y,abs}}^{\mathrm{HII}}=f_{\mathrm{HII}}Q_{\mathrm{H,y}},
\qquad
Q_{\mathrm{H,y,abs}}^{\mathrm{non-HII}}=f_{\mathrm{leak}}Q_{\mathrm{H,y}}.
\tag{28}
$$

The SFR calibration is the statement that, **for its assumed stellar population and complete absorption of the modelled young-star photons**, the fully absorbed H-alpha luminosity per projected area would be \(\overline{\Sigma}_{\mathrm{SFR,\alpha}}/C_\alpha\). Equation (27) therefore identifies the corresponding young-star photon budget as \(Q_{\mathrm{H,y}}=\overline{\Sigma}_{\mathrm{SFR,\alpha}}/(C_\alpha h\nu_\alpha p_\alpha)\) under that calibration. We multiply this *one* production rate by each disjoint destination in equation (28), then apply the Case-B yield:

$$
\begin{aligned}
\mathcal L_\alpha^{\mathrm{HII}}
&=h\nu_\alpha p_\alpha
Q_{\mathrm{H,y,abs}}^{\mathrm{HII}}
=h\nu_\alpha p_\alpha f_{\mathrm{HII}}Q_{\mathrm{H,y}},\\
\mathcal L_\alpha^{\mathrm{non-HII,y}}
&=h\nu_\alpha p_\alpha
Q_{\mathrm{H,y,abs}}^{\mathrm{non-HII}}
=h\nu_\alpha p_\alpha f_{\mathrm{leak}}Q_{\mathrm{H,y}}.
\end{aligned}
$$

Inserting the calibrated expression for \(Q_{\mathrm{H,y}}\) cancels \(h\nu_\alpha p_\alpha\) in both lines and yields

$$
\mathcal L_\alpha^{\mathrm{HII}}
=f_{\mathrm{HII}}\,
\frac{\overline{\Sigma}_{\mathrm{SFR,\alpha}}}{C_\alpha},
\qquad
\mathcal L_\alpha^{\mathrm{non-HII,y}}
=f_{\mathrm{leak}}\,
\frac{\overline{\Sigma}_{\mathrm{SFR,\alpha}}}{C_\alpha}.
\tag{29}
$$

Equation (29) applies the same Case-B calibration to separate destinations of **one** young-photon budget; it does not create another copy of young-star radiation. [Belfiore et al. (2022)](#ref-belfiore22) motivate the leaked-photon and evolved-star channels but do not provide the numerical fractions below. If one chooses the special simplification \(f_{\mathrm{HII}}=1\), then equation (28) forces \(f_{\mathrm{leak}}=f_{\mathrm{unabs}}=0\), and equation (29)'s H II expression reduces to \(\overline{\Sigma}_{\mathrm{SFR,\alpha}}/C_\alpha\). One may drop the printed fraction in that *special case*, as suggested in the discussion. It cannot simultaneously power a young-leakage non-H II component. The data-scaled example intentionally **does not** set \(f_{\mathrm{HII}}=1\).

Other non-H II excitation can be added without reusing those young photons. For hot evolved low-mass stars (HOLMES), define \(q_{\mathrm{H,old}}\) as ionizing photons s\(^{-1}\,M_\odot^{-1}\), \(f_{\mathrm{abs,old}}\leq1\) as the gas-absorbed fraction of that **independent** old-stellar budget, and \(\Sigma_*\) as stellar mass per area. In a recombination-powered old-star subcomponent,

$$
\mathcal L_\alpha^{\mathrm{old}}
=h\nu_\alpha p_\alpha
f_{\mathrm{abs,old}}q_{\mathrm{H,old}}\Sigma_*,
\qquad
\mathcal L_\alpha^{\mathrm{non-HII}}
=\mathcal L_\alpha^{\mathrm{non-HII,y}}
+\mathcal L_\alpha^{\mathrm{old}}
+\mathcal L_\alpha^{\mathrm{mech}}
+\mathcal L_\alpha^{\mathrm{AGN}}.
\tag{30}
$$

The last two labels stand for mechanically and AGN-powered contributions and are **set to zero in the numerical example**. Equation (30)'s old-star factor is a photon-counting consequence of equation (27), with the source population discussed by [Belfiore et al. (2022)](#ref-belfiore22). Crucially, \(f_{\mathrm{abs,old}}\) should **not** be dropped: gas removal can make old-star line emission fade even while \(q_{\mathrm{H,old}}\Sigma_*\) stays roughly fixed. [Belfiore et al. (2017)](#ref-belfiore16) explicitly note the gas-absorption condition.

As a scale check, the [Belfiore et al. (2022)](#ref-belfiore22) HOLMES model has \(q_{\mathrm{H,old}}\) of order \(7\times10^{40}\ {\mathrm{photons\,s^{-1}}}\,M_\odot^{-1}\). With \(\Sigma_*=10^{8.625}\ M_\odot\,{\mathrm{kpc}}^{-2}\), \(h\nu_\alpha=3.027\times10^{-12}\) erg, \(p_\alpha\simeq0.45\), and the deliberately maximal \(f_{\mathrm{abs,old}}=1\), equation (30) gives

$$
\mathcal L_{\alpha,\mathrm{old}}^{\mathrm{max}}
\simeq4.02\times10^{37}
\ {\mathrm{erg\,s^{-1}\,kpc^{-2}}},
\tag{31}
$$

approximately **3.6%** of Table 1's post-peak NSF mean H-alpha at the same mass bin. This is a *conditional ceiling under the adopted old-population ionizing yield*, not a theorem about every evolved population or an upper bound on young-leakage, shocks, or AGN. It makes an old-stars-only explanation of that **mean, common-detection-support** luminosity quantitatively strained under those assumptions. Differences in selected population, source model, gas absorption, and attenuation treatment require direct checking before a source fraction is inferred.

## 4.3 When can non-H II emission fade more slowly?

The report's minimal numerical example omits old, shock, and AGN light and uses only the disjoint young-photon channels of equations (28)–(29). It couples their absorption fractions to the gas solution through the following **declared, testable closure**, not an established empirical law:

$$
\begin{aligned}
f_{\mathrm{HII}}(t)
&=f_{{\mathrm{HII}},0}
\frac{\Sigma_{\mathrm{H_2}}(t)}{\Sigma_{{\mathrm{H_2}},0}},\\
f_{\mathrm{unabs}}(t)
&=f_{{\mathrm{unabs}},0}
+(f_{{\mathrm{unabs}},\infty}-f_{{\mathrm{unabs}},0})
\left[1-\frac{\Sigma_{\mathrm{HI}}(t)}{\Sigma_{{\mathrm{HI}},0}}\right],\\
f_{\mathrm{leak}}(t)
&=1-f_{\mathrm{HII}}(t)-f_{\mathrm{unabs}}(t).
\end{aligned}
\tag{32}
$$

The first line assumes the compact H II absorption fraction falls in proportion to the molecular reservoir in this particular transformed patch. The second allows a bounded increase in photons not absorbed by the modelled gas as the atomic reservoir disappears. The remaining photons are *assigned* to diffuse gas. Equation (32) is admissible only where \(0\leq f_{\mathrm{HII}},f_{\mathrm{unabs}}\leq1\) and \(f_{\mathrm{HII}}+f_{\mathrm{unabs}}\leq1\), so that the derived \(f_{\mathrm{leak}}\) is nonnegative; the numerical trajectory is checked against these bounds. These assumptions also require that diffuse gas and recombination capacity actually persist; equations (26)–(27) provide the independent check that the photon fractions alone cannot guarantee. The limiting fraction \(f_{{\mathrm{unabs}},\infty}\) is a chosen ceiling, not a measured escape fraction.

Logarithmically differentiating equation (29) reveals the precise slower-fading criterion when both fractions are positive. For example, \(d\ln(f_{\mathrm{HII}}\overline{\Sigma}_{\mathrm{SFR,\alpha}}/C_\alpha)/dt=d\ln f_{\mathrm{HII}}/dt+d\ln\overline{\Sigma}_{\mathrm{SFR,\alpha}}/dt\), because \(C_\alpha\) is constant in this experiment. Repeating the same step for the leaked channel gives

$$
\begin{aligned}
\frac{d\ln\mathcal L_\alpha^{\mathrm{HII}}}{dt}
&=\frac{d\ln\overline{\Sigma}_{\mathrm{SFR,\alpha}}}{dt}
+\frac{d\ln f_{\mathrm{HII}}}{dt},\\
\frac{d\ln\mathcal L_\alpha^{\mathrm{non-HII,y}}}{dt}
&=\frac{d\ln\overline{\Sigma}_{\mathrm{SFR,\alpha}}}{dt}
+\frac{d\ln f_{\mathrm{leak}}}{dt}.
\end{aligned}
\tag{33}
$$

If \(f_{\mathrm{HII}}\) falls while \(f_{\mathrm{leak}}\) grows, compact H II emission fades more rapidly. Non-H II emission itself still **fades** only when \(-d\ln\overline{\Sigma}_{\mathrm{SFR,\alpha}}/dt>d\ln f_{\mathrm{leak}}/dt\). During an early redistribution it can instead brighten temporarily even as the total H-alpha falls; in the computed example it peaks and then ends below its initial level. This is a consequence of the specified source partition, *not* a universal multi-hundred-Myr “non-H II lifetime.” An old-star-powered component follows its own \(f_{\mathrm{abs,old}}\), and rapid loss of its gas can make it fade first. Thus the often-assumed inequality \(\tau_{\mathrm{non-HII}}>\tau_{\mathrm{HII}}\) is not imposed as a general law.

# 5. Individual forbidden lines and BPT ratios

## 5.1 Add physical line luminosities before taking ratios

For each physical emission component \(j\in\{{\mathrm{HII}},{\mathrm{non-HII}}\}\), let \(\mathcal L_\ell^j\) be the luminosity surface density in line \(\ell\), and let \(\mathcal R_{\ell/B}^{j}=\mathcal L_\ell^{j}/\mathcal L_B^{j}\) be its intrinsic **linear** ratio to the appropriate Balmer line \(B\). [N II] \(\lambda6583\) and the [S II] \(\lambda\lambda6716,6730\) **sum** use H-alpha; [O III] \(\lambda5007\) uses H-beta. Adding contributions first gives

$$
\mathcal L_\ell
=\mathcal L_\ell^{\mathrm{HII}}
+\mathcal L_\ell^{\mathrm{non-HII}}
=\mathcal R_{\ell/B}^{\mathrm{HII}}\mathcal L_B^{\mathrm{HII}}
+\mathcal R_{\ell/B}^{\mathrm{non-HII}}\mathcal L_B^{\mathrm{non-HII}}.
\tag{34}
$$

This is luminosity conservation for the chosen decomposition, not a relation imported from a photoionization grid. [Belfiore et al. (2022)](#ref-belfiore22) and [Zhang et al. (2017)](#ref-zhang) motivate why physically distinct emission can occupy the same spatial element and alter measured BPT coordinates. The components' line ratios are separate physical inputs; the observed SF/NSF **class means** in Table 1 do not directly measure them.

For completeness, define each intrinsic Balmer decrement \(B_j=\mathcal L_\alpha^j/\mathcal L_\beta^j\). The numerical model chooses

$$
B_{\mathrm{HII}}=B_{\mathrm{non-HII}}=2.86,
\qquad
\mathcal L_\beta^j
=\frac{\mathcal L_\alpha^j}{2.86}.
\tag{35}
$$

This is the familiar approximate **intrinsic** Case-B H-alpha/H-beta value near \(T_e=10^4\) K and low nebular density, based on the Case-B atomic calculations of [Hummer & Storey (1987)](#ref-hs87) and [Storey & Hummer (1995)](#ref-sh95). It simplifies this controlled experiment; it is not an assertion that observed or pipeline-corrected ratios all equal 2.86. Temperature, density, optical depth, collisions, and differential dust attenuation can change component decrements. The live dust routine clips negative inferred reddening but does not force every corrected Balmer ratio to precisely 2.86.

Let \(w_\alpha=\mathcal L_\alpha^{\mathrm{non-HII}}/\mathcal L_\alpha\) and \(w_\beta=\mathcal L_\beta^{\mathrm{non-HII}}/\mathcal L_\beta\). We should work out both Balmer denominators explicitly. For [N II] or [S II], divide the corresponding sum in equation (34) by \(\mathcal L_\alpha=\mathcal L_\alpha^{\mathrm{HII}}+\mathcal L_\alpha^{\mathrm{non-HII}}\):

$$
\begin{aligned}
\frac{\mathcal L_{\mathrm{[N\,II]}}}{\mathcal L_\alpha}
&=\frac{
\mathcal R_{\mathrm{[N\,II]}/\alpha}^{\mathrm{HII}}
\mathcal L_\alpha^{\mathrm{HII}}
+\mathcal R_{\mathrm{[N\,II]}/\alpha}^{\mathrm{non-HII}}
\mathcal L_\alpha^{\mathrm{non-HII}}
}{\mathcal L_\alpha}\\
&=\mathcal R_{\mathrm{[N\,II]}/\alpha}^{\mathrm{HII}}
\frac{\mathcal L_\alpha^{\mathrm{HII}}}{\mathcal L_\alpha}
+\mathcal R_{\mathrm{[N\,II]}/\alpha}^{\mathrm{non-HII}}
\frac{\mathcal L_\alpha^{\mathrm{non-HII}}}{\mathcal L_\alpha}\\
&=(1-w_\alpha)
\mathcal R_{\mathrm{[N\,II]}/\alpha}^{\mathrm{HII}}
+w_\alpha
\mathcal R_{\mathrm{[N\,II]}/\alpha}^{\mathrm{non-HII}}.
\end{aligned}
$$

The identical derivation applies to the [S II] doublet **sum**. For [O III], replace every \(\alpha\) above by \(\beta\); its weight is \(w_\beta\), which need not equal \(w_\alpha\). To relate the weights, equation (35) before the equal-decrement choice says \(\mathcal L_\beta^j=\mathcal L_\alpha^j/B_j\). Substituting \(\mathcal L_\alpha^{\mathrm{HII}}=(1-w_\alpha)\mathcal L_\alpha\) and \(\mathcal L_\alpha^{\mathrm{non-HII}}=w_\alpha\mathcal L_\alpha\) gives

$$
\begin{aligned}
w_\beta
&=\frac{\mathcal L_\alpha^{\mathrm{non-HII}}/B_{\mathrm{non-HII}}}
{\mathcal L_\alpha^{\mathrm{HII}}/B_{\mathrm{HII}}
+\mathcal L_\alpha^{\mathrm{non-HII}}/B_{\mathrm{non-HII}}}\\
&=\frac{w_\alpha\mathcal L_\alpha/B_{\mathrm{non-HII}}}
{(1-w_\alpha)\mathcal L_\alpha/B_{\mathrm{HII}}
+w_\alpha\mathcal L_\alpha/B_{\mathrm{non-HII}}}.
\end{aligned}
$$

The common total \(\mathcal L_\alpha\) cancels. We can now write the general result for either Balmer denominator:

$$
\begin{aligned}
\mathcal R_{\ell/\alpha}
&=(1-w_\alpha)\mathcal R_{\ell/\alpha}^{\mathrm{HII}}
+w_\alpha\mathcal R_{\ell/\alpha}^{\mathrm{non-HII}},\\
\mathcal R_{\ell/\beta}
&=(1-w_\beta)\mathcal R_{\ell/\beta}^{\mathrm{HII}}
+w_\beta\mathcal R_{\ell/\beta}^{\mathrm{non-HII}},\\
w_\beta
&=\frac{w_\alpha/B_{\mathrm{non-HII}}}
{(1-w_\alpha)/B_{\mathrm{HII}}+w_\alpha/B_{\mathrm{non-HII}}}.
\end{aligned}
\tag{36}
$$

For a nontrivial mixture with \(0<w_\alpha<1\), equality \(w_\beta=w_\alpha\) requires equal component decrements \(B_{\mathrm{HII}}=B_{\mathrm{non-HII}}\), as chosen in equation (35). A pure single-component endpoint has equal weights trivially. A ratio is mixed **linearly in luminosity**; the logarithm is taken afterward for a BPT plot. Mixing two already logarithmic coordinates would be mathematically wrong.

## 5.2 Derive differential fading instead of assigning each line a lifetime

Differentiate the definition of \(w_\alpha\) by the quotient rule. We retain the physical component names throughout. Because \(\mathcal L_\alpha=\mathcal L_\alpha^{\mathrm{HII}}+\mathcal L_\alpha^{\mathrm{non-HII}}\), its derivative is the sum of the two component derivatives:

$$
\begin{aligned}
\frac{dw_\alpha}{dt}
&=\frac{\mathcal L_\alpha\,
d\mathcal L_\alpha^{\mathrm{non-HII}}/dt
-\mathcal L_\alpha^{\mathrm{non-HII}}\,
d\mathcal L_\alpha/dt}{\mathcal L_\alpha^2}\\
&=\frac{
\mathcal L_\alpha^{\mathrm{HII}}\,
d\mathcal L_\alpha^{\mathrm{non-HII}}/dt
-\mathcal L_\alpha^{\mathrm{non-HII}}\,
d\mathcal L_\alpha^{\mathrm{HII}}/dt}
{\mathcal L_\alpha^2}\\
&=\frac{
\mathcal L_\alpha^{\mathrm{HII}}
\mathcal L_\alpha^{\mathrm{non-HII}}}
{\mathcal L_\alpha^2}
\left[
\frac{1}{\mathcal L_\alpha^{\mathrm{non-HII}}}
\frac{d\mathcal L_\alpha^{\mathrm{non-HII}}}{dt}
-\frac{1}{\mathcal L_\alpha^{\mathrm{HII}}}
\frac{d\mathcal L_\alpha^{\mathrm{HII}}}{dt}
\right].
\end{aligned}
$$

The prefactor in the last line is \(w_\alpha(1-w_\alpha)\), and each reciprocal-luminosity derivative is a logarithmic derivative. Thus, for positive component luminosities,

$$
\frac{dw_\alpha}{dt}
=w_\alpha(1-w_\alpha)
\left[
\frac{d\ln\mathcal L_\alpha^{\mathrm{non-HII}}}{dt}
-\frac{d\ln\mathcal L_\alpha^{\mathrm{HII}}}{dt}
\right].
\tag{37}
$$

Thus \(w_\alpha\) grows if the **fractional** non-H II H-alpha luminosity declines more slowly than the H II luminosity, even while both endpoints are fainter. Substituting equation (33) shows that, in the two young-photon channels with the same \(\overline{\Sigma}_{\mathrm{SFR,\alpha}}\), their shared stellar history cancels:

$$
\frac{dw_\alpha}{dt}
=w_\alpha(1-w_\alpha)
\left[
\frac{d\ln f_{\mathrm{leak}}}{dt}
-\frac{d\ln f_{\mathrm{HII}}}{dt}
\right].
\tag{38}
$$

This gives a concrete physical meaning to “differential fading” in the example: photon destinations change as compact gas is lost, while one common young-star source diminishes. If old stars or shocks contribute, their separate histories enter equation (37) and equation (38) no longer holds by itself.

Differentiate the first line of equation (36), allowing each component spectrum itself to evolve. The product rule acts on both \(1-w_\alpha\) and \(w_\alpha\):

$$
\begin{aligned}
\frac{d\mathcal R_{\ell/\alpha}}{dt}
&=-\frac{dw_\alpha}{dt}\mathcal R_{\ell/\alpha}^{\mathrm{HII}}
+(1-w_\alpha)\frac{d\mathcal R_{\ell/\alpha}^{\mathrm{HII}}}{dt}\\
&\quad+\frac{dw_\alpha}{dt}\mathcal R_{\ell/\alpha}^{\mathrm{non-HII}}
+w_\alpha\frac{d\mathcal R_{\ell/\alpha}^{\mathrm{non-HII}}}{dt}.
\end{aligned}
$$

Collecting the two terms proportional to \(dw_\alpha/dt\) gives

$$
\frac{d\mathcal R_{\ell/\alpha}}{dt}
=\left(
\mathcal R_{\ell/\alpha}^{\mathrm{non-HII}}
-\mathcal R_{\ell/\alpha}^{\mathrm{HII}}
\right)\frac{dw_\alpha}{dt}
+(1-w_\alpha)\frac{d\mathcal R_{\ell/\alpha}^{\mathrm{HII}}}{dt}
+w_\alpha\frac{d\mathcal R_{\ell/\alpha}^{\mathrm{non-HII}}}{dt}.
\tag{39}
$$

The H-beta-normalized version substitutes \(w_\beta\). The first term is a **changing mixture**; the latter two are **evolving intrinsic spectra** caused, for example, by ionization parameter, stellar age, density, abundance, or shocks. In the deliberately fixed-template numerical experiment only the first term remains. A ratio rises if the increasingly weighted component has a larger intrinsic ratio; there is no rule that all forbidden/Balmer ratios must rise. [Citro et al. (2017)](#ref-citro) provide a quenching photoionization counterexample in which high-ionization lines can fade strongly.

Equation (34) also makes the absolute-luminosity test explicit. With fixed positive templates,

$$
\frac{d\mathcal L_\ell}{dt}
=\mathcal R_{\ell/B}^{\mathrm{HII}}
\frac{d\mathcal L_B^{\mathrm{HII}}}{dt}
+\mathcal R_{\ell/B}^{\mathrm{non-HII}}
\frac{d\mathcal L_B^{\mathrm{non-HII}}}{dt}.
\tag{40}
$$

A growing [N II]/H-alpha ratio does **not** require [N II] photons to become more numerous. If the weighted sum in equation (40) is negative, [N II] itself fades while its fraction relative to H-alpha rises. The numerical calculation checks the five individual line luminosities as well as the three ratios. With evolving templates, add \(\mathcal L_B^j\,d\mathcal R_{\ell/B}^j/dt\) for each component.

An observed NSF classification is **not** equivalent to \(w_\alpha\) crossing a universal threshold. The actual pipeline combines BPT domains, H-alpha EW, velocity width, Balmer detection, map quality, and H II SFR availability. The present model predicts fluxes and ratios; it does not turn them into an SF/NSF/ND count without applying that full observation operator and its noise.

# 6. Data-scaled numerical experiment

## 6.1 Inputs, calibration choices, and distinct branches

The calculation asks whether one **internally consistent choice** of the declared laws can produce (a) rapid loss of atomic supply and a slower molecular decline; (b) a small, temporary local SFR rise if efficiency increases; and (c) declining *absolute* Balmer and forbidden-line luminosities together with rising forbidden/Balmer ratios. It does **not** fit a trajectory to the Virgo stage bins. The common-support pre-peak SF H-alpha mean in Table 1 sets an initial luminosity *scale*; the post-peak NSF row is a comparison scale. It is essential that those rows are not treated as two observations of one patch at known times.

For the no-pulse branch, take \(\Sigma_{{\mathrm{HI}},0}=10\ M_\odot\,{\mathrm{pc}}^{-2}\), \(\tau_{\mathrm{dep,0}}=2.0\) Gyr, \(R=0.4\), \(\eta=0\), \(\gamma_{{\mathrm{strip}},{\mathrm{HI}}}=3.0\ {\mathrm{Gyr}}^{-1}\), and \(\gamma_{{\mathrm{strip}},{\mathrm{H_2}}}=2.22\ {\mathrm{Gyr}}^{-1}\) over an **illustrative** 1-Gyr interval. These are model choices; MAUVE does not independently measure their gas phase values in this bin. The constant coefficients are an exactly solvable effective history, not a measured constant RPS or a universal orbit. The direct molecular coefficient is chosen so that the model fades to the *order of magnitude* of the post-peak NSF H-alpha scale; that agreement is not an independent prediction.

Choose the initial young-photon destinations \(f_{{\mathrm{HII}},0}=0.855\), \(f_{{\mathrm{leak}},0}=0.095\), and \(f_{{\mathrm{unabs}},0}=0.050\), with \(f_{{\mathrm{unabs}},\infty}=0.20\). Thus \(95\%\) of the calibrated initial ionizing-photon budget contributes to the modelled H-alpha, and the compact component supplies \(90\%\) of that H-alpha. Set the initial *total* H-alpha to Table 1's \(1.561418\times10^{40}\ {\mathrm{erg\,s^{-1}\,kpc^{-2}}}\). Equations (2) and (29) then imply, step by step,

$$
\begin{aligned}
\Sigma_{{\mathrm{SFR}},0}
&=C_\alpha\,
\frac{\mathcal L_{\alpha,0}}
{f_{{\mathrm{HII}},0}+f_{{\mathrm{leak}},0}}
=0.08191\ M_\odot\,{\mathrm{yr}}^{-1}\,{\mathrm{kpc}}^{-2},\\
\Sigma_{{\mathrm{H_2}},0}
&=10^3\,\tau_{\mathrm{dep,0}}[{\mathrm{Gyr}}]\,
\Sigma_{{\mathrm{SFR}},0}[M_\odot\,{\mathrm{yr}}^{-1}{\mathrm{kpc}}^{-2}]
=163.82\ M_\odot\,{\mathrm{pc}}^{-2},\\
\tau_{\mathrm{conv}}
&=\frac{\Sigma_{{\mathrm{HI}},0}}
{(1-R+\eta)\Sigma_{{\mathrm{H_2}},0}/\tau_{\mathrm{dep,0}}}
=0.2035\ {\mathrm{Gyr}}.
\end{aligned}
\tag{41}
$$

The factor \(10^3\) only converts area and time units. The final line enforces pre-perturbation molecular balance from equation (18) with zero pre-perturbation stripping. At \(t=0^+\), the chosen environmental removal turns on. Therefore \(\gamma_{{\mathrm{loss}},{\mathrm{HI}}}=7.9146\ {\mathrm{Gyr}}^{-1}\), \(\gamma_{{\mathrm{loss}},{\mathrm{H_2}}}=2.52\ {\mathrm{Gyr}}^{-1}\), and \(\Sigma_{\Phi,0}=49.146\ M_\odot\,{\mathrm{pc}}^{-2}\,{\mathrm{Gyr}}^{-1}\). The inferred high initial H2 column is a **conditional** consequence of imposing a 2-Gyr depletion time on the bright selected H-alpha mean; it is not a CO measurement and may represent a small bright area rather than an entire annulus.

The **separate pulse branch** re-integrates equation (4) with equation (20), \(A_\epsilon=0.25\), \(t_{\mathrm{pk}}=60\) Myr, and \(m=9\), retaining the same gas-loss assumptions. It demonstrates the sign test in equation (9). The line calculation below uses the **no-pulse gas trajectory**, so that the temporary efficiency event is not silently fitted to BPT line values. The pulse shape and amplitude are illustrative and were chosen to produce an offset of the same order as the small NGC4654 facing-side residual, not inferred from that residual.

For the fixed line spectra choose the *linear* intrinsic ratios

**Table 2.** Conditional component spectra used in the numerical calculation; ratios are dimensionless.

| Physical component | [N II]/H-alpha | [S II]/H-alpha | [O III]/H-beta | Intrinsic H-alpha/H-beta |
|:--|--:|--:|--:|--:|
| Compact H II | 0.22 | 0.20 | 0.57 | 2.86 |
| Non-H II, young-leakage powered | 0.60 | 0.45 | 0.77 | 2.86 |

The numbers are **chosen rounded templates**, guided by the observed low- and high-ratio scales in Table 1. They are not directly measured physical-component spectra, outputs of a fitted photoionization model, or independent confirmation of this scenario. The stellar-response kernel is the explicitly illustrative exponential with \(\tau_{\mathrm{ion}}=3\) Myr. No old-star, shock, or AGN component is added to the plotted line luminosities.

Assigning the leakage-powered component a larger [O III]/H-beta ratio implicitly assumes a changed or filtered ionizing spectrum and nebular state; it does **not** follow merely from redirecting photons. [Belfiore et al. (2022)](#ref-belfiore22) find that spectral filtering can alter diffuse line ratios, while evolved stars are needed for some high-[O III] diffuse emission. The present template is therefore a conditional spectral choice. Its later failure for two observed rows below is informative.

## 6.2 Gas loss, a transient efficiency rise, and the meaning of the sign

We now insert the stated numbers into equations (11), (14), and (14a); this is a useful check on the units as well as the physical interpretation. For the no-pulse branch, \(\tau_{\mathrm{conv}}^{-1}=4.9146\ {\mathrm{Gyr}}^{-1}\), hence \(\gamma_{\mathrm{loss,HI}}=4.9146+3.0=7.9146\ {\mathrm{Gyr}}^{-1}\). The molecular consumption contribution is \((1-R+\eta)/\tau_{\mathrm{dep}}=0.6/2=0.30\ {\mathrm{Gyr}}^{-1}\), hence \(\gamma_{\mathrm{loss,H_2}}=0.30+2.22=2.52\ {\mathrm{Gyr}}^{-1}\). At \(t=1\) Gyr, the actual numerical substitution is

$$
\begin{aligned}
\Sigma_{\mathrm{HI}}(1)
&=10e^{-7.9146}
=0.003654\ M_\odot\,{\mathrm{pc}}^{-2},\\
\Sigma_{\mathrm{H_2}}(1)
&=163.820e^{-2.52}
+\frac{10}{0.203475}
\frac{e^{-7.9146}-e^{-2.52}}{2.52-7.9146}\\
&=13.911\ M_\odot\,{\mathrm{pc}}^{-2},\\
\Sigma_{\mathrm{SFR}}(1)
&=10^{-3}\frac{13.911}{2.0}
=0.006955\ M_\odot\,{\mathrm{yr}}^{-1}
\,{\mathrm{kpc}}^{-2}.
\end{aligned}
$$

The first term in the molecular line is surviving initial H2; the second is molecular gas replenished from H I and retained until 1 Gyr. At 300 Myr the same equations give \(\Sigma_{\mathrm{HI}}=0.931\), \(\Sigma_{\mathrm{H_2}}=80.35\ M_\odot\,{\mathrm{pc}}^{-2}\), and \(\Sigma_{\mathrm{SFR}}=0.04017\ M_\odot\,{\mathrm{yr}}^{-1}\,{\mathrm{kpc}}^{-2}\). The molecular and SFR ratios to their initial values at 1 Gyr are both \(0.0849\), as equation (2) requires for constant depletion time. Atomic supply is lost much faster, so \(\tau_\Phi=\Sigma_{\mathrm{H_2}}/\Sigma_\Phi\) becomes large and the positive supply term in equation (6) vanishes. This produces **remaining, suppressed SF**, not instantaneous extinction.

![Figure 1. Gas columns, no-pulse and pulse SFRs, and the pulse branch's instantaneous \(\Gamma_{\mathrm{SFR}}\) and accumulated \(\mathcal E\). Time is an illustrative local model coordinate, not an inferred orbital clock. The gas columns shown for the main line prediction are the no-pulse branch.](assets/20260929_detailed_gas_line_derivation/figure_01_reservoir_and_pulse.png)

The pulse branch starts by declining because environmental loss starts at \(t=0\). Its SFR first has a **positive time derivative** near 22 Myr, then reaches a maximum near 56 Myr. At that sampled maximum, direct insertion into equations (20) and (21a) shows the distinction between efficiency and SFR:

$$
\begin{aligned}
g(0.056)
&=0.25\left(\frac{0.056}{0.060}\right)^9
\exp\!\left[9\left(1-\frac{0.056}{0.060}\right)\right]
=0.24482,\\
\tau_{\mathrm{dep}}(0.056)
&=\frac{2.0}{1+0.24482}
=1.60666\ {\mathrm{Gyr}},\\
\Sigma_{\mathrm{SFR}}^{\mathrm{pulse}}(0.056)
&=10^{-3}\frac{144.1393}{1.60666}
=0.08971\ M_\odot\,{\mathrm{yr}}^{-1}
\,{\mathrm{kpc}}^{-2}.
\end{aligned}
$$

The pulse has left slightly *less* molecular gas than the simultaneous no-pulse branch (\(144.1393\) versus \(144.3217\ M_\odot\,{\mathrm{pc}}^{-2}\)) because higher efficiency also consumes gas. Nevertheless its SFR is \(1.0953\) times its initial value (\(+0.0395\) dex) and \(1.2432\) times the simultaneous no-pulse control. At the nearest 1-Myr grid point \(\Gamma_{\mathrm{SFR}}=+0.267\ {\mathrm{Gyr}}^{-1}\) while \(\mathcal E=-0.0910\). The grid's 56-Myr point lies just past the continuous SFR maximum, where \(\Gamma_{\mathrm{SFR}}=0\); it is already turning downward *instantaneously* but remains above its initial level *integrally*, exactly as equations (7)–(9) distinguish. A 25% peak **efficiency** increase does not imply a 25% rise over the initial **SFR**, because molecular mass has already been removed. The illustrative \(+0.0395\)-dex model rise has a similar order of magnitude to NGC4654's \(+0.0354\)-dex facing-side residual, but they have different baselines; this numerical proximity is not a physical fit or a detection significance.

## 6.3 Photon partition and absolute line fading

For the no-pulse branch, equations (29), (34), and (35) give a direct calculation of **each** line from the recent-SFR response. Let \(\mathcal P(t)=\overline{\Sigma}_{\mathrm{SFR,\alpha}}(t)/C_\alpha\), with units of H-alpha luminosity surface density if all modelled young photons were absorbed. Then

$$
\begin{aligned}
\mathcal L_\alpha
&=\mathcal P(f_{\mathrm{HII}}+f_{\mathrm{leak}}),&
\mathcal L_\beta
&=\frac{\mathcal P}{2.86}(f_{\mathrm{HII}}+f_{\mathrm{leak}}),\\
\mathcal L_{\mathrm{[N\,II]}}
&=\mathcal P(0.22f_{\mathrm{HII}}+0.60f_{\mathrm{leak}}),&
\mathcal L_{\mathrm{[S\,II]}}
&=\mathcal P(0.20f_{\mathrm{HII}}+0.45f_{\mathrm{leak}}),\\
\mathcal L_{\mathrm{[O\,III]}}
&=\frac{\mathcal P}{2.86}
\left(0.57f_{\mathrm{HII}}+0.77f_{\mathrm{leak}}\right).
\end{aligned}
\tag{42}
$$

The final line uses the common H-beta denominator explicitly. To see how the gas history enters an actual line prediction, consider the same 1-Gyr point. The age-response ODE gives \(\overline{\Sigma}_{\mathrm{SFR,\alpha}}=0.00700824\ M_\odot\,{\mathrm{yr}}^{-1}{\mathrm{kpc}}^{-2}\), so \(\mathcal P=0.00700824/C_\alpha=1.40627\times10^{39}\ {\mathrm{erg\,s^{-1}\,kpc^{-2}}}\). Equation (32) then gives \(f_{\mathrm{HII}}=0.072601\), \(f_{\mathrm{unabs}}=0.199945\), and the **remaining** absorbed fraction \(f_{\mathrm{leak}}=1-0.072601-0.199945=0.727454\). We therefore obtain

$$
\begin{aligned}
\mathcal L_\alpha^{\mathrm{HII}}(1)
&=(1.40627\times10^{39})(0.072601)
=1.0210\times10^{38},\\
\mathcal L_\alpha^{\mathrm{non-HII}}(1)
&=(1.40627\times10^{39})(0.727454)
=1.0230\times10^{39},\\
\mathcal L_{\mathrm{[N\,II]}}(1)
&=0.22(1.0210\times10^{38})
+0.60(1.0230\times10^{39})
=6.363\times10^{38},\\
\mathcal L_{\mathrm{[O\,III]}}(1)
&=\frac{0.57(1.0210\times10^{38})
+0.77(1.0230\times10^{39})}{2.86}
=2.958\times10^{38}.
\end{aligned}
$$

The luminosities in the last display are all in \({\mathrm{erg\,s^{-1}\,kpc^{-2}}}\). Dividing **after summation**, \(\mathcal L_\alpha=1.12509\times10^{39}\) and \(\mathcal L_\beta=\mathcal L_\alpha/2.86=3.93388\times10^{38}\), gives \(\mathcal L_{\mathrm{[N\,II]}}/\mathcal L_\alpha=0.5655\) and \(\mathcal L_{\mathrm{[O\,III]}}/\mathcal L_\beta=0.7519\). The [S II] calculation uses \(0.20\mathcal L_\alpha^{\mathrm{HII}}+0.45\mathcal L_\alpha^{\mathrm{non-HII}}\) and gives \(4.808\times10^{38}\), hence \(\mathcal L_{\mathrm{[S\,II]}}/\mathcal L_\alpha=0.4273\). All three forbidden-line **luminosities** are lower than at \(t=0\), even though their Balmer ratios are larger. These five expressions are the actual predictive part of the fixed-template experiment. They also expose what would have to change to improve the model: an evolving photon partition, an evolving intrinsic spectrum, a distinct non-H II power source, or gas absorption and radiative-transfer physics.

**Table 3.** No-pulse H-alpha components. Luminosities are in erg s\(^{-1}\) kpc\(^{-2}\), and \(w_\alpha\) is dimensionless.

| \(t\) (Gyr) | \(\mathcal L_\alpha\) | \(\mathcal L_\alpha^{\mathrm{HII}}\) | \(\mathcal L_\alpha^{\mathrm{non-HII}}\) | \(w_\alpha\) |
|--:|--:|--:|--:|--:|
| 0.0 | \(1.561\times10^{40}\) | \(1.405\times10^{40}\) | \(1.561\times10^{39}\) | 0.100 |
| 0.3 | \(6.611\times10^{39}\) | \(3.406\times10^{39}\) | \(3.205\times10^{39}\) | 0.485 |
| 0.6 | \(3.082\times10^{39}\) | \(7.637\times10^{38}\) | \(2.318\times10^{39}\) | 0.752 |
| 1.0 | \(1.125\times10^{39}\) | \(1.021\times10^{38}\) | \(1.023\times10^{39}\) | 0.909 |

**Table 4.** The other four absolute line luminosities from the same calculation, in erg s\(^{-1}\) kpc\(^{-2}\). [S II] is the doublet sum. Values in both tables are rounded independently from the saved full-precision calculation.

| \(t\) (Gyr) | \(\mathcal L_\beta\) | \(\mathcal L_{\mathrm{[N\,II]}}\) | \(\mathcal L_{\mathrm{[S\,II]}}\) | \(\mathcal L_{\mathrm{[O\,III]}}\) |
|--:|--:|--:|--:|--:|
| 0.0 | \(5.460\times10^{39}\) | \(4.028\times10^{39}\) | \(3.513\times10^{39}\) | \(3.221\times10^{39}\) |
| 0.3 | \(2.311\times10^{39}\) | \(2.672\times10^{39}\) | \(2.123\times10^{39}\) | \(1.542\times10^{39}\) |
| 0.6 | \(1.078\times10^{39}\) | \(1.559\times10^{39}\) | \(1.196\times10^{39}\) | \(7.763\times10^{38}\) |
| 1.0 | \(3.934\times10^{38}\) | \(6.363\times10^{38}\) | \(4.808\times10^{38}\) | \(2.958\times10^{38}\) |

![Figure 2. The model's disjoint destinations for one young-ionizing-photon budget and the resulting H-alpha components. Diffuse, non-H II H-alpha can rise early as the allocation changes; at 1 Gyr it too has faded. The chosen \(f_{\mathrm{unabs}}\) counts photons not absorbed in the modelled hydrogen and is not separately measured escape.](assets/20260929_detailed_gas_line_derivation/figure_02_photon_partition_and_fading.png)

At 1 Gyr, total H-alpha is \(7.21\%\) of its initial value. The compact H II H-alpha is \(0.73\%\) of its initial value, whereas the young-leakage-powered non-H II H-alpha is \(65.5\%\) of its initial value. Its *weight* consequently rises from \(0.100\) to \(0.909\). The non-H II luminosity is **not monotone**: it first brightens because \(f_{\mathrm{leak}}\) rises faster than the young source fades, and it later falls as the SFR loss dominates. This is more precise than assigning the non-H II component an independent fixed e-folding time. The absolute H-alpha, H-beta, [N II], [S II] doublet, and [O III] luminosities in Tables 3–4 all end below their initial values.

There is a further **gas-capacity requirement**, not guaranteed by the two *neutral* reservoirs. Equation (27) applied to the model's final non-H II H-alpha, \(1.023\times10^{39}\ {\mathrm{erg\,s^{-1}\,kpc^{-2}}}\), with \(h\nu_\alpha=3.027\times10^{-12}\) erg, \(p_\alpha=0.45\), and illustrative Case-B \(\alpha_B=2.6\times10^{-13}\ {\mathrm{cm^3\,s^{-1}}}\) at \(T_e\sim10^4\) K, implies

$$
{\mathrm{EM}}_{\mathrm{non-HII}}
\equiv\int n_en_p\,dz
=\frac{\mathcal L_\alpha^{\mathrm{non-HII}}}
{A_{\mathrm{kpc}}h\nu_\alpha\alpha_\alpha^{\mathrm{eff}}}
\simeq 98\ {\mathrm{pc\,cm^{-6}}},
\qquad
\alpha_\alpha^{\mathrm{eff}}=p_\alpha\alpha_B.
\tag{43}
$$

For a **uniform**, fully ionized hydrogen layer of depth 200 pc or 1 kpc, this corresponds respectively to \(n_e\simeq0.70\) or \(0.31\ {\mathrm{cm}}^{-3}\), and ionized-hydrogen columns of about 3.5 or \(7.7\ M_\odot\,{\mathrm{pc}}^{-2}\). Clumping, filling factor, helium, and geometry change these values. The two neutral equations do not track this ionized reservoir or its supply. Hence the plotted large final \(f_{\mathrm{leak}}\simeq0.727\) is **conditional on enough diffuse ionized gas and covering area remaining**; a measured lower emission measure would falsify this photon-partition closure. The approximate recombination coefficients and the Case-B restrictions come from [Hummer & Storey (1987)](#ref-hs87) and [Storey & Hummer (1995)](#ref-sh95); the numerical layer depths are illustrative, not measured in MAUVE.

Because the intrinsic Balmer decrements are equal, equation (36) reduces the three ratios to particularly transparent functions:

$$
\frac{\mathcal L_{\mathrm{[N\,II]}}}{\mathcal L_\alpha}
=0.22+0.38w_\alpha,\qquad
\frac{\mathcal L_{\mathrm{[S\,II]}}}{\mathcal L_\alpha}
=0.20+0.25w_\alpha,\qquad
\frac{\mathcal L_{\mathrm{[O\,III]}}}{\mathcal L_\beta}
=0.57+0.20w_\alpha.
\tag{44}
$$

For \(w_\alpha:0.100\rightarrow0.909\), these become \(0.258\rightarrow0.566\), \(0.225\rightarrow0.427\), and \(0.590\rightarrow0.752\), respectively. The [S II] ratio ends somewhat above Table 1's post-peak NSF mean \(0.417\), but the comparison has no independent test power after choosing the component ratios with the class means in view. No different “forbidden-line fading time” was inserted: the rising ratios result from the shared luminosity equation (42) and an increasing non-H II weight. Figure 3 plots both the declining absolute lines and rising ratios to prevent those two statements being conflated.

![Figure 3. All plotted absolute line luminosities end below their initial values, while the illustrated forbidden-to-Balmer ratios rise. A brief unresolved initial turn in one line is possible under the continuous closure. Dashed and dotted horizontal guides mark the observed pre-peak SF and post-peak NSF mean ratios in Table 1, respectively. The final panel shows the conditional BPT track, not an SF/NSF classification.](assets/20260929_detailed_gas_line_derivation/figure_03_lines_and_ratios.png)

## 6.4 Which observations are reproduced, and which are not

The modeled **end-point** H-alpha \(1.125\times10^{39}\ {\mathrm{erg\,s^{-1}\,kpc^{-2}}}\) is close to Table 1's post-peak NSF mean \(1.111\times10^{39}\), and its final ratios (0.566, 0.427, 0.752) are close to that row's (0.566, 0.417, 0.753). This checks arithmetic and internal plausibility at a MAUVE-relevant surface-luminosity scale. It is **not a blind fit**: the initial H-alpha, the approximate final H-alpha through the chosen loss coefficient, and both line-ratio templates were set with these data in view. In particular, the near equality of final [N II]/H-alpha is not new evidence for the assumed photon fractions.

The example provides a causal bridge for three *kinds* of trend: H I supply disappears early; molecular gas and true SFR remain but decline; and the compact-H II line fraction decreases faster than young-leakage-powered non-H II emission, so absolute lines fade while the weighted BPT ratios can rise. A separate short efficiency pulse shows how a local enhancement can coexist with net gas loss. It does **not** predict the observed *area fractions* of ND, NSF, and SF, their outer versus inner spatial pattern, the survivor-only stage offset \(0.407\), or a unique timescale between catalogue stages. These need radial initial-condition distributions, a pressure/orbit history, line detection and classification rules, and whole-galaxy statistical comparison.

The fixed two-template model also has a direct **failure** visible without an independent fit. Any positive mixture in equation (36) with the same Balmer decrement lies between its component [O III]/H-beta ratios \(0.57\) and \(0.77\). Table 1's *pre-peak NSF* mean is \(0.898\) and its *post-peak SF* mean is \(0.210\), both outside that interval. No photon partition can reproduce either row with Table 2's fixed [O III] templates. The failure points to changing intrinsic spectra, additional excitation, different attenuation/decrements, and/or sample selection, not to a contradiction of the luminosity-conservation algebra. A future quantitative model should confront all four rows and their covariance, not only the chosen two anchors.

# 7. What the calculation can establish

Equations (4)–(17) connect local supply, stripping, molecular survival, and the exact condition for a falling or rising SFR. Equations (22)–(33) connect that gas history to the *available young ionizing photons* and their separate absorption destinations. Equations (34)–(44) connect component luminosities to the individual observed lines and ratios without assuming that a high forbidden/Balmer ratio means a growing forbidden-line luminosity. The numerical example demonstrates that these three layers can coexist at approximately the measured MAUVE line scale, while explicitly rejecting a universal slow-fading non-H II component and a unique RPS interpretation.

The most useful next test is a **forward observational comparison**: assign \(\Sigma_{{\mathrm{HI}},0}(\boldsymbol{x})\), \(\Sigma_{{\mathrm{H_2}},0}(\boldsymbol{x})\), and a spatially and temporally varying \(\gamma_{{\mathrm{strip}},i}(\boldsymbol{x},t)\); calculate all five line maps and continuum; convolve to MAUVE resolution; apply the exact S/N, EW, width, BPT, and H II SFR rules; then compare SF/NSF/ND *occupancy* and survivor-only intensity at equal-galaxy weight. CO and H I maps, shock-sensitive lines and kinematics, EW, stellar-population age constraints, and sensitivity to diffuse-photon propagation would distinguish the possible photon-budget terms. Without that forward operator, a rising non-H II **luminosity fraction** is only a mechanism-level prediction, not a prediction of the measured NSF **area fraction**.

# 8. Limitations and falsifiable consequences

First, the two-reservoir ODE neglects in-plane transport between neighbouring apertures, replenishment after onset, phase dissociation, gas re-accretion, and feedback delay; stripping remains an explicit effective loss across the aperture boundary. Its \(\tau_{\mathrm{conv}}\), \(\gamma_{{\mathrm{strip}},i}\), and \(\tau_{\mathrm{dep}}\) are effective resolved-scale parameters. A constant \(\gamma\) over the example is a solvable exposure, not an orbital-pressure measurement. Molecular abundance is especially degenerate with \(\tau_{\mathrm{dep}}\) in optical-only data. A physical RPS model would compute \(\rho_{\mathrm{ICM}}(t)\), \(v_{\mathrm{rel}}(t)\), orientation, restoring force, and gas response before assigning a local loss rate. The illustrative 1 Gyr is not dated by the three stage bins.

Second, a high line ratio cannot alone identify HOLMES, leakage, shocks, or AGN. The numerical line model omits old and mechanical sources and uses a gas-dependent compact/diffuse partition without solving radiative transfer, **photon exchange between resolved apertures**, or ionized-gas continuity. Its final non-H II luminosity requires the emission measure in equation (43), which remains unverified. Its fixed component spectra ignore metallicity, ionization parameter, density, stellar age, dust, and differential Balmer decrements. The approximate Case-B value 2.86 and exponential age kernel are controlled simplifications. Equation (31)'s evolved-star ceiling uses one published ionizing yield and a maximal gas-absorption fraction; it should not be applied unchanged to another stellar population or to an EW-selected subset. The pre-peak NSF and post-peak SF [O III]/H-beta failures are already evidence that the chosen two fixed spectra cannot describe every measured class.

Third, detection and class selection are absent from the calculation. A region can enter ND as lines cross a sensitivity threshold, or enter NSF because EW, width, or a diagnostic boundary changes, even when it contains some star formation. Conversely, the SF class is not a pure H II photon channel. The Table 1 rows condition on five-line detected support and differ in galaxy membership; their ratios are ratios of mean corrected luminosities. Noise, dust-correction covariance, stellar-continuum subtraction, and whole-galaxy rather than pixel resampling matter for a quantitative comparison. The NGC4654 residual has no established galaxy-level significance or orientation-error propagation, and it cannot select the pulse physics.

The closure has testable consequences rather than a claim of unique identification. If compact absorption truly tracks \(\Sigma_{\mathrm{H_2}}/\Sigma_{{\mathrm{H_2}},0}\) as in equation (32), then at fixed recent SFR the H II H-alpha fraction should fall with molecular gas and the non-H II fraction should grow *only while gas remains available to absorb leaked photons*. If an independently measured diffuse-gas emission measure falls too rapidly to support equation (29), the closure fails. If observed [O III]/H-beta moves beyond the allowed interval of the independently measured component spectra, equation (36) demands evolving spectra, differing Balmer decrements, or another component. If the local enhancement is efficiency driven, \(\Sigma_{\mathrm{SFR}}/\Sigma_{\mathrm{H_2}}\) should rise on that side; if it is content driven, the SFR increase can occur at unchanged \(\tau_{\mathrm{dep}}\). These are discriminants that the current optical residual alone cannot settle.

# Appendix A. Definitions, units, and provenance of every model symbol

Throughout, a subscript 0 denotes a value immediately before the modelled perturbation, and \(\infty\) denotes a chosen limiting value of one closure, not an observed late-stage limit. Surface densities use a common projected area; gas columns are usually printed in \(M_\odot\,{\mathrm{pc}}^{-2}\), whereas SFR and line luminosity surface densities are printed per kpc². The required conversion is stated below equation (1). A subscript \(i\) means either H I or H2; \(j\) means H II or non-H II; \(\ell\) is one specified emission line; \(B\) is its appropriate Balmer denominator.

**Table A1.** Coordinates, reservoirs, and environmental forcing.

| Symbol | Detailed definition, units, and status |
|:--|:--|
| \(\boldsymbol{x}=(r,\varphi)\), \(t\), \(t'\), \(u\), \(a\), \(z\) | Position in the disc, elapsed local-model time (Gyr unless stated), integration dummy times, stellar cohort age, and line-of-sight distance (cm inside recombination integrals). No stage-to-time conversion is inferred. |
| \(\Sigma_{\mathrm{HI}}\), \(\Sigma_{\mathrm{H_2}}\) | Local atomic-associated and molecular-associated **neutral** gas columns, \(M_\odot\,{\mathrm{pc}}^{-2}\), on one helium convention. H2 denotes molecular gas associated with the modelled SFR; neither column tracks the ionized gas needed for equation (43). |
| \(\Sigma_*\) | Stellar-mass surface density, \(M_\odot\,{\mathrm{kpc}}^{-2}\) when used in the old-star photon budget. The example uses the observed stellar-density bin centre, not a pixel-specific stellar map. |
| \(\Sigma_\Phi\), \(\tau_\Phi\) | Supply rate surface density from the atomic-associated reservoir into H2 and molecular supply time \(\Sigma_{\mathrm{H_2}}/\Sigma_\Phi\); \(M_\odot\,{\mathrm{pc}}^{-2}\,{\mathrm{Gyr}}^{-1}\) and Gyr in the numerical calculation. Notation follows [Huang et al. (2026)](#ref-huang); the equality \(\Sigma_\Phi=\Sigma_{\mathrm{HI}}/\tau_{\mathrm{conv}}\) is this report's closure. |
| \(\tau_{\mathrm{conv}}\) | Effective transfer time of the modelled H I-associated reservoir into the molecular one, Gyr. It is fitted here only through chosen initial balance, not directly measured or equated with a microscopic grain-formation time. |
| \(\gamma_{{\mathrm{strip}},i}\) | Environmental fractional removal coefficient for phase \(i\), Gyr\(^{-1}\). The same indexed symbol is used for H I and H2. It is neither \(P_{\mathrm{ram}}\) nor a universal function of it. |
| \(\gamma_{{\mathrm{loss}},i}\) | Composite fractional loss coefficient in the restricted ODE, Gyr\(^{-1}\); sums in equation (10), **not additional parameters**. H I includes phase transfer and stripping; H2 includes locking/feedback and stripping. |
| \(\mathcal G_i(t)\), \(\mathcal G_{\mathrm{H_2}}^{\mathrm{pulse}}(t)\) | Dimensionless integrated loss exposure \(\int_0^t\gamma_{{\mathrm{loss}},i}(u)\,du\); its explicitly evaluated molecular-pulse version appears before equation (21a). Neither is a new removal coefficient. |
| \(I(t)\), \(I_{\mathrm{H_2}}(t)\) | Dimensionless integrating factors \(e^{\gamma_{\mathrm{loss,H_2}}t}\) in the constant-rate derivation and \(e^{\mathcal G_{\mathrm{H_2}}(t)}\) in the time-varying derivation. They are algebraic devices used to integrate the H2 balance. |
| \(P_{\mathrm{ram}}\), \(\rho_{\mathrm{ICM}}\), \(v_{\mathrm{rel}}\) | Incident ram pressure, intracluster-medium mass density, and gas–ICM relative speed; pressure, mass/volume, and speed. Equation (3) is a forcing scale, not a specified law for \(\gamma_{{\mathrm{strip}},i}\). |

**Table A2.** Star formation and the optional local-efficiency event.

| Symbol | Detailed definition, units, and status |
|:--|:--|
| \(\Sigma_{\mathrm{SFR}}\) | True local stellar-mass formation rate per area, \(M_\odot\,{\mathrm{yr}}^{-1}{\mathrm{kpc}}^{-2}\). In equation (2) it is the instantaneous gas-law rate, not yet the H-alpha-inferred rate averaged over young-stellar ages. |
| \(\tau_{\mathrm{dep}}\), \(\tau_{\mathrm{dep,0}}\) | Molecular depletion time \(\Sigma_{\mathrm{H_2}}/\Sigma_{\mathrm{SFR}}\) in consistent units, and its initial value; Gyr. A smaller value means higher molecular efficiency, \(1/\tau_{\mathrm{dep}}\). |
| \(R\), \(\eta\) | Dimensionless effective promptly recycled fraction returned to the modelled cold phase and feedback mass-loading factor. The combination \(1-R+\eta\) is molecular gas removed or locked per formed stellar mass; its exact phase return is a declared approximation. |
| \(\Gamma_{\mathrm{SFR}}\), \(\mathcal E\) | Instantaneous negative logarithmic SFR derivative, Gyr\(^{-1}\), and its dimensionless time integral. Positive \(\Gamma_{\mathrm{SFR}}\) means falling *at that instant*; positive \(\mathcal E\) means below the initial SFR. They are diagnostic definitions, not new environmental coefficients. |
| \(g(t)\), \(A_\epsilon\), \(t_{\mathrm{pk}}\), \(m\) | Dimensionless time profile of a chosen molecular-efficiency pulse, its dimensionless peak amplitude, peak time, and dimensionless shape index \(m>1\) for a smooth zero-slope onset; all defined in equation (20). None is a measured ram-pressure law. |
| \(\Sigma_{\mathrm{H_2}}^{\mathrm{pulse}}\), \(\Sigma_{\mathrm{SFR}}^{\mathrm{pulse}}\) | Molecular column and true instantaneous SFR when the specified pulse changes \(\tau_{\mathrm{dep}}(t)\), with the same units as their no-pulse counterparts. The superscript marks a separate illustrative branch, not another physical gas phase. |

**Table A3.** Young and old ionizing photons, recombination, and gas capacity.

| Symbol | Detailed definition, units, and status |
|:--|:--|
| \(q_{\mathrm{H,y}}(a)\), \(Q_{\mathrm{H,y}}\) | Hydrogen-ionizing photon production rate per unit **formed young-stellar mass** at age \(a\), photons s\(^{-1}M_\odot^{-1}\), and its summed production rate per projected area, photons s\(^{-1}{\mathrm{kpc}}^{-2}\). Population synthesis sets the former; equation (22) sums cohorts. |
| \(K_\alpha(a)\), \(\tau_{\mathrm{ion}}\), \(\overline{\Sigma}_{\mathrm{SFR,\alpha}}\) | Unit-normalized age-response kernel (inverse time), illustrative exponential mean response time (3 Myr), and resulting recent-SFR equivalent, in the same SFR surface units as \(\Sigma_{\mathrm{SFR}}\). An actual synthesis kernel is not claimed to be exponential. |
| \(C_\alpha\) | MAUVE pipeline conversion between fully young-powered H-alpha luminosity and SFR, \(4.983582\times10^{-42}\ M_\odot\,{\mathrm{yr}}^{-1}/({\mathrm{erg\,s}}^{-1})\); the ratio of the corresponding *surface* quantities uses identical projected-area units. The calibration inherits its stellar-population assumptions from [Kennicutt & Evans (2012)](#ref-ke12) and the pipeline's Chabrier factor. |
| \(n_e\), \(n_p\), \(T_e\), \(z\) | Electron and proton densities (cm\(^{-3}\)), electron temperature (K), and distance through emitting gas (cm); their product integrated along \(z\) is the emission measure. |
| \(\alpha_B\), \(\alpha_\alpha^{\mathrm{eff}}\), \(p_\alpha\) | Total hydrogen Case-B recombination coefficient and effective H-alpha coefficient (cm³ s\(^{-1}\)), and their dimensionless ratio, H-alpha photons per absorbed ionization in the stated balance. Adopted values depend on \(T_e\) and density; see [Hummer & Storey (1987)](#ref-hs87) and [Storey & Hummer (1995)](#ref-sh95). |
| \(h\nu_\alpha\), \(j_\alpha\), \(A_{\mathrm{kpc}}\) | H-alpha photon energy (erg), volume emissivity (erg s\(^{-1}\) cm\(^{-3}\)), and surface-unit conversion factor \(9.521\times10^{42}\ {\mathrm{cm^2\,kpc^{-2}}}\). Multiplication by this factor converts a surface flux per cm² into one per kpc² in equations (26)–(27) and (43). |
| \(Q_{\mathrm{H,abs}}\), \({\mathrm{EM}}_{\mathrm{non-HII}}\) | Hydrogen-ionizing photon absorption rate **per kpc²** in photons s\(^{-1}\) kpc\(^{-2}\), and non-H II emission measure \(\int n_en_p dz\), often printed in pc cm\(^{-6}\). The latter is a physical gas-capacity requirement, not solved by the neutral-reservoir ODE. |
| \(f_{\mathrm{HII}}\), \(f_{\mathrm{leak}}\), \(f_{\mathrm{unabs}}\) | Mutually exclusive fractions of young-star photons produced **inside the modelled aperture**: compact H II absorption there, absorption in diffuse non-H II gas there after leakage, and photons not absorbed by that aperture's modelled hydrogen (escape from it, dust absorption, or another destination). Their sum is unity under the local no-import assumption. A fraction can change while the source photon rate declines. |
| \(f_{{\mathrm{HII}},0}\), \(f_{{\mathrm{unabs}},0}\), \(f_{{\mathrm{unabs}},\infty}\) | Initial compact fraction, initial unabsorbed fraction, and the chosen limiting unabsorbed fraction in equation (32); dimensionless closure inputs, not measured escape fractions. The initial leaked fraction follows from unity. |
| \(q_{\mathrm{H,old}}\), \(f_{\mathrm{abs,old}}\) | Old-stellar ionizing photon rate per unit **present stellar mass**, photons s\(^{-1}M_\odot^{-1}\), and fraction of that independent budget absorbed by gas. The old-star numerical ceiling uses the adopted yield and \(f_{\mathrm{abs,old}}=1\); old stars are zero in the plotted trajectory. |

**Table A4.** Line luminosities, ratios, and component labels.

| Symbol | Detailed definition, units, and status |
|:--|:--|
| \(\mathcal L_\ell^j\), \(\mathcal L_\ell\) | Surface luminosity of line \(\ell\) from component \(j\) and summed over components, erg s\(^{-1}\) kpc\(^{-2}\). \(\alpha\) and \(\beta\) denote H-alpha and H-beta; [N II] is \(\lambda6583\), [S II] is the \(\lambda6716+\lambda6730\) doublet, and [O III] is \(\lambda5007\). |
| H II, non-H II, SF, NSF, ND | H II and non-H II are **physical line-component labels** in the model. SF, NSF, and ND are **pipeline observational categories** with line S/N, BPT, EW, width, and SFR criteria. Their labels must not be equated one-to-one. |
| \(\mathcal L_\alpha^{\mathrm{old}}\), \(\mathcal L_\alpha^{\mathrm{mech}}\), \(\mathcal L_\alpha^{\mathrm{AGN}}\) | Old-stellar, mechanically, and AGN-powered non-H II H-alpha luminosities per area. They show possible extra sources in equation (30); only the conditional old-star ceiling is evaluated, and all are set to zero in the plotted two-component calculation. |
| \(\mathcal R_{\ell/B}^{j}\), \(\mathcal R_{\ell/B}\), \(B_j\) | Intrinsic **linear** line-to-Balmer ratio in component \(j\), observed mixture ratio after adding luminosities, and intrinsic H-alpha/H-beta decrement of component \(j\). They are dimensionless. The fixed templates in Table 2 are assumptions. |
| \(w_\alpha\), \(w_\beta\) | Fraction of total H-alpha or H-beta surface luminosity from the physical non-H II component. They are distinct if component Balmer decrements differ. Neither is a measured NSF area fraction. |
| \(\mathcal P(t)\) | Shorthand \(\overline{\Sigma}_{\mathrm{SFR,\alpha}}/C_\alpha\), in erg s\(^{-1}\) kpc\(^{-2}\): hypothetical fully absorbed young-powered H-alpha luminosity density before partition among photon destinations. This is a derived quantity, not a new reservoir. |

# Appendix B. Numerical audit and source-to-equation ledger

The reproducible calculation is [model_predictions.py](assets/20260929_detailed_gas_line_derivation/model_predictions.py), with complete [model parameters](assets/20260929_detailed_gas_line_derivation/model_parameters.json), [1-Myr model trajectories](assets/20260929_detailed_gas_line_derivation/numerical_predictions.csv), [MAUVE anchor values](assets/20260929_detailed_gas_line_derivation/mauve_anchor_table.csv), and [numerical checks](assets/20260929_detailed_gas_line_derivation/numerical_checks.json). It reads the 14 September line-profile scalar export rather than reprocessing large map cubes. The six recorded source fingerprints of the three 9 September stage notebooks, the MAUVE SFR pipeline, master catalogue, and effective-radii table were rechecked unchanged; see [fingerprints](assets/20260929_detailed_gas_line_derivation/source_fingerprints.json). The large FITS inputs to the prior extraction were **not** all rehashed. The 28 September [NGC4654 audit](assets/20260929_detailed_gas_line_derivation/ngc4654_audit.json) reran only the required data and estimator cells of the explicitly named 14 September gradient notebook, whose current SHA-256 is recorded there; plotting/gallery cells, PA uncertainty propagation, and a post-fit continuum \(S/N>25\) cut were not included in that gradient run.

**Table B1.** Verification checks from the saved no-pulse and pulse calculations.

| Check | Computed result | Interpretation |
|:--|--:|:--|
| Exact two-reservoir solution versus tight-tolerance ODE | maximum relative error \(2.24\times10^{-12}\) | Confirms equations (11) and (14) were implemented as written. |
| Three young-photon fractions form an admissible partition | maximum sum error \(1.11\times10^{-16}\); minimum fraction \(0.050\) | Checks both normalization and nonnegativity under equations (28) and (32). |
| H-alpha components add to total | maximum relative error \(0\) | Checks line-component arithmetic. |
| Molecular and SFR decline ratio at 1 Gyr | both \(0.08491\) | Checks the fixed-\(\tau_{\mathrm{dep}}\) law in equation (2). |
| Kernel versus instantaneous SFR at 1 Gyr | relative difference \(0.761\%\) | Only the chosen smooth 3-Myr response, not a universal H-alpha-kernel error. |
| Recombination gas-capacity requirement | \({\mathrm{EM}}\simeq98\ {\mathrm{pc\,cm^{-6}}}\) | Physical condition external to the neutral ODE; not verified by the optical class means alone. |

**Table B2.** Scientific provenance of the main relations.

| Equations | Basis in the literature | What this report adds |
|:--|:--|:--|
| (1)–(2), (4)–(18), including (14a), (15a), (17a) | The MAUVE molecular-supply regulator of [Huang et al. (2026)](#ref-huang), broader gas-regulator framing of [Lilly et al. (2013)](#ref-lilly), and time-dependent stripping context from [Köppen et al. (2018)](#ref-koppen) and [Singh et al. (2019)](#ref-singh). | A restricted atomic reservoir, declared stripping sinks, exact constant- and variable-coefficient solutions, explicit SFR histories, and the local sign test. No paper is credited with these particular numerical coefficients. |
| (19)–(21a) | Resolved molecular-content and efficiency interpretations in [Brown et al. (2023)](#ref-brown), [Lizee et al. (2021)](#ref-lizee), [Vollmer et al. (2012)](#ref-vollmer12), and [Nehlig et al. (2016)](#ref-nehlig). | A chosen smooth depletion-time pulse, its cumulative loss exposure, and the causal H2/SFR solution. It is not a published pressure-efficiency law. |
| (22)–(25), (41) | H-alpha age sensitivity and SFR conversion from [Kennicutt & Evans (2012)](#ref-ke12), plus the live MAUVE Chabrier-adjusted coefficient; spatial-scale caveat from [Kruijssen & Longmore (2014)](#ref-kl14). | Cohort-convolution derivation, an explicitly illustrative exponential response, and data-scale initialization. |
| (26)–(31), (43) | Case-B emissivities and coefficients from [Hummer & Storey (1987)](#ref-hs87) and [Storey & Hummer (1995)](#ref-sh95); diffuse and evolved-star source models from [Belfiore et al. (2022)](#ref-belfiore22). | Disjoint photon fractions, an old-star conditional luminosity ceiling, and the emission-measure capacity check. |
| (32)–(40), (42), (44) | Diffuse-source mixtures and line-ratio shifts discussed by [Belfiore et al. (2022)](#ref-belfiore22) and [Zhang et al. (2017)](#ref-zhang); different quenching-line responses from [Citro et al. (2017)](#ref-citro). | A testable gas-dependent photon-partition closure, exact component-mixture algebra, and MAUVE-scale numerical templates. The fixed ratios were chosen, not derived from those papers. |

# References

<a id="ref-huang"></a>
Huang, R., et al. (2026), *MAUVE-MUSE: When Metallicity Follows or Fights Star Formation—A Mass-Dependent Inversion in Virgo Galaxies*. [Published DOI](https://doi.org/10.1093/mnras/stag1019); [author manuscript](https://arxiv.org/html/2605.31412v1). Equations cited above refer to the manuscript version.

<a id="ref-brown"></a>
Brown, T., et al. (2023), *VERTICO VII: Environmental quenching caused by suppression of molecular gas content and star formation efficiency in Virgo Cluster galaxies*. [Author manuscript](https://arxiv.org/html/2308.10943v1).

<a id="ref-koppen"></a>
Köppen, J., Jáchym, P., Taylor, R., & Palouš, J. (2018), *Ram Pressure Stripping Made Easy: An Analytical Approach*, MNRAS, 479, 4367. [DOI](https://doi.org/10.1093/mnras/sty1610); [author manuscript](https://arxiv.org/abs/1806.05887).

<a id="ref-singh"></a>
Singh, A., Gulati, M., & Bagla, J. S. (2019), *Ram pressure stripping: an analytical approach*, MNRAS, 489, 5582–5593. [DOI](https://doi.org/10.1093/mnras/stz2523); [publisher article](https://academic.oup.com/mnras/article/489/4/5582/5575211).

<a id="ref-lilly"></a>
Lilly, S. J., et al. (2013), *Gas-regulation of galaxies: the evolution of the cosmic sSFR, the metallicity–mass–SFR relation and the stellar content of haloes*, ApJ, 772, 119. [DOI](https://doi.org/10.1088/0004-637X/772/2/119).

<a id="ref-lizee"></a>
Lizée, T., Vollmer, B., Braine, J., & Nehlig, F. (2021), *Gas compression and stellar feedback in the tidally interacting and ram-pressure stripped Virgo spiral galaxy NGC 4654*, A&A, 645, A111. [DOI](https://doi.org/10.1051/0004-6361/202038910); [author manuscript](https://arxiv.org/abs/2011.10531).

<a id="ref-vollmer12"></a>
Vollmer, B., Wong, O. I., Braine, J., Chung, A., & Kenney, J. D. P. (2012), *The influence of the cluster environment on the star formation efficiency of 12 Virgo spiral galaxies*, A&A, 543, A33. [DOI](https://doi.org/10.1051/0004-6361/201118690).

<a id="ref-nehlig"></a>
Nehlig, F., Vollmer, B., & Braine, J. (2016), *Effects of environmental gas compression on the multiphase ISM and star formation*, A&A, 587, A108. [DOI](https://doi.org/10.1051/0004-6361/201527021).

<a id="ref-ke12"></a>
Kennicutt, R. C., & Evans, N. J. (2012), *Star Formation in the Milky Way and Nearby Galaxies*, ARA&A, 50, 531–608. [DOI](https://doi.org/10.1146/annurev-astro-081811-125610); [author manuscript](https://arxiv.org/html/1204.3552).

<a id="ref-kl14"></a>
Kruijssen, J. M. D., & Longmore, S. N. (2014), *An uncertainty principle for star formation. I. Why galactic star formation relations break down below a certain spatial scale*, MNRAS, 439, 3239. [DOI](https://doi.org/10.1093/mnras/stu098).

<a id="ref-hs87"></a>
Hummer, D. G., & Storey, P. J. (1987), *Recombination-line intensities for hydrogenic ions. I. Case B calculations for H I and He II*, MNRAS, 224, 801–820. [DOI](https://doi.org/10.1093/mnras/224.3.801).

<a id="ref-sh95"></a>
Storey, P. J., & Hummer, D. G. (1995), *Recombination line intensities for hydrogenic ions—IV. Total recombination coefficients and machine-readable tables for Z = 1 to 8*, MNRAS, 272, 41–48. [DOI](https://doi.org/10.1093/mnras/272.1.41).

<a id="ref-belfiore22"></a>
Belfiore, F., et al. (2022), *A tale of two DIGs: The relative role of H II regions and low-mass hot evolved stars in powering the diffuse ionised gas in PHANGS-MUSE galaxies*, A&A, 659, A26. [DOI](https://doi.org/10.1051/0004-6361/202141859); [author manuscript](https://arxiv.org/html/2111.14876v3).

<a id="ref-belfiore16"></a>
Belfiore, F., et al. (2017), *SDSS-IV MaNGA—The spatially resolved transition from star formation to quiescence*, MNRAS, 466, 2570–2589. [DOI](https://doi.org/10.1093/mnras/stw3211); [author preprint](https://arxiv.org/abs/1605.07189).

<a id="ref-zhang"></a>
Zhang, K., et al. (2017), *SDSS-IV MaNGA: The impact of diffuse ionized gas on emission-line ratios, interpretation of diagnostic diagrams, and gas metallicity measurements*, MNRAS, 466. [DOI](https://doi.org/10.1093/mnras/stw3308); [author preprint](https://arxiv.org/abs/1612.02000).

<a id="ref-citro"></a>
Citro, A., Pozzetti, L., Quai, S., Moresco, M., Vallini, L., & Cimatti, A. (2017), *A methodology to select galaxies just after the quenching of star formation*, MNRAS. [DOI](https://doi.org/10.1093/mnras/stx932); [author preprint](https://arxiv.org/abs/1704.05462).
