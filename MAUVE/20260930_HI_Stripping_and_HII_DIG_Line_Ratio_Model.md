---
title: "H I stripping, delayed star formation, and H II–DIG line mixing"
subtitle: "Detailed analytical derivation and conditional numerical tests for MAUVE"
author: "Research report for Rongjun Huang"
date: "30 September 2026"
lang: en
---

# 1. What we want this revision to establish

Here we want to examine a simpler physical chain: ram-pressure stripping (RPS) removes atomic gas, the local molecular reservoir receives less replenishment, and star formation subsequently declines as molecular gas is consumed. In the main calculation, we set **direct ram-pressure removal of molecular gas to zero**. We then ask separately whether fading young-star emission can change the mixture of H II-region and diffuse emission enough to explain increasing forbidden-to-Balmer ratios. These are related questions, but the gas continuity equations alone do not specify how ionizing photons are distributed between compact and diffuse gas.

The main result is that H I-only stripping can produce a delayed, eventually approximately exponential star-formation-rate (SFR) decline. However, with a normal-disc depletion time, this decline can be substantially slower than in the previous calculation that also removed H2 directly. Likewise, an increasing non-H II luminosity weight can increase a line ratio when the non-H II component has the higher intrinsic ratio. But if both components are powered by the same massive O- and B-type (OB) stars and their photon-allocation fractions stay constant, common fading changes their luminosities without changing their relative weight. We therefore retain **constant component spectra as the main testable hypothesis**, while making the required evolution of the luminosity weight explicit. Evolving gas coefficients and component spectra are discussed in the appendices.

## 1.1 Observational constraints and their status

The existing classes are star-forming (SF), non-SF-classified but Balmer-detected (NSF), and Balmer-nondetected (ND): ND fails the joint Hα/Hβ detection requirement, and NSF is the remaining Balmer-detected population that fails the full SF selection. Our resolved stage analyses motivate four constraints: outer/low-stellar-density regions increasingly occupy ND; inner/high-stellar-density regions increasingly occupy NSF; the surviving SF population has lower SFR intensity; and detected NSF emission can become fainter while some forbidden-to-Balmer ratios increase. A possible local enhancement is an additional, spatially restricted signal rather than a necessary property of every patch. These are observational categories and population comparisons. **SF is not a pure H II spectrum, NSF is not a pure diffuse ionized gas (DIG) spectrum, and ND is not a direct measurement of absent gas.**

For a numerical reference, Table 1 uses the saved five-line common-detection-support export at \(\log_{10}[\Sigma_*/(M_\odot\,\mathrm{kpc}^{-2})]=8.625\). The underlying analysis imposes \(S/N_{\mathrm{POSTFIT}}>25\). The entries are ratios of galaxy-averaged line surface luminosities on the common support, rather than averages of individual-spaxel ratios. The [S II] numerator is the sum of the 6716 and 6731 Å lines. The source notebook/pipeline fingerprints were checked again on 30 September before regenerating this table from the existing scalar export. This check does not constitute a fresh reduction of the large input FITS maps.

**Table 1.** Representative MAUVE scales. \(\mathcal L_\alpha\) is in \(10^{39}\,\mathrm{erg\,s^{-1}\,kpc^{-2}}\); ratios are linear, dimensionless values.

| Stage | Category | Galaxies | \(\mathcal L_\alpha\) | [N II]/Hα | [S II]/Hα | [O III]/Hβ |
|:--|:--|--:|--:|--:|--:|--:|
| Pre-peak | SF | 7 | 15.614 | 0.254 | 0.222 | 0.596 |
| Close-to-peak | SF | 5 | 9.892 | 0.317 | 0.258 | 0.251 |
| Post-peak | SF | 12 | 7.580 | 0.319 | 0.228 | 0.210 |
| Pre-peak | NSF | 6 | 3.558 | 0.373 | 0.377 | 0.898 |
| Close-to-peak | NSF | 5 | 2.077 | 0.437 | 0.386 | 0.468 |
| Post-peak | NSF | 13 | 1.111 | 0.566 | 0.417 | 0.753 |

Several distinctions matter. In this bin, NSF [N II]/Hα and [S II]/Hα increase from pre- to post-peak, whereas [O III]/Hβ decreases between those endpoints and is non-monotonic across all three stages. Thus “higher BPT ratios,” referring to Baldwin–Phillips–Terlevich diagnostic diagrams, is too broad a description of these particular numbers. The approximately \(0.407\) post-peak surviving-SF attenuation factor from the earlier fit is a separate, matched-profile statistic; it is not the Table 1 SF luminosity ratio and is not an NSF attenuation factor. Its stored 16th–84th percentile whole-galaxy bootstrap interval is approximately \(0.341\)–\(0.504\). We use it only as a scale against which to compare conditional gas-model attenuation.

We retain the pre-peak sample as a **field-like internal-control assumption**. Virgo membership prevents us from treating it as an independently established pristine field population. Further, the stage bins contain different galaxies and different surviving spatial support. They do not follow the same region through time. The model variable \(t\) below is elapsed time under a specified local history, not a measured infall age assigned to a stage.

## 1.2 How the two linked discussions enter this report

The textual contents of [“Summarize stripping trends”](https://chatgpt.com/c/6abccf3b-89e4-83ec-a49a-3c7018e6f05d) and [“分支 · Model Dependence Clarification”](https://chatgpt.com/c/6abbd62b-d27c-83ec-8ff0-fc933bdbc716) were read through the connected conversation interface. The first discussion motivates the H I → H2 → stars supply chain and cautions against interpreting NSF as a physical endmember. The second distinguishes changing luminosity weights from evolving component spectra and proposes tests using H II-dominated regions and several ratios simultaneously. These discussions define the questions; they are not primary scientific evidence. Image attachments in the second conversation were not available through its textual retrieval. Numerical statements here therefore use the local scalar products in Table 1 rather than visual claims about those attachments.

The algebra is written in the style of the November 2025 regulator note: define the quantities first, substitute explicitly, carry out each integration, and distinguish an instantaneous condition from a sustained state. Literature-derived ingredients are cited where introduced. Algebra obtained by combining the stated assumptions is identified as a derivation of this report, rather than attributed to a paper that did not publish it.

# 2. Definitions and the restricted two-reservoir setup

## 2.1 Local quantities, units, and the meaning of gas supply

Consider a patch at disc position \(\boldsymbol{x}=(r,\varphi)\), where \(r\) and \(\varphi\) are galactocentric radius and azimuth. All quantities are per projected area. We suppress \(\boldsymbol{x}\) in the intermediate algebra, but initial conditions and coefficients may differ between patches. The main equations have no spatial-transport term: they describe local effective evolution, not a prediction of wind direction or migration between patches.

Let \(\Sigma_{\mathrm{HI}}(t)\) and \(\Sigma_{\mathrm{H_2}}(t)\) be atomic-associated and molecular-associated cold-gas mass surface densities. Both must use the same helium convention. Let \(\Sigma_{\mathrm{SFR}}(t)\) be the rate of formed stellar mass per area, and \(\tau_{\mathrm{dep}}\) the molecular depletion time. We adopt the linear molecular star-formation law

$$
\Sigma_{\mathrm{SFR}}(t)
\equiv\frac{\Sigma_{\mathrm{H_2}}(t)}{\tau_{\mathrm{dep}}}.
\tag{1}
$$

This is a regulator assumption motivated by the approximately linear resolved molecular relation in normal nearby discs ([Bigiel et al. 2011](https://doi.org/10.1088/2041-8205/730/2/L13)); it is not guaranteed in every environmentally disturbed or individual-cloud-scale region. We use one internally consistent area/time unit system in the analytic equations. Numerically, gas columns are in \(M_\odot\,\mathrm{pc}^{-2}\) and times in Gyr. The right-hand side of equation (1) then initially has units \(M_\odot\,\mathrm{pc}^{-2}\,\mathrm{Gyr}^{-1}\). Its numerical value is multiplied by \(10^{-3}\) when printing \(M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}\).

Next, let \(\tau_{\mathrm{conv}}>0\) denote the effective time over which the modeled atomic reservoir transfers mass into the molecular reservoir. We close the transfer rate by

$$
\Sigma_\Phi(t)\equiv
\frac{\Sigma_{\mathrm{HI}}(t)}{\tau_{\mathrm{conv}}},
\qquad
\tau_\Phi(t)\equiv
\frac{\Sigma_{\mathrm{H_2}}(t)}{\Sigma_\Phi(t)}
=\tau_{\mathrm{conv}}
\frac{\Sigma_{\mathrm{H_2}}(t)}{\Sigma_{\mathrm{HI}}(t)}.
\tag{2}
$$

Here \(\Sigma_\Phi\) is a **gas supply-rate surface density**, consistent with the local regulator notation in [Huang et al. (2026), section 4](https://doi.org/10.1093/mnras/stag1019). It has units of mass per area per time. It is not a timescale. The associated molecular replenishment time is \(\tau_\Phi\), which is generally different from \(\tau_{\mathrm{conv}}\). The former divides the *molecular* reservoir by its supply rate, whereas the latter divides the *atomic* reservoir by the same transfer rate. The linear law \(\Sigma_\Phi=\Sigma_{\mathrm{HI}}/\tau_{\mathrm{conv}}\) is an additional closure adopted here; the published regulator framework does not establish that it must hold with constant \(\tau_{\mathrm{conv}}\).

## 2.2 Continuity equations and assumptions

Let \(\gamma_{\mathrm{strip,HI}}\geq0\) be the effective fractional rate of direct atomic-gas removal by ram pressure, in inverse time. We use the same \(\gamma\) notation for either gas phase; the baseline explicitly sets \(\gamma_{\mathrm{strip,H_2}}=0\). Let \(0\leq R<1\) be the effective prompt returned fraction of formed stellar mass, and \(\lambda\geq0\) the feedback mass-loading factor, defined as feedback gas removal divided by the formed-stellar-mass rate. They are dimensionless. The molecular sink associated with stars and feedback is \((1-R+\lambda)\Sigma_{\mathrm{SFR}}\), following the regulator bookkeeping of [Lilly et al. (2013)](https://doi.org/10.1088/0004-637X/772/2/119).

For the main closed-reservoir interval, we assume no external supply into H I, no resolved transport between patches, no explicit H2 dissociation back into H I, and constant \(\tau_{\mathrm{conv}}\), \(\tau_{\mathrm{dep}}\), \(R\), \(\lambda\), and \(\gamma_{\mathrm{strip,HI}}\). Prompt return is incorporated into the effective molecular sink; this is not a phase-resolved recycling calculation. These assumptions give

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{HI}}}{dt}
&=-\Sigma_\Phi-\gamma_{\mathrm{strip,HI}}\Sigma_{\mathrm{HI}},\\
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
&=\Sigma_\Phi-(1-R+\lambda)\Sigma_{\mathrm{SFR}}.
\end{aligned}
\tag{3}
$$

Equation (3) is our two-phase extension of regulator continuity with a declared linear stripping sink. The literature motivates indirect molecular depletion after substantial inner-disc H I loss ([Fumagalli et al. 2009](https://doi.org/10.1088/0004-637X/697/2/1811)); it does not supply this exact two-phase closure. Direct molecular loss can occur in Virgo ([Boselli et al. 2014](https://doi.org/10.1051/0004-6361/201322313)), so setting its coefficient to zero isolates a mechanism rather than establishing its universal absence.

Adding the equations verifies the mass bookkeeping:

$$
\frac{d}{dt}\left(\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H_2}}\right)
=-\gamma_{\mathrm{strip,HI}}\Sigma_{\mathrm{HI}}
-(1-R+\lambda)\Sigma_{\mathrm{SFR}}.
\tag{4}
$$

The internal transfer cancels. Thus an H I molecule-forming supply term is not also an external gain of total cold gas. To keep the solution readable, define physically identified total decay coefficients

$$
\gamma_{\mathrm{HI}}\equiv
\frac{1}{\tau_{\mathrm{conv}}}+\gamma_{\mathrm{strip,HI}},
\qquad
\gamma_{\mathrm{H_2}}\equiv
\frac{1-R+\lambda}{\tau_{\mathrm{dep}}}.
\tag{5}
$$

The first is total disappearance from the modeled atomic reservoir, including conversion and stripping. The second is molecular consumption plus feedback loss per unit molecular mass. **It contains no direct RPS molecular sink.** Both have inverse-time units; neither is a new physical process beyond equation (3).

# 3. Step-by-step solution for H I, H2, and SFR

## 3.1 Solve the atomic reservoir first

Insert equations (1)–(2) into the first balance in equation (3). Because the coefficients are constant on this interval,

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{HI}}}{dt}
&=-\left(\frac{1}{\tau_{\mathrm{conv}}}
+\gamma_{\mathrm{strip,HI}}\right)\Sigma_{\mathrm{HI}}\\
&=-\gamma_{\mathrm{HI}}\Sigma_{\mathrm{HI}}.
\end{aligned}
$$

For positive initial \(\Sigma_{\mathrm{HI},0}\equiv\Sigma_{\mathrm{HI}}(0)\), divide both sides by \(\Sigma_{\mathrm{HI}}\), and integrate between 0 and \(t\):

$$
\begin{aligned}
\int_{\Sigma_{\mathrm{HI},0}}^{\Sigma_{\mathrm{HI}}(t)}
\frac{d\Sigma_{\mathrm{HI}}}{\Sigma_{\mathrm{HI}}}
&=-\gamma_{\mathrm{HI}}\int_0^t du,\\
\ln\!\left[\frac{\Sigma_{\mathrm{HI}}(t)}{\Sigma_{\mathrm{HI},0}}\right]
&=-\gamma_{\mathrm{HI}}t,\\
\Sigma_{\mathrm{HI}}(t)
&=\Sigma_{\mathrm{HI},0}e^{-\gamma_{\mathrm{HI}}t}.
\end{aligned}
\tag{6}
$$

The integration variable \(u\) has time units. If the initial atomic column is zero, the zero solution follows directly without dividing by it. Equation (2) then gives the explicit declining molecular supply:

$$
\Sigma_\Phi(t)=
\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\tag{7}
$$

## 3.2 Substitute this supply into the molecular balance

Using equation (1) in the molecular sink, and then equation (7) in its source, gives

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
&=\Sigma_\Phi
-\frac{1-R+\lambda}{\tau_{\mathrm{dep}}}\Sigma_{\mathrm{H_2}},\\
\frac{d\Sigma_{\mathrm{H_2}}}{dt}
+\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2}}
&=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\end{aligned}
\tag{8}
$$

We now solve this first-order inhomogeneous equation. Its integrating factor is \(e^{\gamma_{\mathrm{H_2}}t}\), since its derivative is \(\gamma_{\mathrm{H_2}}e^{\gamma_{\mathrm{H_2}}t}\). Multiplying **each term** of equation (8) by this factor yields

$$
\begin{aligned}
e^{\gamma_{\mathrm{H_2}}t}\frac{d\Sigma_{\mathrm{H_2}}}{dt}
+\gamma_{\mathrm{H_2}}e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2}}
&=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})t},\\
\frac{d}{dt}\left[e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2}}(t)\right]
&=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})t}.
\end{aligned}
\tag{9}
$$

The second line follows from the product rule, not from discarding the molecular-loss term. Define \(\Sigma_{\mathrm{H_2},0}\equiv\Sigma_{\mathrm{H_2}}(0)\). Integrate equation (9) from the known initial boundary to \(t\):

$$
\begin{aligned}
e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2}}(t)
-\Sigma_{\mathrm{H_2},0}
&=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\int_0^t e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})u}\,du\\
&=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\frac{e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})t}-1}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}},
\end{aligned}
\tag{10}
$$

where the last line presently assumes \(\gamma_{\mathrm{H_2}}\ne\gamma_{\mathrm{HI}}\). Add the initial column to both sides and multiply by \(e^{-\gamma_{\mathrm{H_2}}t}\). In the supplied term, the product of exponentials becomes \(e^{-\gamma_{\mathrm{HI}}t}\). Thus

$$
\boxed{
\begin{aligned}
\Sigma_{\mathrm{H_2}}(t)
&=\Sigma_{\mathrm{H_2},0}e^{-\gamma_{\mathrm{H_2}}t}\\
&\quad+\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}}
\tag{11}
$$

The first term is surviving initial molecular gas. The second is the accumulated transfer from H I, reduced by molecular consumption after arrival. Its positivity is easiest to see before evaluating the integral:

$$
\Sigma_{\mathrm{H_2}}(t)
=\Sigma_{\mathrm{H_2},0}e^{-\gamma_{\mathrm{H_2}}t}
+\int_0^t\Sigma_\Phi(u)
e^{-\gamma_{\mathrm{H_2}}(t-u)}\,du.
\tag{12}
$$

For a nonnegative supply, every element of the integral is nonnegative. At \(t=0\), the integral vanishes, so the solution returns its prescribed initial condition. Equations (6)–(12) are direct solutions of our stated balances.

If the two rates are equal, the integrand in equation (10) equals one. We must integrate that case directly rather than divide by zero:

$$
\Sigma_{\mathrm{H_2}}(t)
=e^{-\gamma_{\mathrm{H_2}}t}
\left(\Sigma_{\mathrm{H_2},0}
+\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}t\right),
\qquad \gamma_{\mathrm{HI}}=\gamma_{\mathrm{H_2}}.
\tag{13}
$$

## 3.3 Write the SFR solution explicitly

Because \(\tau_{\mathrm{dep}}\) is constant in this solution, equation (1) permits us to divide **both terms** of equation (11) by the same depletion time:

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

The initial SFR is \(\Sigma_{\mathrm{SFR},0}=\Sigma_{\mathrm{H_2},0}/\tau_{\mathrm{dep}}\). Dividing equation (14) by it gives a dimensionless expression whose shape depends on the initial supply-to-molecular-column ratio:

$$
\frac{\Sigma_{\mathrm{SFR}}(t)}{\Sigma_{\mathrm{SFR},0}}
=e^{-\gamma_{\mathrm{H_2}}t}
+\frac{\Sigma_{\mathrm{HI},0}}
{\tau_{\mathrm{conv}}\Sigma_{\mathrm{H_2},0}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\tag{15}
$$

For the numerical baseline, we impose **instantaneous molecular balance at the onset**,

$$
\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
=\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0}
=(1-R+\lambda)\Sigma_{\mathrm{SFR},0}.
\tag{16}
$$

This states \(d\Sigma_{\mathrm{H_2}}/dt=0\) at \(t=0\). It does not state that the closed atomic reservoir is in equilibrium. Indeed, equation (6) already says that H I will decline. Substituting equation (16) into equation (15), collecting the coefficients of each exponential, and reversing the sign of numerator and denominator yields

$$
\begin{aligned}
\frac{\Sigma_{\mathrm{SFR}}(t)}{\Sigma_{\mathrm{SFR},0}}
&=e^{-\gamma_{\mathrm{H_2}}t}
+\gamma_{\mathrm{H_2}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}\\
&=\frac{\gamma_{\mathrm{H_2}}e^{-\gamma_{\mathrm{HI}}t}
-\gamma_{\mathrm{HI}}e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}\\
&=\boxed{\frac{\gamma_{\mathrm{HI}}e^{-\gamma_{\mathrm{H_2}}t}
-\gamma_{\mathrm{H_2}}e^{-\gamma_{\mathrm{HI}}t}}
{\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}}}}.
\end{aligned}
\tag{17}
$$

For equal rates, equations (13) and (16) instead give \((1+\gamma_{\mathrm{H_2}}t)e^{-\gamma_{\mathrm{H_2}}t}\). Both forms preserve positive gas and SFR for nonnegative time.

# 4. What kind of SFR decline follows from H I-only stripping?

## 4.1 A delayed onset, followed by a decline

Differentiate the final expression in equation (17) explicitly:

$$
\begin{aligned}
\frac{d}{dt}\!\left(\frac{\Sigma_{\mathrm{SFR}}}{\Sigma_{\mathrm{SFR},0}}\right)
&=\frac{\gamma_{\mathrm{HI}}\gamma_{\mathrm{H_2}}}
{\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}}}
\left(e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}\right),\\
\left.\frac{d\Sigma_{\mathrm{SFR}}}{dt}\right|_0&=0,\\
\left.\frac{d^2\Sigma_{\mathrm{SFR}}}{dt^2}\right|_0
&=-\gamma_{\mathrm{HI}}\gamma_{\mathrm{H_2}}
\Sigma_{\mathrm{SFR},0}.
\end{aligned}
\tag{18}
$$

For \(\gamma_{\mathrm{HI}}>\gamma_{\mathrm{H_2}}>0\), the parenthesis in the first line is negative at \(t>0\), while its prefactor is positive. Reversing the rate ordering reverses both signs, leaving the derivative negative. For equal positive rates, differentiating \((1+\gamma_{\mathrm{H_2}}t)e^{-\gamma_{\mathrm{H_2}}t}\) gives \(-\gamma_{\mathrm{H_2}}^2t e^{-\gamma_{\mathrm{H_2}}t}<0\). Thus the initially balanced molecular reservoir declines monotonically after its zero-slope onset. Its early expansion is

$$
\frac{\Sigma_{\mathrm{SFR}}(t)}{\Sigma_{\mathrm{SFR},0}}
=1-\frac{\gamma_{\mathrm{HI}}\gamma_{\mathrm{H_2}}}{2}t^2
+O(t^3).
\tag{19}
$$

The \(O(t^3)\) notation in equation (19) denotes omitted Taylor terms starting at third order about \(t=0\), with their appropriate inverse-time coefficients; the expression is a short-time expansion, not a late-time approximation.

This delay is the reservoir response: molecular gas remains available while its replenishment begins to fall. The instantaneous effective SFR decline rate follows directly from equation (8):

$$
\begin{aligned}
\Gamma_{\mathrm{SFR}}(t)
&\equiv-\frac{d\ln\Sigma_{\mathrm{SFR}}}{dt}\\
&=\gamma_{\mathrm{H_2}}
-\frac{\Sigma_\Phi(t)}{\Sigma_{\mathrm{H_2}}(t)}
=\gamma_{\mathrm{H_2}}-\frac{1}{\tau_\Phi(t)}.
\end{aligned}
\tag{20}
$$

Hence \(\Gamma_{\mathrm{SFR}}(0)=0\) under equation (16). A single exponential with a constant decline rate cannot represent the whole onset accurately. With different, unbalanced initial conditions, supply can initially exceed consumption and the SFR can initially rise; monotonicity here depends on the stated balance condition.

## 4.2 When is “exponential-like” an appropriate description?

When atomic disappearance is faster than molecular consumption, factor \(e^{-\gamma_{\mathrm{H_2}}t}\) out of equation (17):

$$
\frac{\Sigma_{\mathrm{SFR}}(t)}{\Sigma_{\mathrm{SFR},0}}
=e^{-\gamma_{\mathrm{H_2}}t}
\left[1+\frac{\gamma_{\mathrm{H_2}}}
{\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}}}
\left(1-e^{-(\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}})t}\right)\right].
\tag{21}
$$

After the second exponential becomes small, the bracket approaches the constant \(\gamma_{\mathrm{HI}}/(\gamma_{\mathrm{HI}}-\gamma_{\mathrm{H_2}})\). The late SFR therefore approaches an exponential with rate \(\gamma_{\mathrm{H_2}}\), determined by molecular consumption, rather than by the H I stripping coefficient alone. If molecular consumption is faster, the late decline instead follows the slower atomic-supply coefficient \(\gamma_{\mathrm{HI}}\). Equal rates retain the polynomial prefactor in equation (13). “Delayed decline with a late exponential asymptote” is the precise statement for this model.

## 4.3 Separate ram-pressure removal from closing the supply system

Our baseline also assumes that external H I replenishment has stopped. Even when \(\gamma_{\mathrm{strip,HI}}=0\), the atomic reservoir is consumed by transfer at rate \(1/\tau_{\mathrm{conv}}\), so the closed no-RPS system eventually fades. We therefore compare the stripped history with a **closed, otherwise identical no-RPS control**, obtained by replacing \(\gamma_{\mathrm{HI}}\) by \(1/\tau_{\mathrm{conv}}\) in equation (17). This isolates the additional attenuation associated with the linear stripping term within the same closed-system assumptions.

For example, subtracting the two initial curvatures from equation (18) gives

$$
\left.\frac{d^2}{dt^2}
\left(\Sigma_{\mathrm{SFR}}^{\mathrm{RPS}}
-\Sigma_{\mathrm{SFR}}^{\mathrm{closed,0}}\right)\right|_0
=-\gamma_{\mathrm{strip,HI}}\gamma_{\mathrm{H_2}}
\Sigma_{\mathrm{SFR},0}.
\tag{22}
$$

The superscripts identify the stripped history and closed control. A pre-perturbation system with continuously maintained atomic supply is a different comparison. Appendix A derives that case. We should not attribute the entire decline relative to a maintained control to H I stripping when external replenishment was also removed by assumption.

## 4.4 A quantitative limit on how fast indirect depletion can act

Because the supply term in equation (12) is nonnegative,

$$
\Sigma_{\mathrm{H_2}}(t)\geq
\Sigma_{\mathrm{H_2},0}e^{-\gamma_{\mathrm{H_2}}t},
\qquad
\frac{\Sigma_{\mathrm{SFR}}(t)}{\Sigma_{\mathrm{SFR},0}}
\geq e^{-(1-R+\lambda)t/\tau_{\mathrm{dep}}}.
\tag{23}
$$

This bound holds even without initial balance, provided the depletion time and molecular sink coefficient stay constant and the supply remains nonnegative. It says that removing H I cannot make molecular gas disappear faster than immediate, complete molecular-supply interruption under the same consumption law. For an attenuation target \(f_{\mathrm{SFR}}\in(0,1)\), a necessary condition is

$$
t\geq\frac{-\ln f_{\mathrm{SFR}}}{\gamma_{\mathrm{H_2}}},
\qquad
\tau_{\mathrm{dep}}\leq
\frac{(1-R+\lambda)t}{-\ln f_{\mathrm{SFR}}}.
\tag{24}
$$

With \(R=0.4\), \(\lambda=0\), \(\tau_{\mathrm{dep}}=2\) Gyr, the one-Gyr lower bound is 0.741. Reaching 0.407 requires at least 3.00 Gyr under this bound, or \(\tau_{\mathrm{dep}}\leq0.667\) Gyr for an assumed one-Gyr interval. These are necessary conditions, not fitted infall times or measured MAUVE depletion times. Residual supply makes attainment slower. Feedback, phase exchange, a changed depletion time, or direct molecular removal would change the bound. A depletion time that grows suppresses instantaneous SFR through equation (1) even without rapidly destroying molecular gas; that is a separate efficiency mechanism.

# 5. From the SFR history to H II and leaked-OB DIG luminosities

## 5.1 The common ionizing source

Let \(q_{\mathrm{H,y}}(a)\) be hydrogen-ionizing photons produced per second per unit formed stellar mass by a young population of age \(a\). The production rate per projected area is the sum of stellar cohorts:

$$
Q_{\mathrm{H,y}}(t)
=\int_0^\infty q_{\mathrm{H,y}}(a)
\Sigma_{\mathrm{SFR}}(t-a)\,da.
\tag{25}
$$

Use a common time unit for \(a\), \(t\), and the mass-formation rate in this integral. For the printed \(M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}\) SFR, integration over ages in years yields photons \(\mathrm{s}^{-1}\,\mathrm{kpc}^{-2}\). The population-synthesis response depends on stellar evolution, metallicity, and the initial mass function (IMF); the short-age sensitivity of Hα is discussed by [Kennicutt & Evans (2012)](https://doi.org/10.1146/annurev-astro-081811-125610). A declining SFR and a declining line luminosity are therefore linked through a response function, not an assumption of instantaneous tracking after any abrupt change.

For ionization-bounded gas in photoionization/recombination balance, define \(\alpha_{\mathrm B}\) as the Case-B recombination coefficient and \(\alpha^{\mathrm{eff}}_\alpha\) as the effective Hα recombination coefficient. If \(Q_{\mathrm{abs}}^j\) is the rate of photons absorbed by hydrogen in component \(j\),

$$
\begin{aligned}
Q_{\mathrm{abs}}^j&=\int_j\alpha_{\mathrm B}n_en_p\,dV,\\
L_\alpha^j&=h\nu_\alpha\int_j\alpha^{\mathrm{eff}}_\alpha n_en_p\,dV,\\
L_\alpha^j&=h\nu_\alpha
\frac{\alpha^{\mathrm{eff}}_\alpha}{\alpha_{\mathrm B}}
Q_{\mathrm{abs}}^j.
\end{aligned}
\tag{26}
$$

Here \(n_e,n_p\) are electron and proton densities; \(dV\) is volume; \(h\nu_\alpha\) is the Hα photon energy; \(L_\alpha^j\) is luminosity, before division by projected area. The last line assumes the coefficient ratio is approximately uniform within the component. This is established recombination physics ([Hummer & Storey 1987](https://doi.org/10.1093/mnras/224.3.801); [Storey & Hummer 1995](https://doi.org/10.1093/mnras/272.1.41)). We denote the corresponding full-hydrogen-absorption young-source surface luminosity by \(\mathcal L_\alpha^{\mathrm{young,max}}\). It is a photon-budget normalization under specified nebular conditions, not a second observed gas component.

## 5.2 Mutually exclusive absorption fractions

Define \(f_{\mathrm{HII}}\) as the fraction of young ionizing photons absorbed in compact H II gas, and \(f_{\mathrm{leak}}\) as the fraction that escapes the compact regions **and is absorbed in diffuse gas**. The remaining \(f_{\mathrm{unabs}}\) accounts for photons that do not ionize either modeled gas component, including escape from the modeled domain or absorption by dust before hydrogen ionization. Then

$$
f_{\mathrm{HII}}+f_{\mathrm{leak}}+f_{\mathrm{unabs}}=1,
\qquad 0\leq f_j\leq1.
\tag{27}
$$

Leakage out of an H II region is not automatically absorption in DIG. These fractions prevent counting the same photon in both components. Under the baseline's equal Hα recombination yield,

$$
\begin{aligned}
\mathcal L_\alpha^{\mathrm{HII}}(t)
&=f_{\mathrm{HII}}(t)\mathcal L_\alpha^{\mathrm{young,max}}(t),\\
\mathcal L_\alpha^{\mathrm{non-HII}}(t)
&=f_{\mathrm{leak}}(t)\mathcal L_\alpha^{\mathrm{young,max}}(t),\\
\mathcal L_\alpha(t)
&=f_{\mathrm{cap}}(t)\mathcal L_\alpha^{\mathrm{young,max}}(t),\quad
f_{\mathrm{cap}}\equiv f_{\mathrm{HII}}+f_{\mathrm{leak}}.
\end{aligned}
\tag{28}
$$

The non-H II component is leaked-OB-photon DIG in this restricted calculation. Treating it as dominant in Hα is motivated by [Belfiore et al. (2022)](https://doi.org/10.1051/0004-6361/202141859). It does not establish that OB leakage explains every NSF forbidden line: a small harder-source contribution can matter disproportionately to [O III]. We retain that distinction in the limitations rather than adding extra emitting components to the main algebra.

For the numerical calculation, we write \(\mathcal L_\alpha^{\mathrm{young,max}}=\Sigma_{\mathrm{SFR}}^{\mathrm{ion}}/C_\alpha\), where \(C_\alpha\) is the unchanged MAUVE Hα calibration and \(\Sigma_{\mathrm{SFR}}^{\mathrm{ion}}\) denotes the SFR filtered by a normalized young-ionizing response. The particular numerical response is

$$
\begin{aligned}
K_\alpha(a)&=\frac{1}{\tau_{\mathrm{ion}}}e^{-a/\tau_{\mathrm{ion}}},
\quad a\geq0,\\
\Sigma_{\mathrm{SFR}}^{\mathrm{ion}}(t)
&=\int_0^\infty K_\alpha(a)\Sigma_{\mathrm{SFR}}(t-a)\,da,\\
\tau_{\mathrm{ion}}\frac{d\Sigma_{\mathrm{SFR}}^{\mathrm{ion}}}{dt}
&=\Sigma_{\mathrm{SFR}}-\Sigma_{\mathrm{SFR}}^{\mathrm{ion}}.
\end{aligned}
\tag{29}
$$

This exponential kernel with \(\tau_{\mathrm{ion}}=3\) Myr is an explicitly illustrative smoothing function, not a population-synthesis prediction. The pre-interval SFR is taken constant to initialize \(\Sigma_{\mathrm{SFR}}^{\mathrm{ion}}(0)=\Sigma_{\mathrm{SFR},0}\). Differentiating the convolution in its causal form gives the last line. Its normalization is unity, so it returns a constant SFR unchanged. The young-source kernel is shared by both components; a separate long-lived DIG source is not created by giving it the name non-H II.

## 5.3 Why OB fading alone leaves the weight unchanged

The non-H II Hα luminosity weight is

$$
\begin{aligned}
w_{\mathrm{non-HII},\alpha}(t)
&\equiv\frac{\mathcal L_\alpha^{\mathrm{non-HII}}}
{\mathcal L_\alpha^{\mathrm{HII}}+\mathcal L_\alpha^{\mathrm{non-HII}}}\\
&=\frac{f_{\mathrm{leak}}\mathcal L_\alpha^{\mathrm{young,max}}}
{(f_{\mathrm{HII}}+f_{\mathrm{leak}})\mathcal L_\alpha^{\mathrm{young,max}}}\\
&=\frac{f_{\mathrm{leak}}}{f_{\mathrm{HII}}+f_{\mathrm{leak}}}.
\end{aligned}
\tag{30}
$$

This cancellation applies while total emission is positive. Consequently, if the fractions stay fixed, the declining common photon source makes both components fainter in the same proportion. It cannot increase their relative weight. The general logarithmic form makes the additional requirement particularly clear:

$$
\frac{d}{dt}\ln\!\left(\frac{w_{\mathrm{non-HII},\alpha}}
{1-w_{\mathrm{non-HII},\alpha}}\right)
=\frac{d\ln f_{\mathrm{leak}}}{dt}
-\frac{d\ln f_{\mathrm{HII}}}{dt}.
\tag{31}
$$

An increasing weight requires preferential loss of compact absorption relative to diffuse absorption, a changed geometry/covering factor, transport of photons from neighboring regions, or another explicitly modeled effect. H I loss can plausibly alter these conditions, but equations (3)–(17) do not establish a quantitative relation between the reservoir columns and these fractions. The main model therefore supplies the SFR history and the mixing algebra, while the photon-allocation history remains a separate constraint to measure or model.

# 6. Derive the observed line ratios before interpreting them

## 6.1 Linear luminosity mixing and normalized weights

Let \(\mathcal L_\ell^j\) denote the emitted, ideally attenuation-corrected surface luminosity of line \(\ell\) in component \(j\in\{\mathrm{HII},\mathrm{non-HII}\}\). Define \(R_{\mathrm{NII}}^j=\mathcal L_{[\mathrm{NII}]6583}^j/\mathcal L_\alpha^j\), with analogous \(R_{\mathrm{SII}}^j\) for the summed [S II] doublet and \(R_{\mathrm{OIII}}^j=\mathcal L_{[\mathrm{OIII}]5007}^j/\mathcal L_\beta^j\). These ratios describe a physical component under specified conditions, not the mean ratio of the observed NSF or SF class.

Fluxes add for unresolved components at the same distance. For [N II], first sum the numerator luminosities, substitute the component ratios, and only then divide by the summed Hα luminosity:

$$
\begin{aligned}
R_{\mathrm{NII}}^{\mathrm{obs}}
&=\frac{\mathcal L_{[\mathrm{NII}]}^{\mathrm{HII}}
+\mathcal L_{[\mathrm{NII}]}^{\mathrm{non-HII}}}
{\mathcal L_\alpha^{\mathrm{HII}}+\mathcal L_\alpha^{\mathrm{non-HII}}}\\
&=\frac{R_{\mathrm{NII}}^{\mathrm{HII}}\mathcal L_\alpha^{\mathrm{HII}}
+R_{\mathrm{NII}}^{\mathrm{non-HII}}\mathcal L_\alpha^{\mathrm{non-HII}}}
{\mathcal L_\alpha^{\mathrm{HII}}+\mathcal L_\alpha^{\mathrm{non-HII}}}\\
&=\boxed{(1-w_{\mathrm{non-HII},\alpha})R_{\mathrm{NII}}^{\mathrm{HII}}
+w_{\mathrm{non-HII},\alpha}R_{\mathrm{NII}}^{\mathrm{non-HII}}}.
\end{aligned}
\tag{32}
$$

The H II weight is \(w_{\mathrm{HII},\alpha}=1-w_{\mathrm{non-HII},\alpha}\); their sum is one because both divide by the same positive total Hα luminosity. They are **Balmer-luminosity fractions**, not gas-mass, area, or numbers-of-region fractions. Equation (32) is exact under additive component luminosities regardless of whether its ratios or weights evolve. An explicit earlier [S II]/Hα mixture appears in [Blanc et al. (2009), equations 7–8](https://doi.org/10.1088/0004-637X/704/1/842), with empirical component spectra and a metallicity scaling. Our derivation here states the normalization directly.

For [S II], the same denominator gives

$$
R_{\mathrm{SII}}^{\mathrm{obs}}
=(1-w_{\mathrm{non-HII},\alpha})R_{\mathrm{SII}}^{\mathrm{HII}}
+w_{\mathrm{non-HII},\alpha}R_{\mathrm{SII}}^{\mathrm{non-HII}}.
\tag{33}
$$

For [O III]/Hβ, the relevant weight instead uses Hβ. Define the component Balmer decrements \(B_j=\mathcal L_\alpha^j/\mathcal L_\beta^j\). Substituting \(\mathcal L_\beta^j=\mathcal L_\alpha^j/B_j\) gives

$$
\begin{aligned}
w_{\mathrm{non-HII},\beta}
&=\frac{w_{\mathrm{non-HII},\alpha}/B_{\mathrm{non-HII}}}
{(1-w_{\mathrm{non-HII},\alpha})/B_{\mathrm{HII}}
+w_{\mathrm{non-HII},\alpha}/B_{\mathrm{non-HII}}},\\
R_{\mathrm{OIII}}^{\mathrm{obs}}
&=(1-w_{\mathrm{non-HII},\beta})R_{\mathrm{OIII}}^{\mathrm{HII}}
+w_{\mathrm{non-HII},\beta}R_{\mathrm{OIII}}^{\mathrm{non-HII}}.
\end{aligned}
\tag{34}
$$

Only when \(B_{\mathrm{HII}}=B_{\mathrm{non-HII}}\) do the two weights coincide. We set both to 2.86 in the numerical baseline, a low-density, approximately \(10^4\) K Case-B example ([Storey & Hummer 1995](https://doi.org/10.1093/mnras/272.1.41)). This value is condition dependent. A pipeline that imposes one intrinsic decrement can make corrected total Hα/Hβ very close to 2.86; that is not independent proof that both hidden components have that decrement. Appendix B treats attenuation and unequal decrements.

## 6.2 What exactly do constant “intrinsic” component ratios assume?

For the main analytical exercise, hold each \(R_k^{\mathrm{HII}}\) and \(R_k^{\mathrm{non-HII}}\) constant **within a restricted comparison of similar local conditions**. The subscript \(k\) indexes the three diagnostic ratios. The spectra can still differ between matched stellar-density/radius/abundance bins. This assumption does not make either endpoint a universal number for every H II region or every DIG region.

In general, each ratio depends on gas abundance \(Z\), ionization parameter \(U\), electron density \(n_e\), electron temperature \(T_e\), ionizing spectral shape, and abundance ratios such as N/O. Declining OB luminosity can change \(U\); selective transmission can change the spectrum reaching DIG. Hence constant spectra are a deliberately restricted empirical hypothesis, not a necessary consequence of gas loss. DIG studies demonstrate such spectral variation and the limitations of solely leaked-H II explanations ([Zhang et al. 2017](https://doi.org/10.1093/mnras/stw3308)). We begin with fixed spectra because they yield a transparent and falsifiable prediction; we permit their evolution only after examining those restrictions and the multi-line tests below.

## 6.3 Increasing ratios do not require brighter forbidden lines

For a ratio with Hα denominator and fixed endpoints, write \(w\equiv w_{\mathrm{non-HII},\alpha}\) within this subsection. Then \(\mathcal L_\ell=R_k^{\mathrm{obs}}\mathcal L_\alpha\), so the product rule gives

$$
\begin{aligned}
\frac{d\ln\mathcal L_\ell}{dt}
&=\frac{d\ln\mathcal L_\alpha}{dt}
+\frac{1}{R_k^{\mathrm{obs}}}\frac{dR_k^{\mathrm{obs}}}{dt},\\
\frac{dR_k^{\mathrm{obs}}}{dt}
&=(R_k^{\mathrm{non-HII}}-R_k^{\mathrm{HII}})\frac{dw}{dt}.
\end{aligned}
\tag{35}
$$

If the non-H II endpoint is higher and its weight grows, the ratio rises. The forbidden line still fades if the fractional Hα decline exceeds the fractional ratio increase. For two comparison endpoints 0 and 1, this is the finite condition

$$
\frac{\mathcal L_{\ell,1}}{\mathcal L_{\ell,0}}
=\frac{\mathcal L_{\alpha,1}}{\mathcal L_{\alpha,0}}
\frac{R_{k,1}^{\mathrm{obs}}}{R_{k,0}^{\mathrm{obs}}}<1.
\tag{36}
$$

An increasing non-H II weight likewise does not imply increasing non-H II luminosity:

$$
\begin{aligned}
\frac{\mathcal L_{\alpha,1}^{\mathrm{non-HII}}}
{\mathcal L_{\alpha,0}^{\mathrm{non-HII}}}
&=\frac{\mathcal L_{\alpha,1}}{\mathcal L_{\alpha,0}}\frac{w_1}{w_0},\\
\frac{\mathcal L_{\alpha,1}^{\mathrm{HII}}}
{\mathcal L_{\alpha,0}^{\mathrm{HII}}}
&=\frac{\mathcal L_{\alpha,1}}{\mathcal L_{\alpha,0}}
\frac{1-w_1}{1-w_0}.
\end{aligned}
\tag{37}
$$

Both fade only if these factors are below one. This provides a quantitative check on any proposed photon-allocation closure; it cannot simply be assumed from a rising DIG fraction.

## 6.4 An analytical differential-fading parameterization

For a line-level illustration, suppose the two Hα luminosities obey constant positive effective fading rates over a finite comparison interval:

$$
\mathcal L_\alpha^{\mathrm{HII}}(t)
=\mathcal L_{\alpha,0}^{\mathrm{HII}}e^{-\gamma_{\alpha,\mathrm{HII}}t},
\qquad
\mathcal L_\alpha^{\mathrm{non-HII}}(t)
=\mathcal L_{\alpha,0}^{\mathrm{non-HII}}e^{-\gamma_{\alpha,\mathrm{non-HII}}t}.
\tag{38}
$$

These \(\gamma_{\alpha,j}\) are **effective luminosity decline rates**, distinct from gas-loss coefficients. Equation (38) is a phenomenological assumption of this report, not a literature-established law that DIG must fade more slowly. It is not the same solution as equation (17). Let \(w_0\) be the initial Hα non-H II weight. Dividing the two component luminosities and using the definition of \(w\) yields

$$
\begin{aligned}
\frac{1-w(t)}{w(t)}
&=\frac{1-w_0}{w_0}
e^{-(\gamma_{\alpha,\mathrm{HII}}-\gamma_{\alpha,\mathrm{non-HII}})t},\\
\boxed{w(t)}
&=\boxed{\left[1+\frac{1-w_0}{w_0}
e^{-(\gamma_{\alpha,\mathrm{HII}}-\gamma_{\alpha,\mathrm{non-HII}})t}\right]^{-1}},\\
\frac{dw}{dt}
&=(\gamma_{\alpha,\mathrm{HII}}-\gamma_{\alpha,\mathrm{non-HII}})w(1-w).
\end{aligned}
\tag{39}
$$

The final line follows by differentiating the inverse bracket, then replacing that bracket using the second line. Positive relative fading, \(\gamma_{\alpha,\mathrm{HII}}>\gamma_{\alpha,\mathrm{non-HII}}\), increases the weight even when both absolute luminosities decline. Substitution of equation (39) into equations (32)–(34) gives explicit analytical ratio histories under fixed spectra.

For an OB-only DIG interpretation, equation (38) must still satisfy equation (28). Over the modeled interval, compute \(f_{\mathrm{HII}}=\mathcal L_\alpha^{\mathrm{HII}}/\mathcal L_\alpha^{\mathrm{young,max}}\) and \(f_{\mathrm{leak}}=\mathcal L_\alpha^{\mathrm{non-HII}}/\mathcal L_\alpha^{\mathrm{young,max}}\), require a sum no larger than one, and reject negative fractions. A slowly fading DIG term cannot persist indefinitely after its OB photon supply disappears. Neither gas retention nor the algebraic mixture generates the missing photons. Consequently, the numerical line interpolation below is reported separately from the H I-only gas prediction; they are not presented as one fitted causal model.

# 7. Quantify changing weights and changing component spectra

## 7.1 Clarify the two proposed interpretations

There are two meanings of “the non-H II ratio becomes more prominent.” If it means that an unchanged non-H II spectrum contributes a larger fraction of the observed luminosity, it is exactly the increasing-weight interpretation in equation (32). It is not a second mechanism. If it means that \(R_k^{\mathrm{non-HII}}\) itself changes, then it is intrinsic spectral evolution. The H II spectrum can remain fixed in either case.

Allow all terms in equation (32) to depend on time, and use \(w=w_{\mathrm{non-HII},\alpha}\) for an Hα-based ratio. Applying the product rule to each term gives

$$
\begin{aligned}
\frac{dR_k^{\mathrm{obs}}}{dt}
&=-R_k^{\mathrm{HII}}\frac{dw}{dt}
+(1-w)\frac{dR_k^{\mathrm{HII}}}{dt}\\
&\quad+R_k^{\mathrm{non-HII}}\frac{dw}{dt}
+w\frac{dR_k^{\mathrm{non-HII}}}{dt}\\
&=\boxed{(R_k^{\mathrm{non-HII}}-R_k^{\mathrm{HII}})\frac{dw}{dt}
+(1-w)\frac{dR_k^{\mathrm{HII}}}{dt}
+w\frac{dR_k^{\mathrm{non-HII}}}{dt}}.
\end{aligned}
\tag{40}
$$

The first term is a changing luminosity weight; the second and third are intrinsic H II and non-H II spectral evolution. For [O III]/Hβ, use \(w_{\mathrm{non-HII},\beta}\) throughout. Equation (40) is our algebraic decomposition, not a claim that the data have separately measured its terms. A single ratio cannot identify three unknown contributions. Even with a fixed H II spectrum, a changing weight and a changing DIG spectrum remain degenerate.

Because our observed stages are snapshots, an exact finite-difference version is more appropriate than treating a stage difference as \(d/dt\). Let \(\Delta R=R_1-R_0\), \(\Delta w=w_1-w_0\), and use bars for arithmetic midpoint values, for example \(\bar w=(w_1+w_0)/2\). The identity \(\Delta(ab)=\bar a\Delta b+\bar b\Delta a\) follows by expanding both midpoint terms; the crossed products cancel. Applying it to both terms of the mixture yields

$$
\boxed{
\Delta R_k^{\mathrm{obs}}
=(\bar R_k^{\mathrm{non-HII}}-\bar R_k^{\mathrm{HII}})\Delta w
+(1-\bar w)\Delta R_k^{\mathrm{HII}}
+\bar w\Delta R_k^{\mathrm{non-HII}}.}
\tag{41}
$$

There is no omitted remainder in this identity. Its decomposition becomes measurable only after obtaining component spectra or weights from additional information. Under both fixed spectra, the last two terms vanish. Under fixed \(w\) and fixed H II spectrum, instead \(\Delta R_k^{\mathrm{non-HII}}=\Delta R_k^{\mathrm{obs}}/w\) for \(w>0\). The same observed change can therefore be assigned to different terms unless independent constraints are available.

## 7.2 Test whether the H II endpoint can be held fixed

The SF-selected sample is useful, but it is not by itself a measurement of a pure H II endpoint. It can contain diffuse emission, and its BPT selection uses some of the very ratios being tested. A stable SF-class mean may result from changes in contamination, surviving support, or abundance distribution. Conversely, a changing mean need not establish temporal evolution of each H II region.

A useful test is to construct an **H II-dominated comparison** using information as independent of the tested BPT ratio as feasible: compact Hα morphology with PSF-aware segmentation, high Hα surface brightness, high EW(Hα), reliable Balmer measurements and stellar-continuum subtraction, and narrow intrinsic line widths. The ordinary MAUVE SF selection, including its finite H II flag, EW(Hα) greater than 6 Å, and intrinsic Hα dispersion below 45 km s\(^{-1}\), remains the parent observational selection. Additional purity cuts should be explored as a grid rather than given a newly invented universal threshold. Brightness, EW, and narrow lines increase plausibility; they do not prove absence of DIG or shocks. Morphological separation has a practical precedent in [Belfiore et al. (2022)](https://doi.org/10.1051/0004-6361/202141859).

Next, compare stages on common support in \(\log_{10}\Sigma_*\), \(r/R_e\), and galaxy mass, where \(R_e\) is effective radius. Control abundance differences with independently constrained quantities when possible. Avoid estimating metallicity from the same N2/O3 lines and then treating it as an independent control for their evolution. Matching \(U\) can test a conditional spectrum at fixed nebular state, but it can also remove the very environmental response we want to diagnose; total and conditional comparisons should be reported separately.

For each diagnostic \(k\), a pre-peak reference surface \(F_k\) can describe the expected **logarithmic** ratio in the H II-dominated subset at these control quantities. Define

$$
\Delta_k^{\mathrm{HII-dom}}
\equiv\log_{10}R_k^{\mathrm{obs}}
-F_k(\log_{10}\Sigma_*,r/R_e,M_*),
\tag{42}
$$

with \(M_*\) the total stellar mass and consistent units implicit in each fitted coordinate. This is a proposed statistical diagnostic, not a fit executed here. Retain equal-galaxy weighting and resample whole galaxies; pixels are not independent environmental realizations. Repeat the comparison as the morphology/EW/brightness purity requirement strengthens. Agreement within a predeclared physically meaningful equivalence interval would support a constant H II endpoint over that range. A non-significant trend alone does not establish constancy. Residual DIG contamination and the selection's use of BPT must enter the interpretation.

## 7.3 Test a single weight against several line ratios

Suppose the two component spectra are independently calibrated in a restricted local-condition bin. Under equal decrements, define the three-element vector of **linear** ratios

$$
\boldsymbol R^{\mathrm{obs}}
=\boldsymbol R^{\mathrm{HII}}
+w\left(\boldsymbol R^{\mathrm{non-HII}}-\boldsymbol R^{\mathrm{HII}}\right),
\qquad 0\leq w\leq1.
\tag{43}
$$

All mixtures then lie on the line segment between the two spectra in linear-ratio space. The usual logarithmic BPT projection is a curved image of this segment, not a linear weighted sum of logarithms. For each diagnostic with nonzero component contrast, inversion gives

$$
\widehat w_k
=\frac{R_k^{\mathrm{obs}}-R_k^{\mathrm{HII}}}
{R_k^{\mathrm{non-HII}}-R_k^{\mathrm{HII}}}.
\tag{44}
$$

Three necessary checks follow. Each inferred weight must lie within [0,1], the several estimates must be consistent after propagating covariance and endpoint uncertainty, and independently separated component spectra must not show incompatible evolution. With unequal decrements, equation (44) applied to O3 yields \(w_\beta\); convert it to the corresponding Hα weight using equation (34) before comparing it to the N2/S2 estimates. Negligible endpoint contrast makes an inferred weight unstable, regardless of measurement precision.

For a covariance-weighted joint estimate, let \(\boldsymbol d=\boldsymbol R^{\mathrm{non-HII}}-\boldsymbol R^{\mathrm{HII}}\) and let \(\boldsymbol C\) be the covariance matrix of measured linear ratios. Minimize

$$
\chi^2(w)=
\left(\boldsymbol R^{\mathrm{obs}}-\boldsymbol R^{\mathrm{HII}}-w\boldsymbol d\right)^T
\boldsymbol C^{-1}
\left(\boldsymbol R^{\mathrm{obs}}-\boldsymbol R^{\mathrm{HII}}-w\boldsymbol d\right).
\tag{45}
$$

Differentiating the quadratic with respect to \(w\) gives

$$
\begin{aligned}
\frac{d\chi^2}{dw}
&=-2\boldsymbol d^T\boldsymbol C^{-1}
\left(\boldsymbol R^{\mathrm{obs}}-\boldsymbol R^{\mathrm{HII}}-w\boldsymbol d\right),\\
\widehat w_{\mathrm{unconstrained}}
&=\frac{\boldsymbol d^T\boldsymbol C^{-1}
(\boldsymbol R^{\mathrm{obs}}-\boldsymbol R^{\mathrm{HII}})}
{\boldsymbol d^T\boldsymbol C^{-1}\boldsymbol d}.
\end{aligned}
\tag{46}
$$

For fixed spectra and a positive-definite covariance matrix, constrain this quadratic solution to [0,1]. An optimum on a boundary does not remove evidence of a poor spectral fit. If the covariance is Gaussian and the spectra are known exactly, the interior formal variance is \((\boldsymbol d^T\boldsymbol C^{-1}\boldsymbol d)^{-1}\). Those assumptions are too strong for the Table 1 summary alone. We therefore do **not** assign a chi-square significance or a formal source-fraction uncertainty to its numerical inversions.

A preferable observational implementation fits the five line fluxes directly with two nonnegative Balmer amplitudes, consistently modeling component decrements and attenuation (Appendix C). Separate data for calibration and testing, ideally by galaxy. A weight derived from [S II]/Hα cannot independently validate [S II]/Hα mixing; it can predict other ratios under the specified spectra. Correlations with Hα brightness or EW should also acknowledge shared denominators and selection. Inconsistent weights reject the particular fixed two-spectrum description under its measurement assumptions; they do not uniquely identify evolving DIG spectra instead of abundance variation, extra ionization sources, PSF mixing, or dust errors.

# 8. Numerical predictions and their relation to MAUVE

## 8.1 Inputs and what is actually calibrated

The executable calculation accompanies this report in `assets/20260930_HI_stripping_line_mixing/model_predictions.py`. The gas calculation is **not a fit to measured H I or CO**, because those local gas columns are not supplied here. Its initial SFR scale comes from Table 1 pre-peak SF Hα, after assuming a captured young-photon fraction of 0.95. The 2-Gyr depletion time is a rounded normal-disc example, rather than a measured value for this high-stellar-density MAUVE bin. For comparison, [Bigiel et al. (2011)](https://doi.org/10.1088/2041-8205/730/2/L13) find a roughly 2.35-Gyr median under their conventions and resolution. The remaining inputs are declared conditions.

**Table 2.** Baseline gas and photon inputs. Numerical \(\Sigma_{\mathrm{SFR}}\) uses \(M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}\); gas columns use \(M_\odot\,\mathrm{pc}^{-2}\).

| Quantity | Value | Origin or role |
|:--|--:|:--|
| \(C_\alpha\) | \(4.9835821\times10^{-42}\) | Unchanged MAUVE `SFR+Z.py`; SFR per Hα luminosity |
| \(f_{\mathrm{cap}}\) | 0.95 | Assumed; not a measured escape fraction |
| \(\Sigma_{\mathrm{SFR},0}\) | 0.0819101 | \(C_\alpha\mathcal L_{\alpha,\mathrm{pre,SF}}/f_{\mathrm{cap}}\) |
| \(\tau_{\mathrm{dep}}\) | 2 Gyr | Constant normal-disc illustration |
| \(\Sigma_{\mathrm{H_2},0}\) | 163.820 | \(1000\tau_{\mathrm{dep}}\Sigma_{\mathrm{SFR},0}\); inferred under assumed law, not CO |
| \(\Sigma_{\mathrm{HI},0}\) | 10 | Illustrative column, not a MAUVE H I measurement |
| \(R,\lambda\) | 0.4, 0 | Effective return and no feedback loss in this example |
| \(\tau_{\mathrm{conv}}\) | 0.203475 Gyr | Imposed initial molecular balance, equation (16) |
| \(\gamma_{\mathrm{strip,HI}}\) | 3 Gyr\(^{-1}\) | Assumed effective interval coefficient |
| \(\gamma_{\mathrm{strip,H_2}}\) | 0 | Baseline restriction |
| \(\gamma_{\mathrm{HI}},\gamma_{\mathrm{H_2}}\) | 7.91460, 0.30000 Gyr\(^{-1}\) | Composite rates from equation (5) |
| \(\tau_{\mathrm{ion}}\) | 3 Myr | Illustrative shared response kernel |
| \(w_{\mathrm{non-HII},\alpha}\) | 0.10 | Fixed-photon-partition null case |
| \(B_{\mathrm{HII}},B_{\mathrm{non-HII}}\) | 2.86, 2.86 | Equal Case-B baseline |

For the null case, \(f_{\mathrm{HII}}=0.855\), \(f_{\mathrm{leak}}=0.095\), and \(f_{\mathrm{unabs}}=0.05\). The combined molecular/helium convention must be made compatible with any future CO calibration; currently the gas normalization follows equation (1) and the adopted Hα calibration. The inferred molecular column is above much of the normal-disc range used to motivate the constant depletion time, reinforcing that it is an illustrative normalization.

**Table 3.** Fixed spectra used to challenge the simple mixture. These rounded values are guided by the observed range and retained from the previous report. They are not independently isolated H II/DIG spectra.

| Ratio | \(R_k^{\mathrm{HII}}\) | \(R_k^{\mathrm{non-HII}}\) |
|:--|--:|--:|
| [N II]6583/Hα | 0.22 | 0.60 |
| [S II] 6716+6731/Hα | 0.20 | 0.45 |
| [O III]5007/Hβ | 0.57 | 0.77 |

The single non-H II template is a two-component mathematical hypothesis. It should not be described as a demonstrated OB-only photoionization spectrum. Its physical realizability must be checked against appropriate photoionization conditions before interpreting fitted weights as OB-powered DIG fractions.

## 8.2 Prediction 1: indirect molecular depletion is real but modest here

At one Gyr, the H I-only solution gives

$$
\begin{aligned}
\Sigma_{\mathrm{HI}}/\Sigma_{\mathrm{HI},0}&\simeq0.000365,\\
\Sigma_{\mathrm{H_2}}/\Sigma_{\mathrm{H_2},0}
=\Sigma_{\mathrm{SFR}}/\Sigma_{\mathrm{SFR},0}&=0.769991,\\
\log_{10}(\Sigma_{\mathrm{SFR}}/\Sigma_{\mathrm{SFR},0})&=-0.113515.
\end{aligned}
\tag{47}
$$

The almost exhausted atomic reservoir does not imply an almost exhausted molecular reservoir. Most of the retained molecular gas survives on its 3.33-Gyr net-consumption timescale \(1/\gamma_{\mathrm{H_2}}\). This reproduces the *direction* of surviving-region SFR suppression, but does not reproduce an attenuation of 0.407 in an assumed one-Gyr interval. The earlier very strong one-Gyr decline depended substantially on nonzero direct molecular removal; we cannot retain that numerical conclusion after setting this term to zero.

The closed no-RPS control gives 0.788502 at the same time. The stripped/control ratio is 0.976523, a further reduction of only about 2.35% within this particular closure. That small difference arises because the assumed atomic conversion time is already very short: the closed control also rapidly loses its molecular supply by exhausting H I. It is not evidence that real RPS has a negligible effect. It identifies a limitation of this initial gas ratio and externally closed comparison. External H I replenishment, different local gas columns, and measured conversion histories need to be constrained before attributing a numerical suppression amplitude to RPS.

![Figure 1. H I-only reservoir response and depletion-time sensitivity. Left: atomic removal is rapid but molecular gas and SFR decline slowly; the closed no-RPS control also fades. Right: sensitivity to constant depletion time, retaining the same initial gas columns and restoring initial molecular balance separately in each curve. Changing depletion time here changes the initial SFR and conversion time; these are normalized sensitivity experiments, not alternative MAUVE fits. The horizontal 0.407 line is a descriptive observational scale, not an assigned stage age.](assets/20260930_HI_stripping_line_mixing/figure01_HI_only_SFR.png)

The faster-depletion curves demonstrate how consumption could produce stronger suppression without direct molecular stripping, but they do not establish that MAUVE actually has those shorter depletion times. Resolved Virgo evidence also allows changing molecular content and efficiency ([Brown et al. 2023](https://arxiv.org/abs/2308.10943)). A model that fixes efficiency isolates only part of that behavior.

## 8.3 Prediction 2: fixed photon fractions fade the lines but preserve all ratios

Combining the same gas history with equation (29) and constant fractions gives an Hα decline factor of 0.770684 at one Gyr; the slight difference from 0.769991 is the short response delay. H II and non-H II Hα fade by exactly the same fractional amount. From Table 3 and \(w=0.10\),

$$
R_{\mathrm{NII}}^{\mathrm{obs}}=0.258,
\qquad R_{\mathrm{SII}}^{\mathrm{obs}}=0.225,
\qquad R_{\mathrm{OIII}}^{\mathrm{obs}}=0.590.
\tag{48}
$$

These ratios remain constant at all modeled times. Their initial scales are close to the pre-peak SF bin, partly because the templates were guided by the observed range. The unchanged ratios are the actual null prediction: **H I-only supply loss plus common OB fading, with fixed photon allocation and fixed spectra, does not itself produce the NSF ratio trend.** To generate a trend, the model needs changed allocation, changed component spectra, or a more complete source/transport description. No parameter in the gas continuity solution uniquely selects those alternatives.

## 8.4 Prediction 3: a conditional NSF differential-fading illustration

We can quantify a changing-weight interpretation without calling it an independent fit. Use the pre- and post-peak **NSF** Hα means and [N II]/Hα ratios in Table 1, with the fixed N2 spectra in Table 3. Equation (44) then returns

$$
w_0=\frac{0.373060-0.22}{0.60-0.22}=0.402789,
\qquad
w_1=\frac{0.566307-0.22}{0.60-0.22}=0.911335.
\tag{49}
$$

The measured Hα ratio of these two category means is \(1.110611/3.557592=0.312181\). Substituting it and the inferred weights into equation (37) yields component decline factors 0.0463482 for H II Hα and 0.706328 for non-H II Hα. Both decline; the H II component declines much more. The non-H II fraction rises despite a roughly 29% reduction in its own luminosity.

For plotting only, take an interpolation interval \(T=1\) Gyr. Solving equation (38) for each effective rate gives

$$
\begin{aligned}
\gamma_{\alpha,\mathrm{HII}}
&=-\frac{1}{T}\ln\!\left[
\frac{(1-w_1)\mathcal L_{\alpha,1}}
{(1-w_0)\mathcal L_{\alpha,0}}\right]
=3.07157\ \mathrm{Gyr}^{-1},\\
\gamma_{\alpha,\mathrm{non-HII}}
&=-\frac{1}{T}\ln\!\left[
\frac{w_1\mathcal L_{\alpha,1}}{w_0\mathcal L_{\alpha,0}}\right]
=0.347676\ \mathrm{Gyr}^{-1}.
\end{aligned}
\tag{50}
$$

The difference is 2.72390 Gyr\(^{-1}\), giving the explicit weight evolution through equation (39). If \(T\) were changed, the fitted rates would scale as \(1/T\); the endpoint component fractions and decline factors would remain unchanged. These are effective interpolation rates of category means, not observed temporal lifetimes of compact regions or DIG. The endpoint Hα and N2 agreement is by construction, using two observables to fix the amplitudes and two rates.

![Figure 2. Fixed partition and conditional differential fading. Left: the H I-only gas solution makes total, compact, and diffuse Hα decline together; the curves coincide when normalized by their own initial luminosity. Right: the NSF Hα/N2-conditioned interpolation makes both component luminosities decline while the non-H II weight rises. Its time coordinate is assumed for presentation and is not a measured pre- to post-peak duration.](assets/20260930_HI_stripping_line_mixing/figure02_photon_partition_and_fading.png)

The two experiments have different normalizations and purposes. The first starts at the pre-peak SF scale and predicts a gas-driven decline. The second interpolates NSF class means and demonstrates mixture arithmetic. The latter is not generated by the first. With the same retained-photon fraction and a gas-shaped OB source decline of 0.7707, a total NSF Hα factor of 0.3122 would require the captured fraction to fall by a factor \(0.3122/0.7707\simeq0.405\), or require a different OB source history. Even that endpoint budget adjustment would not determine the split into compact and diffuse absorption. Treating the illustrative luminosity rates as two independent long-lived stellar sources would be physically incorrect for an OB-leakage baseline.

We can nevertheless check one explicit photon-budget realization. Normalize a hypothetical young source to \(\mathcal L_{\alpha,0}/0.95\) at the NSF initial scale, and assign it the same *normalized* ionizing response as the gas calculation. Dividing the two luminosities from equation (38) by this source, as required by equation (28), gives \((f_{\mathrm{HII}},f_{\mathrm{leak}},f_{\mathrm{unabs}})=(0.56735,0.38265,0.05000)\) initially and \((0.03412,0.35070,0.61518)\) at the assumed endpoint. The script verifies that all three fractions remain admissible throughout the interval. This realizes the interpolation with one OB source, but requires compact absorption to collapse and the fraction unavailable to the two modeled components to increase strongly. These fractions are inferred requirements of the assumed interpolation/source pairing, not measured photon fractions or predictions of H I stripping. In particular, \(f_{\mathrm{unabs}}\) is defined relative to the modeled components and cannot be interpreted as a measured galaxy-wide escape fraction.

## 8.5 Prediction 4: other lines test the conditional N2 interpretation

Once the NSF weights are fixed from N2, [S II]/Hα and [O III]/Hβ are predictions of Table 3 rather than further adjustable fits. Their endpoint values are:

**Table 4.** Conditional interpolation compared with NSF means. Residual is \(R_{\mathrm{pred}}/R_{\mathrm{data}}-1\). It is not an error-normalized significance.

| Stage | Ratio | Data | Predicted | Residual |
|:--|:--|--:|--:|--:|
| Pre-peak | N2 | 0.373060 | 0.373060 | 0% by construction |
| Post-peak | N2 | 0.566307 | 0.566307 | 0% by construction |
| Pre-peak | S2 | 0.377298 | 0.300697 | −20.3% |
| Post-peak | S2 | 0.416836 | 0.427834 | +2.64% |
| Pre-peak | O3 | 0.897783 | 0.650558 | −27.5% |
| Post-peak | O3 | 0.752765 | 0.752267 | −0.066% |

The chosen spectra reproduce the post-peak NSF ratios approximately, partly because their endpoints were guided by that range. They fail the pre-peak NSF spectrum. In particular, the observed pre-peak O3 ratio exceeds both fixed O3 endpoints. No admissible weight can repair that mismatch. Moreover, with \(R_{\mathrm{OIII}}^{\mathrm{non-HII}}>R_{\mathrm{OIII}}^{\mathrm{HII}}\), the inferred increasing weight predicts an increasing O3 ratio, while the measured pre/post NSF means decrease. These failures are scientifically more informative than merely displaying a successful N2 interpolation.

![Figure 3. A multi-ratio test of the same inferred NSF weight. The shaded interval is the range allowed by the chosen fixed component spectra. The N2 agreement is imposed, while S2 and O3 are predictions. The pre-peak O3 point lies outside the allowed interval. These mean-value mismatches do not identify the missing physics or their statistical significance without measurement and endpoint uncertainties.](assets/20260930_HI_stripping_line_mixing/figure03_multi_ratio_falsification.png)

All line luminosities in the interpolation nevertheless fade: [N II] declines by approximately 0.474, [S II] by 0.444, and [O III] by 0.361 relative to their own predicted initial values. The last two factors are model values, not matches to the observed line decline factors. Thus differential fading can demonstrate “fainter forbidden lines but larger ratios” algebraically while still failing the joint spectrum. This is precisely why numerator and denominator luminosities, as well as ratios, must be compared.

## 8.6 Direct single-ratio inversions across all six means

**Table 5.** Equation (44) applied separately to Table 1 with Table 3 spectra. O3 entries are Hβ weights; under the imposed equal decrements they should equal the Hα weights. Values outside [0,1] are inadmissible under these spectra.

| Stage/category | \(\widehat w_{\mathrm{NII}}\) | \(\widehat w_{\mathrm{SII}}\) | \(\widehat w_{\mathrm{OIII}}\) |
|:--|--:|--:|--:|
| Pre-peak SF | 0.088 | 0.090 | 0.128 |
| Close-to-peak SF | 0.255 | 0.231 | **−1.597** |
| Post-peak SF | 0.260 | 0.112 | **−1.800** |
| Pre-peak NSF | 0.403 | 0.709 | **1.639** |
| Close-to-peak NSF | 0.570 | 0.743 | **−0.511** |
| Post-peak NSF | 0.911 | 0.867 | 0.914 |

These are aggregate effective-weight diagnostics of ratios of mean luminosities under common assumed spectra. They are not the mean of individual-region DIG fractions, nor the area fraction labeled NSF. Several O3 ratios fall outside the chosen endpoint interval, and some otherwise admissible ratios imply different weights. Endpoint recalibration, local-condition differences, spectral evolution, another source, or dust/systematic effects must therefore be examined. Since there are no independently isolated component spectra here, this table cannot establish that DIG intrinsic ratios evolve while H II intrinsic ratios remain fixed.

The numerical result supports keeping fixed spectra as a **local starting hypothesis to test**, while rejecting the claim that the particular previous fixed spectra already reproduce all MAUVE behavior. Allowing an arbitrary spectrum and weight in every stage would remove the hypothesis's predictive content. The next useful measurement is an independent component constraint, not unconstrained extra fitting freedom.

# 9. Physical interpretation, limits, and the next decisive measurements

## 9.1 What the restricted model explains

The gas model explains why star formation can decline even if ram pressure directly removes only H I: the molecular source \(\Sigma_\Phi\) falls, the retained molecular reservoir responds with a delay, and a late decline is approximately exponential under constant coefficients. Lower initial molecular columns and stronger atomic-supply losses can make outer regions reach observational detection thresholds earlier than gas-rich inner patches. This is a plausible direction for outer ND growth, rather than a prediction of its occupancy fraction without a distribution of patch properties and an observational selection function.

The luminosity model explains how NSF ratios could increase while every emission line becomes fainter: compact Balmer emission can lose a greater fraction than diffuse Balmer emission, shifting the mixture toward a spectrum with higher low-ionization ratios. A growing NSF *area fraction* can then coexist with fading NSF *luminosity*. To establish this in the data, we must show that the spectra and relative photon allocation behave as required, and forward-apply the actual BPT, EW, dispersion, reliability, and detection criteria. The model does not identify \(w\) with an NSF occupancy fraction.

## 9.2 What is not yet a self-contained quantitative explanation

Two links remain open. First, the assumed normal-disc depletion time gives a weak one-Gyr SFR decline after direct H2 stripping is removed. Second, the reservoir equations do not determine the relative compact/diffuse absorption fractions, and common OB fading alone keeps them fixed. The conditional N2 illustration proves a possible luminosity-mixture behavior, but its fixed O3 spectrum fails part of the representative data. These limitations prevent us from claiming that the revised simple model already fits the surviving-SF attenuation, NSF spectra, and ND/NSF fractions simultaneously.

An OB-leakage DIG baseline also requires sufficient diffuse gas and an adequate continuing source of OB photons, locally or from other regions. Recombination balance is appropriate only while its response is sufficiently rapid; an ionized reservoir without new photons does not generally preserve a Gyr-scale emission floor. Shocks, evolved stars, or an AGN can contribute to NSF and need targeted checks. The physical causes should be examined using spatial morphology, line widths, stellar populations, and additional diagnostic lines, rather than added by default as free components.

## 9.3 A practical sequence for improving the model

First constrain local H I and H2 columns and molecular depletion times on matched resolution and mass conventions. Compare the H I-only prediction with both a closed control and an explicitly replenished control, rather than infer a stripping rate from a stage attenuation alone. The literature indicates that whether H I is lost from the optical disc matters ([Fumagalli et al. 2009](https://doi.org/10.1088/0004-637X/697/2/1811)); losing only the remote outer reservoir need not immediately suppress an inner patch's molecular supply.

Then test H II-dominated component spectra on matched support and independent morphology, as in section 7.2. Calibrate diffuse spectra in noncentral, low-brightness regions while controlling PSF spillover and inspecting shocks/old-star context. Use separate galaxies or held-out regions to avoid fitting and validating the same endpoints. Apply a joint flux model and ask whether one weight reproduces N2, S2, O3, the Balmer lines, and independently estimated diffuse emission. Finally, forward-model the categories using the unchanged MAUVE selection and whole-galaxy uncertainty propagation.

Our preferred analytical starting point remains constant gas coefficients and fixed component spectra within a restricted condition bin. This choice makes the derivation interpretable and the failures explicit. The appendices show what survives when these quantities evolve. Evolution should be introduced through a constrained physical or empirical relation, not by using a new unconstrained endpoint for every observation.

# Appendix A. Evolving gas coefficients and external supply

## A.1 General integrating-factor solution without direct H2 stripping

We now allow \(\tau_{\mathrm{conv}}(t)\), \(\tau_{\mathrm{dep}}(t)\), \(R(t)\), \(\lambda(t)\), and \(\gamma_{\mathrm{strip,HI}}(t)\) to evolve. Their composite rates retain the definitions in equation (5) at each time. Let \(\Sigma_{\mathrm{acc,HI}}(t)\) be an external mass-supply-rate surface density into the atomic reservoir. This differs from \(\Sigma_\Phi\), which transfers mass *between* the modeled reservoirs. Assume no other phase-exchange or transport terms for now:

$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{HI}}}{dt}+\gamma_{\mathrm{HI}}(t)\Sigma_{\mathrm{HI}}
&=\Sigma_{\mathrm{acc,HI}}(t),\\
\frac{d\Sigma_{\mathrm{H_2}}}{dt}+\gamma_{\mathrm{H_2}}(t)\Sigma_{\mathrm{H_2}}
&=\frac{\Sigma_{\mathrm{HI}}(t)}{\tau_{\mathrm{conv}}(t)}.
\end{aligned}
\tag{51}
$$

Define the reservoir survival factor between arrival time \(u\) and observation time \(t\), with \(t\geq u\), by

$$
\mathcal E_j(t,u)
\equiv\exp\!\left[-\int_u^t\gamma_j(v)\,dv\right],
\qquad j\in\{\mathrm{HI},\mathrm{H_2}\}.
\tag{52}
$$

The factor is dimensionless; \(u\) and \(v\) are time-integration variables. To see the integrating-factor solution explicitly, multiply the first line of equation (51) by \(\exp[\int_0^t\gamma_{\mathrm{HI}}(v)dv]\). Its left-hand side becomes the derivative of that factor times \(\Sigma_{\mathrm{HI}}\). Integrating the product from 0 to \(t\), then dividing by the factor at \(t\), gives

$$
\Sigma_{\mathrm{HI}}(t)
=\Sigma_{\mathrm{HI},0}\mathcal E_{\mathrm{HI}}(t,0)
+\int_0^t\Sigma_{\mathrm{acc,HI}}(u)\mathcal E_{\mathrm{HI}}(t,u)\,du.
\tag{53}
$$

Repeat the same multiplication and integration for the molecular line of equation (51):

$$
\Sigma_{\mathrm{H_2}}(t)
=\Sigma_{\mathrm{H_2},0}\mathcal E_{\mathrm{H_2}}(t,0)
+\int_0^t\frac{\Sigma_{\mathrm{HI}}(u)}{\tau_{\mathrm{conv}}(u)}
\mathcal E_{\mathrm{H_2}}(t,u)\,du.
\tag{54}
$$

No new rate is required for evolving parameters. Constant coefficients and zero external source reduce these equations to (6) and (12). Nonnegative source and survival factors preserve positive gas. The explicit SFR is equation (54) divided by the **current** \(\tau_{\mathrm{dep}}(t)\). In particular, the original lower bound generalizes to

$$
\frac{\Sigma_{\mathrm{SFR}}(t)}{\Sigma_{\mathrm{SFR},0}}
\geq\frac{\tau_{\mathrm{dep}}(0)}{\tau_{\mathrm{dep}}(t)}
\exp\!\left[-\int_0^t
\frac{1-R(u)+\lambda(u)}{\tau_{\mathrm{dep}}(u)}\,du\right].
\tag{55}
$$

Thus the retained-reservoir response remains an integral of past supply, but a constant exponential decay rate is no longer guaranteed. A time-increasing depletion time can lower SFR substantially even while molecular gas persists. A time-dependent stripping history need not be compressed into one coefficient unless the interval approximation is adequate. These equations are generalizations derived here from the same balances, not hydrodynamic RPS solutions.

## A.2 A maintained external atomic supply provides a different control

For constant coefficients and a constant \(\Sigma_{\mathrm{acc,HI}}\), define the asymptotic atomic column

$$
\Sigma_{\mathrm{HI},\infty}
\equiv\frac{\Sigma_{\mathrm{acc,HI}}}{\gamma_{\mathrm{HI}}},
\qquad
\Sigma_{\mathrm{H_2},\infty}
\equiv\frac{\Sigma_{\mathrm{HI},\infty}}
{\tau_{\mathrm{conv}}\gamma_{\mathrm{H_2}}}.
\tag{56}
$$

The source integral in equation (53) is
\(\Sigma_{\mathrm{acc,HI}}(1-e^{-\gamma_{\mathrm{HI}}t})/\gamma_{\mathrm{HI}}\). Therefore \(\Sigma_{\mathrm{HI}}(t)=\Sigma_{\mathrm{HI},\infty}+(\Sigma_{\mathrm{HI},0}-\Sigma_{\mathrm{HI},\infty})e^{-\gamma_{\mathrm{HI}}t}\). The molecular source is now a constant term plus one exponential. The constant-source integral contributes \(\Sigma_{\mathrm{H_2},\infty}(1-e^{-\gamma_{\mathrm{H_2}}t})\); the decaying-source integral is the same elementary integral as equation (10). Collecting terms gives

$$
\begin{aligned}
\Sigma_{\mathrm{H_2}}(t)
&=\Sigma_{\mathrm{H_2},\infty}
+(\Sigma_{\mathrm{H_2},0}-\Sigma_{\mathrm{H_2},\infty})e^{-\gamma_{\mathrm{H_2}}t}\\
&\quad+\frac{\Sigma_{\mathrm{HI},0}-\Sigma_{\mathrm{HI},\infty}}
{\tau_{\mathrm{conv}}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}
\tag{57}
$$

For a system initially maintained without stripping, choose \(\Sigma_{\mathrm{acc,HI}}=\Sigma_{\mathrm{HI},0}/\tau_{\mathrm{conv}}\) and the molecular balance in equation (16). If this external supply continues unchanged after switching on atomic stripping, equations (5) and (56) imply

$$
\frac{\Sigma_{\mathrm{SFR},\infty}}{\Sigma_{\mathrm{SFR},0}}
=\frac{\Sigma_{\mathrm{HI},\infty}}{\Sigma_{\mathrm{HI},0}}
=\frac{1}{1+\gamma_{\mathrm{strip,HI}}\tau_{\mathrm{conv}}}.
\tag{58}
$$

Here the source is not shut off, so the SFR tends to a nonzero asymptote rather than zero. The otherwise identical no-RPS control remains constant. With Table 2 coefficients, the asymptotic factor is approximately 0.621. The approach to it is delayed by the retained molecular reservoir. This case isolates RPS against ongoing replenishment and demonstrates why specifying the upstream source is essential to any quantitative interpretation.

## A.3 Local enhancement with no direct molecular stripping

Differentiate equation (1) when \(\tau_{\mathrm{dep}}\) varies:

$$
\begin{aligned}
\frac{d\ln\Sigma_{\mathrm{SFR}}}{dt}
&=\frac{1}{\Sigma_{\mathrm{H_2}}}\frac{d\Sigma_{\mathrm{H_2}}}{dt}
-\frac{d\ln\tau_{\mathrm{dep}}}{dt}\\
&=\frac{\Sigma_\Phi}{\Sigma_{\mathrm{H_2}}}
-\frac{1-R+\lambda}{\tau_{\mathrm{dep}}}
-\frac{d\ln\tau_{\mathrm{dep}}}{dt}.
\end{aligned}
\tag{59}
$$

A local enhancement requires this right-hand side to be positive. A temporarily faster conversion or locally enhanced supply can raise its first term; a decreasing depletion time makes its final term positive. Both provide possible compression responses without direct molecular stripping. By contrast, the initially balanced, constant-coefficient closed solution in section 4 is monotonic and cannot yield a local enhancement on its own. Evidence for increased molecular content in early-stage Virgo outskirts offers motivation for a supply/compression extension ([Brown et al. 2023](https://arxiv.org/abs/2308.10943)); it does not establish a unique pressure-to-\(\tau_{\mathrm{conv}}\) or pressure-to-\(\tau_{\mathrm{dep}}\) law for our patches. A directional SFR residual must also distinguish actual enhancement above its control from weaker suppression. No enhancement amplitude is fitted in this revision.

# Appendix B. Evolving spectra, attenuation, and nonlocal photons

## B.1 What survives when the component ratios evolve?

Equations (32)–(34), derived by adding luminosities, remain valid. What fails is the fixed-segment restriction in equation (43) and the attribution of every ratio change to \(dw/dt\). For a component \(j\), write the conditional spectrum as

$$
R_k^j(t)=\mathcal R_k^j\!\left[Z(t),U_j(t),n_{e,j}(t),T_{e,j}(t),\mathrm{SED}_j(t),\mathrm{N/O}(t)\right].
\tag{60}
$$

Here \(\mathcal R_k^j\) is a photoionization/emission mapping, and SED denotes the shape of the ionizing spectral energy distribution. This equation defines dependencies rather than a new calibration. For variables summarized as \(\theta_p\), with \(p\) indexing the stated physical inputs, the chain rule is

$$
\frac{dR_k^j}{dt}
=\sum_p\frac{\partial\mathcal R_k^j}{\partial\theta_p}
\frac{d\theta_p}{dt}.
\tag{61}
$$

For example, define the dimensionless ionization parameter \(U_j=\Phi_{\mathrm H,j}/(n_{\mathrm H,j}c)\), where \(\Phi_{\mathrm H,j}\) is incident ionizing photon flux, \(n_{\mathrm H,j}\) is hydrogen density, and \(c\) is the speed of light. If geometry, density, and spectral shape stay fixed while the source fades, the incident flux can fall and so can \(U_j\). Fixed component ratios then need not follow. However, the ratio response has to be calculated or measured; there is no universal sign for every diagnostic. Changes in leaked-photon filtering and a harder source can particularly alter O3 ([Belfiore et al. 2022](https://doi.org/10.1051/0004-6361/202141859)).

Equation (40) specifies exactly where these responses enter. Small enough changes in the component spectra would leave the changing-weight interpretation approximately valid, but the relevant comparison is quantitative:

$$
\left|(1-w)\frac{dR_k^{\mathrm{HII}}}{dt}
+w\frac{dR_k^{\mathrm{non-HII}}}{dt}\right|
\ll
\left|(R_k^{\mathrm{non-HII}}-R_k^{\mathrm{HII}})\frac{dw}{dt}\right|.
\tag{62}
$$

This is a sufficient condition for the *net* intrinsic-evolution contribution to be small relative to the weight term for that diagnostic. Small net contribution can still conceal cancellation between larger endpoint changes; independently constraining each is stronger. The analogous finite condition follows from equation (41). We cannot assert either inequality for MAUVE from its class ratios alone.

Allowing an evolving DIG spectrum while retaining a constant H II spectrum is a reasonable restricted next hypothesis **if** the H II-dominated test supports that restriction and independent diffuse data favor evolution. It is not compelled by the present N2 trends. Allowing both endpoints to evolve without external constraints makes a two-component decomposition underidentified. A hierarchical spectrum tied to independently measured local conditions would be more meaningful than a free endpoint in every bin.

## B.2 Different attenuation and Balmer decrements

If \(A_\ell^j\) is component attenuation in magnitudes, the measured attenuated luminosity is

$$
\mathcal L_\ell^{\mathrm{meas}}
=\sum_j\mathcal L_\ell^{j,\mathrm{intrinsic}}
10^{-0.4A_\ell^j}.
\tag{63}
$$

With different attenuation, the relevant observed Balmer weights are wavelength dependent. Correcting the total mixed spectrum with a single screen does not generally reconstruct each intrinsic component separately. Equations (32)–(34) can always be written using each component's *measured* luminosities and measured ratios, but an interpretation in terms of intrinsic templates requires a consistent attenuation model. Nearby line pairs reduce the effect of differential attenuation within a component; they do not eliminate its effect on the relative component weight.

The numerical unequal-decrement identity check uses \(w_\alpha=0.4\), \(B_{\mathrm{HII}}=2.86\), \(B_{\mathrm{non-HII}}=3.10\). It gives \(w_\beta\simeq0.381\), showing explicitly that one cannot reuse the Hα weight unmodified. The value 3.10 is a diagnostic example, not an inferred MAUVE DIG decrement.

## B.3 Retained gas is not automatically a persistent emitter

If the source disappears, recombining gas changes its ionization state. The characteristic hydrogen recombination time under specified density and temperature is \(\tau_{\mathrm{rec}}=(\alpha_{\mathrm B}n_e)^{-1}\). This follows from a recombination loss per ion of approximately \(\alpha_{\mathrm B}n_e\), using the coefficients of [Storey & Hummer (1995)](https://doi.org/10.1093/mnras/272.1.41). Because the density and ionized fraction evolve, it is not a general exponential light-curve law. Retained neutral or molecular gas supplies potential material; it does not supply ionizing energy. The model's late diffuse emission therefore requires continuing OB photons, an explicitly modeled external photon field, or another source.

At resolved scales, the source illuminating a DIG patch can be elsewhere in the disc. A more general young-powered diffuse contribution may be written

$$
\mathcal L_\alpha^{\mathrm{non-HII}}(\boldsymbol x,t)
=h\nu_\alpha p_\alpha
\int G_{\mathrm{abs}}(\boldsymbol x,\boldsymbol x',t)
Q_{\mathrm{H,y}}(\boldsymbol x',t)\,d^2\boldsymbol x',
\tag{64}
$$

where \(p_\alpha=\alpha_\alpha^{\mathrm{eff}}/\alpha_{\mathrm B}\) is the Hα photon yield and \(G_{\mathrm{abs}}\) is a propagation-plus-diffuse-absorption kernel per receiving area. The element \(d^2\boldsymbol x'\) is projected source area in the same area convention as the photon surface density. A source at \(\boldsymbol x'\) contributes photons to the receiving patch at \(\boldsymbol x\). Integrating the kernel over receiving area must not allocate more than the available leaked fraction for each source. This generalization follows the concept of leakage propagation used by [Belfiore et al. (2022)](https://doi.org/10.1051/0004-6361/202141859), but no kernel is fitted here. It preserves OB-powered DIG without requiring local star formation in every emitting NSF patch. It also shows why a local Hα-derived SFR proxy in mixed emission is not automatically the local true \(\Sigma_{\mathrm{SFR}}\) from equation (1).

# Appendix C. A flux-based two-component test

## C.1 Construct one consistent spectral model

For ideally corrected fluxes/surface luminosities in the order Hα, Hβ, [N II]6583, [S II] sum, [O III]5007, define a measured five-element vector \(\boldsymbol y\). Normalize each fixed component spectrum to its own Hα:

$$
\boldsymbol s_j=
\begin{pmatrix}
1\\ B_j^{-1}\\ R_{\mathrm{NII}}^j\\ R_{\mathrm{SII}}^j\\ R_{\mathrm{OIII}}^j/B_j
\end{pmatrix},
\qquad
\boldsymbol S=(\boldsymbol s_{\mathrm{HII}},\boldsymbol s_{\mathrm{non-HII}}).
\tag{65}
$$

Let \(\boldsymbol{\mathcal L}_\alpha=(\mathcal L_\alpha^{\mathrm{HII}},\mathcal L_\alpha^{\mathrm{non-HII}})^T\). The component amplitudes retain direct physical names rather than arbitrary coefficients. Then

$$
\boldsymbol y_{\mathrm{model}}=\boldsymbol S\boldsymbol{\mathcal L}_\alpha,
\qquad \mathcal L_\alpha^{\mathrm{HII}},\mathcal L_\alpha^{\mathrm{non-HII}}\geq0.
\tag{66}
$$

This predicts both Balmer lines and all three forbidden numerators consistently. It avoids combining an independently varied measured Hβ with an imposed Hα/Hβ spectrum. If fitting uncorrected fluxes, attenuate each entry of each spectral column by its own screen according to equation (63), or fit a physically justified shared-screen restriction. More attenuation parameters require additional constraints; they should not be introduced merely to force a solution.

## C.2 Derive the unconstrained amplitudes and assess residuals

Let \(\boldsymbol C_L\) be the five-line measurement covariance matrix. The residual quadratic is

$$
\chi_L^2=
(\boldsymbol y-\boldsymbol S\boldsymbol{\mathcal L}_\alpha)^T
\boldsymbol C_L^{-1}
(\boldsymbol y-\boldsymbol S\boldsymbol{\mathcal L}_\alpha).
\tag{67}
$$

Set its gradient with respect to the two amplitudes to zero:

$$
\begin{aligned}
-2\boldsymbol S^T\boldsymbol C_L^{-1}
(\boldsymbol y-\boldsymbol S\boldsymbol{\mathcal L}_\alpha)&=\boldsymbol0,\\
\boldsymbol S^T\boldsymbol C_L^{-1}\boldsymbol S\boldsymbol{\mathcal L}_\alpha
&=\boldsymbol S^T\boldsymbol C_L^{-1}\boldsymbol y,\\
\widehat{\boldsymbol{\mathcal L}}_\alpha
&=(\boldsymbol S^T\boldsymbol C_L^{-1}\boldsymbol S)^{-1}
\boldsymbol S^T\boldsymbol C_L^{-1}\boldsymbol y.
\end{aligned}
\tag{68}
$$

The last line assumes the two spectral columns are identifiable and the matrix is invertible. If it gives negative amplitudes, solve the nonnegative constrained problem, evaluating its boundary solutions rather than clipping both amplitudes independently. Recover \(w_\alpha\) only after dividing the fitted non-H II amplitude by their positive sum. Compare residuals in the held-out lines or galaxies; a mathematically admissible amplitude does not guarantee a good spectrum.

Flux covariance includes shared continuum-subtraction errors, line-fitting correlations, and any propagated correction uncertainty. Endpoint spectra also have uncertainty and galaxy-to-galaxy variation. These are not included by simply assigning diagonal errors to ratios with a common Balmer denominator. A practical analysis should sample endpoint and attenuation uncertainty together with whole-galaxy resampling or a galaxy-level hierarchical model. Selection on detected lines needs a censoring/detection model when fainter regions enter the inference. The current report executes only a synthetic algebraic-recovery check and the actual mean-value diagnostics, not this full observational fit.

# Appendix D. Definitions of symbols and parameters

These tables summarize quantities introduced in the text. Constants apply only within the specified main interval; appended generalizations use their explicitly time-dependent versions. The \(\gamma\) family always means a physically specified inverse-time rate, with subscripts identifying the process or observable.

**Table D1.** Gas, time, and spatial quantities.

| Symbol | Definition and units |
|:--|:--|
| \(\boldsymbol x=(r,\varphi)\), \(\boldsymbol x'\) | Receiving and source positions in the disc; radius is a length, azimuth is an angle. Suppressing position in an ODE does not assume all patches share one gas column. |
| \(t\), \(u\), \(v\) | Elapsed model time and dummy time-integration variables. Numerical gas calculations use Gyr; not measured stage ages. |
| \(T\), \(a\) | Assumed endpoint interpolation interval and stellar-population age. Age integration must share the SFR time unit. |
| \(O(t^3)\), \(d^2\boldsymbol x'\) | Taylor terms starting at third order near the initial time, and projected source-area element, respectively. |
| \(\Sigma_{\mathrm{HI}}\), \(\Sigma_{\mathrm{H_2}}\) | Atomic-associated and molecular-associated cold-gas surface densities under one helium convention; \(M_\odot\,\mathrm{pc}^{-2}\) numerically. |
| Subscripts \(0,1,\infty\) | Initial state, comparison endpoint, and asymptotic state where defined. Endpoint 1 in line interpolation is a population mean, not the same patch observed again. |
| \(\Sigma_{\mathrm{SFR}}\) | True formed-stellar-mass rate per area. Printed in \(M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}\); not automatically the Hα proxy in a mixed NSF patch. |
| \(\tau_{\mathrm{dep}}\) | Molecular column divided by formed-mass SFR in consistent units. Two Gyr in the baseline; net consumption time also depends on return and feedback. |
| \(\Sigma_\Phi\) | Molecular gas supply-rate surface density, \(\Sigma_{\mathrm{HI}}/\tau_{\mathrm{conv}}\) under our closure; mass per area per time. |
| \(\tau_{\mathrm{conv}}\) | Atomic column divided by its transfer rate to H2; a positive effective conversion time, not necessarily the replenishment time of H2. |
| \(\tau_\Phi\) | Molecular column divided by its current supply rate; \(\Sigma_{\mathrm{H_2}}/\Sigma_\Phi\). Infinite in the zero-supply limit. |
| \(R\), \(\lambda\) | Dimensionless effective prompt return fraction and feedback mass-loading factor per unit formed stellar mass. Return to the modeled molecular balance is an approximation. |
| \(\gamma_{\mathrm{strip,HI}}\), \(\gamma_{\mathrm{strip,H_2}}\) | Direct RPS fractional gas-removal rates; inverse time. The molecular coefficient is explicitly zero in the modeled cases. |
| \(\gamma_{\mathrm{HI}}\) | Total atomic disappearance rate \(1/\tau_{\mathrm{conv}}+\gamma_{\mathrm{strip,HI}}\); conversion is included, so this is not the stripping rate alone. |
| \(\gamma_{\mathrm{H_2}}\) | Molecular consumption/feedback coefficient \((1-R+\lambda)/\tau_{\mathrm{dep}}\); no direct RPS term. |
| \(\Gamma_{\mathrm{SFR}}\) | Instantaneous logarithmic SFR decline rate \(-d\ln\Sigma_{\mathrm{SFR}}/dt\); variable even when gas coefficients are fixed. |
| \(f_{\mathrm{SFR}}\) | Target dimensionless SFR attenuation factor in equation (24), not a photon fraction. |
| \(\Sigma_{\mathrm{acc,HI}}\) | External supply to the atomic reservoir; mass per area per time. It is zero in the main closed calculation. |
| \(\mathcal E_j(t,u)\) | Dimensionless survival factor of reservoir \(j\) for material present at time \(u\), equation (52). |
| \(\Sigma_*\), \(M_*\), \(R_e\) | Stellar mass surface density, galaxy total stellar mass, and effective radius. Matching/control variables, not fitted reservoir coefficients. |

**Table D2.** Ionizing photons and emission luminosities.

| Symbol | Definition and units |
|:--|:--|
| \(q_{\mathrm{H,y}}(a)\) | Young-population hydrogen-ionizing production rate per formed stellar mass at age \(a\); photons s\(^{-1}\,M_\odot^{-1}\). |
| \(Q_{\mathrm{H,y}}\) | Young ionizing production per projected area; photons s\(^{-1}\,\mathrm{kpc}^{-2}\), equation (25). |
| \(Q_{\mathrm{abs}}^j\) | Photons absorbed by hydrogen per second in component volume \(j\), before division by projected area in equation (26). |
| \(n_e,n_p,n_{\mathrm H}\) | Electron, proton, and total hydrogen number densities; cm\(^{-3}\). They are different quantities. |
| \(dV\), \(h\nu_\alpha\) | Volume element and Hα photon energy; cm\(^3\) and erg in equation (26). \(h\) is Planck's constant and \(\nu_\alpha\) the line frequency. |
| \(\alpha_{\mathrm B}\), \(\alpha_\alpha^{\mathrm{eff}}\), \(p_\alpha\) | Case-B and effective Hα recombination coefficients in cm\(^3\) s\(^{-1}\), and their dimensionless Hα photon-yield ratio. Condition dependent. |
| \(L_\alpha^j\), \(\mathcal L_\ell^j\) | Component Hα luminosity (erg s\(^{-1}\)) and component line surface luminosity (erg s\(^{-1}\,\mathrm{kpc}^{-2}\)). Calligraphic quantities are per area. |
| \(\ell\), \(\alpha\), \(\beta\) | Generic emission-line label; \(\alpha\) and \(\beta\) identify Hα and Hβ. In recombination-coefficient symbols, \(\alpha\) instead denotes a coefficient. |
| \(j\) | Component label HII or non-HII in luminosity equations; reservoir label HI or H2 only in \(\mathcal E_j\). Its set is specified in each use. |
| \(\mathcal L_\alpha^{\mathrm{young,max}}\) | Full-hydrogen-absorption Hα normalization of the young source under the chosen recombination yield. Not an extra gas component. |
| \(f_{\mathrm{HII}},f_{\mathrm{leak}},f_{\mathrm{unabs}}\) | Mutually exclusive fractions of young photons absorbed in compact gas, absorbed in diffuse gas after leakage, or unavailable to either modeled hydrogen-emitting component. |
| \(f_{\mathrm{cap}}\) | Total fraction absorbed by the two modeled gas components, \(f_{\mathrm{HII}}+f_{\mathrm{leak}}\). All photon fractions are dimensionless. |
| \(C_\alpha\) | MAUVE formed-mass SFR per Hα luminosity; \(M_\odot\,\mathrm{yr}^{-1}/(\mathrm{erg\,s}^{-1})\). Same calibration acts on luminosity and surface luminosity with consistent area units. |
| \(K_\alpha(a)\), \(\tau_{\mathrm{ion}}\) | Normalized young-ionizing response kernel (inverse time) and its assumed exponential smoothing time. |
| \(\Sigma_{\mathrm{SFR}}^{\mathrm{ion}}\) | SFR filtered with \(K_\alpha\); same units as SFR. A model source-normalization proxy, not an independent observation. |
| \(\gamma_{\alpha,\mathrm{HII}},\gamma_{\alpha,\mathrm{non-HII}}\) | Effective Hα luminosity fading rates of equation (38), inverse time. Not atomic/molecular stripping rates or stellar-population lifetimes. |
| \(\tau_{\mathrm{rec}}\) | Characteristic recombination time \((\alpha_{\mathrm B}n_e)^{-1}\) at specified density and temperature. |
| \(G_{\mathrm{abs}}\) | Young-photon propagation and diffuse-absorption kernel per receiving area; constrained by photon conservation. |

**Table D3.** Ratios, weights, and possible spectral evolution.

| Symbol | Definition and units |
|:--|:--|
| \(R_k^j\), \(R_k^{\mathrm{obs}}\) | Dimensionless linear component and mixture ratios; \(k\) denotes N2=[N II]6583/Hα, S2=[S II] sum/Hα, or O3=[O III]5007/Hβ. |
| \(w_{\mathrm{non-HII},\alpha}\), \(w_{\mathrm{non-HII},\beta}\) | Non-H II fractions of Hα and Hβ luminosity, respectively. \(w_{\mathrm{HII},B}=1-w_{\mathrm{non-HII},B}\) for either Balmer denominator \(B\). |
| \(w\), \(w_0,w_1\), \(\widehat w_k\) | Hα weight abbreviated in stated equal-decrement/Hα contexts, its two endpoints, and a weight inferred from diagnostic \(k\). A hat is an estimator. |
| \(B_j\) | Component Hα/Hβ decrement; dimensionless, with intrinsic or measured convention explicitly stated. |
| \(\Delta\), overbar | Endpoint difference and arithmetic midpoint in equation (41); not temporal derivatives. |
| \(Z\), N/O | Gas abundance coordinate and nitrogen-to-oxygen abundance ratio. Their particular observational calibrations are not specified by the mixture algebra. |
| \(U_j,\Phi_{\mathrm H,j},c\) | Dimensionless incident-photon ionization parameter, ionizing flux in photons cm\(^{-2}\) s\(^{-1}\), and speed of light. |
| \(T_e,\mathrm{SED}_j\) | Electron temperature in K and ionizing spectral shape. They can differ between compact and diffuse gas. |
| \(\mathcal R_k^j,\theta_p\) | Conditional physical mapping for a component spectrum and its inputs; \(p\) indexes those inputs in the chain rule. |
| \(A_\ell^j\) | Component attenuation at a line wavelength, in magnitudes; not a gas-loss rate or an emitting-component amplitude. |

**Table D4.** Statistical notation.

| Symbol | Definition and role |
|:--|:--|
| \(F_k,\Delta_k^{\mathrm{HII-dom}}\) | Proposed pre-peak log-ratio reference surface in matched controls and residual in dex for an H II-dominated subset. Not fitted in this report. |
| \(\boldsymbol R,\boldsymbol d\) | Three linear ratios and the difference of the two fixed spectra in that space. |
| \(\boldsymbol C,\chi^2\) | Ratio covariance and joint one-weight residual quadratic; uncertainty in component spectra needs separate propagation. |
| \(\boldsymbol y,\boldsymbol s_j,\boldsymbol S\) | Five-line data vector, component spectrum normalized to Hα, and matrix with the two spectral columns. |
| \(\boldsymbol{\mathcal L}_\alpha\) | Two nonnegative component Hα surface-luminosity amplitudes in the flux model. |
| \(\boldsymbol C_L,\chi_L^2\) | Five-line covariance and corresponding flux-model residual quadratic. |
| Superscript \(T\), \(^{-1}\) on matrices | Transpose and matrix inverse; \(T\) here does not denote the interpolation time interval. |
| SF, NSF, ND; DIG | Observational classes versus a diffuse ionized gas component. NSF occupancy is not \(w\). |
| RPS, OB, IMF, PSF, AGN, BPT | Ram-pressure stripping; O- and B-type stars; initial mass function; point spread function; active galactic nucleus; Baldwin–Phillips–Terlevich diagnostic diagrams. |
| \(S/N_{\mathrm{POSTFIT}}\), EW, dispersion | Post-fit continuum reliability, line equivalent width, and intrinsic Hα velocity dispersion. Selection variables, not source-fraction measurements. |

# Appendix E. Provenance and executed verification

## E.1 Source-to-equation ledger

| Ingredient | Provenance and boundary |
|:--|:--|
| Local \(\Sigma_\Phi\), \(\tau_\Phi\), and regulator notation | Huang et al. (2026), section 4; the assumed linear HI→H2 transfer is added here. |
| Linear molecular law and constant depletion-time illustration | Bigiel et al. (2011); used as a restricted normal-disc approximation. |
| Consumption/return/feedback bookkeeping | Lilly et al. (2013); two-phase sinks and return-to-phase approximation specified here. |
| H I loss as a route to H2/SFR decline | Fumagalli et al. (2009); neither a universal guarantee nor the exact linear model. |
| Limits of zero direct molecular removal and constant efficiency | Boselli et al. (2014); Brown et al. (2023). |
| Equations (3)–(24), (51)–(59) | Our balances, assumptions, and explicit algebraic solutions/generalizations. |
| Young-ionizing cohort sum and tracer response | Population-synthesis/tracer framework summarized by Kennicutt & Evans (2012); the 3-Myr exponential kernel is our illustrative choice. |
| Case-B photon-to-Hα conversion and decrements | Hummer & Storey (1987); Storey & Hummer (1995); nebular-condition dependent. |
| Photon allocation and OB-leakage DIG motivation | Belfiore et al. (2022); the local fractions are definitions, not measured by that paper for MAUVE. |
| Linear H II/DIG ratio mixture antecedent | Blanc et al. (2009), equations 7–8; our generic Balmer-weight derivation follows additivity. |
| Intrinsic DIG spectral variation and limitations of pure leakage | Zhang et al. (2017); Belfiore et al. (2022). |
| Equations (30)–(50), (60)–(68) | Our mixing identities, conditional parameterization, statistical tests, and numerical evaluations; not a claimed pre-existing unique RPS solution. |

## E.2 Local inputs and reproducibility

The source table is `assets/20260914_resolved_RPS_academic_model/stage_bpt_line_profiles_with_hbeta.csv`. The matched-profile attenuation summary is `assets/20260914_resolved_RPS_academic_model/attenuation_fit_results.json`. Their original extraction is documented by the 14 September report and scripts. On 30 September, the numerical script checked SHA-256 fingerprints for the three stage-analysis notebooks, `further/SFR+Z.py`, the stage wiki FITS, and effective-radius CSV against that extraction's record. Their paths and exact hashes are saved in `source_fingerprints.json`. The scalar line export hash is saved in `numerical_audit.json`. Large emission-map FITS files were not rehashed, and no notebook or FITS map was rerun or changed.

The main numerical command was:

```bash
/opt/miniconda3/envs/ICRAR/bin/python \
  /Users/Igniz/Desktop/ICRAR/MAUVE/assets/20260930_HI_stripping_line_mixing/model_predictions.py
```

The script regenerates the anchor table, gas/null predictions, N2-conditioned NSF interpolation, per-ratio weights, endpoint residuals, three figures, and numerical audit. It does not optimize a fit to the full observations. `eta` in its numerical input dictionary is the conventional code name for the mass-loading factor denoted \(\lambda\) in this report; it is zero in the example.

**Table E1.** Fresh numerical checks executed for this report.

| Check | Result | Meaning |
|:--|:--|:--|
| Analytical HI/H2 versus independent ODE integration | Maximum relative difference \(6.73\times10^{-12}\) | Checks the nondegenerate constant-coefficient solution |
| Equal-rate solution versus independent ODE | \(1.53\times10^{-15}\) | Checks the nondivergent branch of equation (13) |
| Initial molecular/SFR derivative | Zero to printed precision | Verifies the imposed instantaneous balance |
| Positive, monotonic balanced SFR and lower bound | Passed over 0–3 Gyr | Verifies implemented response and equation (23) |
| Synthetic common-weight recovery | Absolute error zero to printed precision | Checks equation (46) with an explicitly synthetic covariance; no MAUVE inference |
| Unequal-Balmer direct luminosity versus weight conversion | Difference \(1.11\times10^{-16}\) | Checks equation (34) |
| Conditional NSF Hα/N2 endpoints | Reproduced to floating-point precision | Construction check, not independent validation |
| Conditional allocation paired with a gas-shaped young source | All fractions in [0,1] and sum to unity over 0–1 Gyr | Finite-interval budget feasibility, not an RPS-derived allocation law |
| Fixed spectra against other MAUVE ratios | Substantial mean-value failures documented | No formal rejection significance without measurement/endpoint uncertainty |

PDF verification uses offline Pandoc MathML conversion, rendered-equation/DOM checks, PDF text and link audits, and visual inspection of all page contact sheets plus selected equation and figure pages. The saved `verification_summary.md` and `pdf_qa/pdf_audit.json` record the final rendering results. Numerical precision checks validate the algebra and implementation; they do not validate the physical closures or establish a fit to MAUVE.

# References

**Belfiore, F., et al. (2022).** *A tale of two DIGs: The relative role of H II regions and low-mass hot evolved stars in powering the diffuse ionised gas in PHANGS–MUSE galaxies.* A&A, 659, A26. [DOI](https://doi.org/10.1051/0004-6361/202141859); [author manuscript](https://arxiv.org/html/2111.14876v3).

**Bigiel, F., et al. (2011).** *A Constant Molecular Gas Depletion Time in Nearby Disk Galaxies.* ApJ Letters, 730, L13. [DOI](https://doi.org/10.1088/2041-8205/730/2/L13); [author manuscript](https://arxiv.org/abs/1102.1720).

**Blanc, G. A., Heiderman, A., Gebhardt, K., Evans, N. J. II, & Adams, J. (2009).** *The Spatially Resolved Star Formation Law from Integral Field Spectroscopy: VIRUS-P Observations of NGC 5194.* ApJ, 704, 842–862. [DOI](https://doi.org/10.1088/0004-637X/704/1/842); [author manuscript](https://arxiv.org/pdf/0908.2810).

**Boselli, A., et al. (2014).** *Cold gas properties of the Herschel Reference Survey. III. Molecular gas stripping in cluster galaxies.* A&A, 564, A67. [DOI](https://doi.org/10.1051/0004-6361/201322313); [author manuscript](https://arxiv.org/abs/1402.0326).

**Brown, T., et al. (2023).** *VERTICO. VII. Environmental Quenching Caused by the Suppression of Molecular Gas Content and Star Formation Efficiency in Virgo Cluster Galaxies.* ApJ, 956, 37. [DOI](https://doi.org/10.3847/1538-4357/acf195); [author manuscript](https://arxiv.org/abs/2308.10943).

**Fumagalli, M., Krumholz, M. R., Prochaska, J. X., Gavazzi, G., & Boselli, A. (2009).** *Molecular hydrogen deficiency in H I-poor galaxies and its implications for star formation.* ApJ, 697, 1811–1821. [DOI](https://doi.org/10.1088/0004-637X/697/2/1811); [author manuscript](https://arxiv.org/abs/0903.3950).

**Huang, R., et al. (2026).** *MAUVE–MUSE: When Metallicity Follows or Fights Star Formation—A Mass-Dependent Inversion in Virgo Galaxies.* MNRAS. [DOI](https://doi.org/10.1093/mnras/stag1019); [author manuscript, section 4](https://arxiv.org/html/2605.31412v1).

**Hummer, D. G., & Storey, P. J. (1987).** *Recombination-line intensities for hydrogenic ions—I. Case B calculations for H I and He II.* MNRAS, 224, 801–820. [DOI](https://doi.org/10.1093/mnras/224.3.801).

**Kennicutt, R. C., Jr., & Evans, N. J. II (2012).** *Star Formation in the Milky Way and Nearby Galaxies.* ARA&A, 50, 531–608. [DOI](https://doi.org/10.1146/annurev-astro-081811-125610); [author manuscript](https://arxiv.org/abs/1204.3552).

**Lilly, S. J., Carollo, C. M., Pipino, A., Renzini, A., & Peng, Y. (2013).** *Gas-regulation of galaxies: the evolution of the cosmic sSFR, the metallicity-mass-SFR relation and the stellar content of haloes.* ApJ, 772, 119. [DOI](https://doi.org/10.1088/0004-637X/772/2/119); [author manuscript](https://arxiv.org/abs/1303.5059).

**Storey, P. J., & Hummer, D. G. (1995).** *Recombination line intensities for hydrogenic ions—IV. Total recombination coefficients and machine-readable tables for Z=1 to 8.* MNRAS, 272, 41–48. [DOI](https://doi.org/10.1093/mnras/272.1.41); [published paper](https://academic.oup.com/mnras/article/272/1/41/967214).

**Zhang, K., et al. (2017).** *SDSS-IV MaNGA: The impact of diffuse ionized gas on emission-line ratios, interpretation of diagnostic diagrams, and gas metallicity measurements.* MNRAS, 466, 3217–3243. [DOI](https://doi.org/10.1093/mnras/stw3308); [author manuscript](https://arxiv.org/abs/1612.02000).
