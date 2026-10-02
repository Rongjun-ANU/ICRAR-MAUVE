"""Build the dated revision from preserved Oct 1 derivations and explicit changes.

Equation labels are resolved in document order; the prior report is read-only.
"""
from pathlib import Path
import re
import json

ROOT=Path('/Users/Igniz/Desktop/ICRAR/MAUVE')
OUT=Path(__file__).resolve().parent
OLD=(ROOT/'20261001_Connected_HI_Stripping_and_HOLMES_Line_Ratio_Model.md').read_text()
AUDIT=json.loads((OUT/'numerical_audit.json').read_text())

def part(start,end):
    text=OLD[OLD.index(start):OLD.index(end)]
    # Preserve old equation references while labels are reordered.
    def refs(m):
        return re.sub(r'\((\d+)\)',lambda n:'[eq:'+n[1]+']',m[0])
    text=re.sub(r'[Ee]quations? \(\d+\)(?:(?:--|, | and )\(\d+\))*',refs,text)
    return text

def time_explicit(text):
    # Fixed-position dynamical quantities retain their time arguments.
    for token in [r'\Sigma_{\mathrm{HI}}',r'\Sigma_{\mathrm{H_2}}',r'\Sigma_{\mathrm{SFR}}',r'\Sigma_\Phi']:
        text=re.sub(re.escape(token)+r'(?![_^]|\()',lambda m:m[0]+'(t)',text)
    text=text.replace(r'\widetilde\Sigma_{\mathrm{HI}}(t)',r'\widetilde\Sigma_{\mathrm{HI}}')
    return text

sections=[]
sections.append(r'''---
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

''')
gas=part('# 2. Definition','# 4. From')
gas=gas[:gas.index('## 3.5')]
gas=time_explicit(gas)
gas=gas.replace('All coefficients are initially constant.',
    r'In Model 0, $\tau_{\mathrm{conv}}$, $\tau_{\mathrm{dep}}$, $R$, $\lambda$, and $\gamma_{\mathrm{strip}}$ are constant in time at fixed $\boldsymbol{x}$. Their possible spatial dependence is retained conceptually. Constant $\tau_{\mathrm{dep}}$ means constant molecular efficiency, not constant SFR.')
gas=gas.replace('If the two rates are equal, the integrand in equation [eq:11] is unity.',
    r'If the HI-reservoir response rate and molecular-consumption rate are equal, $\gamma_{\mathrm{HI}}=\gamma_{\mathrm{H_2}}$, then $1/\tau_{\mathrm{conv}}+\gamma_{\mathrm{strip}}=(1-R+\lambda)/\tau_{\mathrm{dep}}$. The integrand in equation [eq:11] is unity.')
gas=gas.replace('For unequal rates, the integral is',r'For $\gamma_{\mathrm{HI}}\ne\gamma_{\mathrm{H_2}}$, the integral is')
gas=gas.replace('## 3.4 Derive the enhancement peak explicitly','## 3.4 A temporal maximum for an initially over-supplied region')
gas=gas.replace('## 3.2 Solve the molecular reservoir without omitting the integrating-factor algebra','## 3.2 Molecular reservoir and explicit SFR solution')
gas=gas.replace('This is a positive-time peak when',r'This is a positive-time peak of the temporal SFR history, not a measurement of spatial enhancement. It occurs when')
sections.append(gas)
sections.append(r'''
## 3.5 A spatial SFR excess is not a positive temporal derivative

The observable motivated by the NGC4654 gradient discussion is a spatial excess relative to another region or a reference relation. We define that comparison explicitly:

$$
\Delta_{\mathrm{spatial}}\log_{10}\Sigma_{\mathrm{SFR}}(t)
\equiv\log_{10}\frac{\Sigma_{\mathrm{SFR}}^{\mathrm{leading}}(t)}
{\Sigma_{\mathrm{SFR}}^{\mathrm{reference}}(t)}.
\tag{spatial}
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
\tag{initial}
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
\tag{young}
$$

Here $f_{\mathrm{young}}$ is the constant fraction of the young ionizing photon budget absorbed by hydrogen in the modeled region, with $0<f_{\mathrm{young}}\leq1$ in this local approximation. $C_\alpha$ is the fully absorbed Halpha-to-SFR calibration, not its inverse. We use the existing MAUVE value $C_\alpha=4.9835821\times10^{-42}\ M_\odot\,\mathrm{yr^{-1}}/(\mathrm{erg\,s^{-1}})$; division converts $M_\odot\,\mathrm{yr^{-1}\,kpc^{-2}}$ to $\mathrm{erg\,s^{-1}\,kpc^{-2}}$. The numerical benchmark sets $f_{\mathrm{young}}=1$. The calibration depends on the stellar population and IMF and is not universal.

Appendix A derives the normalized response kernel and its exact exponential solution. For the illustrative 3-Myr response and the present gas coefficients, the normalized one-Gyr response is 0.79937, compared with instantaneous $F_{\mathrm{SFR}}=0.79867$: a relative correction of 0.0881%. This supports equation [eq:young] for the smooth Model 0 history. It does not justify an instantaneous Halpha jump at a sudden SFR discontinuity. If the initial gas accumulation was recent on a few-Myr timescale, the actual prehistory must enter the response calculation.

There is also a spatial condition: young emission in an NSF patch may be powered by photons from neighbouring star-forming regions. In that case a local Halpha luminosity need not trace the local SFR. Equation [eq:young] is then a coarse-grained or local-absorption approximation. This limitation is especially relevant to the illustrative gas normalization in section 8; Appendix E gives the nonlocal form.

# 5. A compact, physically normalized HOLMES contribution

Let $\Sigma_*^{\mathrm{old}}$ be the current mass surface density in the old population, including its associated remnants, and let $q_{\mathrm{H,HOLMES}}$ be the production rate of hydrogen-ionizing photons per unit of that current mass. Their product is the emitted photon production per area. Multiplying by $f_{\mathrm{abs,HOLMES}}$ gives the photon rate actually absorbed by hydrogen in the gas assigned to the region:

$$
\begin{aligned}
\mathcal Q_{\mathrm{HOLMES}}&=q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}},\\
\mathcal Q_{\mathrm{abs,HOLMES}}&=f_{\mathrm{abs,HOLMES}}\mathcal Q_{\mathrm{HOLMES}}.
\end{aligned}
\tag{35}
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
\tag{38}
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
\tag{weight}
$$

These are light fractions, not fractions of area or gas mass. Write $w_{\mathrm{HOLMES},0}=w_{\mathrm{HOLMES}}(0)$. Dividing numerator and denominator by the initial total Halpha luminosity and using equation [eq:young] gives

$$
\boxed{w_{\mathrm{HOLMES}}(t)=
\frac{w_{\mathrm{HOLMES},0}}
{w_{\mathrm{HOLMES},0}+(1-w_{\mathrm{HOLMES},0})F_{\mathrm{SFR}}(t)}.}
\tag{45}
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
\tag{weightderiv}
$$

The second line uses constant $f_{\mathrm{young}}$ and $C_\alpha$; the third uses equation [eq:17]. This connects the gas balance directly to the changing source weight. Molecular consumption exceeding replenishment makes SFR decline and the HOLMES fraction rise. An initially over-supplied region has the opposite response until its SFR maximum. Proportional fading of two young components alone would leave their mutual weight unchanged; the independently supplied HOLMES term is what changes this conclusion.

# 7. The forbidden-line ratios and the complete connection

For a forbidden line $\ell$, let $B$ denote its Balmer denominator and define the linear component ratio $R_{\ell/B}^j=\mathcal L_\ell^j/\mathcal L_B^j$, with $j$ equal to young or HOLMES. N2 is [N II]6583/Halpha, S2 is ([S II]6716+[S II]6731)/Halpha, and O3 is [O III]5007/Hbeta. Shock and AGN emission are neglected as a baseline hypothesis, not established absent by these equations.

In Model 0 we adopt $\mathcal L_\alpha^j/\mathcal L_\beta^j=2.86$ for each component. Thus

$$
\frac{\mathcal L_\beta^{\mathrm{HOLMES}}}{\mathcal L_\beta}
=\frac{\mathcal L_\alpha^{\mathrm{HOLMES}}/2.86}{\mathcal L_\alpha/2.86}
=w_{\mathrm{HOLMES}}.
\tag{common}
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
\tag{44}
$$

This luminosity-weighted identity has an HII/DIG antecedent in [Blanc et al. (2009), equations 7--8](#ref-blanc). Their numerical [S II] template is for a single line and is not imported for our doublet sum. Ratios mix linearly; logarithmic BPT coordinates are taken only after addition.

For constant component spectra, differentiate equation [eq:44]:

$$
\frac{dR_{\ell/B}(t)}{dt}
=(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\frac{dw_{\mathrm{HOLMES}}(t)}{dt}.
\tag{ratioderiv}
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
\tag{51}
$$

Here $F_{\mathrm{SFR}}(t)$ is explicitly given by equation [eq:16], or by its equal-rate limit from equation [eq:15]. The second line applies to Hbeta ratios as well because of equation [eq:common]. All luminosity factors and rate coefficients are defined; no arbitrary residual luminosity function is introduced. In the restricted model, Halpha is the observable connection between declining molecular supply and evolving line-ratio weights.

Only after deriving that relation should we identify its observational scope. SF and NSF are selections, not ionizing sources. SF in the source analysis requires a finite HII-selected SFR product, Halpha EW greater than 6 Angstrom, and intrinsic Halpha dispersion below $45\ \mathrm{km\,s^{-1}}$. ND is the joint Balmer non-detection category; NSF is the remaining Balmer-detected category outside SF. Neither component is switched off merely because a region is called SF or NSF. NSF is not synonymous with DIG, LIER, or HOLMES domination, and NSF occupancy is not $w_{\mathrm{HOLMES}}$.

The present model concerns emission within retained gas. Predicting changes in SF/NSF/ND area fractions requires the stellar continuum, EW, noise, widths, masks, and line-detection rules to be modeled and reapplied. Completely stripped outer regions also require an evolving gas absorption capacity. Those calculations are outside Model 0.

''')

# Preserve the measured anchors and the normalization test.
num=part('# 8. Numerical','## 8.3 Gas')
sections.append(num)
sections.append(r'''
## 8.3 The gas response and an elevated initial spatial state

For comparability with the previous report, subtract the fiducial HOLMES luminosity from the pre-peak NSF Halpha scale and apply $C_\alpha$ with $f_{\mathrm{young}}=1$. The resulting effective young SFR is $0.01753\ M_\odot\,\mathrm{yr^{-1}\,kpc^{-2}}$. With $\tau_{\mathrm{dep}}=2$ Gyr, equation [eq:6] assigns $\Sigma_{\mathrm{H_2},0}=35.055\ M_\odot\,\mathrm{pc^{-2}}$.

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

Equation [eq:19] yields $F_{\mathrm{SFR}}(1\ \mathrm{Gyr})=0.79867$, or a 0.0976-dex decline. The initial slope is zero and subsequent slopes are negative. The HI reservoir declines rapidly, but H2 buffers the SFR. The closed comparison with $\gamma_{\mathrm{strip}}=0$, all other quantities unchanged, gives 0.89706. Relative to that declining comparison the extra suppression is 0.0505 dex. The comparison is closed and has no fresh external supply; it is not a permanently maintained field equilibrium.

Now raise both initial gas columns by 1.5 without changing $\tau_{\mathrm{conv}}$, $\tau_{\mathrm{dep}}$, or the stripping coefficient. Linearity makes the entire stripped SFR curve 1.5 times the baseline stripped curve. Its initial SFR is 0.1761 dex above the uncompressed initial reference, but it has the same zero initial fractional slope and subsequent decline. At one Gyr its SFR is $1.5\times0.79867=1.1980$ in units of the original reference SFR, and it remains 0.1256 dex above the contemporaneous closed no-RPS reference, 0.89706.

This example makes the distinction explicit: an elevated spatial SFR and a negative temporal derivative coexist. The factor 1.5 is a chosen initial-condition illustration, not a fitted compression amplitude for NGC4654. The model evolves the accumulated gas but does not calculate the accumulation process.

![Figure 1. Model 0 at fixed conversion and depletion times. Left: SFR in units of the uncompressed initial reference; the orange curve begins with both gas columns multiplied by 1.5 and subsequently declines. Right: normalized atomic columns for the stripped and closed no-stripping cases. The orange atomic fraction coincides with the blue stripped fraction and is omitted. The elevated orange SFR is a spatial-level example, not a burst generated by switching on the HI sink.](assets/20261002_Model0_Derivation/figure01_SFR_response.png)

''')
num=part('## 8.4 How','# 9. What')
num=num.replace('F_\\alpha','F_{\\mathrm{SFR}}')
num=num.replace('equation [eq:45]','equation [eq:45]')
start=num.index('For Figure 2,')
end=num.index('## 8.5',start)
num=num[:start]+r'''Figure 2 uses only the effective young N2 ratio 0.35 and HOLMES N2 ratio 1.5. These positive constant spectra illustrate the mechanism without splitting HII and leakage. The remaining young fraction is varied directly; a hundredfold fading is not claimed to occur within the two-Gyr gas interval or within the validity of every fixed-population approximation.

![Figure 2. Effective young-star plus HOLMES fading at the bright NSF scale and a hypothetical faint scale. Left: increasing HOLMES Halpha weight. Middle: increasing N2 as total Halpha decreases. Right: [N II] still declines. The remaining young fraction decreases from left to right. These are illustrative component spectra, not a fit.](assets/20261002_Model0_Derivation/figure02_fading_and_ratio.png)

'''+num[end:]
num=num.replace('using its **observed Balmer luminosity**','using its **observed Halpha luminosity** to calculate the common Model 0 weight')
num=num.replace('equation [eq:49]', 'equation [eq:inverse]')
num=num.replace('equation [eq:55]', 'equation [eq:55]')
# Main inversion now uses common Halpha weights, including for O3.
num=num.replace(r'\mathcal L_B',r'\mathcal L_\alpha')
num=num.replace(r'\mathcal L_{B,',r'\mathcal L_{\alpha,')
num=num.replace(r'w_{\mathrm{HOLMES},B}(0)',r'w_{\mathrm{HOLMES},0}')
num=num.replace(r'w_{\mathrm{target},B}',r'w_{\mathrm{target}}')
num=num.replace('**line luminosities**','ratio-times-Halpha quantities (the actual forbidden luminosity for N2 and S2, and 2.86 times that luminosity for the idealized O3 model)')
row=AUDIT['line_budget'][2]
num=num.replace('| O3 | 3.0 | 0.8736 | 0.9510 | 0.7528 | +0.1015 |',
    f"| O3 | 3.0 | {row['r_young_initial']:.4f} | {row['ratio_post_predicted']:.4f} | 0.7528 | {row['residual_dex']:+.4f} |")
num=num.replace('O3$_{\\mathrm{HOLMES}}=-4.84$',f"O3$_{{\\mathrm{{HOLMES}}}}={row['r_holmes_required_at_fiducial']:.2f}$")
num=num.replace('constant-spectrum/photon-budget assumptions rather than the gas clock.', 'constant-spectrum/photon-budget assumptions rather than the gas clock. Model 0 uses one common Halpha weight; the actual-Hbeta variant is in Appendix D.')
num=num.replace('Post-peak Balmer luminosities','Post-peak Halpha luminosities')
num=num.replace('The second line follows by subtracting the ratio-times-Halpha', 'The second line follows by subtracting the ratio-times-Halpha')
num=num.replace('Write equation [eq:inverse] at two observed points indexed by 0 and 1.',
    'Using the common weight, write equation [eq:44] as a function of total Halpha luminosity:\n\n$$\nR_{\\ell/B}=R_{\\ell/B}^{\\mathrm{young}}+\\frac{(R_{\\ell/B}^{\\mathrm{HOLMES}}-R_{\\ell/B}^{\\mathrm{young}})\\mathcal L_\\alpha^{\\mathrm{HOLMES}}}{\\mathcal L_\\alpha}.\n\\tag{inverse}\n$$\n\nEvaluate it at two observed points indexed by 0 and 1.')
start=num.index('For a monotonically declining history')
num=num[:start]+r'''Under the instantaneous Model 0 approximation, define $F_{\mathrm{required}}\equiv F_{\alpha,\mathrm{required}}$ as the target value of $F_{\mathrm{SFR}}(t)$. Because the supplied molecular contribution is nonnegative, equation [eq:21] bounds the instantaneous decline. Taking logarithms at the time when that target is reached gives the minimum time step by step:

$$
F_{\mathrm{SFR}}(t)\geq e^{-\gamma_{\mathrm{H_2}}t},\qquad
\ln F_{\mathrm{required}}\geq-\gamma_{\mathrm{H_2}}t,
\qquad t\geq\frac{-\ln F_{\mathrm{required}}}{\gamma_{\mathrm{H_2}}}.
\tag{timebound}
$$

For $F_{\mathrm{required}}=0.30426$ and $\gamma_{\mathrm{H_2}}=0.3\ \mathrm{Gyr^{-1}}$, the no-supply minimum is 3.966 Gyr. Solving the full balanced **instantaneous** Model 0 gives 4.223 Gyr; retaining the Appendix A stellar response gives 4.226 Gyr. The approximation in section 4 therefore leaves the physical tension unchanged. These long extrapolations are diagnostics of the assumed consumption rate, not inferred infall times, and extend beyond the intended 1--2 Gyr fixed-population interval.

For a hypothetical one-Gyr constraint, the same inequality requires

$$
\gamma_{\mathrm{H_2}}\geq-\frac{\ln F_{\mathrm{required}}}{1\ \mathrm{Gyr}}
\quad\Longrightarrow\quad
\tau_{\mathrm{dep}}\leq\frac{(1-R+\lambda)(1\ \mathrm{Gyr})}{-\ln F_{\mathrm{required}}}.
\tag{deplimit}
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

'''
sections.append(num)
sections.append(r'''# 9. Interpretation and limitations

Model 0 establishes a self-consistent conditional sequence. Atomic stripping reduces the future molecular supply; existing H2 buffers the SFR decline; young-star Halpha subsequently fades; and a retained old-star photon budget can become a larger fraction of the remaining Balmer emission. Fixed positive spectra can then produce a rising forbidden-to-Balmer ratio while both lines fade. No changing HII-to-leakage partition is needed for that mechanism.

The model also separates two statements that should not be conflated. A spatially enhanced molecular column produces an elevated SFR at fixed efficiency. Its subsequent derivative can already be negative. The closed reservoir solution predicts that subsequent evolution but does not generate the compression or transport that supplied the initial column. Switching on the HI sink alone does not generate an enhancement from molecular balance.

Its numerical limitations are substantive. The chosen 2-Gyr depletion time cannot yield the required fading within one Gyr even after replenishment stops. The fiducial local HOLMES budget gives only about 1.14% and 3.65% of the pre/post NSF Halpha means. With the assumed spectra the N2 increase is too small and O3 changes in the wrong direction. Allowing freely chosen but fixed spectra still demands a negative O3 HOLMES contribution at this photon budget. Lowering the absorbed fraction does not repair these particular amplitude tests.

These results do not prove that HOLMES are absent or that all NSF regions obey the same spectrum. The points combine different galaxies, the bin-center stellar mass is only representative, and the photon yield and stellar mass conventions carry uncertainties. The model uses a local photon budget and independent emitting templates. If OB and HOLMES photons illuminate the same gas, the ionic fractions and temperature respond jointly, so the forbidden spectrum need not equal the sum of two separately fixed source spectra ([Belfiore et al. 2022, sections 5.4--5.5](#ref-belfiore)). An increasingly important HOLMES spectrum and a declining ionization parameter can act together; their quantitative effects require a physical photoionization calculation, not independent adjustment of three ratios.

For application to MAUVE, the next constraints are concrete: use CO/HI rather than NSF Halpha alone to normalize the gas reservoirs; estimate the old-population photon budget from the fitted stellar populations; test ratios and absolute luminosities against $\mathcal L_\alpha/\Sigma_*^{\mathrm{old}}$ within matched stellar-density/radius support; and assess leakage from nearby HII regions. Galaxy-level uncertainty and spatial support must accompany any comparison among stages. A successful spectral extension should predict N2, S2, O3 and Balmer emission with one consistent gas state. The present report supplies the analytical reference against which such extensions can be tested; it does not claim a fit to MAUVE occupancy fractions or the full data set.

# Appendix A. The finite Halpha response and its short-timescale limit

''')
kernel=part('Halpha traces the ionizing','# 4.') if False else part('Halpha traces the ionizing','## 4.2 Allocate')
kernel=kernel.replace('the equal-gas-rate case','the equal-gas-rate case')
sections.append(kernel)
sections.append(r'''
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
\tag{kernelderiv}
$$

The boundary term comes from the upper limit; differentiating the exponential gives the second term. Normalization ensures that constant SFR is preserved. The exponential kernel is not a 10-Myr top-hat: for $\tau_{\mathrm{ion}}=3$ Myr it assigns $1-e^{-A/\tau_{\mathrm{ion}}}$ of the weight to ages below $A$, giving 90% below 6.91 Myr and 95% below 8.99 Myr.

For equal gas rates, $\gamma_{\mathrm{HI}}=\gamma_{\mathrm{H_2}}=\gamma$, the input is $(1+t/\tau_{\Phi,0})e^{-\gamma t}$. Separating the contribution from its constant prehistory yields a finite integral without a difference of nearly equal gas rates:

$$
F_\alpha(t)=e^{-t/\tau_{\mathrm{ion}}}
\left[1+\frac{1}{\tau_{\mathrm{ion}}}
\int_0^t e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)u}
\left(1+\frac{u}{\tau_{\Phi,0}}\right)du\right].
\tag{equalkernel}
$$

The integral of $u e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)u}$ follows by integration by parts. Writing the result without an additional physical parameter,

$$
\begin{aligned}
\int_0^t u e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)u}du
&=\frac{t e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)t}}{\tau_{\mathrm{ion}}^{-1}-\gamma}
-\frac{e^{(\tau_{\mathrm{ion}}^{-1}-\gamma)t}-1}{(\tau_{\mathrm{ion}}^{-1}-\gamma)^2}.
\end{aligned}
\tag{kernelparts}
$$

If also $\gamma=\tau_{\mathrm{ion}}^{-1}$, the integrand in equation [eq:equalkernel] is simply $1+u/\tau_{\Phi,0}$, giving

$$
F_\alpha(t)=e^{-t/\tau_{\mathrm{ion}}}
\left[1+\frac{t}{\tau_{\mathrm{ion}}}
+\frac{t^2}{2\tau_{\mathrm{ion}}\tau_{\Phi,0}}\right].
\tag{kerneldegenerate}
$$

For smooth SFR evolution, rearranging the ODE and substituting its leading approximation on the derivative side gives

$$
\overline\Sigma_{\mathrm{SFR},\alpha}(t)
\simeq\Sigma_{\mathrm{SFR}}(t)
-\tau_{\mathrm{ion}}\frac{d\Sigma_{\mathrm{SFR}}(t)}{dt}.
\tag{shortkernel}
$$

The leading fractional correction is $-\tau_{\mathrm{ion}}\,d\ln\Sigma_{\mathrm{SFR}}/dt$. It is small when all relevant variation timescales exceed the response time and the prehistory has relaxed. A small first derivative at a single instant is insufficient if higher derivatives or an unresolved jump are large. For the present coefficients, $\gamma_{\mathrm{HI}}\tau_{\mathrm{ion}}=0.0122$ and $\gamma_{\mathrm{H_2}}\tau_{\mathrm{ion}}=0.0009$, and the exact calculation verifies the small one-Gyr correction.

# Appendix B. Variable gas coefficients and sensitivity tests

## B.1 Exact quadrature when the coefficients evolve

''')
sections.append(part('Let $\\Sigma_{\\mathrm{in}}','## A.2 Evolving') if 'Let $\\Sigma_{\\mathrm{in}}' in OLD else part('## A.1 Time-dependent','## A.2 Evolving').split('\n',1)[1])
sections.append(r'''
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
\tag{27}
$$

For piecewise-constant coefficients, solve each interval with its own constants and carry the terminal gas columns into the next interval as initial values. Gas masses remain continuous unless an explicit transport or removal impulse is imposed. If $\tau_{\mathrm{dep}}$ jumps, instantaneous SFR can jump even at continuous gas mass; Halpha must still use the stellar response in Appendix A.

## B.2 Faster conversion as a mathematical test only

Shortening $\tau_{\mathrm{conv}}$ is not the adopted explanation of RPS compression. Its sign of change is not established by Model 0 or by the cited observations. To retain the previous algebra as a sensitivity test, if initially balanced columns are held fixed while $\tau_{\mathrm{conv}}$ is reduced, the new fractional initial slope is

$$
\frac{1}{\Sigma_{\mathrm{SFR},0}}\left.\frac{d\Sigma_{\mathrm{SFR}}(t)}{dt}\right|_{0^+}
=\gamma_{\mathrm{H_2}}
\left(\frac{\tau_{\mathrm{conv,pre}}}{\tau_{\mathrm{conv,post}}}-1\right).
\tag{26}
$$

This follows by substituting the pre-change balance into equation [eq:17]. The labels pre/post in this equation refer to an imposed coefficient change, not the observational infall categories. Reducing the example conversion time by four gives a peak at 0.18368 Gyr with $F_{\mathrm{SFR}}=1.06458$, only 0.0272 dex. This is an alternative mathematical experiment, not one of the main Figure 1 curves.

Adding the two gas balances and integrating also gives a bound at fixed depletion time:

$$
\Sigma_{\mathrm{H_2}}(t)\leq\Sigma_{\mathrm{H_2},0}+\Sigma_{\mathrm{HI},0}
\quad\Longrightarrow\quad
F_{\mathrm{SFR}}(t)\leq1+\frac{\Sigma_{\mathrm{HI},0}}{\Sigma_{\mathrm{H_2},0}}.
\tag{28}
$$

For the illustrative columns this bound is 1.2853, or 0.1090 dex. It limits conversion of the existing local gas, not compression supplied by mass convergence, and not models with changing efficiency. It must not be used to rule out all RPS-induced spatial enhancement.

# Appendix C. Separating compact HII and leaked-OB emission

''')
photon=part('Let $\\mathcal Q_{\\mathrm{OB}}','# 5. Deriving')
photon=photon.replace('derived in the next section','defined in section 5').replace('Appendix A states','Appendix E states')
sections.append(photon)
young_split=part('## 6.3 Different','## 6.4 The').split('\n',1)[1]
sections.append(young_split[:young_split.index('Inserting $')])
sections.append(r'''
A different young effective spectrum can arise if the partition changes or the emitting gas changes. These are additional terms, not a consequence of multiplying both young luminosities by the same fading factor. Compact and leaked emission need not have equal line ratios because radiation filtering, ionization parameter, and gas conditions can differ ([Belfiore et al. 2022](#ref-belfiore)). The main derivation needs only their fixed effective young spectrum.

# Appendix D. Unequal Balmer decrements and the actual Hbeta audit

''')
dec=part('Let $\\mathcal B_j','## 6.3 Different')
dec=dec.replace('The empirical [O III]/Hbeta test uses the actual exported Hbeta denominator.', 'The secondary numerical test in this appendix uses the actual exported Hbeta denominator.')
dec=dec.replace(' A single dust correction applied to a mixture does not necessarily recover the individually corrected source components.','')
sections.append(dec)
sections.append(r'''
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
\tag{continuity}
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
\tag{transport}
$$

Every field in this display depends on $(\boldsymbol{x},t)$; only here the arguments are suppressed to keep the spatial conservation equations readable. Any vertical removal represented by $\gamma_{\mathrm{strip}}$ is already included in that sink and must not be counted again as a boundary loss. The transport term has a minus sign on the right. It adds column where the divergence of the mass flux is negative. Velocity convergence alone is not identical to this condition, since $\boldsymbol\nabla\cdot(\Sigma_i\boldsymbol v_i)=\boldsymbol v_i\cdot\boldsymbol\nabla\Sigma_i+\Sigma_i\boldsymbol\nabla\cdot\boldsymbol v_i$.

Integrating the transport term over a patch converts it into minus the outward mass flux through the boundary. A net inward flux can establish the initial columns in section 3.5 without shortening $\tau_{\mathrm{conv}}$ or $\tau_{\mathrm{dep}}$. Continuing transport requires a specified velocity field or another closure, so the simple closed local ODE is no longer sufficient. Model 0 neglects that later transport; it does not claim that RPS has no compression phase.

## E.2 Nonlocal photons and shared gas

''')
sections.append(part('A local SFR cannot','# Appendix B. Complete'))
sections.append(r'''
# Appendix F. Recombination details and differential fading

## F.1 Recovering the HOLMES normalization from the emission measure

''')
volume=part('In photoionization equilibrium','## 5.3 When')
# Main equation 38 already provides the compact result; omit its duplicate.
a=volume.index('$$\n\\begin{aligned}')
b=volume.index('This is the verified physical content')
volume=volume[:a]+r'''Equations [eq:36] and [eq:37] contain the same volume integral. Solving the former for that integral, substituting into the latter, and dividing by $A_{\mathrm{reg}}$ gives equation [eq:38]. The area must have the same units on both sides before converting to kpc squared.

'''
volume=volume.replace('Solving equation [eq:36] for the volume integral, inserting it into equation [eq:37], and dividing by $A_{\\mathrm{reg}}$ gives','')
sections.append(volume)
sections.append(part('Taking the logarithmic derivative','# 6. The three'))
sections.append('## F.2 Fixed spectra: luminosity can fall while the ratio rises\n\n')
fade=part('The total forbidden-line luminosity','# 7. The complete')
fade=fade.replace('F_\\alpha','F_{\\mathrm{SFR}}')
sections.append(fade)
sections.append('## F.3 Evolving populations, absorption, and intrinsic spectra\n\n')
sections.append(part('Without setting the old-star','## A.3 Nonlocal'))

glossary=part('# Appendix B. Complete','# Appendix C. Provenance')
glossary=glossary.replace('Appendix B.','Appendix G.').replace('Table B','Table G')
glossary=glossary.replace('generalized here to three','generalized here to two effective')
glossary=glossary.replace(r'w_{\mathrm{target},B}',r'w_{\mathrm{target}}')
glossary+=r'''
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

'''
sections.append(glossary)
sections.append(r'''# Appendix H. Provenance, revision coverage, and verification

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

''')
checks=AUDIT['checks']
sections.append(f'''The executed gas-versus-ODE checks have maximum relative difference {max(checks[k]['gas_ode_max_rel'] for k in ['balanced','closed_no_RPS','compression']):.3g}; the filtered-response check has maximum absolute normalized difference {max(checks[k]['filtered_ode_max_abs'] for k in ['balanced','closed_no_RPS','compression']):.3g}. The equal-gas-rate check differs by {checks['equal_rate_max_abs']:.3g}; the scaled-initial-state ODE check differs by {checks['spatial_initial_state_max_abs']:.3g}. Peak, positive-supply, total-mass, and unequal-Balmer-weight identities pass for the implemented examples. These numerical tolerances test the equations, not the physical assumptions.

''')
sections.append(r'''The observational input is `assets/20260914_resolved_RPS_academic_model/stage_bpt_line_profiles_with_hbeta.csv`. Its fingerprint and the six recorded source fingerprints are checked and saved in the new asset directory. This validates reuse of that extraction record, not the present contents of every large FITS map. No full map pipeline, bootstrap, fitted stellar population, new CO/HI analysis, radiation transport calculation, or photoionization grid was executed. No significance is assigned to the cross-sectional amplitude tests. Figure colors identify model cases, not galaxies.

The PDF is generated from the same Markdown and checked for equation numbering, internal links, text boundaries, and page rendering. Detailed acceptance evidence and the revision checklist are saved with the assets. The source report and other user files are preserved.

''')
refs=OLD[OLD.index('# References'):]
refs+=r'''
<span id="ref-brown"></span>
**Brown, T., et al. (2023).** *VERTICO VII: Environmental Quenching Caused by Suppression of Molecular Gas Content and Star Formation Efficiency in Virgo Cluster Galaxies.* ApJ, 956, 37. [DOI](https://doi.org/10.3847/1538-4357/acf195); [primary manuscript](https://arxiv.org/pdf/2308.10943).

<span id="ref-armitage"></span>
**Armitage, P. J. (2022).** *Lecture notes on accretion disk physics.* arXiv:2201.07262, sections II.A and III.A.1. [Primary manuscript](https://arxiv.org/pdf/2201.07262). Equations 9 and 91--97 define surface density and derive its continuity equation; the phase-specific sources and sinks are added explicitly in this report.
'''
ref_blocks=re.findall(r'<span id="ref-[\s\S]*?(?=<span id="ref-|\Z)',refs)
ref_blocks.sort(key=lambda block:re.search(r'\*\*([^*]+)',block)[1])
sections.append('# References\n\n'+'\n'.join(ref_blocks))
text='\n'.join(sections)
text=text.replace('assets/20261001_HI_stripping_HOLMES','assets/20261002_Model0_Derivation')
text=text.replace('Appendix A states the corresponding limitation','Appendix E states the corresponding limitation')
labels=re.findall(r'\\tag\{([^}]+)\}',text)
assert len(labels)==len(set(labels)), [x for x in labels if labels.count(x)>1]
mapping={label:str(i+1) for i,label in enumerate(labels)}
text=re.sub(r'\\tag\{([^}]+)\}',lambda m:r'\tag{'+mapping[m[1]]+'}',text)
text=re.sub(r'\[eq:([^\]]+)\]',lambda m:'('+mapping[m[1]]+')',text)
assert '[eq:' not in text
(ROOT/'20261002_Model0_Gas_Fading_and_HOLMES_Derivation.md').write_text(text)
(OUT/'equation_label_map.json').write_text(json.dumps(mapping,indent=2)+'\n')
print(json.dumps(dict(equations=len(labels),words=len(text.split()),characters=len(text))))
