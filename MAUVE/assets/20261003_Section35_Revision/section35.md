## 3.5 Atomic-gas asymmetry and the leading/trailing SFR contrast

Sections 3.1--3.4 describe the temporal evolution of one region. Here we compare regions whose initial atomic columns differ, while their initial molecular columns and subsequent evolution coefficients are identical. The molecular and SFR differences are then derived from the different atomic supply. The potential NGC4654 gradient motivates this comparison; it is not a fitted constraint or a newly established detection in this report.

### 3.5.1 Physical motivation and the comparison being made

Ram pressure does not couple equally to every gas structure. At a given incident momentum flux, a lower gas mass per exposed area is more easily accelerated; gravitational restoring forces, geometry, and shielding also matter. The relevant column is along the incident wind, which need not equal the observed projected column ([Cramer et al. 2020, section 6.1](https://doi.org/10.3847/1538-4357/abaf54)). Thus the relevant distinction is not strictly HI versus H2, nor volume density alone. Diffuse atomic gas is usually more susceptible to removal and displacement than high-column molecular clouds. The finding that molecular gas is affected less efficiently than HI in Virgo provides motivation for the HI-only baseline ([Boselli et al. 2014](#ref-boselli)). It does not establish that the molecular response is exactly zero.

We therefore impose the spatial perturbation only on the initial atomic reservoir, consistently with placing the direct stripping coefficient $\gamma_{\mathrm{strip}}$ only in the atomic balance. This is a deliberately restrictive approximation. The initial perturbation represents a local atomic excess or deficit already established by compression, displacement, or transport. Subsequent evolution is the closed Model 0 system. We do not add a continuous atomic source, change $\tau_{\mathrm{conv}}$, or prescribe an initial molecular enhancement.

Let $j$ label a comparison region: leading, reference, or trailing. All regions use the same area convention, elapsed time $t$, and constant $\tau_{\mathrm{conv}}$, $\tau_{\mathrm{dep}}$, $R$, $\lambda$, and $\gamma_{\mathrm{strip}}$. Consequently, they share the rates $\gamma_{\mathrm{HI}}$ and $\gamma_{\mathrm{H_2}}$ defined in equation (5). Let $\Sigma_{\mathrm{HI},0}$ and $\Sigma_{\mathrm{H_2},0}$ denote the positive initial columns of the reference. We define the dimensionless atomic-column factor and impose the common initial molecular column:

$$
\begin{aligned}
c_{\mathrm{HI},j}&\equiv
\frac{\Sigma_{\mathrm{HI},j}(0)}{\Sigma_{\mathrm{HI},0}}>0,
\qquad c_{\mathrm{HI,reference}}=1,\\
\Sigma_{\mathrm{H_2},j}(0)&=\Sigma_{\mathrm{H_2},0},
\qquad \Sigma_{\mathrm{SFR},j}(0)
=\frac{\Sigma_{\mathrm{H_2},0}}{\tau_{\mathrm{dep}}}
\equiv\Sigma_{\mathrm{SFR},0}.
\end{aligned}
\tag{26}
$$

An atomic excess has $c_{\mathrm{HI},j}>1$; an atomic deficit has $0<c_{\mathrm{HI},j}<1$. The reference $c_{\mathrm{HI},j}=1$ still experiences the same stripping coefficient as the other regions. It is a matched spatial comparison, not the separate no-stripping comparison in section 8.3. These assumptions isolate the effect of the initial atomic column; they are not asserted to hold between arbitrary observed patches.

### 3.5.2 From the atomic column to the molecular and SFR solutions

The first balance is unchanged: $d\Sigma_{\mathrm{HI},j}/dt=-\gamma_{\mathrm{HI}}\Sigma_{\mathrm{HI},j}$. Applying equation (8) to the initial condition in equation (26), and then using the supply closure in equation (1), gives

$$
\begin{aligned}
\Sigma_{\mathrm{HI},j}(t)
&=c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}e^{-\gamma_{\mathrm{HI}}t},\\
\Sigma_{\Phi,j}(t)
&=\frac{\Sigma_{\mathrm{HI},j}(t)}{\tau_{\mathrm{conv}}}
=\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\end{aligned}
\tag{27}
$$

Thus a larger atomic column produces a larger absolute molecular supply rate at the same conversion time. The factor multiplies the initial condition and the ensuing supply; it is not an extra mass source. Substituting this supply into the molecular balance yields

$$
\frac{d\Sigma_{\mathrm{H_2},j}(t)}{dt}
+\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},j}(t)
=\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\tag{28}
$$

Multiply both sides by $e^{\gamma_{\mathrm{H_2}}t}$. The product rule combines the two terms on the left:

$$
\begin{aligned}
\frac{d}{dt}\left[e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2},j}(t)\right]
&=e^{\gamma_{\mathrm{H_2}}t}\frac{d\Sigma_{\mathrm{H_2},j}(t)}{dt}
+\gamma_{\mathrm{H_2}}e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2},j}(t)\\
&=\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})t}.
\end{aligned}
\tag{29}
$$

Integrate from $0$ to $t$, with $u$ the integration time. The lower limit is the common initial molecular column, since $e^0=1$. For unequal rates this gives

$$
\begin{aligned}
e^{\gamma_{\mathrm{H_2}}t}\Sigma_{\mathrm{H_2},j}(t)-\Sigma_{\mathrm{H_2},0}
&=\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\int_0^t e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})u}\,du\\
&=\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\left[\frac{e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})u}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}\right]_0^t\\
&=\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\frac{e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})t}-1}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}
\tag{30}
$$

Add $\Sigma_{\mathrm{H_2},0}$ to both sides and multiply by $e^{-\gamma_{\mathrm{H_2}}t}$. In the supply term, $e^{-\gamma_{\mathrm{H_2}}t}e^{(\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}})t}=e^{-\gamma_{\mathrm{HI}}t}$. Therefore

$$
\boxed{
\begin{aligned}
\Sigma_{\mathrm{H_2},j}(t)
&=\Sigma_{\mathrm{H_2},0}e^{-\gamma_{\mathrm{H_2}}t}\\
&\quad+\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}}
\tag{31}
$$

Dividing each term by the common molecular depletion time gives the explicit SFR solution:

$$
\boxed{
\begin{aligned}
\Sigma_{\mathrm{SFR},j}(t)
&=\frac{\Sigma_{\mathrm{H_2},0}}{\tau_{\mathrm{dep}}}e^{-\gamma_{\mathrm{H_2}}t}\\
&\quad+\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}
{\tau_{\mathrm{conv}}\tau_{\mathrm{dep}}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}}
\tag{32}
$$

Equation (32) uses consistent physical units; equation (6) supplies the conversion when reporting SFR in $M_\odot\,\mathrm{yr}^{-1}\,\mathrm{kpc}^{-2}$. Setting $c_{\mathrm{HI},j}=1$ recovers equations (13)--(14). The result follows from this report's linear transfer assumption and the regulator bookkeeping described in section 2 ([Lilly et al. 2013](#ref-lilly); [Huang et al. 2026](#ref-huang)); the spatial family and its solution are our derivation, not an empirical law taken from those papers.

If $\gamma_{\mathrm{HI}}=\gamma_{\mathrm{H_2}}\equiv\gamma>0$, the integrand in equation (30) is one, so its integral is $t$. Repeating the last multiplication, rather than dividing by zero, gives

$$
\begin{aligned}
\Sigma_{\mathrm{H_2},j}(t)
&=\left[\Sigma_{\mathrm{H_2},0}
+\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}t\right]e^{-\gamma t},\\
\Sigma_{\mathrm{SFR},j}(t)
&=\left[\frac{\Sigma_{\mathrm{H_2},0}}{\tau_{\mathrm{dep}}}
+\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}
{\tau_{\mathrm{conv}}\tau_{\mathrm{dep}}}t\right]e^{-\gamma t}.
\end{aligned}
\tag{33}
$$

### 3.5.3 Why the atomic ordering produces an SFR ordering

To expose the sign without relying on the ordering of the two rates, separate the surviving initial molecular column from the molecular column supplied after $t=0$ in the reference:

$$
\begin{aligned}
\Sigma_{\mathrm{H_2}}^{\mathrm{init}}(t)
&\equiv\Sigma_{\mathrm{H_2},0}e^{-\gamma_{\mathrm{H_2}}t},\\
\Sigma_{\mathrm{H_2}}^{\mathrm{sup,ref}}(t)
&\equiv\int_0^t\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}u}e^{-\gamma_{\mathrm{H_2}}(t-u)}\,du,\\
\Sigma_{\mathrm{H_2},j}(t)
&=\Sigma_{\mathrm{H_2}}^{\mathrm{init}}(t)
+c_{\mathrm{HI},j}\Sigma_{\mathrm{H_2}}^{\mathrm{sup,ref}}(t).
\end{aligned}
\tag{34}
$$

Both defined quantities are surface densities with the same mass and area convention as $\Sigma_{\mathrm{H_2}}$. The second is the portion of the supplied molecular gas still present at time $t$, not the total mass ever converted: its exponential factor accounts for consumption between supply time $u$ and observation time $t$. These are contributions to the existing molecular reservoir, not extra reservoirs or free parameters. Every factor in its integrand is positive. Hence $\Sigma_{\mathrm{H_2}}^{\mathrm{sup,ref}}(t)>0$ for every finite $t>0$, whether the rates are unequal or equal.

For any two regions $j$ and $k$, the common initial-molecular contribution cancels on subtraction:

$$
\begin{aligned}
\Sigma_{\mathrm{H_2},j}(t)-\Sigma_{\mathrm{H_2},k}(t)
&=(c_{\mathrm{HI},j}-c_{\mathrm{HI},k})
\Sigma_{\mathrm{H_2}}^{\mathrm{sup,ref}}(t),\\
\Sigma_{\mathrm{SFR},j}(t)-\Sigma_{\mathrm{SFR},k}(t)
&=\frac{c_{\mathrm{HI},j}-c_{\mathrm{HI},k}}{\tau_{\mathrm{dep}}}
\Sigma_{\mathrm{H_2}}^{\mathrm{sup,ref}}(t).
\end{aligned}
\tag{35}
$$

The SFR difference therefore has exactly the sign of the atomic-column-factor difference. If the proposed geometry establishes $c_{\mathrm{HI,leading}}>1>c_{\mathrm{HI,trailing}}>0$, then

$$
\boxed{
\begin{aligned}
\Sigma_{\mathrm{H_2,leading}}(t)
&>\Sigma_{\mathrm{H_2,reference}}(t)
>\Sigma_{\mathrm{H_2,trailing}}(t),\\
\Sigma_{\mathrm{SFR,leading}}(t)
&>\Sigma_{\mathrm{SFR,reference}}(t)
>\Sigma_{\mathrm{SFR,trailing}}(t),\qquad t>0.
\end{aligned}}
\tag{36}
$$

At $t=0$ the molecular columns and SFRs are equal by construction. Their differences grow continuously from zero through the changed supply; there is no imposed instantaneous molecular or SFR jump. The strict ordering at later times follows from the chosen atomic ordering and the matched coefficients. It is not a proof that every RPS galaxy must have that atomic geometry. Gas transported downstream can instead accumulate in a trailing region.

The corresponding spatial offset from the contemporaneous reference, in dex, is

$$
\begin{aligned}
\Delta_{\mathrm{spatial},j}\log_{10}\Sigma_{\mathrm{SFR}}(t)
&\equiv\log_{10}\frac{\Sigma_{\mathrm{SFR},j}(t)}
{\Sigma_{\mathrm{SFR,reference}}(t)}\\
&=\log_{10}\left[
\frac{\Sigma_{\mathrm{H_2}}^{\mathrm{init}}(t)
+c_{\mathrm{HI},j}\Sigma_{\mathrm{H_2}}^{\mathrm{sup,ref}}(t)}
{\Sigma_{\mathrm{H_2}}^{\mathrm{init}}(t)
+\Sigma_{\mathrm{H_2}}^{\mathrm{sup,ref}}(t)}\right].
\end{aligned}
\tag{37}
$$

For $t>0$ it is positive for an atomic excess and negative for an atomic deficit. A facing-to-opposite contrast alone can combine an excess on one side with a deficit on the other; it does not determine either offset relative to an independently specified reference.

### 3.5.4 Spatial enhancement versus temporal increase and decline

The initial replenishment time of region $j$, $\tau_{\Phi,j,0}$, follows by dividing the common initial molecular column by its initial supply. Substitution into the molecular balance also gives its initial SFR slope:

$$
\begin{aligned}
\tau_{\Phi,j,0}
&=\frac{\Sigma_{\mathrm{H_2},0}}
{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}/\tau_{\mathrm{conv}}}
=\frac{\tau_{\Phi,0}}{c_{\mathrm{HI},j}},\\
\left.\frac{d\Sigma_{\mathrm{SFR},j}}{dt}\right|_0
&=\frac{1}{\tau_{\mathrm{dep}}}
\left[\frac{c_{\mathrm{HI},j}\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
-\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0}\right]\\
&=\Sigma_{\mathrm{SFR},0}
\left[\frac{c_{\mathrm{HI},j}}{\tau_{\Phi,0}}-\gamma_{\mathrm{H_2}}\right].
\end{aligned}
\tag{38}
$$

Here $\tau_{\Phi,0}$ is the reference replenishment time already defined in section 3.3. An atomic excess shortens $\tau_{\Phi,j,0}$ by increasing the available supply, while $\tau_{\mathrm{conv}}$ itself remains unchanged. Temporal growth requires $c_{\mathrm{HI},j}/\tau_{\Phi,0}>\gamma_{\mathrm{H_2}}$. The weaker condition $c_{\mathrm{HI},j}>1$ guarantees a spatial excess at $t>0$, but does not guarantee a positive temporal derivative when the reference is initially supply-starved.

For the additional special case of an initially balanced reference, the relation between the two statements becomes explicit:

$$
\begin{aligned}
\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
&=\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0},
\qquad \tau_{\Phi,0}^{-1}=\gamma_{\mathrm{H_2}},\\
\left.\frac{d\Sigma_{\mathrm{SFR},j}}{dt}\right|_0
&=(c_{\mathrm{HI},j}-1)\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{SFR},0}.
\end{aligned}
\tag{39}
$$

The leading case then initially rises, the reference has zero initial slope followed by decline, and the trailing case initially declines. Section 3.4 proves that an initially rising solution has one maximum and subsequently declines. Its peak time is obtained from equation (24) by replacing $\tau_{\Phi,0}$ with $\tau_{\Phi,0}/c_{\mathrm{HI},j}$; the same substitution applies to the equal-rate limit stated there. The atomic reservoir and its supply fade, so they cannot maintain the initial rise indefinitely.

For every finite positive $c_{\mathrm{HI},j}$ and positive response rates, equations (31)--(33) approach zero as $t\rightarrow\infty$. All regions eventually decline, while the spatial ordering in equation (36) holds at every finite $t>0$. The absolute SFR differences also approach zero; their ratios need not approach unity. Thus a positive spatial offset can coexist with a negative temporal derivative, including after the leading region has passed its SFR maximum.

### 3.5.5 Physical limits and the relation to the retained numerical example

The atomic factors are imposed local initial conditions, not a prediction of the wind pressure or a calculation of how the asymmetry formed. A higher surface density requires mass accumulation in the adopted patch or a reduction of its physical area; vertical compression alone need not increase the mass per disc area. The transport contribution $-\boldsymbol{\nabla}\cdot(\Sigma_i\boldsymbol v_i)$ in Appendix E can establish local column changes. The present derivation starts after that unresolved episode and neither creates mass during the closed evolution nor enforces a global redistribution budget. Identical initial H2 is an idealization that isolates the atomic pathway, not a necessary result of a finite real compression episode.

[Lee et al. (2017), sections 5.2--5.3](https://doi.org/10.1093/mnras/stw3162), find disturbed CO morphology and kinematics and upstream enhancements in NGC4330, NGC4402, and NGC4522, despite no clear molecular-stripping signature at their sensitivity. [Cramer et al. (2020), sections 5.2 and 6.1--6.3](https://doi.org/10.3847/1538-4357/abaf54), identify compression and displaced molecular material in NGC4402. They interpret dense clouds as initially more resistant than diffuse molecular gas, while allowing subsequent cloud disruption and effective removal. The phase history of the stripped material remains uncertain. These observations support differential coupling and the possibility of a molecular response; they do not establish our particular initial conditions or a universal leading/trailing pattern.

Direct molecular responses therefore remain possible. Their implications can be stated without including them in the atomic-only baseline. If the initial molecular column also changes by a dimensionless factor $c_{\mathrm{H_2},j}$, while subsequent coefficients remain the same, linearity gives

$$
\begin{aligned}
\Sigma_{\mathrm{H_2},j}^{\mathrm{extended}}(t)
&=c_{\mathrm{H_2},j}\Sigma_{\mathrm{H_2}}^{\mathrm{init}}(t)
+c_{\mathrm{HI},j}\Sigma_{\mathrm{H_2}}^{\mathrm{sup,ref}}(t),\\
\Sigma_{\mathrm{SFR},j}^{\mathrm{extended}}(t)
-\Sigma_{\mathrm{SFR},j}^{\mathrm{atomic\ only}}(t)
&=\frac{c_{\mathrm{H_2},j}-1}{\tau_{\mathrm{dep}}}
\Sigma_{\mathrm{H_2}}^{\mathrm{init}}(t).
\end{aligned}
\tag{40}
$$

Here $c_{\mathrm{H_2},j}\equiv\Sigma_{\mathrm{H_2},j}(0)/\Sigma_{\mathrm{H_2},0}$ applies only to this extension; the atomic-only calculation has $c_{\mathrm{H_2},j}=1$. An additional initial molecular enhancement on the leading side, or deficit on the trailing side, therefore amplifies the corresponding SFR offset. Amplification is conditional on the molecular response having the same spatial sign as the atomic response. Molecular displacement, changes in depletion time, and continuing transport need not have that effect and are not described by this initial-condition correction.

The molecular-gas/SFR association is consistent with [Brown et al. (2023), section 3.3 and Figure 5](#ref-brown): their four early-stage RPS galaxies show enhanced outer-disc SFR associated with higher molecular surface density at fixed stellar density, with molecular SFE consistent with the field. This motivates a fixed-efficiency comparison but does not identify whether the molecular excess arose from atomic conversion, molecular transport, or another process. Nor does it establish constant efficiency for every MAUVE region.

The existing section 8.3 and Figure 1 retain a distinct illustrative initial condition in which both gas columns are multiplied by 1.5. That example is the extension in equation (40), not a numerical evaluation of the atomic-only family in equations (26)--(39). Its elevated SFR at $t=0$ is therefore consistent with, but different from, the equal initial SFR imposed here. The remainder of the report's numerical calculations and line-emission model are unchanged.

**References introduced in this section.** Lee, B., et al. (2017), *The effect of ram pressure on the molecular gas of galaxies: three case studies in the Virgo cluster*, MNRAS, 466, 1382--1398, [doi:10.1093/mnras/stw3162](https://doi.org/10.1093/mnras/stw3162), [primary manuscript](https://arxiv.org/pdf/1701.02750). Cramer, W. J., et al. (2020), *ALMA evidence for ram pressure compression and stripping of molecular gas in the Virgo cluster galaxy NGC 4402*, ApJ, 901, [doi:10.3847/1538-4357/abaf54](https://doi.org/10.3847/1538-4357/abaf54), [primary manuscript](https://arxiv.org/pdf/1910.14082).
