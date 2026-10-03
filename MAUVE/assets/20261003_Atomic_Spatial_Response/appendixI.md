# Appendix I. Temporal behaviour and physical limits of the atomic spatial model

## I.1 Equal rates and the condition for an initial SFR increase

Section 3.5 compares positions at a common elapsed time. Its spatial ordering does not by itself fix the temporal derivative at any position. For completeness, if $\gamma_{\mathrm{HI}}=\gamma_{\mathrm{H_2}}\equiv\gamma$, the exponential integrand in section 3.2 is unity. Its time integral is $t$, giving the finite solutions

$$
\begin{aligned}
\Sigma_{\mathrm{H_2}}(x,t)
&=\left[\Sigma_{\mathrm{H_2},0}
+\frac{c_{\mathrm{HI}}(x)\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}t\right]e^{-\gamma t},\\
\Sigma_{\mathrm{SFR}}(x,t)
&=\left[\frac{\Sigma_{\mathrm{H_2},0}}{\tau_{\mathrm{dep}}}
+\frac{c_{\mathrm{HI}}(x)\Sigma_{\mathrm{HI},0}}
{\tau_{\mathrm{conv}}\tau_{\mathrm{dep}}}t\right]e^{-\gamma t}.
\end{aligned}
\tag{81}
$$

For either rate ordering, the initial replenishment time is the common molecular column divided by the local initial supply. With $\tau_{\Phi,0}$ denoting the reference value and $\Sigma_{\mathrm{SFR},0}=\Sigma_{\mathrm{H_2},0}/\tau_{\mathrm{dep}}$, the molecular balance gives

$$
\begin{aligned}
\tau_\Phi(x,0)
&=\frac{\Sigma_{\mathrm{H_2},0}}
{c_{\mathrm{HI}}(x)\Sigma_{\mathrm{HI},0}/\tau_{\mathrm{conv}}}
=\frac{\tau_{\Phi,0}}{c_{\mathrm{HI}}(x)},\\
\left.\frac{d\Sigma_{\mathrm{SFR}}(x,t)}{dt}\right|_{t=0}
&=\frac{1}{\tau_{\mathrm{dep}}}
\left[\frac{c_{\mathrm{HI}}(x)\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
-\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0}\right]\\
&=\Sigma_{\mathrm{SFR},0}
\left[\frac{c_{\mathrm{HI}}(x)}{\tau_{\Phi,0}}-\gamma_{\mathrm{H_2}}\right].
\end{aligned}
\tag{82}
$$

An atomic excess shortens the replenishment time through a larger available supply, while $\tau_{\mathrm{conv}}$ is unchanged. Initial temporal growth requires $c_{\mathrm{HI}}(x)/\tau_{\Phi,0}>\gamma_{\mathrm{H_2}}$. The weaker condition $c_{\mathrm{HI}}(x)>1$ guarantees only a spatial excess over the matched reference at $t>0$.

If the reference is initially balanced, the result simplifies to

$$
\begin{aligned}
\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
&=\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2},0},
\qquad \tau_{\Phi,0}^{-1}=\gamma_{\mathrm{H_2}},\\
\left.\frac{d\Sigma_{\mathrm{SFR}}(x,t)}{dt}\right|_{t=0}
&=[c_{\mathrm{HI}}(x)-1]\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{SFR},0}.
\end{aligned}
\tag{83}
$$

The leading case then initially rises, the reference has zero initial slope followed by decline, and the trailing case initially declines. Section 3.4 proves that any initially rising solution has one maximum and subsequently declines. Its peak time follows from equation (24) after replacing $\tau_{\Phi,0}$ by $\tau_{\Phi,0}/c_{\mathrm{HI}}(x)$; the same substitution applies to the equal-rate limit given there.

With positive response rates and finite positive atomic factors, every gas and SFR solution approaches zero at late times. All regions eventually decline, while the spatial ordering in equation (30) persists at every finite $t>0$. Absolute SFR differences approach zero, although their ratios need not approach unity. Thus a leading-side spatial excess can remain after that region has passed its temporal SFR maximum.

## I.2 Physical interpretation and an additional molecular response

The atomic factors specify local initial conditions, not a prediction of the wind pressure or how the asymmetry formed. A higher surface density requires accumulation in the adopted patch or a reduction of its physical area; vertical compression alone need not increase the mass per disc area. The transport contribution $-\boldsymbol{\nabla}\cdot(\Sigma_i\boldsymbol v_i)$ in Appendix E can establish local column changes. Model 0 begins after that unresolved episode. It neither creates mass during the subsequent closed evolution nor enforces a global redistribution budget. Common initial H2 is an idealization isolating the atomic pathway, not a necessary outcome of a finite compression episode.

The density dependence of ram-pressure coupling motivates this choice of perturbed phase. A lower mass per exposed area is more easily accelerated at fixed incident momentum flux, while gravity, shielding, and geometry also matter. The relevant column is along the incident wind and can differ from an observed projected column ([Cramer et al. 2020, section 6.1](#ref-cramer)). Holding $\gamma_{\mathrm{strip}}$ fixed across the comparison isolates the initial-column effect; it does not solve for the possible column dependence of that coefficient. Likewise, a leading atomic excess and trailing deficit are testable assumptions, not a universal consequence of wind geometry. Downstream transport can instead accumulate gas in a trailing region.

Molecular gas is not immune. [Lee et al. (2017), sections 5.2--5.3](#ref-lee), find disturbed CO morphology and kinematics and upstream enhancements in NGC4330, NGC4402, and NGC4522, despite no clear molecular-stripping signature at their sensitivity. [Cramer et al. (2020), sections 5.2 and 6.1--6.3](#ref-cramer), identify compressed and displaced molecular gas in NGC4402. Their interpretation distinguishes initially resistant dense clouds from diffuse molecular gas, while allowing cloud disruption and effective removal; the phase history of stripped material is uncertain. These observations support differential coupling, not the exact initial conditions imposed here.

An additional initial molecular response can be represented by $c_{\mathrm{H_2}}(x)\equiv\Sigma_{\mathrm{H_2}}(x,0)/\Sigma_{\mathrm{H_2},0}$. Keeping all subsequent coefficients fixed, subtracting the atomic-only solution from this extended solution leaves only the change in the surviving initial molecular term:

$$
\begin{aligned}
\Sigma_{\mathrm{H_2}}^{\mathrm{extended}}(x,t)
-\Sigma_{\mathrm{H_2}}^{\mathrm{atomic\ only}}(x,t)
&=[c_{\mathrm{H_2}}(x)-1]\Sigma_{\mathrm{H_2},0}e^{-\gamma_{\mathrm{H_2}}t},\\
\Sigma_{\mathrm{SFR}}^{\mathrm{extended}}(x,t)
-\Sigma_{\mathrm{SFR}}^{\mathrm{atomic\ only}}(x,t)
&=\frac{[c_{\mathrm{H_2}}(x)-1]\Sigma_{\mathrm{H_2},0}}
{\tau_{\mathrm{dep}}}e^{-\gamma_{\mathrm{H_2}}t}.
\end{aligned}
\tag{84}
$$

The main calculation and Figure 1 set $c_{\mathrm{H_2}}(x)=1$. An additional molecular excess on the leading side or deficit on the trailing side would amplify the corresponding SFR offset. Amplification requires the molecular response to have the same spatial sign as the atomic response. Displacement, changing depletion times, or continuing transport need not do so and require a different extension.

[Brown et al. (2023), section 3.3 and Figure 5](#ref-brown), find enhanced outer-disc SFR associated with larger molecular columns at fixed stellar density in their four early-stage RPS galaxies, with molecular SFE consistent with the field. This motivates the fixed-efficiency comparison but does not distinguish atomic conversion from molecular transport or establish constant efficiency for every MAUVE region. The potential NGC4654 gradient likewise motivates a test; the factors in section 8.3 are illustrative and are not fitted to that galaxy.

