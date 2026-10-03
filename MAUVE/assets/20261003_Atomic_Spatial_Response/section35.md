## 3.5 Atomic-gas asymmetry and the leading/trailing SFR contrast

We now compare positions $x$ with different initial atomic columns, but the same initial molecular column $\Sigma_{\mathrm{H_2},0}$ and the same $\tau_{\mathrm{conv}}$, $\tau_{\mathrm{dep}}$, $R$, $\lambda$, and $\gamma_{\mathrm{strip}}$. Here $x$ denotes a spatial position, written without boldface for this local comparison. The rates in equation (5) are therefore common to all positions. We define $c_{\mathrm{HI}}(x)\equiv\Sigma_{\mathrm{HI}}(x,0)/\Sigma_{\mathrm{HI},0}>0$, where $\Sigma_{\mathrm{HI},0}$ is the reference initial atomic column and $c_{\mathrm{HI}}(x_{\mathrm{reference}})=1$.

This factor represents an atomic excess or deficit established before model onset. Applying the perturbation only to HI is motivated by the greater susceptibility of diffuse gas than dense molecular clouds to ram pressure ([Boselli et al. 2014](#ref-boselli); [Lee et al. 2017](#ref-lee); [Cramer et al. 2020](#ref-cramer)). It is an initial-condition approximation; the physical qualifications are discussed in Appendix I. Equation (8) and the unchanged supply law immediately give

$$
\begin{aligned}
\Sigma_{\mathrm{HI}}(x,t)
&=c_{\mathrm{HI}}(x)\Sigma_{\mathrm{HI},0}e^{-\gamma_{\mathrm{HI}}t},\\
\Sigma_\Phi(x,t)
&=\frac{c_{\mathrm{HI}}(x)\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\end{aligned}
\tag{26}
$$

The molecular equation remains $d\Sigma_{\mathrm{H_2}}(x,t)/dt=\Sigma_\Phi(x,t)-\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2}}(x,t)$. Relative to section 3.2, only the supply amplitude is multiplied by the constant $c_{\mathrm{HI}}(x)$. The integrating-factor solution in equation (13) therefore becomes

$$
\begin{aligned}
\Sigma_{\mathrm{H_2}}(x,t)
&=\Sigma_{\mathrm{H_2},0}e^{-\gamma_{\mathrm{H_2}}t}\\
&\quad+\frac{c_{\mathrm{HI}}(x)\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}
\tag{27}
$$

Only the supplied term is multiplied by $c_{\mathrm{HI}}(x)$; the initial molecular term is common. Dividing by the common depletion time gives

$$
\begin{aligned}
\Sigma_{\mathrm{SFR}}(x,t)
&=\frac{\Sigma_{\mathrm{H_2},0}}{\tau_{\mathrm{dep}}}e^{-\gamma_{\mathrm{H_2}}t}\\
&\quad+\frac{c_{\mathrm{HI}}(x)\Sigma_{\mathrm{HI},0}}
{\tau_{\mathrm{conv}}\tau_{\mathrm{dep}}}
\frac{e^{-\gamma_{\mathrm{HI}}t}-e^{-\gamma_{\mathrm{H_2}}t}}
{\gamma_{\mathrm{H_2}}-\gamma_{\mathrm{HI}}}.
\end{aligned}
\tag{28}
$$

These expressions use consistent physical units, with the numerical conversion in equation (6). For equal response rates, replace the exponential quotient by $t e^{-\gamma_{\mathrm{H_2}}t}$, as in equation (15).

Let $\Sigma_{\mathrm{H_2}}^{\mathrm{reference}}(t)$ be equation (27) evaluated at $c_{\mathrm{HI}}=1$. Subtracting the solutions at any two positions $x_1$ and $x_2$ cancels the common initial term:

$$
\begin{aligned}
&\Sigma_{\mathrm{H_2}}(x_1,t)-\Sigma_{\mathrm{H_2}}(x_2,t)\\
&\quad=[c_{\mathrm{HI}}(x_1)-c_{\mathrm{HI}}(x_2)]
\left[\Sigma_{\mathrm{H_2}}^{\mathrm{reference}}(t)
-\Sigma_{\mathrm{H_2},0}e^{-\gamma_{\mathrm{H_2}}t}\right].
\end{aligned}
\tag{29}
$$

The second bracket is the surviving molecular gas supplied after $t=0$. It is strictly positive for every finite $t>0$, since the numerator and denominator of the exponential quotient have the same sign; its equal-rate limit is also positive. It is not the full reference molecular column: that would incorrectly predict a nonzero difference at $t=0$. Consequently, the molecular difference, and the SFR difference obtained by dividing by $\tau_{\mathrm{dep}}>0$, have the sign of $c_{\mathrm{HI}}(x_1)-c_{\mathrm{HI}}(x_2)$.

If the proposed geometry establishes $c_{\mathrm{HI}}(x_{\mathrm{leading}})>1>c_{\mathrm{HI}}(x_{\mathrm{trailing}})>0$, then

$$
\boxed{
\begin{aligned}
\Sigma_{\mathrm{H_2}}^{\mathrm{leading}}(t)
&>\Sigma_{\mathrm{H_2}}^{\mathrm{reference}}(t)
>\Sigma_{\mathrm{H_2}}^{\mathrm{trailing}}(t),\\
\Sigma_{\mathrm{SFR}}^{\mathrm{leading}}(t)
&>\Sigma_{\mathrm{SFR}}^{\mathrm{reference}}(t)
>\Sigma_{\mathrm{SFR}}^{\mathrm{trailing}}(t),\qquad t>0.
\end{aligned}}
\tag{30}
$$

The superscripts denote evaluation at the corresponding position. All three SFRs are equal at $t=0$; their differences develop through the atomic supply. The reference has the same stripping coefficient as the other positions. This is a spatial ordering, not a statement that the leading SFR is increasing with time or that every RPS geometry has the assumed atomic pattern. Appendix I develops those qualifications; section 8.3 evaluates this atomic-only family numerically.

