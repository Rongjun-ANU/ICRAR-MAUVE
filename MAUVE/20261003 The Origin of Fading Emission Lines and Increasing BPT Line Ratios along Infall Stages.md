# 20261003 The Origin of Fading Emission Lines and Increasing BPT Line Ratios along Infall Stages

## 1. Physical picture and disclaimer

We consider a local region in which RPS removes atomic gas but does not directly remove molecular gas. The loss of HI reduces the subsequent molecular supply. The retained H$_2$ reservoir then supports star formation while it is gradually consumed. Young-star H$\alpha$ emission follows this declining SFR. When old stars and absorbing gas remain, hot low-mass evolved stars (HOLMES) can provide a more slowly varying ionizing contribution, whose fraction of the total Balmer emission increases. A forbidden-line-to-Balmer ratio rises if the HOLMES-powered emitting component has the larger intrinsic ratio.

For simplicity, we set up a restricted spatially resolved analytical model. Its conversion time, molecular depletion time, recycling fraction, feedback loading, and HI stripping coefficient are constant in time at a specified location. The H$\alpha$ budget combines compact HII emission and leaked-OB-powered emission into one effective young-star component, with a constant HOLMES component and adopts the same intrinsic H$\alpha$/H$\beta$ ratio, 2.86, for both. Time $t=0$ denotes the onset of the imposed environmental and supply conditions. It is not an observed time assigned to an infall-stage category.

## 2. Definition of the local two-reservoir model

### 2.1 Key quantities

Following the resolved regulator notation in [Huang et al. (2026), section 4](#ref-huang), we define SFR surface density and molecular GSR surface density as 
$$
\Sigma_{\mathrm{SFR}}(t)\equiv
\frac{\Sigma_{\mathrm{H_2}}(t)}{\tau_{\mathrm{dep}}},
\qquad
\Sigma_\Phi(t)\equiv\frac{\Sigma_{\mathrm{H_2}}(t)}{\tau_{\Phi}}\equiv\frac{\Sigma_{\mathrm{HI}}(t)}{\tau_{\mathrm{conv}}}.
\tag{1}
$$

Here $\tau_{\mathrm{dep}}$ is the molecular depletion time defined relative to the total rate of star formation. The gas supply-rate surface density $\Sigma_\Phi(t)$ feeds the molecular reservoir. The conversion time $\tau_{\mathrm{conv}}$ describes the assumed net effective transfer from the local atomic phase. The second equality is acting as the bridge to connect molecular and actomic phase of Hydrogen gas in our model. Thus, the replenishment timescale of H$_2$ or consumption timescale of HI is 

$$
\tau_\Phi(t)\equiv
\frac{\Sigma_{\mathrm{H_2}}(t)}{\Sigma_\Phi(t)}
=\tau_{\mathrm{conv}}
\frac{\Sigma_{\mathrm{H_2}}(t)}{\Sigma_{\mathrm{HI}}(t)}.
\tag{2}
$$

Consequently, $\Sigma_\Phi(t)$ is a rate, whereas $\tau_\Phi$ and $\tau_{\mathrm{conv}}$ are times. Even when $\tau_{\mathrm{conv}}$ is constant, $\tau_\Phi$ generally evolves. 

Let $0\leq R<1$ be the prompt stellar mass return fraction and $\eta\geq0$ the feedback mass-loading factor, so that the feedback mass loss rate is $\eta\Sigma_{\mathrm{SFR}}(t)$. The net removal associated with star formation and feedback is $(1-R+\eta)\Sigma_{\mathrm{SFR}}(t)$. This is the usual regulator bookkeeping, here assigned effectively to the molecular reservoir (e.g., [Lilly et al. 2013](#ref-lilly); [Huang et al. 2026](#ref-huang)). 

We use only $\gamma_{\mathrm{strip}}$ for the direct RPS loss coefficient and it has units of inverse time and acts only on HI gas. The baseline has no direct molecular stripping, and no external supply to HI after $t=0$ due to the fact of starvation in cluster environment. Again, $\tau_{\mathrm{conv}}$, $\tau_{\mathrm{dep}}$, $R$, $\eta$, and $\gamma_{\mathrm{strip}}$ are constant in time at fixed location $\boldsymbol{x}$. Their possible spatial dependence is retained conceptually. Constant $\tau_{\mathrm{dep}}$ means constant molecular SFE, not constant SFR. A linear molecular law is an empirical first approximation in nearby discs, with substantial environmental and scale-dependent limitations ([Leroy et al. 2013](#ref-leroy)). The HI-only stripping choice is a hypothesis for this calculation, not a claim that molecular stripping never occurs; observations provide counterexamples ([Boselli et al. 2014](#ref-boselli)).

### 2.2 Continuity equations

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

Thus conversion does not destroy gas. The total cold reservoir declines through stripping and net stellar/feedback consumption. Furthermore, we define two physically labelled response rates,

$$
\gamma_{\mathrm{HI}}\equiv\frac{1}{\tau_{\mathrm{conv}}}+\gamma_{\mathrm{strip}},
\qquad
\gamma_{\mathrm{H_2}}\equiv\frac{1-R+\lambda}{\tau_{\mathrm{dep}}}.
\tag{5}
$$

The first is the total fractional removal rate from the HI reservoir, including conversion. The second is the net molecular consumption rate. Both use the same rate notation; $\gamma_{\mathrm{strip}}$ remains the only direct stripping coefficient.

Thus, we can rewrite the equation (3) as
$$
\begin{aligned}
\frac{d\Sigma_{\mathrm{HI}}(t)}{dt}
&=-\gamma_{\mathrm{HI}}\Sigma_{\mathrm{HI}}(t),\\
\frac{d\Sigma_{\mathrm{H_2}}(t)}{dt}
&=\Sigma_\Phi(t)-\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2}}(t).
\end{aligned}
\tag{6}
$$

## 3. Explicit solution and overall SF suppression

### 3.1 Solve the atomic reservoir

For positive initial column $\Sigma_{\mathrm{HI},0}$, the exact solution of atomic gas mass is 
$$
\Sigma_{\mathrm{HI}}(t)=\Sigma_{\mathrm{HI},0}e^{-\gamma_{\mathrm{HI}}t}.
\tag{7}
$$
Immediately, we also have
$$
\Sigma_\Phi(t)=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\tag{8}
$$
### 3.2 Molecular reservoir and explicit SFR solution

First, insert equation (7) and (8) into the second balance:

$$
\frac{d\Sigma_{\mathrm{H_2}}(t)}{dt}
+\gamma_{\mathrm{H_2}}\Sigma_{\mathrm{H_2}}(t)
=\frac{\Sigma_{\mathrm{HI},0}}{\tau_{\mathrm{conv}}}
e^{-\gamma_{\mathrm{HI}}t}.
\tag{9}
$$

Multiplication by $e^{\gamma_{\mathrm{H_2}}t}$ makes the left-hand side a product derivative:

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

The first term is the evolution of surviving initial molecular gas. The second is the evolution of molecular gas supplied after $t=0$, with its subsequent consumption included. 

Finally, dividing every term by the constant depletion time gives the explicit SFR solution:

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

If the HI-reservoir response rate and molecular-consumption rate are equal, $\gamma_{\mathrm{HI}}=\gamma_{\mathrm{H_2}}$, then $1/\tau_{\mathrm{conv}}+\gamma_{\mathrm{strip}}=(1-R+\eta)/\tau_{\mathrm{dep}}$. The integrand in equation (11) is unity, then both H$_2$ and SFR surface densities can be further reduced to:

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

### 3.3 What determines whether the SFR initially rises then decreases or monotonic falls?

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

This explains the delayed decrease. Atomic removal due to RPS immediately reduces future gas supply, but it does not instantly destroy the molecular reservoir or its ongoing star formation.

At late times the slower exponential usually dominates. In the common regime $\gamma_{\mathrm{HI}}>\gamma_{\mathrm{H_2}}$, the molecular response approaches an exponential with timescale $1/\gamma_{\mathrm{H_2}}=\tau_{\mathrm{dep}}/(1-R+\eta)$. It is therefore exponential-like, but generally not one exponential from the onset. For any nonnegative supply, equation (13) also implies

$$
F_{\mathrm{SFR}}(t)\geq e^{-\gamma_{\mathrm{H_2}}t}.
\tag{21}
$$

Increasing HI stripping cannot make retained H$_2$ disappear faster than the no-supply consumption solution (closed-box consumption of H$_2$).

### 3.4 A temporal maximum for an initially over-supplied region

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

Every stationary point has negative curvature. 

Therefore a solution of either $\Sigma_{\mathrm{H_2}}$ or $\Sigma_{\mathrm{SFR}}$ starting with positive slope has one maximum followed by an exponential-like decline; a solution starting with nonpositive slope leads to a monotonic decrease. Finally, we show the indirect route by which atomic deficiency affects molecular gas due to cluster environment, and thus the suppresion of SF (e.g., [Fumagalli et al. 2009](#ref-fumagalli)).

### 3.5 Atomic-gas asymmetry and the leading/trailing SFR contrast 

We now want to explain potential SF enhancement/suppression in the leading/tailing RPS side, e.g., in NGC4654. Here we compare positions $x$ with different initial atomic columns, but the same initial molecular column $\Sigma_{\mathrm{H_2},0}$ and the same $\tau_{\mathrm{conv}}$, $\tau_{\mathrm{dep}}$, $R$, $\eta$, and $\gamma_{\mathrm{strip}}$. Here $x$ denotes a spatial position. The rates in equation (5) are therefore common to all positions. We define $c_{\mathrm{HI}}(x)\equiv\Sigma_{\mathrm{HI}}(x,0)/\Sigma_{\mathrm{HI},0}>0$, where $\Sigma_{\mathrm{HI},0}$ is the reference initial atomic column and $c_{\mathrm{HI}}(x_{\mathrm{reference}})=1$.

This factor represents an atomic excess or deficit established before model onset. Applying the perturbation only to HI is motivated by the greater susceptibility of diffuse atomic gas than dense molecular clouds to ram pressure ([Boselli et al. 2014](#ref-boselli); [Lee et al. 2017](#ref-lee); [Cramer et al. 2020](#ref-cramer)). It is an initial-condition approximation on how RPS act on the compresion of HI gas. Equation (8) and the unchanged supply law immediately give

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

The superscripts denote evaluation at the corresponding position. All three SFRs are equal at $t=0$; their differences develop through the atomic supply. The reference has the same stripping coefficient as the other positions. This is a spatial ordering, not a statement that the leading SFR is increasing with time or that every RPS geometry has the assumed atomic pattern. 

## 4. From instantaneous SFR to young-star H$\alpha$ emission

Let $\mathcal L_\alpha^{\mathrm{young}}$ be the H$\alpha$ luminosity per adopted area powered by the young stellar population. It includes both compact HII emission and emission powered by O/B star photons absorbed outside compact HII regions, i.e. the leaky photons. For the main calculation these are one effective component. 

The conversion from SFR to an ionizing population normally averages the recent star formation history over the lifetimes of massive stars ([Kennicutt & Evans 2012](#ref-ke)). For gas evolution much slower than that stellar response, we use

$$
\boxed{\mathcal L_\alpha^{\mathrm{young}}(t)
\simeq\frac{1}{C_\alpha}\Sigma_{\mathrm{SFR}}(t)
=\mathcal L_{\alpha,0}^{\mathrm{young}}F_{\mathrm{SFR}}(t).}
\tag{31}
$$

$C_\alpha$ is the H$\alpha$-to-SFR calibration, which already absorbs stellar population and the assumed IMF. In MAUVE observation the adopted value is $C_\alpha=4.9835821\times10^{-42}\ M_\odot\,\mathrm{yr^{-1}}/(\mathrm{erg\,s^{-1}})$ with Chabrier IMF. 

## 5. HOLMES's contribution to H$\alpha$ budget

Let $\Sigma_*^{\mathrm{old}}$ be the current mass surface density in the old population, including its associated remnants, and let $q_{\mathrm{H,HOLMES}}$ be the production rate of hydrogen-ionizing photons per unit of that current mass. Their product is the emitted photon production per area:

$$
\begin{aligned}
\mathcal Q_{\mathrm{HOLMES}}&=q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}}.
\end{aligned}
\tag{32}
$$

Here we adopt $q_{\mathrm{H,HOLMES}}=7\times10^{40}\ \mathrm{s^{-1}}M_\odot^{-1}$ from the PEGASE normalization used by [Belfiore et al. (2022), section 3.2](#ref-belfiore).  In ionization equilibrium, one absorbed ionizing photon balances a Case-B recombination. The probability that such a recombination produces H$\alpha$ is $p_\alpha=\alpha_\alpha^{\mathrm{eff}}/\alpha_B$. Each emitted H$\alpha$ photon carries energy $h_{\mathrm P}\nu_\alpha$. Therefore

$$
\begin{aligned}
\mathcal L_\alpha^{\mathrm{HOLMES}}
&=h\nu_\alpha p_\alpha\mathcal Q_{\mathrm{HOLMES}}\\
&=\boxed{\epsilon_\alpha 
q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}}},\\
\epsilon_\alpha&\equiv h\nu_\alpha p_\alpha,\qquad
p_\alpha\equiv\frac{\alpha_\alpha^{\mathrm{eff}}}{\alpha_B}\simeq\frac{1}{2.206}.
\end{aligned}
\tag{33}
$$

$\alpha_B$ excludes recombinations directly to the hydrogen ground state; $\alpha_\alpha^{\mathrm{eff}}$ counts recombinations yielding H$\alpha$ [(Hummer & Storey 1987)](#ref-hs). With $\lambda_\alpha=6562.8\AA$, $\epsilon_\alpha=1.3721\times10^{-12}$ erg per absorbed photon. The adopted numerical conversion $p_\alpha$ is given in [Cid Fernandes et al. (2011), equation 2](#ref-cid), for photoionization by populations older than $10^8$ yr. Therefore, since the old stellar mass ($\Sigma_*^{\mathrm{old}}$) assigned to the region changes little, we consider HOLMES's H$\alpha$ emission as constant (or at least relatively constant compared to the young stellar's contribution).

## 6. How the luminosity weights change

Add the two positive Halpha contributions before defining their weights:

$$
\begin{aligned}
\mathcal L_\alpha(t)&=\mathcal L_\alpha^{\mathrm{young}}(t)+\mathcal L_\alpha^{\mathrm{HOLMES}},\\
w_{\mathrm{HOLMES}}(t)&\equiv
\frac{\mathcal L_\alpha^{\mathrm{HOLMES}}}
{\mathcal L_\alpha^{\mathrm{young}}(t)+\mathcal L_\alpha^{\mathrm{HOLMES}}},\\
w_{\mathrm{young}}(t)&=1-w_{\mathrm{HOLMES}}(t).
\end{aligned}
\tag{34}
$$

These are light fractions, not fractions of area or gas mass. Write $w_{\mathrm{HOLMES},0}=w_{\mathrm{HOLMES}}(0)$. For fixed positive HOLMES luminosity, the quotient rule gives

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
\tag{35}
$$

This connects the gas balance directly to the changing source weight. Molecular consumption exceeding replenishment makes SFR decline and the HOLMES fraction rise. An initially over-supplied region has the opposite response until its SFR maximum. Thus, over the cluster corssing timescale, we expect to see the increase of HOLMES's weight.

## 7. The BPT line ratios and the complete connection

 For a forbidden line $\ell$, let $B$ denote its Balmer denominator and define the linear component ratio $R_{\ell/B}^j=\mathcal L_\ell^j/\mathcal L_B^j$, with $j$ belongs to young or HOLMES. N2 is [N II]6583/Halpha, S2 is ([S II]6716+[S II]6731)/Halpha, and O3 is [O III]5007/Hbeta. Due to the harder ionization of HOLMES, we have 
$$
R_{\ell/B}^{\mathrm{young}}<R_{\ell/B}^{\mathrm{HOLMES}}.
$$
For simplicity, these two terms are set to be constant to indicate the fact that these two sources have representative BPT line ratios; in reality, temperature, ionic fractions, metallicity, N/O, and ionization parameter can change the line ratios ([Byler et al. 2019](#ref-byler)). 

We adopt the intrinsic Balmer ratio $\mathcal L_\alpha^j/\mathcal L_\beta^j=2.86$ for each component, an approximation consistent with the Case-B convention used for the observational dust correction ([Hummer & Storey 1987](#ref-hs)). Thus

$$
\frac{\mathcal L_\beta^{\mathrm{HOLMES}}}{\mathcal L_\beta}
=\frac{\mathcal L_\alpha^{\mathrm{HOLMES}}/2.86}{\mathcal L_\alpha/2.86}
=w_{\mathrm{HOLMES}}.
\tag{37}
$$

Substitute $\mathcal L_\ell^j=R_{\ell/B}^j\mathcal L_B^j$ into the total ratio yields the evolution of BPT line ratio, which show the similar forms of [equations 7 and 8 in Blanc et al. (2009)](#ref-blanc):

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
\tag{38}
$$

This luminosity-weighted identity has an HII/DIG antecedent in . For constant component spectra, differentiate equation (38):

$$
\frac{dR_{\ell/B}(t)}{dt}
=(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\frac{dw_{\mathrm{HOLMES}}(t)}{dt}.
\tag{39}
$$

When young emission fades, the sign of the ratio change is the sign of the spectral contrast. Harder ionization alone does not require every ratio to increase. Both Balmer and forbidden-line luminosities can decline while their ratio rises because the Balmer line fades faster.

Finally, inserting the gas solution, the young-star conversion, and the old-star normalization gives the complete prediction:

$$
\boxed{\begin{aligned}
\mathcal L_\alpha(t)
&=\frac{\Sigma_{\mathrm{SFR},0}}{C_\alpha}F_{\mathrm{SFR}}(t)
+\epsilon_\alpha q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}},\\
R_{\ell/B}(t)
&=R_{\ell/B}^{\mathrm{young}}
+(R_{\ell/B}^{\mathrm{HOLMES}}-R_{\ell/B}^{\mathrm{young}})
\frac{\epsilon_\alpha q_{\mathrm{H,HOLMES}}\Sigma_*^{\mathrm{old}}}
{\mathcal L_\alpha(t)}.
\end{aligned}}
\tag{40}
$$

## References

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


<span id="ref-cramer"></span>
**Cramer, W. J., et al. (2020).** *ALMA evidence for ram pressure compression and stripping of molecular gas in the Virgo cluster galaxy NGC 4402.* ApJ, 901. [DOI](https://doi.org/10.3847/1538-4357/abaf54); [primary manuscript](https://arxiv.org/pdf/1910.14082).

<span id="ref-fumagalli"></span>
**Fumagalli, M., Krumholz, M. R., Prochaska, J. X., Gavazzi, G., & Boselli, A. (2009).** *Molecular hydrogen deficiency in HI-poor galaxies and its implications for star formation.* ApJ, 697, 1811--1821. [DOI](https://doi.org/10.1088/0004-637X/697/2/1811); [primary manuscript](https://arxiv.org/abs/0903.3950).


<span id="ref-huang"></span>
**Huang, R., et al. (2026).** *MAUVE-MUSE: When Metallicity Follows or Fights Star Formation--A Mass-Dependent Inversion in Virgo Galaxies.* MNRAS, 549, stag1019. [DOI](https://doi.org/10.1093/mnras/stag1019); [primary manuscript](https://arxiv.org/html/2605.31412v1).


<span id="ref-hs"></span>
**Hummer, D. G., & Storey, P. J. (1987).** *Recombination-line intensities for hydrogenic ions--I. Case B calculations for HI and HeII.* MNRAS, 224, 801--820. [DOI](https://doi.org/10.1093/mnras/224.3.801).


<span id="ref-ke"></span>
**Kennicutt, R. C., Jr., & Evans, N. J. II (2012).** *Star Formation in the Milky Way and Nearby Galaxies.* ARA&A, 50, 531--608. [DOI](https://doi.org/10.1146/annurev-astro-081811-125610); [primary manuscript](https://arxiv.org/abs/1204.3552).

<span id="ref-lee"></span>
**Lee, B., et al. (2017).** *The effect of ram pressure on the molecular gas of galaxies: three case studies in the Virgo cluster.* MNRAS, 466, 1382--1398. [DOI](https://doi.org/10.1093/mnras/stw3162); [primary manuscript](https://arxiv.org/pdf/1701.02750).

<span id="ref-leroy"></span>
**Leroy, A. K., et al. (2013).** *Molecular Gas and Star Formation in Nearby Disk Galaxies.* AJ, 146, 19. [DOI](https://doi.org/10.1088/0004-6256/146/2/19); [primary manuscript](https://arxiv.org/pdf/1301.2328).

<span id="ref-lilly"></span>
**Lilly, S. J., Carollo, C. M., Pipino, A., Renzini, A., & Peng, Y. (2013).** *Gas Regulation of Galaxies: The Evolution of the Cosmic Specific Star Formation Rate, the Metallicity-Mass-Star-formation Rate Relation, and the Stellar Content of Halos.* ApJ, 772, 119. [DOI](https://doi.org/10.1088/0004-637X/772/2/119); [primary manuscript](https://arxiv.org/abs/1303.5059).