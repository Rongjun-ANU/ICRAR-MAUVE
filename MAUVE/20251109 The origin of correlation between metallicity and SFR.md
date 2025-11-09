# 20251109 The origin of correlation between metallicity and SFR

## Setup of our spatially-resolved bathtub model

All the properties here are written in surface density to highlight the spatially-resolved scenario. We want to highlight that the following model described here is only applicable to the evolution of gas regulator in local universe (present day). Therefore, we don't consider the past history of either star formation or metal enrichment, and we simply assume each spaxel share similar initial properties when they start to evolve under the following model. Also, because we only consider the evolution in relatively shorter period ($\lesssim$ Gyr) compared to the past history ($\sim10$ Gyr), the change of steallar mass surface density ($\Sigma_*$) under this framework is negligible (or at most less than observation uncertainty) and thus we can ignore the deivation of $\Sigma_*$.  

First, we define the SFR surface density ($\Sigma_\mathrm{SFR}[t_m]$) in the following way: 
$$
\begin{split}
\Sigma_\mathrm{SFR}[t] \equiv \frac{\Sigma_\mathrm{g}[t]}{\tau_\mathrm{dep}[t]}
\end{split}
$$
where $\Sigma_\mathrm{g}[t_m]$ is the gas mass surface density and $\tau_\mathrm{dep}$ is the depletion timescale. 

Also, we start with the change of gas mass surface density as the first continuity equation:
$$
\begin{split}
\frac{d}{dt}\Sigma_\mathrm{g}[t] = \Sigma_\Phi[t]-(1-R+\lambda)\Sigma_\mathrm{SFR}[t]
\end{split}
$$
with $R$ is the return fraction from stars, $\lambda$ is mass loading factor, $\Sigma_\Phi[t_m]$ is gas supply rate (GSR) surface density, which denotes the ability to gain supply gas. 

Similarly, we have the change of metal mass surface density as the second continuity equation:
$$
\begin{split}
\frac{d}{dt}\Sigma_Z[t] = y(1-R)\Sigma_\mathrm{SFR}[t]-(1-R+\lambda)\Sigma_\mathrm{SFR}[t]Z[t]+\Sigma_\Phi[t]Z_\Phi \\
\end{split}
$$
with $y$ is the metal yield per stellar and $Z_\Phi$ is the metallicity of supply gas, which is the local background metallicity. Next, we also define the gas-phase metallicity as:
$$
\begin{split}
Z[t]=\frac{\Sigma_Z[t]}{\Sigma_\mathrm{g}[t]}
\end{split}
$$
This yields the change of metallicity: 
$$
\begin{split}
\frac{dZ[t]}{dt} = \frac{1}{\Sigma_\mathrm{g}[t]}(y(1-R)\Sigma_\mathrm{SFR}[t]+(Z_\Phi-Z[t])\Sigma_\Phi[t]) \\
\end{split}
$$

## Correct interpretation of setting $\frac{dZ[t]}{dt}=0$: local minimum/maximum value of metallicity

If we let $\frac{dZ[t]}{dt}=0$ at $t=t_m$, then we can find the local minimum or maximum value of metallicity:
$$
\begin{split}
0 &= \frac{1}{\Sigma_\mathrm{g}[t_m]}(y(1-R)\Sigma_\mathrm{SFR}[t_m]+(Z_\Phi-Z[t_m])\Sigma_\Phi[t_m]) \\
\Leftrightarrow Z[t_m] &= Z_\Phi + y(1-R)\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]} \\
\end{split}
$$
Here we can define the instantaneous target metallicity as 
$$
\begin{split}
Z^\dagger[t]=Z_\Phi + y(1-R)\frac{\Sigma_\mathrm{SFR}[t]}{\Sigma_\Phi[t]} \\
\end{split}
$$
Clearly, only at $t=t_m$ that the realtime metallicity is equal to the instantaneous target (i.e., $Z[t_m] = Z^\dagger[t_m]$), the real-time metallicity will reach local minimum or maximum.  Note that in [Lilly+2013](https://iopscience.iop.org/article/10.1088/0004-637X/772/2/119/meta) and other works, they also take $\frac{dZ[t_m]}{dt}=0$ to obtain the so-called equilibrium metallicity under such condition, and say that the system will rapidly drive the metallicity to the equilibrium value within the timescale 
$$
\begin{split}
\frac{\Sigma_\mathrm{SFR}}{\Sigma_\Phi}\tau_\mathrm{dep}\leq\tau_\mathrm{dep} \\
\end{split}
$$
However, we don't favour this name but instead call it instantaneous target metallicity, because 

1. $\frac{dZ[t_m]}{dt}=0$ doesn't necessiarly mean the relatively steady state; 
2. If we substitute the $Z^\dagger[t_m]$ back to $\frac{dZ[t_m]}{dt}$, we find:

$$
\begin{split}
\frac{dZ[t]}{dt}&=\frac{1}{\tau_\Phi(t)}(Z^\dagger[t]-Z[t]) \\
\Leftrightarrow Z^\dagger[t]&=Z[t]+\frac{dZ[t]}{dt}\tau_\Phi[t] \\
\end{split}
$$

With the gas supply timescale
$$
\begin{split}
\tau_\Phi[t] \equiv \frac{\Sigma_g[t]}{\Sigma_\Phi[t]} \leq\tau_\mathrm{dep}[t] \\
\end{split}
$$
That means $Z^\dagger[t]$ is the target value if the realtime metallicity $Z[t]$ keep the current change rate within the gas supply timescale. 

Now, suepose we observe a set of spaxels, and each spaxel has correspoding properties. Then, given any pair of spaxels, we can study the deviation of each property. Also, mathmatically speaking, for each property, although we observe different values from different spaxels at the same time, it is equivalent to consider them as time variation under the same function expression. Alternatively, we treat $t$ as a parameter to control different values of invididual property and trace the corresponding values in other properties. For example,  let's say we observe spaxel 1 and 2 at the same time, and we denote their properties as $\Rho_1$ and $\Rho_2$, where $\Rho$ can be any properties: $Z$, $\Sigma_\mathrm{SFR}$, etc. If we plugin their values into the continuity ODEs, we can find that they return the corresponding $t_1$ and $t_2$ (i.e., $\Rho_1 = \Rho[t_1]$ and $\Rho_2 = \Rho[t_2]$). However, different values of $t_1$ and $t_2$ do not contradict the fact that we observe $\Rho_1$ and $\Rho_2$ simultaneously, because as mentioned previously, we don't care (know) the past history of each spaxel, so the $t$ only indicates how long have these spaxels evolved under this framework. Hoever, for simplicity, we can solely treat $t$ as an independent variable to trace different properties of each spaxel. 

Here, to understand the correlation between metallicity and SFR, depending on if the metallicity satisfies $t=t_m$ or not, we can separate them into the following two cases:

1. Local minimum/maximum metallicities. In other words, such data have $Z[t_m] = Z^\dagger[t_m]$. Hence, we study the special and limited case: the correlation between $Z^\dagger[t_m]$ and their correponding $\Sigma_\mathrm{SFR}[t_m]$. 

2. Non-local minimum/maximum metallicities. For the rest of the data, which should represent the majority of observational data, they can be treated as deivations from the "nearby" (or the closest) local minimum/maximum moment ($t_m$). To avoid confusion, we introduce a new parameter $\zeta$ to re-express the properties of these data: $\Rho[t\neq t_m] = \Rho[t_m+\zeta]$.  With the new parameter $\zeta = t-t_m$, where $t_m$ is constant in this case (and thus any $\Rho[t_m]$ is also constant), for any given property $\Rho[t]$, we can introduce new notations $\Delta_\Rho[\zeta]$ in $\zeta$ space to indicate individual property's deviation from the local minimum/maximum metallicity moment $t_m$:
   $$
   \begin{split}
   \Delta_\Rho[\zeta] = \Rho[t_m+\zeta]-\Rho[t_m] = \delta \Rho[t] \\
   \end{split}
   $$
   Note that at the moment $\zeta = 0 \Leftrightarrow t = t_m$, we have $\Delta \Rho[\zeta=0] = \delta \Rho[t_m] = 0$.  Hence, in this more common case, the question we want to answer will become for any given $Z[t_m+\zeta]$ and $\Sigma_\mathrm{SFR}[t_m+\zeta]$, how do they compare to the closest local minimum/maximum values $Z[t_m]$ and $\Sigma_\mathrm{SFR}[t_m]$. 

## Correlation for the local minimum/maximum metallicities

In observation, we are normally tracing the Oxygen abundance ($12+\log_{10}(\mathrm{O/H})$) to inidicate the gas-phase metallicity.    Then, for a set of local minimum/maximum observables $12+\log_{10}(\mathrm{O/H})[t_m]$ (i.e., $t_m$ is the variable in this case), if we pick any individual value, then the deviation between it and other observables can be written as 
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})[t_m]) &= \delta\log_{10}(Z^\dagger[t_m])  \\
&= \frac{1}{\ln10}\frac{\delta Z^\dagger[t_m]}{Z^\dagger[t_m]} \\
&= \frac{1}{\ln10}\frac{y(1-R)\delta(\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]})}{Z_\Phi+y(1-R)\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]}} \\
&= \frac{1}{\ln10}\frac{y(1-R)\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]}}{Z_\Phi+y(1-R)\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]}}\frac{\delta(\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]})}{\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]}} \\
&= \frac{y(1-R)\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]}}{Z_\Phi+y(1-R)\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]}}(\delta\log_{10}(\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]})) \\
&= \frac{y(1-R)\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]}}{Z_\Phi+y(1-R)\frac{\Sigma_\mathrm{SFR}[t_m]}{\Sigma_\Phi[t_m]}}(\delta\log_{10}(\Sigma_\mathrm{SFR}[t_m])-\delta\log_{10}(\Sigma_\Phi[t_m])) \\
&= \frac{Z^\dagger[t_m]-Z_\Phi}{Z^\dagger[t_m]}(\delta\log_{10}(\Sigma_\mathrm{SFR}[t_m])-\delta\log_{10}(\Sigma_\Phi[t_m])) \\
\end{split}
$$
Since we always have $Z^\dagger[t_m]\geq Z_\Phi$, then the sign of Oxygen abundance deviation in this case ($\delta(12+\log_{10}(\mathrm{O/H})[t_m])$) is only depends on the bracket on the RHS (i.e., $\delta\log_{10}(\Sigma_\mathrm{SFR}[t_m])-\delta\log_{10}(\Sigma_\Phi[t_m])$).

Note that, if we generalize $t_m$ to $t$, which in fact describes the change of instantaneous target metallicity rather than the realtime metallicity, we still have the same equation:
$$
\begin{split}
\delta\log_{10}(Z^\dagger[t]) = \frac{\delta Z^\dagger[t]}{\ln(10)Z^\dagger[t]} = \frac{Z^\dagger[t]-Z_\Phi}{Z^\dagger[t]}(\delta\log_{10}(\Sigma_\mathrm{SFR}[t])-\delta\log_{10}(\Sigma_\Phi[t])) \\
\end{split}
$$

This tells us a broader context that change of $\log$ scale instantaneous target metallicity is dertermined by the subtraction between change of $\log$ scale SFR surface density and that of GSR surface density. We will get back to this relation later. 

## Correlation for the non-local minimum/maximum metallicities

With our previous definition, we can rewrite the $\frac{dZ[t]}{dt}$ in $\zeta$ space: 
$$
\begin{split}
\frac{dZ[t_m+\zeta]}{d(t_m+\zeta)}&=\frac{1}{\tau_\Phi[t_m+\zeta]}(Z^\dagger[t_m+\zeta]-Z[t_m+\zeta]) \\
\Leftrightarrow \frac{d(Z[t_m]+\Delta_Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]+\Delta_{\tau_\Phi}[\zeta]}((Z^\dagger[t_m]+\Delta_{Z^\dagger}[\zeta])-(Z[t_m]+\Delta_Z[\zeta])) \\
\Leftrightarrow \frac{dZ[t_m]}{dt}+\frac{d(\Delta_Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]+\Delta_{\tau_\Phi}[\zeta]}(\Delta_{Z^\dagger}[\zeta]-\Delta_Z[\zeta]) \\
\Leftrightarrow \frac{d(\Delta_Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]+\Delta_{\tau_\Phi}[\zeta]}(\Delta_{Z^\dagger}[\zeta]-\Delta_Z[\zeta]) \\
\Leftrightarrow \frac{d(\Delta_Z[\zeta])}{d\zeta}&=\frac{1}{\tau_\Phi[t_m]+\Delta_{\tau_\Phi}[\zeta]}(\Delta_{Z^\dagger}[\zeta]-\Delta_Z[\zeta]) \\
\end{split}
$$
Using first-oder Taylor expansion, for $\Delta_{\tau_\Phi}[\zeta]<<\tau_\Phi[t_m]$, we have:
$$
\begin{split}
\frac{1}{\tau_\Phi[t_m] + \Delta_{\tau_\Phi}[\zeta]}=\frac{1}{\tau_\Phi[t_m](1+\frac{\Delta_{\tau_\Phi}[\zeta]}{\tau_\Phi[\zeta]})} \approx \frac{1}{\tau_\Phi[t_m]}(1-\frac{\Delta_{\tau_\Phi}[\zeta]}{\tau_\Phi[t_m]}) \\
\end{split}
$$
Thus, we can express the ODE and neglect the second order deviation to have the final ODE that describes the relation between $\Delta_Z [\zeta]$ and $\Delta_{Z^\dagger} [\zeta]$:
$$
\begin{split}
\frac{d(\Delta_Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]}(1-\frac{\Delta_{\tau_\Phi}[\zeta]}{\tau_\Phi[t_m]})(\Delta_{Z^\dagger}[\zeta]-\Delta_Z[\zeta]) \\
\Leftrightarrow \frac{d(\Delta_Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]}(\Delta_{Z^\dagger}[\zeta]-\Delta_Z[\zeta]) \\
\Leftrightarrow \frac{d(\Delta_Z[\zeta])}{dt}+\frac{\Delta_Z[\zeta]}{\tau_\Phi[t_m]}&=\frac{\Delta_{Z^\dagger}[\zeta]}{\tau_\Phi[t_m]} \\
\end{split}
$$
Now we can take the one-side Laplace transform ($\mathcal{L}\{\dot x\}=sX[s]-x[0]$ ) of the final ODE:
$$
\begin{split}
\mathcal{L}\{\frac{d(\Delta_Z[\zeta])}{dt}\}+\frac{1}{\tau_\Phi[t_m]}\mathcal{L}\{\Delta_Z[\zeta]\}&=\frac{1}{\tau_\Phi[t_m]}\mathcal{L}\{\Delta_{Z^\dagger}[\zeta]\} \\
\Leftrightarrow s\mathcal{Z}[s]-\Delta_Z[\zeta=0]+\frac{1}{\tau_\Phi[t_m]}\mathcal{Z}[s]&=\frac{1}{\tau_\Phi[t_m]}\mathcal{Z}^\dagger[s] \\
\Leftrightarrow (s+\frac{1}{\tau_\Phi[t_m]})\mathcal{Z}[s]&=\frac{1}{\tau_\Phi[t_m]}\mathcal{Z}^\dagger[s] \\
\Leftrightarrow \mathcal{Z}[s]&=\frac{1}{1+s\cdot\tau_\Phi[t_m]}\mathcal{Z}^\dagger[s] \\
\end{split}
$$
where $\mathcal{Z}[s]$ and $\mathcal{Z}^\dagger[s]$ are $\Delta_Z[\zeta]$ and $\Delta_{Z^\dagger}[\zeta]$ in Laplace $s$ domain. 

Then, we can take the inverse Laplace transform
$$
\begin{split}
\mathcal{L}^{-1}\{\mathcal{Z}[s]\}&=\mathcal{L}^{-1}\{\frac{1}{1+s\cdot\tau_\Phi[t_m]}\mathcal{Z}^\dagger[s]\} \\
\Leftrightarrow \Delta_Z[\zeta]&=\frac{1}{\tau_\Phi[t_m]}\int_0^\zeta e^{-(\frac{\zeta-\nu}{\tau_\Phi[t_m]})}\Delta_{Z^\dagger}[\nu]d\nu \\
\end{split}
$$
This convolution relation says that $\Delta_Z[\zeta]$ is low-pass-filtered version of the driver $\Delta_{Z^\dagger}[\zeta]$, or in other words, an exponential moving average of input $\Delta_{Z^\dagger}[\zeta]$, smoothed over a timescale $\tau_\Phi[t_m]$. In astrophyiscs, this equation describes the nature (a process with "memory") that the change of realtime metallicity $\Delta_Z[\zeta]$ does not instantaneously track the changes in target metallicity $\Delta_{Z^\dagger}[\zeta]$, but responds more to recent changes of $\Delta_{Z^\dagger}[\zeta]$, weighted by the timescale of $\tau_\Phi[t_m]$. Such delay or smoothing effect may raise from the physical processes like delay in star formation impact (metal injection) and turbulent mixing (metal diffusion). Note that in practice, due to detection method and measurement uncertainty, our observation is in the regime that $\zeta \geq \tau_\Phi[t_m]$. 

Recall that $\Delta_Z[\zeta] = \delta Z[t]$ and $\Delta_{Z^\dagger}[\zeta] = \delta Z^\dagger[t]$, we can convert the equation above back to $t$ space
$$
\begin{split}
\delta Z[t]&=\frac{1}{\tau_\Phi[t_m]}\int_0^{t-t_m} e^{-(\frac{t-t_m-\nu}{\tau_\Phi[t_m]})}\delta Z^\dagger[\nu]d\nu \\
\end{split}
$$
Despite the gradual response of $\delta Z[t]$ (or$\Delta_Z[\zeta]$) and rapid change of $\delta Z^\dagger[t]$ (or $\Delta_{Z^\dagger}[\zeta]$), since all weights are positive, then $\delta Z[t]$ (or$\Delta_Z[\zeta]$) is always propotional to $\delta Z^\dagger[t]$ (or $\Delta_{Z^\dagger}[\zeta]$). This tells us, in observation, the change of realtime metallicity is propotional to the change of instantaneous target metallicity. Thus, for any observed data at $t= t_m+\zeta$ with $\zeta \geq \tau_\Phi[t_m]$, the Oxygen abundance difference compared to the data at $t_m$ will be
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})[t]) &= \delta\log_{10}(Z[t])  \\
&= \frac{1}{\ln10}\frac{\delta Z[t]}{Z[t]} \\
&\propto \delta Z^\dagger[t] \\
&\propto \delta\log_{10}(\Sigma_\mathrm{SFR}[t])-\delta\log_{10}(\Sigma_\Phi[t]) \\
\end{split}
$$
which again indicates that the change of Oxygen abundance is set by the relative change of $\log_{10}(\Sigma_\mathrm{SFR})$ compared to that of $\log_{10}(\Sigma_\Phi)$. 

## Understand the correlation between metallicity and SFR

Now we have proved that for observational data, there exits such relation:
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})[t]) &\propto \delta\log_{10}(\Sigma_\mathrm{SFR}[t])-\delta\log_{10}(\Sigma_\Phi[t]) \\
\end{split}
$$
Obviously, the sign of $\delta(12+\log_{10}(\mathrm{O/H})[t])$, which reflects the increase or decrease of Oxygen abundance, is controlled by if $\delta\log_{10}(\Sigma_\mathrm{SFR}[t])>\delta\log_{10}(\Sigma_\Phi[t])$, which means if the relative change of SFR surface density is greater than that of GSR surface density.  Hence, this leaves us three cases:

1. Positive correlation between Oxygen abundance and SFR surface density. In the case that the change of GSR surface density is higher than the change of GSR surface density, the increase or decrease of $\Sigma_\mathrm{SFR}$ will be directly reflected by the increase or decrease of $12+\log_{10}(\mathrm{O/H})$, resulting with a positive correlation. 
2. Negative correlation between Oxygen abundance and SFR surface density. In the case that the change of GSR surface density is negligible compared to relative change of GSR surface density, even though there may still be increase or decrease of $\Sigma_\mathrm{SFR}$, the increase or decrease of $\Sigma_\Phi[t]$  will control the decrease or increase of $12+\log_{10}(\mathrm{O/H})$ by the minus sign. Note that to further end up with negative correlation, it requires SFR surface density and GSR surface density increases or decreases simultaneously, otherwise it still leads to positive correlation. 
3. No correlation between Oxygen abundance and SFR surface density. In the case that the change of SFR surface density and change of GSR surface density are roughly balanced, there is nearly no variation on $12+\log_{10}(\mathrm{O/H})$, leading to no correlation. 

More specific, recalling the definitions of $\tau_\mathrm{dep}\equiv\frac{\Sigma_\mathrm{g}[t]}{\Sigma_\mathrm{SFR}[t]}$ and $\tau_\Phi \equiv \frac{\Sigma_g[t]}{\Sigma_\Phi[t]}$, we can further reduce the RHS of the relation to be the comparison of relative change of these two timescales:
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})[t]) &\propto \delta\log_{10}(\tau_\Phi[t])-\delta\log_{10}(\tau_\mathrm{dep}[t]) \\
\end{split}
$$
In our data, with binning by stellar mass surface density ($\Sigma_*$), we found that it tends to show positive correlation in high mass end, while negative correlation in low mass end, with invertion point at $\Sigma_*\sim7.6M_\odot/\mathrm{kpc}^2$. This phenomenon suggests that in the denser region of a galaxy (typically the inner part), the variation of SFR (or depletion timescale) across different spaxels in similar mass range is more severe than that of GSR (or its timescale). This is expected, because there is more stochasticity and burstiness of star formation in this region. On the other hand, in the outskirts of a galaxy (less $\Sigma_*$), the observed negative correlation raises from the fact that the variation of GSR (or gas supply timescale) becomes more prominent, due to the more intense effect of either gas accreting or stripping. Particularly, for our data, in cluster environment, RPS will significantly disturb the ability to get gas supply for low-mass SF regions in outskirts. In general, the positive correlation will emerge from a more close-box-like environment (the ability to obtain gas supply remasin relatively unchanged), while as the system is getting more interactions with or more easily to be affected by the nearby environment to change that ability, the correlation will become less postitive and towards no correlation or even become negative.  

Note that in [WL21](https://doi.org/10.3847/1538-4357/abe413), they conclude that the positive or negative correlation is due to either dominated by time-varying star formation efficiency ($\epsilon[t]\equiv\frac{\Sigma_\mathrm{SFR}[t]}{\Sigma_\mathrm{gas}[t]}$) or time-varying inflow rate, which is $\Sigma_\Phi[t]$ in our resovled scenario. In this case, the relation will become:
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})[t]) &\propto \delta\log_{10}(\epsilon[t])+\delta\log_{10}(\Sigma_\mathrm{gas})-\delta\log_{10}(\tau_\mathrm{dep}[t]) \\
\end{split}
$$
Given the fact that the gas mass surface density is another long-term accumulation so its relative change will be also negligible, the relation above reflects their findings. 
