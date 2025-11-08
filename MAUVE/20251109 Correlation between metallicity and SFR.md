# 20251106 Correlation between metallicity and SFR

## Setup of our spatially-resolved bathtub model

All the properties here are written in surface density to highlight the spatially-resolved secnario. We want to highlight that the following model described here is only applicable to the evolution of gas in local universe (present day). Therefore, we don't consider the past history of either star formation or metal enrichment and simply assume each spaxel share similar initial properties when they start to evolve under the following model. Also, since we only consider the evolution in relatively shorter period ($\lesssim$ Gyr) compared to the past history ($\sim10$ Gyr), so the change of $M_*$ under this framework is negligible (or at most less than observation uncertainty) and thus we can ignore the deivation of $M_*$ compared to the initial value.  

First, as usual, we define the SFR surface density ($\Sigma_\mathrm{SFR}[t_m]$) in the following way: 
$$
\begin{split}
\Sigma_\mathrm{SFR}[t] \equiv \frac{\Sigma_\mathrm{g}[t]}{\tau_\mathrm{dep}[t]}
\end{split}
$$
where $\Sigma_\mathrm{g}[t_m]$ is the gas mass surface density and $\tau_\mathrm{dep}$ is the depletion timescale. Also, we start with the change of gas mass surface density:
$$
\begin{split}
\frac{d}{dt}\Sigma_\mathrm{g}[t] = \Sigma_\Phi[t]-(1-R+\lambda)\Sigma_\mathrm{SFR}[t]
\end{split}
$$
with $R$ is the return fraction from stars, $\lambda$ is mass loading factor, $\Sigma_\Phi[t_m]$ is gas supply rate surface density, which denotes the ability to gain supply gas; 

and the change of metal mass surface density:
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
Then, we have the change of metallicity: 
$$
\begin{split}
\frac{dZ[t]}{dt} = \frac{1}{\Sigma_\mathrm{g}[t]}(y(1-R)\Sigma_\mathrm{SFR}[t]+(Z_\Phi-Z[t])\Sigma_\Phi[t]) \\
\end{split}
$$
## Correct interpretation of setting $\frac{dZ[t]}{dt}=0$: local minimum/maximum metallicities

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
Clearly, only at $t=t_m$ that the realtime metallicity is equal to the instantaneous target (i.e., $Z[t_m] = Z^\dagger[t_m]$), the real-time metallicity will reach local minimum or maximum.  Note that in [Lilly+2013](https://iopscience.iop.org/article/10.1088/0004-637X/772/2/119/meta), they also take $\frac{dZ[t_m]}{dt}=0$ to obtain the so-called equilibrium metallicity under such condition, and say that the system will rapidly drive the metallicity to the equilibrium value within the timescale 
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
\tau_\Phi[t] \equiv \frac{\Sigma_g(t)}{\Sigma_\Phi(t)} \leq\tau_\mathrm{dep} \\
\end{split}
$$
That means $Z^\dagger[t]$ is the target value if the realtime metallicity $Z[t]$ keep the current change rate within the gas supply timescale. 

Now, suepose we observe a set of spaxels, and each spaxel has correspoding properties. Then, given any pair of spaxels, we can study the deviation of each property. Also, mathmatically speaking, for each property, although we observe different values from different spaxels at the same time, it is equivalent to consider them as time variation under the same function expression. Alternatively, we treat $t$ as a parameter to control different values of invididual property and trace the corresponding values in other properties. For example,  let's say we observe spaxel 1 and 2 at the same time, and we denote their properties as $p_1$ and $p_2$, where $p$ can be any properties: $Z$, $\Sigma_\mathrm{SFR}$, etc. If we plugin their values into the continuity ODEs, we can find that they return the corresponding $t_1$ and $t_2$ (i.e., $p_1 = p[t_1]$ and $p_2 = p[t_2]$). However, different values of $t_1$ and $t_2$ do not contradict the fact that we observe $p_1$ and $p_2$ simultaneously, because as mentioned previously, we don't care (know) the past history of each spaxel, so the $t$ only indicates how long have these spaxels evolved under this framework. 

Here, to understand the correlation between metallicity and SFR, depending on if the metallicity satisfies $t=t_m$ or not, we can separate them into the following two cases:

1. Local minimum/maximum metallicities. In other words, such data have $Z[t_m] = Z^\dagger[t_m]$. Hence, we study the special and limited case: the correlation between $Z^\dagger[t_m]$ and their correponding $\Sigma_\mathrm{SFR}[t_m]$. 
2. Non-local minimum/maximum metallicities. For the rest of the data, which should represent the majority of observational data, they can be treated as deivations from the "nearby" (or the closest) local minimum/maximum moment ($t_m$). To avoid confusion, we introduce a new parameter $\zeta$ to re-express the properties of these data: $p[t\neq t_m] = p[t_m+\zeta]$. Hence, in this more common case, the question we interested will become for any given $Z[t_m+\zeta]$ and $\Sigma_\mathrm{SFR}[t_m+\zeta]$, how do they compare to the closest local minimum/maximum values $Z[t_m]$ and $\Sigma_\mathrm{SFR}[t_m]$. 

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
Since we always have $Z^\dagger[t_m]\geq Z_\Phi$, then the sign of Oxygen abundance deviation in this case ($\delta(12+\log_{10}(\mathrm{O/H})[t_m])$) is only depends on the bracket on the RHS ($(\delta\log_{10}(\Sigma_\mathrm{SFR}[t_m])-\delta\log_{10}(\Sigma_\Phi[t_m]))$), which means the relative change of $\log_{10}(\Sigma_\mathrm{SFR}[t_m])$ and $\log_{10}(\Sigma_\mathrm{in}[t_m])$, respectively. We will get back to this relation later. 

Note that, if we generalize $t_m$ to $t$, which in fact studies the change of instantaneous target metallicity rather than the realtime metallicity, we still have the same equation:
$$
\begin{split}
\delta\log_{10}(Z^\dagger[t]) = \frac{\delta Z^\dagger[t]}{\ln(10)Z^\dagger[t]} = \frac{Z^\dagger[t]-Z_\Phi}{Z^\dagger[t]}(\delta\log_{10}(\Sigma_\mathrm{SFR}[t])-\delta\log_{10}(\Sigma_\Phi[t])) \\
\end{split}
$$

## Correlation for the non-local minimum/maximum metallicities

With the new parameter $t = t_m+\zeta$, where $t_m$ is constant in this case (and thus any $p[t_m]$ is also constant), for any given property $p[t]$, we can use $\delta p[\zeta] = p[t_m+\zeta]-p[t_m]$ to indicate the deviation from the local minimum/maximum metallicity moment $t_m$. Note that at the moment $\zeta = 0 \Leftrightarrow t = t_m$, we have $\delta p[\zeta=0] = 0$.  Hence, we can rewrite the $\frac{dZ[t]}{dt}$ in $\zeta$ space: 
$$
\begin{split}
\frac{dZ[t_m+\zeta]}{d(t_m+\zeta)}&=\frac{1}{\tau_\Phi[t_m+\zeta]}(Z^\dagger[t_m+\zeta]-Z[t_m+\zeta]) \\
\Leftrightarrow \frac{d(Z[t_m]+\delta Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]+\delta\tau_\Phi[\zeta]}((Z^\dagger[t_m]+\delta Z^\dagger[\zeta])-(Z[t_m]+\delta Z[\zeta])) \\
\Leftrightarrow \frac{dZ[t_m]}{dt}+\frac{d(\delta Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]+\delta\tau_\Phi[\zeta]}(\delta Z^\dagger[\zeta]-\delta Z[\zeta]) \\
\Leftrightarrow \frac{d(\delta Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]+\delta\tau_\Phi[\zeta]}(\delta Z^\dagger[\zeta]-\delta Z[\zeta]) \\
\Leftrightarrow \frac{d(\delta Z[\zeta])}{d\zeta}&=\frac{1}{\tau_\Phi[t_m]+\delta\tau_\Phi[\zeta]}(\delta Z^\dagger[\zeta]-\delta Z[\zeta]) \\
\end{split}
$$
Using first-oder Taylor expansion, for $\delta\tau_\Phi[\zeta]<<\tau_\Phi[t_m]$, we have:
$$
\begin{split}
\frac{1}{\tau_\Phi[t_m] + \delta\tau_\Phi[\zeta]}=\frac{1}{\tau_\Phi[t_m](1+\frac{\delta\tau_\Phi[\zeta]}{\tau_\Phi[\zeta]})} \approx \frac{1}{\tau_\Phi[t_m]}(1-\frac{\delta\tau_\Phi[\zeta]}{\tau_\Phi[t_m]}) \\
\end{split}
$$
Thus, we can express the ODE and neglect the second order deviation to have the final ODE that denote the relation between $Z[\zeta]$ and $Z^\dagger[\zeta]$:
$$
\begin{split}
\frac{d(\delta Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]}(1-\frac{\delta\tau_\Phi[\zeta]}{\tau_\Phi[t_m]})(\delta Z^\dagger[\zeta]-\delta Z[\zeta]) \\
\Leftrightarrow \frac{d(\delta Z[\zeta])}{dt}&=\frac{1}{\tau_\Phi[t_m]}(\delta Z^\dagger[\zeta]-\delta Z[\zeta]) \\
\Leftrightarrow \frac{d(\delta Z[\zeta])}{dt}+\frac{\delta Z[\zeta]}{\tau_\Phi[t_m]}&=\frac{\delta Z^\dagger[\zeta]}{\tau_\Phi[t_m]} \\
\end{split}
$$
Now we can take the one-side Laplace transform ($\mathcal{L}\{\dot x\}=sX[s]-x[0]$ ) of the final ODE:
$$
\begin{split}
\mathcal{L}\{\frac{d(\delta Z[\zeta])}{dt}\}+\frac{1}{\tau_\Phi[t_m]}\mathcal{L}\{\delta Z[\zeta]\}&=\frac{1}{\tau_\Phi[t_m]}\mathcal{L}\{\delta Z^\dagger[\zeta]\} \\
\Leftrightarrow s\mathcal{Z}[s]-\delta Z[\zeta=0]+\frac{1}{\tau_\Phi[t_m]}\mathcal{Z}[s]&=\frac{1}{\tau_\Phi[t_m]}\mathcal{Z}^\dagger[s] \\
\Leftrightarrow (s+\frac{1}{\tau_\Phi[t_m]})\mathcal{Z}[s]&=\frac{1}{\tau_\Phi[t_m]}\mathcal{Z}^\dagger[s] \\
\Leftrightarrow \mathcal{Z}[s]&=\frac{1}{1+s\cdot\tau_\Phi[t_m]}\mathcal{Z}^\dagger[s] \\
\end{split}
$$
where $\mathcal{Z}[s]$ and $\mathcal{Z}^\dagger[s]$ are $\delta Z[\zeta]$ and $\delta Z^\dagger[\zeta]$ in Laplace $s$ domain. Then, we can take the inverse Laplace transform
$$
\begin{split}
\mathcal{L}^{-1}\{\mathcal{Z}[s]\}&=\mathcal{L}^{-1}\{\frac{1}{1+s\cdot\tau_\Phi[t_m]}\mathcal{Z}^\dagger[s]\} \\
\Leftrightarrow \delta Z[\zeta]&=\frac{1}{\tau_\Phi[t_m]}\int_0^\zeta e^{-(\frac{\zeta-\nu}{\tau_\Phi[t_m]})}\delta Z^\dagger[\nu]d\nu \\
\end{split}
$$
This convolution relation says that $\delta Z[\zeta]$ is low-pass-filtered version of the driver $\delta Z^\dagger[\zeta]$, or in other words, an exponential moving average of input $\delta Z^\dagger[\zeta]$, smoothed over a timescale up to $\tau_\Phi[t_m]$. In astrophyiscs, this equation describes the nature (a process with "memory") that the change of realtime metallicity $\delta Z[\zeta]$ does not instantaneously track the changes in target metallicity $\delta Z^\dagger[\zeta]$, but responds more to recent changes of $\delta Z^\dagger[\zeta]$, weighted by the timescale of $\tau_\Phi[t_m]$. Such delay or smoothing effect may raise from the physical processes like delay in star formation impact (metal injection) and turbulent mixing (metal diffusion).  Note that in practice, due to detection method and measurement uncertainty, our observation is in the regime that $\zeta \geq \tau_\Phi[t_m]$. 

Despite the gradual response of $\delta Z[\zeta]$ and rapid change of $\delta Z^\dagger[\zeta]$, since all weights are positive, $\delta Z[\zeta]$ is always propotional to $\delta Z^\dagger[\zeta]$. Thus, for any data at $t_m+\zeta$, the Oxygen abundance difference compared to the data at $t_m$ will be
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})[\zeta]) &= \delta\log_{10}(Z[\zeta])  \\
&= \frac{1}{\ln10}\frac{\delta Z[\zeta]}{Z[\zeta]} \\
&\propto \delta Z[\zeta] \\
&\propto \delta Z^\dagger[\zeta] \\
&\propto \delta\log_{10}(\Sigma_\mathrm{SFR}[\zeta])-\delta\log_{10}(\Sigma_\Phi[\zeta]) \\
\end{split}
$$
which again indicates that the change of Oxygen abundance difference is proportional to the relative change of $\Sigma_\mathrm{SFR}$ compared to that of $\Sigma_\Phi$. 

