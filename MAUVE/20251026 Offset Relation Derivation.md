# 20251026 Offset Relation Derivation

Here we start with two continuity equations from Lilly+2013 or Wang+2021:
$$
\begin{split}
\frac{dM_\mathrm{g}[t]}{dt} &= \Phi[t] - (1-R)\rm{SFR}[t]-\lambda\mathrm{SFR}[t]  \\
&= \Phi[t] - (1-R+\lambda)\mathrm{SFR}[t]  \\
\label{eq:dMgdt}
\end{split}
$$
Here $M_\rm{g}[t]$ is the gas mass, $\Phi[t]$ is the inflow rate, $R$ is the returned fraction from stars and $\lambda$ is the mass loading factor. So these 3 terms are gas mass changes due to inflow, loss to star formation and loss to outflow. 
$$
\begin{split}
\frac{dM_Z[t]}{dt} = y(1-R)\mathrm{SFR}[t]-Z[t](1-R+\lambda)\mathrm{SFR}[t]+Z_\mathrm{in}\Phi[t]
\label{eq:dMZdt} \\
\end{split}
$$
Here $M_Z[t] = Z[t]\cdot M_\rm{g}[t]$ is the metal mass, y is the metal yield per stellar generation and $Z_\mathrm{in}$ is the metallicity of the inflow gas. Again, these 3 terms are metal mass changes due to gain from metal production, loss to star formation and gain from inflow. 

Now we can expand Eq. 2 by product rule to get
$$
\begin{split}
Z[t]\frac{dM_\mathrm{g}[t]}{dt}+M_\mathrm{g}[t]\frac{dZ[t]}{dt} = y(1-R)\mathrm{SFR}[t]-Z[t](1-R+\lambda)\mathrm{SFR}[t]+Z_\mathrm{in}\Phi[t] \\
\end{split}
$$
And substitute the Eq. 1 to equation above and simplify it to get $\frac{dZ[t]}{dt}$
$$
Z[t](\Phi[t] - (1-R+\lambda)\mathrm{SFR}[t])+M_\mathrm{g}[t]\frac{dZ[t]}{dt} = y(1-R)\mathrm{SFR}[t]-Z[t](1-R+\lambda)\mathrm{SFR}[t]+Z_\mathrm{in}\Phi[t] \\
M_\mathrm{g}[t]\frac{dZ[t]}{dt} = y(1-R)\mathrm{SFR}[t]+(Z_\mathrm{in}-Z[t])\Phi[t] \\
\frac{dZ[t]}{dt} = \frac{1}{M_\mathrm{g}[t]}(y(1-R)\mathrm{SFR}[t]+(Z_\mathrm{in}-Z[t])\Phi[t]) \\
\label{eq:dZdt}
$$
In the equilibrium state, metallicity $Z_\mathrm{eq}$ will be independent on time $t$, then we can force Eq. 4 to be 0 and rearange to get the expression of $Z_\mathrm{eq}$ in a function of $\mathrm{SFR}[t]$ and $\Phi[t]$:
$$
\begin{split}
\frac{dZ[t]}{dt} = \frac{1}{M_\mathrm{g}[t]}(y(1-R)\mathrm{SFR}[t]+(Z_\mathrm{in}-Z_\mathrm{eq})\Phi[t]) = 0 \\
y(1-R)\mathrm{SFR}[t]+(Z_\mathrm{in}-Z_\mathrm{eq})\Phi[t] = 0 \\
Z_\mathrm{eq} = Z_\mathrm{in} + y(1-R)\frac{\mathrm{SFR}[t]}{\Phi[t]} \\
\end{split}
$$
Since the inflow metallicity is much less than equilibrium metallicity ($Z_\mathrm{in}<<Z_\mathrm{eq}$), we can approximate the equation above as 
$$
\begin{split}
Z_\mathrm{eq} \approx y(1-R)\frac{\mathrm{SFR}[t]}{\Phi[t]} \\
\end{split}
$$
Now let consider a small perdubation of $Z_\mathrm{eq}$ (offset from the mean value). Since $y$ and $R$ are all constant here, we have 
$$
\begin{split}
\delta Z_\mathrm{eq} = y(1-R)\delta(\frac{\mathrm{SFR}[t]}{\Phi[t]})  \\
\label{eq:dZ}
\end{split}
$$
 Because we are interested in the offest in $12+\log_{10}(\mathrm{O/H})$, we have
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})_\mathrm{eq}) &= \delta\log_{10}(Z_\mathrm{eq})  \\
&= \frac{1}{\ln10}\frac{\delta Z_\mathrm{eq}}{Z_\mathrm{eq}} \\
&= \frac{1}{\ln10}\frac{y(1-R)\delta(\frac{\mathrm{SFR}[t]}{\Phi[t]})}{y(1-R)\frac{\mathrm{SFR}[t]}{\Phi[t]}} \\
&= \frac{1}{\ln10}\frac{\delta(\mathrm{SFR}[t]/\Phi[t])}{\mathrm{SFR}[t]/\Phi[t]} \\
&= \frac{1}{\ln10}\frac{\frac{\delta(\mathrm{SFR}[t]/\Phi[t])}{\delta\mathrm{SFR}[t]}\delta\mathrm{SFR}[t]+\frac{\delta(\mathrm{SFR}[t]/\Phi[t])}{\delta\Phi[t]}\delta\Phi[t]}{\mathrm{SFR}[t]/\Phi[t]} \\
&= \frac{1}{\ln10}\frac{\frac{1}{\Phi[t]}\delta\mathrm{SFR}[t]-\frac{\mathrm{SFR}[t]}{\Phi^2[t]}\delta\Phi[t]}{\mathrm{SFR}[t]/\Phi[t]} \\
&= \frac{1}{\ln10}(\frac{\delta\mathrm{SFR}[t]}{\mathrm{SFR}[t]}-\frac{\delta\Phi[t]}{\Phi[t]}) \\
&= \frac{1}{\ln10} (\delta\ln(\mathrm{SFR}[t])-\delta\ln(\Phi[t])) \\
&= \delta\log_{10}(\mathrm{SFR}[t])-\delta\log_{10}(\Phi[t])  \\ 
\end{split} \\
$$
Thus, we can see that the sign of $\delta(12+\log_{10}(\mathrm{O/H})_\mathrm{eq})$ is only dertermined by relative change of $\log_{10}(\mathrm{SFR}[t])$ and $\log_{10}(\Phi[t])$. 

More specific, by $\mathrm{SFR}[t]=\epsilon[t]M_\mathrm{g}[t]$, where $\epsilon[t]$ is the star formation efficiency, we have 
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})_\mathrm{eq}) = \delta\log_{10}(\epsilon[t])+\delta\log_{10}(M_\mathrm{g}[t])-\delta\log_{10}(\Phi[t])  \\ 
\end{split}
$$
Therefore, this leaves us 3 cases.

1. In the case that $\delta\log_{10}(\epsilon[t])$ dominants, if $\delta\log_{10}(\epsilon[t])$ increases, both $\delta(12+\log_{10}(\mathrm{O/H})_\mathrm{eq})$ and $\delta\log_{10}(\mathrm{SFR}[t])$ increase immediately, so we can observe the positive correlation between them. 
2. However, if $\delta\log_{10}(\Phi[t])$ dominates, even though $\delta\log_{10}(\mathrm{SFR}[t])$ is still positive, due to the minus sign and relatively large change in $\delta\log_{10}(\Phi[t])$, we can observe negative correlation. 
3. If the change of $\log_{10}(\mathrm{SFR}[t])$ and $\log_{10}(\Phi[t])$ is comparable, that lead to nearly stationary $12+\log_{10}(\mathrm{O/H})$, which show no correlation. In this scenario, it requires $\Phi[t]$ is linearly proportional to $\mathrm{SFR}[t]$. 