## 20251027 Numerical Bathtub

## 1. Model build-up

Let's consider solving the time-dependent regulator model numerically. We start with writing down the set of ODEs and note that I replace all the parameters into surface density to describe the spatially-resolved scenario.

First, as usual, we define the SFR surface density in two ways: 
$$
\begin{split}
\Sigma_\mathrm{SFR}[t] = \epsilon[t]\Sigma_\mathrm{g}[t] = \frac{d}{dt}\Sigma_*[t]
\end{split}
$$
Then, the change of gas mass surface density:
$$
\begin{split}
\frac{d}{dt}\Sigma_\mathrm{g}[t] = \Sigma_\mathrm{in}[t]-(1-R+\lambda)\epsilon[t]\Sigma_\mathrm{g}[t]
\end{split}
$$
Next, the change of metal mass surface density:
$$
\begin{split}
\frac{d}{dt}\Sigma_Z[t] = y(1-R)\epsilon[t]\Sigma_\mathrm{g}[t]-(1-R+\lambda)\epsilon[t]\Sigma_\mathrm{g}[t]Z[t]+\Sigma_\mathrm{in}[t]Z_\mathrm{in}
\end{split}
$$
And, we also define the metallicity as:
$$
\begin{split}
Z[t]=\frac{\Sigma_Z[t]}{\Sigma_\mathrm{g}[t]}
\end{split}
$$
Also, the sSFR:
$$
\begin{split}
\mathrm{sSFR} = \frac{\Sigma_\mathrm{SFR}}{\Sigma_*}
\end{split}
$$
Finally, we can convert metallicity to observable with solar abundace at 8.69 and $Z_\mathrm{surface} = 0.0139$ (Asplund+2021):
$$
\begin{split}
12+\log_{10}(\mathrm{O/H})[t] = 8.69+\log_{10}(\frac{Z[t]}{0.0139})
\end{split}
$$

## 2. Equilibrium state

Taking the product rule of Eq. 3 and substituting Eq.1 to it, we get $\frac{dZ[t]}{dt}$
$$
Z[t](\Sigma_\mathrm{in}[t] - (1-R+\lambda)\Sigma_\mathrm{SFR}[t])+\Sigma_\mathrm{g}[t]\frac{dZ[t]}{dt} = y(1-R)\Sigma_\mathrm{SFR}[t]-Z[t](1-R+\lambda)\Sigma_\mathrm{SFR}[t]+Z_\mathrm{in}\Sigma_\mathrm{in}[t] \\
\Sigma_\mathrm{g}[t]\frac{dZ[t]}{dt} = y(1-R)\Sigma_\mathrm{SFR}[t]+(Z_\mathrm{in}-Z[t])\Sigma_\mathrm{in}[t] \\
\frac{dZ[t]}{dt} = \frac{1}{\Sigma_\mathrm{g}[t]}(y(1-R)\Sigma_\mathrm{SFR}[t]+(Z_\mathrm{in}-Z[t])\Sigma_\mathrm{in}[t]) \\
\label{eq:dZdt}
$$
Now let's consider the equilibrium state that the model rapidly self-regulates to a steady state, which means $\frac{dZ[t]}{dt} = 0$ in a short time scale (compared to the total depletion time). Thus, such metallicity at equilibrium state ($Z_\mathrm{eq}$) are a series of constant that is time-independent in a short period of time (but it is still time-dependent in a much longer period, so it is still written as $Z_\mathrm{eq}[t]$). And therefore it requires
$$
\begin{split}
\frac{dZ[t]}{dt} = \frac{1}{\Sigma_\mathrm{g}[t]}(y(1-R)\Sigma_\mathrm{SFR}[t]+(Z_\mathrm{in}-Z_\mathrm{eq}[t])\Sigma_\mathrm{in}[t]) = 0 \\
y(1-R)\Sigma_\mathrm{SFR}[t]+(Z_\mathrm{in}-Z_\mathrm{eq}[t])\Sigma_\mathrm{in}[t] = 0 \\
Z_\mathrm{eq}[t] = Z_\mathrm{in} + y(1-R)\frac{\Sigma_\mathrm{SFR}[t]}{\Sigma_\mathrm{in}[t]} \\
\end{split}
$$
Again, we are interested in the metallicity that satisfy $Z[t] = Z_\mathrm{in} + y(1-R)\frac{\Sigma_\mathrm{SFR}[t]}{\Sigma_\mathrm{in}[t]}$, at any given $t$. 

Thus, with this set of $Z_\mathrm{eq}$, we can find the mean value ($<Z_\mathrm{eq}>$) and then consider a small perdubation (offset from $<Z_\mathrm{eq}>$). Since $y$ and $R$ are all constant here, we have 
$$
\begin{split}
\delta Z_\mathrm{eq}[t] = y(1-R)\delta(\frac{\Sigma_\mathrm{SFR}[t]}{\Sigma_\mathrm{in}[t]})  \\
\label{eq:dZ}
\end{split}
$$
 Because we are interested in the offest in $12+\log_{10}(\mathrm{O/H})$, we have
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})_\mathrm{eq}[t]) &= \delta\log_{10}(Z_\mathrm{eq}[t])  \\
&= \frac{1}{\ln10}\frac{\delta Z_\mathrm{eq}[t]}{Z_\mathrm{eq}[t]} \\
&= \frac{1}{\ln10}\frac{y(1-R)\delta(\frac{\Sigma_\mathrm{SFR}[t]}{\Sigma_\mathrm{in}[t]})}{Z_\mathrm{in}+y(1-R)\frac{\Sigma_\mathrm{SFR}[t]}{\Sigma_\mathrm{in}[t]}} \\
&= \frac{y(1-R)\frac{\Sigma_\mathrm{SFR}[t]}{\Sigma_\mathrm{in}[t]}}{Z_\mathrm{in}+y(1-R)\frac{\Sigma_\mathrm{SFR}[t]}{\Sigma_\mathrm{in}[t]}}(\delta\log_{10}(\Sigma_\mathrm{SFR}[t])-\delta\log_{10}(\Sigma_\mathrm{in}[t])) \\
&= \frac{Z_\mathrm{eq}[t]-Z_\mathrm{in}}{Z_\mathrm{eq}[t]}(\delta\log_{10}(\Sigma_\mathrm{SFR}[t])-\delta\log_{10}(\Sigma_\mathrm{in}[t])) \\
&\propto \delta\log_{10}(\Sigma_\mathrm{SFR}[t])-\delta\log_{10}(\Sigma_\mathrm{in}[t]) \\
\end{split} \\
$$
Thus, we can see that the sign of $\delta(12+\log_{10}(\mathrm{O/H})_\mathrm{eq})$ is only dertermined by relative change of $\log_{10}(\Sigma_\mathrm{SFR}[t])$ and $\log_{10}(\Sigma_\mathrm{in}[t])$. 

More specific, by $\Sigma_\mathrm{SFR}[t]=\epsilon[t]\Sigma_\mathrm{g}[t]$, where $\epsilon[t]$ is the star formation efficiency, we have 
$$
\begin{split}
\delta(12+\log_{10}(\mathrm{O/H})_\mathrm{eq}) \propto \delta\log_{10}(\epsilon[t])+\delta\log_{10}(\Sigma_\mathrm{g}[t])-\delta\log_{10}(\Sigma_\mathrm{in}[t])  \\ 
\end{split}
$$
Therefore, this leaves us 3 cases. Note that $\delta$ in $\log$ scale of a quantity indicates the relative change of the quantity, and the sign indicates increase or decrease. 

1. In the case that $\delta\log_{10}(\epsilon[t])$ dominants, if $\delta\log_{10}(\epsilon[t])$ increases, both $\delta(12+\log_{10}(\mathrm{O/H})_\mathrm{eq})$ and $\delta\log_{10}(\Sigma_\mathrm{SFR}[t])$ increase immediately, so we can observe the positive correlation between them. 
2. However, if $\delta\log_{10}(\Sigma_\mathrm{in}[t])$ dominates, even though $\delta\log_{10}(\Sigma_\mathrm{SFR}[t])$ is still positive, due to the minus sign and relatively large change in $\delta\log_{10}(\Phi[t])$, we can observe negative correlation. 
3. If the change of $\log_{10}(\Sigma_\mathrm{SFR}[t])$ and $\log_{10}(\Sigma_\mathrm{in}[t])$ is comparable, that lead to nearly stationary $12+\log_{10}(\mathrm{O/H})$, which show no correlation. In this scenario, it requires $\Phi[t]$ is linearly proportional to $\mathrm{SFR}[t]$. 

## 3. Constants and initial values

Here, return fraction, metal yield per stellar, mass loading factor and inflow metallicity are constants:
$$
\left\{
\begin{split}
R &= 0.5 \\
y &= 0.01 \\
\lambda &= 1.5 \\
Z_\mathrm{in} &= 0.1 \\
\end{split}
\right.
$$
Now let's say at $t=0\mathrm{Myr}$,  we have 
$$
\left\{
\begin{split}
\epsilon_0 &= \epsilon[t=0] = \frac{\Sigma_\mathrm{SFR,0}}{\Sigma_\mathrm{g,0}} = 10^{-9}/\mathrm{yr} \\
\Sigma_\mathrm{in,0} &= \Sigma_\mathrm{in}[t=0] = 10^{-3}\mathrm{M_\odot/yr/kpc^2} \\
Z_0 &= Z[t=0] = 0.5 \\
\Sigma_{*,0} &= \Sigma_*[t=0] = 10^7\mathrm{M_\odot/kpc^2} \\
\Sigma_\mathrm{g,0} &= \Sigma_\mathrm{g}[t=0] = 10^6\mathrm{M_\odot/kpc^2} \\
\Sigma_\mathrm{SFR,0} &= \Sigma_\mathrm{SFR}[t=0] = 10^{-3}\mathrm{M_\odot/yr/kpc^2} \\
\mathrm{sSFR}_0 &= \mathrm{sSFR}[t=0] = \frac{\Sigma_\mathrm{SFR,0}}{\Sigma_{*,0}} = 10^{-10}/\mathrm{yr} \\
\Sigma_{Z,0} &= \Sigma_{Z}[t=0] = \Sigma_\mathrm{in,0}\Sigma_{Z,0} = 0.5\times10^6\mathrm{M_\odot/kpc^2} \\
12+\log_{10}(\mathrm{O/H})_0 &= 12+\log_{10}(\mathrm{O/H})[t=0] = 8.69+\log_{10}(\frac{Z_0}{0.0139}) \\
\end{split}
\right.
$$

## 4. Numerical forms of $\epsilon[t]$ and $\Sigma_\mathrm{in}[t]$

We consider the following cases:

### 4.1 Constant $\epsilon[t]$ and constant $\Sigma_\mathrm{in}[t]$

A simple case that both are constant:
$$
\left\{
\begin{split}
\epsilon[t] &= \epsilon_0 = 10^{-9}/\mathrm{yr} \\
\Sigma_\mathrm{in}[t] &= \Sigma_\mathrm{in,0} = 10^{-3}\mathrm{M_\odot/yr/kpc^2} \\
\end{split}
\right.
$$

### 4.2 Sinusoidal variation of $\epsilon[t]$ and constant $\Sigma_\mathrm{in}[t]$

Here we define a constant $A=0.1$ here to denote the amplitide of the sinusoidal variation compared to initial value, and $T_p = 1000$ (in the unit of Myr) to indicate the variation period:
$$
\left\{
\begin{split}
\epsilon[t] &= \epsilon_0 + A\epsilon_0\sin(\frac{2\pi t}{T_p}) \\
\Sigma_\mathrm{in}[t] &= \Sigma_\mathrm{in,0} = 10^{-3}\mathrm{M_\odot/yr/kpc^2} \\
\end{split}
\right.
$$

### 4.3 Constant $\epsilon[t]$ and sinusoidal variation of $\Sigma_\mathrm{in}[t]$

Similarly, we have:
$$
\left\{
\begin{split}
\epsilon[t] &= \epsilon_0 = 10^{-9}/\mathrm{yr} \\
\Sigma_\mathrm{in}[t] &= \Sigma_\mathrm{in,0} + A\Sigma_\mathrm{in,0}\sin(\frac{2\pi t}{T_p}) \\
\end{split}
\right.
$$

## 5. Evolution

We consider a total period of $T = 10\mathrm{Gyr}$, with a timestep of $\Delta t = 0.1\mathrm{Myr}$. Thus, for any given moment $t_i \in [0,T]$, the iteration in the next timestamp ($t_{i+1} = t_i+\Delta t$) will be 
$$
\left\{
\begin{split}
\epsilon[t_{i+1}] &= \epsilon[t_{i+1}] \\
\Sigma_\mathrm{in}[t_{i+1}] &= \Sigma_\mathrm{in}[t_{i+1}] \\
\Sigma_\mathrm{g}[t_{i+1}] &= \Sigma_\mathrm{g}[t_i]+\frac{d}{dt}\Sigma_\mathrm{g}[t_i]\Delta t \\ 
&= \Sigma_\mathrm{g}[t_i]+(\Sigma_\mathrm{in}[t_i]-(1-R+\lambda)\epsilon[t_i]\Sigma_\mathrm{g}[t_i])\Delta t \\
\Sigma_Z[t_{i+1}] &= \Sigma_\mathrm{Z}[t_i]+\frac{d}{dt}\Sigma_\mathrm{Z}[t_i]\Delta t \\ 
&= \Sigma_\mathrm{Z}[t_i]+(y(1-R)\epsilon[t_i]\Sigma_\mathrm{g}[t_i]-(1-R+\lambda)\epsilon[t_i]\Sigma_\mathrm{g}[t_i]Z[t_i]+\Sigma_\mathrm{in}[t_i]Z_\mathrm{in})\Delta t \\ 
\Sigma_\mathrm{SFR}[t_{i+1}] &= \epsilon[t_{i+1}]\Sigma_\mathrm{g}[t_{i+1}] \\
Z[t_{i+1}] &= \frac{\Sigma_Z[t_{i+1}]}{\Sigma_\mathrm{g}[t_{i+1}]} \\
\Sigma_*[t_{i+1}] &= \Sigma_*[t_{i}]+\Sigma_\mathrm{SFR}[t_i]\Delta t \\
\mathrm{sSFR}[t_{i+1}] &= \frac{\Sigma_\mathrm{SFR}[t_{i+1}]}{\Sigma_*[t_{i+1}]} \\
12+\log_{10}(\mathrm{O/H})[t_{i+1}] &= 8.69+\log_{10}(\frac{Z[t_{i+1}]}{0.0139}) \\
\end{split}
\right.
$$

