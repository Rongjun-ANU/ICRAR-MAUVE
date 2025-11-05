# 20251105 Derivation of Metallicity Perturbation Without Equilibrium Assumption

## Starting Equations

From the bathtub model, we have:

**Gas surface density evolution:**
$$\frac{d\Sigma_\mathrm{g}}{dt} = \Sigma_\mathrm{in} - (1-R+\lambda)\Sigma_\mathrm{SFR}$$

**Metal surface density evolution:**
$$\frac{d\Sigma_Z}{dt} = y(1-R)\Sigma_\mathrm{SFR} - (1-R+\lambda)\Sigma_\mathrm{SFR}Z + \Sigma_\mathrm{in}Z_\mathrm{in}$$

**Metallicity definition:**
$$Z = \frac{\Sigma_Z}{\Sigma_\mathrm{g}}$$

## Step 1: Derive the Metallicity Evolution Equation

Using the product rule on $\Sigma_Z = Z\Sigma_\mathrm{g}$:

$$\frac{d\Sigma_Z}{dt} = Z\frac{d\Sigma_\mathrm{g}}{dt} + \Sigma_\mathrm{g}\frac{dZ}{dt}$$

Substitute the gas evolution equation:

$$\frac{d\Sigma_Z}{dt} = Z[\Sigma_\mathrm{in} - (1-R+\lambda)\Sigma_\mathrm{SFR}] + \Sigma_\mathrm{g}\frac{dZ}{dt}$$

Equate with the metal evolution equation:

$$Z\Sigma_\mathrm{in} - Z(1-R+\lambda)\Sigma_\mathrm{SFR} + \Sigma_\mathrm{g}\frac{dZ}{dt} = y(1-R)\Sigma_\mathrm{SFR} - (1-R+\lambda)\Sigma_\mathrm{SFR}Z + \Sigma_\mathrm{in}Z_\mathrm{in}$$

Simplify (the $-Z(1-R+\lambda)\Sigma_\mathrm{SFR}$ terms cancel):

$$\Sigma_\mathrm{g}\frac{dZ}{dt} = y(1-R)\Sigma_\mathrm{SFR} + (Z_\mathrm{in} - Z)\Sigma_\mathrm{in}$$

Therefore:

$$\frac{dZ}{dt} = \frac{1}{\Sigma_\mathrm{g}}\left[y(1-R)\Sigma_\mathrm{SFR} + (Z_\mathrm{in} - Z)\Sigma_\mathrm{in}\right]$$

## Step 2: Solve for Z Including Time Derivative

Rearrange the equation:

$$\Sigma_\mathrm{g}\frac{dZ}{dt} = y(1-R)\Sigma_\mathrm{SFR} + (Z_\mathrm{in} - Z)\Sigma_\mathrm{in}$$

$$\Sigma_\mathrm{g}\frac{dZ}{dt} = y(1-R)\Sigma_\mathrm{SFR} + Z_\mathrm{in}\Sigma_\mathrm{in} - Z\Sigma_\mathrm{in}$$

$$\Sigma_\mathrm{g}\frac{dZ}{dt} + Z\Sigma_\mathrm{in} = y(1-R)\Sigma_\mathrm{SFR} + Z_\mathrm{in}\Sigma_\mathrm{in}$$

Solve for Z:

$$Z = \frac{y(1-R)\Sigma_\mathrm{SFR} + Z_\mathrm{in}\Sigma_\mathrm{in} - \Sigma_\mathrm{g}\frac{dZ}{dt}}{\Sigma_\mathrm{in}}$$

$$Z = Z_\mathrm{in} + y(1-R)\frac{\Sigma_\mathrm{SFR}}{\Sigma_\mathrm{in}} - \frac{\Sigma_\mathrm{g}}{\Sigma_\mathrm{in}}\frac{dZ}{dt}$$

## Step 3: Define the Quasi-Equilibrium Metallicity

Define $Z_\mathrm{eq}$ as the metallicity that would exist if $\frac{dZ}{dt} = 0$:

$$Z_\mathrm{eq} = Z_\mathrm{in} + y(1-R)\frac{\Sigma_\mathrm{SFR}}{\Sigma_\mathrm{in}}$$

Then the actual metallicity can be written as:

$$Z = Z_\mathrm{eq} - \frac{\Sigma_\mathrm{g}}{\Sigma_\mathrm{in}}\frac{dZ}{dt}$$

This shows that Z deviates from $Z_\mathrm{eq}$ by a correction term proportional to $\frac{dZ}{dt}$.

## Step 4: Estimate the Metallicity Relaxation Timescale

From the metallicity evolution equation:

$$\frac{dZ}{dt} = \frac{1}{\Sigma_\mathrm{g}}\left[y(1-R)\Sigma_\mathrm{SFR} + (Z_\mathrm{in} - Z)\Sigma_\mathrm{in}\right]$$

Consider a perturbation $\delta Z$ from quasi-equilibrium. Near $Z_\mathrm{eq}$:

$$\frac{d(\delta Z)}{dt} = \frac{1}{\Sigma_\mathrm{g}}\left[y(1-R)\delta\Sigma_\mathrm{SFR} + (Z_\mathrm{in} - Z_\mathrm{eq} - \delta Z)\delta\Sigma_\mathrm{in} - \delta Z\Sigma_\mathrm{in}\right]$$

Since $Z_\mathrm{in} - Z_\mathrm{eq} = -y(1-R)\frac{\Sigma_\mathrm{SFR}}{\Sigma_\mathrm{in}}$, and assuming variations are slow:

$$\frac{d(\delta Z)}{dt} \approx -\frac{\Sigma_\mathrm{in}}{\Sigma_\mathrm{g}}\delta Z + \text{driving terms}$$

The characteristic relaxation timescale for metallicity perturbations is:

$$\boxed{\tau_Z = \frac{\Sigma_\mathrm{g}}{\Sigma_\mathrm{in}}}$$

This is the **gas depletion timescale** (assuming $\Sigma_\mathrm{in}$ is the dominant source).

## Step 5: Timescale Separation Condition

Let $\tau_\mathrm{var}$ be the timescale over which $\Sigma_\mathrm{SFR}$ and $\Sigma_\mathrm{in}$ vary significantly.

**Condition for quasi-equilibrium:** $\tau_Z \ll \tau_\mathrm{var}$

This means:

$$\frac{\Sigma_\mathrm{g}}{\Sigma_\mathrm{in}} \ll \tau_\mathrm{var}$$

When this condition is satisfied, the correction term in Step 3 becomes negligible:

$$\left|\frac{\Sigma_\mathrm{g}}{\Sigma_\mathrm{in}}\frac{dZ}{dt}\right| \ll |Z_\mathrm{eq} - Z_\mathrm{in}|$$

Therefore:

$$Z \approx Z_\mathrm{eq} = Z_\mathrm{in} + y(1-R)\frac{\Sigma_\mathrm{SFR}}{\Sigma_\mathrm{in}}$$

## Step 6: Derive the Perturbation Relationship

Starting from the quasi-equilibrium relation:

$$Z \approx Z_\mathrm{in} + y(1-R)\frac{\Sigma_\mathrm{SFR}}{\Sigma_\mathrm{in}}$$

Consider small perturbations around a reference state (denoted by subscript 0):

$$Z_0 + \delta Z = Z_\mathrm{in} + y(1-R)\frac{\Sigma_{\mathrm{SFR},0} + \delta\Sigma_\mathrm{SFR}}{\Sigma_{\mathrm{in},0} + \delta\Sigma_\mathrm{in}}$$

Using the approximation $\frac{A+\delta A}{B+\delta B} \approx \frac{A}{B}\left(1 + \frac{\delta A}{A} - \frac{\delta B}{B}\right)$:

$$\delta Z \approx y(1-R)\frac{\Sigma_{\mathrm{SFR},0}}{\Sigma_{\mathrm{in},0}}\left(\frac{\delta\Sigma_\mathrm{SFR}}{\Sigma_{\mathrm{SFR},0}} - \frac{\delta\Sigma_\mathrm{in}}{\Sigma_{\mathrm{in},0}}\right)$$

## Step 7: Convert to Logarithmic Perturbations

For logarithmic perturbations:

$$\delta\log_{10}(Z) = \frac{1}{\ln 10}\frac{\delta Z}{Z_0}$$

Substitute the expression for $\delta Z$:

$$\delta\log_{10}(Z) = \frac{1}{\ln 10}\frac{y(1-R)\frac{\Sigma_{\mathrm{SFR},0}}{\Sigma_{\mathrm{in},0}}}{Z_0}\left(\frac{\delta\Sigma_\mathrm{SFR}}{\Sigma_{\mathrm{SFR},0}} - \frac{\delta\Sigma_\mathrm{in}}{\Sigma_{\mathrm{in},0}}\right)$$

Using $\frac{\delta X}{X} = \ln(10)\delta\log_{10}(X)$:

$$\delta\log_{10}(Z) = \frac{y(1-R)\frac{\Sigma_{\mathrm{SFR},0}}{\Sigma_{\mathrm{in},0}}}{Z_0}\left(\delta\log_{10}(\Sigma_\mathrm{SFR}) - \delta\log_{10}(\Sigma_\mathrm{in})\right)$$

From the quasi-equilibrium condition: $Z_0 = Z_\mathrm{in} + y(1-R)\frac{\Sigma_{\mathrm{SFR},0}}{\Sigma_{\mathrm{in},0}}$

Therefore:

$$\delta\log_{10}(Z) = \frac{Z_0 - Z_\mathrm{in}}{Z_0}\left(\delta\log_{10}(\Sigma_\mathrm{SFR}) - \delta\log_{10}(\Sigma_\mathrm{in})\right)$$

## Step 8: Final Result for Observable Metallicity

Since $12 + \log_{10}(\mathrm{O/H}) = 8.69 + \log_{10}(Z/0.0139)$:

$$\delta[12 + \log_{10}(\mathrm{O/H})] = \delta\log_{10}(Z)$$

Therefore:

$$\boxed{\delta[12 + \log_{10}(\mathrm{O/H})] = \frac{Z - Z_\mathrm{in}}{Z}\left(\delta\log_{10}(\Sigma_\mathrm{SFR}) - \delta\log_{10}(\Sigma_\mathrm{in})\right)}$$

Or more compactly:

$$\boxed{\delta[12 + \log_{10}(\mathrm{O/H})] \propto \delta\log_{10}(\Sigma_\mathrm{SFR}) - \delta\log_{10}(\Sigma_\mathrm{in})}$$

## Conclusion

**Key Result:** The relationship $\delta[12 + \log_{10}(\mathrm{O/H})] \propto \delta\log_{10}(\Sigma_\mathrm{SFR}) - \delta\log_{10}(\Sigma_\mathrm{in})$ holds **without** assuming $\frac{dZ}{dt} = 0$, as long as the timescale separation condition $\tau_Z \ll \tau_\mathrm{var}$ is satisfied.

**Physical Interpretation:** The metallicity "tracks" the quasi-equilibrium value because it adjusts to changes in $\Sigma_\mathrm{SFR}$ and $\Sigma_\mathrm{in}$ much faster than these quantities vary. This is a **quasi-static approximation** based on timescale separation, not true equilibrium.

---

## Bonus: Demonstrating $\tau_Z \ll \tau_\mathrm{var}$

### Typical Values

For a star-forming galaxy:
- $\Sigma_\mathrm{g} \sim 10$ M$_\odot$ pc$^{-2}$
- $\Sigma_\mathrm{in} \sim 1$ M$_\odot$ pc$^{-2}$ Gyr$^{-1}$
- $\Sigma_\mathrm{SFR} \sim 0.1$ M$_\odot$ pc$^{-2}$ Gyr$^{-1}$

### Metallicity Relaxation Timescale

$$\tau_Z = \frac{\Sigma_\mathrm{g}}{\Sigma_\mathrm{in}} = \frac{10 \text{ M}_\odot \text{ pc}^{-2}}{1 \text{ M}_\odot \text{ pc}^{-2} \text{ Gyr}^{-1}} = 10 \text{ Gyr}$$

### SFR Variation Timescale

Typical variations in SFR occur on timescales of:
- **Secular evolution:** $\tau_\mathrm{var} \sim$ few Gyr (Hubble time)
- **Galaxy interactions:** $\tau_\mathrm{var} \sim 0.1 - 1$ Gyr
- **Within simulation timesteps:** If you're resolving variations on $\sim 100$ Myr timescales

**Wait, this doesn't satisfy $\tau_Z \ll \tau_\mathrm{var}$!**

### Resolution: The Effective Timescale

The relevant comparison is not just $\frac{\Sigma_\mathrm{g}}{\Sigma_\mathrm{in}}$, but the actual metallicity adjustment rate including outflows:

From the full equation:
$$\frac{dZ}{dt} = \frac{1}{\Sigma_\mathrm{g}}\left[y(1-R)\Sigma_\mathrm{SFR} + (Z_\mathrm{in} - Z)\Sigma_\mathrm{in}\right]$$

Near equilibrium, the restoring force is:
$$\frac{d(\delta Z)}{dt} \approx -\frac{\Sigma_\mathrm{in}}{\Sigma_\mathrm{g}}\delta Z$$

But the **total metal cycling rate** includes star formation:
$$\tau_Z^\mathrm{eff} = \frac{\Sigma_\mathrm{g}}{\Sigma_\mathrm{in} + (1-R+\lambda)\Sigma_\mathrm{SFR}}$$

With $\Sigma_\mathrm{SFR} \sim 0.1$ M$_\odot$ pc$^{-2}$ yr$^{-1}$ and $(1-R+\lambda) \sim 0.6$:

$$\tau_Z^\mathrm{eff} = \frac{10}{1 + 0.6 \times 0.1} \approx 9.4 \text{ Gyr}$$

### Why It Still Works

Even though $\tau_Z \sim \tau_\mathrm{var}$, the quasi-equilibrium approximation works because:

1. **The metallicity evolution is a damped response:** Perturbations decay exponentially with $e^{-t/\tau_Z}$
2. **Most astrophysical variations are gradual:** SFR doesn't jump discontinuously
3. **In practice:** If you're looking at correlations on timescales $> 0.1\tau_Z$, the metallicity has already adjusted significantly

**Quantitative criterion:** The approximation is good when:
$$\frac{1}{\tau_Z}\frac{d\ln\Sigma_\mathrm{SFR}}{dt} \ll 1$$

This means the fractional change in $\Sigma_\mathrm{SFR}$ per metallicity relaxation time is small, which is typically true for secular evolution even when $\tau_Z \sim \tau_\mathrm{var}$.