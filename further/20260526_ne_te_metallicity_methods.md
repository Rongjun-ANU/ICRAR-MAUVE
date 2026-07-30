# 20260526 Notes: NE_SII and Te-based metallicity products in SFR+Z.py

This note documents the current implementation in `SFR+Z.py` for the PyNeb
electron-density product and the three Te-based oxygen-abundance products.  The
goal is to make the CANFAR products auditable: every finite value should be
traceable to a specific detection gate, line ratio, PyNeb call, temperature
assumption, and error-propagation step.

The three metallicity branches documented here are:

1. `O_H_BR24_DIRECT_{HII,SF}`: Brazzini et al. 2024 multi-zone direct-style
   method using measured Te([N II]) and measured Te([S III]).
2. `O_H_NII_OII7325_{HII,SF}`: MAUVE semi-direct ionic-sum method using
   measured Te([N II]), red [O II] 7325, and an inferred high-zone temperature.
3. `O_H_NII_K25_{HII,SF}`: Kreckel et al. 2025 / Mendez-Delgado et al. 2023
   single-temperature Te([N II]) metallicity proxy.

All new abundance-like FITS maps use `BUNIT = 1`, not `12+log(O/H)`, so that
CARTA/CASA treats them as dimensionless images rather than warning about an
unknown FITS unit.

## 1. Region masks and detection gates

The script first applies the nGIST spatial mask when
`v3tk_v7.6.8/{GALID}/{GALID}_mask.fits` exists.  In the current MAUVE v3tk
products, the script interprets:

```text
MASK == 0  -> keep / usable spatial pixel
MASK == 1  -> masked spatial pixel
```

The mask is applied to `V_STARS2`, `SIGMA_STARS2`, and all raw gas-line flux
and flux-error maps.  Therefore later checks such as `np.isfinite(V_STARS2)`
implicitly respect the spatial mask.

For Te products, an emission line is considered detected only if the raw
gas-line maps satisfy:

```math
\mathrm{detected}(F_\lambda) =
\mathrm{finite}(F_\lambda)
\land \mathrm{finite}(\sigma_\lambda)
\land \sigma_\lambda > 0
\land F_\lambda / \sigma_\lambda \ge 3
\land F_\lambda \ge 20 .
```

The two numeric thresholds are the current script-level values:

```text
cut   = 3
noise = 20
```

The Te products are computed twice:

- `_HII`: conservative Classification-2 mask, where SF means HII only.
- `_SF`: Classification-1 mask, where the existing script treats SF as
  HII+Composite in the [N II] BPT branch and HII in the [S II] BPT branch.

So if an HDU exists as `_HII`, the matching `_SF` HDU should also exist.

## 2. Dust correction used by all Te products

### 2.1 Balmer decrement

The observed Balmer decrement is:

```math
BD = \frac{F_{\mathrm{H}\alpha}}{F_{\mathrm{H}\beta}} .
```

If `BD` is non-finite in a spatially usable pixel, the script sets it to the
intrinsic value.  If the corrected `BD` is below the intrinsic value, the script
also floors it to the intrinsic value:

```math
BD_{\mathrm{used}} = \max(BD, 2.86).
```

Thus the pipeline does not create negative gas reddening from
`Halpha/Hbeta < 2.86`.

### 2.2 Gas E(B-V)

The local default extinction curve is Calzetti et al. 2000, matching the
original MAUVE SFR+Z products.  The code still contains an optional Milky Way
CCM89 branch for explicit tests, but it is not the default.  The current default
is:

```python
extinction_law = "calzetti"
mw_rv = 3.1  # used only if extinction_law is set to "mw"/"ccm89"
```

For each line, the script computes:

```math
k(\lambda) = \frac{A(\lambda)}{E(B-V)} .
```

Then the Balmer-decrement reddening is:

```math
E(B-V)_{\mathrm{BD}} =
\frac{2.5}{k(\mathrm{H}\beta)-k(\mathrm{H}\alpha)}
\log_{10}
\left(
\frac{BD_{\mathrm{used}}}{2.86}
\right).
```

This follows the standard Balmer-decrement form used in nebular work.  Brazzini
et al. 2024 and Kreckel et al. 2025 use an intrinsic Halpha/Hbeta ratio of 2.86
for their PHANGS-MUSE dust corrections, with an O'Donnell/CCM-family Milky Way
curve.  The local script deliberately keeps the original MAUVE convention and
uses the Calzetti et al. 2000 curve by default for consistency with the previous
SFR+Z products.

### 2.3 Corrected fluxes and errors

Each line flux is dereddened as:

```math
I_\lambda =
F_\lambda\,10^{0.4\,k(\lambda)\,E(B-V)_{\mathrm{BD}}}.
```

The propagated uncertainty includes both the line-flux error and the
Balmer-decrement reddening error:

```math
\sigma(I_\lambda)^2 =
\left[
\sigma(F_\lambda)\,10^{0.4\,k(\lambda)\,E(B-V)}
\right]^2
+
\left[
|I_\lambda|\,\ln(10)\,0.4\,|k(\lambda)|\,\sigma(E(B-V))
\right]^2 .
```

The script now writes corrected flux and corrected error maps for every
available gas line, with global, `_HII`, and `_SF` variants.

## 3. Electron density from [S II] 6716/6731

### 3.1 Required lines

The density map requires both [S II] lines to pass the detection gate:

```text
SII6716_FLUX and SII6716_FLUX_ERR
SII6730_FLUX and SII6730_FLUX_ERR
```

The valid-density mask is:

```math
M_{n_e} =
\mathrm{detected}([\mathrm{S\,II}]6716)
\land
\mathrm{detected}([\mathrm{S\,II}]6731).
```

### 3.2 Density ratio

The density-sensitive ratio is computed from dust-corrected fluxes:

```math
R_{\mathrm{SII}} =
\frac{I(6716)}{I(6731)} .
```

The line-ratio uncertainty is the standard independent-error ratio formula:

```math
\sigma(R_{\mathrm{SII}}) =
|R_{\mathrm{SII}}|
\sqrt{
\left(\frac{\sigma_{6716}}{I_{6716}}\right)^2
+
\left(\frac{\sigma_{6731}}{I_{6731}}\right)^2
}.
```

### 3.3 PyNeb inversion

The physical density is obtained by inverting the PyNeb emissivity ratio:

```math
R_{\mathrm{SII}} =
\frac{\epsilon_{6716}(T_e,n_e)}
{\epsilon_{6731}(T_e,n_e)} .
```

The current product assumes:

```text
Te = 10000 K
atom = PyNeb Atom("S", 2)
ratio = I(6716)/I(6731)
```

This gives the main HDU:

```text
NE_SII
```

with units:

```text
cm-3
```

### 3.4 Lookup table and exact high-density fallback

For speed, the script does not call `getTemDen` independently for every normal
pixel.  It creates or loads a reusable PyNeb lookup table in the current working
directory:

```text
pyneb_sii_6716_6731_density_lookup_te8000_13000_step1000_ne20_100000_n4096.npz
pyneb_sii_6716_6731_density_lookup_te8000_13000_step1000_ne20_100000_n4096.png
```

The cache stores:

```text
temperature_grid = 8000, 9000, ..., 13000 K
density_grid     = geomspace(20, 100000, 4096) cm-3
ratio_grid_2d    = I(6716)/I(6731) from PyNeb emissivities
```

The density products currently use the 10000 K row of this table.  The grid is
logarithmic in density.  With 4096 points from 20 to 100000 cm-3, the fractional
density spacing is about:

```math
\exp\left[
\frac{\ln(100000/20)}{4095}
\right] - 1
\approx 0.0021 ,
```

or about 0.2 percent per lookup step.  This is much finer than the astrophysical
uncertainty implied by the line-ratio errors.

The lookup is deliberately conservative:

- Ratios at or above the low-density limit are assigned the low-density floor.
- Ratios outside the high-density side of the cached range are not silently
  extrapolated.
- Any pixel whose interpolated density is non-finite or at least
  `1e4 cm-3` is recomputed with the exact PyNeb `getTemDen` call.

The exact-fallback threshold is therefore:

```text
NE_SII >= 1e4 cm-3
```

The flag map is:

```text
NE_SII_EXACT_HIDEN
```

where 1 means that the exact high-density PyNeb fallback was used.

### 3.5 Low-density floor at 20 cm-3

Brazzini et al. 2024 use [S II] 6716/6731 for density, and for regions below
the low-density limit they impose:

```math
n_e = 20\ \mathrm{cm}^{-3}.
```

The local implementation follows this logic.  A pixel is flagged as fixed to
20 cm-3 if:

```math
R_{\mathrm{SII}} \ge R_{\mathrm{low-density}}
```

or if PyNeb returns a finite density below 20 cm-3.  The final density is:

```math
n_e =
\begin{cases}
20\ \mathrm{cm}^{-3}, & \mathrm{if\ below\ the\ low-density\ limit},\\
n_{e,\mathrm{PyNeb}}, & \mathrm{otherwise}.
\end{cases}
```

The flag map is:

```text
NE_SII_FIXED20
```

where 1 means the final `NE_SII` value was fixed to 20 cm-3.

The uncertainty is:

```text
NE_SII_ERR
```

computed by re-inverting the ratio at `R_SII - sigma(R_SII)` and
`R_SII + sigma(R_SII)`, then taking half the difference:

```math
\sigma(n_e) =
\frac{1}{2}
\left|
n_e(R+\sigma_R) - n_e(R-\sigma_R)
\right| .
```

For `NE_SII_FIXED20` pixels, the script sets:

```math
\sigma(n_e) = 0 .
```

This encodes the adopted floor, not a measurement that the physical density is
known exactly.

### 3.6 NE_SII_ALL

`NE_SII` is finite only when the [S II] doublet passes the detection gate.
For some downstream operations it is useful to have a density everywhere in the
spatially usable galaxy footprint.  The script therefore also writes:

```text
NE_SII_ALL
```

with:

```math
n_{e,\mathrm{all}} =
\begin{cases}
n_e([\mathrm{S\,II}]), & \mathrm{if}\ NE\_SII\ \mathrm{is\ finite},\\
20\ \mathrm{cm}^{-3}, & \mathrm{if\ the\ spatial\ pixel\ is\ usable},\\
\mathrm{NaN}, & \mathrm{if\ spatially\ masked}.
\end{cases}
```

In code, "spatially usable" is `np.isfinite(V_STARS2)` after the spatial mask
has been applied.

## 4. Common Te line ratios

All three metallicity methods use dust-corrected line intensities.

### 4.1 Te([N II])

The required [N II] temperature ratio is:

```math
R_{\mathrm{NII}} =
\frac{I(5755)}
{I(6548)+I(6583)} .
```

The script evaluates:

```python
PyNeb Atom("N", 2).getTemDen(
    R_NII,
    den=NE_SII,
    to_eval="L(5755)/(L(6548)+L(6584))",
)
```

and writes:

```text
TE_NII_HII, TE_NII_HII_ERR
TE_NII_SF,  TE_NII_SF_ERR
```

Only temperatures in the range 1000 to 30000 K are retained.

The error is a finite-difference propagation in the line ratio and density:

```math
\sigma(T_e)^2 =
\sigma(T_e; R)^2 + \sigma(T_e; n_e)^2 .
```

### 4.2 Te([S III])

The measured [S III] ratio is:

```math
R_{\mathrm{SIII}} =
\frac{I(6312)}
{I(9069)+I(9532)} .
```

Because the current products use 9069 but not a robust measured 9532 map, the
script adopts the theoretical relation:

```math
I(9532) = 2.5\,I(9069).
```

Thus:

```math
R_{\mathrm{SIII}} =
\frac{I(6312)}
{3.5\,I(9069)} .
```

The script evaluates:

```python
PyNeb Atom("S", 3).getTemDen(
    R_SIII,
    den=NE_SII,
    to_eval="L(6312)/(L(9069)+2.5*L(9069))",
)
```

and writes:

```text
TE_SIII_BR24_HII, TE_SIII_BR24_HII_ERR
TE_SIII_BR24_SF,  TE_SIII_BR24_SF_ERR
```

### 4.3 Inferred high-zone Te([O III])

MUSE does not normally cover the full [O III] auroral/nebular pair needed for a
direct Te([O III]) measurement in the same way as classical optical direct
methods.  Brazzini et al. 2024 therefore infer the high-zone temperature from
the [O III]-[S III] temperature relation, and state that the [O III]-[S III]
relation is preferred over [O III]-[N II] because it is tighter.

The local code uses the implemented Brazzini-chain relations:

```math
T_e([\mathrm{S\,III}]) =
1.22\,T_e([\mathrm{N\,II}]) - 2000\ \mathrm{K},
```

and:

```math
T_e([\mathrm{O\,III}]) =
0.80\,T_e([\mathrm{S\,III}]) + 2000\ \mathrm{K}.
```

For the full Brazzini-style branch, `TE_OIII_BR24` is inferred from measured
`TE_SIII_BR24`.  For the NII+OII7325 branch, `TE_OIII_NII_CHAIN` is inferred
from measured `TE_NII` through the two-step chain above.

The relation-error terms currently include the coefficient-error terms in the
code plus intrinsic scatters of:

```text
724 K   for Te([S III]) from Te([N II])
1270 K  for Te([O III]) from Te([S III])
```

## 5. Method 1: Brazzini+2024 multi-zone direct-style O/H

### 5.1 Output names

```text
O_PLUS_BR24_HII,       O_PLUS_BR24_HII_ERR
O_DPLUS_BR24_HII,      O_DPLUS_BR24_HII_ERR
O_H_BR24_DIRECT_HII,   O_H_BR24_DIRECT_HII_ERR

O_PLUS_BR24_SF,        O_PLUS_BR24_SF_ERR
O_DPLUS_BR24_SF,       O_DPLUS_BR24_SF_ERR
O_H_BR24_DIRECT_SF,    O_H_BR24_DIRECT_SF_ERR

BR24_DIRECT_VALID_HII
BR24_DIRECT_VALID_SF
```

### 5.2 Required detections

For a pixel to enter this branch, all of the following must hold:

```text
region mask is true: HII or SF
Hbeta detected
Halpha detected
SII6716 and SII6731 detected, giving finite NE_SII
NII5755 detected
NII6548+NII6583 detected
SIII6312 detected
SIII9069 detected
OII7318+7319+7329+7330 detected
OIII4959+5007 detected
```

Equivalently:

```math
M_{\mathrm{BR24}} =
M_{\mathrm{region}}
\land M_{\mathrm{Balmer}}
\land M_{n_e}
\land M_{T_{\mathrm{NII}}}
\land M_{T_{\mathrm{SIII}}}
\land M_{\mathrm{OII7325}}
\land M_{\mathrm{OIII}} .
```

This is the strictest branch.  It is normal for many galaxies to have zero valid
spaxels because it requires both faint auroral systems and red [O II].

### 5.3 Temperatures used

The branch uses the multi-zone assignment described by Brazzini et al. 2024:

```text
O+    uses measured Te([N II])
O++   uses Te([O III]) inferred from measured Te([S III])
```

The local implementation therefore sets:

```math
T_{\mathrm{low}} = T_e([\mathrm{N\,II}])
```

and:

```math
T_{\mathrm{high}} =
T_e([\mathrm{O\,III}])
= 0.80\,T_e([\mathrm{S\,III}]) + 2000\ \mathrm{K}.
```

### 5.4 Ionic abundances

The O+ intensity is:

```math
I([\mathrm{O\,II}]7325) =
I(7318)+I(7319)+I(7329)+I(7330).
```

The O++ intensity is:

```math
I([\mathrm{O\,III}]) =
I(4959)+I(5007).
```

PyNeb receives each line ratio with Hbeta normalized to 100:

```math
R_{\mathrm{OII}} =
100\frac{I([\mathrm{O\,II}]7325)}{I(\mathrm{H}\beta)}
```

and:

```math
R_{\mathrm{OIII}} =
100\frac{I([\mathrm{O\,III}]4959+5007)}{I(\mathrm{H}\beta)} .
```

Then:

```python
O+/H+ = PyNeb Atom("O", 2).getIonAbundance(
    int_ratio=R_OII,
    tem=Te_NII,
    den=NE_SII,
    to_eval="L(7318)+L(7319)+L(7329)+L(7330)",
    Hbeta=100,
)
```

and:

```python
O++/H+ = PyNeb Atom("O", 3).getIonAbundance(
    int_ratio=R_OIII,
    tem=Te_OIII_BR24,
    den=NE_SII,
    to_eval="L(4959)+L(5007)",
    Hbeta=100,
)
```

The total oxygen abundance is:

```math
\frac{O}{H} =
\frac{O^+}{H^+}
+
\frac{O^{++}}{H^+},
```

and the map value is:

```math
12+\log(O/H) =
12 + \log_{10}
\left[
\frac{O^+}{H^+}
+
\frac{O^{++}}{H^+}
\right].
```

This is the closest local analogue to the Brazzini et al. 2024 PHANGS-MUSE
direct-method abundance workflow.  The important difference is error handling:
Brazzini et al. use Monte Carlo line-flux realizations, while the current local
script uses deterministic finite-difference propagation.

## 6. Method 2: MAUVE NII+OII7325 semi-direct ionic-sum O/H

### 6.1 Output names

```text
O_PLUS_NII_OII7325_HII,       O_PLUS_NII_OII7325_HII_ERR
O_DPLUS_NII_OII7325_HII,      O_DPLUS_NII_OII7325_HII_ERR
O_H_NII_OII7325_HII,          O_H_NII_OII7325_HII_ERR

O_PLUS_NII_OII7325_SF,        O_PLUS_NII_OII7325_SF_ERR
O_DPLUS_NII_OII7325_SF,       O_DPLUS_NII_OII7325_SF_ERR
O_H_NII_OII7325_SF,           O_H_NII_OII7325_SF_ERR

NII_OII7325_VALID_HII
NII_OII7325_VALID_SF
```

### 6.2 Required detections

This branch requires:

```text
region mask is true: HII or SF
Hbeta detected
Halpha detected
SII6716 and SII6731 detected, giving finite NE_SII
NII5755 detected
NII6548+NII6583 detected
OII7318+7319+7329+7330 detected
OIII4959+5007 detected
```

It does not require [S III] 6312 or [S III] 9069.  Therefore it should usually
have more finite pixels than the full Brazzini-style branch.

### 6.3 Temperatures used

The low-zone temperature is measured:

```math
T_{\mathrm{low}} =
T_e([\mathrm{N\,II}]).
```

The high-zone temperature is inferred through the local Brazzini-chain
implementation:

```math
T_e([\mathrm{S\,III}])_{\mathrm{chain}} =
1.22\,T_e([\mathrm{N\,II}]) - 2000\ \mathrm{K},
```

then:

```math
T_e([\mathrm{O\,III}])_{\mathrm{chain}} =
0.80\,T_e([\mathrm{S\,III}])_{\mathrm{chain}}
+ 2000\ \mathrm{K}.
```

Combining those two equations:

```math
T_e([\mathrm{O\,III}])_{\mathrm{chain}} =
0.976\,T_e([\mathrm{N\,II}]) + 400\ \mathrm{K}.
```

The code keeps the two-step form so that the uncertainty includes both relation
error terms.

### 6.4 Ionic abundances and total abundance

The O+ abundance is computed from red [O II] 7325 with measured Te([N II]):

```math
\frac{O^+}{H^+}
=
f_{\mathrm{PyNeb}}
\left(
\frac{I(7318)+I(7319)+I(7329)+I(7330)}
{I(\mathrm{H}\beta)},
T_e([\mathrm{N\,II}]),
n_e
\right).
```

The O++ abundance is computed from [O III] 4959+5007 with the inferred
high-zone temperature:

```math
\frac{O^{++}}{H^+}
=
f_{\mathrm{PyNeb}}
\left(
\frac{I(4959)+I(5007)}
{I(\mathrm{H}\beta)},
T_e([\mathrm{O\,III}])_{\mathrm{chain}},
n_e
\right).
```

The total abundance is the same ionic sum:

```math
12+\log(O/H) =
12 + \log_{10}
\left[
\frac{O^+}{H^+}
+
\frac{O^{++}}{H^+}
\right].
```

This method is best described as a MAUVE-tailored semi-direct ionic-sum
experiment.  It is not a standard named PHANGS product.  Its motivation is that
Brazzini et al. 2024 already use Te([N II]) for O+ and red [O II] lines for the
MUSE-accessible O+ abundance, while older work such as Pilyugin & Thuan 2007
also discusses using [O II] 7320,7330 when the blue [O II] 3727 doublet is not
available.  The main caveat is that the high-zone temperature is inferred from
Te([N II]) rather than anchored by measured Te([S III]).

## 7. Method 3: Kreckel+2025 / Mendez-Delgado+2023 NII-only proxy

### 7.1 Output names

```text
O_H_NII_K25_HII,  O_H_NII_K25_HII_ERR
O_H_NII_K25_SF,   O_H_NII_K25_SF_ERR

NII_K25_VALID_HII
NII_K25_VALID_SF
```

### 7.2 Required detections

This branch requires only the common Te inputs plus [N II]:

```text
region mask is true: HII or SF
Hbeta detected
Halpha detected
SII6716 and SII6731 detected, giving finite NE_SII
NII5755 detected
NII6548+NII6583 detected
```

Then the script further requires the measured temperature to fall in the
calibration range used by the implementation:

```math
8000\ \mathrm{K}
\le
T_e([\mathrm{N\,II}])
\le
13000\ \mathrm{K}.
```

This extra temperature-range gate is why the integrated
`O_H_NII_K25` value can be `nan` even when many spaxels have finite
`TE_NII`: the integrated Te([N II]) summary must also be finite and within
8000-13000 K.

### 7.3 Metallicity equation

The implemented calibration is:

```math
12+\log(O/H) =
(-1.19)\,
\frac{T_e([\mathrm{N\,II}])}{10^4\ \mathrm{K}}
+ 9.68 .
```

The uncertainty includes the Te uncertainty and the quoted slope/intercept
calibration uncertainties:

```math
\sigma_Z =
\sqrt{
\left(
\frac{1.19\,\sigma(T_e)}{10^4}
\right)^2
+
\left[
0.14\,
\frac{T_e([\mathrm{N\,II}])}{10^4}
\right]^2
+
0.15^2
}.
```

This branch does not compute O+ and O++ separately.  It is a one-temperature
oxygen-abundance proxy.  Kreckel et al. 2025 motivate this because [N II] 5755
is the most commonly detected PHANGS-MUSE auroral line, and they show that
Te([N II])-based gradients agree well with other metallicity-gradient measures.
The calibration basis is tied to the Mendez-Delgado et al. 2023 temperature-
metallicity relation, so its absolute scale should not be interpreted as the
same thing as a classical collisionally excited-line ionic sum.

## 8. Error propagation summary

The local implementation is deterministic, not Monte Carlo.  For a ratio:

```math
R = \frac{A}{B},
```

it uses:

```math
\sigma(R) =
|R|
\sqrt{
\left(\frac{\sigma_A}{A}\right)^2
+
\left(\frac{\sigma_B}{B}\right)^2
}.
```

For PyNeb quantities such as density, temperature, and ionic abundance, the
script evaluates the quantity at the central value and at plus/minus one sigma
in each driving variable, then uses:

```math
\sigma_Q =
\frac{1}{2}|Q_+ - Q_-|.
```

Independent error contributions are added in quadrature.

For ionic-sum oxygen abundance:

```math
y =
\frac{O^+}{H^+}
+
\frac{O^{++}}{H^+},
```

and:

```math
Z = 12+\log_{10}(y),
```

the abundance uncertainty is:

```math
\sigma_Z =
\frac{
\sqrt{
\sigma(O^+/H^+)^2
+
\sigma(O^{++}/H^+)^2
}
}
{\ln(10)\,y}.
```

This is fast and transparent, but it does not include all covariance terms and
does not reproduce the full Monte Carlo procedure used by Brazzini et al. 2024.
The local errors should therefore be treated as first-order propagated errors.

## 9. Integrated values printed to the log

For each of the HII and SF masks, the script prints integrated Te summaries:

```text
NE_SII
Te(NII)
O_H_BR24_DIRECT
O_H_NII_OII7325
O_H_NII_K25
```

The integrated calculation is not the median of the spaxel maps.  It does:

1. Selects the pixels that pass the relevant valid mask.
2. Sums raw line fluxes over those pixels.
3. Computes one integrated Balmer decrement from the summed Halpha and Hbeta.
4. Computes one integrated E(B-V).
5. Dust-corrects the integrated line fluxes.
6. Repeats the same PyNeb density, temperature, ionic-abundance, or K25 proxy
   equations on the integrated line ratios.

This is why a map can have many finite pixels while an integrated value is
`nan`: the integrated branch still requires the summed-line ratios, integrated
density, integrated temperature, and method-specific range gates to be valid.
For `O_H_NII_K25`, the most common extra gate is:

```math
8000\ \mathrm{K}
\le
T_e([\mathrm{N\,II}])_{\mathrm{integrated}}
\le
13000\ \mathrm{K}.
```

## 10. Practical interpretation of the valid-pixel counters

The log block:

```text
PyNeb Te-method valid-pixel counts:
  NE_SII finite pixels: ...
  NE_SII fixed at 20 cm^-3: ...
  NE_SII exact high-density fallback pixels: ...
  NE_SII_ALL finite pixels: ...
  Te(NII) HII pixels: ...
  Te(NII) SF pixels: ...
  Brazzini+2024 direct HII pixels: ...
  Brazzini+2024 direct SF pixels: ...
  NII+OII7325 semi-direct HII pixels: ...
  NII+OII7325 semi-direct SF pixels: ...
  Kreckel+2025 NII-only HII pixels: ...
  Kreckel+2025 NII-only SF pixels: ...
```

should be read as follows:

- `NE_SII finite pixels`: pixels where the [S II] doublet passed detection and
  the final density is finite.
- `NE_SII fixed at 20 cm^-3`: subset of finite `NE_SII` pixels below the
  [S II] low-density limit or below the adopted 20 cm-3 floor.
- `NE_SII exact high-density fallback pixels`: pixels where interpolation was
  replaced by exact PyNeb `getTemDen`.
- `NE_SII_ALL finite pixels`: spatially usable pixels after the mask, using
  measured density when available and 20 cm-3 otherwise.
- `Te(NII)`: pixels where the common Te gate and [N II] temperature gate pass.
- `Brazzini+2024 direct`: pixels where [N II], [S III], [O II] 7325, and
  [O III] requirements all pass.
- `NII+OII7325 semi-direct`: pixels where [N II], [O II] 7325, and [O III]
  requirements pass, without requiring [S III].
- `Kreckel+2025 NII-only`: pixels where Te([N II]) is finite and in the
  8000-13000 K calibration range.

## 11. References

1. Brazzini, M. et al. 2024, "Metallicity calibrations based on auroral lines
   from PHANGS-MUSE data", A&A, 691, A173.
   DOI: https://doi.org/10.1051/0004-6361/202451007
   arXiv: https://arxiv.org/abs/2410.00106

2. Kreckel, K. et al. 2025, "Temperature based radial metallicity gradients in
   nearby galaxies", A&A, 703, A42.
   DOI: https://doi.org/10.1051/0004-6361/202556017
   arXiv: https://arxiv.org/abs/2507.20744

3. Mendez-Delgado, J. E. et al. 2023, "Temperature inhomogeneities cause the
   abundance discrepancy in H II regions", Nature, 618, 249-251.
   DOI: https://doi.org/10.1038/s41586-023-05956-2
   arXiv: https://arxiv.org/abs/2305.11578

4. Luridiana, V., Morisset, C., and Shaw, R. A. 2015, "PyNeb: a new tool for
   analyzing emission lines. I. Code description and validation of results",
   A&A, 573, A42.
   DOI: https://doi.org/10.1051/0004-6361/201323152
   ADS: https://ui.adsabs.harvard.edu/abs/2015A%26A...573A..42L/abstract

5. Morisset, C. et al. 2020, "Atomic Data Assessment with PyNeb", Atoms, 8, 66.
   DOI: https://doi.org/10.3390/atoms8040066
   ADS: https://ui.adsabs.harvard.edu/abs/2020Atoms...8...66M/abstract

6. Calzetti, D. et al. 2000, "The dust content and opacity of actively star-
   forming galaxies", ApJ, 533, 682.
   DOI: https://doi.org/10.1086/308692
   ADS: https://ui.adsabs.harvard.edu/abs/2000ApJ...533..682C/abstract

7. Cardelli, J. A., Clayton, G. C., and Mathis, J. S. 1989, "The relationship
   between infrared, optical, and ultraviolet extinction", ApJ, 345, 245.
   DOI: https://doi.org/10.1086/167900
   ADS: https://ui.adsabs.harvard.edu/abs/1989ApJ...345..245C/abstract

8. Pilyugin, L. S. and Thuan, T. X. 2007, "The Oxygen Abundance of Nearby
   Galaxies from Sloan Digital Sky Survey Spectra", ApJ, 669, 299.
   DOI: https://doi.org/10.1086/521188
   arXiv: https://arxiv.org/abs/0707.2856
