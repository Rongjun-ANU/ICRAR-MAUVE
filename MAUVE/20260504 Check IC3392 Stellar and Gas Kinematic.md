# 20260504 Check IC3392 Stellar and Gas Kinematic

## 1. IC3392 stellar and gas kinematic comparison across three runs

This notebook repeats the `KIN` map checks from `20260430_check_KIN.ipynb`, but only for IC3392 and with three runs. Here we compare both stellar `KIN` maps and the gas kinematic tracers `HA6562` and `NII6583`:

- `old_v721_7000`: original `v3_v7.2.1` run, old `CONFIG`, fitting to $7000\AA$.
- `v768_7000`: new `v3tk_v7.6.8_7000` run, v7.6.8 code and v3tk cube, old-style `CONFIG`, fitting to $7000\AA$.
- `v768_full`: original `v3tk_v7.6.8` run, v7.6.8 code and v3tk cube, wider/new `CONFIG`, fitting to 9100 A.

The key comparisons are:

- `version_proxy`: `v768_7000 - old_v721_7000`; this is the best proxy for the nGIST-version effect, with the v3 to v3tk cube change assumed negligible (less than 0.4%).
- ``CONFIG`_effect`: `v768_full - v768_7000`; this isolates the newer `CONFIG`/range/template/mask choices within v7.6.8 and the v3tk cube.
- `total_effect`: `v768_full - old_v721_7000`; this is the original old-to-new comparison.

Here we show the maps of 3 runs and 3 comparisons first in stellar velocity & dispersion, and then `HA6562` and `NII6583` as well:

![image-20260504184206232](assets/image-20260504184206232.png)

![image-20260504184932846](assets/image-20260504184932846.png)

![image-20260504185038258](assets/image-20260504185038258.png)

![image-20260504185046659](assets/image-20260504185046659.png)

![image-20260504185054465](assets/image-20260504185054465.png)

![image-20260504185059564](assets/image-20260504185059564.png)

And if we plot the comparisons in pairwise scatter, we confirm that the primary differences are more contributed by `CONFIG` than the `nGIST` version, especially in gas kinematic.  

![image-20260504190321661](assets/image-20260504190321661.png)

![image-20260504190346221](assets/image-20260504190346221.png)

![image-20260504190357519](assets/image-20260504190357519.png)

![image-20260504190404314](assets/image-20260504190404314.png)

![image-20260504190410803](assets/image-20260504190410803.png)

![image-20260504190416731](assets/image-20260504190416731.png)

## 2. IC3392 gas minus stellar LOS velocity residuals across three runs

Here we subtracts the stellar LOS velocity map from the gas LOS velocity maps for IC3392:

```
delta V = V_gas - V_stars
```

Again, by visual check, both runs with fitting upto $7000\AA$ look similar, while fitting to $9100\AA$ makes more changes.

![image-20260504190604102](assets/image-20260504190604102.png)

Then i also show the differences between residual maps across 3 comparisons. Similarly, `CONFIG` makes more differences. 

![image-20260504191944404](assets/image-20260504191944404.png)

## 3. IC3392 stellar and gas kinematic comparison across four runs

Now we extend the 3-run IC3392 kinematic checks with a fourth run:

- `old_v721_7000`: original `v3_v7.2.1`, old `CONFIG`, fitting to $7000\AA$.
- `v768_7000`: `v3tk_v7.6.8_7000`, v7.6.8 code with old-style `CONFIG`, fitting to $7000\AA$.
- `v768_full`: original `v3tk_v7.6.8`, wide/new `CONFIG` with the full red-line GAS file.
- `v768_uptoSII`: same broad setup as `v768_full`, but `GAS` module uses the shorter emission-line file up to [S II]$\lambda\lambda6716,30$, which is same as fitting to $7000\AA$.

Here the primary test is `extra_red_lines_effect = v768_full - v768_uptoSII`. This isolates the effect of including the red/faint emission lines (which tie the amplitudes and/or kinematic to strong lines) in the `GAS` fit while leaving the stellar kinematic setup unchanged.

The stellar kinematic will be exactly the same in these two runs, and the gas kinematic are still very similar visually, at least in high S/N regions.

![image-20260504204536263](assets/image-20260504204536263.png)

Then, we also show the pairwise scatter and residual plots for `v768_full` vs `v768_uptoSII`,  `v768_full` vs `old_v721_7000`, and `v768_uptoSII` vs `old_v721_7000`. Clearly, we can say that fitting red/faint-line does not perturb the results with new $9100\AA$ `CONFIG`. 

![image-20260504205620376](assets/image-20260504205620376.png)

![image-20260504205631327](assets/image-20260504205631327.png)

![image-20260504205640780](assets/image-20260504205640780.png)