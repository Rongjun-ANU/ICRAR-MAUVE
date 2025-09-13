# 20250515 Dust and [Lee et al. 2025](https://iopscience.iop.org/article/10.3847/1538-3881/adb285#ajadb285s2)

## Dust

As mentioned yesterday, the `dust` parameters actually converge to $A_V = 0.897$ and $\delta = -1.089$. Now I explore how dust attenuation curves vary with $\delta$:

![image-20250515195820411](/Users/maclaptop29/Desktop/ICRAR/MAUVE/20250515 Dust and Lee et al 2025.assets/image-20250515195820411.png)

Here, all gray curves are dust attenuation curves with $\delta$ from -1.2 to 0.2 (steeper $\rightarrow$ smaller $\delta$), the black curve is $\delta = -1.089$, the blue curve is $\delta = 0$, and cyan one is from [Calzetti et al. 2000](https://iopscience.iop.org/article/10.1086/308692/fulltext/) ($\delta = 0$ but no bump at $2175\AA$). Clearly, since we are only interested in optical band so there is no significant difference for choosing difference . Hence, for simplicity, from now on I will fix $\delta = 0$ (main-sequence galaxy). 

## Compare stellar mass with [Lee et al. 2025](https://iopscience.iop.org/article/10.3847/1538-3881/adb285#ajadb285s2)

In Figure 13 of [Lee et al. 2025](https://iopscience.iop.org/article/10.3847/1538-3881/adb285#ajadb285s2), they compare stellar masses obtained from `pPXF` with best-fit models with linear combinations of specific luminosity (at $1.63 \mu m$, one of the channels in simulated SPHEREx data) and SDSS *g* − *r* color for different SPS templates:

![img](https://content.cld.iop.org/journals/1538-3881/169/3/185/revision1/ajadb285f13_hr.jpg)

However, in legacy data, I only have h-band ($1.66\mu m$) from 2MASS, so I take it as the $\lambda$ here. By running my previous pipeline, I get $\log L_h/L_\odot = 9.43$ (h-band magnitude is 11.427, which is the same as the data in [NED](https://ned.ipac.caltech.edu/byname?objname=ic3392&hconst=67.8&omegam=0.308&omegav=0.692&wmap=4&corr_z=1)). Now with $g-r = 12.658-11.958 = 0.7$, I can estimate the stellar mass by model 2 and compare with those obtained by `pPXF`: 

|                | model 2         | `pPXF` |
| :------------- | :-------------- | :----- |
| E-MILES        | 9.269$\pm$0.132 | 9.355  |
| BC03 (GALAXEV) | 9.015$\pm$0.126 | 9.314  |
| FSPS           | 9.267$\pm$0.121 | 9.385  |

As mentioned in [Lee et al. 2025](https://iopscience.iop.org/article/10.3847/1538-3881/adb285#ajadb285s2) that "BC03 and CB19 yield mass-to-light ratios on average  $\sim$ 0.2−0.3 dex lower than those from E-MILES and FSPS",  we do observe that BC03 is $\sim0.25$ dex lower than E-MILES and FSPS with my legacy data using model 2. On the other hand, we can see that results using E-MILES and FSPS obtained from `pPXF` are still within uncertainty ranges of  those from model 2. Again, BC03 systematically underestimates $M_*/L$ than E-MILES and FSPS, especially using $M_*/L$ vs color relation. 