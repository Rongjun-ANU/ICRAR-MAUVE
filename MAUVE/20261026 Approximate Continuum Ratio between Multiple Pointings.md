# 20261026 Approximate Continuum Ratio between Multiple Pointings

Although i dont have the pointing information, we use the `VARIANCE` maps at different wavelengths to roughly segment each galaxy into miltiple pointings. This works well if the pointing pattern is clear but fails when pointings are messy like NGC4264, NGC4321 and NGC4383. Anyway, then i take the continuum flux ratios of these pairs of neighouring pointings to inspect if we still have some large scale sine-like variations along wavelength, so that we may increase `MDEG` in `ppxf` fitting to mitigate the pointing issue in our final mass-to-light ratio/stellar mass map.

So here i show how i determine the pointing and overlapping regions of each galaxy, and then show the continuum ratio vs wavelength of each pair of neighbouring pointings (both 7000 and 8900 runs), and finally the SFH QC maps (7000 run only).

Below is what i think there are clear sine-like residual between two neighbouring pointings and so we may need to redo the `SFH` module in 7000 run with `MDEG` increase from 12 to 16:

1. NGC4216: R6/R7, the top two pointings. 
2. NGC4351
3. NGC4380: R3/R6, R5/R6
4. NGC4450: R1/R4, R2/R5, R4/R5
5. NGC4457: R2/R4, R3/R4
6. NGC4501: R4/R8, R7/R8
7. NGC4535: R1/R4, R4/R5
8. NGC4654: R3/R4, R5/R6, R7/R8
9. NGC4689: R2/R3, R3/R6, R5/R6 show clear sine-like variations, and they are indeed fixed by refitting with `MDEG=16`. 

Also some galaxies I am unsure: 

1. NGC4192, 
2. NGC4254, 
3. NGC4298 (SFH QC maps seem fine, but ratio vs wavelength show sine-like variation in R1/R2, R3/R5 and R5/R6), 
4. NGC4383, 
5. NGC4567_8: the top left R3, maybe worth to try
6. NGC4569:  bottom right R6, maybe worth to try
7. NGC4698

And for the resf of galaxies, I am happy with their current `nGIST` products: 

1. NGC4064, 
2. NGC4189,  
3. NGC4222, 
4. NGC4293, 
5. NGC4294, 
6. NGC4302, 
7. NGC4321, 
8. NGC4330, 
9. NGC4388, 
10. NGC4394, 
11. NGC4396 (although R2/R3 seems a bit sine-like vacation, the QC maps look fine for me), 
12. NGC4402, 
13. NGC4405, 
14. NGC4419 (although seems pointing residual in QC maps, the continuum ratio does not show sine-like variation so increasing `MDEG` does not expect to fix this issue), 
15. NGC4424, 
16. NGC4522,
17. NGC4580, 
18. NGC4607.

## NGC4064 (fine)

![image-20261026142337141](assets/image-20261026142337141.png)

![image-20261026142343189](assets/image-20261026142343189.png)

![image-20261026142417952](assets/image-20261026142417952.png)

## NGC4189 (fine)

![image-20261026142433525](assets/image-20261026142433525.png)

![image-20261026142509870](assets/image-20261026142509870.png)

![image-20261026142524555](assets/image-20261026142524555.png)

## NGC4192 (maybe)

![image-20261026142619182](assets/image-20261026142619182.png)

![image-20261026142624597](assets/image-20261026142624597.png)

![image-20261026144427125](assets/image-20261026144427125.png)

## NGC4216 (needed)

![image-20261026142918230](assets/image-20261026142918230.png)

![image-20261026142926880](assets/image-20261026142926880.png)

![image-20261026144523917](assets/image-20261026144523917.png)

## NGC4222 (fine)

![image-20261026143520838](assets/image-20261026143520838.png)

![image-20261026143528429](assets/image-20261026143528429.png)

![image-20261026143742782](assets/image-20261026143742782.png)

## NGC4254 (maybe)

![image-20261026143836957](assets/image-20261026143836957.png)

![image-20261026143846504](assets/image-20261026143846504.png)

![image-20261026143939260](assets/image-20261026143939260.png)

## NGC4293 (fine)

![image-20261026144020391](assets/image-20261026144020391.png)

![image-20261026144031369](assets/image-20261026144031369.png)

![image-20261026144049786](assets/image-20261026144049786.png)

## NGC4294 (fine)

![image-20261026144642144](assets/image-20261026144642144.png)

![image-20261026144650674](assets/image-20261026144650674.png)

![image-20261026144709866](assets/image-20261026144709866.png)

## NGC4298 (maybe)

![image-20261026144852151](assets/image-20261026144852151.png)

![image-20261026144927805](assets/image-20261026144927805.png)

![image-20261026144946911](assets/image-20261026144946911.png)

## NGC4302 (fine)

![image-20261026145400875](assets/image-20261026145400875.png)

![image-20261026145407190](assets/image-20261026145407190.png)

![image-20261026145417074](assets/image-20261026145417074.png)

## NGC4321 (fine)

![image-20261026145505541](assets/image-20261026145505541.png)

![image-20261026145511200](assets/image-20261026145511200.png)

![image-20261026145844730](assets/image-20261026145844730.png)

## NGC4330 (fine)

![image-20261026150021972](assets/image-20261026150021972.png)

![image-20261026150028545](assets/image-20261026150028545.png)

![image-20261026150103343](assets/image-20261026150103343.png)

## NGC4351 (needed)

![image-20261026152837971](assets/image-20261026152837971.png)

![image-20261026152848472](assets/image-20261026152848472.png)

![image-20261026154547600](assets/image-20261026154547600.png)

## NGC4380 (needed)

![image-20261026154644561](assets/image-20261026154644561.png)

![image-20261026154649507](assets/image-20261026154649507.png)

![image-20261026154746871](assets/image-20261026154746871.png)

## NGC4383 (maybe)

![image-20261026154833102](assets/image-20261026154833102.png)

![image-20261026154915994](assets/image-20261026154915994.png)

![image-20261026155019640](assets/image-20261026155019640.png)

## NGC4388 (fine)

![image-20261026155113846](assets/image-20261026155113846.png)

![image-20261026155117386](assets/image-20261026155117386.png)

![image-20261026155157445](assets/image-20261026155157445.png)

## NGC4394 (fine)

![image-20261026155621942](assets/image-20261026155621942.png)

![image-20261026155626275](assets/image-20261026155626275.png)

![image-20261026155643006](assets/image-20261026155643006.png)

## NGC4396 (fine)

![image-20261026160547733](assets/image-20261026160547733.png)

![image-20261026160553526](assets/image-20261026160553526.png)

![image-20261026160600658](assets/image-20261026160600658.png)

## NGC4402 (fine)

![image-20261026160937647](assets/image-20261026160937647.png)

![image-20261026160945642](assets/image-20261026160945642.png)

![image-20261026161125834](assets/image-20261026161125834.png)

## NGC4405 (fine)

![image-20261026161200478](assets/image-20261026161200478.png)

![image-20261026161207132](assets/image-20261026161207132.png)

![image-20261026161543687](assets/image-20261026161543687.png)

## NGC4419 (fine)

![image-20261026161635914](assets/image-20261026161635914.png)

![image-20261026161641349](assets/image-20261026161641349.png)

![image-20261026161716331](assets/image-20261026161716331.png)

## NGC4424 (fine)

![image-20261026162013840](assets/image-20261026162013840.png)

![image-20261026162019873](assets/image-20261026162019873.png)

![image-20261026162139267](assets/image-20261026162139267.png)

## NGC4450 (needed)

![image-20261026162247835](assets/image-20261026162247835.png)

![image-20261026162427583](assets/image-20261026162427583.png)

![image-20261026162435223](assets/image-20261026162435223.png)

## NGC4457 (needed)

![image-20261026162535517](assets/image-20261026162535517.png)

![image-20261026162539073](assets/image-20261026162539073.png)

![image-20261026162639770](assets/image-20261026162639770.png)

## NGC4501 (needed)

![image-20261026162749539](assets/image-20261026162749539.png)

![image-20261026162758108](assets/image-20261026162758108.png)

![image-20261026162819007](assets/image-20261026162819007.png)

## NGC4522 (fine)

![image-20261026163026542](assets/image-20261026163026542.png)

![image-20261026163033006](assets/image-20261026163033006.png)

![image-20261026163051980](assets/image-20261026163051980.png)

## NGC4535 (needed)

![image-20261026163323380](assets/image-20261026163323380.png)

![image-20261026163328113](assets/image-20261026163328113.png)

![image-20261026163336003](assets/image-20261026163336003.png)

## NGC4567_8 (maybe)

![image-20261026163658386](assets/image-20261026163658386.png)

![image-20261026163707107](assets/image-20261026163707107.png)

![image-20261026164013038](assets/image-20261026164013038.png)

## NGC4569 (maybe)

![image-20261026164311816](assets/image-20261026164311816.png)

![image-20261026164318017](assets/image-20261026164318017.png)

![image-20261026164607352](assets/image-20261026164607352.png)

## NGC4580 (fine)

![image-20261026164941583](assets/image-20261026164941583.png)

![image-20261026164946592](assets/image-20261026164946592.png)

![image-20261026165015005](assets/image-20261026165015005.png)

## NGC4607 (fine)

![image-20261026165116101](assets/image-20261026165116101.png)

![image-20261026165151641](assets/image-20261026165151641.png)

![image-20261026165236023](assets/image-20261026165236023.png)

## NGC4654 (needed)

![image-20261026165323640](assets/image-20261026165323640.png)

![image-20261026165337465](assets/image-20261026165337465.png)

![image-20261026165457348](assets/image-20261026165457348.png)

## NGC4689

![image-20261026165601602](assets/image-20261026165601602.png)

![image-20261026165617220](assets/image-20261026165617220.png)

![image-20261026165636057](assets/image-20261026165636057.png)

![image-20261026170301556](assets/image-20261026170301556.png)

## NGC4698 (maybe)

![image-20261026170228082](assets/image-20261026170228082.png)

![image-20261026170234453](assets/image-20261026170234453.png)

![image-20261026170404655](assets/image-20261026170404655.png)