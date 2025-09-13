# 20250701 NII and SII BPT Comparison.md

## 1. 4 catalogs

I have the QC for each line: $\rm FLUX/FLUX_{ERR} \geq 3$ and $\rm FLUX \geq 20\times10\ erg/s/cm^2$. 

Now I have divided all non-nan spaxels into 4 catalogs (here I use [NII] BPT as an example):

1. `Unknown` (in color blue): not detected in Balmer lines, i.e., either H$\alpha$ or H$\beta$ or both does not pass the QC. 
2. `Unconstrainted` (in color red): detected in Balmer lines, but error bars are too large to locate the spaxel's position on the [NII] BPT diagram. 
3. `Upper` (in color orange): detected in Balmer lines, constrainted on [NII] BPT diagram, and above the Kewley+2001 curve (red curve). 
4. `SF` (in color green): detected in Balmer lines, constrainted on [NII] BPT diagram, and below the Kewley+2001 curve (red curve). 

## 2. Comparison

Now I apply this classification for both [NII] and [SII] BPT diagrams to see any difference. As expected, most of galaxies have similar distributions across two BPT diagrams, but some of them are different. Here are three cases:

1. Similar distributions: IC3392, NGC4064, NGC4192, NGC4298, NGC4330, NGC4457, NGC4624
2. Differences in `Upper` (`Upper` regions appear in [NII] but become `Unconstrainted` in [SII] BPT diagram): NGC4293 (the biggest differences), NGC4419, NGC4501, NGC4698
3. Differences in `SF` (extended `SF` regions in [NII] compared to [SII] BPT diagram): NGC4383 (the biggest differences), NGC4396, NGC4522

Meanwhile, I am still trying to understand the physics behind the ionization of these two lines, but here I show the details of each galaxy. 

### 2.1 IC3392

![image-20250701170726302](assets/image-20250701170726302.png)

### 2.2 NGC4064

![image-20250701170741292](assets/image-20250701170741292.png)

### 2.3 NGC4192

![image-20250701171615604](assets/image-20250701171615604.png)

### 2.4 NGC4293

![image-20250701171717418](assets/image-20250701171717418.png)

### 2.5 NGC4298

![image-20250701171824755](assets/image-20250701171824755.png)

### 2.6 NGC4330

![image-20250701171937595](assets/image-20250701171937595.png)

### 2.7 NGC4383

![image-20250701172009035](assets/image-20250701172009035.png)

### 2.8 NGC4396

![image-20250701172057538](assets/image-20250701172057538.png)

### 2.9 NGC4419

![image-20250701172200263](assets/image-20250701172200263.png)

### 2.10 NGC4457

![image-20250701172233243](assets/image-20250701172233243.png)

### 2.11 NGC4501

![image-20250701172401570](assets/image-20250701172401570.png)

### 2.12 NGC4522

![image-20250701172444784](assets/image-20250701172444784.png)

### 2.13 NGC4694

![image-20250701172514634](assets/image-20250701172514634.png)

### 2.14 NGC4698

![image-20250701172525775](assets/image-20250701172525775.png)

## 3. Together

To put these two BPT selections together, the current approach I use is to make sure each spaxel is classified on `both` BPT diagrams. That means:

1. `Unknown`: still not detected in Balmer lines
2. `SF` (in color green): `SF` on both [NII] and [SII] BPT diagram simulanously. 
3. `Upper` (in color orange): detected in Balmer lines, constrainted on both [NII] and [SII] BPT diagram simulanously, but not `SF` on either  [NII] or [SII] BPT diagram or both. 
4. `Unconstrainted` (in color red): detected in Balmer lines, but `Unconstrainted` on either  [NII] or [SII] BPT diagram or both. 

Again, below are the details of applying the `both` selection on each galaxy.

### 3.1 IC3392

![image-20250701170224223](assets/image-20250701170224223.png)

### 3.2 NGC4064

![image-20250701170250190](assets/image-20250701170250190.png)

### 3.3 NGC4192

![image-20250701171644766](assets/image-20250701171644766.png)

### 3.4 NGC4293

![image-20250701171730256](assets/image-20250701171730256.png)

### 3.5 NGC4298

![image-20250701171839898](assets/image-20250701171839898.png)

### 3.6 NGC4330

![image-20250701171951467](assets/image-20250701171951467.png)

### 3.7 NGC4383

![image-20250701172029595](assets/image-20250701172029595.png)

### 3.8 NGC4396

![image-20250701172134425](assets/image-20250701172134425.png)

### 3.9 NGC4419

![image-20250701172217362](assets/image-20250701172217362.png)

### 3.10 NGC4457

![image-20250701172337815](assets/image-20250701172337815.png)

### 3.11 NGC4501

![image-20250701172427781](assets/image-20250701172427781.png)

### 3.12 NGC4522

![image-20250701172456647](assets/image-20250701172456647.png)

### 3.13 NGC4694

![image-20250701172540204](assets/image-20250701172540204.png)

### 3.14 NGC4698

![image-20250701172559328](assets/image-20250701172559328.png)
