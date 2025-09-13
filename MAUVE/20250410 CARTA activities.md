# CARTA activities

## IC3392

https://cloud.datacentral.org.au/remote.php/webdav/mauve/teamdata/cubes/red_v3/IC3392_DATACUBE_FINAL_WCS_Pall_mad_red_v3.fits

### 1. Load data in CARTA and check extensions

I see three extensions: DATA, STAT, and DQ. These match the description on MAUVE wiki.  

![image-20250410154406662](/Users/maclaptop29/Library/Application Support/typora-user-images/image-20250410154406662.png)

### 2. Skim through the Animator tool

Just a bit confused with the terminology, `frame rate`. But I guess the animator go through the z-direction of datacube  so i think that should refer to wavelength space?

![image-20250410154823191](/Users/maclaptop29/Library/Application Support/typora-user-images/image-20250410154823191.png)

### 3. Spectral Profiler

I stop the animator at frame 1490 (i.e., 6562.5 $\AA$). I open the spectral profiler and choose a datapoint roughly the galactic centre. 

![image-20250410230122404](/Users/maclaptop29/Desktop/MAUVE/20250410 CARTA activities.assets/image-20250410230122404.png)

I do see some emission lines. I can try to identify 5 of them: 

N II doublet, H$\alpha$ and S II doublet?

![image-20250410230933819](/Users/maclaptop29/Desktop/MAUVE/20250410 CARTA activities.assets/image-20250410230933819.png)

Also a question here, is there any figures or tables that I can use as a reference to identify them. 

Then I move the data point to lefter region and clearly the intensity drops. 

![image-20250410231349713](/Users/maclaptop29/Desktop/MAUVE/20250410 CARTA activities.assets/image-20250410231349713.png)

### 4. Moments

I guess the moment tool is integrating the image over the selected wavelength range. Here i choose stacking 6500 - 6800 $\AA$ and subtracting the continuum (intensity below 400) to get a H$\alpha$ map, if that peak at 6600$\AA$ is the redshifted H$\alpha$ line (z = 0.0056). 

![image-20250411000121243](/Users/maclaptop29/Desktop/MAUVE/20250410 CARTA activities.assets/image-20250411000121243.png)

