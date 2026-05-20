from matplotlib import pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
plt.rcParams['figure.figsize'] = [30, 30]
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams.update({'font.size': 8})
from matplotlib.backends.backend_pdf import PdfPages


##################################################
data_folder='/arc/projects/mauve/products/v0.6'
galaxyid = 'NGC4192'
gas_bin_sub = '_gas_BIN_maps'
stars_sub ='_KIN_maps'
sfh_sub ='_SFH_maps'
ls_sub = '_LS_ADAPTED_maps'
version ='v0.6'
#################################################

gas_bin_image = pyfits.open(data_folder+'/'+galaxyid+'/'+galaxyid+gas_bin_sub+'.fits')
stars_bin_image = pyfits.open(data_folder+'/'+galaxyid+'/'+galaxyid+stars_sub+'.fits')
sfh_bin_image = pyfits.open(data_folder+'/'+galaxyid+'/'+galaxyid+sfh_sub+'.fits')
ls_bin_image = pyfits.open(data_folder+'/'+galaxyid+'/'+galaxyid+ls_sub+'.fits')

out=galaxyid+'_'+version+'_QC.pdf'

id= ['HB4861_FLUX',
     'HB4861_FLUX_ERR',
     'HB4861_VEL',
     'HB4861_SIGMA',
     'HB4861_SIGMA_ERR',
     'OIII5006_FLUX',
     'OIII5006_FLUX_ERR',
     'OIII5006_VEL',
     'OIII5006_SIGMA',
     'OIII5006_SIGMA_ERR',
     'OI6300_FLUX',
     'OI6300_FLUX_ERR',
     'OI6300_VEL',
     'OI6300_SIGMA',
     'OI6300_SIGMA_ERR',
     'OI6363_FLUX',
     'OI6363_FLUX_ERR',
     'OI6363_VEL',
     'OI6363_SIGMA',
     'OI6363_SIGMA_ERR',
     'HA6562_FLUX',
     'HA6562_FLUX_ERR',
     'HA6562_VEL',
     'HA6562_SIGMA',
     'HA6562_SIGMA_ERR',
     'NII6583_FLUX',
     'NII6583_FLUX_ERR',
     'NII6583_VEL',
     'NII6583_SIGMA',
     'NII6583_SIGMA_ERR',
     'SII6716_FLUX',
     'SII6716_FLUX_ERR',
     'SII6716_VEL',
     'SII6716_SIGMA',
     'SII6716_SIGMA_ERR',
     'Ha/Hb',
     'NII/Ha',
     'SII/NII',
     'Ha-Hb vel',
     'Ha/Hb sigma',
     'V',
     'SIGMA',
     'H3',
     'H4',
     'FORM_ERR_SIGMA',
     'Vs-Vha',
     'Sigma_s - Sigma_ha',
     'SII6716/SII6730',
     'VHa-VNII',
     'SHa-SNII',
     'AGE',
     'METAL',
     'AGE',
     'METAL',
     'ALPHA',
     'NAD',
     'EBV',
    ]



with PdfPages(out) as pdf:
    for i in range(0, len(id)):
        # Create a new figure every 6 plots (2x3 grid)
        if i % 5 == 0:
            fig, axs = plt.subplots(2, 3, figsize=(20, 12))
            axs = axs.flatten()  # Flatten to make it easier to index each subplot
            
        # Select the subplot based on the position within the 2x3 grid
        ax = axs[i % 5]

        # Clear any unused subplot for a clean layout
        if (i % 5) == 4 or i == len(id) - 1:
            for j in range((i % 5) + 1, 6):
                fig.delaxes(axs[j])

        count = (i % 5) + 1
        #######
        #######
        if (i<35):
            # Retrieve data for the current plot
            to_plot = gas_bin_image[id[i]].data
            mean = np.nanmean(to_plot) 
            std  = np.nanstd(to_plot)  

            # Plot based on the count (1 to 5) for the different visualization settings
            if count == 1:
                if (~np.isnan(mean)):
                    img=ax.imshow(to_plot, origin='lower', norm=LogNorm(vmin=2), cmap='plasma')
            elif count == 2:
                if (~np.isnan(mean)):
                    img=ax.imshow(to_plot, origin='lower', norm=LogNorm())
            elif count == 3:
                if (~np.isnan(mean)):
                    img=ax.imshow(to_plot, origin='lower', vmin=mean - 1 * std, vmax=mean + 1 * std, cmap='PiYG_r')
            elif count == 4:
                if (~np.isnan(mean)):
                    img=ax.imshow(to_plot, origin='lower', vmin=50, vmax=100, cmap='magma')
            elif count == 5:
                if (~np.isnan(mean)):
                    img=ax.imshow(to_plot, origin='lower', vmin=1, vmax=30, cmap='cividis')
                    
                
                
        ######
        ######
        if ((i>=35) & (i<40)):
            if count==1:
                to_plot= gas_bin_image['HA6562_FLUX'].data / gas_bin_image['HB4861_FLUX'].data
                img=ax.imshow(to_plot,origin='lower',norm=LogNorm(vmin=2.5,vmax=8),cmap='plasma')

            elif count==2:
                to_plot= gas_bin_image['NII6583_FLUX'].data / gas_bin_image['HA6562_FLUX'].data
                img=ax.imshow(to_plot,origin='lower',norm=LogNorm(vmin=0.1,vmax=1))
        
            elif count==3:
                to_plot= gas_bin_image['SII6716_FLUX'].data / gas_bin_image['NII6583_FLUX'].data
                mean=np.nanmedian(to_plot)
                vmin=mean/2.
                vmax=mean*2.
                img=ax.imshow(to_plot,origin='lower',norm=LogNorm(vmin=vmin,vmax=vmax),cmap='magma')
            elif count==4:
                to_plot= gas_bin_image['HA6562_VEL'].data - gas_bin_image['HB4861_VEL'].data
                img=ax.imshow(to_plot,origin='lower',vmin=-3,vmax=+3,cmap='PiYG_r')
            elif count==5:
                to_plot= gas_bin_image['HA6562_SIGMA'].data / gas_bin_image['HB4861_SIGMA'].data
                img=ax.imshow(to_plot,origin='lower',vmin=0.6,vmax=1.2,cmap='PiYG')
    
       #######
       #######
        if ((i>=40) & (i<45)):
            to_plot= stars_bin_image[id[i]].data
            mean=np.nanmean(to_plot)
            std=np.nanstd(to_plot)
            if count==1:
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='RdBu_r')
            elif count==2:
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='magma')
            elif count==3:
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='PRGn_r')
            elif count==4:
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='RdYlBu')
            elif count==5:
                img=ax.imshow(to_plot,origin='lower',vmin=10,vmax=50,cmap='cividis')
            
        if ((i>=45) & (i<50)):
            if count==1:
                to_plot= gas_bin_image['HA6562_VEL'].data - stars_bin_image['V'].data
                img=ax.imshow(to_plot,origin='lower',vmin=-20,vmax=+20,cmap='PuOr_r')

            elif count==2:
                to_plot= gas_bin_image['HA6562_SIGMA'].data - stars_bin_image['SIGMA'].data
                img=ax.imshow(to_plot,origin='lower',vmin=-30, vmax=30, cmap='PRGn')
        
            elif count==3:
                to_plot= gas_bin_image['SII6716_FLUX'].data / gas_bin_image['SII6730_FLUX'].data
                img=ax.imshow(to_plot,origin='lower',vmin=0.6, vmax=1.7)

            elif count==4:
                to_plot= gas_bin_image['HA6562_VEL'].data - gas_bin_image['NII6583_VEL'].data
                img=ax.imshow(to_plot,origin='lower',vmin=-20,vmax=+20,cmap='BrBG_r')
        
            elif count==5:
                to_plot= gas_bin_image['HA6562_SIGMA'].data - gas_bin_image['NII6583_SIGMA'].data
                img=ax.imshow(to_plot,origin='lower',vmin=-20,vmax=+20,cmap='BrBG_r')
        
       #######
       #######
        if ((i>=50) & (i<55)):
            if count==1:
                to_plot= sfh_bin_image[id[i]].data
                mean=np.nanmean(to_plot)
                std=np.nanstd(to_plot)
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='magma')
            elif count==2:
                to_plot= sfh_bin_image[id[i]].data
                mean=np.nanmean(to_plot)
                std=np.nanstd(to_plot)
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='plasma')
            elif count==3:
                to_plot= ls_bin_image[id[i]].data
                mean=np.nanmean(to_plot)
                std=np.nanstd(to_plot)
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='cividis')
            elif count==4:
                to_plot= ls_bin_image[id[i]].data
                mean=np.nanmean(to_plot)
                std=np.nanstd(to_plot)
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='magma')
            elif count==5:
                to_plot= ls_bin_image[id[i]].data
                mean=np.nanmean(to_plot)
                std=np.nanstd(to_plot)
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='plasma')
             
        
       #######
       #######
        if ((i>=55) & (i<60)):
            if count==1:
                to_plot= ls_bin_image[id[i]].data
                mean=np.nanmedian(to_plot)
                std=np.nanstd(to_plot)
                img=ax.imshow(to_plot,origin='lower',vmin=0,vmax=mean+0.5*std,cmap='magma')
        
            elif count==2:
                to_plot= sfh_bin_image[id[i]].data
                mean=np.nanmedian(to_plot)
                std=np.nanstd(to_plot)
                img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='magma')
       
        
        
        
        
        
        
        
        
        
        
        
        
        # Set title or other subplot details as needed
        ax.set_title(f"{id[i]}")
        cbar = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.04)
        # After every 6 plots, save the page and close the figure
        if (i + 1) % 5 == 0 or i == len(id) - 1:
            plt.tight_layout()
            pdf.savefig(fig)  # Save the current page in the PDF
            plt.close(fig)  # Close the figure to free memory



