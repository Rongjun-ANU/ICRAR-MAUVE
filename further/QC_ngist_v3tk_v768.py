#!/usr/bin/env python
import sys
from matplotlib import pyplot as plt
import numpy as np
import astropy.io.fits as pyfits
plt.rcParams['figure.figsize'] = [30, 30]
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams.update({'font.size': 8})
from matplotlib.backends.backend_pdf import PdfPages


def calzetti_k(w_um):
    """Return k(lambda) = A(lambda) / E(B-V) for Calzetti (2000)."""
    w = np.asarray(w_um, dtype=float)
    Rv = 4.05
    k = np.empty_like(w, dtype=float)

    short = (w >= 0.12) & (w < 0.63)
    long = (w >= 0.63) & (w <= 2.2)

    k[short] = 2.659 * (
        -2.156 + 1.509 / w[short] - 0.198 / w[short]**2 + 0.011 / w[short]**3
    ) + Rv
    k[long] = 2.659 * (-1.857 + 1.040 / w[long]) + Rv
    k[~(short | long)] = np.nan
    return k.item() if k.ndim == 1 and k.size == 1 else k


def mass_correction_map(sfh_hdul, mass_hdul, snr_threshold=25.0, lam_eff_um=0.623):
    """Map the stellar-mass correction term A(r)/2.5 + log10(M/L_R)."""
    ebv = np.asarray(sfh_hdul["EBV"].data, dtype=float)
    snr_postfit = np.asarray(sfh_hdul["SNR_POSTFIT"].data, dtype=float)
    ml_r = np.asarray(mass_hdul["ML_R"].data, dtype=float)

    ml_r_positive = np.where(ml_r > 0, ml_r, np.nan)
    correction = calzetti_k(lam_eff_um) * ebv / 2.5 + np.log10(ml_r_positive)
    valid = (snr_postfit >= snr_threshold) & np.isfinite(correction) & (ml_r > 0)
    return np.where(valid, correction, np.nan)


def open_mass_product(data_folder, galaxyid):
    path = f"{data_folder}/{galaxyid}/{galaxyid}_spatial_binning_maps_further.fits"
    try:
        return pyfits.open(path)
    except FileNotFoundError:
        print(f"Warning: mass product not found, skipping mass-correction QC map: {path}")
        return None


##################################################
data_folder='/Users/Igniz/Desktop/ICRAR/further/v3tk_v7.6.8'
# On Pawsey these products are stored in:
# /scratch/pawsey1308/mauve/products/v3tk_v7.6.8
galaxyid = 'IC3392'
if len(sys.argv) > 1:
    galaxyid = sys.argv[1]
gas_bin_sub = '_gas_bin_maps'
gas_spaxel_sub = '_gas_spaxel_maps'
stars_sub ='_kin_maps'
sfh_sub ='_sfh_maps'
version ='v3tk_v7.6.8'
#################################################

if __name__ == "__main__":
    gas_bin_image = pyfits.open(data_folder+'/'+galaxyid+'/'+galaxyid+gas_bin_sub+'.fits')
    gas_spaxel_image = pyfits.open(data_folder+'/'+galaxyid+'/'+galaxyid+gas_spaxel_sub+'.fits')
    stars_bin_image = pyfits.open(data_folder+'/'+galaxyid+'/'+galaxyid+stars_sub+'.fits')
    sfh_bin_image = pyfits.open(data_folder+'/'+galaxyid+'/'+galaxyid+sfh_sub+'.fits')
    mass_bin_image = open_mass_product(data_folder, galaxyid)

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
         'EBV',
         'A_R/2.5+logML_R',
        ]

    plot_id = []
    for j in range(0, 35, 5):
        plot_id.extend([(k, 'BIN', id[k]) for k in range(j, j+5)])
        plot_id.extend([(k, 'SPAXELS', id[k]) for k in range(j, j+5)])
    plot_id.extend([(j, 'BIN', id[j]) for j in range(35, 40)])
    plot_id.extend([(j, 'SPAXELS', id[j]) for j in range(35, 40)])
    plot_id.extend([(j, 'BIN', id[j]) for j in range(40, 45)])
    plot_id.extend([(j, 'BIN', id[j]) for j in range(45, 50)])
    plot_id.extend([(j, 'SPAXELS', id[j]) for j in range(45, 50)])
    plot_id.extend([(j, 'BIN', id[j]) for j in range(50, len(id))])

    gas_limits = {}
    for k in range(0, 35):
        if id[k].endswith('_VEL') or id[k].endswith('_SIGMA') or id[k].endswith('_SIGMA_ERR'):
            bin_data = gas_bin_image[id[k]].data
            spaxel_data = gas_spaxel_image[id[k]].data
            values = np.concatenate([bin_data[np.isfinite(bin_data)], spaxel_data[np.isfinite(spaxel_data)]])
            if len(values) == 0:
                continue
            if id[k].endswith('_VEL'):
                mean = np.nanmean(values)
                std = np.nanstd(values)
                gas_limits[id[k]] = (mean - 1 * std, mean + 1 * std)
            elif id[k].endswith('_SIGMA'):
                vmin, vmax = np.nanpercentile(values, [5, 95])
                gas_limits[id[k]] = (max(0, vmin), vmax)
            elif id[k].endswith('_SIGMA_ERR'):
                vmax = min(500, np.nanpercentile(values, 80))
                gas_limits[id[k]] = (1, max(30, vmax))


    with PdfPages(out) as pdf:
        for i in range(0, len(plot_id)):
            # Create a new figure every 5 plots (2x3 grid)
            if i % 5 == 0:
                fig, axs = plt.subplots(2, 3, figsize=(20, 12))
                axs = axs.flatten()  # Flatten to make it easier to index each subplot
            
            # Select the subplot based on the position within the 2x3 grid
            ax = axs[i % 5]

            # Clear any unused subplot for a clean layout
            if (i % 5) == 4 or i == len(plot_id) - 1:
                for j in range((i % 5) + 1, 6):
                    fig.delaxes(axs[j])

            id_i, level, id_name = plot_id[i]
            count = (i % 5) + 1
            img = None
            #######
            #######
            if (id_i<35):
                gas_image = gas_bin_image
                if level == 'SPAXELS':
                    gas_image = gas_spaxel_image
                # Retrieve data for the current plot
                to_plot = gas_image[id_name].data
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
                        vmin, vmax = gas_limits.get(id_name, (mean - 1 * std, mean + 1 * std))
                        img=ax.imshow(to_plot, origin='lower', vmin=vmin, vmax=vmax, cmap='PiYG_r')
                elif count == 4:
                    if (~np.isnan(mean)):
                        vmin, vmax = gas_limits.get(id_name, (50, 100))
                        img=ax.imshow(to_plot, origin='lower', vmin=vmin, vmax=vmax, cmap='magma')
                elif count == 5:
                    if (~np.isnan(mean)):
                        vmin, vmax = gas_limits.get(id_name, (1, 30))
                        img=ax.imshow(to_plot, origin='lower', vmin=vmin, vmax=vmax, cmap='cividis')
                    
                
                
            ######
            ######
            if ((id_i>=35) & (id_i<40)):
                gas_image = gas_bin_image
                if level == 'SPAXELS':
                    gas_image = gas_spaxel_image
                if count==1:
                    to_plot= gas_image['HA6562_FLUX'].data / gas_image['HB4861_FLUX'].data
                    img=ax.imshow(to_plot,origin='lower',norm=LogNorm(vmin=2.5,vmax=8),cmap='plasma')

                elif count==2:
                    to_plot= gas_image['NII6583_FLUX'].data / gas_image['HA6562_FLUX'].data
                    img=ax.imshow(to_plot,origin='lower',norm=LogNorm(vmin=0.1,vmax=1))
        
                elif count==3:
                    to_plot= gas_image['SII6716_FLUX'].data / gas_image['NII6583_FLUX'].data
                    mean=np.nanmedian(to_plot)
                    vmin=mean/2.
                    vmax=mean*2.
                    img=ax.imshow(to_plot,origin='lower',norm=LogNorm(vmin=vmin,vmax=vmax),cmap='magma')
                elif count==4:
                    to_plot= gas_image['HA6562_VEL'].data - gas_image['HB4861_VEL'].data
                    img=ax.imshow(to_plot,origin='lower',vmin=-3,vmax=+3,cmap='PiYG_r')
                elif count==5:
                    to_plot= gas_image['HA6562_SIGMA'].data / gas_image['HB4861_SIGMA'].data
                    img=ax.imshow(to_plot,origin='lower',vmin=0.6,vmax=1.2,cmap='PiYG')
    
           #######
           #######
            if ((id_i>=40) & (id_i<45)):
                to_plot= stars_bin_image[id_name].data
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
            
            if ((id_i>=45) & (id_i<50)):
                gas_image = gas_bin_image
                if level == 'SPAXELS':
                    gas_image = gas_spaxel_image
                if count==1:
                    if level == 'SPAXELS':
                        to_plot= gas_image['HA6562_VEL'].data - gas_image['V_STARS2'].data
                    else:
                        to_plot= gas_image['HA6562_VEL'].data - stars_bin_image['V'].data
                    img=ax.imshow(to_plot,origin='lower',vmin=-20,vmax=+20,cmap='PuOr_r')

                elif count==2:
                    if level == 'SPAXELS':
                        to_plot= gas_image['HA6562_SIGMA'].data - gas_image['SIGMA_STARS2'].data
                    else:
                        to_plot= gas_image['HA6562_SIGMA'].data - stars_bin_image['SIGMA'].data
                    img=ax.imshow(to_plot,origin='lower',vmin=-30, vmax=30, cmap='PRGn')
        
                elif count==3:
                    to_plot= gas_image['SII6716_FLUX'].data / gas_image['SII6730_FLUX'].data
                    img=ax.imshow(to_plot,origin='lower',vmin=0.6, vmax=1.7)

                elif count==4:
                    to_plot= gas_image['HA6562_VEL'].data - gas_image['NII6583_VEL'].data
                    img=ax.imshow(to_plot,origin='lower',vmin=-20,vmax=+20,cmap='BrBG_r')
        
                elif count==5:
                    to_plot= gas_image['HA6562_SIGMA'].data - gas_image['NII6583_SIGMA'].data
                    img=ax.imshow(to_plot,origin='lower',vmin=-20,vmax=+20,cmap='BrBG_r')
        
           #######
           #######
            if ((id_i>=50) & (id_i<55)):
                if count==1:
                    to_plot= sfh_bin_image[id_name].data
                    mean=np.nanmean(to_plot)
                    std=np.nanstd(to_plot)
                    img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='magma')
                elif count==2:
                    to_plot= sfh_bin_image[id_name].data
                    mean=np.nanmean(to_plot)
                    std=np.nanstd(to_plot)
                    img=ax.imshow(to_plot,origin='lower',vmin=mean-1*std,vmax=mean+1*std,cmap='plasma')
                elif count==3:
                    to_plot= sfh_bin_image[id_name].data
                    img=ax.imshow(to_plot,origin='lower',vmin=0,vmax=0.3,cmap='magma')
                elif count==4:
                    if mass_bin_image is not None:
                        to_plot = mass_correction_map(sfh_bin_image, mass_bin_image)
                        finite = to_plot[np.isfinite(to_plot)]
                        if len(finite) > 0:
                            vmin, vmax = np.nanpercentile(finite, [5, 95])
                            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
                                median = np.nanmedian(finite)
                                vmin, vmax = median - 0.2, median + 0.2
                            img=ax.imshow(to_plot,origin='lower',vmin=vmin,vmax=vmax,cmap='viridis')
                    if img is None:
                        ax.text(
                            0.5,
                            0.5,
                            "ML_R unavailable",
                            transform=ax.transAxes,
                            ha="center",
                            va="center",
                        )
       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            # Set title or other subplot details as needed
            title = id_name
            if id_i < 40 or ((id_i>=45) & (id_i<50)):
                title = level+' '+id_name
            ax.set_title(f"{title}")
            if img is not None:
                cbar = fig.colorbar(img, ax=ax, fraction=0.046, pad=0.04)
            # After every 5 plots, save the page and close the figure
            if (i + 1) % 5 == 0 or i == len(plot_id) - 1:
                plt.tight_layout()
                pdf.savefig(fig)  # Save the current page in the PDF
                plt.close(fig)  # Close the figure to free memory
