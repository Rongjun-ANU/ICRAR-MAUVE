import json

with open('20260617_Moran_MAUVE_WHAN.ipynb', 'r') as f:
    nb = json.load(f)

# Step 1: Remove O_H from the notebook
for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        
        # Remove oh check in physics_ready_galaxies
        if 'and fits_has_hdus(gas_path, [\'O_H_O3N2_M13_HII\', \'O_H_O3N2_M13_SF\'])' in source:
            source = source.replace('and fits_has_hdus(gas_path, [\'O_H_O3N2_M13_HII\', \'O_H_O3N2_M13_SF\'])', '')
            source = source.replace('mass/SFR/OH+EW products', 'mass/SFR+EW products')

        # Remove O_H loading from load_maps
        oh_loading_block = """        if 'O_H_O3N2_M14_HII' in h_gas:
            oh_map = h_gas['O_H_O3N2_M14_HII'].data
            indicator_name = 'O3N2-M14'
        elif 'O_H_O3N2_M13_HII' in h_gas:
            print(f"Warning: O3N2-M14 not available for {gal}; using O3N2-M13 instead.")
            oh_map = h_gas['O_H_O3N2_M13_HII'].data
            indicator_name = 'O3N2-M13'
        else:
            available = [hdu.name for hdu in h_gas]
            raise KeyError(f"No O3N2 metallicity extension (M14/M13) found for {gal}. Available: {available}")"""
        
        if oh_loading_block in source:
            source = source.replace(oh_loading_block, "")
            source = source.replace('return sigM, sigSFR, oh_map, indicator_name, wcs, gas_header, binid_map, ew_ha, ha_corr, nii_corr', 
                                  'return sigM, sigSFR, wcs, gas_header, binid_map, ew_ha, ha_corr, nii_corr')
        
        # Remove O_H unpacking in the plotting loop
        if 'sigM, sigSFR, oh_o3n2_map, indicator_name, wcs_full, gas_header, binid_map, ew_ha, ha_corr, nii_corr = load_maps(gal)' in source:
            source = source.replace('sigM, sigSFR, oh_o3n2_map, indicator_name, wcs_full, gas_header, binid_map, ew_ha, ha_corr, nii_corr = load_maps(gal)', 
                                  'sigM, sigSFR, wcs_full, gas_header, binid_map, ew_ha, ha_corr, nii_corr = load_maps(gal)')

        # Split back to lines and preserve newlines correctly
        lines = []
        import re
        parts = re.split('(\n)', source)
        for i in range(0, len(parts)-1, 2):
            lines.append(parts[i] + parts[i+1])
        if parts[-1]:
            lines.append(parts[-1])
            
        cell['source'] = lines

# Step 2: Add the new cell for spatial WHAN maps
new_cell_source = """import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

for gal in galaxies:
    sigM, sigSFR, wcs_full, gas_header, binid_map, ew_ha, ha_corr, nii_corr = load_maps(gal)
    
    # Calculate WHAN diagram components
    valid_whan = (
        np.isfinite(ew_ha) & (ew_ha > 0) & 
        np.isfinite(ha_corr) & (ha_corr > 0) & 
        np.isfinite(nii_corr) & (nii_corr > 0) & 
        np.isfinite(sigSFR) & np.isfinite(sigM)
    )
    
    log_nii_ha = np.full_like(ew_ha, np.nan)
    log_nii_ha[valid_whan] = np.log10(nii_corr[valid_whan] / ha_corr[valid_whan])
    
    # WHAN classification map
    whan_class = np.full(ew_ha.shape, np.nan)
    
    # Passive: -1, SF: 0, LINERs: 1, Seyferts: 2
    whan_class[valid_whan & (ew_ha < 1.0)] = -1
    whan_class[valid_whan & (log_nii_ha < -0.4) & (ew_ha >= 1.0)] = 0
    whan_class[valid_whan & (log_nii_ha >= -0.4) & (ew_ha >= 1.0) & (ew_ha < 6.0)] = 1
    whan_class[valid_whan & (log_nii_ha >= -0.4) & (ew_ha >= 6.0)] = 2
    
    # Setup coordinates
    wcs_celestial = WCS(gas_header).celestial
    y_size, x_size = sigSFR.shape
    x_coords = np.arange(x_size)
    y_coords = np.arange(y_size)
    xx, yy = np.meshgrid(x_coords, y_coords)
    ra, dec = wcs_celestial.pixel_to_world_values(xx, yy)
    extent = [ra.max(), ra.min(), dec.min(), dec.max()]
    
    # Scale bar
    pix_scale_arcsec = 0.2
    distance_mpc = 16.5
    kpc_arcsec = 206265.0 / (distance_mpc * 1000.0)
    scalebar_length_deg = kpc_arcsec / 3600.0
    
    def add_scalebar(ax, length_deg, label='1 kpc'):
        x0, x1 = ax.get_xlim()
        x_range = abs(x1 - x0)
        if x_range == 0:
            return
        frac_len = length_deg / x_range
        x_start = 0.05
        y_start = 0.95
        x_end = x_start + frac_len
        line, = ax.plot([x_start, x_end], [y_start, y_start], transform=ax.transAxes,
                        color='k', lw=3, clip_on=False)
        line.set_in_layout(False)
        txt = ax.text((x_start + x_end) / 2, y_start - 0.03, label,
                      transform=ax.transAxes, ha='center', va='top', color='k')
        txt.set_in_layout(False)
        
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Set up ListedColormap
    cmap_whan = mcolors.ListedColormap(['red', 'blue', 'green', 'orange'])
    bounds = [-1.5, -0.5, 0.5, 1.5, 2.5]
    norm_whan = mcolors.BoundaryNorm(bounds, cmap_whan.N)
    
    im = ax.imshow(
        np.ma.masked_invalid(whan_class),
        origin="lower", cmap=cmap_whan, norm=norm_whan,
        extent=extent, aspect="equal"
    )
    
    ax.xaxis.set_major_locator(mticker.MaxNLocator(3))
    ax.yaxis.set_major_locator(mticker.MaxNLocator(5))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
    
    ax.set_xlabel("RA (deg)", fontsize=16)
    ax.set_ylabel("Dec (deg)", fontsize=16)
    
    add_scalebar(ax, scalebar_length_deg, label="1 kpc")
    
    # Custom colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, ticks=[-1, 0, 1, 2])
    cbar.ax.set_yticklabels(['Passive', 'SF', 'LINERs', 'Seyferts'], fontsize=14)
    
    ax.set_title(f"{gal} WHAN Classification", fontsize=18)
    
    fig.savefig(f'Moran_MAUVE/{gal}_WHAN_map.png', dpi=600, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f'Moran_MAUVE/{gal}_WHAN_map.pdf', bbox_inches="tight", pad_inches=0.02)
    plt.show()"""

new_cell = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [line + '\n' for line in new_cell_source.split('\n')]
}
# Remove the last newline
new_cell['source'][-1] = new_cell['source'][-1].rstrip('\n')

nb['cells'].append(new_cell)

with open('20260617_Moran_MAUVE_WHAN.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

print("Modification complete.")
