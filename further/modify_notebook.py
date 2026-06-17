import json

with open('20260617_Moran_MAUVE_WHAN.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        
        # Replacement 1
        old_1 = "        ha_sigma[sigma_ok] = np.sqrt(sigma2_intrinsic[sigma_ok])\n        if 'O_H_O3N2_M14_HII' in h_gas:\n"
        new_1 = """        ha_sigma[sigma_ok] = np.sqrt(sigma2_intrinsic[sigma_ok])
        
        # Load NII and Ha corrected fluxes for WHAN diagram
        ha_corr = h_gas['HA6562_FLUX_CORR'].data
        nii_corr = h_gas['NII6583_FLUX_CORR'].data
        
        if 'O_H_O3N2_M14_HII' in h_gas:\n"""
        
        # Replacement 2
        old_2 = "    sigSFR = np.where(ew_mask, sigSFR, np.nan)\n    oh_map = np.where(ew_mask, oh_map, np.nan)\n    return sigM, sigSFR, oh_map, indicator_name, wcs, gas_header, binid_map\n"
        new_2 = """    sigSFR = np.where(ew_mask, sigSFR, np.nan)
    oh_map = np.where(ew_mask, oh_map, np.nan)
    return sigM, sigSFR, oh_map, indicator_name, wcs, gas_header, binid_map, ew_ha, ha_corr, nii_corr\n"""
        
        # Replacement 3
        old_3_start = "for gal in galaxies:\n    sigM, sigSFR, oh_o3n2_map, indicator_name, wcs_full, gas_header, binid_map = load_maps(gal)\n    # Valid pixels across maps"
        
        if old_3_start in source:
            start_idx = source.find("for gal in galaxies:\n    sigM, sigSFR, oh_o3n2_map, indicator_name, wcs_full, gas_header, binid_map = load_maps(gal)\n")
            end_idx = source.find("    plt.show()") + len("    plt.show()")
            
            new_3 = """for gal in galaxies:
    sigM, sigSFR, oh_o3n2_map, indicator_name, wcs_full, gas_header, binid_map, ew_ha, ha_corr, nii_corr = load_maps(gal)
    
    # Calculate WHAN diagram components
    valid_whan = (
        np.isfinite(ew_ha) & (ew_ha > 0) & 
        np.isfinite(ha_corr) & (ha_corr > 0) & 
        np.isfinite(nii_corr) & (nii_corr > 0) & 
        np.isfinite(sigSFR) & np.isfinite(sigM)
    )
    
    ew_ha_valid = ew_ha[valid_whan]
    log_nii_ha = np.log10(nii_corr[valid_whan] / ha_corr[valid_whan])
    log_ssfr = sigSFR[valid_whan] - sigM[valid_whan]
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Increase tick label font sizes and tick length/width specifically for WHAN
    ax.tick_params(axis='both', which='major', labelsize=18, size=6, width=1.2, direction='in')
    ax.tick_params(axis='both', which='minor', size=4, width=1.0, direction='in')
    
    scatter = ax.scatter(
        log_nii_ha, 
        ew_ha_valid, 
        c=log_ssfr, 
        cmap='viridis',
        s=10, 
        alpha=0.6,
        marker='d',
        vmin=-12.0,
        vmax=-9.0
    )
    
    ax.set_yscale('log')
    ax.set_ylim(0.2, 400)
    ax.set_xlim(-1.3, 0.45)
    
    ax.set_xlabel(r'$\\log_{10}([\\mathrm{NII}]/\\mathrm{H}\\alpha)$', fontsize=22)
    ax.set_ylabel(r'$\\mathrm{EW}_{\\mathrm{H}\\alpha} (\\AA)$', fontsize=22)
    
    # Draw WHAN boundaries
    ax.axhline(1.0, linestyle=':', color='blue', label='EW = 1Å', zorder=10)
    ax.plot([-0.4, 0.45], [6.0, 6.0], linestyle='-.', color='blue', label='EW = 6Å', zorder=10)
    ax.plot([-0.4, -0.4], [1.0, 400], linestyle='-', color='blue', label='S06', zorder=10)
    
    ax.text(-0.85, 60, 'SF', fontsize=18, fontweight='bold', ha='center')
    ax.text(0.15, 40, 'Seyferts', fontsize=18, fontweight='bold', ha='center')
    ax.text(0.15, 3.5, 'LINERs', fontsize=18, fontweight='bold', ha='center')
    ax.text(-0.05, 0.35, 'Passive', fontsize=18, fontweight='bold', ha='center')
    
    cbar = fig.colorbar(scatter, ax=ax, pad=0.05)
    cbar.set_label(r'$\\log_{10}(\\mathrm{sSFR} / \\mathrm{yr}^{-1})$', fontsize=22)
    cbar.ax.tick_params(labelsize=18)
    
    ax.legend(loc='lower left', fontsize=12)
    
    fig.savefig(f'Moran_MAUVE/{gal}_WHAN_diagram.png', dpi=600, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f'Moran_MAUVE/{gal}_WHAN_diagram.pdf', bbox_inches="tight", pad_inches=0.02)
    plt.show()"""
            source = source[:start_idx] + new_3 + source[end_idx:]
            
        source = source.replace(old_1, new_1)
        source = source.replace(old_2, new_2)
        
        # Split back to lines and preserve newlines correctly
        lines = []
        import re
        parts = re.split('(\n)', source)
        for i in range(0, len(parts)-1, 2):
            lines.append(parts[i] + parts[i+1])
        if parts[-1]:
            lines.append(parts[-1])
            
        cell['source'] = lines

with open('20260617_Moran_MAUVE_WHAN.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

print("Modification complete.")
