import json

with open('20260617_Moran_MAUVE_WHAN.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        
        if 'ew_mask = np.isfinite(ew_ha)' in source:
            old_mask = """    ew_ha = align_map_to_shape(ew_ha, sigSFR.shape, label=f"{gal} {EW_HDU}")
    ew_mask = np.isfinite(ew_ha) & (ew_ha > EW_MIN)
    ew_mask = ew_mask & np.isfinite(ha_sigma) & (ha_sigma < 45.0)
    sigSFR = np.where(ew_mask, sigSFR, np.nan)
    oh_map = np.where(ew_mask, oh_map, np.nan)
    return sigM, sigSFR, oh_map, indicator_name, wcs, gas_header, binid_map, ew_ha, ha_corr, nii_corr"""
            
            new_mask = """    ew_ha = align_map_to_shape(ew_ha, sigSFR.shape, label=f"{gal} {EW_HDU}")
    # Removed EW and velocity dispersion selections for HII spaxels
    # ew_mask = np.isfinite(ew_ha) & (ew_ha > EW_MIN)
    # ew_mask = ew_mask & np.isfinite(ha_sigma) & (ha_sigma < 45.0)
    # sigSFR = np.where(ew_mask, sigSFR, np.nan)
    # oh_map = np.where(ew_mask, oh_map, np.nan)
    return sigM, sigSFR, oh_map, indicator_name, wcs, gas_header, binid_map, ew_ha, ha_corr, nii_corr"""
            
            source = source.replace(old_mask, new_mask)
            
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
