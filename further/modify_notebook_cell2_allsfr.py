import json

with open('20260617_Moran_MAUVE_WHAN.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        
        old_str = "sigSFR = h_gas['LOGSFR_SURFACE_DENSITY_HII'].data"
        new_str = "sigSFR = h_gas['LOGSFR_SURFACE_DENSITY'].data"
        
        if old_str in source:
            source = source.replace(old_str, new_str)
            
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
