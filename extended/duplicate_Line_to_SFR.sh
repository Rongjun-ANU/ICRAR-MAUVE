#!/usr/bin/env bash
# duplicate_Line_to_SFR_NGC4298.sh  –  clone & retarget Line_to_SFR_NGC4298.ipynb
# Usage:  bash duplicate_Line_to_SFR_NGC4298.sh

set -euo pipefail

template="Line_to_SFR_IC3392.ipynb"   # master file
orig_tag="IC3392"                     # string to replace

# List of destination galaxy IDs (edit here if needed)
galaxies=(  
  NGC4064  
  NGC4192  
  NGC4293  
  NGC4298
  NGC4330
  NGC4383 
  NGC4396  
  NGC4419  
  NGC4457  
  NGC4501
  NGC4522 
  NGC4694  
  NGC4698
)

for gal in "${galaxies[@]}"; do
    new_file="Line_to_SFR_${gal}.ipynb"

    # 1. Copy the notebook
    cp -- "$template" "$new_file"

    # 2. In-place substitution (Perl is cross-platform & JSON-safe)
    #    -pi  : edit file in place
    #    -e   : execute the given Perl code
    perl -pi -e "s/${orig_tag}/${gal}/g" "$new_file"

    echo "✓  Created $new_file"
done


