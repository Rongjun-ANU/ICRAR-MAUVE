"""Rerun previously audited scalar extraction, without source/product writes."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

OUT = Path(__file__).resolve().parent
OLD = OUT.parent / '20260911_sf_nsf_physical_model'
ROOT = Path('/Users/Igniz/Desktop/ICRAR/further')
records = []
for item in json.loads((OLD / 'input_fingerprints.json').read_text()):
    path = Path(item['path'])
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    records.append(dict(path=str(path), sha256=digest, bytes=path.stat().st_size,
                        matches_11_September=digest == item['sha256']))
for name in ['inspect_live.py', 'inspect_halpha.py']:
    source = (OLD / name).read_text().replace('/private/tmp/mauve_20260911', str(OUT))
    adapted = OUT / name
    adapted.write_text(source)
    subprocess.run([sys.executable, str(adapted)], check=True)
# Reuse only observational estimates/figures; exclude the old illustrative model.
source = (OLD / 'reproduce_analysis.py').read_text().split('# Exact two-reservoir solution')[0]
source = source.replace('20260911', '20260914')
adapted = OUT / 'observational_bootstrap.py'
adapted.write_text(source)
subprocess.run([sys.executable, str(adapted)], check=True)
(OUT / 'input_fingerprints.json').write_text(json.dumps(records, indent=2))
print('CURRENT_SCALAR_EXTRACTION_AND_GALAXY_BOOTSTRAP_PASS')
