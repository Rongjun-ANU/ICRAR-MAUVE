from pathlib import Path
import re
import json
import hashlib
import unicodedata
from difflib import SequenceMatcher
from bs4 import BeautifulSoup
import fitz

ROOT=Path('/Users/Igniz/Desktop/ICRAR/MAUVE')
ASSETS=ROOT/'assets/20260913_jwst_close_revision'
data=json.loads((ASSETS/'content.json').read_text())
source=json.loads((ASSETS/'source.json').read_text())
pdf=ROOT/(data['stem']+'.pdf')
md=ROOT/(data['stem']+'.md')
doc=fitz.open(pdf)
original=ROOT/'_JWST_Cycle5__MAUVE_NIRCam_MIRI_Imaging.pdf'
assert hashlib.sha256(original.read_bytes()).hexdigest()==source['source_sha256']
assert len(doc)==12

def plain(t):
    return BeautifulSoup(t,'html.parser').get_text()
def norm(t):
    return re.sub(r'\s+', '', unicodedata.normalize('NFKC',t))

spans=[]
actual=[]
bounds=[]
for n,page in enumerate(doc,1):
    for block in page.get_text('dict')['blocks']:
        for line in block.get('lines',[]):
            for span in line['spans']:
                if span['bbox'][1]>738: continue
                actual.append(span['text'])
                spans.append((n,span))
                x0,y0,x1,y1=span['bbox']
                if not (x0>=71.5 and x1<=540.5 and y0>=68 and y1<=733): bounds.append((n,span))
assert not bounds, bounds
expected=''.join(plain(b['text']) for b in data['blocks'])
assert norm(expected)==norm(''.join(actual)), 'PDF text differs from Markdown content model'
assert '\ufffd' not in ''.join(actual)
colors=sorted(set(s['color'] for _,s in spans))
assert colors==[0,int('147D3B',16)], colors

# Ensure the entire new Goal 3 and all revised duplication prose are physically green.
for n,s in spans:
    if 'Young stellar feedback in an externally disturbed ISM.' in s['text']:
        assert s['color']==int('147D3B',16)

expected_fig_hash=sorted(f['sha256'] for f in source['figures'])
actual_fig_hash=[]
for page in doc:
    for img in page.get_images(full=True):
        b=doc.extract_image(img[0])['image']
        actual_fig_hash.append(hashlib.sha256(b).hexdigest())
assert sorted(actual_fig_hash)==expected_fig_hash, actual_fig_hash

overlap=['NGC 4294','NGC 4351','NGC 4388','NGC 4402','NGC 4548','NGC 4569','NGC 4579','NGC 4606']
dup=next(plain(b['text']) for b in data['blocks'] if 'J-Virgo GO 7763 [81] overlaps eight targets:' in b['text'])
assert re.findall(r'NGC \d+',dup)==overlap
assert 'Goal #3. Precise distance' not in md.read_text()
assert 'orange' not in md.read_text().lower()
assert md.read_text().count('![Original Cycle 5 Figure')==3

# Quantify word preservation of Goals 1 and 2, excluding their figure captions.
preservation={}
for g,end in [(1,'Fig. 2:'),(2,'The sharp, sensitive PAH maps')]:
    start='Goal #%s.'%g
    old=source['text'][source['text'].index(start):source['text'].index(end,source['text'].index(start))]
    new=plain(next(b['text'] for b in data['blocks'] if plain(b['text']).startswith(start)))
    a,b=plain(old).split(),new.split()
    matching=sum(m.size for m in SequenceMatcher(None,a,b,autojunk=False).get_matching_blocks())
    preservation['goal_%s_opening_original_words_retained_fraction'%g]=round(matching/len(a),4)

result={'status':'PASS','pages':len(doc),'source_unchanged_sha256':source['source_sha256'],
        'original_figures_preserved_byte_for_byte':3,'matching_md_pdf_content':True,
        'text_colors':[hex(c) for c in colors],'green_text_spans':sum(s['color']!=0 for _,s in spans),
        'margin_overflows':len(bounds),'references':87,'overlap':overlap,
        **preservation,'output_sha256':{'md':hashlib.sha256(md.read_bytes()).hexdigest(),
                                     'pdf':hashlib.sha256(pdf.read_bytes()).hexdigest()},
        'limits':'Close editorial revision. No new ETC/APT, footprint validation, or scientific re-audit of unchanged Cycle 5 claims.'}
(ASSETS/'verification.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
