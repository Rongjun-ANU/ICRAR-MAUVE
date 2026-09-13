from pathlib import Path
import fitz
import re
import json
import hashlib

ROOT = Path('/Users/Igniz/Desktop/ICRAR/MAUVE')
ASSETS = ROOT / 'assets/20260913_jwst_close_revision'
SOURCE = ROOT / '_JWST_Cycle5__MAUVE_NIRCam_MIRI_Imaging.pdf'
doc = fitz.open(SOURCE)
pages = []
for page in doc:
    text = page.get_text()
    text = re.sub(r'\n\d+\s*$', '', text)
    pages.append(text)
raw = '\n'.join(pages)
print('Line-end hyphens:', re.findall(r'\S+-\n\S+', raw))
# Preserve lexical hyphens while removing typographic line-end hyphenation.
keep = {'mass-matched', 'UV-optical', 'star-forming', 'line-based',
        'ICM-affected', 'off-target', 'on-target', 'mid-IR', 'peer-review',
        'medium-and', 'high-level', 'continuum-subtracted', 'small-scale'}
def unhyphen(m):
    a, b = m.group(1), m.group(2)
    return a + ('-' if a + '-' + b in keep else '') + b
raw = re.sub(r'([A-Za-z]+)-\n([A-Za-z]+)', unhyphen, raw)
raw = raw.replace('\ufb01', 'fi').replace('\ufb02', 'fl')
text = re.sub(r'\s+', ' ', raw).strip()
for a, b in [('Hii', 'H II'), ('HII', 'H II'), ('Hi,', 'H I,'), ('Hi ', 'H I '), ('Hi.', 'H I.'),
             ('M⊙clusters', 'M⊙ clusters'), ('103 M⊙','10<sup>3</sup> M⊙'),
             ('pc−2', 'pc<sup>−2</sup>'), ('cm2/arcsec2', 'cm<sup>2</sup>/arcsec<sup>2</sup>'),
             ('10−16', '10<sup>−16</sup>'), ('10−17', '10<sup>−17</sup>'),
             ('AV ', 'A<sub>V</sub> '), ('E(B−V)', 'E(B−V)'),
             ('∼<', '≲'), ('∼>', '≳')]:
    text = text.replace(a, b)
# Extract the three embedded original figure files, without redrawing or altering them.
figures = []
for page_index, figure_ids in [(2, [1, 2]), (4, [3])]:
    imgs = sorted(doc[page_index].get_image_info(xrefs=True), key=lambda im: im['bbox'][1])
    for number, item in zip(figure_ids, imgs):
        image = doc.extract_image(item['xref'])
        path = ASSETS / ('original_figure_%d.' % number + image['ext'])
        path.write_bytes(image['image'])
        figures.append({'number':number, 'path':str(path), 'source_page':page_index+1,
                        'xref':item['xref'], 'width':image['width'], 'height':image['height'],
                        'sha256':hashlib.sha256(image['image']).hexdigest()})
(ASSETS / 'source.json').write_text(json.dumps({'text':text, 'figures':figures,
    'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest()}, ensure_ascii=False, indent=2))
print('Extracted original figures:', figures)
