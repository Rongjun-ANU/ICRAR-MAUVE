"""Check source preservation, report links and PDF content; build QA contact sheets."""
from pathlib import Path
import hashlib
import json
import re
from urllib.parse import unquote
import html
import pdfplumber
from pypdf import PdfReader
from PIL import Image, ImageDraw

ASSET = Path(__file__).resolve().parent
ROOT = ASSET.parents[1]
SOURCE = ROOT / '20261011_MAUVE_Paper_II_Findings_Synthesis_and_Paper_Plan.md'
PDF = SOURCE.with_suffix('.pdf')
SCRATCH = Path('/private/tmp/mauve_second_paper_20261011')
SCRATCH.mkdir(exist_ok=True)
md = SOURCE.read_text()
notebooks = json.loads((ASSET/'notebook_inventory.json').read_text())
reports = json.loads((ASSET/'report_sources.json').read_text())
for row in notebooks:
    p = ROOT.parent/'further'/row['name']
    assert hashlib.sha256(p.read_bytes()).hexdigest() == row['sha256'], str(p)
for row in reports:
    p = Path(row['path'])
    assert hashlib.sha256(p.read_bytes()).hexdigest() == row['sha256'], str(p)
links = re.findall(r'\]\(([^)]+)\)', md)
local = []
for target in links:
    if target.startswith(('https://','http://','#')):
        continue
    p = Path(unquote(target))
    p = p if p.is_absolute() else ROOT/p
    assert p.exists(), str(p)
    local.append(str(p))
assert '<!--' not in md
reader = PdfReader(PDF)
page_text = [p.extract_text() or '' for p in reader.pages]
text = '\n'.join(page_text)
assert all(len(t) > 100 for t in page_text)
assert '\ufffd' not in text
for row in notebooks + reports:
    assert row['id'] in text, row['id']
for needle in ['0.00004043', '-0.243', '0.494', '-1.204', '0.223',
               '7000 Angstrom', 'Appendix C. Verified external context']:
    assert needle in text, needle
# Paragraphs can cross pages, where headers intervene. Check every short
# segment separately and report unmatched segments for a manual review.
plain_lines = []
for line in md.splitlines():
    if not line.strip() or line.startswith('|'):
        continue
    line = re.sub(r'!?\[([^]]+)\]\([^)]+\)', r'\1', line)
    line = re.sub(r'^[#>\- ]+', '', line)
    line = line.replace('**','').replace('*','').replace('`','')
    plain_lines.append(line)
normalize = lambda s: re.sub(r'\s+', '', html.unescape(s))
norm_pdf = normalize(text)
segments = [normalize(line)[i:i+90] for line in plain_lines
            for i in range(0,len(normalize(line)),90) if len(normalize(line)[i:i+90]) > 12]
unmatched = [s for s in segments if s not in norm_pdf]
coverage = 1-len(unmatched)/len(segments)
assert coverage > .96, (coverage, unmatched[:8])
outside = []
cropped_text = []
orphan_headings = []
headings = [normalize(re.sub(r'^#+ ', '', line)) for line in md.splitlines() if line.startswith('##')]
with pdfplumber.open(PDF) as pdf:
    for pi,page in enumerate(pdf.pages,1):
        body = page.crop((0,44,page.width,page.height-44)).extract_text() or ''
        cropped_text.append(body)
        tail = normalize(body.splitlines()[-1]) if body.splitlines() else ''
        if len(tail) > 15 and any(h.endswith(tail) for h in headings):
            orphan_headings.append([pi,tail])
        for ch in page.chars:
            if ch['x0'] < -0.1 or ch['x1'] > page.width+.1 or ch['top'] < -.1 or ch['bottom'] > page.height+.1:
                outside.append([pi,ch['text'],ch['x0'],ch['top']])
assert not outside, outside[:10]
assert not orphan_headings, orphan_headings
norm_pdf = normalize(''.join(cropped_text))
unmatched = [s for s in segments if s not in norm_pdf]
coverage = 1-len(unmatched)/len(segments)
assert coverage > .99, (coverage, unmatched[:8])
images = sorted(SCRATCH.glob('report_page-*.png'))
assert len(images) == len(reader.pages), (len(images),len(reader.pages))
contact_files = []
for start in range(0,len(images),9):
    chosen=images[start:start+9]
    cellw,cellh=500,730
    canvas=Image.new('RGB',(3*cellw,((len(chosen)+2)//3)*cellh),'#D9E1E6')
    draw=ImageDraw.Draw(canvas)
    for j,path in enumerate(chosen):
        im=Image.open(path).convert('RGB'); im.thumbnail((490,700))
        x=(j%3)*cellw+(cellw-im.width)//2
        y=(j//3)*cellh+23
        canvas.paste(im,(x,y));draw.text(((j%3)*cellw+12,(j//3)*cellh+5),f'Page {start+j+1}',fill='black')
    out=SCRATCH/f'qa_contact_{start//9+1}.jpg'
    canvas.save(out,quality=90);contact_files.append(str(out))
record = dict(status='PASS', pages=len(reader.pages), notebook_hashes_unchanged=len(notebooks),
              report_hashes_unchanged=len(reports), local_links_checked=len(local),
              body_text_segment_coverage=coverage, unmatched_segments=unmatched,
              off_page_characters=len(outside), orphan_headings=len(orphan_headings), rendered_pages=len(images),
              contact_sheets=contact_files, markdown_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              pdf_sha256=hashlib.sha256(PDF.read_bytes()).hexdigest())
(ASSET/'verification.json').write_text(json.dumps(record,indent=2))
(SCRATCH/'report_pdf_text.txt').write_text(text)
print(json.dumps({k:v for k,v in record.items() if k not in ['unmatched_segments','contact_sheets']},indent=2))
print('UNMATCHED_SEGMENTS',len(unmatched),unmatched[:8])
