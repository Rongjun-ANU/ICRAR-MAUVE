from pathlib import Path
import json
import re
from html import escape

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, PageBreak, Image, Table, TableStyle, KeepTogether

ROOT = Path('/Users/Igniz/Desktop/ICRAR/MAUVE')
ASSETS = ROOT / 'assets/20260913_jwst_close_revision'
data = json.loads((ASSETS / 'content.json').read_text())
fontroot = Path('/Users/Igniz/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype')
for name, file in [('Serif','LiberationSerif-Regular.ttf'), ('Serif-Bold','LiberationSerif-Bold.ttf'),
                   ('Serif-Italic','LiberationSerif-Italic.ttf'), ('Serif-BoldItalic','LiberationSerif-BoldItalic.ttf'),
                   ('Math','DejaVuSans.ttf')]:
    pdfmetrics.registerFont(TTFont(name, str(fontroot/file)))
pdfmetrics.registerFontFamily('Serif', normal='Serif', bold='Serif-Bold', italic='Serif-Italic', boldItalic='Serif-BoldItalic')

body = ParagraphStyle('body',fontName='Serif', fontSize=12, leading=14.5,
                      spaceAfter=6, alignment=TA_JUSTIFY, allowWidows=0, allowOrphans=0)
heading = ParagraphStyle('heading', parent=body, fontName='Serif-Bold', fontSize=14,
                         leading=16.5, spaceBefore=8, spaceAfter=8, keepWithNext=True, alignment=TA_LEFT)
caption = ParagraphStyle('caption', parent=body, fontSize=11.5, leading=13.5, spaceAfter=8)
cover = ParagraphStyle('cover',parent=body,spaceAfter=15)
cover_note = ParagraphStyle('cover_note',parent=body, fontSize=11, spaceAfter=24, alignment=TA_LEFT)
reference = ParagraphStyle('reference',parent=body, fontSize=11.5, leading=14.5, spaceAfter=8,
                           leftIndent=23, firstLineIndent=-23, alignment=TA_LEFT)
disclosure = ParagraphStyle('disclosure',parent=body,fontSize=10,leading=12.5,spaceBefore=12,alignment=TA_LEFT)

def xml(t):
    t = re.sub(r'<span style="color: (#[A-Fa-f0-9]+);">', r'<font color="\1">', t)
    t = t.replace('</span>','</font>')
    # Retain the intentionally authored inline tags; escape literal comparison operators.
    parts = re.split(r'(<(?:/?(?:b|i|sup|sub|font|a)\b)[^>]*>)', t)
    for i, part in enumerate(parts):
        if part.startswith('<') and re.match(r'</?(?:b|i|sup|sub|font|a)\b', part): continue
        part = escape(part, quote=False)
        # Use a dedicated mathematical font only for glyphs absent from the text family.
        for char in set(part):
            if ord(char)>127 and ord(char) not in pdfmetrics.getFont('Serif').face.charWidths:
                assert ord(char) in pdfmetrics.getFont('Math').face.charWidths, repr(char)
                part=part.replace(char, '<font name="Math">'+char+'</font>')
        parts[i] = part
    return ''.join(parts)

class Proposal(BaseDocTemplate):
    def afterFlowable(self, flowable):
        if isinstance(flowable, Paragraph):
            plain=flowable.getPlainText()
            if flowable.style.name=='heading' or plain.startswith('Goal #'):
                key='section-'+str(self.seq.nextf('section'))
                self.canv.bookmarkPage(key)
                title=plain.split('. We will')[0] if plain.startswith('Goal #') else plain
                self.canv.addOutlineEntry(title,key,level=0,closed=False)

story=[]
for b in data['blocks']:
    k,t=b['kind'],b['text']
    if k=='pagebreak': story.append(PageBreak())
    elif k=='h2':
        if t.startswith('Special Requirements'): story.append(PageBreak())
        story.append(Paragraph(xml(t),heading))
    elif k=='figure':
        if b['number']==1:
            story.append(PageBreak())
            w=201.22
            im=Image(b['path'],width=w,height=w*b['height']/b['width'])
            table=Table([[im, Paragraph(xml(t),caption)]],colWidths=[213.22,254.78])
            table.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'MIDDLE'),
                ('LEFTPADDING',(0,0),(-1,-1),0),('RIGHTPADDING',(0,0),(-1,-1),0),
                ('TOPPADDING',(0,0),(-1,-1),0),('BOTTOMPADDING',(0,0),(-1,-1),0)]))
            story.append(table)
            story.append(Spacer(1,10))
        else:
            im=Image(b['path'],width=468,height=468*b['height']/b['width'])
            story.append(KeepTogether([Spacer(1,3),im,Spacer(1,3),Paragraph(xml(t),caption)]))
            if b['number']==2: story.append(PageBreak())
    else:
        sty={'cover':cover,'cover_note':cover_note,'reference':reference,'disclosure':disclosure}.get(k,body)
        story.append(Paragraph(xml(t),sty))

def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont('Serif',11)
    canvas.setFillColor(colors.black)
    canvas.drawCentredString(306,39,str(doc.page-1))
    canvas.restoreState()

path=ROOT/(data['stem']+'.pdf')
doc=Proposal(str(path), pagesize=(612,792),leftMargin=72,rightMargin=72,
             topMargin=72,bottomMargin=61,title='The Evolving ISM and Stellar Populations in Virgo Cluster Galaxies',
             author='MAUVE proposal working draft',subject='13 September 2026 close revision; updates in green')
frame=Frame(72,61,468,659,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0,id='main')
doc.addPageTemplates([PageTemplate(id='proposal',frames=[frame],onPage=footer)])
doc.build(story)
print('Rendered',path)
