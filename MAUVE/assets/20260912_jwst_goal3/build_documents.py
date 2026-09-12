"""Render the two canonical Markdown documents with bundled ReportLab Python.

Run build_figures.py with ICRAR Python first. PDF text is parsed from the
Markdown source; equations are rendered by the figure builder.
"""
from pathlib import Path
import html
import json
import re
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import BaseDocTemplate, Frame, PageTemplate, Paragraph, Spacer, Image, LongTable, TableStyle, KeepTogether

ASSETS = Path(__file__).resolve().parent
ROOT = ASSETS.parent.parent
SCRATCH = Path('/private/tmp/mauve_goal3_20260912')
DATA = json.loads((SCRATCH/'build_data.json').read_text())
W,H=A4
MARGIN=47
WIDTH=W-2*MARGIN
INK=colors.HexColor('#203945')
TEAL=colors.HexColor('#245e73')
MUTED=colors.HexColor('#5f7380')

styles={
 'title':ParagraphStyle('Title',fontName='Helvetica-Bold',fontSize=23,leading=27,textColor=TEAL,spaceAfter=16),
 'body':ParagraphStyle('Body',fontName='Helvetica',fontSize=10.3,leading=14.4,textColor=INK,spaceAfter=8),
 'h2':ParagraphStyle('Heading2',fontName='Helvetica-Bold',fontSize=15,leading=18.4,textColor=TEAL,spaceBefore=15,spaceAfter=9,keepWithNext=True),
 'h3':ParagraphStyle('Heading3',fontName='Helvetica-Bold',fontSize=11.9,leading=15,textColor=INK,spaceBefore=10,spaceAfter=7,keepWithNext=True),
 'caption':ParagraphStyle('Caption',fontName='Helvetica-Oblique',fontSize=9,leading=12.2,textColor=MUTED,spaceAfter=11),
 'cell':ParagraphStyle('Cell',fontName='Helvetica',fontSize=8.5,leading=11.5,textColor=INK,spaceAfter=0),
 'th':ParagraphStyle('TableHead',fontName='Helvetica-Bold',fontSize=8.5,leading=11.5,textColor=colors.white),
 'ref':ParagraphStyle('Reference',fontName='Helvetica',fontSize=9.0,leading=12.3,textColor=INK,spaceAfter=6),
}

def inline(s):
 s=html.escape(s,quote=False)
 s=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',lambda m:f'<a href="{m[2]}" color="#24677e">{m[1]}</a>',s)
 s=re.sub(r'`([^`]+)`',r'<font name="Courier">\1</font>',s)
 s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s)
 s=re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)',r'<i>\1</i>',s)
 return s

def para(s,style='body'):
 return Paragraph(inline(s),styles[style])

class Doc(BaseDocTemplate):
 def __init__(self,path,label):
  self.label=label;self.bookmark_index=0
  super().__init__(str(path),pagesize=A4,leftMargin=MARGIN,rightMargin=MARGIN,topMargin=53,bottomMargin=45,
                   title=label,author='MAUVE internal working documents')
  frame=Frame(MARGIN,45,WIDTH,H-98,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
  self.addPageTemplates(PageTemplate(id='main',frames=frame,onPage=self.decorate))
 def decorate(self,c,doc):
  c.saveState();c.setFont('Helvetica',8.1);c.setFillColor(MUTED)
  c.drawString(MARGIN,H-30,self.label.upper());c.setStrokeColor(colors.HexColor('#d4e0e5'))
  c.line(MARGIN,H-37,W-MARGIN,H-37)
  c.drawString(MARGIN,25,'12 September 2026  |  MAUVE-JWST  |  Internal working document')
  c.drawRightString(W-MARGIN,25,str(doc.page));c.restoreState()
 def afterFlowable(self,f):
  if isinstance(f,Paragraph) and f.style.name in ['Heading2','Heading3']:
   self.bookmark_index+=1;key=f'section-{self.bookmark_index}';self.canv.bookmarkPage(key)
   self.canv.addOutlineEntry(f.getPlainText(),key,level=0 if f.style.name=='Heading2' else 1,closed=False)

def render(md):
 text=md.read_text();lines=text.splitlines();story=[];i=0;eqn=0;tablecount=0;figurecount=0
 eqpaths=DATA['equations'].get(md.stem,[])
 while i<len(lines):
  line=lines[i].strip()
  if not line:i+=1;continue
  if line=='$$':
   j=i+1
   while lines[j].strip()!='$$':j+=1
   path=Path(eqpaths[eqn]);eqn+=1;pw,ph=PILImage.open(path).size
   sw=min(WIDTH,pw/2.7);sh=sw*ph/pw
   story.extend([Spacer(1,3),Image(str(path),width=sw,height=sh,hAlign='CENTER'),Spacer(1,10)])
   i=j+1;continue
  if line.startswith('!['):
   m=re.match(r'!\[([^\]]*)\]\(([^)]+)\)',line);path=md.parent/m[2];pw,ph=PILImage.open(path).size
   sw=min(WIDTH,WIDTH*.99);sh=sw*ph/pw
   if sh>265:sw*=265/sh;sh=265
   contents=[Image(str(path),width=sw,height=sh,hAlign='CENTER'),Spacer(1,6)];figurecount+=1
   j=i+1
   while j<len(lines) and not lines[j].strip():j+=1
   if j<len(lines) and lines[j].startswith('*Figure '):contents.append(para(lines[j].strip().strip('*'),'caption'));i=j+1
   else:i+=1
   story.append(KeepTogether(contents));continue
  if line.startswith('|'):
   rows=[]
   while i<len(lines) and lines[i].strip().startswith('|'):
    row=[x.strip() for x in lines[i].strip().strip('|').split('|')]
    if not all(re.fullmatch(r':?-+:?',x) for x in row):rows.append(row)
    i+=1
   n=len(rows[0]);fractions={2:[.27,.73],3:[.27,.36,.37],4:[.28,.18,.18,.36]}.get(n,[1/n]*n)
   if rows[0][0]=='Conditional explanation':fractions=[.23,.41,.36]
   if rows[0][0]=='Evidence':fractions=[.28,.35,.37]
   tabledata=[[para(cell,'th' if ri==0 else 'cell') for cell in row] for ri,row in enumerate(rows)]
   table=LongTable(tabledata,colWidths=[WIDTH*f for f in fractions],repeatRows=1,hAlign='LEFT')
   table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),TEAL),('VALIGN',(0,0),(-1,-1),'TOP'),
      ('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),
      ('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),
      ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#eff3f5'),colors.white]),
      ('LINEBELOW',(0,-1),(-1,-1),.4,colors.HexColor('#d4e0e5'))]))
   story.extend([table,Spacer(1,10)]);tablecount+=1;continue
  if line.startswith('# '):story.append(para(line[2:],'title'));i+=1;continue
  if line.startswith('## '):story.append(para(line[3:],'h2'));i+=1;continue
  if line.startswith('### '):story.append(para(line[4:],'h3'));i+=1;continue
  if line.startswith('- '):story.append(para(line[2:],'ref' if line.startswith('- **[R') else 'body'));i+=1;continue
  paragraph=[line];i+=1
  while i<len(lines) and lines[i].strip() and not re.match(r'^(#|\||!\[|\$\$|- )',lines[i]):
   paragraph.append(lines[i].strip());i+=1
  story.append(para(' '.join(paragraph)))
 label='Goal 3: scientific justification' if 'justification' in md.stem else 'Cycle 6: revised Treasury proposal'
 out=md.with_suffix('.pdf');doc=Doc(out,label);doc.build(story)
 return {'source':str(md),'pdf':str(out),'words':len(text.split()),'equations':eqn,'figures':figurecount,'tables':tablecount,'bytes':out.stat().st_size}

if __name__=='__main__':
 results=[render(p) for p in sorted(ROOT.glob('20260912_MAUVE_JWST_*.md'))]
 (ASSETS/'build_record.json').write_text(json.dumps(results,indent=2))
 print(json.dumps(results,indent=2))
