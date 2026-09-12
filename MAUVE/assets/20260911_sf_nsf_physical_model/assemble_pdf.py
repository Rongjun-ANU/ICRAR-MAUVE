"""Add the dated report cover and page furniture to the offline MathML PDF."""
from io import BytesIO
from pathlib import Path
import argparse
import json
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
from pypdf import PdfReader, PdfWriter
from pypdf.generic import NameObject, TextStringObject

def paragraph(c, text, top, size=12, leading=19, color='#273e47'):
    p=Paragraph(text, ParagraphStyle('cover',fontName='AtlasSans',fontSize=size,
                leading=leading,textColor=HexColor(color)))
    _,h=p.wrap(A4[0]-105,1000)
    p.drawOn(c,50,top-h)

def cover():
    pdfmetrics.registerFont(TTFont('AtlasSans','/Library/Fonts/PlusJakartaSans-Regular.ttf'))
    pdfmetrics.registerFont(TTFont('AtlasBold','/Library/Fonts/PlusJakartaSans-Bold.ttf'))
    s=BytesIO(); c=canvas.Canvas(s,pagesize=A4);w,h=A4
    c.setFillColor(HexColor('#f8faf9'));c.rect(0,0,w,h,fill=1,stroke=0)
    c.setFillColor(HexColor('#176e72'));c.rect(50,h-90,42,4,fill=1,stroke=0)
    c.setFont('AtlasBold',10);c.drawString(50,h-116,'MAUVE / ICRAR   -   DEEP RESEARCH')
    c.setFillColor(HexColor('#173f4b'));c.setFont('AtlasBold',35)
    for y,t in [(185,'Remaining SF'),(234,'and Rising NSF')]:c.drawString(48,h-y,t)
    paragraph(c,'An observational audit and<br/>an analytical physical scenario',h-267,19,28)
    paragraph(c,'Outside-in gas loss, fading young stars, and continued ionization '
                'of the gas that remains. A framework linking occupancy, intensity '
                'and excitation to the actual MAUVE selection.',h-360,12,20)
    c.setStrokeColor(HexColor('#becfd1'));c.setLineWidth(.65);c.line(50,323,w-50,323)
    c.setFillColor(HexColor('#176e72'));c.setFont('AtlasBold',11)
    c.drawString(50,294,'26 SYSTEMS  |  40 EQUATIONS  |  6 FIGURES')
    paragraph(c,'Live notebook and FITS-derived summaries<br/>26 literature references, including contrary evidence<br/>'
                'Closed-form predictions and a practical fitting programme',269,11,21)
    c.setFillColor(HexColor('#173f4b'));c.setFont('AtlasBold',10)
    c.drawString(50,118,'11 SEPTEMBER 2026')
    c.setFont('AtlasSans',9);c.drawString(50,98,'Prepared for Rongjun | MAUVE observational framework')
    c.setFillColor(HexColor('#62747a'));c.setFont('AtlasSans',8)
    c.drawString(50,59,'An explanatory model with numerical checks; physical parameters have not been fitted.')
    c.showPage();c.save();s.seek(0);return s

def furniture(n,total,w,h):
    s=BytesIO();c=canvas.Canvas(s,pagesize=(w,h));c.setFillColor(HexColor('#52686e'))
    c.setFont('Helvetica',7.2);c.drawString(51,h-29,'REMAINING SF AND RISING NSF')
    c.drawRightString(w-51,h-29,'MAUVE / ICRAR  |  11 SEP 2026')
    c.setStrokeColor(HexColor('#d2dddf'));c.setLineWidth(.4);c.line(51,37,w-51,37)
    c.drawString(51,24,'Observations - gas reservoirs - ionization - falsifiable predictions')
    c.drawRightString(w-51,24,f'{n} / {total}');c.save();s.seek(0);return PdfReader(s).pages[0]

def main():
    p=argparse.ArgumentParser();p.add_argument('body_pdf',type=Path);p.add_argument('output_pdf',type=Path);a=p.parse_args()
    w=PdfWriter();w.append(PdfReader(cover()));w.append(PdfReader(a.body_pdf),import_outline=True)
    def normalize(node):
        while node is not None:
            obj=node.get_object();title=str(obj.get('/Title',''));mid=len(title)//2
            if title and len(title)%2==0 and title[:mid]==title[mid:]:obj[NameObject('/Title')]=TextStringObject(title[:mid])
            if obj.get('/First') is not None:normalize(obj['/First'])
            node=obj.get('/Next')
    normalize(w.get_outline_root().get('/First'))
    total=len(w.pages)
    for i,page in enumerate(w.pages[1:],2):page.merge_page(furniture(i,total,float(page.mediabox.width),float(page.mediabox.height)))
    w.add_metadata({'/Title':'Remaining Star Formation and Rising NSF Occupancy',
      '/Author':'Research synthesis prepared for Rongjun / MAUVE',
      '/Subject':'Observational audit and analytical scenario, 11 September 2026',
      '/Keywords':'MAUVE, ram pressure, star formation, SF, NSF, ND, DIG, gas regulation'})
    with a.output_pdf.open('wb') as f:w.write(f)
    r=PdfReader(a.output_pdf);assert len(r.pages)==total
    assert all(len(p.extract_text().strip())>50 for p in r.pages)
    print(json.dumps({'output':str(a.output_pdf),'pages':total,'bytes':a.output_pdf.stat().st_size}))

if __name__=='__main__':main()
