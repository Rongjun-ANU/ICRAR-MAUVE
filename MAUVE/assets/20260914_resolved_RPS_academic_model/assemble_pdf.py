"""Academic title page and restrained page furniture for the offline report."""
from io import BytesIO
from pathlib import Path
import argparse
import json
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.colors import HexColor
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
from pypdf import PdfReader, PdfWriter
from pypdf.generic import NameObject, TextStringObject

def para(c,text,top,size=14,leading=21,bold=False):
    style=ParagraphStyle('cover',fontName='Times-Bold' if bold else 'Times-Roman',
                         fontSize=size,leading=leading,textColor=HexColor('#21343d'))
    p=Paragraph(text,style);_,height=p.wrap(A4[0]-112,800);p.drawOn(c,56,top-height)

def cover():
    s=BytesIO();c=canvas.Canvas(s,pagesize=A4);w,h=A4
    c.setFillColor(HexColor('#174e5a'));c.rect(56,h-91,w-112,2,fill=1,stroke=0)
    para(c,'MAUVE-MUSE / ICRAR',h-112,12,18,True)
    para(c,'Resolved gas loss,<br/>star formation, and<br/>continuing ionization<br/>in Virgo galaxies',h-162,29,37,True)
    para(c,'Observational constraints, analytical derivation,<br/>and a partial fit to MAUVE-MUSE',h-345,16,24)
    para(c,'Research report<br/>14 September 2026<br/>Prepared for Rongjun Huang',h-450,12,21)
    c.setStrokeColor(HexColor('#b8c9ce'));c.line(56,225,w-56,225)
    para(c,'26 galaxy products &nbsp; | &nbsp; 57 numbered equations<br/>7 research figures &nbsp; | &nbsp; 26 literature references',202,11,20)
    para(c,'Physical notation follows the local molecular-gas regulator.<br/>Observations, model assumptions, and fitted constraints are distinguished.',128,10,16)
    c.showPage();c.save();s.seek(0);return s

def furniture(n,total,w,h):
    s=BytesIO();c=canvas.Canvas(s,pagesize=(w,h))
    c.setFillColor(HexColor('#536b73'));c.setFont('Helvetica',7)
    c.drawString(51,h-29,'RESOLVED GAS LOSS, STAR FORMATION, AND IONIZATION')
    c.drawRightString(w-51,h-29,'14 SEP 2026')
    c.setStrokeColor(HexColor('#ccd8da'));c.setLineWidth(.4);c.line(51,37,w-51,37)
    c.drawString(51,24,'MAUVE-MUSE / ICRAR')
    c.drawRightString(w-51,24,f'{n} / {total}')
    c.save();s.seek(0);return PdfReader(s).pages[0]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('body');ap.add_argument('output');a=ap.parse_args()
    writer=PdfWriter();writer.append(PdfReader(cover()));writer.append(PdfReader(a.body),import_outline=True)
    def norm(node):
        while node is not None:
            obj=node.get_object();title=str(obj.get('/Title',''));half=len(title)//2
            if title and len(title)%2==0 and title[:half]==title[half:]:
                obj[NameObject('/Title')]=TextStringObject(title[:half])
            if obj.get('/First') is not None:norm(obj['/First'])
            node=obj.get('/Next')
    norm(writer.get_outline_root().get('/First'))
    total=len(writer.pages)
    for i,page in enumerate(writer.pages[1:],2):
        page.merge_page(furniture(i,total,float(page.mediabox.width),float(page.mediabox.height)))
    writer.add_metadata({'/Title':'Resolved gas loss, star formation, and continuing ionization in Virgo galaxies',
                         '/Author':'Research synthesis prepared for Rongjun Huang / MAUVE',
                         '/Subject':'Observations, analytical derivation and partial fit, 14 September 2026'})
    with open(a.output,'wb') as f:writer.write(f)
    print(json.dumps(dict(pages=total,output=a.output,bytes=Path(a.output).stat().st_size)))

if __name__=='__main__':main()
