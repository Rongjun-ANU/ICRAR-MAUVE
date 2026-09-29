"""Page numbers for the plain, author-style 29 September report."""
from io import BytesIO
from pathlib import Path
import argparse
import json
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from pypdf import PdfReader, PdfWriter
from pypdf.generic import NameObject, TextStringObject

def furniture(n,total,w,h):
    s=BytesIO();c=canvas.Canvas(s,pagesize=(w,h))
    c.setFillColor(HexColor('#444444'));c.setFont('Times-Roman',8)
    c.drawString(51,h-29,'MAUVE-MUSE / local gas and emission-line response')
    c.drawRightString(w-51,h-29,'29 September 2026')
    c.setStrokeColor(HexColor('#999999'));c.setLineWidth(.35);c.line(51,37,w-51,37)
    c.drawString(51,24,'Derivation and numerical example')
    c.drawRightString(w-51,24,f'{n} / {total}')
    c.save();s.seek(0);return PdfReader(s).pages[0]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('body');ap.add_argument('output');a=ap.parse_args()
    writer=PdfWriter();writer.append(PdfReader(a.body),import_outline=True)
    def norm(node):
        while node is not None:
            obj=node.get_object();title=str(obj.get('/Title',''));half=len(title)//2
            if title and len(title)%2==0 and title[:half]==title[half:]:
                obj[NameObject('/Title')]=TextStringObject(title[:half])
            if obj.get('/First') is not None:norm(obj['/First'])
            node=obj.get('/Next')
    norm(writer.get_outline_root().get('/First'))
    total=len(writer.pages)
    for i,page in enumerate(writer.pages,1):
        page.merge_page(furniture(i,total,float(page.mediabox.width),float(page.mediabox.height)))
    writer.add_metadata({'/Title':'The local response of gas, star formation, and emission lines during Virgo infall',
                         '/Author':'Research synthesis prepared for Rongjun Huang / MAUVE',
                         '/Subject':'Spatially resolved derivation and MAUVE-scale numerical example, 29 September 2026'})
    with open(a.output,'wb') as f:writer.write(f)
    print(json.dumps(dict(pages=total,output=a.output,bytes=Path(a.output).stat().st_size)))

if __name__=='__main__':main()
