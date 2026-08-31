"""Add a ReportLab cover and page furniture to the verified MathML body."""
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


def text_box(c, text, x, top, width, size=12, leading=18, color="#273e47"):
    style = ParagraphStyle("cover",fontName="AtlasSans",fontSize=size,leading=leading,
                           textColor=HexColor(color),spaceAfter=0)
    p = Paragraph(text,style)
    _, height = p.wrap(width,1000)
    p.drawOn(c,x,top-height)
    return top-height


def make_cover():
    pdfmetrics.registerFont(TTFont("AtlasSans","/Library/Fonts/PlusJakartaSans-Regular.ttf"))
    pdfmetrics.registerFont(TTFont("AtlasBold","/Library/Fonts/PlusJakartaSans-Bold.ttf"))
    stream = BytesIO()
    c = canvas.Canvas(stream,pagesize=A4)
    w,h = A4
    c.setFillColor(HexColor("#f8faf9")); c.rect(0,0,w,h,fill=1,stroke=0)
    c.setFillColor(HexColor("#176e72")); c.rect(50,h-90,42,4,fill=1,stroke=0)
    c.setFont("AtlasBold",10); c.drawString(50,h-116,"MAUVE / ICRAR   -   DEEP RESEARCH")
    c.setFillColor(HexColor("#173f4b")); c.setFont("AtlasBold",44)
    c.drawString(48,h-190,"Molecular")
    c.drawString(48,h-245,"Gas Supply")
    text_box(c,"Definitions, chemical rates,<br/>and observational inference",
             50,h-275,w-100,size=19,leading=27)
    text_box(c,"From the resolved bathtub model in Paper I to a practical programme "
               "for distinguishing molecular replenishment, H I-to-H2 conversion "
               "and cloud-boundary accretion.",50,h-363,w-105,size=12,leading=19)
    c.setStrokeColor(HexColor("#becfd1")); c.setLineWidth(.65)
    c.line(50,325,w-50,325)
    c.setFillColor(HexColor("#176e72")); c.setFont("AtlasBold",11)
    c.drawString(50,294,"MUSE  |  ALMA  |  MeerKAT  |  UVIT  |  HST")
    text_box(c,"Step-by-step conservation laws and estimators<br/>"
               "Resolution-aware inference and H I alternatives<br/>"
               "A ranked programme of additional observations",
             50,269,w-105,size=11,leading=21)
    c.setFillColor(HexColor("#173f4b")); c.setFont("AtlasBold",10)
    c.drawString(50,118,"31 AUGUST 2026")
    c.setFont("AtlasSans",9); c.drawString(50,98,"Prepared for Rongjun | MAUVE observational framework")
    c.setFillColor(HexColor("#62747a")); c.setFont("AtlasSans",8)
    c.drawString(50,59,"Literature, derivations and illustrative calculations; no measured MAUVE rates.")
    c.showPage(); c.save(); stream.seek(0)
    return stream


def furniture(number,total,width,height):
    stream = BytesIO(); c = canvas.Canvas(stream,pagesize=(width,height))
    c.setFillColor(HexColor("#52686e")); c.setFont("Helvetica",7.2)
    c.drawString(51,height-29,"MOLECULAR GAS SUPPLY")
    c.drawRightString(width-51,height-29,"MAUVE / ICRAR  |  31 AUG 2026")
    c.setStrokeColor(HexColor("#d2dddf")); c.setLineWidth(.4); c.line(51,37,width-51,37)
    c.drawString(51,24,"Reservoirs - chemistry - transport - observations")
    c.drawRightString(width-51,24,f"{number} / {total}")
    c.save(); stream.seek(0)
    return PdfReader(stream).pages[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("body_pdf",type=Path)
    parser.add_argument("output_pdf",type=Path)
    args = parser.parse_args()
    writer = PdfWriter()
    writer.append(PdfReader(make_cover()))
    writer.append(PdfReader(args.body_pdf),import_outline=True)
    def normalize_outline(first):
        node = first
        while node is not None:
            obj = node.get_object(); title = str(obj.get("/Title","")); mid=len(title)//2
            if title and len(title)%2==0 and title[:mid]==title[mid:]:
                obj[NameObject("/Title")] = TextStringObject(title[:mid])
            if obj.get("/First") is not None: normalize_outline(obj["/First"])
            node=obj.get("/Next")
    normalize_outline(writer.get_outline_root().get("/First"))
    total=len(writer.pages)
    for i,page in enumerate(writer.pages[1:],start=2):
        page.merge_page(furniture(i,total,float(page.mediabox.width),float(page.mediabox.height)))
    writer.add_metadata({"/Title":"Molecular Gas Supply: Definitions and Observational Inference",
        "/Author":"Research synthesis prepared for Rongjun / MAUVE",
        "/Subject":"Resolved molecular supply, HI-H2 conversion, GMC accretion and MAUVE observing strategy",
        "/Keywords":"MAUVE, MUSE, ALMA, MeerKAT, UVIT, HST, H2 formation, photodissociation, molecular supply"})
    with args.output_pdf.open("wb") as stream: writer.write(stream)
    result=PdfReader(args.output_pdf)
    assert len(result.pages)==total
    assert all(len(page.extract_text().strip())>50 for page in result.pages)
    print(json.dumps({"output":str(args.output_pdf),"pages":total,"bytes":args.output_pdf.stat().st_size},indent=2))


if __name__=="__main__": main()
