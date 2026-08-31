"""Add a ReportLab cover and unobtrusive running furniture to the math-safe PDF."""
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


def text_box(c, text, x, top, width, size=12, leading=17, color="#273e47"):
    style = ParagraphStyle(
        "cover", fontName="AtlasSans", fontSize=size, leading=leading,
        textColor=HexColor(color), spaceAfter=0,
    )
    p = Paragraph(text, style)
    _, height = p.wrap(width, 1000)
    p.drawOn(c, x, top - height)
    return top - height


def make_cover():
    pdfmetrics.registerFont(TTFont("AtlasSans", "/Library/Fonts/PlusJakartaSans-Regular.ttf"))
    pdfmetrics.registerFont(TTFont("AtlasBold", "/Library/Fonts/PlusJakartaSans-Bold.ttf"))
    stream = BytesIO()
    c = canvas.Canvas(stream, pagesize=A4)
    w, h = A4
    c.setFillColor(HexColor("#f8faf9"))
    c.rect(0, 0, w, h, fill=1, stroke=0)
    c.setFillColor(HexColor("#176e72"))
    c.rect(50, h - 90, 42, 4, fill=1, stroke=0)
    c.setFont("AtlasBold", 10)
    c.drawString(50, h - 116, "MAUVE / ICRAR   -   RESEARCH REPORT")
    c.setFillColor(HexColor("#173f4b"))
    c.setFont("AtlasBold", 45)
    c.drawString(48, h - 190, "Nebular")
    c.drawString(48, h - 245, "Metallicity")
    text_box(c, "Atomic physics, calibration methods,<br/>and a MUSE observing guide",
             50, h - 275, w - 100, size=19, leading=27)
    text_box(c, "From abundances and ionizing photons to level populations, line fluxes, "
             "and the inverse problem of measuring gas-phase oxygen.",
             50, h - 370, w - 105, size=12, leading=19)
    c.setStrokeColor(HexColor("#becfd1"))
    c.setLineWidth(0.65)
    c.line(50, 330, w - 50, 330)
    windows = [("7000 A", "4800-7000"), ("8900 A", "4800-8900"), ("9300 A", "4800-9300")]
    for j, (label, interval) in enumerate(windows):
        x = 50 + j * 166
        c.setFillColor(HexColor("#176e72"))
        c.setFont("AtlasBold", 18)
        c.drawString(x, 295, label)
        c.setFillColor(HexColor("#4e6167"))
        c.setFont("AtlasSans", 9)
        c.drawString(x, 275, interval + " Angstrom")
    text_box(c, "Direct and semi-direct abundances - strong-line calibrations - photoionization "
             "inference - density and temperature diagnostics - ionization sensitivity - "
             "other elements and wavelengths.",
             50, 234, w - 105, size=10.5, leading=17)
    c.setFont("AtlasBold", 10)
    c.setFillColor(HexColor("#173f4b"))
    c.drawString(50, 118, "30 AUGUST 2026")
    c.setFont("AtlasSans", 9)
    c.drawString(50, 98, "Prepared for Rongjun | Oxygen abundance by default")
    c.setFillColor(HexColor("#62747a"))
    c.setFont("AtlasSans", 8)
    c.drawString(50, 59, "Literature and physics synthesis. No MAUVE science products were modified.")
    c.showPage()
    c.save()
    stream.seek(0)
    return stream


def furniture(page_number, page_count, width, height):
    stream = BytesIO()
    c = canvas.Canvas(stream, pagesize=(width, height))
    c.setFillColor(HexColor("#52686e"))
    c.setFont("Helvetica", 7.2)
    c.drawString(51, height - 29, "NEBULAR METALLICITY")
    c.drawRightString(width - 51, height - 29, "MAUVE / ICRAR  |  30 AUG 2026")
    c.setStrokeColor(HexColor("#d2dddf"))
    c.setLineWidth(0.4)
    c.line(51, 37, width - 51, 37)
    c.drawString(51, 24, "Atomic physics - calibration atlas - MUSE")
    c.drawRightString(width - 51, 24, f"{page_number} / {page_count}")
    c.save()
    stream.seek(0)
    return PdfReader(stream).pages[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("body_pdf", type=Path)
    parser.add_argument("output_pdf", type=Path)
    args = parser.parse_args()
    writer = PdfWriter()
    writer.append(PdfReader(make_cover()))
    writer.append(PdfReader(args.body_pdf), import_outline=True)
    # Skia may repeat a heading in an outline title when it crosses layout fragments.
    # Fix only exact doubled strings, leaving all destinations unchanged.
    def normalize_outline(first):
        node = first
        while node is not None:
            obj = node.get_object()
            title = str(obj.get("/Title", ""))
            mid = len(title) // 2
            if title and len(title) % 2 == 0 and title[:mid] == title[mid:]:
                obj[NameObject("/Title")] = TextStringObject(title[:mid])
            if obj.get("/First") is not None:
                normalize_outline(obj["/First"])
            node = obj.get("/Next")
    normalize_outline(writer.get_outline_root().get("/First"))
    total = len(writer.pages)
    for index, page in enumerate(writer.pages[1:], start=2):
        page.merge_page(furniture(index, total, float(page.mediabox.width), float(page.mediabox.height)))
    writer.add_metadata({
        "/Title": "Nebular Metallicity: Atomic Physics, Calibration Atlas and MUSE Guide",
        "/Author": "Research report prepared for Rongjun / MAUVE",
        "/Subject": "Gas-phase abundances; atomic physics; MUSE 7000/8900/9300 Angstrom",
        "/Keywords": "O/H, metallicity, direct method, strong lines, photoionization, density, temperature",
    })
    with args.output_pdf.open("wb") as stream:
        writer.write(stream)
    check = PdfReader(args.output_pdf)
    assert len(check.pages) == total
    assert all(len(p.extract_text().strip()) > 50 for p in check.pages)
    print(json.dumps({"output": str(args.output_pdf), "pages": total,
                      "bytes": args.output_pdf.stat().st_size}, indent=2))


if __name__ == "__main__":
    main()
