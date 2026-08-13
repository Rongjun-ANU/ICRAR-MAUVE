#!/usr/bin/env python3
"""Render the MAUVE analytical-environment Markdown review as an archival PDF."""

from __future__ import annotations

import html
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    NextPageTemplate,
    PageBreak,
    PageTemplate,
    Paragraph,
    Preformatted,
    Spacer,
    Table,
    TableStyle,
)
from reportlab.platypus.tableofcontents import TableOfContents


ASSET_DIR = Path(__file__).resolve().parent
REPORT_DIR = ASSET_DIR.parent.parent
SOURCE = REPORT_DIR / (
    "20260813 Analytical Models of Galaxy Evolution in Cluster Environments "
    "and Ram Pressure Stripping.md"
)
OUTPUT = SOURCE.with_suffix(".pdf")

PAGE_W, PAGE_H = A4
MARGIN_L = 19 * mm
MARGIN_R = 17 * mm
MARGIN_T = 18 * mm
MARGIN_B = 18 * mm
CONTENT_W = PAGE_W - MARGIN_L - MARGIN_R
CONTENT_H = PAGE_H - MARGIN_T - MARGIN_B

NAVY = colors.HexColor("#16324F")
BLUE = colors.HexColor("#2F6690")
GREEN = colors.HexColor("#4A7C59")
RED = colors.HexColor("#A23E48")
GOLD = colors.HexColor("#D5A021")
INK = colors.HexColor("#20252A")
MUTED = colors.HexColor("#5C646B")
PAPER = colors.HexColor("#F7F4EE")
GRID = colors.HexColor("#D8D4CB")


def register_fonts() -> None:
    root = Path("/System/Library/Fonts/Supplemental")
    pdfmetrics.registerFont(TTFont("MAUVEArial", root / "Arial.ttf"))
    pdfmetrics.registerFont(TTFont("MAUVEArial-Bold", root / "Arial Bold.ttf"))
    pdfmetrics.registerFont(TTFont("MAUVEArial-Italic", root / "Arial Italic.ttf"))
    pdfmetrics.registerFont(TTFont("MAUVEArial-BoldItalic", root / "Arial Bold Italic.ttf"))
    pdfmetrics.registerFontFamily(
        "MAUVEArial",
        normal="MAUVEArial",
        bold="MAUVEArial-Bold",
        italic="MAUVEArial-Italic",
        boldItalic="MAUVEArial-BoldItalic",
    )


def normalize(text: str) -> str:
    replacements = {
        "\u2010": "-", "\u2011": "-", "\u2012": "-", "\u2013": "-", "\u2014": "-",
        "\u2212": "-", "\u00a0": " ", "\u2009": " ", "\u202f": " ",
        "\u2018": "'", "\u2019": "'", "\u201c": '"', "\u201d": '"',
        "\u2265": ">=", "\u2264": "<=", "\u2248": "approximately ", "\u2273": ">=",
        "\u00d7": "x", "\u03b1": "alpha", "\u03b2": "beta", "\u03c1": "rho",
        "\u03c3": "sigma", "\u03c4": "tau", "\u03bb": "lambda", "\u0394": "Delta",
        "\u03a3": "Sigma", "\u03a6": "Phi", "\u03c0": "pi", "\u221a": "sqrt",
        "\u2299": "sun", "\u2192": "->", "\u201e": '"', "\u2026": "...",
        "\u2003": " ", "\u2002": " ", "\ufeff": "",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def inline_markup(text: str) -> str:
    text = normalize(text)
    placeholders: list[tuple[str, str, str]] = []

    def hold_link(match: re.Match[str]) -> str:
        key = f"@@LINK{len(placeholders)}@@"
        placeholders.append((key, match.group(1), match.group(2)))
        return key

    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", hold_link, text)
    text = html.escape(text, quote=False)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", text)
    text = re.sub(r"`([^`]+)`", r"<font name='Courier'>\1</font>", text)
    for key, label, target in placeholders:
        label_markup = html.escape(normalize(label), quote=False)
        target_clean = html.escape(normalize(target), quote=True)
        replacement = f'<link href="{target_clean}" color="#2F6690">{label_markup}</link>'
        text = text.replace(key, replacement)
    return text


def make_styles() -> dict[str, ParagraphStyle]:
    sample = getSampleStyleSheet()
    return {
        "body": ParagraphStyle(
            "Body", parent=sample["BodyText"], fontName="MAUVEArial", fontSize=9.25,
            leading=12.7, textColor=INK, spaceAfter=4.5 * mm, allowWidows=0, allowOrphans=0,
        ),
        "small": ParagraphStyle(
            "Small", parent=sample["BodyText"], fontName="MAUVEArial", fontSize=7.6,
            leading=9.6, textColor=INK, spaceAfter=2 * mm,
        ),
        "caption": ParagraphStyle(
            "Caption", parent=sample["BodyText"], fontName="MAUVEArial-Italic", fontSize=8.1,
            leading=10.5, textColor=MUTED, alignment=TA_CENTER, spaceBefore=1.4 * mm,
            spaceAfter=4.2 * mm,
        ),
        "bullet": ParagraphStyle(
            "Bullet", parent=sample["BodyText"], fontName="MAUVEArial", fontSize=9.1,
            leading=12.3, leftIndent=5.5 * mm, firstLineIndent=-3.5 * mm, textColor=INK,
            spaceAfter=1.7 * mm,
        ),
        "quote": ParagraphStyle(
            "Quote", parent=sample["BodyText"], fontName="MAUVEArial-Italic", fontSize=8.9,
            leading=12.2, leftIndent=7 * mm, rightIndent=4 * mm, borderColor=BLUE,
            borderWidth=1.5, borderPadding=3 * mm, backColor=colors.HexColor("#EEF3F6"),
            textColor=INK, spaceAfter=4 * mm,
        ),
        "code": ParagraphStyle(
            "Code", parent=sample["Code"], fontName="Courier", fontSize=7.6, leading=10.2,
            leftIndent=3 * mm, rightIndent=3 * mm, borderColor=GRID, borderWidth=0.5,
            borderPadding=3 * mm, backColor=colors.HexColor("#F1F2F2"), textColor=INK,
            spaceBefore=1 * mm, spaceAfter=4 * mm,
        ),
        "h1": ParagraphStyle(
            "H1", parent=sample["Heading1"], fontName="MAUVEArial-Bold", fontSize=17,
            leading=21, textColor=NAVY, spaceBefore=8 * mm, spaceAfter=4 * mm,
            keepWithNext=True,
        ),
        "h2": ParagraphStyle(
            "H2", parent=sample["Heading2"], fontName="MAUVEArial-Bold", fontSize=13.2,
            leading=16.2, textColor=BLUE, spaceBefore=5.5 * mm, spaceAfter=3 * mm,
            keepWithNext=True,
        ),
        "h3": ParagraphStyle(
            "H3", parent=sample["Heading3"], fontName="MAUVEArial-Bold", fontSize=10.8,
            leading=13.5, textColor=GREEN, spaceBefore=4.2 * mm, spaceAfter=2.2 * mm,
            keepWithNext=True,
        ),
        "h4": ParagraphStyle(
            "H4", parent=sample["Heading4"], fontName="MAUVEArial-Bold", fontSize=9.6,
            leading=12, textColor=INK, spaceBefore=3 * mm, spaceAfter=1.8 * mm,
            keepWithNext=True,
        ),
        "title": ParagraphStyle(
            "Title", parent=sample["Title"], fontName="MAUVEArial-Bold", fontSize=26,
            leading=31, textColor=NAVY, alignment=TA_LEFT, spaceAfter=7 * mm,
        ),
        "subtitle": ParagraphStyle(
            "Subtitle", parent=sample["BodyText"], fontName="MAUVEArial", fontSize=13,
            leading=17, textColor=BLUE, spaceAfter=6 * mm,
        ),
        "meta": ParagraphStyle(
            "Meta", parent=sample["BodyText"], fontName="MAUVEArial", fontSize=9.6,
            leading=13, textColor=MUTED, spaceAfter=2 * mm,
        ),
        "toc_h": ParagraphStyle(
            "TOCHeading", parent=sample["Heading1"], fontName="MAUVEArial-Bold", fontSize=18,
            leading=22, textColor=NAVY, spaceAfter=5 * mm,
        ),
        "toc1": ParagraphStyle(
            "TOC1", fontName="MAUVEArial", fontSize=9.2, leading=12.2,
            leftIndent=0, firstLineIndent=0, textColor=INK, spaceBefore=1.4 * mm,
        ),
        "toc2": ParagraphStyle(
            "TOC2", fontName="MAUVEArial", fontSize=8.5, leading=11.2,
            leftIndent=5 * mm, firstLineIndent=0, textColor=MUTED, spaceBefore=0.8 * mm,
        ),
    }


class ReviewDocTemplate(BaseDocTemplate):
    def __init__(self, filename: str, styles: dict[str, ParagraphStyle]):
        super().__init__(
            filename, pagesize=A4, leftMargin=MARGIN_L, rightMargin=MARGIN_R,
            topMargin=MARGIN_T, bottomMargin=MARGIN_B,
            title="Analytical Models of Galaxy Evolution in Cluster Environments and Ram Pressure Stripping",
            author="MAUVE literature review",
            subject="Analytic, semi-analytic, simulation-calibrated, and individual models of environmental galaxy evolution",
        )
        frame = Frame(MARGIN_L, MARGIN_B, CONTENT_W, CONTENT_H, id="body")
        self.addPageTemplates([
            PageTemplate(id="Title", frames=frame, onPage=self.draw_title_page),
            PageTemplate(id="Body", frames=frame, onPage=self.draw_body_page),
        ])
        self.styles = styles
        self._heading_count = 0

    def beforeDocument(self) -> None:
        self._heading_count = 0

    def draw_title_page(self, canvas, doc) -> None:
        canvas.saveState()
        canvas.setFillColor(NAVY)
        canvas.rect(0, PAGE_H - 16 * mm, PAGE_W, 16 * mm, fill=1, stroke=0)
        canvas.setFillColor(GOLD)
        canvas.rect(0, PAGE_H - 17.5 * mm, PAGE_W, 1.5 * mm, fill=1, stroke=0)
        canvas.restoreState()

    def draw_body_page(self, canvas, doc) -> None:
        canvas.saveState()
        canvas.setStrokeColor(GRID)
        canvas.setLineWidth(0.45)
        canvas.line(MARGIN_L, PAGE_H - 12.5 * mm, PAGE_W - MARGIN_R, PAGE_H - 12.5 * mm)
        canvas.setFont("MAUVEArial", 7.2)
        canvas.setFillColor(MUTED)
        canvas.drawString(MARGIN_L, PAGE_H - 9.5 * mm, "MAUVE | Analytical cluster-environment models")
        canvas.drawRightString(PAGE_W - MARGIN_R, 9 * mm, f"{doc.page}")
        canvas.setFont("MAUVEArial", 6.7)
        canvas.drawString(MARGIN_L, 9 * mm, "Search and verification date: 2026-08-13")
        canvas.restoreState()

    def afterFlowable(self, flowable) -> None:
        if not isinstance(flowable, Paragraph):
            return
        style_name = flowable.style.name
        level_map = {"H1": 0, "H2": 0, "H3": 1}
        if style_name not in level_map:
            return
        text = flowable.getPlainText()
        if text == "Contents":
            return
        self._heading_count += 1
        anchor = f"heading-{self._heading_count}"
        self.canv.bookmarkPage(anchor)
        self.canv.addOutlineEntry(text, anchor, level=level_map[style_name], closed=False)
        self.notify("TOCEntry", (level_map[style_name], text, self.page, anchor))


def image_flowable(path: Path, max_w: float = CONTENT_W, max_h: float = 160 * mm) -> Image:
    img = Image(str(path))
    scale = min(max_w / img.imageWidth, max_h / img.imageHeight)
    img.drawWidth = img.imageWidth * scale
    img.drawHeight = img.imageHeight * scale
    img.hAlign = "CENTER"
    return img


def table_flowable(rows: list[list[str]], styles: dict[str, ParagraphStyle]) -> Table:
    ncols = max(len(row) for row in rows)
    normalized_rows = [row + [""] * (ncols - len(row)) for row in rows]
    if ncols == 4:
        widths = [0.19 * CONTENT_W, 0.28 * CONTENT_W, 0.29 * CONTENT_W, 0.24 * CONTENT_W]
    elif ncols == 3:
        widths = [0.25 * CONTENT_W, 0.39 * CONTENT_W, 0.36 * CONTENT_W]
    else:
        widths = [CONTENT_W / ncols] * ncols
    data = []
    for ridx, row in enumerate(normalized_rows):
        row_style = ParagraphStyle(
            f"TableCell{ridx}", parent=styles["small"], fontName=("MAUVEArial-Bold" if ridx == 0 else "MAUVEArial"),
            fontSize=7.0 if ncols >= 4 else 7.5, leading=8.7 if ncols >= 4 else 9.3,
            textColor=(colors.white if ridx == 0 else INK), spaceAfter=0,
        )
        data.append([Paragraph(inline_markup(cell), row_style) for cell in row])
    table = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT", splitByRow=1)
    commands = [
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("GRID", (0, 0), (-1, -1), 0.35, GRID),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 3.2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3.2),
        ("TOPPADDING", (0, 0), (-1, -1), 3.4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3.4),
    ]
    for ridx in range(1, len(data)):
        if ridx % 2 == 0:
            commands.append(("BACKGROUND", (0, ridx), (-1, ridx), colors.HexColor("#F2F4F4")))
    table.setStyle(TableStyle(commands))
    return table


def parse_markdown(source: str, styles: dict[str, ParagraphStyle]) -> list:
    lines = normalize(source).splitlines()
    flowables = []
    paragraph: list[str] = []
    code: list[str] = []
    table_lines: list[str] = []
    in_code = False
    seen_title = False
    appendix_started = False

    def flush_paragraph() -> None:
        nonlocal paragraph
        if paragraph:
            text = " ".join(item.strip() for item in paragraph).strip()
            if text:
                flowables.append(Paragraph(inline_markup(text), styles["body"]))
            paragraph = []

    def flush_table() -> None:
        nonlocal table_lines
        if not table_lines:
            return
        parsed = []
        for line in table_lines:
            cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
            if cells and all(re.fullmatch(r":?-{3,}:?", cell) for cell in cells):
                continue
            parsed.append(cells)
        if parsed:
            flowables.extend([table_flowable(parsed, styles), Spacer(1, 4 * mm)])
        table_lines = []

    for index, line in enumerate(lines):
        if line.strip().startswith("```"):
            flush_paragraph()
            flush_table()
            if in_code:
                flowables.append(Preformatted("\n".join(code), styles["code"]))
                code = []
                in_code = False
            else:
                in_code = True
            continue
        if in_code:
            code.append(line)
            continue

        if line.startswith("|") and line.rstrip().endswith("|"):
            flush_paragraph()
            table_lines.append(line)
            continue
        if table_lines:
            flush_table()

        image_match = re.fullmatch(r"!\[([^]]*)\]\(([^)]+)\)", line.strip())
        if image_match:
            flush_paragraph()
            image_path = (REPORT_DIR / image_match.group(2)).resolve()
            flowables.extend([Spacer(1, 2 * mm), image_flowable(image_path), Spacer(1, 1 * mm)])
            continue

        heading = re.match(r"^(#{1,4})\s+(.+)$", line)
        if heading:
            flush_paragraph()
            level = len(heading.group(1))
            title = heading.group(2).strip()
            if level == 1 and not seen_title:
                seen_title = True
                continue
            if level == 1 and title.startswith("Appendix A") and not appendix_started:
                appendix_started = True
                flowables.append(PageBreak())
            elif level == 1:
                flowables.append(PageBreak())
            flowables.append(Paragraph(inline_markup(title), styles[f"h{level}"]))
            continue

        if not line.strip():
            flush_paragraph()
            continue

        list_match = re.match(r"^\s*(?:[-*]|(\d+)\.)\s+(.+)$", line)
        if list_match:
            flush_paragraph()
            marker = f"{list_match.group(1)}." if list_match.group(1) else "-"
            flowables.append(Paragraph(f"{marker}&nbsp;&nbsp;{inline_markup(list_match.group(2))}", styles["bullet"]))
            continue

        if line.startswith("> "):
            flush_paragraph()
            flowables.append(Paragraph(inline_markup(line[2:]), styles["quote"]))
            continue

        if line.startswith("*") and line.endswith("*") and not line.startswith("**"):
            flush_paragraph()
            flowables.append(Paragraph(inline_markup(line), styles["caption"]))
            continue

        paragraph.append(line)

    flush_paragraph()
    flush_table()
    return flowables


def title_story(styles: dict[str, ParagraphStyle]) -> list:
    hero = ASSET_DIR / "figure_02_model_stack.png"
    return [
        Spacer(1, 13 * mm),
        Paragraph(
            "Analytical Models of Galaxy Evolution in Cluster Environments and Ram Pressure Stripping",
            styles["title"],
        ),
        Paragraph("A systematic-style scoping review with a MAUVE implementation roadmap", styles["subtitle"]),
        Spacer(1, 2 * mm),
        Paragraph("Search and verification date: 2026-08-13 (Australia/Perth)", styles["meta"]),
        Paragraph("Evidence base: 80 DOI-verified papers plus one DOI-less foundational ICM model", styles["meta"]),
        Paragraph("Primary application: resolved nearby-cluster galaxies, especially Virgo and MAUVE", styles["meta"]),
        Spacer(1, 8 * mm),
        image_flowable(hero, max_w=CONTENT_W, max_h=92 * mm),
        Spacer(1, 5 * mm),
        Paragraph(
            "Core result: analytical models constrain susceptibility, pressure histories, and statistical clocks; "
            "a defensible history for one galaxy requires a probabilistic multi-observable forward model.",
            styles["quote"],
        ),
        NextPageTemplate("Body"),
        PageBreak(),
    ]


def toc_story(styles: dict[str, ParagraphStyle]) -> list:
    toc = TableOfContents()
    toc.levelStyles = [styles["toc1"], styles["toc2"], styles["toc2"]]
    return [Paragraph("Contents", styles["toc_h"]), toc, PageBreak()]


def main() -> None:
    register_fonts()
    styles = make_styles()
    source = SOURCE.read_text(encoding="utf-8")
    body = re.sub(r"^# .+?\n", "", source, count=1)
    story = title_story(styles) + toc_story(styles) + parse_markdown(body, styles)
    doc = ReviewDocTemplate(str(OUTPUT), styles)
    doc.multiBuild(story)
    print(f"created={OUTPUT}")


if __name__ == "__main__":
    main()
