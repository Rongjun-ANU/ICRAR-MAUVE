"""Render the companion Markdown report without changing scientific content."""
from pathlib import Path
import html
import re
from urllib.parse import unquote

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, CondPageBreak, Frame, Image, KeepTogether,
    NextPageTemplate, PageBreak, PageTemplate, Paragraph, Spacer, Table, TableStyle,
)
from reportlab.platypus.tableofcontents import TableOfContents

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / '20261011_MAUVE_Paper_II_Findings_Synthesis_and_Paper_Plan.md'
TARGET = SOURCE.with_suffix('.pdf')
FONT = Path('/System/Library/Fonts/Supplemental')
for name, filename in [('Arial', 'Arial.ttf'), ('Arial-Bold', 'Arial Bold.ttf'),
                       ('Arial-Italic', 'Arial Italic.ttf'), ('Arial-BoldItalic', 'Arial Bold Italic.ttf'),
                       ('Georgia', 'Georgia.ttf'), ('Georgia-Bold', 'Georgia Bold.ttf')]:
    pdfmetrics.registerFont(TTFont(name, str(FONT / filename)))
pdfmetrics.registerFontFamily('Arial', normal='Arial', bold='Arial-Bold',
                             italic='Arial-Italic', boldItalic='Arial-BoldItalic')
pdfmetrics.registerFontFamily('Georgia', normal='Georgia', bold='Georgia-Bold')

INK = colors.HexColor('#1C2933')
BLUE = colors.HexColor('#245873')
MUTED = colors.HexColor('#526574')
PALE = colors.HexColor('#F0F5F7')
PAGE_W, PAGE_H = A4
MARGIN = 47
WIDTH = PAGE_W - 2 * MARGIN
styles = {
    'body': ParagraphStyle('body', fontName='Arial', fontSize=9.7, leading=14.1,
                           textColor=INK, spaceAfter=8, splitLongWords=True,
                           allowWidows=0, allowOrphans=0),
    'h2': ParagraphStyle('h2', fontName='Georgia-Bold', fontSize=16, leading=20,
                         textColor=BLUE, spaceBefore=16, spaceAfter=10, keepWithNext=True),
    'h3': ParagraphStyle('h3', fontName='Arial-Bold', fontSize=11.7, leading=15.2,
                         textColor=INK, spaceBefore=13, spaceAfter=7, keepWithNext=True),
    'title': ParagraphStyle('title', fontName='Georgia-Bold', fontSize=28, leading=34,
                            textColor=BLUE, spaceAfter=25),
    'kicker': ParagraphStyle('kicker', fontName='Arial-Bold', fontSize=10, leading=14,
                             textColor=MUTED, spaceAfter=18),
    'caption': ParagraphStyle('caption', fontName='Arial', fontSize=8.5, leading=12,
                              textColor=MUTED, spaceBefore=6, spaceAfter=13),
    'table': ParagraphStyle('table', fontName='Arial', fontSize=8.6, leading=11.8,
                            textColor=INK, splitLongWords=True),
    'tablehead': ParagraphStyle('tablehead', fontName='Arial-Bold', fontSize=8.6, leading=11.8,
                                textColor=colors.white, splitLongWords=True),
    'quote': ParagraphStyle('quote', fontName='Arial-Bold', fontSize=10, leading=14,
                            textColor=BLUE, backColor=PALE, borderPadding=10,
                            leftIndent=10, rightIndent=10, spaceBefore=7, spaceAfter=15),
    'inventory': ParagraphStyle('inventory', fontName='Arial', fontSize=8.4, leading=11.5,
                                textColor=INK, leftIndent=10, firstLineIndent=-10,
                                spaceAfter=7, splitLongWords=True, allowWidows=0),
    'list': ParagraphStyle('list', fontName='Arial', fontSize=9.7, leading=14.1,
                           textColor=INK, leftIndent=12, firstLineIndent=-12,
                           spaceAfter=8, splitLongWords=True, allowWidows=0),
}

TOKEN = re.compile(r'(`[^`]+`|\[[^\]]+\]\([^\)]+\)|\*\*.+?\*\*|(?<!\*)\*[^*]+\*(?!\*))')


def inline(text):
    """Translate only the Markdown constructs actually used by the report."""
    out = []
    cursor = 0
    for match in TOKEN.finditer(text):
        out.append(html.escape(text[cursor:match.start()]))
        part = match.group(0)
        if part.startswith('`'):
            out.append('<font color="#35596D">' + html.escape(part[1:-1]) + '</font>')
        elif part.startswith('['):
            label, target = re.match(r'\[([^\]]+)\]\(([^\)]+)\)', part).groups()
            if not target.startswith(('https://', 'http://', '#')):
                target_path = Path(unquote(target))
                if not target_path.is_absolute():
                    target_path = SOURCE.parent / target_path
                target = target_path.resolve().as_uri()
            out.append('<link href="' + html.escape(target, quote=True) + '" color="#245873">'
                       + inline(label) + '</link>')
        elif part.startswith('**'):
            out.append('<b>' + inline(part[2:-2]) + '</b>')
        else:
            out.append('<i>' + inline(part[1:-1]) + '</i>')
        cursor = match.end()
    out.append(html.escape(text[cursor:]))
    return ''.join(out)


def para(text, style='body'):
    return Paragraph(inline(text), styles[style])


class ReportDoc(BaseDocTemplate):
    def __init__(self, filename):
        super().__init__(str(filename), pagesize=A4, leftMargin=MARGIN,
                         rightMargin=MARGIN, topMargin=49, bottomMargin=45,
                         title='MAUVE Paper II: Findings Synthesis and Paper Plan',
                         author='Research synthesis for Rongjun Huang',
                         subject='Source-based review of the further workspace, 11 October 2026')
        frame = Frame(MARGIN, 45, WIDTH, PAGE_H-94, leftPadding=0, rightPadding=0,
                      topPadding=0, bottomPadding=0)
        self.addPageTemplates([
            PageTemplate(id='cover', frames=frame, onPage=self.cover),
            PageTemplate(id='body', frames=frame, onPage=self.header),
        ])

    def cover(self, canvas, doc):
        canvas.saveState()
        canvas.setFillColor(BLUE)
        canvas.rect(MARGIN, PAGE_H-38, 78, 3, fill=1, stroke=0)
        canvas.setFont('Arial', 8)
        canvas.setFillColor(MUTED)
        canvas.drawString(MARGIN, 27, 'MAUVE  /  RESEARCH SYNTHESIS  /  11 OCTOBER 2026')
        canvas.restoreState()

    def header(self, canvas, doc):
        canvas.saveState()
        canvas.setFont('Arial-Bold', 8)
        canvas.setFillColor(MUTED)
        canvas.drawString(MARGIN, PAGE_H-29, 'MAUVE  /  PAPER II')
        canvas.setFont('Arial', 8)
        canvas.drawRightString(PAGE_W-MARGIN, PAGE_H-29, 'Findings synthesis and paper plan')
        canvas.setStrokeColor(colors.HexColor('#D4E0E6'))
        canvas.line(MARGIN, PAGE_H-36, PAGE_W-MARGIN, PAGE_H-36)
        canvas.setFont('Arial', 8)
        canvas.drawString(MARGIN, 26, '11 October 2026  |  Current 8900 Angstrom evidence')
        canvas.drawRightString(PAGE_W-MARGIN, 26, str(doc.page))
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if isinstance(flowable, Paragraph) and hasattr(flowable, 'section_key'):
            self.canv.bookmarkPage(flowable.section_key)
            self.canv.addOutlineEntry(flowable.getPlainText(), flowable.section_key, 0, False)
            self.notify('TOCEntry', (0, flowable.getPlainText(), self.page, flowable.section_key))


def make_table(lines):
    raw = [[cell.strip() for cell in row.strip().strip('|').split('|')] for row in lines]
    raw = [row for row in raw if not all(re.fullmatch(r'[:\- ]+', cell) for cell in row)]
    ncols = len(raw[0])
    weights = {2:[.28,.72], 3:[.25,.36,.39], 4:[.40,.20,.20,.20], 5:[.24,.19,.19,.19,.19]}[ncols]
    content = [[para(cell, 'tablehead' if i == 0 else 'table') for cell in row] for i,row in enumerate(raw)]
    table = Table(content, colWidths=[WIDTH*w for w in weights], repeatRows=1, hAlign='LEFT')
    table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), BLUE),
        ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, PALE]),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('LEFTPADDING', (0,0), (-1,-1), 7), ('RIGHTPADDING', (0,0), (-1,-1), 7),
        ('TOPPADDING', (0,0), (-1,-1), 7), ('BOTTOMPADDING', (0,0), (-1,-1), 7),
        ('LINEBELOW', (0,-1), (-1,-1), .4, colors.HexColor('#D4E0E6')),
    ]))
    return [table, Spacer(1,10)]


def build():
    lines = SOURCE.read_text().splitlines()
    intro_end = next(i for i,line in enumerate(lines) if line.startswith('## '))
    intro = '\n'.join(lines[1:intro_end]).strip().split('\n\n')
    story = [Spacer(1,40), para('MAUVE-MUSE / SECOND PAPER', 'kicker'),
             para(lines[0][2:], 'title')]
    story.extend(para(item) for item in intro)
    story.extend([Spacer(1,25), para('A complete notebook register, a source-traceable synthesis, '
                                  'and a plan for the 7000 Angstrom transition.', 'kicker'),
                  NextPageTemplate('body'), PageBreak(), para('Contents', 'h2')])
    toc = TableOfContents()
    toc.levelStyles = [ParagraphStyle('toc0', fontName='Arial', fontSize=10, leading=15,
                                     leftIndent=0, firstLineIndent=0, textColor=INK,
                                     spaceBefore=7, rightIndent=22)]
    story.extend([toc, PageBreak()])
    i = intro_end
    appendix = False
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.startswith('## '):
            heading = line[3:]
            if heading.startswith('Appendix'):
                appendix = True
                story.append(PageBreak())
            else:
                story.append(CondPageBreak(120))
            p = para(heading, 'h2')
            p.section_key = 'section-' + str(i)
            story.append(p)
            i += 1
        elif line.startswith('### '):
            story.append(para(line[4:], 'h3'))
            i += 1
        elif line.startswith('|'):
            table_lines = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                table_lines.append(lines[i]); i += 1
            group = make_table(table_lines)
            while story and isinstance(story[-1], Paragraph) and story[-1].getKeepWithNext():
                group.insert(0, story.pop())
            story.append(KeepTogether(group))
        elif line.startswith('!['):
            caption, target = re.match(r'!\[(.*)\]\((.*)\)',line).groups()
            path = SOURCE.parent / unquote(target)
            with PILImage.open(path) as im:
                w,h = im.size
            graphic = Image(str(path), width=WIDTH, height=WIDTH*h/w)
            story.append(KeepTogether([Spacer(1,6), graphic, para(caption, 'caption')]))
            i += 1
        elif line.startswith('> '):
            story.append(para(line[2:], 'quote')); i += 1
        elif line.startswith('- '):
            style = 'inventory' if appendix else 'list'
            story.append(para('- ' + line[2:], style)); i += 1
        elif re.match(r'^\d+\. ', line):
            story.append(para(line, 'list')); i += 1
        else:
            paragraph = [line]; i += 1
            while i < len(lines) and lines[i].strip() and not re.match(r'^(#|\||!\[|> |\- |\d+\. )', lines[i]):
                paragraph.append(lines[i].strip()); i += 1
            p = para(' '.join(paragraph))
            lookahead = i
            while lookahead < len(lines) and not lines[lookahead].strip():
                lookahead += 1
            if lookahead < len(lines) and lines[lookahead].startswith('|') and len(' '.join(paragraph)) < 350:
                p.keepWithNext = True
            story.append(p)
    ReportDoc(TARGET).multiBuild(story)
    print('PDF_WRITTEN', TARGET)


if __name__ == '__main__':
    build()
