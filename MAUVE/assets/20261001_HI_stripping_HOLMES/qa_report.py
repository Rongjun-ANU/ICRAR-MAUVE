"""Structural PDF checks and contact sheets; visual acceptance remains manual."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re

import fitz
from PIL import Image, ImageDraw, ImageFont


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf", type=Path)
    parser.add_argument("qa_dir", type=Path)
    args = parser.parse_args()
    args.qa_dir.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(args.pdf)
    pages, outside, bad_links, replacements = [], [], [], []
    link_counts = {"internal": 0, "external": 0}
    for i, page in enumerate(doc):
        txt = page.get_text()
        compact = re.sub(r"\s+", " ", txt).strip()
        pages.append({"page": i + 1, "characters": len(txt),
                      "start": compact[:220], "end": compact[-180:]})
        for word in page.get_text("words"):
            if word[0] < -0.5 or word[1] < -0.5 or word[2] > page.rect.width + 0.5 or word[3] > page.rect.height + 0.5:
                outside.append({"page": i + 1, "word": word[:5]})
        if "\ufffd" in txt:
            replacements.append(i + 1)
        for link in page.get_links():
            if link["kind"] in (fitz.LINK_GOTO, fitz.LINK_NAMED):
                link_counts["internal"] += 1
                if not 0 <= link.get("page", -1) < len(doc):
                    bad_links.append({"page": i + 1, "link": link})
            elif link["kind"] == fitz.LINK_URI:
                link_counts["external"] += 1
                if link.get("uri", "").startswith(("file:", "/private/tmp")):
                    bad_links.append({"page": i + 1, "link": link})
            elif link["kind"] == fitz.LINK_LAUNCH:
                bad_links.append({"page": i + 1, "link": repr(link)})
    audit = {
        "pdf": str(args.pdf.resolve()),
        "sha256": hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
        "pages": len(doc),
        "outline_entries": len(doc.get_toc()),
        "link_counts": link_counts,
        "invalid_links": bad_links,
        "out_of_page_words": outside,
        "replacement_character_pages": replacements,
        "short_pages_under_250_chars": [p["page"] for p in pages if p["characters"] < 250],
        "page_text_summaries": pages,
    }
    (args.qa_dir / "pdf_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({k: v for k, v in audit.items() if k != "page_text_summaries"}, indent=2))
    print("\nPAGE MAP")
    for p in pages:
        print(f'{p["page"]:02d} ({p["characters"]:4d} chars) {p["start"][:160]}')
    rendered = sorted(args.qa_dir.glob("page-*.png"))
    if rendered:
        assert len(rendered) == len(doc), (len(rendered), len(doc))
        columns, rows, thumb_w, label_h, gutter = 3, 4, 350, 22, 12
        thumb_h = round(thumb_w * 297 / 210)
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 15)
        sheet_paths = []
        for start in range(0, len(rendered), columns * rows):
            subset = rendered[start:start + columns * rows]
            used_rows = math.ceil(len(subset) / columns)
            sheet = Image.new("RGB", (columns * (thumb_w + gutter) + gutter,
                                      used_rows * (thumb_h + label_h + gutter) + gutter), "#dfe5e7")
            draw = ImageDraw.Draw(sheet)
            for j, image_path in enumerate(subset):
                with Image.open(image_path) as im:
                    im = im.convert("RGB")
                    im.thumbnail((thumb_w, thumb_h), Image.Resampling.LANCZOS)
                    x = gutter + (j % columns) * (thumb_w + gutter)
                    y = gutter + (j // columns) * (thumb_h + label_h + gutter)
                    sheet.paste(im, (x, y + label_h))
                    draw.text((x, y), f"PAGE {start + j + 1:02d}", font=font, fill="#173f4b")
            output = args.qa_dir / f"contact-{start // (columns * rows) + 1:02d}.png"
            sheet.save(output)
            sheet_paths.append(str(output))
        print("\nCONTACT SHEETS\n" + "\n".join(sheet_paths))
    assert not outside, "Text outside the PDF page."
    assert not replacements, "Missing/replacement glyphs in extracted text."
    assert not bad_links, "Broken internal or temporary-file links."
    assert all(p["characters"] > 50 for p in pages), "Empty or near-empty page."


if __name__ == "__main__":
    main()
