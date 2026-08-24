#!/usr/bin/env python3
"""Apply A4 geometry and a page-number footer without rewriting DOCX XML."""

from __future__ import annotations

import tempfile
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "tmp/pdfs/20260823 RPS Full Derivation.docx"
OUTPUT = ROOT / "tmp/pdfs/20260823 Finite-Thickness Arbitrary-Orientation Ram Pressure Stripping - Full Derivation.docx"
RID = "rIdFooter20260823"


def replace_once(text: str, old: str, new: str, label: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"Expected one {label} marker; found {text.count(old)}")
    return text.replace(old, new, 1)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="rps_docx_") as temp_name:
        temp = Path(temp_name)
        with zipfile.ZipFile(SOURCE) as archive:
            archive.extractall(temp)

        document_path = temp / "word/document.xml"
        document = document_path.read_text(encoding="utf-8")
        document = replace_once(
            document,
            "<w:sectPr>",
            f'<w:sectPr><w:footerReference w:type="default" r:id="{RID}"/>',
            "section start",
        )
        document = replace_once(
            document,
            "</w:sectPr>",
            '<w:pgSz w:w="11906" w:h="16838"/>'
            '<w:pgMar w:top="1134" w:right="1020" w:bottom="1134" w:left="1134" '
            'w:header="567" w:footer="567" w:gutter="0"/>'
            "</w:sectPr>",
            "section end",
        )
        document_path.write_text(document, encoding="utf-8")

        rels_path = temp / "word/_rels/document.xml.rels"
        rels = rels_path.read_text(encoding="utf-8")
        rels = replace_once(
            rels,
            "</Relationships>",
            f'<Relationship Id="{RID}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" '
            'Target="footer1.xml"/></Relationships>',
            "relationship end",
        )
        rels_path.write_text(rels, encoding="utf-8")

        types_path = temp / "[Content_Types].xml"
        types = types_path.read_text(encoding="utf-8")
        types = replace_once(
            types,
            "</Types>",
            '<Override PartName="/word/footer1.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/>'
            "</Types>",
            "content-types end",
        )
        types_path.write_text(types, encoding="utf-8")

        footer = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"
       xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <w:p>
    <w:pPr>
      <w:pBdr><w:top w:val="single" w:sz="4" w:space="4" w:color="D8D4CB"/></w:pBdr>
      <w:jc w:val="center"/>
    </w:pPr>
    <w:r><w:rPr><w:color w:val="5C646B"/><w:sz w:val="16"/></w:rPr><w:t xml:space="preserve">MAUVE | RPS derivation | 2026-08-23 | </w:t></w:r>
    <w:fldSimple w:instr=" PAGE "><w:r><w:rPr><w:color w:val="5C646B"/><w:sz w:val="16"/></w:rPr><w:t>1</w:t></w:r></w:fldSimple>
  </w:p>
</w:ftr>
'''
        (temp / "word/footer1.xml").write_text(footer, encoding="utf-8")

        with zipfile.ZipFile(OUTPUT, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(temp.rglob("*")):
                if path.is_file():
                    archive.write(path, path.relative_to(temp))

    print(OUTPUT)


if __name__ == "__main__":
    main()
