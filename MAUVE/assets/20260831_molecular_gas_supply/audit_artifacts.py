"""Verify report structure and build all-page contact sheets for visual QA."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import unicodedata
from PIL import Image, ImageOps, ImageDraw, ImageFont
from pypdf import PdfReader

ASSET = Path(__file__).resolve().parent
MAUVE = ASSET.parent.parent
STEM = "20260831_Molecular_Gas_Supply_Definitions_and_Observational_Inference"


def norm(text):
    text = unicodedata.normalize("NFKD",text)
    return "".join(c.lower() for c in text if c.isalnum())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("render_dir",type=Path)
    args=parser.parse_args()
    md=MAUVE/f"{STEM}.md"; pdf=MAUVE/f"{STEM}.pdf"
    canonical=ASSET/"report-source.md"
    source=md.read_text()
    assert md.read_bytes()==canonical.read_bytes(), "Delivered MD differs from canonical source"
    assert not [c for c in source if ord(c)<32 and c not in "\n\t"], "Unexpected control characters"
    tags=re.findall(r"\\tag\{(\d+)\}",source)
    assert tags==[str(i) for i in range(1,61)], "Equation numbering incomplete or duplicated"
    assert source.count("\n$$\n")==120, "Display-math delimiters are not paired"
    assert not re.search(r"\\rm\b|\\hbox\b",source), "Unsupported legacy math command"
    figures=re.findall(r"!\[[^\]]*\]\(([^)]+)\)",source)
    assert len(figures)==2 and all((MAUVE/f).is_file() for f in figures), "Missing figure"
    doc=PdfReader(pdf)
    texts=[page.extract_text() or "" for page in doc.pages]
    assert all(len(text.strip())>50 for text in texts), "Blank or near-blank page"
    pdftext="\n".join(texts)
    assert not re.search(r"\\frac|\\mathrm|\\tag\{|\\Sigma|\\begin",pdftext), "Unrendered TeX in PDF"
    assert "\ufffd" not in pdftext, "PDF text extraction replacement glyph"
    for number in range(1,61):
        assert f"({number})" in pdftext, f"Missing visible equation number {number}"
    headings=[]
    for line in source.splitlines():
        if re.match(r"^#{1,2} ",line):
            heading=re.sub(r"^#+ ","",line).replace("$_2$","2")
            headings.append(heading)
    missing=[h for h in headings if norm(h) not in norm(pdftext)]
    assert not missing, f"Missing PDF headings: {missing}"
    dom=json.loads((ASSET/"verification_dom.json").read_text())
    for key in ["mathErrors","unconvertedMath","brokenAnchors","brokenImages","overflow","nonWebLinks"]:
        assert not dom[key], f"DOM check failed: {key}"
    assert dom["displayMath"]==60 and dom["equationNumbers"]==[f"({i})" for i in range(1,61)]
    links=[]
    for page in doc.pages:
        for annotation in page.get("/Annots",[]):
            obj=annotation.get_object()
            if "/A" in obj and obj["/A"].get("/URI"):
                links.append(str(obj["/A"]["/URI"]))
    assert links, "No live PDF web links"
    renders=sorted(args.render_dir.glob("page-*.png"))
    assert len(renders)==len(doc.pages), "Not all PDF pages were rendered"
    sheets=[]
    label_font=ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc",18)
    for start in range(0,len(renders),9):
        canvas=Image.new("RGB",(1260,1854),"#dce4e6")
        draw=ImageDraw.Draw(canvas)
        for local,path in enumerate(renders[start:start+9]):
            thumb=ImageOps.contain(Image.open(path).convert("RGB"),(400,566))
            x=10+(local%3)*420; y=28+(local//3)*618
            canvas.paste(thumb,(x,y))
            draw.text((x,y+572),f"Page {start+local+1}",font=label_font,fill="#203f4b")
        output=args.render_dir/f"contact-{start//9+1}.png"
        canvas.save(output); sheets.append(str(output))
    visual_path=ASSET/"verification_visual.json"
    visual=json.loads(visual_path.read_text()) if visual_path.is_file() else {}
    visual_current=visual.get("status")=="passed" and visual.get("pdf_sha256")==sha(pdf)
    visual_summary=(f"Passed: all-page contact sheets and {len(visual.get('high_resolution_pages',[]))} selected 1600-pixel pages inspected; see verification_visual.json"
                    if visual_current else "Pending: no passing visual record matches this PDF hash")
    result={"status":"passed","pdf_pages":len(doc.pages),"source_words":len(source.split()),
        "display_equations":60,"equation_labels":60,"section_headings_checked":len(headings),
        "tables":dom["tables"],"figures":len(figures),"pdf_web_link_annotations":len(links),
        "unique_pdf_web_targets":len(set(links)),"md_equals_canonical":True,
        "markdown_sha256":sha(md),"pdf_sha256":sha(pdf),
        "figure_sha256":{name:sha(MAUVE/name) for name in figures},
        "rendered_pages":len(renders),"contact_sheets":sheets,
        "visual_qa":visual_summary,"visual_record_matches_pdf":visual_current,
        "limits":"Structural checks do not establish scientific closure validity or all URL reachability.",
        "page_text_inventory":[{"page":i+1,"characters":len(t),"start":" ".join(t.split())[:160]} for i,t in enumerate(texts)]}
    (ASSET/"verification_artifacts.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__": main()
