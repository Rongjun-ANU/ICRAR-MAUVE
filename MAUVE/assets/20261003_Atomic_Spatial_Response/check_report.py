"""Structural and unchanged-content verification for the focused revision."""
from pathlib import Path
import hashlib,json,re
OUT=Path(__file__).resolve().parent
ROOT=OUT.parent.parent
report=ROOT/'20261002_Model0_Gas_Fading_and_HOLMES_Derivation.md'
before=(OUT/'report_before_simplification.md').read_text()
after=report.read_text()
pattern=re.compile(r'\b[Ee]quations?\s+\(\d+\)(?:(?:--|\s+and\s+|,\s*)\(\d+\))*')
def shift_old(text):
    def shifted(n):return str(int(n)-10) if int(n)>=41 else n
    text=re.sub(r'\\tag\{(\d+)\}',lambda m:r'\tag{'+shifted(m[1])+'}',text)
    return pattern.sub(lambda m:re.sub(r'\((\d+)\)',lambda n:'('+shifted(n[1])+')',m[0]),text)
def between(text,start,end):return text[text.index(start):text.index(end,text.index(start))]
unchanged=[]
for start,end in [('# 5.','# 8.'),('## 8.1 ','## 8.3 '),('## 8.4 ','# 9.'),
                  ('# Appendix A.','# Appendix B.'),('# Appendix C.','# Appendix G.')]:
    assert shift_old(between(before,start,end))==between(after,start,end),(start,end)
    unchanged.append(start+' through before '+end)
assert before[:before.index('## 3.5 ')]==after[:after.index('## 3.5 ')]
assert re.findall(r'\\tag\{(\d+)\}',after)==[str(n) for n in range(1,85)]
section=between(after,'## 3.5 ','# 4.')
assert section.count('\\tag')==5 and '### ' not in section
assert not re.search(r'(?<![A-Za-z])j(?![A-Za-z])',section)
assert not any(s in after for s in ['References introduced in this section',
    'Now raise both initial gas columns','the spatial example uses 1.5 for each','1.1980','0.1256'])
assert after.index('<span id="ref-lee">')>after.index('# References')
assert after.index('<span id="ref-cramer">')>after.index('# References')
# Verify all equation citations resolve; range endpoints are checked too.
for match in pattern.finditer(after):
    for n in re.findall(r'\((\d+)\)',match[0]):assert 1<=int(n)<=84,(match[0],n)
image_paths=re.findall(r'^!\[.*\]\(([^)]+)\)$',after,re.MULTILINE)
assert len(image_paths)==3
assert all((ROOT/p).exists() for p in image_paths)
assert '20261003_Atomic_Spatial_Response' in image_paths[0]
prior_hashes=json.loads((OUT/'before_hashes.json').read_text())
for path,digest in prior_hashes.items():
    if path.startswith('assets/'):
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
audit=json.loads((OUT/'spatial_audit.json').read_text())
for r in audit['cases']:
    assert f"{r['sfr_fraction_1gyr']:.5f}" in after
    if r['region']!='reference':assert f"{r['offset_dex_1gyr']:+.5f}" in after
assert '0.09839' in after and '1.00686' in after
result=dict(unchanged_before_section35=True,unchanged_except_renumbering=unchanged,
    existing_assets_unchanged=True,section35_equations=5,total_equations=84,
    local_j_absent=True,subsubsections_removed=True,references_moved=True,
    figure_and_numeric_table_match_audit=True,md_sha256=hashlib.sha256(report.read_bytes()).hexdigest())
(OUT/'report_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
