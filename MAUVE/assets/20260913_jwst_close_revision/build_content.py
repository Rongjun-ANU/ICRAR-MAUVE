from pathlib import Path
import json
import re

ROOT = Path('/Users/Igniz/Desktop/ICRAR/MAUVE')
ASSETS = ROOT / 'assets/20260913_jwst_close_revision'
STEM = '20260913_MAUVE_JWST_proposal_close_revision'
source = json.loads((ASSETS / 'source.json').read_text())
text = source['text']
GREEN = '#147D3B'
blocks = []
changes = []

def green(t):
    return '<span style="color: %s;">%s</span>' % (GREEN, t)

def s(start, end):
    a = text.index(start)
    b = text.index(end, a + len(start))
    return text[a:b].strip()

def replace(t, old, new):
    assert old in t, old
    changes.append({'old':old, 'new':new})
    return t.replace(old, green(new))

def add(t, kind='p', **kw):
    blocks.append({'kind':kind, 'text':t, **kw})

def new(t, kind='p', **kw):
    add(green(t), kind, **kw)

def fig(n, caption):
    item = next(i for i in source['figures'] if i['number']==n)
    blocks.append({'kind':'figure', 'text':caption, **item})

def bold_start(t, stop):
    k = t.index(stop)+len(stop)
    return '<b>'+t[:k]+'</b>'+t[k:]

add('This page is for your information and will be removed before submission', 'cover_note')
add('<b>Title:</b> The Evolving ISM and Stellar Populations in Virgo Cluster Galaxies', 'cover')
abstract = s('Abstract:', 'Scientific Category:')
abstract = replace(abstract, 'infall stages into the cluster', 'infall stages into the Virgo galaxy cluster')
abstract = replace(abstract, 'galaxy evolution in cluster environments', 'galaxy evolution in galaxy cluster environments')
abstract = replace(abstract, 'determine precise relative distances to all targets, constraining their orbital histories in Virgo as well as the assembly of the Virgo Cluster itself.', 'determine how environmental processing changes the emergence of young stellar clusters from their natal clouds, connecting local stellar feedback to galaxy-scale gas loss.')
abstract = abstract.replace('Abstract:', '<b>Abstract:</b>', 1)
add(abstract, 'cover')
add('<b>Scientific Category:</b> Nearby Galaxies to Cosmic Noon', 'cover')
add('<b>Alternate Category:</b> Gas, Dust and the ISM', 'cover')
add('<b>Science Keywords:</b> Disk galaxies, Galaxy clusters, Galaxy environments, Galaxy evolution, Interstellar dust, Star formation, Stellar populations', 'cover')
new('<b>Science / Parallel / Charged Time:</b> Cycle 5 reference: 55.7 h / 7.6 h / 156.3 h. The revised allocation requires updated ETC and APT calculations.', 'cover')
add('<b>Sample Information:</b> <a href="https://docs.google.com/spreadsheets/d/1Csa-zoVhpEk_VxK_cwWB_K5jGOpy93r0ora8CrhHP2I/edit?usp=sharing">'+green('Google Spreadsheet')+'</a>', 'cover')
new('<b>13 September 2026 revision:</b> Green text marks additions and replacements relative to Cycle 5. Original figures and section order are retained. The original technical estimates are retained for drafting; TRGB depth and halo coverage no longer drive the request.', 'cover')
add('', 'pagebreak')

add('Scientific Justification (required for all)', 'h2')
add(s('Galaxy evolution is governed', 'Over the past decades,'))
intro = s('Over the past decades,', 'We propose a JWST Treasury program')
intro = replace(intro, 'cluster satellites', 'satellite galaxies')
intro = bold_start(intro, 'and important questions are unanswered:') if False else intro
q = intro.index('How does star formation')
intro = intro[:q]+'<b>'+intro[q:]+'</b>'
add(intro)
intro2 = s('We propose a JWST Treasury program', 'This program is uniquely timed')
intro2 = replace(intro2, 'Virgo Cluster', 'Virgo galaxy cluster')
intro2 = replace(intro2, 'We will also get precise distances for all targets, mapping out their 3D distribution and connecting their infall with Virgo’s global assembly, a feat only possible with JWST in Virgo [d=16.2 Mpc; 4].', 'We will also determine how young stellar clusters emerge from their natal clouds as the surrounding ISM is compressed or stripped, linking local feedback to gas loss in Virgo [d=16.2 Mpc; 4].')
intro2 = intro2.replace('All this science will be achieved with one integrated program', '<b>All this science will be achieved with one integrated program</b>')
add(intro2)
intro3 = s('This program is uniquely timed', 'The 40 target galaxies')
add(intro3)
sample = s('The 40 target galaxies', 'Fig. 1:')
sample = replace(sample, 'A stellar mass-matched field control sample also exists', 'A stellar mass-matched nearby comparison sample also exists')
sample = replace(sample, 'enabling rigorous distinction between environmental and internal processes across the full range of galaxy conditions.', 'enabling comparisons across the full range of galaxy conditions. We will select controls by their gas content and disturbance, accounting for Virgo members already present in PHANGS.')
add(sample)
new('J-Virgo [81] overlaps eight targets and supplies valuable stellar-continuum imaging; our recombination-line, PAH, and MIRI data reveal how those populations interact with their ISM.')

cap1 = s('Fig. 1:', 'Goal #1.')
cap1 = replace(cap1, 'for field galaxies (from PHANGS;', 'for nearby comparison galaxies (from PHANGS;')
cap1 = replace(cap1, 'provide an ideal control sample for separating internal versus environmental effects.', 'provide candidate controls for separating internal versus environmental effects, after checking their environmental disturbance.')
fig(1, cap1)
g1 = s('Goal #1.', 'Fig. 2:')
g1 = replace(g1, 'with the cluster age distribution', 'with the stellar cluster age distribution')
g1 = bold_start(g1, 'IR SFR tracers.')
add(g1, 'goal')
fig(2, s('Fig. 2:', 'We will also derive recent SFR'))
add(s('We will also derive recent SFR', 'Goal #2.'))
g2 = s('Goal #2.', 'The sharp, sensitive PAH maps')
g2 = bold_start(g2, 'at high resolution.')
add(g2, 'goal')
g2b = s('The sharp, sensitive PAH maps', 'When low column density gas')
g2b = replace(g2b, 'the field control sample', 'the selected comparison sample')
add(g2b)
out = s('When low column density gas', 'Fig. 3:')
out = replace(out, 'stripped from cluster galaxies', 'stripped from galaxies in Virgo')
add(out)
fig(3, s('Fig. 3:', 'systems [11, 71].'))
add(s('systems [11, 71].', 'As the ISM gets stripped'))
dust = s('As the ISM gets stripped', 'Goal #3.')
dust = replace(dust, 'targeting more distant cluster lensing fields', 'targeting fields lensed by more distant galaxy clusters')
add(dust)

new('<b>Goal #3. Young stellar feedback in an externally disturbed ISM.</b> We will determine how environmental gas compression and removal change the emergence of young stellar clusters from their natal clouds. In Virgo, stellar feedback acts on an ISM already being reshaped by the ICM. Compression may retain obscuring material, while stripping may expose young stars or remove gas made more diffuse by their feedback [20, 35, 86]. <b>Do young stellar clusters remain embedded for longer in compressed regions, or emerge earlier where gas is being removed?</b>', 'goal')
new('Our NIRCam and MIRI imaging, combined with HST, will uncover embedded young stellar clusters and connect them to the gas and dust around them, building on recent JWST studies of their emergence [33, 87]. We will compare embedded and exposed populations of similar age and mass across compressed, stripped, and less disturbed regions. The recent SFHs from Goal #1 and the evolving ISM mapped in Goal #2 will place these differences in the broader history of each galaxy. Together, these measurements will establish whether galaxy-scale gas loss is accompanied by a change in how young stellar populations emerge from their birth material, linking local stellar feedback to environmental quenching.')

treasury = s('Why Treasury:', 'Technical Justification')
treasury = treasury.replace('Why Treasury:', '<b>Why Treasury:</b>', 1)
treasury = replace(treasury, 'host systems of other transient events [57].', 'host systems of other transient events [57]. TRGB distances from adequate archival or incidental imaging, including J-Virgo [81], will add three-dimensional context as a supporting science product.')
treasury = replace(treasury, 'local clusters', 'local galaxy clusters')
treasury = replace(treasury, 'higher-z clusters', 'higher-z galaxy clusters')
add(treasury)
add('Technical Justification (required for GO, DD and Survey only)', 'h2')
obs = s('Observation Design:', 'We use NIRCam’s F090W')
add(obs.replace('Observation Design:', '<b>Observation Design:</b>'), 'p')
continuum = s('We use NIRCam’s F090W', 'We use NIRCam’s F335M')
continuum = replace(continuum, 'M⊙ clusters', 'M⊙ stellar clusters')
continuum = replace(continuum, 'less embedded clusters', 'less embedded stellar clusters')
continuum = replace(continuum, 'young clusters', 'young stellar clusters')
continuum = replace(continuum, 'These setups directly support Goal #1.', 'These setups directly support Goals #1 and #3. The F090W/F150W depths will be set by recovery of stellar populations, using adequate archival imaging first; no dedicated TRGB depth is required.')
add(continuum)
add(s('We use NIRCam’s F335M', 'We use NIRCam’s F187N'))
add(s('We use NIRCam’s F187N', 'We use NIRCam’s F090W and F150W filters to detect RGB'))
new('For Goal #3, we use the same stellar-continuum, recombination-line, and PAH images to identify embedded and exposed young stellar populations and characterize their surroundings [33, 87]. The line-based comparison will use targets and regions with adequate Paα or Brα transmission. Artificial-source tests will establish the common stellar-mass and completeness limits, while nearby comparison images will be matched to the Virgo resolution and depth. We place the NIRCam modules to cover the star-forming disk and environmentally disturbed regions; neither dedicated halo coverage nor reaching below the TRGB is a requirement.')
data = s('Data Processing and Dissemination:', 'Special Requirements (if any)')
data = data.replace('Data Processing and Dissemination:', '<b>Data Processing and Dissemination:</b>')
data = replace(data, 'star cluster / H II region catalogs and value-added measurements', 'stellar cluster / H II region catalogs, young stellar emergence classifications, and value-added measurements')
data = replace(data, 'published online at CADC', 'delivered to MAST with a mirror at CADC')
add(data)

add('Special Requirements (if any)', 'h2')
special = s('We constrain the position angle', 'We request the NIRCam and MIRI observations')
special = replace(special, 'All observations remain schedulable under these PA constraints.', 'The revised pointings and PA constraints will be checked for schedulability in APT.')
special = replace(special, 'guarantee coverage of both the galaxy disk and part of the outer halo (for distance measurements).', 'cover the star-forming disk and environmentally disturbed regions, using adequate archival coverage where available.')
add(special)
add(s('We request the NIRCam and MIRI observations', 'Justify Coordinated Parallel Observations'))
add('Justify Coordinated Parallel Observations (if any, GO/DD only)', 'h2')
add(s('We require coordinated parallels', 'Justify Duplications'))
add('Justify Duplications (if any, GO/DD only)', 'h2')
new('The Cycle 5 archive check identified adequate imaging for 12/40 targets, mostly with MIRI, and omitted those duplicated observations. We will update this assessment filter by filter for the revised footprint and science requirements, using adequate archival data in place of new exposures.')
new('(1) <b>J-Virgo GO 7763 [81] overlaps eight targets:</b> NGC 4294, NGC 4351, NGC 4388, NGC 4402, NGC 4548, NGC 4569, NGC 4579, and NGC 4606. Its NIRCam F115W/F150W/F277W imaging provides valuable stellar-continuum and distance information, including disk coverage. It does not provide the targeted infrared recombination lines, 3.3 μm PAH band, or MIRI imaging required here. We will use its adequate continuum data and restrict additional continuum exposures to uncovered regions or demonstrated depth requirements for Goals #1 and #3.')
new('(2) <b>F150W and F187N from GO 3707 [PHANGS; 39] for eight targets:</b> We will use existing F150W wherever it supports the stellar-population measurements. Any additional Paα imaging will be justified by the required footprint and demonstrated line sensitivity after continuum subtraction. TRGB depth and halo coverage are no longer reasons for repeating these observations.')
new('(3) <b>F300M and F335M from GO 2107 [37], 3707 [39], and 4793 [67]:</b> We will use adequate archival data and add the missing disk and disturbed-region coverage. Exposure pairings will be set by stellar-population and PAH requirements, with no TRGB-driven integration.')
new('(4) <b>HST/ACS F814W parallel imaging from GO 18103 [72]:</b> These data will add stellar-population constraints and incidental TRGB measurements where feasible; halo coverage and TRGB precision no longer drive this program.')

add('', 'pagebreak')
add('References', 'h2')
refs = text[text.index('References [1]')+len('References '):]
refitems = re.findall(r'\[(\d+)\]\s+(.*?)(?=\s+\[\d+\]\s+|$)', refs)
assert len(refitems)==85, len(refitems)
for n, body in refitems:
    add('['+n+'] '+body, 'reference')
new('[86] Cramer, W. J., Kenney, J. D. P., Cortes, J. R., et al. 2020, ApJ, 901, 95. <a href="https://doi.org/10.3847/1538-4357/abaf54">doi:10.3847/1538-4357/abaf54</a>', 'reference')
new('[87] Pedrini, A., Adamo, A., Calzetti, D., et al. 2026, Nature Astronomy, 10, 1179. <a href="https://doi.org/10.1038/s41550-026-02857-y">doi:10.1038/s41550-026-02857-y</a>', 'reference')
new('Draft preparation: OpenAI Codex desktop, GPT-6-based assistant, 13 September 2026; used to adapt the supplied Cycle 5 text, check the added references and J-Virgo overlap, and prepare this marked revision. Final scientific and submission review rests with the proposing team.', 'disclosure')

# Consistent, legible mathematical notation; these are extraction repairs, not prose revisions.
for b in blocks:
    b['text'] = b['text'].replace('3.3−11.3µm', '3.3−11.3 μm').replace('µm', 'μm')
    b['text'] = b['text'].replace('∼100', '∼100').replace(' ≲10', ' ≲10')
    b['text'] = b['text'].replace('1−2', '1–2').replace('9−30', '9–30').replace('20−30', '20–30')
    b['text'] = b['text'].replace('M⊙pc', 'M⊙ pc')

md = []
for b in blocks:
    k, t = b['kind'], b['text']
    if k == 'pagebreak': md.append('<!-- pagebreak -->')
    elif k == 'h2': md.append('## '+t)
    elif k == 'figure':
        rel = Path(b['path']).relative_to(ROOT)
        md.append('![Original Cycle 5 Figure %d](%s)\n\n%s' % (b['number'], rel.as_posix(), t))
    else: md.append(t)
(ROOT / (STEM+'.md')).write_text('\n\n'.join(md)+'\n')
(ASSETS / 'content.json').write_text(json.dumps({'blocks':blocks, 'changes':changes,
    'stem':STEM, 'green':GREEN, 'source_sha256':source['source_sha256']}, ensure_ascii=False, indent=2))
print('Wrote Markdown:', ROOT / (STEM+'.md'))
print('Blocks:', len(blocks), 'References:',len(refitems)+2)
