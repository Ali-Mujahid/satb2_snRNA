"""One-time, bounded migration of the existing communication page to strict DEG-only.

No biological input or statistical fit is changed.
"""
from pathlib import Path
import hashlib
import re

ROOT = Path(__file__).resolve().parent.parent
PAGE = ROOT / 'site/communication/index.html'


def migrate():
    text = PAGE.read_text()
    if 'id="ligandDirection"' in text:
        raise SystemExit("This one-time migration has already been applied.")
    text = text.replace('SATB2 · Ligand–receptor explorer', 'SATB2 · DEG-only ligand–receptor explorer')
    text = text.replace('CA1 ↔ inhibitory communication candidates', 'DEG-only CA1 ↔ inhibitory candidates')
    text = text.replace('LIANA consensus with CellPhoneDB scoring · saved results · interactive figures', 'Exact five DESeq2 endpoints · ligand AND receptor DEGs · signed gene effects')
    text = re.sub(r'<label for="scope">.*?</select>', '<label for="scope">Fixed DEG endpoint scope</label><select id="scope" disabled><option value="IN_clusters">CA1 DL / SL ↔ inhibitory C1, C6, C8</option></select>', text, count=1)
    text = re.sub(r'<label for="tier">.*?<p id="tier-note" class="small"></p>', '<label for="tier">Candidate selection</label><select id="tier"><option value="DEG_only">DEGs only · all matching endpoint subunits</option><option value="replicated_specific">DEGs only + replicated/context-specific shortlist</option></select><p id="tier-note" class="small"></p>', text, count=1)
    text = text.replace('step="0.01" value="0.05"><label class="check"><input id="bothDEG" type="checkbox"> Gene-level DEG support at both ends</label>', 'step="0.01" value="1"><p class="small">1 includes every DEG-only candidate; 0.05 retains eight pooled records. This native rank is not genotype FDR.</p><label class="check"><input id="bothDEG" type="checkbox" checked disabled> Both ends must be DEGs · always enforced</label>')
    signed = '''<details open><summary>Gene up/down in cKO</summary>
<p>Directions use the original gene DESeq2 log2(cKO / WT), separately for the source ligand and target receptor. They do not indicate net interaction activity.</p>
<label for="ligandDirection">Ligand DEG direction</label><select id="ligandDirection"><option value="any">Any direction</option><option value="up">Up in cKO</option><option value="down">Down in cKO</option><option value="mixed">Mixed complex subunits</option></select>
<label for="receptorDirection">Receptor DEG direction</label><select id="receptorDirection"><option value="any">Any direction</option><option value="up">Up in cKO</option><option value="down">Down in cKO</option><option value="mixed">Mixed complex subunits</option></select>
<label for="directionPair">Ligand → receptor gene effects</label><select id="directionPair"><option value="any">All combinations</option><option value="up→up">Up → Up</option><option value="up→down">Up → Down</option><option value="down→up">Down → Up</option><option value="down→down">Down → Down</option><option value="mixed">Mixed subunit directions</option></select>
<label class="check"><input id="showDEGDirections" type="checkbox" checked> Show signed gene annotations</label>
<div class="triple"><div><label for="degUpColor">Up in cKO</label><input id="degUpColor" type="color" value="#c63d4f"></div><div><label for="degDownColor">Down in cKO</label><input id="degDownColor" type="color" value="#397eb5"></div><div><label for="degMixedColor">Mixed</label><input id="degMixedColor" type="color" value="#8b6bb1"></div></div>
<p>Dot rims: left = ligand, right = receptor. Gene-circle outlines and ↑/↓ labels show each gene’s sign. Numeric score colors and detection sizes retain their original meanings.</p></details>
'''
    text = text.replace('<details open><summary>Plot &amp; expression scale</summary>', signed + '<details open><summary>Plot &amp; expression scale</summary>', 1)
    text = text.replace('C1–C20', 'C1 / C6 / C8')
    a = text.index('<div class="notice">')
    b = text.index('<div id="workspace" hidden>', a)
    text = text[:a] + '''<div class="notice"><strong>Corrected: DEGs at both matching endpoints.</strong> Every ligand and receptor subunit must pass its own saved DESeq2 adjusted P &lt; 0.05. Only CA1 deep/superficial programs and inhibitory C1, C6, C8 are included. Up/down describes each gene’s cKO/WT RNA effect. The earlier one-endpoint rule is no longer available here. Native LIANA scores retain their original full-gene scoring background; no fits or ranks were recalculated.</div>
<div class="toolbar"><a href="RESULTS.html">Corrected findings &amp; exact endpoint audit</a><span>·</span><a href="data/deg_only/index.html">DEG-only tables &amp; up/down gene files</a><span>·</span><a href="https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/atlas-v3.6.0/SATB2_DEG_ONLY_LR_interactive_offline.html" download>Corrected offline HTML</a></div>
<div id="load-status" class="loading" role="status">Loading the three DEG-only projections…</div>''' + text[b:]
    a = text.index('<section id="findings-view"')
    b = text.index('<section id="downloads-view"', a)
    findings = '''<section id="findings-view" hidden><article class="card prose"><h2>Corrected results: ligand AND receptor DEGs</h2>
<p>The previous default accepted a DEG at either endpoint. It therefore included non-DEG partners and populations outside your five requested endpoints. This corrected view requires <b>every ligand AND receptor subunit to be a significant, unambiguous DEG in its exact source or target population</b>. Other inhibitory clusters and inhibitory total are excluded.</p>
<p>The pooled dataset contains <b>16 DEG-only population × ordered LR records</b>: six CA1→inhibitory and ten inhibitory→CA1. WT has 15 and cKO 16. With magnitude rank ≤0.05, these become pooled eight, WT nine and cKO seven. The optional replicated/context-specific shortlist contains two records. These rank thresholds describe native context priority, not genotype significance.</p>
<button id="show-lead" class="primary">Inspect C6 Cntn3 → Ptprg</button><button id="show-strong">Show the two context-specific DEG candidates</button>
<h3>Exact original DESeq2 endpoints</h3><div class="table-wrap"><table><thead><tr><th>Population / original endpoint</th><th>DEGs</th><th>Up in cKO</th><th>Down in cKO</th></tr></thead><tbody><tr><td>CA1 deep program · P0244</td><td>742</td><td>360</td><td>382</td></tr><tr><td>CA1 superficial program · P0247</td><td>1,825</td><td>928</td><td>897</td></tr><tr><td>Inhibitory C1 · P0140</td><td>27</td><td>11</td><td>16</td></tr><tr><td>Inhibitory C6 · P0145</td><td>22</td><td>4</td><td>18</td></tr><tr><td>Inhibitory C8 · P0147</td><td>1</td><td>1</td><td>0</td></tr></tbody></table></div>
<p>Deep and superficial refer to molecular expression programs, not measured anatomical layers. The complete endpoint names, hashes and unchanged seven-column DEG files are in the <a href="data/deg_only/index.html">corrected downloads</a>.</p>
<h3>Signed gene effects</h3><p><b>Cntn3 is down in inhibitory C6; Ptprg is up in both CA1 programs.</b> Their expression gate is supported by all six WT and six cKO mice. The two opposite gene signs do not determine net signaling activity.</p>
<div class="table-wrap"><table><thead><tr><th>Component &amp; population</th><th>log2(cKO / WT)</th><th>Adjusted P</th></tr></thead><tbody><tr><td>Cntn3 · inhibitory C6 · down</td><td>−0.53835</td><td>0.0037006</td></tr><tr><td>Ptprg · CA1 deep program · up</td><td>+0.39666</td><td>0.016489</td></tr><tr><td>Ptprg · CA1 superficial program · up</td><td>+0.49044</td><td>2.9271 × 10⁻⁸</td></tr></tbody></table></div>
<p>Every result exposes the source ligand and target receptor’s original gene log2FC, adjusted P and up/down direction. The same gene effects remain attached when viewing WT-only or cKO-only native scores. Mixed-sign complex subunits would be labeled mixed; all qualifying complexes in this dataset are single genes.</p>
<h3>Why C8 has no eligible interaction</h3><p>C8’s only DEG is <b>Calb1, up in cKO</b> (log2FC +1.93207; adjusted P 0.0154178). Calb1 is absent from mouseconsensus ligand/receptor subunits. C8 remains visible with zero qualifying edges and its mixed-identity review flag. This does not establish absence of biological communication.</p>
<h3>How this correction was calculated</h3><p>The five saved DESeq2 significant-gene lists were checked against their original stable-ID audits. Only matching endpoint symbols with unambiguous identities and adjusted P &lt;0.05 qualify. Every complex subunit must qualify at both ends. The 47 saved rows across pooled/WT/cKO were then selected from completed native LIANA results, preserving every native scalar score.</p>
<p>This is a <b>DEG-only projection of completed LIANA scores</b>. Full-gene normalization and the original population/background scoring context are retained. DESeq2, LIANA and its rank aggregation were not refitted on a truncated DEG matrix. Previous analyses and broad results remain preserved as historical files.</p>
<h3>How to read the native scores</h3><ul><li><b>magnitude_rank:</b> lower is higher consensus priority; not a genotype P value or FDR.</li><li><b>lr_logfc:</b> native source/target group-versus-rest specificity, not cKO/WT fold change.</li><li><b>lr_gmeans:</b> geometric mean of endpoint means on the native log-normalized expression scale.</li><li><b>CellPhoneDB P:</b> cell-label specificity from 1,000 permutations inside a scoring condition; not a differential genotype test. Saved zero is below empirical resolution.</li><li><b>Up/down:</b> original gene DESeq2 log2(cKO/WT), with each gene’s adjusted P &lt;0.05. No interaction-activity genotype test was fitted.</li></ul>
<h3>Completed methods and limits</h3><p>Frozen Seurat original_full counts/memberships were used for LIANA 2.0.0 with mouseconsensus and its native CellPhoneDB, Connectome, logFC, NATMI and SingleCellSignalR methods. Original log1p(CP10k) normalization uses the full source library. Every complex subunit must be detected in ≥10% of the scored group; minimum group size five; 1,000 cell-label permutations, seed 1337. Six independent WT and six cKO mice supply gene DESeq2 replication.</p>
<p>Nuclear RNA and curated interaction priors prioritize hypotheses; they do not measure protein availability, localization, contact, activation/inhibition or synapses. Mirrored Nrxn/Nlgn resource orientations are not independent physiological routes. WT/cKO native score or record-count differences are descriptive.</p><p><a href="RESULTS.html">All corrected findings and the complete 16-record pooled table</a> · <a href="https://liana.readthedocs.io/en/stable/">LIANA method documentation</a></p></article></section>
'''
    text = text[:a] + findings + text[b:]
    a = text.index('<section id="downloads-view"')
    b = text.index('</div></main></div>', a)
    text = text[:a] + '''<section id="downloads-view" hidden><div class="card"><h2>Corrected DEG-only results and figures</h2><p>Chart PDF, SVG and PNG buttons save the displayed filters, colors, sizes and signed gene annotations. Seven-column LR exports retain the requested native score columns; extended exports include both genes’ original DESeq2 log2FC, adjusted P and up/down direction.</p>
<div class="download-grid"><div><h3>Portable corrected atlas</h3><a href="https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/atlas-v3.6.0/SATB2_DEG_ONLY_LR_interactive_offline.html" download>DEG-only standalone HTML · works offline</a><a href="https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/atlas-v3.6.0/SATB2_DEG_ONLY_LR_results_20261023.zip">Corrected tables, signed gene files, figures and source · ZIP</a><a href="data/DEG_ONLY_MANIFEST.json" download>DEG-only manifest, rules and hashes</a><button id="download-settings-copy">Download this view’s settings</button></div>
<div><h3>Exact five endpoints</h3><a href="data/deg_only/index.html">All DEG-only tables, original DEGs and up/down subsets</a><a href="https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.6.0/communication_data/deg_only_v36/ENDPOINT_GENE_DIRECTIONS.tsv">Complete endpoint gene direction audit</a><a href="https://github.com/Ali-Mujahid/satb2_snRNA/tree/atlas-v3.6.0/analysis/communication">Native analysis source and configuration</a><a href="https://github.com/Ali-Mujahid/satb2_snRNA/blob/atlas-v3.6.0/scripts/audit_deg_only_communication.py">Reproduce the corrected DEG-only projection</a><a href="RESULTS.html">Corrected detailed findings</a></div></div>
<p>The standalone HTML embeds only the exact five populations and 47 DEG-only records across three conditions, with matching gene evidence. No server is needed for exploration or plot export. External source/table links require internet; the ZIP supplies corrected tables and figures locally.</p>
<details><summary>Historical broader results · superseded for the DEG-only question</summary><p>These archives include non-DEG partners under the earlier OR rule. They are preserved for provenance and must not be used as the corrected DEG-only result.</p><a href="https://github.com/Ali-Mujahid/satb2_snRNA/releases/tag/atlas-v3.5.1">Historical v3.5.1 explorer</a> · <a href="https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/atlas-v3.5.0/SATB2_CA1_IN_LIANA_20261022.zip">Historical full native analysis package</a> · <a href="data/full_result_tables/index.html">Historical complete native score tables</a></details></div></section>
<footer>Atlas v3.6.0 · corrected DEG-only projection, 23 October 2026. Original biological inputs, native LIANA scores and DESeq2 fits are unchanged. PDF plots reproduce the displayed figure as a high-resolution raster; SVG preserves vector geometry. <span id="provenance"></span></footer>
''' + text[b:]
    text = text.replace('data/data.js?', 'data/deg_only_data.js?')
    def version(match):
        src = match.group(1)
        path = PAGE.parent / src.split('?')[0]
        return 'src="' + src.split('?')[0] + '?v=' + hashlib.sha256(path.read_bytes()).hexdigest()[:12] + '"'
    text = re.sub(r'src="([^"]+\.js\?v=[^"]+)"', version, text)
    PAGE.write_text(text)


if __name__ == '__main__':
    migrate()
