"""Package reviewed aggregate results for the existing atlas; never recount BAMs."""
from pathlib import Path
import base64, collections, csv, gzip, hashlib, json, re, shutil

BASE = Path(__file__).resolve().parent
REPO = BASE.parent / 'satb2_snRNA_github'
DEST = REPO / 'analysis/grm7_isoforms_20261126'
READER = DEST / 'reader'
shutil.copytree(BASE / 'deliverables', READER, dirs_exist_ok=True,
                ignore=shutil.ignore_patterns('plot_data.json', 'all_counts_long.tsv'))
archive = BASE / 'SATB2_Grm7_12mice_results_20261126.zip'
shutil.copy2(archive, DEST / archive.name)
h = (READER / 'index.html').read_text()
grm1_path = REPO / 'analysis/grm1_isoforms_20261126/reader/index.html'
h1 = grm1_path.read_text()
data_pattern = r'(<script id="data" type="application/json">)([\s\S]*?)(</script>)'
original_data = re.search(data_pattern, h)[2]
data = json.loads(original_data)
packed = json.dumps(dict(encoding='gzip-base64', payload=base64.b64encode(
    gzip.compress(original_data.encode(), mtime=0)).decode()), separators=(',', ':'))
h = re.sub(data_pattern, lambda m: m[1] + packed + m[3], h, count=1)
library = re.findall(r'<script[^>]*>[\s\S]*?</script>', h1)[1]
h = re.sub(r'<script>/\*\* @license[\s\S]*?</script>', lambda _: library, h, count=1)
h = h.replace("const D=JSON.parse(document.getElementById('data').textContent)",
    "const packed=JSON.parse(document.getElementById('data').textContent),D=packed.encoding==='gzip-base64'?JSON.parse(new TextDecoder().decode(fflate.gunzipSync(Uint8Array.from(atob(packed.payload),c=>c.charCodeAt(0))))):packed")
h = h.replace('<main>', '<main><p><a href="https://ali-mujahid.github.io/satb2_snRNA/">← SATB2 atlas</a> · <a href="https://ali-mujahid.github.io/satb2_snRNA/grm1_isoforms/index.html">Grm1 α/β evidence</a> · <a href="../' + archive.name + '" download>Download complete offline Grm7 results (ZIP)</a></p>', 1)
h = h.replace("IN_C1:'IN C1 · PV basket-like'", "IN_C2:'IN C2 · parent SST / O-LM-like',IN_C1:'IN C1 · PV basket-like'")
handler = """document.querySelectorAll('a[download]').forEach(a=>a.addEventListener('click',async e=>{if(new URL(a.href).origin===location.origin)return;e.preventDefault();try{const r=await fetch(a.href);if(!r.ok)throw Error('Download returned '+r.status);saveBlob(await r.blob(),new URL(a.href).pathname.split('/').pop());}catch(x){$('error').textContent=x.message;}}));"""
h = h.replace('scheme();draw();', handler + 'scheme();draw();')
(READER / 'index.html').write_text(h)
# Reader payload and all reviewed biological records remain identical.
assert json.loads(gzip.decompress(base64.b64decode(json.loads(re.search(data_pattern,h)[2])['payload']))) == data
for p in (BASE / 'deliverables').glob('*.tsv'):
    if (READER/p.name).exists():
        assert p.read_bytes() == (READER/p.name).read_bytes()

# Add a traceable interpretation table to Grm1 without modifying its counts.
G1 = BASE.parent / 'SATB2_GRM1_ISOFORMS_20261126/deliverables'
def read(name):
    return list(csv.DictReader((G1/name).open(), delimiter='\t'))
gene, family, junction = (read(n) for n in ['gene_UMI_counts_per_mouse.tsv',
    'alpha_beta_evidence_per_mouse.tsv','alpha_beta_specific_junction_counts_per_mouse.tsv'])
summary = []
for pop in ['IN_C2','IN_C6']:
    counts = {k:sum(int(r['count']) for r in gene if r['population']==pop and r['metric']==k)
              for k in ['gene','intronic_evidence','exonic_gene','spliced_gene']}
    counts['alpha_diagnostic'] = sum(int(r['count']) for r in family if r['population']==pop and r['metric']=='family_diagnostic_feature' and r['feature']=='alpha')
    counts['alpha_conservative'] = sum(int(r['count']) for r in family if r['population']==pop and r['metric']=='family_splice_supported' and r['feature']=='alpha')
    counts['alpha_terminal_junction'] = sum(int(r['count']) for r in junction if r['population']==pop and r['feature']=='junction_10689903_10719222')
    summary.append(dict(population=pop, **counts))
assert [r['alpha_conservative'] for r in summary] == [15,1]
assert [r['alpha_terminal_junction'] for r in summary] == [2,0]
table = grm1_path.parent / 'C2_C6_Grm1_evidence_explanation.tsv'
with table.open('w') as f:
    w=csv.DictWriter(f, fieldnames=list(summary[0]), delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
headings = [('gene','All Grm1 gene UMIs'),('intronic_evidence','UMIs with intronic evidence'),('exonic_gene','UMIs with exonic evidence'),('spliced_gene','UMIs with compatible anchored splice evidence'),('alpha_diagnostic','α-diagnostic feature UMIs'),('alpha_conservative','Conservative α-associated RNA UMIs'),('alpha_terminal_junction','Direct α-discriminating terminal-junction UMIs')]
rows=''.join('<tr><th style="text-align:left">'+label+'</th>'+''.join('<td style="padding:5px 22px;text-align:right">'+format(r[key],',')+'</td>' for r in summary)+'</tr>' for key,label in headings)
card='''<section class="card" id="c2-c6-interpretation"><h2>Why only 15 α-associated UMIs in C2 and one in C6?</h2><p>Gene-level Grm1 abundance and evidence for a processed α or β transcript answer different questions. Totals below sum the twelve mice. Most gene UMIs include intronic evidence; shared exons also cannot identify an isoform. These categories overlap and must not be added.</p><table><thead><tr><th>Evidence</th><th>IN C2</th><th>IN C6</th></tr></thead><tbody>'''+rows+'''</tbody></table><p><strong>The 15 and 1 are conservative RNA-support counts, not α abundance estimates or numbers of α-positive nuclei.</strong> Assignment requires family-specific feature evidence, compatibility with the α model, a compatible anchored splice and no incompatible anchored splice. The splice need not be the α-discriminating tail junction. Direct tail-junction support is only two in C2 and zero in C6. Zero β support and sparse α support do not establish absence, an α/β ratio or a genotype-dependent isoform switch.</p><p><strong>Validation priority:</strong> use α-skip and β-exon-inclusion junction RT-PCR on independently identified C2/C6-like cells, sequence the amplicons, then quantify with validated junction-specific RT-qPCR or RT-ddPCR across mice. Check assays against all six reference RNA models; a region coding only in α can be transcribed as β UTR and is not automatically α-specific. Include DNase, no-RT, no-template and cloned α/β specificity controls. Targeted long-read cDNA can resolve connected transcript structures. Confirm protein with validated C-terminal α/β antibodies and Grm1-knockout controls in identified cells.</p><p><a href="C2_C6_Grm1_evidence_explanation.tsv" download>Download evidence table</a> · <a href="https://pubmed.ncbi.nlm.nih.gov/9779943/">Ferraguti et al., 1998: isoform-specific protein localization</a> · <a href="https://www.frontiersin.org/journals/molecular-neuroscience/articles/10.3389/fnmol.2018.00439/full">Naito et al., 2018: splice-variant RNA and protein assays</a></p></section>'''
if 'id="c2-c6-interpretation"' not in h1:
    h1=h1.replace('<section class="card"><h2>Per-animal evidence</h2>',card+'\n<section class="card"><h2>Per-animal evidence</h2>',1)
grm1_path.write_text(h1)

# Pages is already near its size ceiling; use a tiny loader and external aggregate assets.
page=REPO/'site/grm7_isoforms/index.html';page.parent.mkdir(exist_ok=True)
page.write_text((REPO/'site/grm1_isoforms/index.html').read_text().replace('Grm1','Grm7').replace('grm1','grm7').replace('alpha/beta','A/B'))
main=REPO/'site/index.html';text=main.read_text()
link='<a href="grm7_isoforms/index.html" class="result-link">Grm7 A/B &amp; junction evidence</a>'
assert link not in text
text=text.replace('<a href="grm1_isoforms/index.html"',link+'<a href="grm1_isoforms/index.html"',1);main.write_text(text)
size=sum(p.stat().st_size for p in (REPO/'site').rglob('*') if p.is_file())
assert size <= 950_000_000, size
receipt=dict(status='PACKAGED',pages_bytes=size,original_count_tables_unchanged=True,
    archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
    embedded_data_compression_verified=True,C2_C6_Grm1_evidence=summary)
(BASE/'PACKAGING_REVIEW.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt,indent=2))
