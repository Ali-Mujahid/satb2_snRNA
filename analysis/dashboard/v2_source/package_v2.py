"""Package a new portable viewer from unchanged native exports and v2 frontend."""
from pathlib import Path
import argparse, datetime, hashlib, json, re, shutil, subprocess, zipfile

BASE = Path(__file__).resolve().parent.parent
SRC = BASE / 'offline_dashboard_v2_source'
OUT = BASE / 'reports/offline_dashboard_v2'

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''): h.update(b)
    return h.hexdigest()

def clone_tree(source, target):
    for p in source.rglob('*'):
        q = target / p.relative_to(source)
        if p.is_dir(): q.mkdir(parents=True, exist_ok=True); continue
        q.parent.mkdir(parents=True, exist_ok=True)
        if q.exists() and sha(q) == sha(p): continue
        if q.exists(): q.unlink()
        r = subprocess.run(['cp', '-c', str(p), str(q)], capture_output=True)
        if r.returncode: shutil.copy2(p, q)

README = '''# SATB2 offline explorer v2

Open index.html in the complete unzipped folder. The saved new-control examples
are in assets/v2_plot_examples/index.html. Keep its data, vendor, assets
and saved_results folders together. No internet, installation or server is
needed. This is a separate delivery; the previous analysis and v1 stay intact.

The green Saved DEG lists & figure books link opens the full searchable result
index. It contains per-total, per-cluster and reviewed-category source tables,
948 PRIMARY endpoint statuses (852 tested, 96 unsupported), 112 sensitivity
endpoints, full mouse coverage, and native UMAP / representative marker /
provided-panel figures. Gene-symbol DESeq2 tables keep seven standard fields:
gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj. Stable IDs are separate
audits. FDR05 files use within-endpoint BH <0.05 with no extra effect cutoff.

Use category shortcuts for all nuclei, excitatory, CA1/ProS candidates, CA1
sublayer-associated expression programs, and inhibitory. Pure anatomical CA1
and physical superficial/deep layers were not established by these data.

Marker expression: all 878 already-produced marker genes are selectable.
Search, paste gene symbols, select matches or use representative / supplied
marker presets. Choose 60, 90, 120 or 180 genes per plot page. Wide plots scroll
horizontally; PDF exports save every selected gene page at a matching page
size. Mean normalized expression, saved gene-scaled z, detected fraction and
mean UMIs are available. Set source-default, automatic or manual color limits;
limits clip color only, not the saved values. Adjust dot area and negative-z
attenuation. Negative z is smaller/fainter by default. Missing produced z,
including cKO-only z, stays unavailable; the UI offers mean expression instead.

UMAPs: select Harmony or unintegrated saved coordinates and reviewed / edited,
MapMyCells, program, genotype, mouse, QC or gene coloring. The gene example
button chooses an available representative gene. All 878 gene overlays use
cluster summaries, not individual-nucleus RNA measurements. Choose all/WT/cKO
summary values independently of which nuclei are displayed. UMAP-only genotype
and CA1-program filters do not change dotplot cohorts or saved DE memberships.

Inhibitory spacing is a display-only centroid displacement. At 1x the original
coordinates are exact. Other settings preserve within-cluster offsets but
alter between-cluster distances and can create overlap. That display has no
biological distance interpretation; PDF metadata and the plot flag it.

Labels & exclusions: type a label, choose an evidence-backed provisional
suggestion or a source reference/program label, and adjust colors. Export
annotation CSV and settings JSON. Source labels, confidence and nucleus
membership remain preserved. Seventy original clusters fail the native fine
identity criterion (61 retained, nine excluded audit-only); all 70 keep that
uncertainty after suggestions. Descriptive '-like' labels are not validated
identities. CA1 expression-program options keep physical layer undetermined.

Mouse boxplots: search all 32,285 source genes, including genes outside the
878-gene marker panel. Select total, cluster or saved category per tool. One
point is one independent mouse, not one nucleus. Raw counts use the completed
all-gene integer mouse arrays. Normalized values divide by the exact saved
endpoint's recorded size factor; log1p is optional. Missing/non-model factors
remain NA. Raw mode can show non-model mice with an explicit option. Boxes show
median and IQR; whiskers use observed values within 1.5xIQR. No new P values are
calculated. CSV exports include all selected genes/mice and source provenance.
Full-gene display does not imply the gene was tested in a DESeq2 model.

Saved DESeq2 fits: select total, original cluster, saved category or a custom
label union. A fit is reused only for exact saved nucleus membership, count
branch and endpoint. An edited total or union without a matching fit displays
NO MATCHING FIT; it does not borrow results or normalize with another fit.
Unchanged clusters retain their own fits when other groups are excluded.
Export requests document selections but cannot run DESeq2 in an offline page.

All/EX baseline exclusions follow the prior shared <100-nucleus rule. IN
pseudobulk uses every available nucleus, with at least three independent WT
and three cKO mice; fewer-than-five-nucleus contributions stay exploratory.
Fine label confidence is separate from model eligibility. Both tools' models
use one R DESeq2 engine, so identical core effects are reproducibility of this
shared calculation, not independent statistical confirmation.

Validation covers source identities, all-gene mouse count totals, recorded fit
factors, exact membership guards, new controls, simulated DOM drawing and
actual PDF exports. Live browser/CSS rendering was not verified: automation
previously rejected local-file navigation and that restriction was respected.
See DASHBOARD_REVIEW.md and BUNDLE_MANIFEST.json for evidence and integrity.
'''

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--zip', action='store_true'); args=ap.parse_args()
    OUT.mkdir(exist_ok=True)
    for name in ['index.html','app.js','styles.css','v2_features.js']:
        shutil.copy2(SRC / name, OUT / name)
    clone_tree(SRC / 'vendor', OUT / 'vendor')
    complete = BASE / 'reports/saved_results_v2/COMPLETE.json'
    assert json.loads(complete.read_text())['status'] == 'COMPLETE'
    clone_tree(BASE / 'reports/saved_results_v2', OUT / 'saved_results')
    ann = BASE / 'reviews/annotation_suggestions_v2'
    data = json.loads((ann / 'editor_label_suggestions.json').read_text())
    (OUT / 'data/annotation_suggestions.js').write_text('window.SATB2_ANNOTATION_SUGGESTIONS='+json.dumps(data, separators=(',',':'), ensure_ascii=False)+';\n')
    shutil.copy2(ann / 'review.md', OUT / 'data/annotation_review.md')
    shutil.copy2(ann / 'cluster_label_suggestions.tsv', OUT / 'data/cluster_label_suggestions.tsv')
    shutil.copy2(ann / 'REVIEW_PROVENANCE.json', OUT / 'data/ANNOTATION_REVIEW_PROVENANCE.json')
    (OUT / 'README.md').write_text(README)
    for name in ['OFFLINE_DELIVERY_COMPLETE.json','DE_CATALOG_COMPLETE.json']:
        # Catalog completion is still source-bound; prior delivery receipts must
        # not imply they cover this revised frontend.
        pass
    html=(OUT/'index.html').read_text()
    for ref in re.findall(r'(?:src|href)="([^"]+)"',html):
        if ref.startswith(('#','data:')): continue
        assert not re.match(r'^(?:https?:)?//',ref),ref
        assert (OUT/ref).is_file(),ref
    files=[p for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='BUNDLE_MANIFEST.json']
    receipt={'schema':'satb2-offline-bundle-2','generated_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'entrypoint':'index.html','file_count':len(files),'bytes':sum(p.stat().st_size for p in files),'scientific_fitting_performed':False,'raw_nucleus_count_aggregation_performed':False,'previous_analysis_files_modified':False,'files':[{'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files]}
    (OUT/'BUNDLE_MANIFEST.json').write_text(json.dumps(receipt,indent=2)+'\n')
    if args.zip:
        z=OUT.with_suffix('.zip')
        with zipfile.ZipFile(z,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as f:
            for p in sorted(OUT.rglob('*')):
                if p.is_file(): f.write(p,Path(OUT.name)/p.relative_to(OUT))
        print(json.dumps({'zip':str(z),'bytes':z.stat().st_size,'sha256':sha(z)}))
    print(json.dumps({'folder':str(OUT),'files':len(files),'bytes':receipt['bytes']}))

if __name__=='__main__': main()
