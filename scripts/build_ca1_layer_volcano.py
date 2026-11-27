"""Package existing Seurat C2+C3 and C1+C5 fits; do not fit statistical models."""
from pathlib import Path
import csv, gzip, hashlib, json, math

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT.parent
OUT = ROOT / 'analysis/ca1_layer_volcano_20261127'
SITE = ROOT / 'site/ca1_layer_volcano'
TAG = 'atlas-v3.20-ca1-layer-volcano-20261127'
OUT.mkdir(exist_ok=True, parents=True)
SITE.mkdir(exist_ok=True, parents=True)
def sha(b): return hashlib.sha256(b).hexdigest()
def read(p):
    with p.open() as f: return list(csv.DictReader(f, delimiter='\t'))
def numeric(v):
    try: x = float(v)
    except (ValueError, TypeError): return None
    return x if math.isfinite(x) else None

sources = {
    'superficial': (WORK/'SATB2_CA1_EPHYS_PLANNING_20261125/superficial_C2_C3_DESeq2.tsv', WORK/'SATB2_CA1_EPHYS_PLANNING_20261125/superficial_mouse_coverage.tsv', [2, 3], (1775,938,837)),
    'deep': (WORK/'SATB2_CA1_IN_EXPANDED_20261102/pseudobulk/audit/CA1_deep_C1_C5/ALL_RESULTS_STABLE_ID_AUDIT.tsv', WORK/'SATB2_CA1_IN_EXPANDED_20261102/pseudobulk/audit/CA1_deep_C1_C5/ALL12_POPULATION_COVERAGE.tsv', [1, 5], (742,360,382)),
}
fits = {}
columns = ['gene_id','gene_symbol','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']
manifest = {'schema':'satb2-ca1-layer-volcano-v1','contrast':'cKO / WT','tool':'Seurat','cohort':'original_full','models_refitted':False,'tag':TAG,'fits':{}}
for key,(path,coverage,clusters,expected) in sources.items():
    original = read(path)
    rows=[]
    for r in original:
        symbol=r.get('gene_symbol',r.get('source_gene_symbol','')) or ''
        if symbol in {'NA','NaN'}: symbol=''
        rows.append([r['gene_id'],symbol,*[numeric(r.get(c)) for c in columns[2:]]])
    assert len(set(r[0] for r in rows)) == len(rows)
    sig=[r for r in rows if r[7] is not None and r[7]<.05]
    actual=(len(sig),sum(r[3]>0 for r in sig),sum(r[3]<0 for r in sig))
    assert actual==expected,(key,actual)
    mice=read(coverage)
    assert len(mice)==12 and sum(m['condition']=='WT' for m in mice)==6 and sum(m['condition']=='cKO' for m in mice)==6
    meta={'id':key,'label':('Superficial' if key=='superficial' else 'Deep')+' · '+' + '.join('C'+str(k) for k in clusters),'clusters':clusters,'n_WT':6,'n_cKO':6,'n_nuclei':sum(int(m['n_nuclei']) for m in mice),'n_rows':len(rows),'DEGs':actual[0],'up':actual[1],'down':actual[2],'source_path':str(path),'source_sha256':sha(path.read_bytes()),'coverage_sha256':sha(coverage.read_bytes()),'fit_date':'2026-11-25' if key=='superficial' else 'archived original fit','C8_included':False}
    fits[key]={'metadata':meta,'columns':columns,'rows':rows,'mice':mice}
    raw=json.dumps(fits[key],separators=(',',':'),allow_nan=False).encode()
    packed=gzip.compress(raw,mtime=0)
    name=key+'.json.gz'
    (OUT/name).write_bytes(packed)
    manifest['fits'][key]={**meta,'url':f'https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/{TAG}/analysis/ca1_layer_volcano_20261127/{name}','bytes':len(packed),'sha256':sha(packed),'payload_sha256':sha(raw)}
    # Preserve complete native numeric text, including missing values, for audit/rerun use.
    (OUT/(key+'_source.tsv.gz')).write_bytes(gzip.compress(path.read_bytes(),mtime=0))
    (OUT/(key+'_mouse_coverage.tsv')).write_bytes(coverage.read_bytes())
    print(key,len(rows),'genes;',actual,'DEGs/up/down;',len(packed),'compressed bytes')
(OUT/'PROVENANCE.json').write_text(json.dumps(manifest,indent=2)+'\n')
(SITE/'manifest.js').write_text('window.CA1_LAYER_VOLCANO_MANIFEST='+json.dumps(manifest,separators=(',',':'))+';\n')
