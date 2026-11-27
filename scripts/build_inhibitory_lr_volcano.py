"""Expose existing original IN C1/C6 fits and the inhibitory ends of reviewed LR rows."""
from pathlib import Path
import csv,gzip,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT.parent
OUT=ROOT/'analysis/inhibitory_lr_volcano_20261127'
SITE=ROOT/'site/ca1_layer_volcano'
TAG='atlas-v3.22-inhibitory-lr-volcano-20261127'
AUDIT=WORK/'SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0/results/deseq2/original_full/seurat'
OUT.mkdir(exist_ok=True,parents=True)
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(b):return hashlib.sha256(b).hexdigest()
def num(v):
    try:x=float(v)
    except (ValueError,TypeError):return None
    return x if math.isfinite(x) else None
manifest=json.loads((ROOT/'analysis/ca1_layer_volcano_20261127/PROVENANCE.json').read_text())
manifest['extension_tag']=TAG
columns=['gene_id','gene_symbol','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']
fits={}
for key,cluster,pid,expected in [('in_c1',1,140,(27,11,16,12769)),('in_c6',6,145,(22,4,18,8989))]:
    directory=AUDIT/f'P{pid:04d}_seurat_Inhibitory_CLUSTER/PRIMARY/AUDIT'
    source=directory/'ALL_RESULTS_STABLE_ID_AUDIT.tsv'
    coverage=directory/'MICE_ALL12_AUDIT.tsv'
    rows=[]
    for r in read(source):
        symbol=r['source_gene_symbol']
        rows.append([r['gene_id'],'' if symbol in {'NA','NaN'} else symbol,*[num(r[c]) for c in columns[2:]]])
    assert len({r[0] for r in rows})==len(rows)
    sig=[r for r in rows if r[7] is not None and r[7]<.05]
    count=(len(sig),sum(r[3]>0 for r in sig),sum(r[3]<0 for r in sig),len(rows))
    assert count==expected,(key,count)
    mice=read(coverage)
    assert len(mice)==12 and all(r['included_in_model']=='TRUE' for r in mice)
    assert sum(r['condition']=='WT' for r in mice)==sum(r['condition']=='cKO' for r in mice)==6
    metadata={'id':key,'label':f'Inhibitory C{cluster} · '+('PV reference' if cluster==1 else 'original SST cluster'),'compartment':'Inhibitory','clusters':[cluster],'n_WT':6,'n_cKO':6,'n_nuclei':sum(int(r['n_nuclei']) for r in mice),'n_rows':len(rows),'DEGs':count[0],'up':count[1],'down':count[2],'source_path':str(source),'source_sha256':sha(source.read_bytes()),'coverage_sha256':sha(coverage.read_bytes()),'fit_date':'archived original full-cohort fit','export_stem':f'Seurat_IN_C{cluster}_original_full_volcano','annotation_note':'Original parent cluster; the C6.0/C6.1/C6.2 subdivision and mouse-subset fits are not substituted.'}
    fits[key]={'metadata':metadata,'columns':columns,'rows':rows,'mice':mice}
    raw=json.dumps(fits[key],separators=(',',':'),allow_nan=False).encode();packed=gzip.compress(raw,mtime=0)
    (OUT/(key+'.json.gz')).write_bytes(packed)
    (OUT/(key+'_source.tsv.gz')).write_bytes(gzip.compress(source.read_bytes(),mtime=0))
    (OUT/(key+'_mouse_audit.tsv')).write_bytes(coverage.read_bytes())
    (OUT/(key+'_model_settings.tsv')).write_bytes((directory/'MODEL_SETTINGS.tsv').read_bytes())
    manifest['fits'][key]={**metadata,'url':f'https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/{TAG}/analysis/inhibitory_lr_volcano_20261127/{key}.json.gz','bytes':len(packed),'sha256':sha(packed),'payload_sha256':sha(raw)}
    print(key,count,'nuclei',metadata['n_nuclei'])

original=json.loads((ROOT/'analysis/ca1_volcano_lr_size_20261127/lr_data.json').read_text())
assert len(original['rows'])==14
views=[]
for r in original['rows']:
    key=r['IN_cluster'].lower()
    fit=fits[key]
    match=[dict(zip(columns,v)) for v in fit['rows'] if v[1]==r['IN_gene']]
    assert len(match)==1
    g=match[0]
    assert g['padj']<.05
    assert math.isclose(g['log2FoldChange'],float(r['IN_log2FC']),rel_tol=1e-7,abs_tol=1e-8)
    assert math.isclose(g['padj'],float(r['IN_FDR']),rel_tol=1e-7,abs_tol=1e-280)
    common={**r,'CA1_gene_id':r['gene_id'],'IN_gene_id':g['gene_id'],'IN_log2FC_exact':g['log2FoldChange'],'IN_FDR_exact':g['padj']}
    views.append({**common,'view_compartment':'CA1','plot_gene':r['CA1_gene'],'plot_role':r['CA1_role'],'plot_log2FC_exact':r['CA1_log2FC_exact'],'plot_FDR_exact':r['CA1_FDR_exact'],'partner_population':r['IN_cluster'],'partner_gene':r['IN_gene'],'partner_log2FC':g['log2FoldChange']})
    views.append({**common,'row_id':key+'_'+r['row_id'],'layer':key,'gene_id':g['gene_id'],'view_compartment':'Inhibitory','plot_gene':r['IN_gene'],'plot_role':r['IN_role'],'plot_log2FC_exact':g['log2FoldChange'],'plot_FDR_exact':g['padj'],'partner_population':r['CA1_sublayer'],'partner_gene':r['CA1_gene'],'partner_log2FC':r['CA1_log2FC_exact']})
summary={key:{'rows':len([r for r in views if r['layer']==key]),'plot_genes':sorted({r['plot_gene'] for r in views if r['layer']==key})} for key in manifest['fits']}
assert summary['in_c1']=={'rows':7,'plot_genes':['Fstl5','Grm7','Nlgn1','Ptprz1']}
assert summary['in_c6']=={'rows':7,'plot_genes':['Cntn3','Robo1','Slit2']}
lr={**{k:v for k,v in original.items() if k not in {'rows','summary'}},'schema':'satb2-ca1-in-volcano-reviewed-lr-v2','original_population_pair_rows':14,'view_rows':28,'summary':summary,'rows':views}
(OUT/'LR_VIEW_DATA.json').write_text(json.dumps(lr,indent=2)+'\n')
with (OUT/'IN_CA1_LR_endpoint_mapping.tsv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(views[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(views)
(OUT/'PROVENANCE.json').write_text(json.dumps(manifest,indent=2)+'\n')
(SITE/'manifest.js').write_text('window.CA1_LAYER_VOLCANO_MANIFEST='+json.dumps(manifest,separators=(',',':'))+';\n')
(SITE/'lr_data.js').write_text('window.CA1_LAYER_LR='+json.dumps(lr,separators=(',',':'))+';\n')
print(json.dumps(summary))
