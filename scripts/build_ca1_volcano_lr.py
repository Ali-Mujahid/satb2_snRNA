"""Join the existing reviewed LR rows to exact CA1 fit stable IDs; no rescoring."""
from pathlib import Path
import csv,gzip,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'analysis/synaptic_evidence_20261125/results/CA1_IN_directed_pair_roles.tsv'
OUT=ROOT/'analysis/ca1_volcano_lr_size_20261127'
OUT.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
with SRC.open() as f:original=list(csv.DictReader(f,delimiter='\t'))
assert len(original)==14
rows=[];summary={}
for key,layer in [('superficial','CA1_Superficial'),('deep','CA1_Deep')]:
    p=ROOT/f'analysis/ca1_layer_volcano_20261127/{key}.json.gz'
    fit=json.loads(gzip.decompress(p.read_bytes()))
    named={}
    for r in fit['rows']:named.setdefault(r[1],[]).append(dict(zip(fit['columns'],r)))
    chosen=[r for r in original if r['CA1_sublayer']==layer]
    for i,r in enumerate(chosen):
        assert len(named[r['CA1_gene']])==1
        gene=named[r['CA1_gene']][0]
        assert gene['padj']<.05 and float(r['IN_FDR'])<.05
        assert math.isclose(gene['log2FoldChange'],float(r['CA1_log2FC']),rel_tol=1e-7,abs_tol=1e-8)
        assert math.isclose(gene['padj'],float(r['CA1_FDR']),rel_tol=1e-7,abs_tol=1e-280)
        rows.append({**r,'row_id':key+'_'+str(i+1),'layer':key,'gene_id':gene['gene_id'],'CA1_log2FC_exact':gene['log2FoldChange'],'CA1_FDR_exact':gene['padj'],'resource_direction':'ca1_to_in' if r['source'].startswith('CA1_') else 'in_to_ca1','magnitude_rank':float(r['magnitude_rank'])})
    summary[key]={'rows':len(chosen),'CA1_genes':sorted(set(r['CA1_gene'] for r in chosen))}
assert summary['superficial']['rows']==9 and len(summary['superficial']['CA1_genes'])==7
assert summary['deep']['rows']==5 and len(summary['deep']['CA1_genes'])==5
data={'schema':'satb2-ca1-volcano-reviewed-lr-v1','source_sha256':sha(SRC),'source_path':str(SRC.relative_to(ROOT)),'criteria':'Both matching endpoints are DEGs at adjusted P <0.05 in the saved original-full 6WT/6cKO fits; superficial C2+C3, deep C1+C5. No new LIANA analysis.','summary':summary,'rows':rows}
(OUT/'lr_data.json').write_text(json.dumps(data,indent=2)+'\n')
(ROOT/'site/ca1_layer_volcano/lr_data.js').write_text('window.CA1_LAYER_LR='+json.dumps(data,separators=(',',':'))+';\n')
with (OUT/'CA1_IN_LR_highlight_mapping.tsv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
(OUT/'PROVENANCE.json').write_text(json.dumps({'models_changed':False,'source':str(SRC.relative_to(ROOT)),'source_sha256':sha(SRC),'gene_ID_mapping':'unique source symbol matched to exact fit gene_id; native fold changes/FDR agree with archived LR evidence within its numeric export precision','summary':summary},indent=2)+'\n')
print(json.dumps(summary))
