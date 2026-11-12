from pathlib import Path
import pandas as pd,json,math
r=Path('/Users/mujahid357/Documents/codex_new/SATB2_INHIBITORY_SUBSET_SENSITIVITY_20261112/reviewed_results');orig=Path('/Users/mujahid357/Documents/codex_new/SATB2_INHIBITORY_AUDIT_20261111/runs/audit_v2/results/results');out=Path('analysis/inhibitory_mouse_subsets_20261112');cols=['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'];s=pd.read_csv(r/'all_planned_scenarios.tsv',sep='\t');g=pd.read_csv(r/'all_significant_gene_calls.tsv.gz',sep='\t');summary={};hits={};baseline={}
def clean(x):return None if pd.isna(x) else (int(x) if isinstance(x,int) else float(x) if isinstance(x,float) else x)
for _,x in s.iterrows():
 d={k:clean(x[k]) for k in ['status','reason','selected_WT','selected_cKO','excluded_WT','excluded_cKO','n_WT','n_cKO','min_nuclei','n_genes','n_FDR_available','n_deg','n_up','n_down','n_shared_original','n_gained','n_lost','n_gained_original_FDR_na','n_lost_subset_FDR_na']};summary.setdefault(x.population,{})[x.scenario]=d
for (pop,scenario),x in g.groupby(['population','scenario'],sort=False):
 rows=[]
 for _,v in x.iterrows():rows.append([str(v.gene_id),'' if pd.isna(v.gene) else str(v.gene),*[clean(v[c]) for c in cols],clean(v.original_padj),clean(v.original_log2FoldChange)])
 hits.setdefault(pop,{})[scenario]=rows
for pop in ['TOTAL']+[f'C{x}' for x in range(1,21)]:
 f=orig/pop/'original/all_results_with_stable_ids.tsv'
 if f.exists():
  d=pd.read_csv(f,sep='\t');d=d[d.padj.lt(.05)];baseline[pop]=[[str(x.gene_id),'' if pd.isna(x.gene) else str(x.gene),*[clean(x[c]) for c in cols]] for _,x in d.iterrows()]
 else:baseline[pop]=None
for pop,scenarios in summary.items():
 for sc,x in scenarios.items():
  if x['status']=='FITTED':
   hh=hits.get(pop,{}).get(sc,[]);assert len(hh)==x['n_deg'],(pop,sc,len(hh),x['n_deg']);bb=baseline[pop]
   if bb is not None:
    ai={z[0] for z in bb};bi={z[0] for z in hh};assert len(ai&bi)==x['n_shared_original'] and len(bi-ai)==x['n_gained'] and len(ai-bi)==x['n_lost'],(pop,sc)
payload={'version':'20261112-v1','stats':['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'],'summary':summary,'hits':hits,'baseline':baseline,'global':json.loads((r/'verification.json').read_text())}
(out/'subset_payload_v1.json').write_text(json.dumps(payload,separators=(',',':'),allow_nan=False))
print('Payload bytes', (out/'subset_payload_v1.json').stat().st_size,'models',sum(len(x) for x in summary.values()),'gene records',len(g))
