from pathlib import Path
import pandas as pd, json, math
src=Path('/Users/mujahid357/Documents/codex_new/SATB2_INHIBITORY_AUDIT_20261111')
root=src/'runs/audit_v2/results/results'; dest=Path('site/inhibitory_sensitivity');dest.mkdir(parents=True,exist_ok=True)
status=pd.read_csv(root/'model_status.tsv',sep='\t'); samples=pd.read_csv(src/'figures/normalized_counts_both_models.tsv',sep='\t')
records={};summ=[]
for label in ['TOTAL']+[f'C{i}' for i in range(1,21)]:
 files=[root/label/x/'all_results_with_stable_ids.tsv' for x in ['original','main_partition']]
 if not all(x.exists() for x in files):summ.append(dict(population=label,status='NOT_TESTED'));continue
 a,b=[pd.read_csv(x,sep='\t').set_index('gene_id') for x in files];u=a.index.union(b.index);a=a.reindex(u);b=b.reindex(u)
 ag=a.padj.lt(.05);bg=b.padj.lt(.05);sel=ag|bg; data=[]
 for geneid in u[sel]:
  x=a.loc[geneid];y=b.loc[geneid];f=bool(x.padj<.05);g=bool(y.padj<.05)
  category='shared' if f and g else ('full_only' if f else 'restricted_only')
  def val(v):return None if pd.isna(v) or not math.isfinite(float(v)) else float(v)
  data.append(dict(id=geneid,gene=str(x.gene) if pd.notna(x.gene) else (str(y.gene) if pd.notna(y.gene) else ''),category=category,full={col:val(x[col]) for col in ['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']},restricted={col:val(y[col]) for col in ['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']}))
 records[label]=data
 counts={k:sum(x['category']==k for x in data) for k in ['shared','full_only','restricted_only']}
 missing_fullonly=sum(x['category']=='full_only' and x['restricted']['padj'] is None for x in data);missing_new=sum(x['category']=='restricted_only' and x['full']['padj'] is None for x in data)
 sm=status[(status.population==label)&(status.model.isin(['original','main_partition']))];summ.append(dict(population=label,status='FITTED',full=int(ag.sum()),restricted=int(bg.sum()),**counts,full_other_na=missing_fullonly,restricted_other_na=missing_new,full_mice=[int(sm[sm.model=='original'].WT.iloc[0]),int(sm[sm.model=='original'].cKO.iloc[0])],restricted_mice=[int(sm[sm.model=='main_partition'].WT.iloc[0]),int(sm[sm.model=='main_partition'].cKO.iloc[0])]))
assert len(summ)==21 and next(x for x in summ if x['population']=='C17')['restricted_only']==2
p={'summary':summ,'genes':records,'samples':samples.to_dict('records')}
(dest/'data.js').write_text('window.SENSITIVITY_DATA='+json.dumps(p,separators=(',',':'),allow_nan=False)+';\n')
pd.DataFrame(summ).to_csv(src/'verification/deg_overlap_by_population.tsv',sep='\t',index=False)
print('Wrote',dest/'data.js','bytes', (dest/'data.js').stat().st_size)
