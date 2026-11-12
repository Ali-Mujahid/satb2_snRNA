from pathlib import Path
import pandas as pd,numpy as np,json,math,zipfile,hashlib
root=Path('/Users/mujahid357/Documents/codex_new/SATB2_INHIBITORY_SELECTED_VOLCANO_20261112');r=root/'runs/selected_v1/results/results';pkg=Path('/Users/mujahid357/Documents/codex_new/SATB2_SEURAT_INHIBITORY_DESEQ2_PACKAGE_20261109');out=Path('analysis/inhibitory_6wt3cko_volcano_20261112');data=out/'data';pops=['TOTAL']+[f'C{i}' for i in range(1,21)];ix=pd.read_csv(pkg/'04_population_indexes/INHIBITORY_POPULATION_INDEX.tsv',sep='\t');summ=[];mouse=[];test=[]
def native(v):
 if isinstance(v,np.generic):v=v.item()
 if pd.isna(v):return None
 if isinstance(v,(float,int)):return float(v) if math.isfinite(float(v)) else None
 return v
for pop in pops:
 row=pd.read_csv(r/pop/'COMPLETE.tsv',sep='\t').iloc[0].to_dict();cover=pd.read_csv(r/pop/'mouse_coverage_all12.tsv',sep='\t');cover.insert(0,'population',pop);mouse.append(cover)
 iv=ix[ix.analysis_kind=='total'] if pop=='TOTAL' else ix[(ix.analysis_kind=='cluster')&(ix.analysis_label.astype(str)==pop[1:])];assert len(iv)==1;f=pkg/'06_primary_results_and_audits'/iv.iloc[0].population_key/'AUDIT/ALL_RESULTS_STABLE_ID_AUDIT.tsv';orig=pd.read_csv(f,sep='\t').set_index('gene_id') if f.exists() else None
 for k in ['full_deg','full_up','full_down','overlap_total','overlap_up','overlap_down','overlap_opposite','selected_only','full_only','selected_only_full_fdr_na','full_only_selected_fdr_na']:row[k]=None
 if orig is not None:
  old=orig[orig.padj.lt(.05)];row.update(full_deg=len(old),full_up=int(old.log2FoldChange.gt(0).sum()),full_down=int(old.log2FoldChange.lt(0).sum()))
 if row['status']=='FITTED':
  x=pd.read_csv(r/pop/'selected_all_results.tsv.gz',sep='\t').set_index('gene_id');assert x.index.is_unique and orig is not None;new=x[x.padj.lt(.05)];ids=set(new.index)&set(old.index);newonly=set(new.index)-set(old.index);oldonly=set(old.index)-set(new.index)
  sameup=sum(new.loc[i,'log2FoldChange']>0 and old.loc[i,'log2FoldChange']>0 for i in ids);samedown=sum(new.loc[i,'log2FoldChange']<0 and old.loc[i,'log2FoldChange']<0 for i in ids);opp=len(ids)-sameup-samedown
  row.update(overlap_total=len(ids),overlap_up=sameup,overlap_down=samedown,overlap_opposite=opp,selected_only=len(newonly),full_only=len(oldonly),selected_only_full_fdr_na=int(orig.reindex(list(newonly)).padj.isna().sum()),full_only_selected_fdr_na=int(x.reindex(list(oldonly)).padj.isna().sum()))
  joined=x.join(orig[['log2FoldChange','padj']],rsuffix='_full',how='left');cols=['gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj','log2FoldChange_full','padj_full'];rows=[]
  for gid,v in joined[cols].iterrows():rows.append([gid,*[native(v[c]) for c in cols]])
  missing=[]
  for gid,v in old.loc[list(set(old.index)-set(x.index))].iterrows():missing.append([gid,native(v.get('display_symbol')),None,None,None,None,None,None,native(v['log2FoldChange']),native(v['padj'])])
  d={'population':pop,'fields':['gene_id',*cols],'rows':rows,'full_only_missing':missing};(data/f'{pop}.json').write_text(json.dumps(d,separators=(',',':'),allow_nan=False));test.append({'population':pop,'genes':len(x),'deg':len(new),'overlap':len(ids)})
 summ.append({str(k):native(v) for k,v in row.items()})
assert sum(x['n_deg'] or 0 for x in summ if x['population']!='TOTAL')==198
assert next(x for x in summ if x['population']=='TOTAL')['n_deg']==60
assert {x['population'] for x in summ if x['status']=='NOT_TESTED'}=={'C18','C20'}
cover=pd.concat(mouse,ignore_index=True);cover.to_csv(out/'nuclei_per_mouse_all12.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(summ).to_csv(out/'cluster_summary.tsv',sep='\t',index=False,na_rep='NA')
(out/'summary.json').write_text(json.dumps({'version':'20261112-v1','summary':summ,'mouse':cover.where(pd.notna(cover),None).to_dict('records')},separators=(',',':'),allow_nan=False,default=native))
with zipfile.ZipFile(out/'all_cluster_DESeq2_results.zip','w',zipfile.ZIP_DEFLATED) as z:
 z.write(out/'cluster_summary.tsv','cluster_summary.tsv');z.write(out/'nuclei_per_mouse_all12.tsv','nuclei_per_mouse_all12.tsv')
 for pop in pops:
  for name in ['COMPLETE.tsv','selected_standard_7_columns_all_genes.tsv.gz','selected_standard_7_columns_DEGs.tsv','selected_all_results.tsv.gz']:
   f=r/pop/name
   if f.exists():z.write(f,f'{pop}/{name}')
(out/'data_validation.json').write_text(json.dumps({'status':'PASS','cluster_gene_calls':198,'total_degs':60,'not_tested':['C18','C20'],'matched_saved_search':True,'populations':test},indent=2))
print('PASS',len(summ),'populations;',sum(x['n_deg'] or 0 for x in summ if x['population']!='TOTAL'),'cluster gene calls; bytes',sum(f.stat().st_size for f in data.glob('*.json')))
