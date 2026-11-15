#!/usr/bin/env python3
"""Recover marker summaries from same original UMI layer, grouped by frozen Seurat nuclei."""
import csv,json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
root=Path('/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013')
h5=root/'runs/20261013_v1_0_5/results/scanpy/original_full/Inhibitory/scanpy.h5ad'
meta=Path('/Users/mujahid357/Documents/codex_new/SATB2_SEURAT_INHIBITORY_DESEQ2_PACKAGE_20261109/01_nucleus_metadata/rendered_nucleus_metadata.tsv')
out=Path(__file__).with_name('data.json')
A=ad.read_h5ad(h5)
M=pd.read_csv(meta,sep='\t',dtype={'cluster_id':str}).set_index('nucleus_id',verify_integrity=True)
assert A.n_obs==8911 and len(M)==8911 and A.obs_names.is_unique and set(A.obs_names)==set(M.index)
M=M.loc[A.obs_names]
assert (A.obs['nCount_RNA'].astype(float).to_numpy()==M['original_nCount_RNA'].astype(float).to_numpy()).all()
genes=['Elfn1','Chodl','Chrna2','Ctnnd2','Htr3a','Cdh7','Sst','Reln','Trpc6','Pcdh11x','Grin3a']; ids=['2','6','10','17']
data=json.loads(out.read_text()); existing={(x['cluster'],x['gene_id']):x for x in data['markers']}; denom=M['original_nCount_RNA'].to_numpy(dtype=float)
assert np.all(denom>0)
for gene in genes:
 pos=np.flatnonzero(A.var['gene_symbol'].to_numpy()==gene)
 print(gene, 'matches',len(pos))
 if len(pos)!=1: continue
 p=pos[0]; gene_id=str(A.var_names[p]); v=np.asarray(A.layers['counts'][:,p].toarray()).ravel(); norm=np.log1p(v*10000/denom)
 group=M['cluster_id'].to_numpy(); level=sorted(set(group),key=int); means=np.array([norm[group==c].mean() for c in level]); mean=means.mean(); sd=means.std(ddof=0)
 for c in ids:
  mask=(group==c); n=int(mask.sum()); det=int((v[mask]>0).sum()); mn=float(norm[mask].mean()); z=float((mn-mean)/sd) if sd else None
  item=dict(cluster='C'+c,gene=gene,gene_id=gene_id,n_nuclei=n,n_detected=det,fraction_detected=det/n,mean_log1p_normalized=mn,scaled_expression_z=z)
  if gene=='Cdh7':
   orig=existing[(item['cluster'],gene_id)]
   assert abs(orig['fraction_detected']-item['fraction_detected'])<1e-9
   assert abs(orig['mean_log1p_normalized']-mn)<1e-6,(c,orig['mean_log1p_normalized'],mn)
   assert abs(orig['scaled_expression_z']-z)<1e-6,(c,orig['scaled_expression_z'],z)
  elif (item['cluster'],gene_id) not in existing:data['markers'].append(item)
  else:
   orig=existing[(item['cluster'],gene_id)];assert abs(orig['mean_log1p_normalized']-mn)<1e-6 and abs(orig['scaled_expression_z']-z)<1e-6
  print(c, n, det, round(det/n,4), round(mn,5), round(z,5))
data['markers'].sort(key=lambda x:(x['gene'].lower(),x['gene_id'],int(x['cluster'][1:])))
data['provenance']['supplemental_marker_source']=str(h5)
data['provenance']['supplemental_markers']='Elfn1, Chodl, Chrna2, Ctnnd2, Htr3a recovered from original-full UMI counts with 8,911/8,911 IDs mapped to frozen Seurat clusters; log1p normalized to 10,000 original full-library UMIs; z computed over all Seurat inhibitory cluster means, ddof 0. Counts are from the shared source count layer, not Scanpy cluster assignments. Validated against saved Seurat Cdh7, Sst, Reln, Trpc6, Pcdh11x and Grin3a summaries.'
data['provenance']['missing']='Markers absent from the available saved and reconstructed summaries are unavailable, not necessarily zero; no imputation.'
out.write_text(json.dumps(data,separators=(',',':'),allow_nan=False)+'\n')
print('total markers',len(data['markers']))
