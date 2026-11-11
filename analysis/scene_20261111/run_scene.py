import os, time, json, gzip
os.environ['MPLBACKEND']='Agg'
os.environ['OMP_NUM_THREADS']='4'
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
import matplotlib.pyplot as plt
from pathlib import Path
from scene import SCENE, joint_umap
import torch
P=Path('/Users/mujahid357/Documents/codex_new/SATB2_SCENE_SEURAT_20261111'); P.mkdir(exist_ok=True)
base=Path('/Users/mujahid357/Documents/codex_new')
raw=base/'SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/shared/original_full/branch_counts.h5ad'
meta=base/'SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0/results/plots/seurat/original_full/CA1_ProS/original_frozen/rendered_nucleus_metadata.tsv.gz'
m=pd.read_csv(meta,sep='\t',low_memory=False); m=m[m.keep_for_analysis.astype(str).str.lower().isin(['true','1'])].copy(); assert m.nucleus_id.is_unique
A=ad.read_h5ad(raw); ix=A.obs_names.get_indexer(m.nucleus_id); assert (ix>=0).all(); A=A[ix,:].copy(); assert A.n_obs==len(m); m.index=A.obs_names
for c in ['sample','condition','cluster_id','layer_program_label']: A.obs[c]=m[c].astype(str)
assert np.all(A.X.data>=0) and np.all(A.X.data==np.floor(A.X.data))
print('input',A.shape,A.obs['sample'].nunique(),A.obs.condition.value_counts().to_dict(),flush=True)
# Descriptive mean log1p library-normalized counts per Seurat cluster; no genotype DEG test.
lib=np.asarray(A.X.sum(axis=1)).ravel(); N=sp.diags(10000/np.maximum(lib,1)) @ A.X; N.data=np.log1p(N.data)
clusters=sorted(A.obs.cluster_id.unique(),key=lambda x:int(x)); means=np.vstack([np.asarray(N[np.asarray(A.obs.cluster_id==c)].mean(axis=0)).ravel() for c in clusters]); detect=np.vstack([np.asarray((A.X[np.asarray(A.obs.cluster_id==c)]>0).mean(axis=0)).ravel() for c in clusters]); sym=A.var.gene_symbol.fillna('').astype(str).to_numpy(); unique=pd.Series(sym).value_counts(); ok=np.array([bool(x) and unique[x]==1 for x in sym]); contrast=means-means.mean(axis=0); markeridx=[]
for j,c in enumerate(clusters):
 rank=np.argsort(-(contrast[j]*np.sqrt(detect[j]))); selected=[i for i in rank if ok[i] and detect[j,i]>=.05][:5]; markeridx.extend(selected)
known=['Satb2','Ano3','Wfs1','Calb1','Ndst4','Kcnip4','Col11a1','Pde10a','Asic2','Ntm','Cntnap5a','Dcc','Gpc5','Rora','Pvalb','Gad1','Slc17a7']
markeridx=list(dict.fromkeys(markeridx+[i for gene in known for i in np.where(sym==gene)[0] if ok[i]]))
# Select 2500 most detected features plus requested and descriptive markers for tractable model fit.
frequency=np.asarray((A.X>0).sum(axis=0)).ravel(); idx=np.unique(np.r_[np.argsort(-frequency)[:2500],markeridx]); idx=idx[frequency[idx]>0]; sub=A[:,idx].copy(); sub.var['symbol']=sym[idx]; sub.obs['sample']=sub.obs['sample'].astype('category'); del A
pd.DataFrame({'gene_id':sub.var_names,'gene_symbol':sub.var.symbol,'detected_nuclei':frequency[idx]}).to_csv(P/'model_features.tsv',sep='\t',index=False)
pd.DataFrame({'gene_id':[str(N) for N in [A for A in []]]}) if False else None
pd.DataFrame(means[:,markeridx],index=['C'+c for c in clusters],columns=sym[markeridx]).to_csv(P/'cluster_mean_log1p_cp10k.tsv',sep='\t'); pd.DataFrame(detect[:,markeridx],index=['C'+c for c in clusters],columns=sym[markeridx]).to_csv(P/'cluster_detected_fraction.tsv',sep='\t')
sub.write_h5ad(P/'seurat_ca1_scene_input.h5ad',compression='gzip'); print('model input',sub.shape,'markers',len(markeridx),flush=True)
model=SCENE(latent_dim=16,likelihood='ZIP'); t=time.time(); model.fit(sub,layer=None,batch_keys=['sample'],epochs=600,seed=42,validate=True,val_interval=20,print_every=10,cell_batch_size=512,loader_workers=0,device='cpu'); print('fit seconds',time.time()-t,flush=True)
pd.DataFrame(model.training_history['train']).to_csv(P/'training_history.tsv',sep='\t',index=False)
model.save(P/'scene_model'); sub.write_h5ad(P/'seurat_ca1_scene_fit.h5ad',compression='gzip'); print('fit saved',flush=True)
joint_umap(sub,seed=42); sub.write_h5ad(P/'seurat_ca1_scene_umap.h5ad',compression='gzip'); print('umap saved',flush=True)
