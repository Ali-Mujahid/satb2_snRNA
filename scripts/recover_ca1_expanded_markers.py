import anndata as ad,pandas as pd,numpy as np,sys
from pathlib import Path
base=Path('/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0/results')
a=ad.read_h5ad('/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/runs/20261013_v1_0_5/results/shared/original_full/branch_counts.h5ad',backed='r')
m=pd.read_csv(base/'plots/seurat/original_full/CA1_ProS/original_frozen/rendered_nucleus_metadata.tsv.gz',sep='\t',usecols=['nucleus_id','source_cluster_id','condition']);m=m[m.condition=='WT'].set_index('nucleus_id');assert m.index.is_unique
idx=a.obs.nucleus_id.isin(m.index).to_numpy();obs=a.obs.loc[idx,['nucleus_id','nCount_RNA']].copy().reset_index(drop=True);obs=obs.join(m,on='nucleus_id',validate='one_to_one');assert len(obs)==len(m)
print('WT matched',len(obs),flush=True)
sup='Asic2 Ntm Cntnap5a Galntl6 Dcc Kcnq3 AC129186.1 Epha3 Fstl4 Erc2 Grm1 Nav2 Kcnc2 Tenm4 Nlgn1 March1 Robo2 Gpc5 Lingo2 Mdga1 Ak5 Fgf14 Pcdh7 Ryr2 Ptprk Ralyl Kcnab1 Thsd7b Adcy2 Sorbs2 Arl15 Xkr4 Ppp3ca Ano3 Fam19a1 Rora Alcam Phactr1 Nrg2 Plcb1'.split();deep='Ndst4 Kcnip4 Col11a1 Pde10a Fras1 Epha10 Ntrk2 Pex5l Astn2 Pde4d Sorcs3 Zfp804a Brinp1 Cacna2d3 2010300C02Rik Man1a R3hdm1 Spock3 Prkca Ptprd Opcml Zfpm2 Nr3c2 Epha6 Slc24a2 Malat1 Nr4a1 Srgap1 Stac Arhgap26 Adgrb3 Enox1 Plekhg1 Ap2b1 Hdac9 Rab3c Plekha7 Zbtb16 Pde1a 4930578G10Rik'.split(); genes=sup+deep
sym=a.var.gene_symbol.astype(str);indices=[i for i,g in enumerate(sym) if g in genes and (sym==g).sum()==1];print('uniquely mapped',len(indices),flush=True)
X=a.X[idx,:][:,indices].tocsr();lib=obs.nCount_RNA.to_numpy();assert (lib>0).all();norm=X.multiply((1e4/lib)[:,None]).tocsr();norm.data=np.log1p(norm.data)
rows=[]
for cl in range(1,12):
 mask=(obs.source_cluster_id.to_numpy()==cl);z=norm[mask];r=X[mask]
 for j,col in enumerate(indices):rows.append({'cluster_id':str(cl),'gene_symbol':sym.iloc[col],'mean_log1p_normalized':float(z[:,j].mean()),'fraction_detected':float((r[:,j]>0).sum()/mask.sum())})
pd.DataFrame(rows).to_csv(Path(__file__).parent/'recovered_wt.tsv',sep='\t',index=False);print('complete',flush=True)
