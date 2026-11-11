import json,anndata as ad,pandas as pd,numpy as np
from pathlib import Path
P=Path('/Users/mujahid357/Documents/codex_new/SATB2_SCENE_SEURAT_20261111'); D=P/'atlas';a=ad.read_h5ad(P/'seurat_ca1_scene_umap.h5ad'); xy=a.obsm['SCENE_umap']; assert len(xy)==27703 and np.isfinite(xy).all()
layers=['Superficial-associated expression program enriched','Deep-associated expression program enriched','Mixed or weak layer-associated expression programs','Layer program not assigned: conflicting identity signatures']; vals=a.obs.layer_program_label.astype(str);assert set(vals).issubset(layers)
rows=[[round(float(x),2),round(float(y),2),int(c),layers.index(l),str(g)] for (x,y),c,l,g in zip(xy,a.obs.cluster_id.astype(str),vals,a.obs.condition.astype(str))];
sel=pd.read_csv(P/'grouped_30_marker_selection.tsv',sep='\t'); data={'points':rows,'layers':['Superficial-associated','Deep-associated','Mixed / weak','Unassigned / conflict'],'markers':{},'preset':sel.gene_symbol.tolist(),'groups':dict(zip(sel.gene_symbol,sel.marker_group))}
for kind in ['cluster','layer']:
 m=pd.read_csv(P/f'{kind}_mean_log1p_cp10k.tsv',sep='\t',index_col=0);d=pd.read_csv(P/f'{kind}_detected_fraction.tsv',sep='\t',index_col=0);assert m.shape==d.shape and m.index.equals(d.index) and m.columns.equals(d.columns)
 v=m if kind=='layer' else ((m-m.mean())/m.std(ddof=0).replace(0,np.nan)).fillna(0).clip(-2,2);data['markers'][kind]={'rows':m.index.tolist(),'genes':m.columns.tolist(),'mean':np.round(m.values,5).tolist(),'fraction':np.round(d.values,5).tolist(),'color':np.round(v.values,5).tolist()}
syms=a.var.symbol.astype(str).to_numpy();gp=a.varm['SCENE_umap']; data['gene_positions']={g:[round(float(gp[ix[0],0]),2),round(float(gp[ix[0],1]),2)] for g in data['markers']['cluster']['genes'] if len(ix:=np.flatnonzero(syms==g))==1}
with open(D/'data.js','w') as f:f.write('window.SCENE_DATA='+json.dumps(data,separators=(',',':'),allow_nan=False)+';\n')
print('points',len(rows),'markers',len(data['markers']['cluster']['genes']),'bytes',(D/'data.js').stat().st_size)
