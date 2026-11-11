from pathlib import Path
import json, numpy as np,pandas as pd,anndata as ad,imageio.v2 as imageio
from sklearn.decomposition import PCA
from PIL import Image,ImageDraw,ImageFont
P=Path('/Users/mujahid357/Documents/codex_new/SATB2_SCENE_SEURAT_20261111');a=ad.read_h5ad(P/'seurat_ca1_scene_umap.h5ad');X=a.obsm['SCENE'];Y=a.varm['SCENE'];p=PCA(n_components=3,svd_solver='full').fit(np.vstack([X,Y]));C=p.transform(X);G=p.transform(Y);ratio=float(p.explained_variance_ratio_.sum());names=a.var.symbol.astype(str).to_numpy();selected=['Satb2','Ano3','Ndst4','Kcnip4','Ntm','Asic2']; pos={};edges=[]
for gene in selected:
 j=np.flatnonzero(names==gene);assert len(j)==1;j=int(j[0]);pos[gene]=G[j];counts=np.asarray(a.X[:,j].toarray()).ravel(); candidates=np.flatnonzero(counts>0);assert len(candidates)>2; distances=np.sum((X[candidates]-Y[j])**2,axis=1); chosen=candidates[np.argsort(distances)[:3]]
 for i in chosen:edges.append({'gene':gene,'cell_index':int(i),'nucleus_id':str(a.obs_names[i]),'umi_count':int(counts[i]),'distance_16d':float(np.sqrt(np.sum((X[i]-Y[j])**2)))})
pd.DataFrame(edges).to_csv(P/'scene_v3_observed_cell_gene_edges.tsv',sep='\t',index=False)
markers=pd.read_csv(P/'grouped_30_marker_selection.tsv',sep='\t');M=pd.read_csv(P/'cluster_mean_log1p_cp10k.tsv',sep='\t',index_col=0)[markers.gene_symbol];Z=((M-M.mean())/M.std(ddof=0).replace(0,np.nan)).fillna(0).clip(-2,2).to_numpy()
W,H=1088,1360;bg='#050a10'; cluster=a.obs.cluster_id.astype(str).to_numpy();cat=sorted(set(cluster),key=int);layer=a.obs.layer_program_label.astype(str).to_numpy();L=['Superficial-associated expression program enriched','Deep-associated expression program enriched','Mixed or weak layer-associated expression programs','Layer program not assigned: conflicting identity signatures'];cc=['#4477aa','#8fb8d5','#ee9966','#f4c29e','#378f5c','#8bd3a2','#d95858','#f4a6b4','#9568ac','#cbb1da','#876b5c'];lc=['#468eae','#bd6847','#9da7ad','#9274aa']; f='/System/Library/Fonts/Supplemental/Arial.ttf';fb='/System/Library/Fonts/Supplemental/Arial Bold.ttf';F=lambda n:ImageFont.truetype(f,n);B=lambda n:ImageFont.truetype(fb,n);white='#eef4f7';muted='#9fb6c3';teal='#78d1d0';rng=np.random.default_rng(42);order=rng.permutation(len(C));r=float(np.linalg.norm(np.vstack([C,np.stack(list(pos.values()))]),axis=1).max());scale=420/r;center=(W/2,695)
def project(X,theta):
 xx=X[:,0]*np.cos(theta)-X[:,1]*np.sin(theta); yy=(X[:,0]*np.sin(theta)+X[:,1]*np.cos(theta))*.42+X[:,2]*.9075;return np.column_stack((center[0]+xx*scale,center[1]-yy*scale))
def mix(a,b,t):return tuple(round(x*(1-t)+y*t) for x,y in zip(a,b))
def color(v):
 t=np.clip((v+2)/4,0,1);return mix((55,126,184),(247,247,247),t*2) if t<.5 else mix((247,247,247),(185,38,67),(t-.5)*2)
def header(d,sub):
 d.text((55,44),'SCENE  |  Joint embedding',font=B(41),fill=white);d.text((58,105),'P60 CA1 / ProS · Seurat clusters · WT and Satb2 cKO',font=F(23),fill=muted);d.text((58,138),'27,703 nuclei · 2,518 genes · 12 mice · saved final model',font=F(20),fill=muted);d.line((58,182,1030,182),fill='#3a5561',width=2);d.text((58,204),sub,font=B(25),fill='#d5e7ed')
def footer(d,phase):
 d.line((58,1180,1030,1180),fill='#3a5561',width=2);d.text((58,1199),'FINAL FIT: 600 / 600 epochs',font=B(20),fill='#d0e2e8');d.text((58,1230),'Camera rotates around saved fit; no biological movement is inferred.',font=F(19),fill=muted);d.text((58,1274),phase,font=B(20),fill=teal)
def heat(k):
 im=Image.new('RGB',(W,H),bg);d=ImageDraw.Draw(im);header(d,f'01  |  pheatmap: CA1 cluster markers ({k+1}/2)');x0,y0=174,355;cw,ch=52,49
 for j in range(15):
  jj=k*15+j
  if jj%3==0:d.line((x0+j*cw-2,y0-14,x0+j*cw-2,y0+11*ch),fill='#526875',width=2);d.text((x0+j*cw+12,y0-43),str(markers.marker_group.iloc[jj]),font=B(18),fill='#c6d9df')
  for i in range(11):d.rounded_rectangle((x0+j*cw,y0+i*ch,x0+j*cw+cw-3,y0+i*ch+ch-4),3,fill=color(Z[i,jj]))
  text=Image.new('RGBA',(195,25));td=ImageDraw.Draw(text);td.text((0,1),markers.gene_symbol.iloc[jj],font=F(20),fill='#d5e5ed');text=text.rotate(55,expand=True,resample=Image.Resampling.BICUBIC);im.paste(text,(x0+j*cw,y0+11*ch+10),text)
 for i,n in enumerate(M.index):d.text((110,y0+i*ch+12),n,font=F(23),fill=white)
 d.text((78,1103),'Red / blue = higher / lower gene-wise z-score of mean expression.',font=F(20),fill=muted);footer(d,'01  PHEATMAP     02  COUNT EDGES     03  3D EMBEDDING');return np.asarray(im)
heats=[heat(0),heat(1)];fps=20;seconds=40;writer=imageio.get_writer(str(P/'SCENE_CA1_interactions_v3.mp4'),fps=fps,codec='libx264',quality=6,macro_block_size=16)
for frame in range(seconds*fps):
 t=frame/fps
 if t<10:writer.append_data(heats[int(t//5)]);continue
 im=Image.new('RGB',(W,H),bg);d=ImageDraw.Draw(im)
 if t<19:sub='02  |  Nucleus–gene count edges'; mode='cluster'
 elif t<29:sub='03  |  CA1 clusters and marker positions';mode='cluster'
 else:sub='04  |  Layer programs and marker positions';mode='layer'
 header(d,sub);d.text((58,248),f'3D PCA of the final 16D SCENE embedding · {ratio*100:.1f}% variance retained',font=F(19),fill=muted)
 if mode=='cluster':
  for k,c in enumerate(cat):
   x=62+(k%6)*158;y=285+(k//6)*34;d.ellipse((x,y,x+12,y+12),fill=cc[k]);d.text((x+19,y-4),'C'+c,font=F(19),fill=white)
 else:
  for k,n in enumerate(['Superficial-associated','Deep-associated','Mixed / weak','Unassigned / conflict']):
   x=62+(k%2)*465;y=285+(k//2)*34;d.ellipse((x,y,x+12,y+12),fill=lc[k]);d.text((x+20,y-4),n,font=F(19),fill=white)
 theta=.2+(t-10)*.16;xy=project(C,theta)
 if mode=='cluster':
  for k,c in enumerate(cat):d.point([tuple(map(int,q)) for q in xy[order[cluster[order]==c]]],fill=cc[k])
 else:
  for k,l in enumerate(L):d.point([tuple(map(int,q)) for q in xy[order[layer[order]==l]]],fill=lc[k])
 gg=project(np.stack([pos[g] for g in selected]),theta)
 if t<19:
  for e in edges:
   gx,gy=gg[selected.index(e['gene'])];cx,cy=xy[e['cell_index']];d.line((gx,gy,cx,cy),fill='#6d5661',width=1);d.ellipse((cx-3,cy-3,cx+3,cy+3),fill='#e8e0bf')
  d.text((60,1095),'18 example edges: each nucleus has a nonzero UMI count for the gene.',font=F(19),fill='#d4dfe3');d.text((60,1122),'These are expression edges, not ligand–receptor interactions.',font=F(19),fill='#d4dfe3')
 else:d.text((60,1110),'Diamond = gene coordinate; label leaders point to its projected position.',font=F(19),fill='#cbdbe2')
 # Fixed, non-overlapping label lanes; connecting leaders end at actual projected points.
 for i,(gene,(x,y)) in enumerate(zip(selected,gg)):
  left=i<3; tx=60 if left else 860;ty=500+(i%3)*108;endx=tx+145 if left else tx-14;d.line((x,y,endx,ty+13),fill='#d5a0a4',width=2);d.polygon([(x,y-9),(x+9,y),(x,y+9),(x-9,y)],fill='#ffe6dd',outline='#ba9094');d.text((tx,ty),gene,font=B(24),fill='#f1c4bd');d.text((tx,ty+28),'gene position',font=F(16),fill=muted)
 footer(d,'02  COUNT EDGES     →     03  CLUSTERS     →     04  LAYERS');writer.append_data(np.asarray(im))
 if frame%(fps*5)==0:print('seconds',t,flush=True)
writer.close();(P/'SCENE_v3_receipt.json').write_text(json.dumps({'fps':fps,'duration_s':seconds,'nuclei':len(C),'genes':len(Y),'pca_variance':ratio,'gene_labels':selected,'edges':len(edges),'edge_definition':'3 observed nonzero-UMI nuclei per gene selected by nearest 16D distance; expression evidence only','heatmap':'30 candidates, 11 clusters, gene-wise z-scores of mean log1p(CP10k)','recorded_training_trajectory':False},indent=2));print('done',flush=True)
