from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad,imageio.v2 as imageio
from sklearn.decomposition import PCA
from PIL import Image,ImageDraw,ImageFont
P=Path('/Users/mujahid357/Documents/codex_new/SATB2_SCENE_SEURAT_20261111');a=ad.read_h5ad(P/'seurat_ca1_scene_umap.h5ad');X=a.obsm['SCENE'];Y=a.varm['SCENE'];pc=PCA(n_components=3,svd_solver='full').fit(np.vstack([X,Y]));C=pc.transform(X);G=pc.transform(Y);variance=float(pc.explained_variance_ratio_.sum());sym=a.var.symbol.astype(str).to_numpy();L=['Superficial-associated expression program enriched','Deep-associated expression program enriched'];short=['Superficial-associated','Deep-associated'];colors=['#53b3e5','#f0a168'];ms=[['Satb2','Ano3','Ntm','Asic2'],['Ndst4','Kcnip4','Pde10a','Col11a1']];genes=['Satb2','Ano3','Ntm','Asic2','Cntnap5a','Gpc5','Rora','Calb1','Ndst4','Kcnip4','Pde10a','Col11a1','Ntrk2','Zfp804a','Prkca','Nr3c2'];gidx={g:int(np.flatnonzero(sym==g)[0]) for g in genes};assert all(np.sum(sym==g)==1 for g in genes)
layer=a.obs.layer_program_label.astype(str).to_numpy(); sample=a.obs['sample'].astype(str).to_numpy(); mice=sorted(set(sample));assert len(mice)==12;rng=np.random.default_rng(42);chosen=[]
for l in L:
 ix=np.concatenate([rng.choice(np.flatnonzero((layer==l)&(sample==mouse)),4,replace=False) for mouse in mice]);chosen.append(ix)
counts=np.asarray(a.X[np.concatenate(chosen)][:,[gidx[g] for g in genes]].toarray()); assert counts.shape==(96,16) and np.all(counts>=0) and np.all(counts==counts.astype(int));pd.DataFrame(counts,columns=genes).to_csv(P/'scene_v4_heatmap_counts.tsv',sep='\t',index=False); ci=np.concatenate(chosen); obs=a.obs.iloc[ci]; pd.DataFrame({'cell_index':ci,'nucleus_id':a.obs_names[ci],'mouse':obs['sample'].astype(str).to_numpy(),'genotype':obs['condition'].astype(str).to_numpy(),'layer_program':obs.layer_program_label.astype(str).to_numpy(),'seurat_cluster':obs.cluster_id.astype(str).to_numpy()}).to_csv(P/'scene_v4_heatmap_metadata.tsv',sep='\t',index=False)
edges=[]
for k in range(2):
 for gene in ms[k]:
  j=gidx[gene]; raw=np.asarray(a.X[chosen[k],j].toarray()).ravel(); eligible=chosen[k][raw>0];assert len(eligible)>=3;distance=np.sum((X[eligible]-Y[j])**2,axis=1); picked=eligible[np.argsort(distance)[:3]]
  for i in picked:edges.append({'layer':short[k],'gene':gene,'cell_index':int(i),'nucleus_id':str(a.obs_names[i]),'umi_count':int(a.X[i,j]),'distance_16d':float(np.linalg.norm(X[i]-Y[j]))})
pd.DataFrame(edges).to_csv(P/'scene_v4_observed_count_edges.tsv',sep='\t',index=False);assert len(edges)==24
W,H=1280,1360;bg='#04090e';white='#f0f7fa';muted='#a2b8c4';accent='#73d4d3';gold='#ffd0b8'; font='/System/Library/Fonts/Supplemental/Arial.ttf';bold='/System/Library/Fonts/Supplemental/Arial Bold.ttf';F=lambda n:ImageFont.truetype(font,n);B=lambda n:ImageFont.truetype(bold,n);allpts=np.vstack([C,np.stack([G[gidx[g]] for g in sum(ms,[])])]);radius=float(np.linalg.norm(allpts,axis=1).max());scale=440/radius;center=(W/2,725);order=np.random.default_rng(10).permutation(len(C))
def project(arr,theta):
 xx=arr[:,0]*np.cos(theta)-arr[:,1]*np.sin(theta);yy=(arr[:,0]*np.sin(theta)+arr[:,1]*np.cos(theta))*.42+arr[:,2]*.9075;return np.column_stack([center[0]+xx*scale,center[1]-yy*scale])
def blend(a,b,t):return tuple(int(x*(1-t)+y*t) for x,y in zip(a,b))
def heat(c):return blend((12,31,42),(128,224,214),min(1,np.log1p(c)/np.log(11)))
def top(d,subtitle):
 d.text((55,40),'SCENE  |  CA1 nuclei and genes',font=B(40),fill=white);d.text((57,101),'P60 · Seurat · WT and Satb2 cKO · 12 mice',font=F(24),fill=muted);d.text((57,137),'27,703 CA1 / ProS nuclei · 2,518 selected genes',font=F(20),fill=muted);d.line((57,180,1222,180),fill='#3b5361',width=2);d.text((57,203),subtitle,font=B(25),fill=white)
def foot(d,label):
 d.line((57,1194,1222,1194),fill='#3b5361',width=2);d.text((57,1210),'Fixed final fit · rotating camera, not a biological trajectory',font=F(19),fill=muted);d.text((57,1242),'3D PCA of 16D SCENE coordinates · 29.0% variance retained',font=F(19),fill=muted);d.text((57,1284),label,font=B(19),fill=accent)
def make_heat():
 im=Image.new('RGB',(W,H),bg);d=ImageDraw.Draw(im);top(d,'01  |  Count matrix excerpt · 96 nuclei × 16 genes');x0,y0=256,400;cw,ch=48,7
 for j,g in enumerate(genes):
  d.line((x0+j*cw+23,y0-14,x0+j*cw+23,y0-3),fill='#ceaaa5',width=1)
  label=Image.new('RGBA',(160,23));q=ImageDraw.Draw(label);q.text((0,2),g,font=F(19),fill='#e4c4be');label=label.rotate(55,expand=True,resample=Image.Resampling.BICUBIC);im.paste(label,(x0+j*cw-3,y0-156),label)
  for i in range(96):
   x=x0+j*cw;y=y0+i*ch;d.rectangle((x,y,x+cw-3,y+ch-1),fill=heat(counts[i,j]))
 for i in range(96):
  y=y0+i*ch+2;d.ellipse((216,y-2,220,y+2),fill=colors[i//48])
 d.line((x0-10,y0+48*ch-3,x0+16*cw,y0+48*ch-3),fill='#77949d',width=2)
 d.text((57,1094),'Each row: one randomly sampled nucleus · 4 per mouse and group',font=F(20),fill=muted);d.text((57,1122),'Cell color: blue superficial-associated, orange deep-associated',font=F(20),fill=muted);d.text((57,1150),'Teal intensity: log1p(raw UMI count), capped at 10 for display',font=F(20),fill=muted);foot(d,'01  COUNTS  →  02  SUPERFICIAL  →  03  DEEP');return np.asarray(im)
heatframe=make_heat();fps=20;seconds=38;writer=imageio.get_writer(str(P/'SCENE_CA1_layers_clear_v4.mp4'),fps=fps,codec='libx264',quality=6,macro_block_size=16)
for frame in range(seconds*fps):
 t=frame/fps
 if t<8:writer.append_data(heatframe);continue
 k=0 if t<23 else 1;start=8 if k==0 else 23;im=Image.new('RGB',(W,H),bg);d=ImageDraw.Draw(im);top(d,('02  |  Superficial-associated nuclei' if k==0 else '03  |  Deep-associated nuclei'));d.text((58,252),f'{np.count_nonzero(layer==L[k]):,} nuclei · same joint embedding and scale',font=F(21),fill=muted);d.ellipse((61,290,74,303),fill=colors[k]);d.text((85,284),short[k]+' nuclei',font=F(20),fill=white);d.ellipse((437,290,450,303),outline=gold,width=3);d.text((461,284),'Gene position',font=F(20),fill=white);d.line((700,298,735,298),fill=accent,width=3);d.text((747,284),'UMI edge',font=F(20),fill=white)
 theta=.25+(t-start)*.20;xy=project(C,theta);idx=order[layer[order]==L[k]];d.point([tuple(map(int,q)) for q in xy[idx]],fill=colors[k]);gg=project(np.stack([G[gidx[g]] for g in ms[k]]),theta)
 # These lines link only real, observed nonzero raw counts; not signaling or anatomical connections.
 for e in [e for e in edges if e['layer']==short[k]]:
  gx,gy=gg[ms[k].index(e['gene'])];cx,cy=xy[e['cell_index']];d.line((gx,gy,cx,cy),fill=accent,width=3);d.ellipse((cx-5,cy-5,cx+5,cy+5),fill='#e2f7e7',outline='#46bbb9',width=2)
 for j,(gene,(x,y)) in enumerate(zip(ms[k],gg)):
  side=j<2;tx=58 if side else 1100;ty=500+(j%2)*160;ex=tx+150 if side else tx-15;d.line((x,y,ex,ty+15),fill='#d9a8a2',width=2);d.ellipse((x-10,y-10,x+10,y+10),fill='#fff0dc',outline=gold,width=3);d.text((tx,ty),gene,font=B(27),fill=gold);d.text((tx,ty+34),'gene position',font=F(18),fill=muted)
 d.text((58,1100),'12 cyan edges: detected nucleus–gene counts in sampled cells.',font=F(20),fill=white);d.text((58,1129),'Not ligand–receptor interactions or physical synapses.',font=F(20),fill=muted);foot(d,'01  COUNTS     02  SUPERFICIAL     03  DEEP');writer.append_data(np.asarray(im))
 if frame%(fps*5)==0:print('seconds',t,flush=True)
writer.close();(P/'SCENE_v4_receipt.json').write_text(json.dumps({'duration_s':seconds,'fps':fps,'nuclei_total':len(C),'pca_variance_retained':variance,'sampling':'random seed 42; 4 nuclei per mouse per layer, all 12 mice','heatmap':'96 nuclei × 16 genes; raw integer UMI count, log1p color display capped at 10','edges':'3 sampled nuclei per gene with nonzero raw UMI, selected by minimum native 16D distance; 24 total','motion':'rotating camera on fixed 3D PCA of final 16D SCENE fit; same camera and scale for both layers','markers':ms,'training_frames_recorded':False},indent=2));print('done',flush=True)
