from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, imageio.v2 as imageio
from sklearn.decomposition import PCA
from PIL import Image,ImageDraw,ImageFont
P=Path('/Users/mujahid357/Documents/codex_new/SATB2_SCENE_SEURAT_20261111'); a=ad.read_h5ad(P/'seurat_ca1_scene_umap.h5ad'); cells=a.obsm['SCENE'];genes=a.varm['SCENE']; pca=PCA(n_components=3,svd_solver='full').fit(np.vstack((cells,genes))); C=pca.transform(cells);G=pca.transform(genes); variance=float(sum(pca.explained_variance_ratio_)); syms=a.var.symbol.astype(str).to_numpy(); labelgenes=['Satb2','Ano3','Ndst4','Kcnip4','Ntm','Asic2']; positions={s:G[np.flatnonzero(syms==s)[0]] for s in labelgenes if sum(syms==s)==1}; assert len(positions)==6
labs=a.obs.cluster_id.astype(str).to_numpy();groups=sorted(set(labs),key=int); layers=a.obs.layer_program_label.astype(str).to_numpy();L=['Superficial-associated expression program enriched','Deep-associated expression program enriched','Mixed or weak layer-associated expression programs','Layer program not assigned: conflicting identity signatures']; Ccolors=['#4477aa','#8fb8d5','#ee9966','#f4c29e','#378f5c','#8bd3a2','#d95858','#f4a6b4','#9568ac','#cbb1da','#876b5c'];Lcolors=['#468eae','#bd6847','#9da7ad','#9274aa']; W,H=1088,1360; bg='#071017'; normal='/System/Library/Fonts/Supplemental/Arial.ttf';bold='/System/Library/Fonts/Supplemental/Arial Bold.ttf'; fbig=ImageFont.truetype(bold,42);fmid=ImageFont.truetype(normal,24);fsmall=ImageFont.truetype(normal,20);fmark=ImageFont.truetype(bold,23)
sel=pd.read_csv(P/'grouped_30_marker_selection.tsv',sep='\t');m=pd.read_csv(P/'cluster_mean_log1p_cp10k.tsv',sep='\t',index_col=0)[sel.gene_symbol];z=((m-m.mean())/m.std(ddof=0).replace(0,np.nan)).fillna(0).clip(-2,2).to_numpy();h=pd.read_csv(P/'training_history.tsv',sep='\t'); loss=h.iloc[:,-1].to_numpy() if 'loss' not in h.columns else h['loss'].to_numpy(); loss=np.asarray(loss,dtype=float); print('history columns',h.columns.tolist(),'PCA variance',variance,flush=True)
# A common orthographic scale for nuclei and selected genes; only camera rotates.
allxyz=np.vstack([C,np.stack(list(positions.values()))]); radius=float(np.max(np.linalg.norm(allxyz,axis=1))); scale=420/radius; center=np.array([W/2,710]);rng=np.random.default_rng(42); order=rng.permutation(len(C));
def rotate(data,angle):
 x=data[:,0]*np.cos(angle)-data[:,1]*np.sin(angle);y=(data[:,0]*np.sin(angle)+data[:,1]*np.cos(angle))*.42+data[:,2]*.9075;return np.column_stack([center[0]+x*scale,center[1]-y*scale])
def rgb(hex):return tuple(bytes.fromhex(hex[1:]));
def blend(a,b,t):return tuple(int(x*(1-t)+y*t) for x,y in zip(a,b))
def heatcolor(v):
 t=(v+2)/4;return blend((55,126,184),(247,247,247),t*2) if t<.5 else blend((247,247,247),(185,38,67),(t-.5)*2)
def header(draw,subtitle):
 draw.text((58,46),'SCENE  |  CA1 / ProS',font=fbig,fill='#e8f1f5');draw.text((60,108),subtitle,font=fmid,fill='#b4c7d0');draw.text((60,142),'Seurat annotations · 27,703 nuclei · 12 mice · 2,518 genes',font=fsmall,fill='#91a8b4');draw.line((60,185,1020,185),fill='#385160',width=2)
def footer(draw,phase):
 draw.line((60,1185,1020,1185),fill='#385160',width=2);draw.text((60,1200),'FINAL FIT  600 / 600 epochs',font=fsmall,fill='#c2d5df');draw.text((60,1232),'Animated camera only · no inferred cell movement or developmental time',font=fsmall,fill='#91a8b4');draw.text((60,1270),phase,font=fsmall,fill='#72c2c4')
def heat_frame():
 im=Image.new('RGB',(W,H),bg);d=ImageDraw.Draw(im);header(d,'01  |  Measured marker expression across CA1 clusters');
 x0,y0=130,375;w=27;ht=43
 for j in range(30):
  if j%3==0:d.line((x0+j*w-2,y0-23,x0+j*w-2,y0+len(m)*ht),fill='#5a6e78',width=1);d.text((x0+j*w+9,y0-50),sel.marker_group.iloc[j],font=fsmall,fill='#bdd0d8')
  for i in range(len(m)):
   d.rounded_rectangle((x0+j*w,y0+i*ht,x0+j*w+w-3,y0+i*ht+ht-4),radius=3,fill=heatcolor(z[i,j]))
 for i,n in enumerate(m.index):d.text((75,y0+i*ht+8),n,font=fsmall,fill='#d2dfe6')
 for j,gene in enumerate(sel.gene_symbol):
  label=Image.new('RGBA',(175,22),(0,0,0,0));ld=ImageDraw.Draw(label);ld.text((0,1),gene,font=ImageFont.truetype(normal,15),fill='#bdd0d8');label=label.rotate(55,expand=True,resample=Image.Resampling.BICUBIC);im.paste(label,(x0+j*w-3,y0+len(m)*ht+8),label)
 d.text((130,1080),'30 cluster-grouped candidates · red: higher relative mean expression',font=fsmall,fill='#b4c7d0');d.text((130,1110),'Blue: lower relative mean expression · log1p(CP10k), gene-wise z-score',font=fsmall,fill='#b4c7d0');footer(d,'01  MARKERS    →    02  CLUSTERS    →    03  LAYER PROGRAMS');return np.asarray(im)
heat=heat_frame();fps=20;dur=28;writer=imageio.get_writer(str(P/'SCENE_CA1_dynamic_final.mp4'),fps=fps,codec='libx264',quality=6,macro_block_size=16)
for frame in range(fps*dur):
 t=frame/fps
 if t<4:writer.append_data(heat);continue
 im=Image.new('RGB',(W,H),bg);d=ImageDraw.Draw(im);phase='02  |  Seurat cluster labels' if t<13 else ('03  |  CA1 layer-associated programs' if t<21 else '04  |  Selected gene coordinates');header(d,phase);d.text((60,208),f'Rotating 3D PCA of final 16D SCENE fit · {variance*100:.1f}% variance retained',font=fsmall,fill='#92adb8');angle=.4+(t-4)*.19;xy=rotate(C,angle)
 if t<13:
  for k,g in enumerate(groups):
   idx=order[labs[order]==g];pts=[tuple(map(int,q)) for q in xy[idx]];d.point(pts,fill=Ccolors[k])
  for k,g in enumerate(groups):
   x=65+(k%6)*155;y=250+(k//6)*33;d.ellipse((x,y,x+11,y+11),fill=Ccolors[k]);d.text((x+17,y-4),'C'+g,font=fsmall,fill='#d8e5eb')
 else:
  for k,l in enumerate(L):
   idx=order[layers[order]==l];d.point([tuple(map(int,q)) for q in xy[idx]],fill=Lcolors[k])
  for k,n in enumerate(['Superficial-associated','Deep-associated','Mixed / weak','Unassigned / conflict']):
   x=65+(k%2)*440;y=250+(k//2)*33;d.ellipse((x,y,x+11,y+11),fill=Lcolors[k]);d.text((x+18,y-4),n,font=fsmall,fill='#d8e5eb')
 if t>=21:
  Gxy=rotate(np.stack(list(positions.values())),angle); offsets=[(-165,-90),(70,-75),(-170,65),(75,45),(-160,120),(80,110)]
  for (name,_),(x,y),(dx,dy) in zip(positions.items(),Gxy,offsets):
   d.ellipse((x-7,y-7,x+7,y+7),fill='#ffffff',outline='#ffd8a6',width=2);d.line((x,y,x+dx,y+dy),fill='#f4c39c',width=2);d.text((x+dx,y+dy-15),name,font=fmark,fill='#ffd4b6')
  d.text((60,1110),'Gene labels mark projected gene coordinates, not gene expression.',font=fsmall,fill='#bed0d9')
 footer(d,'03D VIEW: PCA PROJECTION  ·  marker dots are descriptive');writer.append_data(np.asarray(im))
 if frame%(fps*5)==0:print('seconds',t,flush=True)
writer.close();(P/'dynamic_video_receipt.json').write_text(__import__('json').dumps({'fps':fps,'seconds':dur,'n_nuclei':len(C),'n_genes':len(G),'pca_variance':variance,'markers':labelgenes,'camera':'rotation of final fitted 3D PCA projection; no training time interpolation','colors':'saved cluster and layer annotations','seed':42},indent=2));print('done',flush=True)
