import os,sys,json,hashlib
from pathlib import Path
os.environ['OMP_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='4'
import numpy as np,pandas as pd,anndata as ad,imageio.v2 as imageio
from PIL import Image,ImageDraw,ImageFont
fit=Path(sys.argv[1]);out=Path(sys.argv[2]);assert not out.exists();out.mkdir(parents=True);A=ad.read_h5ad(fit/'fitted.h5ad'); ck=sorted((fit/'checkpoints').glob('epoch_*.npz'));assert len(ck)==242 and A.shape==(27703,2518);history=pd.read_csv(fit/'training.tsv',sep='\t');assert len(history)==1200
obs=A.obs;cluster=obs.cluster_id.astype(str).to_numpy();layer=obs.layer_program_label.astype(str).to_numpy();sample=obs['sample'].astype(str).to_numpy();symbols=A.var.symbol.astype(str).to_numpy();L=['Superficial-associated expression program enriched','Deep-associated expression program enriched','Mixed or weak layer-associated expression programs','Layer program not assigned: conflicting identity signatures'];colors=['#53b3e5','#f0a168','#9ca8ae','#bba0d3'];cat=sorted(set(cluster),key=int);cc=['#4477aa','#8fb8d5','#ee9966','#f4c29e','#378f5c','#8bd3a2','#d95858','#f4a6b4','#9568ac','#cbb1da','#876b5c'];markers=['Satb2','Ano3','Ntm','Asic2','Cntnap5a','Gpc5','Rora','Calb1','Ndst4','Kcnip4','Pde10a','Col11a1','Ntrk2','Zfp804a','Prkca','Nr3c2'];assert all(np.sum(symbols==g)==1 for g in markers);gi={g:int(np.flatnonzero(symbols==g)[0]) for g in markers};rng=np.random.default_rng(42); samples=sorted(set(sample));assert len(samples)==12;chosen=[np.concatenate([rng.choice(np.flatnonzero((layer==L[k])&(sample==m)),4,replace=False) for m in samples]) for k in range(2)];idx=np.concatenate(chosen);matrix=A.X[idx][:,[gi[x] for x in markers]].toarray();assert matrix.shape==(96,16) and np.all(matrix>=0);pd.DataFrame(matrix,columns=markers).to_csv(out/'heatmap_counts.tsv',sep='\t',index=False)
meta=obs.iloc[idx];pd.DataFrame({'cell_index':idx,'sample':meta['sample'].astype(str).to_numpy(),'condition':meta.condition.astype(str).to_numpy(),'layer':meta.layer_program_label.astype(str).to_numpy(),'cluster':meta.cluster_id.astype(str).to_numpy()}).to_csv(out/'heatmap_metadata.tsv',sep='\t',index=False)
last=np.load(ck[-1]);finalc=last['cells'];finalg=last['genes'];assert np.array_equal(finalc,A.obsm['SCENE']) and np.array_equal(finalg,A.varm['SCENE']);layer_genes=[['Satb2','Ano3','Ntm','Asic2'],['Ndst4','Kcnip4','Pde10a','Col11a1']];edges=[]
for k in range(2):
 for g in layer_genes[k]:
  vals=np.asarray(A.X[chosen[k],gi[g]].toarray()).ravel();ix=chosen[k][vals>0];assert len(ix)>=3;rank=np.argsort(np.linalg.norm(finalc[ix]-finalg[gi[g]],axis=1))[:3]
  for j in ix[rank]:edges.append({'layer':['Superficial-associated','Deep-associated'][k],'gene':g,'cell_index':int(j),'umi_count':int(A.X[j,gi[g]])})
pd.DataFrame(edges).to_csv(out/'observed_expression_edges.tsv',sep='\t',index=False)
maxnorm=0.
for f in ck:
 z=np.load(f);c=z['cells'];g=z['genes'];assert c.shape==(27703,3) and g.shape==(2518,3) and np.isfinite(c).all() and np.isfinite(g).all();maxnorm=max(maxnorm,float(np.linalg.norm(c,axis=1).max()),float(np.linalg.norm(g,axis=1).max()))
W,H=1280,1360;cen=(W/2,715);scale=440/maxnorm;bg='#04090e';white='#edf5f8';muted='#9eb4c2';cyan='#70d2d2';gold='#f6c7bd';fn='/System/Library/Fonts/Supplemental/Arial.ttf';bold='/System/Library/Fonts/Supplemental/Arial Bold.ttf';F=lambda n:ImageFont.truetype(fn,n);B=lambda n:ImageFont.truetype(bold,n);order=np.random.default_rng(23).permutation(A.n_obs)
def project(v,theta=.35):
 x=v[:,0]*np.cos(theta)-v[:,1]*np.sin(theta);y=(v[:,0]*np.sin(theta)+v[:,1]*np.cos(theta))*.42+v[:,2]*.9075;return np.column_stack([cen[0]+x*scale,cen[1]-y*scale])
def blend(a,b,t):return tuple(int(x*(1-t)+y*t) for x,y in zip(a,b))
def heat(c):return blend((12,31,42),(128,224,214),min(1,np.log1p(c)/np.log(11)))
def base(sub):
 im=Image.new('RGB',(W,H),bg);d=ImageDraw.Draw(im);d.text((55,42),'SCENE  |  Native 3D CA1 embedding',font=B(39),fill=white);d.text((58,104),'P60 Seurat CA1 / ProS · 6 WT + 6 cKO mice',font=F(24),fill=muted);d.text((58,139),'27,703 nuclei · 2,518 genes · real training checkpoints',font=F(20),fill=muted);d.line((58,180,1222,180),fill='#3a5561',width=2);d.text((58,204),sub,font=B(25),fill=white);return im,d
def footer(d,text,epoch=None):
 d.line((58,1185,1222,1185),fill='#3a5561',width=2);d.text((58,1204),text,font=F(20),fill=muted);d.text((58,1235),'Model optimization is not cell movement or developmental time.',font=F(19),fill=muted);d.text((58,1275),'01  COUNTS    →    02  FIT    →    03  LAYERS',font=B(19),fill=cyan)
 if epoch is not None:d.text((945,1275),f'{epoch} / 1200',font=B(20),fill=white)
def make_heat():
 im,d=base('01  |  Raw-count matrix excerpt · 96 nuclei × 16 genes');x0,y0=256,400;cw,ch=48,7
 for j,g in enumerate(markers):
  d.line((x0+j*cw+23,y0-14,x0+j*cw+23,y0-3),fill='#d2aba6',width=1); txt=Image.new('RGBA',(160,23));ImageDraw.Draw(txt).text((0,1),g,font=F(19),fill='#e4c4be');txt=txt.rotate(55,expand=True,resample=Image.Resampling.BICUBIC);im.paste(txt,(x0+j*cw-3,y0-156),txt)
  for i in range(96):d.rectangle((x0+j*cw,y0+i*ch,x0+j*cw+cw-3,y0+i*ch+ch-1),fill=heat(matrix[i,j]))
 for i in range(96):y=y0+i*ch+2;d.ellipse((216,y-2,220,y+2),fill=colors[i//48])
 d.line((x0-10,y0+48*ch-3,x0+16*cw,y0+48*ch-3),fill='#77949d',width=2);d.text((58,1100),'Rows: four sampled nuclei per mouse and layer group.',font=F(20),fill=muted);d.text((58,1128),'Teal: log1p(raw UMI count), capped at 10 for display.',font=F(20),fill=muted);d.text((58,1156),'Blue: superficial-associated · orange: deep-associated.',font=F(20),fill=muted);footer(d,'Count matrix uses measured values, not inferred expression.');return np.asarray(im)
heatframe=make_heat();fps=20;writer=imageio.get_writer(str(out/'SCENE_native3d_training.mp4'),fps=fps,codec='libx264',quality=6,macro_block_size=16)
for _ in range(100):writer.append_data(heatframe)
# Exactly recorded native checkpoints; each is held for two frames, without interpolation.
for q,file in enumerate(ck):
 z=np.load(file);c=z['cells'];g=z['genes'];epoch=int(file.stem.split('_')[1]);im,d=base('02  |  Native 3D model training · actual saved checkpoint');d.text((58,248),'Nuclei: existing Seurat clusters  ·  rose points: genes',font=F(20),fill=muted)
 for k,name in enumerate(cat):
  x=58+(k%6)*180;y=285+(k//6)*32;d.ellipse((x,y,x+11,y+11),fill=cc[k]);d.text((x+18,y-3),'C'+name,font=F(19),fill=white)
 xy=project(c);gp=project(g)
 for k,name in enumerate(cat):d.point([tuple(map(int,v)) for v in xy[order[cluster[order]==name]]],fill=cc[k])
 d.point([tuple(map(int,v)) for v in gp],fill='#eeb9b3')
 # Plot only actually recorded losses, never a simulated training history.
 d.text((58,1070),'OBSERVED TRAINING LOSS (log scale)',font=F(17),fill=muted)
 upto=history[history.epoch<=epoch]; vals=np.log(np.maximum(upto.loss.to_numpy(dtype=float),1e-12)); allvals=np.log(np.maximum(history.loss.to_numpy(dtype=float),1e-12))
 if len(vals)>1:
  lo,hi=float(allvals.min()),float(allvals.max()); line=[(410+int(500*(int(e)-1)/1199),1152-int(78*(float(v)-lo)/max(hi-lo,1e-9))) for e,v in zip(upto.epoch,vals)];d.line(line,fill=cyan,width=3);d.ellipse((line[-1][0]-3,line[-1][1]-3,line[-1][0]+3,line[-1][1]+3),fill=cyan)
 footer(d,'True saved training coordinates at epoch '+str(epoch),epoch);arr=np.asarray(im);writer.append_data(arr);writer.append_data(arr)
 if q%30==0:print('rendered checkpoint',epoch,flush=True)
# Final coordinates rotate as one rigid native 3D space. Separate layer views use identical camera and scale.
for k in range(2):
 for n in range(140):
  theta=.35+n*.012;xy=project(finalc,theta);gp=project(finalg,theta);im,d=base('03  |  '+['Superficial-associated','Deep-associated'][k]+' CA1 nuclei');d.text((58,250),f'{np.count_nonzero(layer==L[k]):,} nuclei  ·  common native 3D coordinates and camera',font=F(20),fill=muted);d.ellipse((60,291,73,304),fill=colors[k]);d.text((85,285),'Nuclei',font=F(20),fill=white);d.ellipse((310,291,323,304),fill=gold);d.text((335,285),'Gene positions',font=F(20),fill=white);d.line((580,298,615,298),fill=cyan,width=3);d.text((626,285),'Measured UMI edge',font=F(20),fill=white);d.point([tuple(map(int,v)) for v in xy[order[layer[order]==L[k]]]],fill=colors[k]);d.point([tuple(map(int,v)) for v in gp],fill='#704b4e')
  for e in [v for v in edges if v['layer']==['Superficial-associated','Deep-associated'][k]]:
   x,y=gp[gi[e['gene']]];xx,yy=xy[e['cell_index']];d.line((x,y,xx,yy),fill=cyan,width=3);d.ellipse((xx-4,yy-4,xx+4,yy+4),fill='#e4fff1',outline=cyan,width=2)
  for j,g in enumerate(layer_genes[k]):
   x,y=gp[gi[g]];side=j<2;tx=58 if side else 1100;ty=500+(j%2)*160;ex=tx+150 if side else tx-15;d.line((x,y,ex,ty+15),fill='#c89e9a',width=2);d.ellipse((x-10,y-10,x+10,y+10),fill='#fff0dd',outline=gold,width=3);d.text((tx,ty),g,font=B(27),fill=gold);d.text((tx,ty+34),'gene position',font=F(18),fill=muted)
  d.text((58,1100),'12 example edges: nonzero counts in sampled nuclei.',font=F(20),fill=white);d.text((58,1128),'Not ligand–receptor interactions or physical synapses.',font=F(20),fill=muted);footer(d,'Final fitted native 3D coordinates · camera rotates; positions fixed.',1200);writer.append_data(np.asarray(im))
writer.close(); (out/'COMPLETE.json').write_text(json.dumps({'status':'COMPLETE','nuclei':A.n_obs,'genes':A.n_vars,'training_snapshots':len(ck),'epochs':1200,'dimension':3,'heatmap_shape':[96,16],'expression_edges':len(edges),'movie_fps':fps,'movie_frames':100+len(ck)*2+280,'recorded_checkpoint_frames':'two frames per actual saved checkpoint; no interpolation','layer_display':'same final native space, camera scale and sample labels'},indent=2));print('MOVIE COMPLETE',flush=True)
