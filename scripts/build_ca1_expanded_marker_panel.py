import json,hashlib
from pathlib import Path
import pandas as pd
sup='Asic2 Ntm Cntnap5a Galntl6 Dcc Kcnq3 AC129186.1 Epha3 Fstl4 Erc2 Grm1 Nav2 Kcnc2 Tenm4 Nlgn1 March1 Robo2 Gpc5 Lingo2 Mdga1 Ak5 Fgf14 Pcdh7 Ryr2 Ptprk Ralyl Kcnab1 Thsd7b Adcy2 Sorbs2 Arl15 Xkr4 Ppp3ca Ano3 Fam19a1 Rora Alcam Phactr1 Nrg2 Plcb1'.split()
deep='Ndst4 Kcnip4 Col11a1 Pde10a Fras1 Epha10 Ntrk2 Pex5l Astn2 Pde4d Sorcs3 Zfp804a Brinp1 Cacna2d3 2010300C02Rik Man1a R3hdm1 Spock3 Prkca Ptprd Opcml Zfpm2 Nr3c2 Epha6 Slc24a2 Malat1 Nr4a1 Srgap1 Stac Arhgap26 Adgrb3 Enox1 Plekhg1 Ap2b1 Hdac9 Rab3c Plekha7 Zbtb16 Pde1a 4930578G10Rik'.split()
src=Path('/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0/results/marker_summaries/gene_cluster_summary_WT.tsv')
out=Path(__file__).parent
d=pd.read_csv(out/'recovered_wt.tsv',sep='\t',dtype={'cluster_id':str})
records=[]; votes={str(x):[] for x in [4,6,7,9,10]}
for side,genes in [('superficial',sup),('deep',deep)]:
 for gene in genes:
  x=d[d.gene_symbol==gene]; ids=x.groupby('cluster_id').size(); valid=len(x)==11 and ids.eq(1).all()
  rec={'gene':gene,'panel':side,'status':'unmapped_or_ambiguous','votes':{}}
  if valid:
   x=x.set_index('cluster_id'); a=float(x.loc[['2','3','8'],'mean_log1p_normalized'].mean());b=float(x.loc[['1','5'],'mean_log1p_normalized'].mean());det=float(x.loc[['1','2','3','5','8'],'fraction_detected'].max());rec.update(sup_anchor=a,deep_anchor=b,anchor_detection=det)
   eligible=(a-b)*(1 if side=='superficial' else -1)>=.1 and det>=.05
   rec['status']='included' if eligible else 'excluded: discordant or weak anchor separation/detection'
   for cl in votes:
    mean=float(x.loc[cl,'mean_log1p_normalized']);frac=float(x.loc[cl,'fraction_detected']);v='superficial' if abs(mean-a)<abs(mean-b) else 'deep' if abs(mean-b)<abs(mean-a) else 'tie';rec['votes'][cl]={'mean':mean,'fraction':frac,'vote':v if eligible else 'excluded'}
    if eligible:votes[cl].append(v)
  records.append(rec)
summary=[]
for cl,vs in votes.items():
 a=vs.count('superficial');b=vs.count('deep');n=a+b;winner='superficial' if a>b else 'deep' if b>a else 'mixed';summary.append({'cluster':cl,'superficial':a,'deep':b,'tie':vs.count('tie'),'winner':winner,'fraction':max(a,b)/n if n else None,'60':winner if n and max(a,b)/n>=.6 and a!=b else 'mixed','70':winner if n and max(a,b)/n>=.7 and a!=b else 'mixed'})
result={'markers':records,'summary':summary,'source_sha256':hashlib.sha256((out/'recovered_wt.tsv').read_bytes()).hexdigest(),'method':'WT pooled cluster means from original count matrix; verified against 385 saved summaries; anchors are equal-weighted means of original superficial C2/C3/C8 and deep C1/C5. Exact unique symbols; anchor separation >=0.1 log1p CP10k in expected direction and >=5% detection in an anchor. Votes compare distance to anchors. Not independent validation.'}
(out/'panel.json').write_text(json.dumps(result,indent=2))
print(json.dumps({'summary':summary,'status':pd.Series([x['status'] for x in records]).value_counts().to_dict()},indent=2))
