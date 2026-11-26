"""Source-bound SVG plots and portable exports; no model fitting."""
import csv,json,math,html,sys,subprocess
from pathlib import Path
OUT=Path(sys.argv[1]);FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
ann=json.loads((OUT/'annotation.json').read_text());tx=ann['transcripts'];features=ann['features'];samples=json.loads((OUT/'COMPLETE.json').read_text())['mice']
rows=list(csv.DictReader((OUT/'all_counts_long.tsv').open(),delimiter='\t'));lookup={(r['population'],r['sample'],r['metric'],r['feature']):r for r in rows}
POP={'CA1_superficial':'Superficial CA1 · C2+C3','CA1_deep':'Deep CA1 · C1+C5','IN_total':'All inhibitory neurons','IN_C2':'Inhibitory C2 · SST/O-LM-like','IN_C6':'Inhibitory C6 · parent SST'};main=list(POP)
C={'alpha':'#176e83','beta':'#ad6437','noncoding':'#7d8e96'}
def T(x,y,s,z=14,c='#20353f',a='start',w=400):return f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{z}" fill="{c}" text-anchor="{a}" font-weight="{w}">{html.escape(str(s))}</text>'
def R(x,y,w,h,c):return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{c}"/>'
def L(x,y,a,b,c='#a5bbc3',w=1):return f'<path d="M{x},{y}L{a},{b}" fill="none" stroke="{c}" stroke-width="{w}"/>'
def start(w,h,title,sub):return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',R(0,0,w,h,'white'),T(32,38,title,25,w=700),T(32,67,sub,13)]
def save(n,p):(FIG/(n+'.svg')).write_text(''.join(p)+'</svg>')
def v(pop,s,m,f,norm='count'):
 r=lookup.get((pop,s,m,f));return float(r[norm]) if r and r[norm] not in ['None',''] else 0
def num(x):return str(round(x)) if x>=10 or x==0 else f'{x:.1f}'
def heatmap(n,title,pop,metric,fs,norm='count'):
 h=180+31*len(fs);p=start(1320,h,title,POP[pop]+' · '+('raw corrected barcode/UMI counts' if norm=='count' else 'UMIs per 1,000 nuclei'))
 x0=410;dx=69;vals=[[v(pop,s,metric,f,norm) for s in samples] for f,label in fs];mx=max([1]+[x for r in vals for x in r])
 for j,s in enumerate(samples):
  x=x0+j*dx;p.extend([R(x,85,dx-2,6,'#285e7d' if s.startswith('wt') else '#b56438'),T(x+dx/2,111,s.split('_')[1],12,a='middle')])
 for i,(f,label) in enumerate(fs):
  y=125+31*i;p.append(T(x0-14,y+20,label,13,a='end'))
  for j,x in enumerate(vals[i]):
   q=math.log1p(x)/math.log1p(mx);c='#'+''.join(f'{round(a+(b-a)*q):02x}' for a,b in zip([243,247,248],[22,99,124]));p.extend([R(x0+j*dx,y,dx-2,29,c),T(x0+j*dx+dx/2,y+20,num(x),11,'white' if q>.62 else '#20353f',a='middle')])
 p.extend([T(32,h-25,'WT: first six · cKO: last six. Color = log(1 + displayed value); feature rows overlap and cannot be summed.',12)])
 save(n,p)

# Transcript order is left-to-right 5' to 3' on a negative-strand gene.
bins=sorted({(f['start'],f['end']) for f in features if f['kind']=='exonic_bin'});bounds=sorted({z for s,e in bins for z in [s,e+1]});segments=[];x=350
for a,b in zip(bounds,bounds[1:]):
 width=24+20*math.log10(b-a+1) if any(s<=a<=e for s,e in bins) else 25;segments.append((a,b-1,x,x+width));x+=width
scale=895/(x-350)
def pos(p):
 for a,b,x1,x2 in segments:
  if a<=p<=b:return 1245-(x1-350+(p-a)/(b-a+1)*(x2-x1))*scale
 return 350
p=start(1320,876,'Grm1 alpha and beta in the actual mapping reference','mm10 · chr10 negative strand · 5′ → 3′ transcript order · six RNA models; three encode proteins')
for i,t in enumerate(tx):
 y=145+112*i;c=C[t['protein_family']];p.extend([T(32,y,t['name']+' · '+t['protein_family'],20,c,w=700),T(32,y+24,('No annotated CDS' if t['protein_family']=='noncoding' else str(t['protein_length'])+' aa')+' · '+str(len(t['exons']))+' exons',13),L(pos(t['exons'][-1][1]),y,pos(t['exons'][0][0]),y,c,2)])
 for j,(s,e) in enumerate(reversed(t['exons'])):
  left=pos(e);right=pos(s);p.extend([R(left,y-9,max(2,right-left),18,c),T((left+right)/2,y-17,'E'+str(j+1),10,a='middle')])
 for s,e in t['CDS']:p.append(R(pos(e),y-14,max(2,pos(s)-pos(e)),28,c))
 p.append(T(350,y+42,'C-terminal 30 aa: '+t['C_terminal_30aa'],12,c))
p.extend([T(32,829,'Grm1-202 and Grm1-203 encode identical beta proteins. Grm1-204 / -205 / -206 are noncoding competing RNA models.',14,w=700),T(32,855,'Thick bars: CDS · thin bars: untranslated sequence. Counts are RNA evidence; protein levels were not measured.',12)])
save('01_Grm1_transcript_schematic',p)

p=start(1320,850,'Grm1 RNA in each of the twelve mice','Frozen Seurat populations · raw counts and per-nucleus normalization · descriptive, not DESeq2')
mx=max(v(pop,s,'gene','gene','per_1000_nuclei') for pop in main for s in samples)
for k,pop in enumerate(main):
 y=118+132*k;p.extend([T(32,y,POP[pop],15,w=700),T(32,y+28,'UMIs / 1,000 nuclei',11),L(360,y+75,1260,y+75)])
 for j,s in enumerate(samples):
  a=v(pop,s,'gene','gene','per_1000_nuclei');z=y+75-a/mx*62;x=390+j*75;c='#285e7d' if s.startswith('wt') else '#b56438';p.extend([L(x,y+75,x,z,c,2),f'<circle cx="{x}" cy="{z}" r="4" fill="{c}"/>',T(x,z-9,num(a),10,a='middle'),T(x,y+95,s.split('_')[1],10,a='middle'),T(x,y+110,'raw '+str(round(v(pop,s,'gene','gene'))),9,a='middle')])
p.extend([T(32,815,'Blue: six WT · rust: six cKO. Shared y scale. Raw counts are printed below each mouse; sample/library sizes differ.',12),T(32,837,'Nuclear gene UMIs include unprocessed and shared RNA and cannot be read as mature isoform or membrane receptor abundance.',12)])
save('02_Grm1_gene_counts_per_mouse',p)

for pop in main:
 heatmap('03_'+pop+'_alpha_beta_splice','Grm1 alpha/beta-specific splice evidence',pop,'family_splice_supported',[('alpha','Alpha · Grm1-201'),('beta','Beta · Grm1-202 / 203')])
 heatmap('04_'+pop+'_alpha_beta_diagnostic','Grm1 alpha/beta diagnostic feature support',pop,'family_diagnostic_feature',[('alpha','Alpha-associated RNA features'),('beta','Beta-associated RNA features')])
 heatmap('05_'+pop+'_transcript_splice','Grm1 individual transcript-model evidence',pop,'splice_supported_unique',[(t['name'],t['name']+' · '+t['protein_family']) for t in tx])
 for kind,prefix in [('exon','06'),('junction','07')]:
  fs=[f for f in features if f['kind']==kind];heatmap(prefix+'_'+pop+'_'+kind+'_counts','Grm1 '+kind+' feature counts',pop,'feature_UMIs',[(f['id'],str(i+1).zfill(2)+' · '+str(f['start'])+'–'+str(f['end'])) for i,f in enumerate(fs)],'per_1000_nuclei')

p=start(1320,520,'One gene, different mature RNAs, different C-terminal protein tails','Sequence mapping from the matching GTF/FASTA; the schematic does not measure receptor protein')
for i,fam in enumerate(['alpha','beta']):
 y=160+150*i;c=C[fam];p.extend([T(32,y,'mGluR1 '+fam,23,c,w=700),R(355,y-25,430,45,'#7c929c'),T(570,y+3,'Shared receptor sequence · 1–886 aa',15,'white',a='middle')])
 tail=313 if fam=='alpha' else 20;w=400 if fam=='alpha' else 85;p.extend([R(785,y-25,w,45,c),T(785+w/2,y+4,str(tail)+' aa',14,'white',a='middle'),T(355,y+48,'Grm1-201 · long tail · 1199 aa' if fam=='alpha' else 'Grm1-202 / -203 · identical short tail · 906 aa',16,c)])
p.extend([T(32,435,'Alternative splicing changes the mature coding sequence. Total Grm1 counts do not identify alpha/beta balance.',16,w=700),T(32,470,'Ferraguti et al. 1998 Figure 4B measures hippocampal synaptosomal protein with antibodies; it does not assign individual cell types.',13)])
save('08_Grm1_alpha_beta_protein_schematic',p)

figs=sorted(FIG.glob('*.svg'));subprocess.run(['/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node',str(Path(__file__).with_name('export_figures.cjs')),str(FIG)],check=True)
(OUT/'FIGURES_COMPLETE.json').write_text(json.dumps(dict(status='complete',n_figures=len(figs),formats=['SVG','PNG','PDF'],figures=[f.stem for f in figs]),indent=2)+'\n')
print('Exported',len(figs),'Grm1 figures.',flush=True)
