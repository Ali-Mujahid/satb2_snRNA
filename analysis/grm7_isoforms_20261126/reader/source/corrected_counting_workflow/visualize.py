"""Editable vector figures from source-bound mouse-level count tables."""
import csv,json,math,sys,html,subprocess,statistics
from pathlib import Path
OUT=Path(sys.argv[1]);FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
ann=json.loads((OUT/'annotation.json').read_text());counts=list(csv.DictReader((OUT/'all_counts_long.tsv').open(),delimiter='\t'));cov=list(csv.DictReader((OUT/'population_mouse_coverage.tsv').open(),delimiter='\t'))
samples=json.loads((OUT/'COMPLETE.json').read_text())['mice']
colors=['#146b80','#c06a32','#7e69a5','#71985d','#b65268'];tx=ann['transcripts'];features=ann['features'];lookup={(r['population'],r['sample'],r['metric'],r['feature']):r for r in counts}
POP={'CA1_superficial':'Superficial CA1 · C2 + C3','CA1_deep':'Deep CA1 · C1 + C5','IN_total':'Inhibitory neurons · all 20 clusters','IN_C1':'Inhibitory C1 · PV basket-like','IN_C6':'Inhibitory C6 · parent SST'}
main=['CA1_superficial','CA1_deep','IN_total'];five=main+['IN_C1','IN_C6']
def esc(x):return html.escape(str(x),quote=True)
def text(x,y,s,size=14,fill='#20353f',anchor='start',weight='400'):
 return f'<text x="{x}" y="{y}" fill="{fill}" font-family="Arial,sans-serif" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}">{esc(s)}</text>'
def rect(x,y,w,h,fill,stroke='none'):return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke}"/>'
def line(x,y,a,b,stroke='#c9d4d8',width=1,dash=''):return f'<path d="M{x},{y} L{a},{b}" fill="none" stroke="{stroke}" stroke-width="{width}"'+(f' stroke-dasharray="{dash}"' if dash else '')+'/>'
def circle(x,y,r,c):return f'<circle cx="{x}" cy="{y}" r="{r}" fill="{c}"/>'
def value(pop,sample,metric,feature,norm='count'):
 r=lookup.get((pop,sample,metric,feature));return float(r[norm]) if r and r[norm] not in ['None',''] else 0

def number(v):return str(round(v)) if v>=10 or v==0 else f'{v:.1f}'
def start(w,h,title,sub):return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',rect(0,0,w,h,'white'),text(32,37,title,25,weight='700'),text(32,65,sub,13)]
def save(name,parts):
 (FIG/(name+'.svg')).write_text(''.join(parts)+'</svg>')
def heatmap(name,title,sub,pop,metric,fspec,norm='per_1000_nuclei'):
 n=len(fspec);h=180+n*31;w=1260;x0=390;dx=65;dy=31
 parts=start(w,h,title,sub+' · '+POP[pop]);vals=[[value(pop,s,metric,key,norm) for s in samples] for key,label in fspec];mx=max([1]+[v for r in vals for v in r])
 for j,s in enumerate(samples):
  color='#285e7d' if s.startswith('wt') else '#b56438'
  parts.extend([rect(x0+j*dx,83,dx-2,6,color),text(x0+j*dx+dx/2,109,s.split('_')[1],13,anchor='middle')])
 parts.extend([text(x0+dx*3,78,'WT',12,anchor='middle',weight='700'),text(x0+dx*9,78,'cKO',12,anchor='middle',weight='700')])
 for i,(key,label) in enumerate(fspec):
  y=120+i*dy;parts.append(text(x0-13,y+21,label,12,anchor='end'))
  for j,v in enumerate(vals[i]):
   q=math.log1p(v)/math.log1p(mx);a=[242,246,247];b=[18,91,115];c='#'+''.join(f'{round(x+(z-x)*q):02x}' for x,z in zip(a,b))
   parts.extend([rect(x0+j*dx,y,dx-2,dy-2,c),text(x0+j*dx+dx/2,y+21,number(v),11,fill='white' if q>.58 else '#20353f',anchor='middle')])
 parts.append(line(x0+6*dx-1,83,x0+6*dx-1,120+n*dy,'white',4))
 parts.append(text(32,h-24,'Values: '+('raw corrected CB/UB support' if norm=='count' else 'corrected CB/UB support per 1,000 retained nuclei')+'. Color uses log(1 + value); shared features are not additive.',12))
 save(name,parts)

# Transcript structure with introns compressed uniformly; exons retain their relative lengths.
bins=sorted({(f['start'],f['end']) for f in features if f['kind']=='exonic_bin'});bounds=sorted({v for s,e in bins for v in [s,e+1]});segments=[];x=360
for a,b in zip(bounds,bounds[1:]):
 exonic=any(s<=a<=e for s,e in bins);width=26+26*math.log10(b-a+1) if exonic else 25
 segments.append((a,b-1,x,x+width));x+=width
scale=(1240-360)/(x-360)
def pos(p):
 for a,b,x1,x2 in segments:
  if a<=p<=b:return 360+(x1-360+(p-a)/(b-a+1)*(x2-x1))*scale
 return 1240
parts=start(1320,660,'Grm7 transcript architecture in the actual mapping reference','mm10-2020-A · GENCODE M23 · chr6, positive strand · five transcripts; introns compressed')
for i,t in enumerate(tx):
 y=128+i*88;c=colors[i];parts.extend([text(32,y,t['name'],19,c,weight='700'),text(32,y+23,t['biotype'].replace('_',' '),12),text(32,y+43,t['tail_assignment']+(' · partial CDS' if 'cds_start_NF' in t['tags'] else ' · '+str(t['protein_length'])+' aa'),12),line(pos(t['exons'][0][0]),y,pos(t['exons'][-1][1]),y,c,2)])
 for j,(s,e) in enumerate(t['exons']):
  a=pos(s);b=pos(e);parts.append(rect(a,y-10,max(2,b-a),20,c));parts.append(text((a+b)/2,y-17,'E'+str(j+1),10,anchor='middle'))
 for s,e in t['CDS']:parts.append(rect(pos(s),y-15,max(2,pos(e)-pos(s)),30,c))
 parts.append(text(365,y+39,'C-tail: '+t['C_terminal_30aa'],10,c))
parts.extend([text(32,603,'Thick blocks: CDS · thin blocks: untranslated sequence · lines: introns. Partial CDS and NMD annotations do not establish functional protein isoforms.',12),text(32,627,'Transcript support is assessed from diagnostic features and anchored splice junctions; shared fragments remain unassigned.',12)])
save('01_Grm7_transcript_schematic',parts)
# Main per-animal gene abundance; same normalized axis across all panels.
parts=start(1320,850,'Grm7 abundance in each of the twelve mice','BAM-derived corrected CB/UB support · mice are biological replicates · counts are descriptive')
maxv=max(value(p,s,'gene','gene','per_1000_nuclei') for p in five for s in samples);plotw=1030;x0=230;dx=80
for k,pop in enumerate(five):
 y0=115+k*132;parts.append(text(32,y0,POP[pop],14,weight='700'))
 parts.append(line(x0,y0+78,x0+12*dx,y0+78))
 for j,s in enumerate(samples):
  v=value(pop,s,'gene','gene','per_1000_nuclei');y=y0+78-v/maxv*62;c='#285e7d' if s.startswith('wt') else '#b56438';n=int(value(pop,s,'gene','gene'));parts.extend([line(x0+j*dx+30,y0+78,x0+j*dx+30,y,c,2),circle(x0+j*dx+30,y,4,c),text(x0+j*dx+30,y-9,number(v),10,anchor='middle'),text(x0+j*dx+30,y0+97,s.split('_')[1],10,anchor='middle'),text(x0+j*dx+30,y0+112,'raw '+str(n),9,anchor='middle')])
 parts.append(text(32,y0+30,'UMIs / 1,000 nuclei',11))
parts.extend([text(32,810,'Blue: six WT mice · rust: six cKO mice. Raw counts printed below each animal. Cell numbers and library sizes differ; these are not DESeq2 effects.',12),text(32,833,'Original Cell Ranger gene counts are preserved separately in original_matrix_vs_BAM_gene_counts.tsv.',12)])
save('02_Grm7_gene_counts_per_mouse',parts)
for pop in five:
 heatmap('03_'+pop+'_isoform_splice_support','Grm7 transcript-specific splice evidence','Conservative uniquely compatible, diagnostic and splice-supported UMIs',pop,'splice_supported_unique',[(t['name'],t['name']+' · '+t['tail_assignment']) for t in tx],norm='count')
 heatmap('04_'+pop+'_diagnostic_features','Grm7 isoform-diagnostic feature support','Diagnostic exon segments or unique junctions; exon-only support may reflect pre-mRNA',pop,'diagnostic_feature',[(t['name'],t['name']) for t in tx],norm='count')
for pop in main:
 exons=[f for f in features if f['kind']=='exon'];junc=[f for f in features if f['kind']=='junction']
 heatmap('05_'+pop+'_exon_counts','Grm7 exon feature counts','Each row is a distinct annotated genomic exon interval',pop,'feature_UMIs',[(f['id'],'X'+str(i+1).zfill(2)+' · '+str(f['start'])+'–'+str(f['end'])) for i,f in enumerate(exons)])
 heatmap('06_'+pop+'_junction_counts','Grm7 splice-junction counts','Annotated introns; at least 8 contiguous aligned bases on both sides',pop,'feature_UMIs',[(f['id'],'J'+str(i+1).zfill(2)+' · '+str(f['start'])+'–'+str(f['end'])) for i,f in enumerate(junc)])
# Ambiguity explicitly visualized rather than forcing shared reads into named transcripts.
parts=start(1320,860,'What the reads can and cannot identify','Raw Grm7 UMIs partitioned into conservative isoform assignments and unresolved evidence')
cat=[t['name'] for t in tx]+['shared_exonic','intronic_or_incompatible'];cc=colors+['#aabac0','#e2e7e9'];x0=230;dx=78
for k,pop in enumerate(main):
 y0=120+k*213;parts.append(text(32,y0,POP[pop],16,weight='700'))
 for j,s in enumerate(samples):
  vals=[value(pop,s,'composition',c) for c in cat];total=sum(vals);y=y0+25
  for v,c in zip(vals,cc):
   h=v/total*135 if total else 0;parts.append(rect(x0+j*dx,y,47,h,c));y+=h
  parts.extend([text(x0+j*dx+24,y0+181,s.split('_')[1],11,anchor='middle'),text(x0+j*dx+24,y0+198,str(round(total)),10,anchor='middle')])
for i,(c,label) in enumerate(zip(cc,[t['name'] for t in tx]+['Exonic, unassigned','Intronic / incompatible'])):
 x=32+(i%4)*305;y=790+(i//4)*26;parts.extend([rect(x,y-11,12,12,c),text(x+20,y,label,12)])
save('07_Grm7_isoform_ambiguity',parts)
# Anatomical protein schematic: validated localization distinguished from hypothetical isoform effects.
parts=start(1320,690,'Grm7 RNA isoforms and the presynaptic mGlu7 receptor','A gene-level increase does not identify the protein isoform, membrane abundance or physiological effect')
parts.extend([rect(80,160,1160,24,'#dae4e6'),text(85,140,'Extracellular / synaptic space',16),text(85,214,'Presynaptic terminal cytoplasm',16),text(75,270,'Shared receptor architecture',19,weight='700')])
for j in range(7):
 x=440+j*40;parts.append(rect(x,149,15,48,'#396b7c'));parts.append(text(x+7,132,str(j+1),11,anchor='middle'))
parts.extend([f'<ellipse cx="330" cy="100" rx="80" ry="42" fill="#d5e8eb" stroke="#396b7c" stroke-width="2"/>',text(330,96,'Glutamate-binding',14,anchor='middle'),text(330,114,'domain',14,anchor='middle'),line(410,100,445,151,'#396b7c',3),text(845,176,'Seven transmembrane helices',16)])
for k,t in enumerate([t for t in tx if t['tail_assignment'] in ['mGlu7A-like','mGlu7B-like']]):
 c=colors[tx.index(t)];y=335+k*110;parts.extend([line(700,197,700,y,c,3),line(700,y,1020,y,c,3),text(80,y,t['name']+' · '+t['tail_assignment'],20,c,weight='700'),text(80,y+27,t['C_terminal_30aa'],15,c),text(735,y-13,'Alternative intracellular C-terminal tail',14,c)])
parts.extend([text(80,580,'Testable question: does SATB2 loss alter total Grm7 RNA, transcript processing, or the balance of specific C-terminal variants?',17,weight='700'),text(80,613,'Validate isoform-specific junctions by RT-PCR / targeted sequencing, then assess terminal protein localization and release physiology.',14),text(80,644,'Schematic is conceptual. Nuclear RNA does not establish presynaptic protein trafficking or confirm APP–GRM7 / PENK–GRM7 binding.',13)])
save('08_Grm7_presynaptic_isoform_schematic',parts)
# Portable reader: real SVG files, no invented estimates or external plotting dependencies.
figs=sorted(FIG.glob('*.svg'))
page='<!doctype html><meta charset="utf-8"><title>Grm7 twelve-mouse isoform evidence</title><style>body{font:16px Arial;color:#20353f;max-width:1380px;margin:30px auto;padding:20px}figure{margin:30px 0;border-top:1px solid #ccd8dc;padding-top:20px}img{width:100%;height:auto}a{color:#146b80;margin-right:15px}p{max-width:100ch}</style><h1>Grm7 isoform and transcript-feature evidence</h1><p>All twelve mice · frozen Seurat labels · superficial C2+C3 · deep C1+C5 · inhibitory neurons. Five transcripts are annotated in the mapping GTF. Diagnostic and splice-supported counts are distinguished from shared or unprocessed RNA. A zero means no support under these filters, not proven biological absence.</p>'
for f in figs:
 page+='<figure><h2>'+esc(f.stem.replace('_',' '))+'</h2><a download href="figures/'+f.name+'">SVG</a><a download href="figures/'+f.stem+'.png">PNG</a><a download href="figures/'+f.stem+'.pdf">PDF</a><img src="figures/'+f.name+'"></figure>'
page+='<h2>Count tables and definitions</h2>'
for f in sorted(OUT.glob('*.tsv')):page+='<p><a download href="'+f.name+'">'+esc(f.name)+'</a></p>'
(OUT/'index.html').write_text(page)
subprocess.run(['/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node',str(Path(__file__).with_name('export_figures.cjs')),str(FIG)],check=True)
(OUT/'FIGURES_COMPLETE.json').write_text(json.dumps(dict(status='complete',n_figures=len(figs),formats=['SVG','PNG','PDF'],figures=[f.stem for f in figs]),indent=2)+'\n')
print('Exported',len(figs),'figures in SVG, PNG and PDF.',flush=True)
