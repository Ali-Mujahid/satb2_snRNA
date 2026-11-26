"""Review-only exports from completed count tables; no BAM processing or refitting."""
import csv,json,html,shutil,hashlib,subprocess,math,gzip,collections
from pathlib import Path
BASE=Path(__file__).resolve().parent;SRC=BASE/'run_v2/results';OUT=BASE/'deliverables';FIG=OUT/'figures';FIG.mkdir(parents=True,exist_ok=True)
assert json.loads((SRC/'COMPLETE.json').read_text())['n_mice']==12
assert (SRC/'FIGURES_COMPLETE.json').exists()
for p in SRC.glob('*.tsv'):shutil.copy2(p,OUT/p.name)
for name in ['annotation.json','COMPLETE.json','FIGURES_COMPLETE.json','Grm7.reference.gtf']:
 shutil.copy2(SRC/name,OUT/name)
for p in (SRC/'figures').glob('*'):shutil.copy2(p,FIG/p.name)
# Review correction: rendered nCount_RNA is the HVG subtotal, not full library size.
# Keep the frozen native run untouched and correct only reviewed exports.
roster=list(csv.DictReader((BASE/'inputs/nucleus_roster.tsv').open(),delimiter='\t'))
full={}
for lineage in ['CA1_ProS','Inhibitory']:
 p=BASE.parent/'SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0/results/plots/seurat/original_full'/lineage/'original_frozen/rendered_nucleus_metadata.tsv.gz'
 with gzip.open(p,'rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   n=int(float(r['original_nCount_RNA']));assert n==int(float(r['representation_full_library_UMIs']))
   full[r['sample'],r['barcode']]=n
den=collections.Counter();hvg=collections.Counter();nn=collections.Counter()
for r in roster:
 for group in r['groups'].split(';'):
  key=(group,r['sample']);den[key]+=full[r['sample'],r['barcode']];hvg[key]+=int(r['total_gene_UMIs']);nn[key]+=1
audit=[dict(population=p,sample=s,nuclei=nn[p,s],full_source_library_UMIs=v,previous_HVG_subtotal=hvg[p,s]) for (p,s),v in sorted(den.items())]
with (OUT/'library_denominator_audit.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(audit[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(audit)
for p in OUT.glob('*.tsv'):
 with p.open() as f:
  reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames
  if not fields or 'total_gene_UMIs' not in fields or 'population' not in fields or 'sample' not in fields:continue
  records=list(reader)
 for r in records:
  r['total_gene_UMIs']=den[r['population'],r['sample']]
  if 'per_million_library_UMIs' in r:r['per_million_library_UMIs']=int(r['count'])/r['total_gene_UMIs']*1e6
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)
ann=json.loads((SRC/'annotation.json').read_text());samples=json.loads((SRC/'COMPLETE.json').read_text())['mice'];features={f['id']:f for f in ann['features']}
def read(name):return list(csv.DictReader((OUT/name).open(),delimiter='\t'))
coverage=read('population_mouse_coverage.tsv');present={(r['population'],r['sample']):r for r in coverage};pops=sorted({r['population'] for r in coverage})
grid=[]
for p in pops:
 for s in samples:
  r=present.get((p,s));grid.append(dict(population=p,sample=s,genotype='WT' if s.startswith('wt') else 'cKO',nuclei=int(r['nuclei']) if r else 0,status='SAMPLED' if r else 'NOT_SAMPLED',Grm7_locus_UMIs=int(r['Grm7_locus_UMIs']) if r else 'NA'))
def write(name,rows):
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
write('all_mouse_population_coverage.tsv',grid)
allrows=[]
for name in ['gene_UMI_counts_per_mouse.tsv','isoform_evidence_per_mouse.tsv','exon_UMI_counts_per_mouse.tsv','junction_UMI_counts_per_mouse.tsv','CDS_UMI_counts_per_mouse.tsv','UTR_UMI_counts_per_mouse.tsv']:
 for row in read(name):
  r=dict(row);r['count']=int(r['count']);r['nuclei']=int(r['nuclei'])
  for k in ['per_1000_nuclei','per_million_library_UMIs','per_1000_Grm7_locus_UMIs']:r[k]=float(r[k]) if r[k] not in ['None',''] else None
  if name.startswith('exon_'):r['kind']='exon'
  elif name.startswith('junction_'):r['kind']='junction'
  elif name.startswith('CDS_'):r['kind']='CDS'
  elif name.startswith('UTR_'):r['kind']='UTR'
  else:r['kind']=r['metric']
  allrows.append(r)
events=[('A_tail_splice','junction_111495900_111566017','A-type terminal exon skip; Grm7-202'),('B_tail_exon_entry','junction_111495900_111501509','B-type 92-nt exon entry; shared by Grm7-201 and Grm7-205'),('B_201_terminal_exit','junction_111501602_111566017','Long B-type model exit; Grm7-201'),('B_205_terminal_exit','junction_111501602_111502146','Partial B-tail model exit; Grm7-205')]
eventrows=[]
for r in allrows:
 for event,key,desc in events:
  if r['feature']==key and r['kind']=='junction':
   eventrows.append(dict(r,event=event,definition=desc))
write('terminal_splice_event_counts_per_mouse.tsv',eventrows)
for r in eventrows:
 z=dict(r);z['source_feature']=z['feature'];z['feature']=z['event'];z['kind']='terminal_event';allrows.append(z)
# Fix the protein cartoon to show complete receptor forms only; partial CDS stays explicit.
def T(x,y,s,z=15,c='#20353f',anchor='start',weight=400):return f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{z}" fill="{c}" text-anchor="{anchor}" font-weight="{weight}">{html.escape(str(s))}</text>'
def R(x,y,w,h,c):return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{c}"/>'
def L(x,y,a,b,c='#627f89',w=2):return f'<path d="M{x},{y}L{a},{b}" fill="none" stroke="{c}" stroke-width="{w}"/>'
def save_svg(name,w,h,parts):
 (FIG/(name+'.svg')).write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">'+R(0,0,w,h,'white')+''.join(parts)+'</svg>')
parts=[T(35,40,'Grm7 terminal splicing: the key A/B distinction',25,weight=700),T(35,72,'Coordinate-defined RNA events; arrows show transcript order, not synaptic transmission',14)]
for k,(name,color) in enumerate([('A-type · Grm7-202','#146b80'),('B-type · Grm7-201','#bd6934'),('Partial B-tail · Grm7-205','#8e75ac')]):
 y=155+k*140;parts.extend([T(35,y,name,20,color,weight=700),R(365,y-20,140,40,color),R(1020,y-20,175,40,color),T(435,y-30,'Shared exon',13,anchor='middle'),T(1105,y-30,'Terminal exon',13,anchor='middle')])
 if k:
  parts.extend([R(710,y-20,85,40,color),T(752,y-30,'92 nt',13,anchor='middle'),L(505,y,710,y,color),L(795,y,1020,y,color),T(605,y+30,'B exon-entry junction',13,anchor='middle'),T(910,y+30,'Model-specific exit',13,anchor='middle')])
 else:parts.extend([f'<path d="M505,{y} Q760,{y-85} 1020,{y}" fill="none" stroke="{color}" stroke-width="3"/>',T(760,y-50,'A-specific junction skips the 92-nt exon',14,color,anchor='middle')])
parts.extend([T(35,582,'A-type and B-type junction support can be counted without claiming that every complete transcript is reconstructed.',16,weight=700),T(35,613,'Grm7-205 contains a B-like tail but has an incomplete 5′ coding sequence. Its support is not proof of a complete receptor protein.',14),T(35,645,'Introns and exon lengths are schematic. See terminal_splice_event_counts_per_mouse.tsv for exact coordinates and all twelve mice.',13)])
save_svg('09_Grm7_A_B_terminal_splicing',1280,675,parts)
parts=[T(35,40,'Grm7 isoforms and presynaptic mGlu7',25,weight=700),T(35,71,'Two complete receptor models in this GTF; transcript and protein abundance are separate measurements',14),R(80,185,1120,22,'#dbe7e9'),T(85,130,'Extracellular / synaptic space',16),T(85,245,'Presynaptic terminal cytoplasm',16)]
for j in range(7):parts.extend([R(490+j*39,169,14,54,'#366c80'),T(497+j*39,160,j+1,11,anchor='middle')])
parts.extend([f'<ellipse cx="345" cy="120" rx="80" ry="40" fill="#d6e9eb" stroke="#366c80" stroke-width="2"/>',T(345,115,'Glutamate-binding',13,anchor='middle'),T(345,134,'domain',13,anchor='middle'),L(425,130,495,180,'#366c80',3),T(820,201,'Seven transmembrane helices',16)])
for k,name in enumerate(['Grm7-202','Grm7-201']):
 t=next(t for t in ann['transcripts'] if t['name']==name);y=355+k*115;c=['#146b80','#bd6934'][k]
 parts.extend([T(80,y,name+' · '+t['tail_assignment']+' · '+str(t['protein_length'])+' aa',20,c,weight=700),T(80,y+30,t['C_terminal_30aa'],15,c),L(730,223,730,y,c,3),L(730,y,1130,y,c,3),T(770,y-16,'Alternative intracellular C-terminal tail',14,c)])
parts.extend([T(80,565,'The partial Grm7-205 transcript also contains a B-like tail; it is not drawn as a third complete receptor.',14),T(80,600,'RNA support does not establish membrane localization, protein abundance or the direction of a physiological change.',15,weight=700),T(80,632,'Presynaptic localization: Shigemoto et al., J Neurosci 1997. A/B splice forms: Corti et al., Eur J Neurosci 1998.',13)])
save_svg('08_Grm7_presynaptic_isoform_schematic',1280,660,parts)
# Keep A/B colors consistent across the reviewed schematic and composition plot.
for name in ['01_Grm7_transcript_schematic','07_Grm7_isoform_ambiguity']:
 p=FIG/(name+'.svg');s=p.read_text()
 swaps={'#146b80':'#bd6934','#c06a32':'#146b80','#7e69a5':'#6f7e83','#71985d':'#709455','#b65268':'#8e75ac'}
 for i,(a,b) in enumerate(swaps.items()):s=s.replace(a,'COLOR_TOKEN_'+str(i))
 for i,(a,b) in enumerate(swaps.items()):s=s.replace('COLOR_TOKEN_'+str(i),b)
 p.write_text(s)
# Compact portable interactive reader; all normalizations are exported in the source tables.
DATA=dict(samples=samples,annotation=ann,coverage=grid,rows=allrows,populations=pops,summary=read('isoform_detection_summary.tsv'))
(OUT/'plot_data.json').write_text(json.dumps(DATA,separators=(',',':')))
template=(BASE/'explorer_template.html').read_text();jspdf=(BASE.parent/'satb2_snRNA_github/site/vendor/jspdf.umd.min.js').read_text().replace('</script','<\\/script')
gallery=''.join('<p>'+html.escape(f.stem)+' · '+''.join('<a download href="figures/'+f.stem+'.'+ext+'">'+ext.upper()+'</a> ' for ext in ['svg','png','pdf'])+'</p>' for f in sorted(FIG.glob('*.svg')))
template=template.replace('<div id="gallery"></div>','<div id="gallery">'+gallery+'</div>')
page=template.replace('DATA_JSON_HERE',json.dumps(DATA,separators=(',',':')).replace('<','\\u003c')).replace('JSPDF_CODE_HERE',jspdf)
(OUT/'index.html').write_text(page)
subprocess.run(['/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node',str(BASE/'export_reviewed_figures.cjs'),str(FIG)],check=True)
print('Built reviewed interactive explorer with',len(allrows),'plotting rows and',len(grid),'mouse/population coverage records.')
