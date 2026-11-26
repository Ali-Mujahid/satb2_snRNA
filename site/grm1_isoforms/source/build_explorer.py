"""Review-only interactive export of the completed immutable Grm1 run."""
from pathlib import Path
import csv,json,shutil,html,re,subprocess,gzip,base64
BASE=Path(__file__).resolve().parent;SRC=BASE/'run_v2/results';OUT=BASE/'deliverables';OUT.mkdir(exist_ok=True)
assert json.loads((SRC/'COMPLETE.json').read_text())['n_mice']==12
assert (SRC/'FIGURES_COMPLETE.json').exists()
for p in SRC.glob('*.tsv'):shutil.copy2(p,OUT/p.name)
for name in ['COMPLETE.json','FIGURES_COMPLETE.json','annotation.json','Grm1.reference.gtf']:shutil.copy2(SRC/name,OUT/name)
shutil.copytree(SRC/'figures',OUT/'figures',dirs_exist_ok=True)
for prefix,message in [('03','Diagnostic RNA features plus compatible splicing; not necessarily the alpha/beta-specific terminal junction. Zero does not prove absence.'),('04','Exon-only diagnostic signal may be unprocessed RNA. This is not a count of full-length mature transcripts or receptor protein.')]:
 for p in (OUT/'figures').glob(prefix+'_*.svg'):
  s=p.read_text().replace('Grm1 alpha/beta-specific splice evidence','Grm1 conservative alpha/beta RNA support');h=int(re.search(r'height="(\d+)"',s)[1]);s=s.replace('</svg>',f'<text x="32" y="{h-7}" font-family="Arial,sans-serif" font-size="11" fill="#20353f">{html.escape(message)}</text></svg>');p.write_text(s)
def read(n):return list(csv.DictReader((OUT/n).open(),delimiter='\t'))
def write(n,rs):
 with (OUT/n).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
ann=json.loads((OUT/'annotation.json').read_text());samples=json.loads((OUT/'COMPLETE.json').read_text())['mice'];cov=read('population_mouse_coverage.tsv');present={(r['population'],r['sample']):r for r in cov};pops=sorted({r['population'] for r in cov})
grid=[dict(population=p,sample=s,genotype='WT' if s.startswith('wt') else 'cKO',nuclei=int(present[p,s]['nuclei']) if (p,s) in present else 0,status='SAMPLED' if (p,s) in present else 'NOT_SAMPLED') for p in pops for s in samples];write('all_mouse_population_coverage.tsv',grid)
rows=[]
for file,kind in [('gene_UMI_counts_per_mouse.tsv',None),('alpha_beta_evidence_per_mouse.tsv',None),('isoform_evidence_per_mouse.tsv',None),('exon_UMI_counts_per_mouse.tsv','exon'),('junction_UMI_counts_per_mouse.tsv','junction'),('CDS_UMI_counts_per_mouse.tsv','CDS'),('UTR_UMI_counts_per_mouse.tsv','UTR')]:
 for r in read(file):
  if r['metric']=='family_composition':continue
  z={k:r[k] for k in ['sample','genotype','population','feature','nuclei','total_gene_UMIs','Grm1_locus_UMIs','count','per_1000_nuclei','per_million_library_UMIs','per_1000_Grm1_locus_UMIs']}
  for k in ['nuclei','total_gene_UMIs','Grm1_locus_UMIs','count']:z[k]=int(z[k])
  for k in ['per_1000_nuclei','per_million_library_UMIs','per_1000_Grm1_locus_UMIs']:z[k]=None if z[k] in ['None',''] else float(z[k])
  z['kind']=kind or r['metric'];rows.append(z)
fams={t['name']:t['protein_family'] for t in ann['transcripts']};events=[]
for f in ann['features']:
 members={fams[t] for t in f['transcripts']}
 if f['kind']=='junction' and len(members)==1 and next(iter(members)) in ['alpha','beta']:
  family=next(iter(members));events.append(dict(feature=f['id'],protein_family=family,start=f['start'],end=f['end'],transcripts=';'.join(f['transcripts'])))
eventmap={e['feature']:e for e in events};eventrows=[]
for r in rows:
 if r['kind']=='junction' and r['feature'] in eventmap:
  z=dict(r,**{k:v for k,v in eventmap[r['feature']].items() if k!='feature'});eventrows.append(z)
rows.extend(dict(r,kind='family_junction') for r in eventrows)
write('alpha_beta_specific_junction_counts_per_mouse.tsv',eventrows);write('alpha_beta_specific_junction_definitions.tsv',events)
DATA=dict(samples=samples,annotation=ann,coverage=grid,populations=pops,rows=rows,summary=read('isoform_detection_summary.tsv'),family_summary=read('alpha_beta_detection_summary.tsv'),events=events)
(OUT/'plot_data.json').write_text(json.dumps(DATA,separators=(',',':')))
page=(BASE/'explorer_template.html').read_text();jspdf=(BASE.parent/'satb2_snRNA_github/site/vendor/jspdf.umd.min.js').read_text().replace('</script','<\\/script')
jspdf+=';\n'+(BASE.parent/'satb2_snRNA_github/site/vendor/fflate.min.js').read_text().replace('</script','<\\/script')
gallery=''.join('<p>'+html.escape(p.stem)+' · '+''.join('<a download href="figures/'+p.stem+'.'+ext+'">'+ext.upper()+'</a> ' for ext in ['svg','png','pdf'])+'</p>' for p in sorted((OUT/'figures').glob('*.svg')))
packed={'encoding':'gzip-base64','payload':base64.b64encode(gzip.compress(json.dumps(DATA,separators=(',',':')).encode(),mtime=0)).decode()}
page=page.replace('<div id="gallery"></div>','<div id="gallery">'+gallery+'</div>').replace('DATA_JSON_HERE',json.dumps(packed,separators=(',',':'))).replace('JSPDF_CODE_HERE',jspdf)
(OUT/'index.html').write_text(page)
subprocess.run(['/Users/mujahid357/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node',str(BASE/'workflow/export_figures.cjs'),str(OUT/'figures')],check=True)
print('Exported',len(rows),'interactive plotting rows,',len(events),'alpha/beta-specific junction definitions.')
