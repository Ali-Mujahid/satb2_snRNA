# Executed on LiSC through a declared read-only SSH task. No remote output files.
import sys,json,base64,subprocess,gzip,hashlib,collections,os
from pathlib import Path

def read_reference(cfg):
 p=Path(cfg['gtf']);h=hashlib.sha256();lines=[]
 with p.open('rb') as f:
  for line in f:
   h.update(line)
   if GENE.encode() in line:lines.append(line.decode())
 text=''.join(lines);ann=annotation(text)
 assert Path(cfg['fasta']+'.fai').is_file()
 regions=sorted({ann['gene']['chrom']+':'+str(s)+'-'+str(e) for t in ann['transcripts'] for s,e in t['CDS']})
 fa=subprocess.check_output([cfg['samtools_executable'],'faidx',cfg['fasta'],*regions],text=True)
 seq={};key=None
 for line in fa.splitlines():
  if line.startswith('>'):key=line[1:];seq[key]=''
  else:seq[key]+=line.strip()
 for t in ann['transcripts']:
  t['coding_sequence']=''.join(seq[ann['gene']['chrom']+':'+str(s)+'-'+str(e)] for s,e in t['CDS'])
 return dict(annotation=ann,gtf_excerpt=text,gtf_sha256=h.hexdigest(),reference=json.loads(Path(cfg['reference_json']).read_text()),reference_path=cfg['reference_json'],gtf_path=cfg['gtf'],fasta_path=cfg['fasta'])

def count_sample(cfg):
 sample=cfg['sample'];s=cfg['sample_paths'];ann=cfg['annotation'];roster=cfg['roster'];sam=cfg['samtools_executable']
 subprocess.run([sam,'quickcheck','-v',s['bam']],check=True,capture_output=True)
 header=subprocess.check_output([sam,'view','-H',s['bam']],text=True)
 rg={x.split('SM:')[1].split('\t')[0] for x in header.splitlines() if x.startswith('@RG') and 'SM:' in x}
 assert rg=={sample},('Sample mismatch',sample,rg)
 assert any(x.startswith('@SQ\tSN:chr6\tLN:149736546') for x in header.splitlines()),'Reference mismatch'
 with gzip.open(s['called_barcodes'],'rt') as f:called={x.strip() for x in f}
 assert set(roster).issubset(called),'Frozen nuclei not present in this BAM call set'
 counts=collections.Counter();read_features=collections.Counter();molecules={};novel=set();diag=collections.Counter()
 region=ann['gene']['chrom']+':'+str(ann['gene']['start'])+'-'+str(ann['gene']['end'])
 cmd=[sam,'view','-q','255','-F','2308','-X',s['bam'],s['index'],region]
 # Pipe stderr is bounded by SAMtools in successful indexed locus extraction.
 p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 for line in p.stdout:
  diag['locus_records']+=1;x=line.rstrip().split('\t');flag=int(x[1])
  tags={v.split(':',2)[0]:v.split(':',2)[2] for v in x[11:] if len(v.split(':',2))==3}
  cb=tags.get('CB');ub=tags.get('UB');xf=int(tags.get('xf','0'))
  if cb not in roster or not ub:continue
  diag['rostered_records']+=1
  if tags.get('GX')!=GENE or int(tags.get('mm','0'))==1 or xf&(2|4|32):continue
  diag['accepted_records']+=1
  blocks,junctions=cigar_parts(x[5],int(x[3]));mask,hits,diagnostic,spliced=read_evidence(blocks,junctions,ann)
  for start,end,l,r in junctions:
   if l<8 or r<8:continue
   if not any(f['kind']=='junction' and (f['start'],f['end'])==(start,end) for f in ann['features']):
    fid='novel_junction_'+str(start)+'_'+str(end);hits.append(fid);novel.add((fid,start,end))
  key=(cb,ub)
  m=molecules.setdefault(key,dict(features=set(),mask=31,has_exonic=False,diagnostic=0,spliced=False,has_intronic=False,incompatible_splice=False))
  m['features'].update(hits);m['diagnostic']|=diagnostic;m['spliced']|=bool(mask and spliced);m['incompatible_splice']|=bool(spliced and not mask)
  if mask:m['mask']&=mask;m['has_exonic']=True
  else:m['has_intronic']=True
  for group in roster[cb]['groups']:
   counts[(group,'records','gene')]+=1
   for fid in hits:read_features[(group,fid)]+=1
 stderr=p.stderr.read();status=p.wait();p.stdout.close();p.stderr.close()
 if status:raise RuntimeError('samtools '+str(status)+' '+stderr[-2000:])
 txnames={t['bit']:t['name'] for t in ann['transcripts']}
 for (cb,ub),m in molecules.items():
  candidate=m['mask'] if m['has_exonic'] else 0
  single=candidate in txnames
  diagnostic=single and bool(m['diagnostic']&candidate)
  strict=diagnostic and m['spliced'] and not m['incompatible_splice']
  category=txnames[candidate] if strict else 'shared_exonic' if m['has_exonic'] and candidate else 'intronic_or_incompatible'
  for group in roster[cb]['groups']:
   counts[(group,'gene','gene')]+=1
   if m['has_exonic']:counts[(group,'exonic_gene','gene')]+=1
   if m['spliced']:counts[(group,'spliced_gene','gene')]+=1
   if m['has_intronic']:counts[(group,'intronic_evidence','gene')]+=1
   counts[(group,'composition',category)]+=1
   names=';'.join(t['name'] for t in ann['transcripts'] if candidate&t['bit']) or 'unresolved'
   counts[(group,'equivalence_class',names)]+=1
   for t in ann['transcripts']:
    if candidate&t['bit']:counts[(group,'compatible_upper_bound',t['name'])]+=1
    if m['diagnostic']&t['bit']:counts[(group,'diagnostic_feature',t['name'])]+=1
    if strict and t['bit']==candidate:counts[(group,'splice_supported_unique',t['name'])]+=1
   for fid in m['features']:counts[(group,'feature_UMIs',fid)]+=1
 for (group,fid),n in read_features.items():counts[(group,'feature_records',fid)]=n
 groupmeta={}
 for cb,meta in roster.items():
  for g in meta['groups']:
   x=groupmeta.setdefault(g,dict(nuclei=0,total_gene_UMIs=0));x['nuclei']+=1;x['total_gene_UMIs']+=meta['total_gene_UMIs']
 stat=Path(s['bam']).stat()
 return dict(sample=sample,genotype=cfg['genotype'],groups=groupmeta,counts=[dict(group=g,metric=k,feature=f,count=n) for (g,k,f),n in sorted(counts.items())],novel_junctions=sorted(novel),diagnostics=dict(diag),provenance=dict(bam=s['bam'],bam_bytes=stat.st_size,bam_mtime_ns=stat.st_mtime_ns,index=s['index'],index_sha256=hashlib.sha256(Path(s['index']).read_bytes()).hexdigest(),header_sha256=hashlib.sha256(header.encode()).hexdigest(),header_RG_samples=sorted(rg),PG=[x for x in header.splitlines() if x.startswith('@PG')],samtools_command=cmd,called_barcodes_sha256=hashlib.sha256(Path(s['called_barcodes']).read_bytes()).hexdigest(),samtools_version=subprocess.check_output([sam,'--version'],text=True).splitlines()[0]))
