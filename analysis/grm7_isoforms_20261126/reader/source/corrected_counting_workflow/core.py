"""Conservative transcript-feature evidence; all coordinates are 1-based inclusive."""
import re,collections
GENE='ENSMUSG00000056755'
CIGAR=re.compile(r'(\d+)([MIDNSHP=X])')
def attrs(s):return dict(re.findall(r'(\w+) "([^"]*)"',s))
def annotation(text):
 tx={};gene=None
 for line in text.splitlines():
  if line.startswith('#'):continue
  x=line.split('\t')
  if len(x)!=9:continue
  a=attrs(x[8])
  if a.get('gene_id')!=GENE:continue
  if x[2]=='gene':gene=dict(chrom=x[0],start=int(x[3]),end=int(x[4]),strand=x[6])
  tid=a.get('transcript_id')
  if not tid:continue
  t=tx.setdefault(tid,dict(id=tid,name=a.get('transcript_name',tid),biotype=a.get('transcript_type'),exons=[],CDS=[],UTR=[],CDS_phases=[],tags=[]))
  if x[2] in ['exon','CDS','UTR']:t[x[2] if x[2]!='exon' else 'exons'].append((int(x[3]),int(x[4])))
  if x[2]=='CDS':t['CDS_phases'].append((int(x[3]),int(x[7])))
  if x[2]=='transcript':t['tags']=re.findall(r'tag "([^"]+)"',x[8])
 assert gene and gene['strand']=='+' and len(tx)==5
 transcripts=sorted(tx.values(),key=lambda t:t['name'])
 for i,t in enumerate(transcripts):
  t['bit']=1<<i;t['exons']=sorted(set(t['exons']));t['CDS']=sorted(set(t['CDS']));t['UTR']=sorted(set(t['UTR']));t['junctions']=[(e[1]+1,f[0]-1) for e,f in zip(t['exons'],t['exons'][1:]) if f[0]>e[1]+1]
 features=[]
 def add(kind,start,end,members,detail=''):
  fid=kind+'_'+str(start)+'_'+str(end)
  features.append(dict(id=fid,kind=kind,start=start,end=end,transcripts=members,detail=detail))
 for i,(s,e) in enumerate(sorted({e for t in transcripts for e in t['exons']}),1):
  members=[t['name'] for t in transcripts if (s,e) in t['exons']]
  detail=';'.join(t['name']+':E'+str(t['exons'].index((s,e))+1) for t in transcripts if (s,e) in t['exons'])
  add('exon',s,e,members,detail)
 for s,e in sorted({j for t in transcripts for j in t['junctions']}):add('junction',s,e,[t['name'] for t in transcripts if (s,e) in t['junctions']])
 for kind in ['CDS','UTR']:
  for s,e in sorted({r for t in transcripts for r in t[kind]}):add(kind,s,e,[t['name'] for t in transcripts if (s,e) in t[kind]])
 bounds=sorted({v for t in transcripts for s,e in t['exons'] for v in [s,e+1]})
 for s,n in zip(bounds,bounds[1:]):
  members=[t['name'] for t in transcripts if any(a<=s and n-1<=b for a,b in t['exons'])]
  if members:add('exonic_bin',s,n-1,members)
 return dict(gene=gene,transcripts=transcripts,features=features)
def cigar_parts(cigar,pos):
 ops=[(int(n),op) for n,op in CIGAR.findall(cigar)]
 if ''.join(str(n)+op for n,op in ops)!=cigar:raise ValueError('Invalid CIGAR')
 blocks=[];junctions=[];p=pos
 for i,(n,op) in enumerate(ops):
  if op in 'M=X':blocks.append((p,p+n-1));p+=n
  elif op=='N':
   l=ops[i-1] if i else (0,'');r=ops[i+1] if i+1<len(ops) else (0,'')
   junctions.append((p,p+n-1,l[0] if l[1] in 'M=X' else 0,r[0] if r[1] in 'M=X' else 0));p+=n
  elif op=='D':p+=n
 return blocks,junctions

def overlap(blocks,s,e):return sum(max(0,min(b,e)-max(a,s)+1) for a,b in blocks)
def compatible(blocks,junctions,t):
 return all(any(s<=a and b<=e for s,e in t['exons']) for a,b in blocks) and all(any(s==a and e==b for a,b in t['junctions']) for s,e,_,_ in junctions)
def read_evidence(blocks,junctions,ann):
 mask=sum(t['bit'] for t in ann['transcripts'] if compatible(blocks,junctions,t))
 hits=[];diagnostic=0
 for f in ann['features']:
  if f['kind']=='junction':yes=any((s,e)==(f['start'],f['end']) and l>=8 and r>=8 for s,e,l,r in junctions)
  else:yes=overlap(blocks,f['start'],f['end'])>0
  if yes:hits.append(f['id'])
  if len(f['transcripts'])==1 and ((f['kind']=='junction' and yes) or (f['kind']=='exonic_bin' and overlap(blocks,f['start'],f['end'])>=20)):
   diagnostic|=next(t['bit'] for t in ann['transcripts'] if t['name']==f['transcripts'][0])
 return mask,hits,diagnostic,any(l>=8 and r>=8 for _,_,l,r in junctions)
