import unittest,sys,tempfile,json,gzip,importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'workflow';sys.path.insert(0,str(ROOT));import core
spec=importlib.util.spec_from_file_location('remote',ROOT/'remote_counter.py');remote=importlib.util.module_from_spec(spec);spec.loader.exec_module(remote)
for name in ['GENE','annotation','cigar_parts','read_evidence','family_diagnostics']:setattr(remote,name,getattr(core,name))

class CounterTest(unittest.TestCase):
 def test_actual_reference_has_three_coding_and_three_noncoding_models(self):
  p=ROOT.parent/'inputs/Grm1.reference_inspection.gtf';a=core.annotation(p.read_text())
  self.assertEqual(a['gene']['strand'],'-');self.assertEqual(len(a['transcripts']),6)
  self.assertEqual(sum(bool(t['CDS']) for t in a['transcripts']),3)
  self.assertEqual({t['name'] for t in a['transcripts'] if not t['CDS']},{'Grm1-204','Grm1-205','Grm1-206'})
 def test_shared_beta_transcripts_count_once_as_beta_family(self):
  ts=[dict(name='alpha_tx',bit=1,protein_family='alpha',exons=[[100,149],[250,299]],junctions=[[150,249]]),dict(name='beta1',bit=2,protein_family='beta',exons=[[100,149],[350,399]],junctions=[[150,349]]),dict(name='beta2',bit=4,protein_family='beta',exons=[[100,149],[350,399],[450,499]],junctions=[[150,349],[400,449]])]
  fs=[dict(id='shared',kind='exonic_bin',start=100,end=149,transcripts=['alpha_tx','beta1','beta2']),dict(id='uniqueA',kind='exonic_bin',start=250,end=299,transcripts=['alpha_tx']),dict(id='beta_bin',kind='exonic_bin',start=350,end=399,transcripts=['beta1','beta2']),dict(id='JA',kind='junction',start=150,end=249,transcripts=['alpha_tx']),dict(id='JB',kind='junction',start=150,end=349,transcripts=['beta1','beta2'])]
  ann=dict(gene=dict(chrom='chr10',start=100,end=499,strand='-'),transcripts=ts,features=fs);sam=[]
  def row(ub,pos,cigar,cb='AAAA-1',gx=core.GENE,mm=0,xf=1):
   sam.append('\t'.join(['r'+str(len(sam)),'16','chr10',str(pos),'255',cigar,'*','0','0','A'*40,'I'*40,'CB:Z:'+cb,'UB:Z:'+ub,'GX:Z:'+gx,'mm:i:'+str(mm),'xf:i:'+str(xf)]))
  row('A',130,'20M100N20M');row('A',130,'20M100N20M');row('A',130,'20M100N20M',cb='CCCC-1')
  row('B',130,'20M200N20M');row('B',130,'20M200N20M')
  row('shared',110,'30M');row('intron',160,'20M');row('badsplice',260,'20M');row('badsplice',130,'20M90N20M');row('short',143,'7M200N7M')
  row('other',130,'20M100N20M',gx='OTHER');row('mm',130,'20M100N20M',mm=1);row('xf',130,'20M100N20M',xf=2)
  with tempfile.TemporaryDirectory() as td:
   d=Path(td);(d/'bam').write_bytes(b'bam');(d/'index').write_bytes(b'index')
   with gzip.open(d/'called.gz','wt') as f:f.write('AAAA-1\nCCCC-1\n')
   script=d/'samtools';header='@HD\tSO:coordinate\n@SQ\tSN:chr10\tLN:123456\n@RG\tID:r\tSM:wt_1\n'
   script.write_text('#!'+sys.executable+'\nimport sys\nif "--version" in sys.argv:print("samtools synthetic")\nelif "-H" in sys.argv:print('+repr(header)+',end="")\nelif "quickcheck" not in sys.argv:print('+repr('\n'.join(sam)+'\n')+',end="")\n');script.chmod(0o700)
   r=remote.count_sample(dict(sample='wt_1',sample_paths=dict(bam=str(d/'bam'),index=str(d/'index'),called_barcodes=str(d/'called.gz')),annotation=json.loads(json.dumps(ann)),roster={cb:dict(groups=['G'],total_gene_UMIs=100) for cb in ['AAAA-1','CCCC-1']},genotype='WT',samtools_executable=str(script),chromosome_length=123456))
  c={(x['metric'],x['feature']):x['count'] for x in r['counts']}
  self.assertEqual(c['gene','gene'],7)
  self.assertEqual(c['family_splice_supported','alpha'],2)
  self.assertEqual(c['family_splice_supported','beta'],1)
  self.assertEqual(c.get(('splice_supported_unique','beta1'),0),0)
  self.assertEqual(c.get(('splice_supported_unique','beta2'),0),0)
  self.assertEqual(c['feature_UMIs','JB'],1)
  self.assertEqual(sum(v for (m,f),v in c.items() if m=='family_composition'),7)
  self.assertEqual(sum(v for (m,f),v in c.items() if m=='composition'),7)

 def test_negative_strand_reference_translation_input(self):
  gid=core.GENE;lines=['chr10\ttest\tgene\t100\t208\t.\t-\t.\tgene_id "'+gid+'";']
  for i in range(3):
   att='gene_id "'+gid+'"; transcript_id "T'+str(i)+'"; transcript_name "Grm1-20'+str(i+1)+'"; transcript_type "protein_coding";'
   lines.extend(['chr10\ttest\ttranscript\t100\t208\t.\t-\t.\t'+att,'chr10\ttest\texon\t200\t208\t.\t-\t.\t'+att,'chr10\ttest\tCDS\t200\t208\t.\t-\t0\t'+att])
  with tempfile.TemporaryDirectory() as td:
   p=Path(td);(p/'a.gtf').write_text('\n'.join(lines)+'\n');(p/'fa.fai').write_text('chr10\t123456\n');(p/'reference.json').write_text('{}')
   tool=p/'samtools';tool.write_text('#!'+sys.executable+'\nprint(">chr10:200-208\\nAAATTTCAT")\n');tool.chmod(0o700)
   r=remote.read_reference(dict(gtf=str(p/'a.gtf'),fasta=str(p/'fa'),reference_json=str(p/'reference.json'),samtools_executable=str(tool)))
  self.assertTrue(all(t['coding_sequence']=='ATGAAATTT' for t in r['annotation']['transcripts']))

 def test_json_junction_coordinates_and_anchor_boundaries(self):
  t=json.loads(json.dumps(dict(exons=[(100,149),(250,299)],junctions=[(150,249)])))
  blocks,j=core.cigar_parts('8M100N8M',142);self.assertTrue(core.compatible(blocks,j,t))
  self.assertEqual(j,[(150,249,8,8)])
  self.assertFalse(core.compatible(blocks,[(150,248,8,8)],t))
if __name__=='__main__':unittest.main()
