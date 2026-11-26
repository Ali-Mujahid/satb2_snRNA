import unittest,sys,tempfile,json,gzip,importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'workflow'
sys.path.insert(0,str(ROOT));import core
spec=importlib.util.spec_from_file_location('remote',ROOT/'remote_counter.py');remote=importlib.util.module_from_spec(spec);spec.loader.exec_module(remote)
remote.GENE=core.GENE;remote.annotation=core.annotation;remote.cigar_parts=core.cigar_parts;remote.read_evidence=core.read_evidence
class CounterTest(unittest.TestCase):
 def test_umi_dedup_compatibility_and_negative_controls(self):
  ts=[dict(name='A',bit=1,exons=[[100,149],[250,299]],junctions=[[150,249]]),dict(name='B',bit=2,exons=[[100,149],[350,399]],junctions=[[150,349]])]
  ann=dict(gene=dict(chrom='chr6',start=100,end=399),transcripts=ts,features=[dict(id='shared',kind='exonic_bin',start=100,end=149,transcripts=['A','B']),dict(id='uniqueA',kind='exonic_bin',start=250,end=299,transcripts=['A']),dict(id='uniqueB',kind='exonic_bin',start=350,end=399,transcripts=['B']),dict(id='JA',kind='junction',start=150,end=249,transcripts=['A']),dict(id='JB',kind='junction',start=150,end=349,transcripts=['B'])])
  sam=[]
  def row(ub,pos,cigar,cb='AAAA-1',gx=core.GENE,mm=0,xf=1):
   sam.append('\t'.join(['read'+str(len(sam)), '0','chr6',str(pos),'255',cigar,'*','0','0','A'*40,'I'*40,'CB:Z:'+cb,'UB:Z:'+ub,'GX:Z:'+gx,'mm:i:'+str(mm),'xf:i:'+str(xf)]))
  row('U1',130,'20M100N20M');row('U1',130,'20M100N20M');row('U1',130,'20M100N20M',cb='CCCC-1')
  row('U2',110,'30M');row('U3',160,'20M');row('U4',260,'20M');row('U4',130,'20M90N20M');row('short',143,'7M100N7M');row('B',130,'20M200N20M')
  row('other',130,'20M100N20M',gx='OTHER');row('mm',130,'20M100N20M',mm=1);row('xf',130,'20M100N20M',xf=2)
  with tempfile.TemporaryDirectory() as td:
   d=Path(td);(d/'bam.bam').write_bytes(b'bam');(d/'bam.bai').write_bytes(b'index')
   with gzip.open(d/'called.gz','wt') as f:f.write('AAAA-1\nCCCC-1\n')
   script=d/'samtools';header='@HD\tSO:coordinate\n@SQ\tSN:chr6\tLN:149736546\n@RG\tID:r\tSM:wt_1\n'
   script.write_text('#!'+sys.executable+'\nimport sys\nif "--version" in sys.argv:print("samtools synthetic")\nelif "-H" in sys.argv:print('+repr(header)+',end="")\nelif "quickcheck" not in sys.argv:print('+repr('\n'.join(sam)+'\n')+',end="")\n');script.chmod(0o700)
   r=remote.count_sample(dict(sample='wt_1',sample_paths=dict(bam=str(d/'bam.bam'),index=str(d/'bam.bai'),called_barcodes=str(d/'called.gz')),annotation=json.loads(json.dumps(ann)),roster={cb:dict(groups=['G'],total_gene_UMIs=100) for cb in ['AAAA-1','CCCC-1']},genotype='WT',samtools_executable=str(script)))
  c={(x['metric'],x['feature']):x['count'] for x in r['counts']}
  self.assertEqual(c['gene','gene'],7);self.assertEqual(c['splice_supported_unique','A'],2);self.assertEqual(c['splice_supported_unique','B'],1)
  self.assertEqual(c['feature_UMIs','JA'],2);self.assertEqual(c['feature_UMIs','JB'],1)
  self.assertEqual(sum(v for (m,f),v in c.items() if m=='composition'),7)
if __name__=='__main__':unittest.main()
