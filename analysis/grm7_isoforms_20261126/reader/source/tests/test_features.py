import unittest,sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow'))
from core import *
class FeatureTests(unittest.TestCase):
 def test_reference_identity(self):
  a=annotation((Path(__file__).resolve().parents[1]/'inputs/Grm7.reference_inspection.gtf').read_text())
  self.assertEqual([t['name'] for t in a['transcripts']],['Grm7-201','Grm7-202','Grm7-203','Grm7-204','Grm7-205'])
  self.assertEqual(sum(t['biotype']=='protein_coding' for t in a['transcripts']),4)
  for f in a['features']:
   self.assertLessEqual(f['start'],f['end']);self.assertTrue(f['transcripts'])
 def test_cigar_and_anchor_off_by_one(self):
  b,j=cigar_parts('20M100N20M',100)
  self.assertEqual(b,[(100,119),(220,239)]);self.assertEqual(j,[(120,219,20,20)])
  self.assertEqual(overlap(b,119,119),1);self.assertEqual(overlap(b,120,219),0)
 def test_shared_and_diagnostic(self):
  ts=[dict(name='A',bit=1,exons=[(100,149),(250,299)],junctions=[(150,249)]),dict(name='B',bit=2,exons=[(100,149),(350,399)],junctions=[(150,349)])]
  a=dict(transcripts=ts,features=[dict(id='shared',kind='exonic_bin',start=100,end=149,transcripts=['A','B']),dict(id='uniqueA',kind='exonic_bin',start=250,end=299,transcripts=['A']),dict(id='JA',kind='junction',start=150,end=249,transcripts=['A'])])
  m,h,d,s=read_evidence([(110,140)],[],a);self.assertEqual(m,3);self.assertEqual(d,0);self.assertFalse(s)
  a=json.loads(json.dumps(a));m,h,d,s=read_evidence(*cigar_parts('20M100N20M',130),a);self.assertEqual(m,1);self.assertEqual(d,1);self.assertTrue(s);self.assertIn('JA',h)
  m,h,d,s=read_evidence(*cigar_parts('7M100N7M',143),a);self.assertEqual(m,1);self.assertFalse(s);self.assertNotIn('JA',h);self.assertEqual(d,0)
  m,h,d,s=read_evidence([(150,170)],[],a);self.assertEqual(m,0);self.assertFalse(s)
 def test_bad_cigar_rejected(self):
  with self.assertRaises(ValueError):cigar_parts('10Mgarbage',1)
if __name__=='__main__':unittest.main()
