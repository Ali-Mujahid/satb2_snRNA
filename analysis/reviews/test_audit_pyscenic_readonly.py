#!/usr/bin/env python3
"""Synthetic arithmetic/streaming checks for the independent read-only audit."""
import csv,gzip,importlib.util,tempfile,unittest
from pathlib import Path
P=Path(__file__).with_name('audit_pyscenic_readonly.py');spec=importlib.util.spec_from_file_location('qa',P);qa=importlib.util.module_from_spec(spec);spec.loader.exec_module(qa)
class AuditTest(unittest.TestCase):
 def test_BH_uses_exported_P_not_refits(self):
  t=[dict(tested='TRUE',pvalue=p,association_fdr=q,padj=q,fdr_family='family') for p,q in [(0.01,.03),(.04,.06),(.20,.20)]]
  r=qa.bh_qa(t);self.assertEqual(r['n_tested'],3);self.assertEqual(r['n_FDR_lt05'],1);self.assertLess(r['max_abs_BH_QA_error'],1e-15)
  t[0]['association_fdr']=.01
  with self.assertRaises(ValueError):qa.bh_qa(t)
 def test_Untested_FDR_stays_missing(self):
  r=qa.bh_qa([dict(tested='FALSE',pvalue='NA',association_fdr='NA',padj='NA',fdr_family='family')]);self.assertEqual(r['n_tested'],0)
  with self.assertRaises(ValueError):qa.bh_qa([dict(tested='FALSE',pvalue='NA',association_fdr=0,padj=0,fdr_family='family')])
 def test_Actual_nucleus_streaming_mean_and_zero(self):
  with tempfile.TemporaryDirectory(prefix='satb2_readonly_QA_SYNTHETIC_') as td:
   p=Path(td)/'scores.tsv.gz';axis=['n1','n2','n3','n4'];meta={n:{'endpoint':'e','sample':'w' if i<2 else 'k','condition':'WT' if i<2 else 'cKO'} for i,n in enumerate(axis)}
   with gzip.open(p,'wt') as f:f.write('nucleus_id\tR1\tR2\nn1\t0\t0.2\nn2\t0.5\t0.4\nn3\t0.1\t0.3\nn4\t0.3\t0.5\n')
   means={('e','R1','w'):{'aucell_mean':.25},('e','R2','w'):{'aucell_mean':.3},('e','R1','k'):{'aucell_mean':.2},('e','R2','k'):{'aucell_mean':.4}}
   r,m,g=qa.stream_auc(p,axis,meta,{'R1','R2'},means);self.assertEqual(r['n_nuclei'],4);self.assertEqual(r['n_zero_AUC'],1);self.assertLess(r['max_abs_mouse_mean_error'],1e-15)
   with self.assertRaisesRegex(ValueError,'axis/order'):qa.stream_auc(p,axis[::-1],meta,{'R1','R2'},means)
   means['e','R1','w']['aucell_mean']=.4
   with self.assertRaisesRegex(ValueError,'does not equal'):qa.stream_auc(p,axis,meta,{'R1','R2'},means)
 def test_Empty_regulon_axis_retains_all_nuclei(self):
  with tempfile.TemporaryDirectory(prefix='satb2_readonly_QA_SYNTHETIC_') as td:
   p=Path(td)/'empty.tsv';p.write_text('nucleus_id\nn1\nn2\n');meta={n:dict(endpoint='e',sample='s',condition='WT') for n in ('n1','n2')}
   r,means,groups=qa.stream_auc(p,['n1','n2'],meta,set(),{});self.assertEqual(r['n_nuclei'],2);self.assertEqual(r['n_scored_regulons'],0);self.assertEqual(means,{})
if __name__=='__main__':unittest.main()
