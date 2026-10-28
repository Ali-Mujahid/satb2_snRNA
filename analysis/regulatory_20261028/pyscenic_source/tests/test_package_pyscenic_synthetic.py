#!/usr/bin/env python3
"""Synthetic producer→R→transport validation, never study analysis."""
import csv,gzip,hashlib,importlib.util,json,shutil,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
RWRAPPER=ROOT.parent/'environment/run_R.sh'

def rec(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def js(p,x):p.write_text(json.dumps(x))
def tsv(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 opener=gzip.open if p.suffix=='.gz' else open
 with opener(p,'wt',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def module(path):
 sp=importlib.util.spec_from_file_location('synthetic_package',path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

def fixture(root,plots=False,empty=False):
 src=root/'workflow/scripts';src.mkdir(parents=True)
 for name in ('analyze_regulon_mouse.R','package_pyscenic_results.py'):shutil.copyfile(ROOT/'scripts'/name,src/name)
 (src/'run_pyscenic_native.py').write_text('# SYNTHETIC producer receipt fixture, not an executable workflow\n')
 m=module(src/'package_pyscenic_results.py');run=root/'run';results=run/'results';native=results/'pyscenic_native';native.mkdir(parents=True)
 inputs=root/'synthetic_inputs';inputs.mkdir();samples=[('wt_fake'+str(i),'WT') for i in range(1,7)]+[('ko_fake'+str(i),'cKO') for i in range(1,7)]
 samplefile=inputs/'samples.tsv';tsv(samplefile,[{'sample':s,'condition':c} for s,c in samples],['sample','condition'])
 coverage=[];audit=[];cells_bylineage={'CA1':[],'IN':[]};means={}
 for ep,info in m.ENDPOINTS.items():
  n=info[2];lin=m.lineage(ep);audit.append(dict(endpoint=ep,population_key=info[1],retained_nuclei=n,membership_sha256='a'*64,native_population_membership_sha256='a'*64))
  for j,(s,c) in enumerate(samples):
   count=n//12+(j<n%12);partition='23HJT3LT4_lanes6_7' if c=='WT' or j<9 else 'other';val=.2+[.001,.003,.006,.009,.014,.021][j%6]+(.05 if c=='cKO' else 0)
   coverage.append(dict(endpoint=ep,sample=s,condition=c,n_nuclei=count,sequencing_partition=partition));means[ep,s]=val
   for k in range(count):cells_bylineage[lin].append(dict(nucleus_id=f'SYNTHETIC_{ep}_{s}_{k}',endpoint=ep,sample=s,condition=c,lineage=lin,sequencing_partition=partition))
 covfile=inputs/'coverage.tsv';tsv(covfile,coverage,list(coverage[0]));auditfile=inputs/'membership.tsv';tsv(auditfile,audit,list(audit[0]))
 cfg={'schema':'satb2-pyscenic/v1','analysis_id':'SYNTHETIC_TEST_ONLY','run_root':str(run),'samples':rec(samplefile),'seed':1337,'min_regulon_targets':10,'mouse_associations':{'plots':plots},'lineages':{},'references':{},'additional_hash_inputs':{'frozen_membership_audit':rec(auditfile),'mouse_coverage':rec(covfile)}}
 dummy=inputs/'COMPLETE.json';js(dummy,{'status':'SYNTHETIC_TEST_ONLY'});cfg['prepared_complete']=rec(dummy)
 for lin in ('CA1','IN'):
  cellfile=inputs/(lin+'_cells.tsv');tsv(cellfile,cells_bylineage[lin],list(cells_bylineage[lin][0]));g=inputs/(lin+'_genes.tsv');g.write_text('SYNTHETIC_FULL_BACKGROUND\n');c=inputs/(lin+'_counts.mtx');c.write_text('SYNTHETIC_NOT_EXPRESSION\n')
  cfg['lineages'][lin]={'cells':rec(cellfile),'genes':rec(g),'counts':rec(c)}
 for name in ('tf_names','motif_annotations','ranking1','ranking2'):
  f=inputs/(name+'.fake');f.write_text('SYNTHETIC_REFERENCE_ONLY\n');cfg['references'][name]=rec(f)
 cfg['references']['rankings']=[cfg['references'].pop('ranking1'),cfg['references'].pop('ranking2')]
 config=root/'config.json';js(config,cfg)
 cov=[];scores=[];targets=[]
 basecols=['lineage','endpoint','endpoint_id','source_endpoint_id','regulon_id','native_name','TF','n_targets','n_background_genes','eligible']
 for lin in ('CA1','IN'):
  cand=[dict(lineage=lin,TF='Satb2',official_TF_candidate=False,assayed_symbol=True,eligible_RNA=False,reason='SYNTHETIC_NOT_IN_OFFICIAL_LIST',forced_into_TF_list=False)]
  tsv(native/lin/'prepared/TF_CANDIDATE_COVERAGE.tsv',cand,list(cand[0]))
  for tf,nt in [('FakeTF',10),('FakeLow',9)]:
   reg=lin+'::'+tf+'(+)'
   for j in range(nt):targets.append(dict(lineage=lin,regulon_id=reg,TF=tf,target=tf if j==0 else 'SYNTHETIC_TARGET'+str(j),weight=1.0/(j+1),self_target=j==0))
   for ep,info in m.ENDPOINTS.items():
    if m.lineage(ep)!=lin:continue
    b=dict(lineage=lin,endpoint=ep,endpoint_id=info[0],source_endpoint_id=info[1],regulon_id=reg,native_name=tf+'(+)',TF=tf,n_targets=nt,n_background_genes=100,eligible=nt>=10)
    cov.append({**b,'reason':'ELIGIBLE' if nt>=10 else 'FEWER_THAN10_MEASURED_TARGETS'})
    if nt>=10:
     for s,c in samples:scores.append({**b,'sample':s,'condition':c,'n_nuclei':next(x['n_nuclei'] for x in coverage if x['endpoint']==ep and x['sample']==s),'aucell_mean':means[ep,s],'score_unit':'native_AUCell_AUC'})
  af=native/lin/'aucell/nucleus_AUC.tsv.gz';meta=native/lin/'aucell/nucleus_axis_metadata.tsv.gz';tsv(meta,cells_bylineage[lin],list(cells_bylineage[lin][0]))
  tsv(af,[{'nucleus_id':c['nucleus_id'],lin+'::FakeTF(+)':means[c['endpoint'],c['sample']]} for c in cells_bylineage[lin]],['nucleus_id',lin+'::FakeTF(+)'])
 if empty:
  cov=[];scores=[];targets=[]
  for lin in ('CA1','IN'):
   tsv(native/lin/'aucell/nucleus_AUC.tsv.gz',[{'nucleus_id':x['nucleus_id']} for x in cells_bylineage[lin]],['nucleus_id'])
 tsv(native/'REGULON_COVERAGE.tsv',cov,basecols+['reason']);tsv(native/'REGULON_ACTIVITY_PER_MOUSE.tsv',scores,basecols+['sample','condition','n_nuclei','aucell_mean','score_unit']);tsv(native/'REGULON_TARGETS.tsv',targets,['lineage','regulon_id','TF','target','weight','self_target'])
 statuses=[dict(lineage=lin,status='EMPTY_NO_SUPPORTED_REGULONS' if empty else 'COMPLETE_WITH_REGULONS',n_regulons=0 if empty else 2,n_nuclei=len(cells_bylineage[lin]),n_genes=100) for lin in ('CA1','IN')];tsv(native/'LINEAGE_STATUS.tsv',statuses,list(statuses[0]))
 nc=[]
 for ep,info in m.ENDPOINTS.items():
  lin=m.lineage(ep);nc.append(dict(lineage=lin,endpoint=ep,endpoint_id=info[0],source_endpoint_id=info[1],regulon_id=lin+'::Satb2(+)_NO_SELF_TARGET',native_name='Satb2(+)',TF='Satb2',n_targets=0,n_background_genes=100,eligible=False,reason='NO_SUPPORTED_NATIVE_SATB2_REGULON',network_variant='SATB2_native_regulon_self_target_excluded'))
 tsv(native/'SATB2_NO_SELF_TARGET_COVERAGE.tsv',nc,basecols+['reason','network_variant']);tsv(native/'SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv',[],basecols+['sample','condition','n_nuclei','aucell_mean','score_unit','network_variant'])
 js(native/'PROVENANCE.json',{'schema':'SYNTHETIC_NATIVE_PROVENANCE','config':rec(config),'source':rec(src/'run_pyscenic_native.py')})
 checks={k:True for k in ['positive_motif_regulons_only','full_genotype_neutral_RNA_background','no_DEG_selection','native_nucleus_AUCell','exact_nucleus_axes_verified','all12_mouse_roster','min_measured_targets_10','Satb2_self_target_sensitivity_separate','inputs_config_source_SHA256_unchanged']}
 seal={'schema':'satb2-pyscenic-native-complete/v1','status':'COMPLETE','analysis_id':cfg['analysis_id'],'config':rec(config),'source':rec(src/'run_pyscenic_native.py'),'checks':checks,'output_files':{str(f.relative_to(native)):rec(f) for f in native.rglob('*') if f.is_file()}}
 js(native/'COMPLETE.json',seal)
 out=results/'regulon_mouse';cp=subprocess.run([str(RWRAPPER),str(src/'analyze_regulon_mouse.R'),'--config',str(config),'--input-dir',str(native),'--outdir',str(out)],check=False,capture_output=True,text=True)
 if cp.returncode:raise RuntimeError(cp.stderr+cp.stdout)
 return m,config,results

def reseal(folder,relative):
 p=folder/'COMPLETE.json';d=json.loads(p.read_text());d.get('output_files',d.get('outputs'))[relative]=rec(folder/relative);js(p,d)

class TestSynthetic(unittest.TestCase):
 def test_completed_transport_preserves_means_and_missing(self):
  with tempfile.TemporaryDirectory(prefix='satb2_pyscenic_SYNTHETIC_') as td:
   m,c,r=fixture(Path(td),plots=True);final=m.package(c,r);payload=json.loads((r/'atlas/data.json').read_text())
   self.assertEqual(final['status'],'COMPLETE');self.assertTrue(final['motif_regulon_run']);self.assertEqual(len(payload['methods']),2);self.assertEqual(len(payload['features']),4);self.assertEqual(len(payload['values']),120)
   self.assertEqual(sum(x['value'] is None for x in payload['values']),60);self.assertEqual(len(payload['effects']),10);self.assertEqual(sum(x['tested'] for x in payload['effects']),5);self.assertEqual(len(payload['embeddings']),0)
   self.assertEqual(len(list((r/'regulon_mouse/figures').glob('*.pdf'))),9)
   self.assertEqual({x['family'] for x in payload['methods']},{'motif_regulon'});self.assertTrue(all('binding' in x['edge_basis'] for x in payload['edges']))
 def test_empty_lineages_never_impute_activity(self):
  with tempfile.TemporaryDirectory(prefix='satb2_pyscenic_SYNTHETIC_') as td:
   m,c,r=fixture(Path(td),empty=True);final=m.package(c,r);payload=json.loads((r/'atlas/data.json').read_text())
   self.assertEqual(final['status'],'COMPLETE');self.assertEqual(len(payload['features']),0);self.assertEqual(len(payload['values']),0);self.assertEqual(len(payload['effects']),0);self.assertTrue(all(x['status']=='not_available' for x in payload['methods']))
 def test_changed_native_hash_rejected(self):
  with tempfile.TemporaryDirectory(prefix='satb2_pyscenic_SYNTHETIC_') as td:
   m,c,r=fixture(Path(td));f=r/'pyscenic_native/REGULON_ACTIVITY_PER_MOUSE.tsv';f.write_text(f.read_text()+'\n')
   with self.assertRaisesRegex(ValueError,'Native output changed'):m.package(c,r)
 def test_changed_native_axis_even_if_resealed_rejected(self):
  with tempfile.TemporaryDirectory(prefix='satb2_pyscenic_SYNTHETIC_') as td:
   m,c,r=fixture(Path(td));n=r/'pyscenic_native';f=n/'CA1/aucell/nucleus_AUC.tsv.gz'
   with gzip.open(f,'rt') as h:s=h.read()
   with gzip.open(f,'wt') as h:h.write(s.replace('SYNTHETIC_CA1_Deep_wt_fake1_0','SYNTHETIC_WRONG_AXIS',1))
   reseal(n,str(f.relative_to(n)))
   # Mouse seal binds all native inputs: update only its native file hash to
   # isolate the adapter's independently checked nucleus-axis gate.
   mm=r/'regulon_mouse/COMPLETE.json';d=json.loads(mm.read_text());d['inputs']['native/CA1/aucell/nucleus_AUC.tsv.gz']=rec(f);d['inputs']['native_complete']=rec(n/'COMPLETE.json');js(mm,d)
   with self.assertRaisesRegex(ValueError,'nucleus order'):m.package(c,r)
 def test_changed_R_mean_even_if_resealed_rejected(self):
  with tempfile.TemporaryDirectory(prefix='satb2_pyscenic_SYNTHETIC_') as td:
   m,c,r=fixture(Path(td));d=r/'regulon_mouse';f=d/'REGULON_MOUSE_SCORES.tsv';s=f.read_text();f.write_text(s.replace('0.201','0.299',1));reseal(d,f.name)
   with self.assertRaisesRegex(ValueError,'R changed original native AUCell'):m.package(c,r)
 def test_changed_bound_reference_rejected(self):
  with tempfile.TemporaryDirectory(prefix='satb2_pyscenic_SYNTHETIC_') as td:
   m,c,r=fixture(Path(td));cfg=json.loads(c.read_text());Path(cfg['references']['tf_names']['path']).write_text('CHANGED_REFERENCE\n')
   with self.assertRaisesRegex(ValueError,'Bound input changed'):m.package(c,r)

if __name__=='__main__':unittest.main()
