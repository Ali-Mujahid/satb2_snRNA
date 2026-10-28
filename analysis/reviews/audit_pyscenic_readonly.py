#!/usr/bin/env python3
"""Independent, streaming arithmetic QA of completed native pySCENIC results.

This is NOT workflow source and never executes inference, AUCell, models, FDR
adjustment for interpretation, RDS/pickle loading, network requests, or native
edits. Recomputed BH and means are QA comparisons to exported native quantities.
Only completed seals are inspected. Output must be outside the native run.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math,sys
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

ENDPOINTS={"CA1_Deep":("CA1_DL",8750),"CA1_Superficial":("CA1_SL",11706),"IN_C1":("IN_C1",1718),"IN_C6":("IN_C6",501),"IN_C8":("IN_C8",468)}
PARTITION="23HJT3LT4_lanes6_7"

def need(ok,message):
 if not ok:raise ValueError(message)
def load(path):return json.loads(Path(path).read_text())
def open_text(path):return gzip.open(path,"rt",encoding="utf-8",newline="") if str(path).endswith('.gz') else open(path,"rt",encoding="utf-8",newline="")
def rows(path):
 with open_text(path) as f:
  r=csv.DictReader(f,delimiter='\t');need(r.fieldnames and len(r.fieldnames)==len(set(r.fieldnames)),f'Invalid header {path}')
  for x in r:need(None not in x,f'Malformed TSV row {path}');yield x

def flag(x):
 if x in [True,'TRUE','True','true','1',1]:return True
 need(x in [False,'FALSE','False','false','0',0],f'Invalid boolean {x}');return False

def number(x,required=True):
 if x is None or x in ('','NA','NaN','nan','null'):
  need(not required,'Required native number missing');return None
 n=float(x);need(math.isfinite(n),'Nonfinite native number');return n

def integer(x):
 n=number(x);need(n>=0 and n.is_integer(),'Invalid native count');return int(n)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def check_receipt(record,path=None):
 p=Path(path or record['path']);need(p.is_file() and not p.is_symlink(),f'Missing/nonregular sealed file: {p}')
 need(p.stat().st_size==record['bytes'] and sha(p)==record['sha256'],f'Sealed bytes changed: {p}')
 return {'path':str(p),'bytes':p.stat().st_size,'sha256':record['sha256']}

def bh_qa(table):
 """Recompute arithmetic BH only; native exported P values are not refitted."""
 usable=[(i,number(r['pvalue'])) for i,r in enumerate(table) if flag(r['tested'])]
 ordered=sorted(usable,key=lambda z:z[1]);m=len(ordered);expected={};last=1.
 for rank in range(m,0,-1):
  idx,p=ordered[rank-1];last=min(last,p*m/rank);expected[idx]=last
 err=0.
 for i,r in enumerate(table):
  q=number(r['association_fdr'],False);alias=number(r.get('padj'),False)
  if i in expected:
   need(q is not None and 0<=q<=1,'Tested native association lacks valid FDR')
   err=max(err,abs(q-expected[i]));need(alias is None or abs(alias-q)<1e-14,'Native padj/FDR alias differs')
  else:need(q is None and alias is None,'Untested native association has manufactured FDR')
 need(err<1e-12,f'Exported FDR differs from arithmetic BH: {err}')
 return {'n_rows':len(table),'n_tested':m,'n_FDR_lt05':sum(number(table[i]['association_fdr'])<.05 for i in expected),'max_abs_BH_QA_error':err,'fdr_families':sorted({r['fdr_family'] for r in table})}

def metadata(path):
 axis=[];mapping={}
 for r in rows(path):
  n=r['nucleus_id'];need(n not in mapping,'Duplicated nucleus metadata ID');axis.append(n)
  mapping[n]={k:r[k] for k in ('endpoint','sample','condition')}
 return axis,mapping

def stream_auc(path,axis,meta,expected_regulons,expected_means=None):
 """Stream actual nucleus scores; retain only 60 x regulon running sums/counts."""
 groups=Counter((meta[n]['endpoint'],meta[n]['sample']) for n in axis)
 sums={key:[0.]*len(expected_regulons) for key in groups};corrections={key:[0.]*len(expected_regulons) for key in groups}
 zero=[0]*len(expected_regulons);minimum=[math.inf]*len(expected_regulons);maximum=[-math.inf]*len(expected_regulons)
 with open_text(path) as f:
  reader=csv.reader(f,delimiter='\t');header=next(reader);need(header[0]=='nucleus_id' and len(header)==len(set(header)),'Invalid AUC header')
  regs=header[1:];need(set(regs)==set(expected_regulons),'Native AUC regulon axis differs from eligibility')
  seen=0
  for i,row in enumerate(reader):
   need(i<len(axis) and len(row)==len(header) and row[0]==axis[i],'Actual nucleus AUC axis/order differs')
   group=meta[row[0]]['endpoint'],meta[row[0]]['sample'];total=sums[group];comp=corrections[group]
   for j,value in enumerate(row[1:]):
    v=number(value);need(0<=v<=1,'Native AUC is outside[0,1]');zero[j]+=v==0;minimum[j]=min(minimum[j],v);maximum[j]=max(maximum[j],v)
    # Kahan sum controls rounding without holding nucleus matrices in memory.
    y=v-comp[j];t=total[j]+y;comp[j]=(t-total[j])-y;total[j]=t
   seen+=1
 need(seen==len(axis),'Missing actual nucleus AUC rows')
 maximum_error=0.;n_compared=0;means={}
 for (ep,sample),total in sums.items():
  for j,reg in enumerate(regs):
   value=total[j]/groups[ep,sample];means[ep,reg,sample]=value
   if expected_means is not None:
    key=ep,reg,sample;need(key in expected_means,'Missing exported actual-nucleus mean')
    got=number(expected_means[key]['aucell_mean']);err=abs(value-got);maximum_error=max(maximum_error,err);n_compared+=1
    need(err<=max(2e-13,2e-12*abs(got)),f'Native mouse AUC mean does not equal all exact nuclei: {key} error={err}')
 if expected_means is not None:need(n_compared==len(expected_means),'Exported means contain unknown axes')
 return {'n_nuclei':seen,'n_scored_regulons':len(regs),'n_mouse_means_compared':n_compared,'max_abs_mouse_mean_error':maximum_error,
  'n_zero_AUC':sum(zero),'n_constant_regulons':sum(a==b for a,b in zip(minimum,maximum)),
  'matrix':str(path),'zero_AUC_preserved':True},means,groups

def completed_stages(native,cfg_sha):
 result={}
 for lin in ('CA1','IN'):
  for stage in ('prepared','grn','ctx','aucell'):
   folder=native/lin/stage;seal=folder/'COMPLETE.json'
   if not seal.is_file():continue
   d=load(seal);need(d['schema']=='satb2-pyscenic-stage/v1' and d['status']=='COMPLETE','Invalid completed checkpoint')
   need(d['lineage']==lin and d['stage']==stage and d['config']['sha256']==cfg_sha,'Checkpoint identity/config differs')
   check_receipt(d['config']);check_receipt(d['source'])
   for rel,rec in d['output_files'].items():check_receipt(rec,folder/rel)
   result[lin+'/'+stage]={'completion_sha256':sha(seal),'summary':d['summary'],'n_sealed_files':len(d['output_files'])}
 return result

def positive_context_audit(native,lin):
 ctx=native/lin/'ctx';nmodules=nregs=nmotifs=0
 for r in rows(ctx/'CANDIDATE_MODULES.tsv'):
  context=json.loads(r['context']);need('activating' in context and 'repressing' not in context,'Nonpositive candidate module');need(integer(r['n_genes'])>=20,'Native candidate module below min20');nmodules+=1
 for r in rows(ctx/'NATIVE_REGULON_COVERAGE.tsv'):
  context=json.loads(r['native_context']);need('activating' in context and 'repressing' not in context and r['native_name'].endswith('(+)'),'Nonpositive final native regulon');nregs+=1
 # ctx's native significance gates are NES/motif annotation, not genotype P.
 for r in rows(ctx/'motif_enrichment_flat.tsv.gz'):
  need(number(r['NES'])>=3.-1e-12,'Native retained motif NES below3')
  context=json.loads(r['Context']);need('activating' in context and 'repressing' not in context,'Nonpositive native motif context')
  need(r.get('Annotation','') not in ('','nan','NaN','NA'),'Missing native motif annotation');nmotifs+=1
 return {'n_positive_candidate_modules':nmodules,'n_native_motif_rows':nmotifs,'n_positive_native_regulons':nregs,
  'native_motif_reference_support_not_direct_binding':True,'genotype_P_not_used_for_discovery':True}

def configured_receipts(cfg):
 out={}
 def walk(x):
  if isinstance(x,dict):
   if {'path','bytes','sha256'}<=set(x):
    key=(x['path'],x['bytes'],x['sha256'])
    if key not in out:out[key]=check_receipt(x)
   else:
    for v in x.values():walk(v)
  elif isinstance(x,list):
   for v in x:walk(v)
 walk(cfg);return list(out.values())

def audit(config,results,allow_partial=False):
 cfg=load(config);native=results/'pyscenic_native';mouse=results/'regulon_mouse';cfg_sha=sha(config)
 need(cfg['schema']=='satb2-pyscenic/v1','Wrong immutable config schema');need(results.resolve()==Path(cfg['run_root']).resolve()/'results','Wrong native run root')
 out={'schema':'satb2-independent-pyscenic-readonly-audit/v1','date_UTC':datetime.now(timezone.utc).isoformat(),'config':check_receipt({'bytes':Path(config).stat().st_size,'sha256':cfg_sha},config),
  'results':str(results),'status':'PARTIAL_CHECKPOINT_AUDIT','scope':'Read-only hashes, exact axes, counts, all-nucleus mean equality and native-statistic arithmetic; no inference or fit',
  'checkpoints':completed_stages(native,cfg_sha),'checklist':['All23,143 exact nuclei retained and raw native AUC axes aligned','Every AUC mean equals all exact nuclei in endpoint/mouse','Six WT +six cKO independent mice, no nucleus pseudo-replication','Native positive motif contexts and unchanged target weights,>=10 measured targets','Joint both-lineage primary BH and separate6WT/3cKO +Satb2 sensitivity families','Exact native Satb2 support, own-gene sensitivity and no forced proxy','All model/evidence limits retained; no direct-binding or causal claim']}
 if not (results/'FINAL_VALIDATION.json').is_file():
  need(allow_partial,'FINAL_VALIDATION is not yet complete; use --allow-partial only for sealed checkpoints')
  return out
 input_receipts=configured_receipts(cfg)
 final=load(results/'FINAL_VALIDATION.json');need(final['configuration_sha256']==cfg_sha,'Final config differs');need(final['schema']=='satb2-pyscenic-regulatory-final-validation/v1' and final['status']=='COMPLETE','Native final validation is not complete')
 for rec in final['native_file_receipts']:check_receipt(rec,results/rec['path'])
 for rec in final['package_file_receipts']:check_receipt(rec,results/rec['path'])
 roster={r['sample']:r['condition'] for r in rows(cfg['samples']['path'])};need(len(roster)==12 and Counter(roster.values())=={'WT':6,'cKO':6},'Wrong independent-mouse roster')
 frozen_cov={(r['endpoint'],r['sample']):r for r in rows(cfg['additional_hash_inputs']['mouse_coverage']['path'])};need(len(frozen_cov)==60,'Wrong frozen endpoint-mouse coverage')
 cov=list(rows(native/'REGULON_COVERAGE.tsv'));regdefs={};maincov={}
 for r in cov:
  key=r['endpoint'],r['regulon_id'];need(key not in maincov,'Duplicate regulon coverage');maincov[key]=r
  lin=r['lineage'];need(r['regulon_id'].startswith(lin+'::') and r['native_name'].endswith('(+)'),'Lineage/native regulon identity changed')
  need(flag(r['eligible'])==(integer(r['n_targets'])>=10),'Target minimum10 changed')
  if r['regulon_id'] in regdefs:need(all(r[k]==regdefs[r['regulon_id']][k] for k in ('lineage','TF','native_name','n_targets','n_background_genes')),'Regulon definition varies across matching endpoints')
  else:regdefs[r['regulon_id']]=r
 means={}
 for r in rows(native/'REGULON_ACTIVITY_PER_MOUSE.tsv'):
  key=r['endpoint'],r['regulon_id'],r['sample'];need(key not in means,'Duplicate native mouse means');means[key]={k:r[k] for k in ('lineage','aucell_mean','n_nuclei','condition')}
 expected={(ep,reg,s) for (ep,reg),c in maincov.items() if flag(c['eligible']) for s in roster};need(set(means)==expected,'Eligible regulons lack exact all12 mice')
 targets=defaultdict(dict);n_target_edges=0
 for r in rows(native/'REGULON_TARGETS.tsv'):
  reg=r['regulon_id'];need(reg in regdefs and r['target'] not in targets[reg],'Unknown/duplicate native target')
  need(number(r['weight'])>0 and flag(r['self_target'])==(r['TF']==r['target']),'Invalid native weight/self flag')
  need(r['edge_origin']==('native_module_added_TF_self_gene' if flag(r['self_target']) else 'native_GRNBoost2_then_cisTarget_leading_edge_union'),'Native target support type changed')
  targets[reg][r['target']]=number(r['weight']);n_target_edges+=1
 need(all(len(targets[reg])==integer(c['n_targets']) for reg,c in regdefs.items()),'Target count differs from native set')
 source_axes={};axis_reports={};all_nuclei=set()
 for lin in ('CA1','IN'):
  original_axis,original_meta=metadata(cfg['lineages'][lin]['cells']['path']);native_axis,native_meta=metadata(native/lin/'aucell/nucleus_axis_metadata.tsv.gz')
  need(original_axis==native_axis and original_meta==native_meta,'Actual scoring metadata differs from full frozen lineage axis')
  need(not all_nuclei.intersection(native_axis),'Nucleus duplicated across lineages');all_nuclei.update(native_axis);source_axes[lin]=(native_axis,native_meta)
  regs={reg for reg,c in regdefs.items() if c['lineage']==lin and flag(c['eligible'])};epmeans={key:r for key,r in means.items() if r['lineage']==lin}
  report,_,groups=stream_auc(native/lin/'aucell/nucleus_AUC.tsv.gz',native_axis,native_meta,regs,epmeans);axis_reports[lin]=report
  for key,n in groups.items():need(n==integer(frozen_cov[key]['n_nuclei']),'Actual raw AUC mouse nucleus count differs')
  for key,r in epmeans.items():need(integer(r['n_nuclei'])==integer(frozen_cov[key[0],key[2]]['n_nuclei']) and r['condition']==roster[key[2]],'Native mean membership changed')
  out.setdefault('motif_context',{})[lin]=positive_context_audit(native,lin)
 need(len(all_nuclei)==23143,'Not all23,143 selected nuclei retained')
 primary=list(rows(mouse/'REGULON_ASSOCIATIONS.tsv'));same=list(rows(mouse/'REGULON_SAME_PARTITION_ASSOCIATIONS.tsv'))
 primary_map={(r['endpoint'],r['regulon_id']):r for r in primary};need(set(primary_map)==set(maincov) and len(primary_map)==len(primary),'Native primary association axes differ')
 for r in primary:
  ep,reg=r['endpoint'],r['regulon_id'];eligible=flag(maincov[ep,reg]['eligible'])
  need(not flag(r['tested']) or eligible and (integer(r['n_WT']),integer(r['n_cKO']))==(6,6),'Invalid primary model coverage')
  if eligible:
   wt=[number(means[ep,reg,s]['aucell_mean']) for s,c in roster.items() if c=='WT'];ko=[number(means[ep,reg,s]['aucell_mean']) for s,c in roster.items() if c=='cKO'];delta=math.fsum(ko)/6-math.fsum(wt)/6
   need(abs(delta-number(r['effect_cKO_minus_WT']))<2e-12,'Native continuous effect differs from independent-mouse mean difference')
 samemap={(r['endpoint'],r['regulon_id']):r for r in same};need(len(samemap)==len(same) and set(samemap)==set(primary_map),'Same-partition model axes differ')
 for r in same:
  ep=r['endpoint'];selected=[s for s in roster if frozen_cov[ep,s]['sequencing_partition']==PARTITION]
  need(Counter(roster[s] for s in selected)=={'WT':6,'cKO':3},'Prespecified partition does not select6WT/3cKO')
  need(not flag(r['tested']) or (integer(r['n_WT']),integer(r['n_cKO']))==(6,3),'Same-partition model uses wrong roster')
  if flag(maincov[ep,r['regulon_id']]['eligible']):
   wt=[number(means[ep,r['regulon_id'],s]['aucell_mean']) for s in selected if roster[s]=='WT'];ko=[number(means[ep,r['regulon_id'],s]['aucell_mean']) for s in selected if roster[s]=='cKO'];delta=math.fsum(ko)/3-math.fsum(wt)/6
   need(abs(delta-number(r['effect_cKO_minus_WT']))<2e-12,'Same-partition effect differs from selected mice')
 out['association_FDR_arithmetic']={'primary':bh_qa(primary),'same_partition':bh_qa(same)}
 if primary:need(out['association_FDR_arithmetic']['primary']['fdr_families']!=out['association_FDR_arithmetic']['same_partition']['fdr_families'],'Primary and same-partition families conflated')
 loo=list(rows(mouse/'REGULON_LOO_MODELS.tsv'));boot=list(rows(mouse/'REGULON_MOUSE_BOOTSTRAP.tsv'));tested_keys={key for key,r in primary_map.items() if flag(r['tested'])}
 loolookup=defaultdict(list)
 for r in loo:loolookup[r['endpoint'],r['regulon_id']].append(r)
 need(set(loolookup)==tested_keys and all(len(v)==12 and {r['omitted_mouse'] for r in v}==set(roster) for v in loolookup.values()),'LOO must omit each actual mouse once')
 for r in loo:need((integer(r['n_WT']),integer(r['n_cKO']))==((5,6) if roster[r['omitted_mouse']]=='WT' else (6,5)),'LOO cohort counts differ from omitted mouse')
 need(len(boot)==len(tested_keys) and all(integer(r['bootstrap_iterations'])==200 and r['resampling_unit']=='independent_mouse_with_replacement_within_each_genotype' for r in boot),'Bootstrap unit or iteration count changed')
 tf_list=set(Path(cfg['references']['tf_names']['path']).read_text().split());sat_audit={}
 nocov=list(rows(native/'SATB2_NO_SELF_TARGET_COVERAGE.tsv'));need(len(nocov)==5 and {r['endpoint'] for r in nocov}==set(ENDPOINTS),'No-self unavailable states not explicit for all5 endpoints')
 nomeanrows=list(rows(native/'SATB2_NO_SELF_TARGET_ACTIVITY_PER_MOUSE.tsv'));nomeans={(r['endpoint'],r['regulon_id'],r['sample']):r for r in nomeanrows};need(len(nomeans)==len(nomeanrows),'Duplicate no-self mouse means')
 for lin in ('CA1','IN'):
  cands=list(rows(native/lin/'prepared/TF_CANDIDATE_COVERAGE.tsv'));need(all(not flag(r['forced_into_TF_list']) for r in cands),'A TF candidate was forced')
  sat=[r for r in cands if r['TF']=='Satb2'];need(len(sat)==1 and flag(sat[0]['official_TF_candidate'])==('Satb2' in tf_list),'Satb2 official TF coverage differs')
  nr=[r for r in nocov if r['lineage']==lin];sreg=lin+'::Satb2(+)';has=sreg in regdefs;eligible=any(flag(r['eligible']) for r in nr);sat_audit[lin]={'candidate':sat[0],'native_positive_motif_regulon_present':has,'main_targets':len(targets[sreg]) if has else None,'self_excluded_available':eligible,'self_excluded_reasons':sorted({r['reason'] for r in nr})}
  for r in nr:
   rid=lin+'::Satb2(+)_NO_SELF_TARGET';need(r['regulon_id']==rid and r['TF']=='Satb2','Wrong no-self regulon identity')
   if flag(r['eligible']):need(has and 'Satb2' in targets[sreg] and len(targets[sreg])-1==integer(r['n_targets'])>=10 and integer(r['n_background_genes'])==integer(regdefs[sreg]['n_background_genes']),'No-self sensitivity changes beyond own gene')
   else:need(not any(key[0]==r['endpoint'] for key in nomeans),'Unavailable self-excluded activity has proxy values')
  if eligible:
   axis,meta=source_axes[lin];sub={key:r for key,r in nomeans.items() if r['lineage']==lin};need(len(sub)==len(nr)*12,'Supported no-self sensitivity lacks all12 mice')
   report,_,_=stream_auc(native/lin/'aucell/SATB2_NO_SELF_TARGET_nucleus_AUC.tsv.gz',axis,meta,{lin+'::Satb2(+)_NO_SELF_TARGET'},sub);sat_audit[lin]['nucleus_AUC_mean_QA']=report
 noassoc=list(rows(mouse/'SATB2_NO_SELF_TARGET_ASSOCIATIONS.tsv'));nosame=list(rows(mouse/'SATB2_NO_SELF_TARGET_SAME_PARTITION_ASSOCIATIONS.tsv'))
 nocovkeys={(r['endpoint'],r['regulon_id']) for r in nocov}
 need(len(noassoc)==len(nosame)==5 and {(r['endpoint'],r['regulon_id']) for r in noassoc}==nocovkeys and {(r['endpoint'],r['regulon_id']) for r in nosame}==nocovkeys,'No-self model axes incomplete')
 for r in noassoc:
  if flag(r['tested']):
   need((integer(r['n_WT']),integer(r['n_cKO']))==(6,6),'No-self primary roster differs')
   wt=[number(nomeans[r['endpoint'],r['regulon_id'],s]['aucell_mean']) for s,c in roster.items() if c=='WT'];ko=[number(nomeans[r['endpoint'],r['regulon_id'],s]['aucell_mean']) for s,c in roster.items() if c=='cKO']
   need(abs(math.fsum(ko)/6-math.fsum(wt)/6-number(r['effect_cKO_minus_WT']))<2e-12,'No-self primary effect differs from native means')
 out['association_FDR_arithmetic']['Satb2_no_self']=bh_qa(noassoc);out['association_FDR_arithmetic']['Satb2_no_self_same_partition']=bh_qa(nosame)
 need(out['association_FDR_arithmetic']['Satb2_no_self']['fdr_families']!=out['association_FDR_arithmetic']['primary']['fdr_families'],'Satb2 sensitivity entered primary BH')
 out.update(status='COMPLETE_READONLY_QA_PASS',n_exact_nuclei=len(all_nuclei),n_native_positive_regulons=len(regdefs),n_score_eligible_regulons=sum(flag(r['eligible']) for r in regdefs.values()),
  n_native_target_edges=n_target_edges,n_native_mouse_means=len(means),n_primary_endpoint_regulon_models=len(primary),n_primary_tested=len(tested_keys),n_LOO_rows=len(loo),n_mouse_bootstrap_summaries=len(boot),
  nucleus_AUC_QA=axis_reports,Satb2_status=sat_audit,all_native_and_package_SHA256_verified=True,n_original_bound_input_files_verified=len(input_receipts),
  primary_significant_by_endpoint=dict(Counter(r['endpoint'] for r in primary if flag(r['tested']) and number(r['association_fdr'])<.05)),
  limitations=['Positive coexpression plus motif-resource support is not measured TF binding or proof of activation/causality','Pooled discovery and mouse comparisons reuse the same cohort; not independent validation','AUC background and target sets vary by lineage; absolute cross-lineage AUC not calibrated','Genotype partly confounds technical partition; same-partition has only3 cKO mice','C8 remains mixed identity; cKO933 has13 selected nuclei; all are scored','Native TF selfgene addition retained, not proof of direct selfbinding','No new gene-level DESeq2/LR fit; original published results preserved'])
 return out

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',required=True,type=Path);parser.add_argument('--results',required=True,type=Path);parser.add_argument('--output',required=True,type=Path);parser.add_argument('--allow-partial',action='store_true');args=parser.parse_args()
 config,results,output=args.config.resolve(),args.results.resolve(),args.output.resolve()
 need(config.is_absolute() and results.is_absolute() and output.is_absolute(),'Use absolute paths');need(not output.is_relative_to(results.parent),'Audit output must be outside entire native run');need(not output.exists(),'Refusing audit overwrite')
 report=audit(config,results,args.allow_partial);output.parent.mkdir(parents=True,exist_ok=True)
 with output.open('x') as f:json.dump(report,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')
 print(json.dumps({'status':report['status'],'output':str(output),'n_completed_checkpoints':len(report['checkpoints']),'n_exact_nuclei':report.get('n_exact_nuclei')},indent=2))
if __name__=='__main__':main()
