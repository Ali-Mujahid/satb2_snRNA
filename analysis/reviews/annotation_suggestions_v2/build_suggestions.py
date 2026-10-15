#!/usr/bin/env python3
"""Provisional display labels from completed aggregate/reference evidence; no new RNA fit."""
import ast,csv,json,hashlib,re,collections,math
from pathlib import Path
BASE=Path('/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision')
RESULTS=BASE/'runs/20261014_v2_0_0/results';OUT=BASE/'reviews/annotation_suggestions_v2'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def rows(p):
 with Path(p).open() as f:yield from csv.DictReader(f,delimiter='\t')

def top(x,n=3):return sorted(json.loads(x or '{}').items(),key=lambda z:(-z[1],z[0]))[:n]

def num(x):return None if x in ['', 'NA','nan','NaN'] else float(x)

def clean(x):return re.sub(r'^\d+\s+','',x).replace(' Glut','').replace(' Gaba','').replace(' NN','')

sources=[RESULTS/'annotation_final/cluster_annotations.tsv',RESULTS/'comparison/LABEL_RELIABILITY_REVIEW.tsv',RESULTS/'marker_summaries/gene_cluster_summary.tsv',RESULTS/'marker_summaries/gene_cluster_summary_WT.tsv',RESULTS/'marker_summaries/representative_markers.tsv',BASE/'source/annotation/annotate_clusters.py']
hashes={str(p):sha(p) for p in sources}
for folder,status in [('annotation_final','ANNOTATION_FINAL_COMPLETE_REVIEW_REQUIRED'),('marker_summaries','COMPLETE'),('comparison','EXECUTED_REVIEW_REQUIRED')]:
 p=RESULTS/folder/'COMPLETE.json';d=json.loads(p.read_text());assert d['status']==status;hashes[str(p)]=sha(p)
# Read dictionary literals only; never import or execute a scientific module.
tree=ast.parse((BASE/'source/annotation/annotate_clusters.py').read_text());coarse={}
for node in tree.body:
 if isinstance(node,ast.Assign) and any(isinstance(x,ast.Name) and x.id=='COARSE' for x in node.targets):coarse=ast.literal_eval(node.value)
assert coarse
allrows=[r for r in rows(RESULTS/'comparison/LABEL_RELIABILITY_REVIEW.tsv') if r['branch']=='original_full' and r['cluster_scope']=='original_frozen']
uncertain=[r for r in allrows if r['fine_label_evidence_supported']=='False'];keys={(r['tool'],r['population'],r['cluster_id']) for r in uncertain}
pooled={};wt={}
for filename,target in [('gene_cluster_summary.tsv',pooled),('gene_cluster_summary_WT.tsv',wt)]:
 for r in rows(RESULTS/'marker_summaries'/filename):
  k=(r['tool'],r['population'],r['cluster_id'])
  if r['branch']=='original_full' and r['cluster_scope']=='original_frozen' and k in keys:
   assert r['gene_symbol'] not in target.setdefault(k,{}), 'Duplicate panel symbol: do not collapse stableIDs'
   target[k][r['gene_symbol']]={'gene_id':r['gene_id'],**{name:num(r[name]) for name in ['fraction_detected','mean_log1p_normalized','n_mice_expression_supported','pct_difference_vs_rest','mean_log1p_difference_vs_rest','n_nuclei']}}
reps=collections.defaultdict(list)
for r in rows(RESULTS/'marker_summaries/representative_markers.tsv'):
 k=(r['tool'],r['population'],r['cluster_id'])
 if r['branch']=='original_full' and k in keys:reps[k].append({'gene_id':r['gene_id'],'gene_symbol':r['gene_symbol'],'rank':int(float(r['representative_rank'])),'criteria_pass':r['representative_criteria_pass']=='True'})

suggestions=[]
for r in uncertain:
 k=(r['tool'],r['population'],r['cluster_id']);n=int(float(r['n_nuclei']));ref=top(r['atlas_fractions']);broads=top(r['atlas_broad_fractions']);canonical=[g for g in r['full_cohort_canonical_genes'].split(';') if g];conf=json.loads(r['final_identity_conflicting_genes'] or '{}');status=r['identity_review_status'];audit=n==0;rare=0<n<100
 support={name:[g for g in genes if pooled.get(k,{}).get(g,{}).get('fraction_detected') is not None and pooled[k][g]['fraction_detected']>=.2] for name,genes in coarse.items()}
 support={name:genes for name,genes in support.items() if len(genes)>=2}
 enriched={name:[g for g in genes if pooled[k][g].get('pct_difference_vs_rest') is not None and pooled[k][g]['pct_difference_vs_rest']>=.1] for name,genes in support.items()}
 enriched={name:genes for name,genes in enriched.items() if len(genes)>=2}
 parent={'Excitatory':'glutamatergic-like','Inhibitory':'GABAergic-like','Astrocyte':'astrocyte-like','Oligodendrocyte':'oligodendrocyte-like','OPC':'OPC-like','Microglia':'microglial-like','Vascular':'vascular-associated','Mixed/uncertain':'mixed-lineage RNA'}[r['broad_label']]
 category='reference_mixture';prefix=[];evidence=[]
 if audit:
  parent='reference-associated '+clean(ref[0][0])+'-like' if ref else 'unresolved reference-associated RNA'
  proposed=parent+' (excluded/audit only)';category='audit_only_reference';confidence='LOW_REFERENCE_ONLY_NO_RETAINED_NUCLEI'
 elif 'CONFLICT' in status:
  # Preserve a supported parent and express additional RNA rather than declaring a new cell class.
  if len(canonical)<2:
   refname=ref[0][0] if ref else ''
   if any(t in refname for t in ['Endo','Peri','VLMC','SMC']) and 'Vascular' in enriched:parent='vascular-like'
   elif any(t in refname for t in ['Epend','CHOR']) and 'Ependymal' in enriched:parent='ependymal-like'
   elif 'Neuron' in support:parent='neuronal-like with nonmatching reference'
   else:parent='mixed-lineage reference-associated RNA'
  if r['population']=='CA1_ProS':parent='CA1/ProS-reference-like neuronal RNA'
  others=[]
  for lineage in conf:
   term={'Microglia':'myeloid','Other_immune':'immune','Astrocyte':'astroglial','Oligodendrocyte':'oligodendroglial','OPC':'OPC','Vascular':'vascular','Ependymal':'ependymal'}.get(lineage,lineage)
   if term not in others:others.append(term)
  if not others:others=['nonmatching lineage']
  if parent=='vascular-like' and others==['vascular']:proposed='vascular-like (endothelial reference; mixed original assignment)'
  elif r['broad_label']=='Microglia' and set(conf).issubset({'Other_immune'}):proposed='microglial-like (additional broad immune-lineage RNA)'
  else:proposed=parent+' with '+ '/'.join(others)+'-associated RNA (mixed signature)'
  category='mixed_RNA_signature';confidence='LOW_FINE_IDENTITY;PROVISIONAL_PARENT_ONLY'
 else:
  names=[clean(name) for name,f in ref if f>=.10][:2]
  if r['population']=='CA1_ProS':
   proposed='CA1/ProS-reference-like glutamatergic (weak local marker contrast)';category='weak_local_state';confidence='PROVISIONAL_PARENT_SUPPORTED;FINE_STATE_UNRESOLVED'
  elif r['broad_label']=='Inhibitory':
   if names and all('LSX' in x for x in names):proposed='LSX-reference-associated GABAergic-like (mixed subclasses)'
   elif names:proposed='/'.join(names)+'-reference-like GABAergic (mixed subclasses)' if 'MIXED' in status else names[0]+'-like GABAergic (rare/limited support)'
   else:proposed='GABAergic-like (subclass unresolved)'
   category='rare_reference_resemblance' if rare else 'reference_mixture';confidence='LOW_FINE_IDENTITY;PROVISIONAL_GABAERGIC_PARENT'
  elif r['broad_label']=='Excitatory':
   proposed=('/'.join(names)+'-reference-like glutamatergic (mixed subclasses)') if names else 'glutamatergic-like (subclass unresolved)';confidence='PROVISIONAL_GLUTAMATERGIC_PARENT;FINE_IDENTITY_UNRESOLVED'
  elif r['broad_label']=='Vascular':
   proposed=('/'.join(names)+'-reference-like vascular cluster (mixed reference)') if names else 'vascular-associated RNA (mixed reference)';confidence='PROVISIONAL_VASCULAR_PARENT;FINE_IDENTITY_UNRESOLVED'
  else:
   canonical_lineages=[x for x in enriched if x!='Neuron']
   if canonical_lineages:proposed='mixed '+ '/'.join(canonical_lineages)+'-associated RNA (reference unresolved)'
   elif 'Neuron' in support:
    parent='neuronal-like';proposed='neuronal-like (mixed excitatory/GABAergic reference)'
   else:proposed='mixed neuronal/non-neuronal reference-like RNA (identity unresolved)'
   confidence='PROVISIONAL_PARENT_RESEMBLANCE;FINE_IDENTITY_UNRESOLVED';category='unresolved_parent'
 if rare:proposed+='; rare cluster';confidence+=';RARE_COVERAGE'
 caveats=['Provisional resemblance/display suggestion only; native identity, masks, memberships, categories and models are unchanged.','Reference probabilities and cluster marker summaries do not validate anatomical location, subtype, morphology or protein.','Cluster-average detection cannot distinguish a cell mixture, same-nucleus mixed signal, ambient RNA or multiplets. Detection alone is not lineage specificity; enriched sets require detection>=0.20 and detected-fraction difference>=0.10 versus the existing rest population.']
 if audit:caveats.append('No retained nuclei under the shared main mask; current aggregate marker values are unavailable. Reference-only option stays audit-only/excluded.')
 if rare:caveats.append('Fewer than100 nuclei does not establish absence of a type. Rare IN/CA1 remain included; model eligibility is separate.')
 if r['n_supported_representative_genes']=='0':caveats.append('No representative candidates pass the native specificity/mouse-support rule; a parent-like label does not resolve local state.')
 program=r['layer_program_label'] if r['population']=='CA1_ProS' else r['expression_program_label']
 if r['population']=='CA1_ProS':
  program_display=('mixed/weak superficial–deep-associated RNA programme' if 'Mixed' in program or 'weak' in program else 'unassigned layer-associated programme (conflicting RNA)' if 'not assigned' in program.lower() else program+' (provisional RNA programme)')
  caveats.append('Identity suggestion is separate from the original-full WT-defined layer expression programme. Physical layer stays UNDETERMINED; no programme/model regrouping occurs.')
 elif r['population']=='Inhibitory':
  program_display=program+' (expression programme only)' if program not in ['Not assigned',''] else 'expression programme unresolved'
  caveats.append('Workbook morphology/transmitter programme headings are not validated cell identities; no cholinergic/basket/bistratified/OLM identity is forced.')
 else:program_display=program
 selected=[]
 for gene in list(dict.fromkeys(canonical+(['Pvalb','Sst','Vip','Sncg','Ndnf','Lamp5','Lhx6'] if r['population']=='Inhibitory' else [])+[g for genes in conf.values() for g in genes]+[x['gene_symbol'] for x in reps[k][:4]]+[g for genes in support.values() for g in genes]))[:18]:
  z=pooled.get(k,{}).get(gene,{});w=wt.get(k,{}).get(gene,{})
  if z:selected.append({'gene_id':z['gene_id'],'gene_symbol':gene,'pooled':z,'WT':w,'native_canonical':gene in canonical,'native_conflict_group':[x for x,v in conf.items() if gene in v]})
 key='::'.join(k);suggestions.append({'key':key,'tool':r['tool'],'population':r['population'],'original_cluster_id':r['cluster_id'],'source_branch':'original_full','cluster_scope':'original_frozen','native_label':r['reviewed_label'],'native_identity_status':status,'native_confidence':r['label_confidence'],'native_fine_label_evidence_supported':False,'n_retained_original_nuclei':n,'original_cluster_n_nuclei_before_filter':int(float(r['original_cluster_n_nuclei'])),'keep_for_analysis':r['keep_for_analysis']=='True','exclusion_reason':r['exclusion_reason'],'parent_label':parent,'suggested_label':proposed,'suggestion_category':category,'suggestion_confidence':confidence,'remains_uncertain':True,'native_mixed_or_conflicting_status':status!='REFERENCE_SUPPORTED','reference_top_subclasses':[{'label':x,'fraction':v} for x,v in ref],'reference_top_broad':[{'label':x,'fraction':v} for x,v in broads],'canonical_genes':canonical,'additional_detected_parent_marker_sets':support,'specificity_qualified_marker_sets':enriched,'conflicting_gene_groups':conf,'native_supported_representative_count':int(float(r['n_supported_representative_genes'])),'representative_candidates':reps[k][:7],'marker_evidence':selected,'native_programme_label':program,'programme_display_suggestion':program_display,'physical_layer_assignment':'UNDETERMINED','caveats':caveats,'assignment_is_provisional':True,'native_membership_or_statistical_fit_changed':False})

for p,h in hashes.items():assert sha(p)==h
by_pop=collections.Counter((x['tool'],x['population']) for x in suggestions);explicit=sum(x['native_mixed_or_conflicting_status'] for x in suggestions);retained=sum(x['keep_for_analysis'] for x in suggestions);programme_options=[]
for r in allrows:
 if r['population']!='CA1_ProS':continue
 label=r['layer_program_label'];unresolved=('Mixed' in label or 'weak' in label or 'not assigned' in label.lower());key='::'.join([r['tool'],r['population'],r['cluster_id']])
 text='mixed/weak superficial–deep-associated RNA programme' if 'Mixed' in label or 'weak' in label else 'unassigned layer-associated programme; conflicting identity RNA' if 'not assigned' in label.lower() else label+' (RNA-associated, not anatomical)'
 programme_options.append({'key':key,'tool':r['tool'],'population':r['population'],'original_cluster_id':r['cluster_id'],'native_programme_label':label,'suggested_programme_display_label':text,'assignment_basis':'original-full WT-only expression, carried unchanged to paired branches','superficial_WT_expression_z':num(r['superficial_program_WT_expression_z']),'deep_WT_expression_z':num(r['deep_program_WT_expression_z']),'n_WT_evidence_nuclei':int(float(r['n_WT_program_evidence_nuclei'])),'n_WT_evidence_mice':int(float(r['n_WT_program_evidence_mice'])),'programme_uncertain_before':unresolved,'programme_uncertain_after':unresolved,'physical_layer_assignment':'UNDETERMINED','membership_changed':False,'identity_is_separate_from_programme':True,'confidence':'PROVISIONAL_RNA_PROGRAMME_ONLY','caveat':'Leading WT-score0.3/gap0.5 rules are operational and threshold sensitive; no deep/superficial physical-layer identity or programme regrouping.'})
summary={'n_original_clusters_reviewed':len(allrows),'n_native_fine_supported':len(allrows)-len(uncertain),'n_uncertain_before':len(uncertain),'n_uncertain_after_provisional_suggestions':len(uncertain),'n_explicit_mixed_low_or_conflicting_before':explicit,'n_explicit_mixed_low_or_conflicting_after':explicit,'n_affected_retained_main_clusters':retained,'n_affected_excluded_audit_only_clusters':len(uncertain)-retained,'n_suggestions':len(suggestions),'n_promoted_to_validated_identity':0,'n_CA1_programme_options':len(programme_options),'n_CA1_mixed_weak_or_conflicting_programmes_before':sum(x['programme_uncertain_before'] for x in programme_options),'n_CA1_mixed_weak_or_conflicting_programmes_after':sum(x['programme_uncertain_after'] for x in programme_options),'by_tool_population':[{'tool':t,'population':p,'n_uncertain_before':v,'n_uncertain_after':v} for (t,p),v in sorted(by_pop.items())]}
categories=[{'id':x,'label':x.replace('_',' '),'validation_status':'PROVISIONAL_ONLY'} for x in sorted({x['suggestion_category'] for x in suggestions})]
payload={'schema_version':'satb2-annotation-suggestions-2','scope':'Original Harmony0.4 clusters; nonbinding display/edit suggestions only','source_run':str(RESULTS.parent),'count_summary':summary,'categories':categories,'suggestions':suggestions,'programme_suggestions':programme_options,'programme_by_key':{x['key']:i for i,x in enumerate(programme_options)},'by_key':{x['key']:i for i,x in enumerate(suggestions)},'instructions_for_editor':'Offer explicit optional display-label choices; never overwrite native annotations or imply changed membership/statistics. Carry chosen aliases by original frozen cluster across paired branches. Every proposal retains uncertainty.'}
OUT.mkdir(exist_ok=True)
(OUT/'editor_label_suggestions.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
fields=['key','tool','population','original_cluster_id','n_retained_original_nuclei','original_cluster_n_nuclei_before_filter','keep_for_analysis','exclusion_reason','native_label','native_identity_status','native_confidence','parent_label','suggested_label','suggestion_category','suggestion_confidence','remains_uncertain','native_mixed_or_conflicting_status','native_supported_representative_count','native_programme_label','programme_display_suggestion','physical_layer_assignment','canonical_genes','conflicting_gene_groups','reference_top_subclasses','caveats']
with (OUT/'cluster_label_suggestions.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for x in suggestions:w.writerow({k:json.dumps(x[k],ensure_ascii=False) if isinstance(x[k],(list,dict)) else x[k] for k in fields})
prov={'status':'COMPLETED_READ_ONLY_PROVISIONAL_ANNOTATION_SUGGESTION_REVIEW','source_code_sha256':sha(__file__),'input_sha256':hashes,'count_summary':summary,'output_sha256':{p.name:sha(p) for p in [OUT/'editor_label_suggestions.json',OUT/'cluster_label_suggestions.tsv']},'scientific_fitting_performed':False,'raw_nucleus_RNA_extraction_or_aggregation_performed':False,'native_outputs_modified':False,'native_confidence_or_model_eligibility_upgraded':False}
(OUT/'REVIEW_PROVENANCE.json').write_text(json.dumps(prov,indent=2)+'\n');print(json.dumps(summary,indent=2))
