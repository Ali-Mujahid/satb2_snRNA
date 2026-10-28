from pathlib import Path
import collections,csv,gzip,hashlib,io,itertools,json,math,zipfile
import numpy as np
import pyarrow as pa
import pyarrow.feather as feather

P=Path('/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026')
R=P/'runs/20261028_pyscenic_v1'
config_path=R/'config/pyscenic.v1.json'
cfg=json.loads(config_path.read_text())
EP={'CA1_Deep':8750,'CA1_Superficial':11706,'IN_C1':1718,'IN_C6':501,'IN_C8':468}
def hrec(p):
 p=Path(p);h=hashlib.sha256()
 with p.open('rb') as stream:
  for b in iter(lambda:stream.read(1024*1024),b''):h.update(b)
 return {'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()}
def verify(r):
 a=hrec(r['path']);assert a['bytes']==r['bytes'] and a['sha256']==r['sha256'];return a
def rows(p):
 p=Path(p);op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt') as stream:yield from csv.DictReader(stream,delimiter='\t')
def flag(r,k):
 assert r[k] in ['True','False'];return r[k]=='True'
official=Path(cfg['references']['tf_names']['path']).read_text().splitlines()
assert len(official)==len(set(official))==1860 and 'Satb2' in official
reference_receipts=[verify(cfg['references']['tf_names']),verify(cfg['references']['motif_annotations'])]
dbs=[]
for record in cfg['references']['rankings']:
 reference_receipts.append(verify(record))
 path=record['path']
 with pa.memory_map(path,'r') as source:
  arrow=pa.ipc.open_file(source);schema=arrow.schema
  index=[n for n in schema.names if str(schema.field(n).type) not in ['int16','int32','int64']]
  assert index==['motifs']
  genes=[n for n in schema.names if n!='motifs']
 motifs=feather.read_table(path,columns=['motifs'],memory_map=True).column('motifs').to_pylist()
 dbs.append({'name':record['name'],'genes':set(genes),'motifs':set(motifs),'n_genes':len(genes),'n_motifs':len(motifs),'Satb2_gene_present':'Satb2' in genes,'gene_types':list(set(str(schema.field(n).type) for n in genes))})
assert dbs[0]['motifs']==dbs[1]['motifs']
ann_TFs=set();ann_ranked_TFs=set();ann_motifs=set();sat_all=[];sat_ranked=[];ann_count=0
for r in rows(cfg['references']['motif_annotations']['path']):
 ann_count+=1
 # Same declared native annotation gate, inspected through numeric metadata only.
 if float(r['motif_similarity_qvalue'])<=.001 and float(r['orthologous_identity'])>=0:
  ann_TFs.add(r['gene_name']);ann_motifs.add(r['#motif_id'])
  if r['#motif_id'] in dbs[0]['motifs']:ann_ranked_TFs.add(r['gene_name'])
  if r['gene_name']=='Satb2':
   sat_all.append(r)
   if r['#motif_id'] in dbs[0]['motifs']:sat_ranked.append(r)
assert dbs[0]['motifs']<=ann_motifs
results={};output_receipts=[]
for lin in ['CA1','IN']:
 out=R/'results/pyscenic_native'/lin/'prepared'
 marker=json.loads((out/'COMPLETE.json').read_text())
 assert marker['status']=='COMPLETE' and marker['stage']=='prepared' and marker['lineage']==lin
 assert marker['config']==hrec(config_path)
 verify(marker['source'])
 assert marker['source']['sha256']=='91a2d534c59fa843f067fc66038ae7eb2039bc9327515aac5d396072d33a6e90'
 for r in marker['output_files'].values():output_receipts.append(verify(r))
 audit=list(rows(out/'gene_background_audit.tsv.gz'));candidate=list(rows(out/'TF_CANDIDATE_COVERAGE.tsv'))
 names=[r['gene_symbol'] for r in rows(out/'expression_genes.tsv')]
 original_gene=list(rows(cfg['lineages'][lin]['genes']['path']))
 assert len(audit)==len(original_gene)==32245
 assert all(all(a[k]==o[k] for k in ['feature_id','gene_symbol','gene_ids','n_source_features','duplicate_symbol','symbol_missing']) for a,o in zip(audit,original_gene))
 expected_names=[r['gene_symbol'] for r in audit if flag(r,'eligible_full_RNA_background')]
 assert names==expected_names and len(names)==len(set(names))
 assert [r['TF'] for r in candidate]==official
 assert all(flag(r,'official_TF_candidate') and not flag(r,'forced_into_TF_list') for r in candidate)
 assert {r['TF'] for r in candidate if flag(r,'eligible_RNA')}==set(official)&set(names)
 detcols=[k for k in audit[0] if k.startswith('detected_fraction__')]
 bad_logic=[]
 for r in audit:
  detection=max(float(r[k]) for k in detcols)>=.05
  eligible=flag(r,'valid_unambiguous_symbol') and detection and flag(r,'normalized_nonconstant')
  if detection!=flag(r,'genotype_blind_detection_eligible') or eligible!=flag(r,'eligible_full_RNA_background'):bad_logic.append(r['feature_id'])
 assert not bad_logic
 # Metadata comparison is streamed: no raw expression matrix or full dense array loaded.
 counts=collections.Counter();mouse_counts=collections.Counter();umi_sums=collections.Counter();n=0;metadata_errors=[];sample_genotype={}
 selected=['nucleus_id','barcode','sample','sample_id','condition','endpoint','lineage','sequencing_partition','library_umi_full_genes']
 for i,(a,b) in enumerate(itertools.zip_longest(rows(cfg['lineages'][lin]['cells']['path']),rows(out/'nucleus_metadata.tsv.gz'))):
  assert a is not None and b is not None
  if any(a[k]!=b[k] for k in selected):metadata_errors.append(i)
  counts[b['endpoint']]+=1;mouse_counts[b['endpoint'],b['sample']]+=1;umi_sums[b['endpoint']]+=int(b['library_umi_full_genes']);sample_genotype[b['sample']]=b['condition'];n+=1
 assert not metadata_errors
 assert all(counts[ep]==EP[ep] for ep in counts) and len(sample_genotype)==12
 assert all(len([s for (ep,s) in mouse_counts if ep==e])==12 for e in counts)
 assert collections.Counter(sample_genotype.values())=={'WT':6,'cKO':6}
 with zipfile.ZipFile(out/'expression_cells_by_genes.npz') as z:
  shape=np.load(io.BytesIO(z.read('shape.npy'))).tolist();fmt=np.load(io.BytesIO(z.read('format.npy'))).item().decode()
  assert shape==[n,len(names)] and fmt=='csc'
 window=list(rows(out/'AUCell_detection_window_audit.tsv.gz'))
 assert len(window)==n and all(r['rank_window_includes_undetected_RNA']=='False' for r in window)
 sat=next(r for r in audit if r['gene_symbol']=='Satb2')
 detection=[]
 for col in detcols:
  ep=col.replace('detected_fraction__','');value=float(sat[col]);detected=round(value*EP[ep]);assert abs(detected/EP[ep]-value)<1e-15
  detection.append({'endpoint':ep,'n_nuclei':EP[ep],'Satb2_detected_nuclei':detected,'detected_fraction':value,'threshold_fraction':.05,'minimum_detected_nuclei_to_pass':math.ceil(EP[ep]*.05),'passes':value>=.05})
 passed_valid=[r for r in audit if flag(r,'valid_unambiguous_symbol')]
 passed_detection=[r for r in passed_valid if flag(r,'genotype_blind_detection_eligible')]
 mapping=[{'ranking_database':d['name'],'eligible_RNA_genes_mapped':len(set(names)&d['genes']),'eligible_RNA_genes_total':len(names),'eligible_RNA_fraction_mapped':len(set(names)&d['genes'])/len(names),'eligible_TF_sources_mapped_as_genes':len(set(names)&set(official)&d['genes'])} for d in dbs]
 results[lin]={'status':'PREPARATION_AUDIT_PASS','native_marker':hrec(out/'COMPLETE.json'),'gene_background':{'all_prepared_features':len(audit),'unambiguous_valid':len(passed_valid),'after_any_endpoint_5percent_detection':len(passed_detection),'after_nonconstant_normalized_RNA':len(names),'reason_counts':dict(collections.Counter(r['reason'] for r in audit))},'TF_candidate_attrition':{'official_unique_candidates':len(candidate),'assayed':sum(flag(r,'assayed_symbol') for r in candidate),'eligible_for_GRN':sum(flag(r,'eligible_RNA') for r in candidate),'excluded':sum(not flag(r,'eligible_RNA') for r in candidate),'reason_counts':dict(collections.Counter(r['reason'] for r in candidate)),'eligible_sources_with_any_native_filtered_annotation':len(set(names)&set(official)&ann_TFs),'eligible_sources_with_annotation_for_ranked_motif':len(set(names)&set(official)&ann_ranked_TFs),'manual_TF_additions':0},'Satb2':{'official_TF_list_member':True,'valid_unambiguous_symbol':flag(sat,'valid_unambiguous_symbol'),'normalized_nonconstant':flag(sat,'normalized_nonconstant'),'eligible':flag(sat,'eligible_full_RNA_background'),'reason':sat['reason'],'per_endpoint_detection':detection},'nucleus_axes':{'n_nuclei':n,'exact_input_metadata_order_and_selected_fields_match':True,'mismatches':metadata_errors,'n_mouse_genotypes':dict(collections.Counter(sample_genotype.values())),'all_endpoints_all12_mice':True,'endpoint_counts':dict(counts),'full_original_library_UMI_sums_by_endpoint':dict(umi_sums),'normalization_full_library_fields_match_original_preparation':True,'sparse_matrix_shape':shape,'sparse_format':fmt,'matrix_shape_header_only_inspection':True,'n_downsampled_nuclei':0,'all_nuclei_detected_RNA_exceeds_AUCell_rank_window':True},'reference_gene_mapping':mapping,'genotype_blind_selection_logic_reconciled_for_all_features':True,'DEG_or_HVG_gene_selection':False}
compat={'reference_receipts':reference_receipts,'official_TF_list_unique_count':len(official),'official_TF_list_contains_Satb2':True,'ranking_resources':[ {k:v for k,v in d.items() if k not in ['genes','motifs']} for d in dbs],'shared_ranked_motifs':len(dbs[0]['motifs']&dbs[1]['motifs']),'annotation_raw_rows':ann_count,'native_filter_unique_TFs':len(ann_TFs),'native_filter_unique_motifs':len(ann_motifs),'native_filter_TFs_with_ranked_motif':len(ann_ranked_TFs),'all_ranked_motifs_have_compatible_annotation':True,'Satb2_filtered_annotation_rows':len(sat_all),'Satb2_ranked_annotation_rows':len(sat_ranked),'Satb2_unique_ranked_annotated_motifs':len({r['#motif_id'] for r in sat_ranked}),'Satb2_ranked_description_classes':dict(collections.Counter('direct_annotation' if 'directly annotated' in r['description'] else 'orthology_or_motif_similarity' for r in sat_ranked)),'scope':'Arrow schema and motif index column only; no full ranking matrix loaded, no scientific inference'}
report={'schema':'satb2-pyscenic-preparation-readonly-audit/v1','audit_date':'2026-10-28','status':'PASS_PREPARATION_ONLY','native_run':str(R),'config':hrec(config_path),'lineages':results,'reference_compatibility':compat,'verified_output_receipts':output_receipts,'verified_native_output_count':len(output_receipts),'scientific_models_run_by_auditor':False,'native_files_modified':False,'limitations':['This is preparation-only review, not proof that GRN or motif pruning has completed or produced valid regulons.','Full original raw count matrices were not reloaded or normalized by this auditor. Library denominators and exact nucleus axes were checked against frozen native metadata; original all-gene reconciliation and tested executed code provide the upstream count evidence.','Sparse expression matrix inspection used the shape/format headers only; file SHA256 seals were independently verified.','Reference file identities and gene/motif annotation compatibility were checked; this is not an independent genome alignment/build audit.','The fixed five-percent gate is a coverage rule; absence of an IN Satb2 regulon cannot establish absence of SATB2 function or protein.','Native motif annotations can include orthology and motif similarity; these are not measured binding evidence in the studied nuclei.']}
out=P/'reviews/20261028_pyscenic_preparation_review.json';out.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
print(json.dumps({'status':report['status'],'outputs_verified':len(output_receipts),'lineages':results,'reference_compatibility':{k:v for k,v in compat.items() if k!='reference_receipts'},'report_path':str(out)},indent=2))
