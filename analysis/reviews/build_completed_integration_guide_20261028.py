from pathlib import Path
import csv,json,math,hashlib,collections
P=Path('/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026')
H=P/'runs/20261026_hdwgcna_v1_0_2_export_recovery/results'
T=P/'runs/20261026_collectri_v1/results'
L=Path('/Users/mujahid357/Documents/codex_new/satb2_snRNA_github/communication_data/deg_only_v36')
O=P/'reviews'
sources={}
def receipt(p):
 p=Path(p);return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):
 p=Path(p);sources[str(p)]=receipt(p)
 with p.open() as h:return list(csv.DictReader(h,delimiter='\t'))
def write(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with Path(p).open('w') as h:
  w=csv.DictWriter(h,delimiter='\t',fieldnames=keys,lineterminator='\n');w.writeheader();w.writerows(rows)
 return receipt(p)
def sign(s):return 1 if float(s)>0 else -1 if float(s)<0 else 0
cfg_path=P/'runs/20261026_collectri_v1/config/collectri.v1.json';cfg=json.loads(cfg_path.read_text());sources[str(cfg_path)]=receipt(cfg_path)
selected=[('IN_C6','CA1_Deep','Cntn3','Ptprg'),('IN_C6','CA1_Superficial','Cntn3','Ptprg'),('CA1_Deep','IN_C6','Slit1','Robo1'),('CA1_Deep','IN_C1','Ptn','Ptprz1'),('CA1_Superficial','IN_C1','Ptn','Ptprz1')]
pairs=[r for r in read(L/'interaction_evidence_signed.tsv') if (r['source'],r['target'],r['ligand_complex'],r['receptor_complex']) in selected]
assert len(pairs)==15 and collections.Counter(r['condition'] for r in pairs)=={'pooled':5,'WT':5,'cKO':5}
assert all(r['both_sides_deg_supported']=='True' and r['ligand_all_subunits_deg']=='True' and r['receptor_all_subunits_deg']=='True' for r in pairs)
components={}
for r in pairs:
 for side,gfield,epfield in [('ligand','ligand_complex','source'),('receptor','receptor_complex','target')]:
  k=(r[epfield],r[gfield]);components.setdefault(k,{'role':side,'pair_keys':set()})['pair_keys'].add('|'.join([r['source'],r['target'],r['ligand_complex'],r['receptor_complex']]))
assert len(components)==8
for ep in ['CA1_Deep','CA1_Superficial']:components[(ep,'Satb2')]={'role':'SATB2_gene_context_not_LR_component','pair_keys':set()}
de={}
for ep in sorted({k[0] for k in components}):
 rec=cfg['existing_deseq2'][ep];assert receipt(rec['path'])['sha256']==rec['sha256'];rr=read(rec['path']);idx=collections.defaultdict(list)
 for r in rr:idx[r['gene']].append(r)
 for e,g in components:
  if e==ep:assert len(idx[g])==1;de[e,g]=idx[g][0]
module={};ranks={};lo={}
for lin in ['CA1','IN']:
 module[lin]={r['gene_symbol']:r for r in read(H/f'hdwgcna/{lin}/pooled/module_membership.tsv')}
 ranks[lin]={r['gene_symbol']:r for r in read(H/f'hdwgcna/{lin}/pooled/all_gene_hub_ranks.tsv')}
 for r in read(H/f'hdwgcna/{lin}/pooled/module_genotype_leave_one_mouse_out.tsv'):lo.setdefault((r['endpoint'],r['module']),[]).append(r)
primary={(r['endpoint'],r['module']):r for r in read(H/'hdwgcna/pooled_module_genotype_associations.tsv')}
part={(r['endpoint'],r['module']):r for r in read(H/'hdwgcna/pooled_module_genotype_same_partition.tsv')}
link={(r['endpoint'],r['gene_symbol']):r for r in read(H/'enrichment/STRICT_LR_AND_SATB2_MODULE_LINKS.tsv')}
tfa={(r['endpoint'],r['TF']):r for r in read(T/'tf_mouse/TF_ASSOCIATIONS.tsv')}
tfp={(r['endpoint'],r['TF']):r for r in read(T/'tf_mouse/TF_SAME_PARTITION_ASSOCIATIONS.tsv')}
sourceRNA={};prior={}
for ep in sorted({k[0] for k in components}):
 prior[ep]=read(T/f'tf_activity/{ep}/mouse_activity_used_signed_prior_edges.tsv')
 sourceRNA[ep]=collections.defaultdict(list)
 for r in read(T/f'tf_activity/{ep}/all_gene_activity_background_audit.tsv'):sourceRNA[ep][r['gene_symbol']].append(r)
component_rows=[];prior_rows=[]
for ep,g in components:
 lin='CA1' if ep.startswith('CA1') else 'IN';m=module[lin][g];mod=m['module'];a=primary[ep,mod];s=part[ep,mod];d=de[ep,g];native_link=link[ep,g]
 assert native_link['module']==mod and float(native_link['DESeq2_log2FC'])==float(d['log2FoldChange']) and float(native_link['DESeq2_padj'])==float(d['padj'])
 ll=lo[ep,mod];assert len(ll)==12
 sat_edges=[r for r in prior[ep] if r['source']=='Satb2' and r['target']==g]
 genuine=[r for r in prior[ep] if r['target']==g]
 # Match only actual saved used prior links, exact target expression endpoint.
 for edge in genuine:
  tf=edge['source'];activity=tfa[ep,tf];sensitivity=tfp[ep,tf];rna=sourceRNA[ep][tf]
  eligibleRNA=len(rna)==1 and rna[0]['eligible_activity_background']=='True'
  absentRNA=len(rna)==1 and rna[0]['all_zero_all12']=='True'
  prior_rows.append({'endpoint':ep,'target_gene':g,'role':components[ep,g]['role'],'TF_prior_source':tf,'native_prior_weight':edge['weight'],'prior_resources':edge['resources'],'prior_reference_PMIDs':edge['references'],'prior_self_target':edge['self_target'],'TF_mouse_score_effect_cKO_minus_WT':activity['effect_cKO_minus_WT'],'TF_primary_association_FDR':activity['association_fdr'],'TF_same_partition_effect':sensitivity['effect_cKO_minus_WT'],'TF_same_partition_FDR':sensitivity['association_fdr'],'TF_source_RNA_unique_assayed_row':len(rna)==1,'TF_source_RNA_all_zero_all12':absentRNA,'TF_source_RNA_meets_background_gate':eligibleRNA,'TF_source_RNA_mice_count_ge5':rna[0]['n_mice_raw_count_ge5'] if len(rna)==1 else 'NA','target_gene_DESeq2_log2FC':d['log2FoldChange'],'target_gene_DESeq2_padj':d['padj'],'evidence_boundary':'saved_curated_prior_and_same_dataset_target_signature_not_binding_or_causal_chain','source_file':str(T/f'tf_activity/{ep}/mouse_activity_used_signed_prior_edges.tsv')})
 component_rows.append({'endpoint':ep,'lineage':lin,'gene_symbol':g,'role':components[ep,g]['role'],'ordered_pair_keys':';'.join(sorted(components[ep,g]['pair_keys'])),'DESeq2_log2FC_cKO_over_WT':d['log2FoldChange'],'DESeq2_padj':d['padj'],'DE_direction':'higher_in_cKO' if sign(d['log2FoldChange'])>0 else 'lower_in_cKO','module_network':'pooled','module_id':lin+'/pooled/'+mod,'module':mod,'native_nucleus_kME':m['kME'],'native_within_module_hub_rank':ranks[lin][g]['hub_rank'],'module_primary_effect_cKO_minus_WT':a['effect_cKO_minus_WT'],'module_primary_FDR':a['padj'],'module_same_partition_effect':s['effect_cKO_minus_WT'],'module_same_partition_FDR':s['padj'],'module_primary_n_WT':a['n_WT'],'module_primary_n_cKO':a['n_cKO'],'module_same_partition_n_WT':s['n_WT'],'module_same_partition_n_cKO':s['n_cKO'],'module_LOO_same_sign_n':sum(sign(r['effect_cKO_minus_WT'])==sign(a['effect_cKO_minus_WT']) for r in ll),'module_LOO_tested_n':len(ll),'gene_and_module_effect_sign_agree':sign(d['log2FoldChange'])==sign(a['effect_cKO_minus_WT']),'module_primary_FDR_lt05':float(a['padj'])<.05,'module_same_partition_FDR_lt05':float(s['padj'])<.05,'native_Satb2_used_prior_link_to_gene':len(sat_edges)>0,'native_Satb2_used_prior_self_link':g=='Satb2' and len(sat_edges)>0,'n_saved_used_prior_sources_for_gene':len(genuine),'n_saved_prior_sources_primary_TF_FDR_lt05':sum(tfa[ep,r['source']]['association_fdr']!='NA' and float(tfa[ep,r['source']]['association_fdr'])<.05 for r in genuine),'source_gene_DE_file':cfg['existing_deseq2'][ep]['path'],'source_module_membership_file':str(H/f'hdwgcna/{lin}/pooled/module_membership.tsv'),'source_module_effect_file':str(H/'hdwgcna/pooled_module_genotype_associations.tsv')})
# Check strict LR's component DE values against the original endpoint-specific seven-column files.
for r in pairs:
 for side,ef,gf in [('ligand','source','ligand_complex'),('receptor','target','receptor_complex')]:
  d=de[r[ef],r[gf]];assert float(r[side+'_deseq2_log2fc'])==float(d['log2FoldChange']);assert float(r[side+'_deseq2_padj'])==float(d['padj'])
selected_pairs=[{k:r[k] for k in ['condition','direction','source','target','ligand_complex','receptor_complex','ligand_deg_direction','receptor_deg_direction','ligand_deseq2_log2fc','ligand_deseq2_padj','receptor_deseq2_log2fc','receptor_deseq2_padj','magnitude_rank','specificity_rank','lr_logfc','lr_gmeans','ligand_props','receptor_props','sample_support_WT','sample_support_cKO','support_ge3_both_genotypes']} for r in pairs]
outputs=[]
for name,r in [('components',component_rows),('directional_pairs',selected_pairs),('saved_prior_overlap',prior_rows)]:outputs.append(write(O/('20261028_LR_HD_TF_'+name+'.tsv'),r))
selected_component_keys=set(components)-{('CA1_Deep','Satb2'),('CA1_Superficial','Satb2')}
sig_overlap=[r for r in prior_rows if r['role']!='SATB2_gene_context_not_LR_component' and r['TF_primary_association_FDR']!='NA' and float(r['TF_primary_association_FDR'])<.05]
assert len(sig_overlap)==1 and sig_overlap[0]['TF_prior_source']=='Hoxa5' and sig_overlap[0]['target_gene']=='Ptn' and sig_overlap[0]['endpoint']=='CA1_Superficial' and sig_overlap[0]['TF_source_RNA_all_zero_all12'] is True
assert sum(r['native_Satb2_used_prior_link_to_gene'] for r in component_rows if r['role']!='SATB2_gene_context_not_LR_component')==0
summary={'schema':'satb2-completed-LR-HD-TF-gene-join/v1','date':'2026-10-28','scope':'Five prespecified ordered LR population pairs plus two CA1 Satb2 gene-context rows; all evidence from already completed analyses only','new_scientific_fits':0,'native_files_modified':False,'selection_changes':False,'LIANA_reranking':False,'future_motif_results_included':False,'counts':{'ordered_LR_pairs':5,'condition_specific_LR_rows':15,'unique_LR_component_endpoint_gene_rows':8,'Satb2_gene_context_rows':2,'all_component_rows':10,'all_genuine_used_prior_links_including_Satb2_context':len(prior_rows),'primary_FDR_TF_overlap_with_LR_components':len(sig_overlap),'Satb2_used_prior_links_to_selected_LR_components':0},'input_receipts':list(sources.values()),'outputs':outputs,'primary_FDR_TF_overlap':sig_overlap,'boundaries':['The same measured RNA supplies gene DE, module structure, TF target signatures and LR qualification; convergence is not independent validation.','Native module directions are PC1 score effects, not directions of every member gene and not communication activity effects.','Original LIANA lr_logfc is group-versus-rest, not genotype logFC; rank is not genotype FDR.','No significant curated TF overlap involving these leads has both expressed source RNA and matched-partition FDR support; Hoxa5→Ptn is a prior join with source RNA absent.','No direct Satb2 prior link to these selected LR components is present in the saved used network; absence from a curated resource does not exclude biology.','Only exact pooled-module identity is used; colors across lineages or consensus fits are not interchangeable.'],'primary_literature_context':[{'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6605312/','purpose':'PTN–PTPRZ has demonstrated oligodendrocyte-lineage context and inhibitory-ligand biology; supports an identity/context caveat, not proof of contamination in this cohort.'}]}
for path,old in sources.items():assert receipt(path)==old
(O/'20261028_LR_HD_TF_sources.json').write_text(json.dumps(summary,indent=2,sort_keys=True,allow_nan=False)+'\n')
print(json.dumps({'counts':summary['counts'],'components':component_rows,'significant_overlap':sig_overlap,'outputs':outputs},indent=2))
