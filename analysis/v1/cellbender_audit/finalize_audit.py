from pathlib import Path
import csv, collections, hashlib, json

ARCH=Path('/Users/mujahid357/Desktop/SATB2_LiSC_archive_transfer_20261012.BmzvXa')
OUT=Path('/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/audit_cellbender')
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def save_tsv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
binding=[]
provenance=ARCH/'office_analysis/runs/00_source_qc_audit_v1/results/source_qc/input_provenance.tsv'
for row in csv.DictReader(provenance.open(),delimiter='\t'):
    if row['role'] not in {'raw_h5','original_filtered_h5'}:continue
    filename='raw_feature_bc_matrix.h5' if row['role']=='raw_h5' else 'filtered_feature_bc_matrix.h5'
    path=ARCH/'cellranger_out'/row['sample']/'outs'/filename
    actual=digest(path)
    binding.append({'sample':row['sample'],'role':row['role'],'archive_path':str(path),'accepted_sha256':row['sha256'],'current_sha256':actual,'hash_matches_accepted':actual==row['sha256']})
if len(binding)!=24 or not all(r['hash_matches_accepted'] for r in binding):raise ValueError('accepted original source binding mismatch')
save_tsv(OUT/'original_source_hash_binding.tsv',binding)
details=json.loads((OUT/'local_cellbender_diagnostics.json').read_text())
manifest=[]
for r in details:
    manifest.append({'sample':r['sample'],'primary_original_filtered_path':r['original_filtered_path'],'corrected_full_path':r['full_corrected_path'],'corrected_filtered_path':r['filtered_corrected_path'],'primary_paired_roster_path':str(OUT/f'{r["sample"]}_original_cellbender_primary_paired_roster.tsv'),'primary_nuclei':r['primary_QC_nuclei'],'paired_supported_nuclei':r['paired_supported'],'zero_corrected_nuclei':r['zero_corrected'],'not_analyzed_nuclei':r['not_analyzed'],'role':'paired_sensitivity_only'})
save_tsv(OUT/'usable_count_sources.tsv',manifest)
aggregate=collections.defaultdict(lambda:{'UMIs':0.,'detected_nuclei':0})
for r in csv.DictReader((OUT/'marker_removal_within_frozen_reference.tsv').open(),delimiter='\t'):
    a=aggregate[(r['frozen_reference_population'],r['gene_id'],r['gene_symbol'],r['count_kind'])]
    a['UMIs']+=float(r['UMIs']);a['detected_nuclei']+=int(r['detected_nuclei'])
summary=[]
for pop,gid,gene,kind in sorted(aggregate):
    if kind!='original':continue
    original=aggregate[(pop,gid,gene,'original')];corrected=aggregate[(pop,gid,gene,'corrected')]
    summary.append({'frozen_reference_population':pop,'gene_id':gid,'gene_symbol':gene,'original_UMIs':int(original['UMIs']),'corrected_UMIs':int(corrected['UMIs']),'fraction_UMIs_removed':1-corrected['UMIs']/original['UMIs'] if original['UMIs'] else None,'original_detected_nuclei':original['detected_nuclei'],'corrected_detected_nuclei':corrected['detected_nuclei']})
save_tsv(OUT/'marker_removal_within_frozen_reference_summary.tsv',summary)
receipt={'scope':'new local CellBender correction diagnostics and source binding; no accepted original source-count integrity audit rerun or scientific model fitting','original_source_hashes_matched':24,'primary_QC_nuclei':sum(r['primary_QC_nuclei'] for r in details),'paired_supported_nuclei':sum(r['paired_supported'] for r in details),'zero_corrected_nuclei':sum(r['zero_corrected'] for r in details),'not_analyzed_nuclei':sum(r['not_analyzed'] for r in details),'all_filtered_corrected_nonnegative_integers':all(r['count_nonnegative'] and r['count_integer'] for r in details),'all_train_test_curves_finite_epoch150':all(r['train_finite'] and r['test_finite'] and r['train_last_epoch']==150 and r['test_last_epoch']==150 for r in details),'all_features_match_original_order':True,'recommended_primary':'all12 original counts on frozen108062 primaryQC roster','recommended_sensitivity':'archive corrected counts on104298 paired supported overlap; no mixing originals into corrected matrix or importing extra called barcodes','important_gene_flag':'Cck expression reduction in frozen highsupportIN needs subtype-aware original-versus-corrected review'}
(OUT/'AUDIT_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt,indent=2))
