from pathlib import Path
import csv, collections, json

OUT=Path('/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/audit_cellbender')
META=Path('/Users/mujahid357/Desktop/SATB2_LiSC_archive_transfer_20261012.BmzvXa/office_analysis/runs/01_original_reference_ingest_v1_0_1/results/original_inputs/SELECTED_NUCLEUS_ACTIVE_METADATA.tsv')
metadata={r['nucleus_id']:r for r in csv.DictReader(META.open(),delimiter='\t')}
aggregate=collections.defaultdict(lambda:{'primary_nuclei':0,'paired_supported':0,'zero_corrected':0,'not_analyzed':0})
membership=[]
def yes(v): return v=='TRUE'
for path in sorted(OUT.glob('*_original_cellbender_primary_paired_roster.tsv')):
    for r in csv.DictReader(path.open(),delimiter='\t'):
        m=metadata[r['nucleus_id']]
        if m['sample']!=r['sample'] or m['barcode']!=r['barcode']:raise ValueError('source identity mismatch')
        populations=['ALL']
        if yes(m['REF_IN_HS']):populations.append('ReferenceHighSupport_IN')
        if yes(m['REF_EX_HS']):populations.append('ReferenceHighSupport_EX')
        if yes(m['REF_CA1_PROS_HS']):populations.append('ReferenceHighSupport_CA1_ProS')
        for level,label in [('broad',m['reference_broad_prediction']),('subclass',m['reference_subclass_name'])]:
            populations.append('Reference_'+level+':'+label)
        for pop in populations:
            key=(r['sample'],pop)
            a=aggregate[key];a['primary_nuclei']+=1;a['paired_supported']+=yes(r['paired_sensitivity_usable']);a['zero_corrected']+=not yes(r['cellbender_nonzero_full']);a['not_analyzed']+=not yes(r['cellbender_analyzed'])
        if yes(m['REF_IN_HS']) or yes(m['REF_CA1_PROS_HS']):
            membership.append({'sample':r['sample'],'barcode':r['barcode'],'nucleus_id':r['nucleus_id'],'ReferenceHighSupport_IN':m['REF_IN_HS'],'ReferenceHighSupport_CA1_ProS':m['REF_CA1_PROS_HS'],'paired_sensitivity_usable':r['paired_sensitivity_usable']})
records=[]
for (sample,pop),a in sorted(aggregate.items()):records.append({'sample':sample,'frozen_reference_population':pop,**a,'zero_fraction':a['zero_corrected']/a['primary_nuclei']})
with (OUT/'frozen_reference_coverage.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
with (OUT/'identity_marker_diagnostic_membership.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(membership[0]),delimiter='\t');w.writeheader();w.writerows(membership)
for pop in ['ALL','ReferenceHighSupport_IN','ReferenceHighSupport_EX','ReferenceHighSupport_CA1_ProS']:
    p=[r for r in records if r['frozen_reference_population']==pop]
    print(pop,{k:sum(r[k] for r in p) for k in ['primary_nuclei','paired_supported','zero_corrected','not_analyzed']})
