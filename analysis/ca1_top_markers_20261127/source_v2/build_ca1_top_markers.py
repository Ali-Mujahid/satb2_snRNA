"""Select two Seurat CA1 candidates per frozen cluster from saved marker audits.

No marker test, differential model, embedding or annotation is rerun. C7 and
its exploratory genes are excluded at the user request. Columns are ordered
by pooled peak expression, from the first displayed cluster to the last.
"""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import os

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(os.environ.get('SATB2_MARKER_SUMMARIES', str(
    ROOT.parent / 'SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0/results/marker_summaries')))
OUT = ROOT / 'site/ca1_sublayers/top_markers'
AUDIT = ROOT / 'analysis/ca1_top_markers_20261127'
OUT.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)


def read(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return list(csv.DictReader(f, delimiter='\t', lineterminator='\n'))


def scope(rows):
    return [r for r in rows if r['tool'] == 'seurat' and r['branch'] == 'original_full'
            and r['population'] == 'CA1_ProS' and r['cluster_scope'] == 'original_frozen']


rep_path = SOURCE / 'representative_markers.tsv'
audit_path = SOURCE / 'representative_candidate_audit.tsv.gz'
registry_path = SOURCE / 'representative_cluster_registry.tsv'
shard_path = ROOT / 'site/data/shards/markers__seurat__original_full__CA1_ProS.json.gz'
representatives = scope(read(rep_path))
candidates = scope(read(audit_path))
registry = scope(read(registry_path))
native = json.loads(gzip.decompress(shard_path.read_bytes()))
clusters = sorted([r for r in registry if r['cluster_id'] != '7'], key=lambda r: int(r['cluster_id']))
assert [int(r['cluster_id']) for r in clusters] == [i for i in range(1, 12) if i != 7]
selected = []
for c in clusters:
    cid = c['cluster_id']
    ranked = sorted([r for r in representatives if r['cluster_id'] == cid],
                    key=lambda r: int(r['representative_rank']))
    supported = bool(ranked)
    if not supported:
        ranked = sorted([r for r in candidates if r['cluster_id'] == cid],
                        key=lambda r: (-float(r['pct_difference_vs_rest']),
                                       -float(r['mean_log1p_difference_vs_rest']), r['gene_id']))
    for rank, r in enumerate(ranked[:2], 1):
        selected.append(dict(source_cluster=cid, source_rank=rank, gene_id=r['gene_id'],
            gene_symbol=r['gene_symbol'], supported=supported,
            selection_status='saved_supported_representative' if supported else 'exploratory_no_supported_marker',
            fraction_detected=float(r['fraction_detected']),
            detection_difference_vs_rest=float(r['pct_difference_vs_rest']),
            mean_log1p_difference_vs_rest=float(r['mean_log1p_difference_vs_rest']),
            mice_expression_supported=int(r['n_mice_expression_supported']),
            mice_eligible=int(r['n_mice_eligible_for_expression_review'])))

assert len(selected) == 20 and len({s['gene_id'] for s in selected}) == 20
selected_ids = {s['gene_id'] for s in selected}
manifest = json.loads((ROOT / 'site/data/manifest.json').read_text())
assert selected_ids <= {g['gene_id'] for g in manifest['genes']}
fields = ['cluster_id', 'gene_id', 'gene_symbol', 'n_nuclei', 'sum_counts', 'n_detected',
          'mean_log1p_normalized', 'fraction_detected', 'mean_source_counts', 'scaled_expression_z',
          'scale_definition', 'scale_reference_branch', 'scale_reference_mean', 'scale_reference_sd']
rows = {}
for condition, key in [('ALL', 'rows'), ('WT', 'WT_rows'), ('cKO', 'cKO_rows')]:
    rows[condition] = [{k: r.get(k) for k in fields} for r in native[key] if r['gene_id'] in selected_ids and r['cluster_id'] != '7']
    assert len(rows[condition]) == 10 * 20
    assert len({(r['cluster_id'], r['gene_id']) for r in rows[condition]}) == 10 * 20

# A common display scale uses the existing pooled reference, without altering any saved field.
reference = {r['gene_id']: r for r in rows['ALL']}
for values in rows.values():
    for row in values:
        ref = reference[row['gene_id']]
        row['display_scaled_expression_z'] = ((row['mean_log1p_normalized'] - ref['scale_reference_mean']) / ref['scale_reference_sd']) if ref['scale_reference_sd'] else None
        row['display_scale_definition'] = 'derived_from_saved_pooled_original_11_cluster_reference; no_refit_or_row_rescaling'
fields += ['display_scaled_expression_z', 'display_scale_definition']

# Order columns by the pooled peak-expression cluster, from the top row down.
for marker in selected:
    values = [r for r in rows['ALL'] if r['gene_id'] == marker['gene_id']]
    peak = sorted(values, key=lambda r: (-r['mean_log1p_normalized'], int(r['cluster_id'])))[0]
    marker['display_peak_cluster'] = peak['cluster_id']
    marker['display_peak_mean_log1p'] = peak['mean_log1p_normalized']

cluster_data = []
for c in clusters:
    cid = c['cluster_id']
    program = 'Superficial-associated' if cid in ['2', '3'] else 'Deep-associated' if cid in ['1', '5'] else 'CA1-prosubiculum' if cid == '8' else 'Identity-conflicting' if cid == '11' else 'Mixed / weak'
    cluster_data.append(dict(id=cid, n=int(c['n_nuclei']), program=program,
                             supported_representatives=int(c['n_supported_representative_genes'])))
provenance = dict(tool='Seurat', branch='original_full', cluster_scope='original_frozen',
    population='CA1_ProS', n_nuclei_reference_population=27703, n_nuclei_displayed=26733, excluded_clusters=['7'], n_mice=12, genotypes='6 WT + 6 cKO',
    selection='Ranks 1 and 2 from the saved supported representative table for C1-C6 and C8-C11. C7 and its exploratory genes are excluded at the user request.',
    column_order='Pooled peak mean-expression cluster, ascending by row ID; original selection order breaks ties. The order is fixed across genotype views and row hiding.',
    universe='Saved curated marker panel, not genome-wide marker discovery',
    condition_selection='Gene selection fixed from pooled full-cohort summaries; WT/cKO controls change display only.',
    dot_area='fraction of nuclei with nonzero source RNA counts',
    color='display_scaled_expression_z = (saved condition mean - saved pooled reference mean) / saved pooled reference SD; original 11-cluster equal-weight reference; no refit or rescaling when rows are hidden. Original scaled_expression_z is retained separately, including unavailable cKO values.',
    expression='mean log1p(10000 × source counts / complete nucleus source library)',
    no_models_refitted=True,
    inputs=[dict(file=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for p in [rep_path, audit_path, registry_path, shard_path]])
data = dict(schema='satb2-ca1-top2-seurat-v2', clusters=cluster_data, markers=selected,
            rows=rows, provenance=provenance)
(OUT / 'data.js').write_text('window.SATB2_CA1_TOP_MARKERS=' + json.dumps(data, separators=(',', ':'), allow_nan=False) + ';\n')
(OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
with (OUT / 'marker_selection.tsv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(selected[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(selected)
with (OUT / 'all_expression_values.tsv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=['condition'] + fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for condition, values in rows.items():
        writer.writerows(dict(condition=condition, **r) for r in values)
(AUDIT / 'selection_and_source_audit.json').write_text(json.dumps(dict(
    provenance=provenance, selection=selected, clusters=cluster_data,
    expected_matrix_cells_per_condition=200, conditions=list(rows),
    excluded_cluster='7', source_values_copied_without_changes=True), indent=2) + '\n')
print(json.dumps(dict(markers=len(selected), supported=sum(s['supported'] for s in selected),
    exploratory=sum(not s['supported'] for s in selected), clusters=len(clusters),
    source_cells_per_condition={k: len(v) for k, v in rows.items()}), indent=2))
