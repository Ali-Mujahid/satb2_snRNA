#!/usr/bin/env python3
"""Read-only integrity and semantic validation of the additive v3.4 marker data.

Uses the Python standard library. This checks saved bytes and descriptive
arithmetic only; it never opens native count objects or performs scientific fits.
"""
from pathlib import Path
import gzip
import hashlib
import json
import math
import re

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / 'site'
TAG = 'atlas-v3.4.0'
PREFIX = 'https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/' + TAG + '/'
BRANCHES = {'original_full', 'paired_original', 'paired_cellbender'}
SAFE = 2**53 - 1
ALLOWED_ALIASES = {
    'Nov': ('Ccn3', 'ENSMUSG00000037362', 'https://www.ncbi.nlm.nih.gov/gene/18133'),
    'Fam19a1': ('Tafa1', 'ENSMUSG00000059187', 'https://www.ncbi.nlm.nih.gov/gene/320265'),
    'Fam19a2': ('Tafa2', 'ENSMUSG00000044071', 'https://www.ncbi.nlm.nih.gov/gene/268354'),
    'Neph1': ('Kirrel', 'ENSMUSG00000041734', 'https://www.ncbi.nlm.nih.gov/gene/170643'),
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def file_sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def close(a, b):
    return finite(a) and finite(b) and math.isclose(a, b, rel_tol=2e-12, abs_tol=1e-10)


def load_json(path):
    return json.loads(path.read_bytes())


def catalog(name, assignment, receipt):
    path = SITE / 'data' / (name + '.json')
    script = SITE / 'data' / (name.lower() + '.js')
    for file in [path, script]:
        assert file_sha(file) == receipt['new_catalog_sha256'][str(file.relative_to(ROOT))], str(file)
    value = load_json(path)
    js = script.read_text().strip()
    start = 'window.' + assignment + '='
    assert js.startswith(start) and js.endswith(';'), name
    assert json.loads(js[len(start):-1]) == value, name + ' JSON/JS disagreement'
    assert value['release_tag'] == TAG and value['scientific_model_fitting_performed'] is False
    return value


def read_payload(ident, desc, folder, receipt_hashes):
    name = ident + '.json.gz'
    assert re.fullmatch(r'[-A-Za-z0-9_.]+\.json\.gz', name), name
    assert desc['path'] == PREFIX + folder + '/' + name, ('Unpinned extension path', ident)
    assert desc['encoding'] == 'gzip-bytes'
    path = ROOT / folder / name
    assert path.resolve().parent == (ROOT / folder).resolve() and not path.resolve().is_relative_to(SITE.resolve())
    packed = path.read_bytes()
    assert sha(packed) == desc['sha256'] == receipt_hashes[name], ident
    raw = gzip.decompress(packed)
    assert sha(raw) == desc['payload_sha256'], ident
    assert len(packed) == desc['gzip_bytes'] and len(raw) == desc['uncompressed_bytes'], ident
    return json.loads(raw)


def verify_panels(merged_genes):
    panel_path = SITE / 'data/USER_MARKER_PANELS.json'
    panels = load_json(panel_path)
    resolved = load_json(ROOT / 'analysis/dashboard/v34_user_panels/RESOLVED_PANELS.json')
    assert panels == resolved, 'Published requested panels differ from reviewed resolution'
    assert panels['schema'] == 'satb2-user-marker-panel-resolved-1'
    requested_file = ROOT / 'analysis/dashboard/v34_user_panels/REQUESTED_PANELS.json'
    assert file_sha(requested_file) == panels['request_sha256']
    requested = {p['id']: p for p in load_json(requested_file)['panels']}
    assert len(panels['panels']) == len(requested) == 3
    alias_names = set()
    counts = []
    for panel in panels['panels']:
        words = [w.strip() for w in requested[panel['id']]['requested_text'].split(',') if w.strip()]
        assert words == panel['requested_symbols']
        assert len(words) == len(panel['entries']) == panel['requested_entry_count']
        ids = []
        unresolved = []
        for i, entry in enumerate(panel['entries'], start=1):
            symbol = entry['requested_symbol']
            assert symbol == words[i-1] and entry['position'] == i
            gid, source_symbol, status = entry['gene_id'], entry['source_gene_symbol'], entry['status']
            if status == 'duplicate':
                first = entry['first_requested_position']
                assert type(first) is int and 0 < first < i and words[first-1] == symbol
                original = panel['entries'][first-1]
                assert original['gene_id'] == gid and original['source_gene_symbol'] == source_symbol
                status = entry['resolution_status']
                assert status in ['exact','verified_alias','unresolved']
                if gid is not None:
                    assert entry['duplicate_selected_gene_id'] is True and gid in ids
            if status == 'unresolved':
                assert gid is None and source_symbol is None
                if symbol not in unresolved:
                    unresolved.append(symbol)
                continue
            assert gid in merged_genes and merged_genes[gid] == source_symbol
            if status == 'exact':
                assert symbol == source_symbol
            else:
                assert status == 'verified_alias' and symbol in ALLOWED_ALIASES
                expected_symbol, expected_id, url = ALLOWED_ALIASES[symbol]
                assert (source_symbol, gid, entry['alias_evidence_url']) == (expected_symbol, expected_id, url), entry
                alias_names.add(symbol)
            if gid not in ids:
                ids.append(gid)
        assert ids == panel['gene_ids'] and len(ids) == panel['selected_gene_count']
        assert unresolved == panel['unresolved_symbols']
        assert len(set(ids)) == len(ids)
        counts.append(dict(id=panel['id'], requested=panel['requested_entry_count'],
                           selected=len(ids), unresolved_unique=len(unresolved)))
    assert alias_names == set(ALLOWED_ALIASES), alias_names
    return counts


def main():
    receipt_path = ROOT / 'analysis/dashboard/v34_user_panels/export/COMPLETE.json'
    receipt = load_json(receipt_path)
    assert receipt['status'] == 'COMPLETE' and receipt['release_tag'] == TAG
    for key in ['scientific_model_fitting_performed', 'previous_scientific_results_modified',
                'original_889_gene_marker_shards_modified', 'original_2670_expression_files_modified']:
        assert receipt[key] is False, key
    assert receipt['original_source_axis_reused'] is True
    assert receipt['new_descriptive_native_count_aggregation_performed'] is True
    assert receipt['n_new_genes'] == 166 and receipt['n_datasets'] == 24
    assert receipt['n_new_gene_vector_shards'] == 498 and receipt['n_new_marker_shards'] == 24
    producer = ROOT / 'scripts/extend_user_marker_panels_v34.py'
    assert file_sha(producer) == receipt['producer_sha256']
    exact_new_file = ROOT / 'analysis/dashboard/v34_user_panels/EXPORT_NEW_GENES.json'
    assert file_sha(exact_new_file) == receipt['new_genes_input_sha256']
    sources_file = ROOT / 'analysis/provenance/PROCESSED_OBJECTS_RELEASE_MANIFEST.json'
    assert file_sha(sources_file) == receipt['source_objects_manifest_sha256']
    sources = {r['branch']:r for r in load_json(sources_file)['files'] if r['kind'] == 'shared_counts'}
    for relative, expected in receipt['existing_catalog_sha256'].items():
        assert file_sha(ROOT / relative) == expected, ('Frozen catalog changed', relative)
    for relative, expected in receipt['native_output_sha256'].items():
        path = receipt_path.parent / relative
        assert path.resolve().is_relative_to(receipt_path.parent.resolve())
        assert file_sha(path) == expected, relative
    assert file_sha(ROOT / 'analysis/v2/executed_workflow_snapshot/plots/summarize_aggregates.py') == receipt['native_summary_code_sha256']
    assert file_sha(ROOT / 'analysis/data/annotation_initial/cluster_annotations.tsv') == receipt['initial_annotation_sha256']
    audit_branches = {r['branch']:r for r in receipt['branch_validation']}
    assert set(audit_branches) == BRANCHES
    integer_checks = 0
    for branch, audit in audit_branches.items():
        assert audit['source_sha256'] == sources[branch]['sha256']
        assert audit['n_source_genes'] == 32285 and audit['n_new_genes'] == 166
        for key in ['all_source_libraries_exact', 'original_source_axis_reused_exactly', 'decoded_new_vectors_equal_native_counts']:
            assert audit[key] is True
        assert audit['old_count_and_detection_agreement'] == 'EXACT'
        assert audit['independent_all_new_gene_integer_sum_agreement'] == 'EXACT'
        assert audit['max_abs_native_log_sum_error'] < 1e-10
        integer_checks += audit['n_independent_new_gene_atom_mouse_integer_checks']
    assert integer_checks == 974088

    manifest = load_json(SITE / 'data/manifest.json')
    old_vectors = load_json(SITE / 'data/GENE_EXPRESSION_CATALOG.json')
    marker = catalog('MARKER_EXTENSION_CATALOG', 'SATB2_MARKER_EXTENSION_CATALOG', receipt)
    vectors = catalog('GENE_EXPRESSION_EXTENSION_CATALOG', 'SATB2_GENE_EXPRESSION_EXTENSION_CATALOG', receipt)
    assert marker['schema'] == 'satb2-user-marker-extension-catalog-1'
    assert vectors['schema'] == 'satb2-nucleus-expression-extension-catalog-1'
    assert marker['genes'] == vectors['genes'] == load_json(exact_new_file)
    gene_symbols = {g['gene_id']:g['gene_symbol'] for g in marker['genes']}
    assert len(marker['genes']) == len(gene_symbols) == len(set(gene_symbols.values())) == 166
    old_genes = {g['gene_id']:g['gene_symbol'] for g in manifest['genes']}
    assert len(old_genes) == 889 and not set(old_genes).intersection(gene_symbols)
    merged_genes = old_genes | gene_symbols
    assert len(merged_genes) == 1055
    all_source_genes = {g['gene_id']:g['gene_symbol'] for g in load_json(SITE / 'data/MOUSE_CATALOG.json')['genes']}
    assert all(all_source_genes[gid] == sym for gid, sym in gene_symbols.items())
    old_ids = set(manifest['shards']) | set(old_vectors['shards'])
    for file in ['DE_CATALOG.json', 'MOUSE_CATALOG.json']:
        old_ids.update(load_json(SITE / 'data' / file)['shards'])
    assert not set(marker['shards']).intersection(old_ids)
    assert not set(vectors['shards']).intersection(old_ids | set(marker['shards']))
    assert len(marker['shards']) == len(marker['datasets']) == 24
    assert len(vectors['shards']) == 498 and set(vectors['branches']) == BRANCHES
    assert vectors['normalization'] == old_vectors['normalization'] == receipt['normalization']
    assert vectors['zero_policy'] == old_vectors['zero_policy'] == receipt['zero_policy']
    assert vectors['existing_axis_catalog_sha256'] == file_sha(SITE / 'data/GENE_EXPRESSION_CATALOG.json')
    assert marker['cKO_z_policy'] == receipt['cKO_z_policy']
    panel_counts = verify_panels(merged_genes)

    vector_files = set()
    sparse_values = 0
    vector_bytes = 0
    for branch, item in vectors['branches'].items():
        base = old_vectors['branches'][branch]
        assert item['axis_shard'] == base['axis_shard'] == 'axis__' + branch
        assert item['n_source_nuclei'] == base['n_source_nuclei'] == audit_branches[branch]['n_source_nuclei']
        for key in ['source_count_object_sha256', 'source_count_object_release_asset']:
            assert item[key] == base[key]
        assert item['source_count_object_sha256'] == sources[branch]['sha256']
        assert item['source_count_object_release_asset'] == sources[branch]['asset_name']
        axis_desc = old_vectors['shards'][item['axis_shard']]
        axis_file = ROOT / 'gene_expression' / (item['axis_shard'] + '.json.gz')
        packed = axis_file.read_bytes(); raw = gzip.decompress(packed)
        assert sha(packed) == axis_desc['sha256'] and sha(raw) == axis_desc['payload_sha256']
        axis = json.loads(raw)
        n = item['n_source_nuclei']; libraries = axis['library_counts']
        assert axis['branch'] == branch and len(axis['nucleus_ids']) == len(libraries) == n
        assert len(set(axis['nucleus_ids'])) == n
        assert all(type(count) is int and 0 < count <= SAFE for count in libraries)
        assert set(item['genes']) == set(gene_symbols)
        for gid, ident in item['genes'].items():
            assert ident == 'gene_v34__' + branch + '__' + gid
            desc = vectors['shards'][ident]
            payload = read_payload(ident, desc, 'gene_expression_v34', receipt['new_gene_vector_sha256'])
            assert payload['schema'] == 'satb2-nucleus-gene-counts-1'
            assert payload['branch'] == branch and payload['gene_id'] == gid and payload['gene_symbol'] == gene_symbols[gid]
            assert payload['n_source_nuclei'] == n and payload['axis_payload_sha256'] == axis_desc['payload_sha256']
            assert len(payload['indices']) == len(payload['counts'])
            last = -1
            for i, count in zip(payload['indices'], payload['counts']):
                assert type(i) is int and last < i < n
                assert type(count) is int and 0 < count <= libraries[i] <= SAFE
                last = i
            sparse_values += len(payload['counts'])
            vector_files.add(ident + '.json.gz'); vector_bytes += desc['gzip_bytes']
    assert vector_files == {p.name for p in (ROOT / 'gene_expression_v34').glob('*.json.gz')}
    assert vector_files == set(receipt['new_gene_vector_sha256'])

    source_datasets = {d['id']:d for d in manifest['datasets']}
    marker_files = set()
    marker_bytes = 0
    summary_rows = 0
    references = {}
    scaled_rows = []
    additive_fields = ['n_nuclei', 'sum_counts', 'n_detected', 'sum_log1p_normalized',
                       'n_mice_with_nuclei', 'n_mice_eligible_for_expression_review', 'n_mice_expression_supported']
    assert {d['id'] for d in marker['datasets']} == set(source_datasets)
    for ds in marker['datasets']:
        base = source_datasets[ds['id']]
        for key in ['tool', 'branch', 'population', 'cluster_scope']:
            assert ds[key] == base[key]
        assert ds['base_markers_shard'] == base['markers_shard']
        ident = ds['markers_shard']
        assert ident == 'markers_v34__' + ds['id']
        desc = marker['shards'][ident]
        payload = read_payload(ident, desc, 'marker_extensions_v34', receipt['new_marker_shard_sha256'])
        assert payload['schema'] == 'satb2-user-marker-extension-1'
        assert payload['dataset_id'] == ds['id']
        assert payload['gene_ids'] == list(gene_symbols)
        assert payload['base_marker_shard'] == base['markers_shard']
        old_desc = manifest['shards'][base['markers_shard']]
        assert payload['base_marker_payload_sha256'] == old_desc['payload_sha256']
        old_packed = (SITE / old_desc['path']).read_bytes(); old_raw = gzip.decompress(old_packed)
        assert sha(old_packed) == old_desc['sha256'] and sha(old_raw) == old_desc['payload_sha256']
        old_payload = json.loads(old_raw)
        scope_maps = {}
        for scope in ['rows', 'WT_rows', 'cKO_rows']:
            rows = payload[scope]
            old_denominators = {r['cluster_id']:r['n_nuclei'] for r in old_payload[scope]}
            pairs = {(r['cluster_id'],r['gene_id']):r for r in rows}
            assert len(pairs) == len(rows) == len(old_denominators) * 166
            assert set(pairs) == {(cid,gid) for cid in old_denominators for gid in gene_symbols}
            scope_maps[scope] = pairs
            for (cid,gid), row in pairs.items():
                assert row['gene_symbol'] == gene_symbols[gid]
                assert row['n_nuclei'] == old_denominators[cid]
                for field in ['n_nuclei','sum_counts','n_detected']:
                    assert finite(row[field]) and 0 <= row[field] <= SAFE and row[field] == math.floor(row[field])
                assert row['n_detected'] <= row['n_nuclei']
                assert finite(row['sum_log1p_normalized']) and row['sum_log1p_normalized'] >= 0
                n = row['n_nuclei']
                for numerator, mean in [('sum_counts','mean_source_counts'),('n_detected','fraction_detected'),
                                        ('sum_log1p_normalized','mean_log1p_normalized')]:
                    if n:
                        assert close(row[mean],row[numerator]/n), (ds['id'],scope,cid,gid,mean)
                    else:
                        assert row[numerator] == 0 and row[mean] is None
                if scope == 'cKO_rows':
                    assert row['scaled_expression_z'] is None
                    assert row['scale_reference_branch'] is None and row['scale_reference_mean'] is None and row['scale_reference_sd'] is None
                else:
                    key = (ds['tool'],ds['population'],scope,gid,ds['branch'])
                    if row['keep_for_plot'] and finite(row['mean_log1p_normalized']):
                        references.setdefault(key,[]).append(row['mean_log1p_normalized'])
                    scaled_rows.append((ds['tool'],ds['population'],scope,gid,ds['branch'],row['scale_reference_branch'],
                                        row['scale_reference_mean'],row['scale_reference_sd'],row['mean_log1p_normalized'],row['scaled_expression_z']))
            summary_rows += len(rows)
        for key, row in scope_maps['rows'].items():
            wt, cko = scope_maps['WT_rows'][key], scope_maps['cKO_rows'][key]
            for field in additive_fields:
                assert close(row[field],wt[field]+cko[field]), (ds['id'],key,field)
        marker_files.add(ident + '.json.gz'); marker_bytes += desc['gzip_bytes']
    assert marker_files == {p.name for p in (ROOT / 'marker_extensions_v34').glob('*.json.gz')}
    assert marker_files == set(receipt['new_marker_shard_sha256'])
    scales = {}
    for key, values in references.items():
        mean = math.fsum(values)/len(values)
        sd = math.sqrt(math.fsum((v-mean)**2 for v in values)/len(values))
        scales[key] = (mean,sd)
    for tool,pop,scope,gid,branch,refbranch,mean,sd,expr,z in scaled_rows:
        assert refbranch == ('paired_original' if branch == 'paired_cellbender' else branch)
        expected_mean,expected_sd = scales[(tool,pop,scope,gid,refbranch)]
        assert close(mean,expected_mean) and close(sd,expected_sd)
        if finite(expr):
            assert close(z, (expr-mean)/sd if sd > 0 else 0)
        else:
            assert z is None
    print(json.dumps(dict(status='PASS', mode='READ_ONLY_SAVED_BYTES_AND_DESCRIPTIVE_ARITHMETIC',
        new_genes=166,total_marker_genes=1055,marker_shards=24,gene_vector_shards=498,
        marker_summary_rows=summary_rows,scaled_rows_checked=len(scaled_rows),sparse_nucleus_counts_checked=sparse_values,
        independent_native_integer_checks_recorded=integer_checks,
        marker_bytes_outside_pages=marker_bytes,gene_vector_bytes_outside_pages=vector_bytes,
        original_catalog_hashes_unchanged=True,old_axes_reused=True,
        exact_pinned_tag=TAG,panels=panel_counts,scientific_model_fitting_performed=False)))


if __name__ == '__main__':
    main()
