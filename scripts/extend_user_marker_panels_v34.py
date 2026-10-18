#!/usr/bin/env python3
"""Append exact native marker summaries and nucleus counts without changing old data.

This is descriptive export from frozen full-gene integer counts. It never fits a
model or changes QC, nuclei, clusters, embeddings, labels or differential results.
The existing 889-gene marker data and 2,670 nucleus-vector files remain untouched.
New files use a separate immutable atlas-v3.4.0 transport namespace.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

BRANCHES = ['original_full', 'paired_original', 'paired_cellbender']
AUDIT_GENES = ['Ccnd2', 'Pvalb', 'Sst']
TAG = 'atlas-v3.4.0'
BASE_URL = 'https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/' + TAG + '/'


def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def compact(value):
    return json.dumps(value, separators=(',', ':'), ensure_ascii=False, allow_nan=False).encode()


def emit(folder, ident, payload):
    raw = compact(payload)
    packed = gzip.compress(raw, compresslevel=6, mtime=0)
    target = folder / (ident + '.json.gz')
    if target.exists():
        raise FileExistsError(target)
    target.write_bytes(packed)
    return dict(path=BASE_URL + folder.name + '/' + target.name,
                encoding='gzip-bytes', sha256=hashlib.sha256(packed).hexdigest(),
                payload_sha256=hashlib.sha256(raw).hexdigest(),
                gzip_bytes=len(packed), uncompressed_bytes=len(raw))


def write_catalog(site, basename, assignment, catalog):
    target = site / 'data' / (basename + '.json')
    script = site / 'data' / (basename.lower() + '.js')
    for file in [target, script]:
        if file.exists():
            raise FileExistsError(file)
    target.write_bytes(compact(catalog) + b'\n')
    script.write_text('window.' + assignment + '=' + compact(catalog).decode() + ';\n')
    return [target, script]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', type=Path, required=True)
    ap.add_argument('--native-results', type=Path, required=True)
    ap.add_argument('--source-count-root', type=Path, required=True)
    ap.add_argument('--new-genes', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    repo, out = args.repo.resolve(), args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'native').mkdir()
    marker_dir, vector_dir = repo / 'marker_extensions_v34', repo / 'gene_expression_v34'
    for folder in [marker_dir, vector_dir]:
        if folder.exists():
            raise FileExistsError(folder)
        folder.mkdir()

    # Reuse the native descriptive formulas and the established row/missingness contract.
    sys.path.insert(0, str(repo / 'analysis/v2/executed_workflow_snapshot/plots'))
    from summarize_aggregates import summarize, add_mask_and_contrasts, add_reference_scaling
    from plot_common import KEYS, read_table
    sys.path.insert(0, str(repo / 'scripts'))
    from extend_screenshot_markers import marker_row, cko_row

    manifest_file = repo / 'site/data/manifest.json'
    manifest = json.loads(manifest_file.read_text())
    old_gene_file = repo / 'site/data/GENE_EXPRESSION_CATALOG.json'
    old_gene_catalog = json.loads(old_gene_file.read_text())
    source_manifest_file = repo / 'analysis/provenance/PROCESSED_OBJECTS_RELEASE_MANIFEST.json'
    sources = {x['branch']: x for x in json.loads(source_manifest_file.read_text())['files']
               if x['kind'] == 'shared_counts'}
    genes = json.loads(args.new_genes.read_text())
    ids, symbols = [g['gene_id'] for g in genes], [g['gene_symbol'] for g in genes]
    ids_set = set(ids)
    if not genes or len(set(ids)) != len(ids) or len(set(symbols)) != len(symbols):
        raise ValueError('New exact gene list must be nonempty with unique IDs and symbols')
    if set(ids) & {g['gene_id'] for g in manifest['genes']}:
        raise ValueError('New exact gene list overlaps frozen 889-gene marker registry')
    annotations_file = repo / 'analysis/data/annotation_initial/cluster_annotations.tsv'
    annotations = read_table(annotations_file)
    catalogs = dict(schema='satb2-user-marker-extension-catalog-1', release_tag=TAG,
        genes=genes, datasets=[], shards={},
        scientific_model_fitting_performed=False,
        source_scope='Existing full-gene integer counts and frozen original atomic memberships',
        cKO_z_policy='Unavailable as in original atlas; use expression mean for cKO')
    vectors = dict(schema='satb2-nucleus-expression-extension-catalog-1', release_tag=TAG,
        genes=genes, branches={}, shards={}, normalization=old_gene_catalog['normalization'],
        zero_policy=old_gene_catalog['zero_policy'], scientific_model_fitting_performed=False,
        existing_axis_catalog_sha256=sha_file(old_gene_file))
    pooled, wt, branch_checks = [], [], []
    native_mouse_path = out / 'native/new_gene_cluster_mouse_summary.tsv.gz'
    source_invariants = {str(p.relative_to(repo)): sha_file(p) for p in [manifest_file, old_gene_file]}
    header = True
    with gzip.open(native_mouse_path, 'wt') as mouse_handle:
        for branch in BRANCHES:
            prepared = args.native_results / 'prepared' / branch
            receipt_file = prepared / 'COMPLETE.json'
            receipt = json.loads(receipt_file.read_text())
            for name in ['POPULATION_INDEX.tsv', 'POPULATION_MEMBERSHIPS.tsv.gz', 'cluster_mouse_counts.tsv']:
                if sha_file(prepared / name) != receipt['output_sha256'][name]:
                    raise ValueError('Native prepared input hash mismatch: ' + branch + '/' + name)
            source = sources[branch]
            source_path = args.source_count_root / branch / 'branch_counts.h5ad'
            if sha_file(source_path) != source['sha256']:
                raise ValueError('Full integer-count source hash mismatch: ' + branch)
            data = ad.read_h5ad(source_path, backed='r')
            try:
                if data.uns['count_branch'] != branch or not data.obs_names.is_unique or not data.var_names.is_unique:
                    raise ValueError('Native branch or axes mismatch')
                new_columns = data.var_names.get_indexer(ids)
                if (new_columns < 0).any() or any(str(data.var.iloc[j].gene_symbol) != g['gene_symbol']
                                                 for j, g in zip(new_columns, genes)):
                    raise ValueError('New exact gene ID/symbol mapping is absent or inconsistent')
                audit_ids = []
                for symbol in AUDIT_GENES:
                    hits = data.var_names[data.var.gene_symbol.astype(str).eq(symbol)].tolist()
                    if len(hits) != 1:
                        raise ValueError('Nonunique native audit symbol')
                    audit_ids += hits
                all_ids, all_symbols = ids + audit_ids, symbols + AUDIT_GENES
                columns = data.var_names.get_indexer(all_ids)
                n = data.n_obs
                selected = np.empty((n, len(columns)), dtype=np.float64)
                libraries = np.empty(n, dtype=np.int64)
                for start in range(0, n, 8192):
                    stop = min(n, start + 8192)
                    full = sparse.csr_matrix(data.X[start:stop, :])
                    if not np.isfinite(full.data).all() or (full.data < 0).any() or not np.equal(full.data, np.floor(full.data)).all():
                        raise ValueError('Invalid native integer counts')
                    sums = np.asarray(full.sum(axis=1)).ravel()
                    if (sums <= 0).any() or (sums > 2**53-1).any():
                        raise ValueError('Unsafe complete source library')
                    libraries[start:stop] = sums.astype(np.int64)
                    selected[start:stop] = full[:, columns].toarray()
                if not np.array_equal(libraries, data.obs.branch_UMIs_verified.to_numpy(dtype=np.int64)):
                    raise ValueError('Full native library disagrees with frozen verified library')
                old_branch = old_gene_catalog['branches'][branch]
                old_axis_desc = old_gene_catalog['shards'][old_branch['axis_shard']]
                axis_file = repo / 'gene_expression' / (old_branch['axis_shard'] + '.json.gz')
                axis_packed = axis_file.read_bytes()
                axis_raw = gzip.decompress(axis_packed)
                if hashlib.sha256(axis_packed).hexdigest() != old_axis_desc['sha256'] or hashlib.sha256(axis_raw).hexdigest() != old_axis_desc['payload_sha256']:
                    raise ValueError('Frozen source axis transport hash mismatch')
                axis = json.loads(axis_raw)
                if axis['nucleus_ids'] != data.obs_names.tolist() or not np.array_equal(libraries, axis['library_counts']):
                    raise ValueError('New vectors cannot reuse the unchanged source axis')
                gene_shards = {}
                for j, gene in enumerate(genes):
                    ix = np.flatnonzero(selected[:, j] > 0)
                    vals = selected[ix, j].astype(np.int64)
                    ident = 'gene_v34__' + branch + '__' + gene['gene_id']
                    payload = dict(schema='satb2-nucleus-gene-counts-1', branch=branch,
                        gene_id=gene['gene_id'], gene_symbol=gene['gene_symbol'], n_source_nuclei=n,
                        axis_payload_sha256=old_axis_desc['payload_sha256'],
                        indices=ix.tolist(), counts=vals.tolist())
                    desc = emit(vector_dir, ident, payload)
                    decoded = json.loads(gzip.decompress((vector_dir / (ident + '.json.gz')).read_bytes()))
                    if decoded['indices'] != ix.tolist() or decoded['counts'] != vals.tolist():
                        raise ValueError('Decoded new count vector differs from native integer values')
                    vectors['shards'][ident] = desc
                    gene_shards[gene['gene_id']] = ident
                vectors['branches'][branch] = dict(axis_shard=old_branch['axis_shard'], genes=gene_shards,
                    n_source_nuclei=n, source_count_object_sha256=source['sha256'],
                    source_count_object_release_asset=source['asset_name'])

                index, memberships, cluster_counts = (read_table(prepared / name) for name in
                    ['POPULATION_INDEX.tsv', 'POPULATION_MEMBERSHIPS.tsv.gz', 'cluster_mouse_counts.tsv'])
                if memberships.duplicated(['population_key', 'nucleus_id']).any():
                    raise ValueError('Duplicate frozen atomic memberships')
                oi = data.obs_names.get_indexer(memberships.nucleus_id)
                pi = pd.Index(index.population_key).get_indexer(memberships.population_key)
                if (oi < 0).any() or (pi < 0).any():
                    raise ValueError('Frozen atom membership absent from native counts')
                for atom in index.itertuples(index=False):
                    member_ids = sorted(memberships.loc[memberships.population_key.eq(atom.population_key), 'nucleus_id'])
                    digest = hashlib.sha256(('\n'.join(member_ids) + '\n' if member_ids else '').encode()).hexdigest()
                    if len(member_ids) != atom.n_nuclei or digest != atom.membership_sha256:
                        raise ValueError('Frozen atom membership count or digest changed')
                assignment = sparse.csr_matrix((np.ones(len(oi)), (oi, pi)), shape=(n, len(index)))
                normalized = np.log1p(selected * (10000.0 / libraries[:, None]))
                obs = data.obs[['sample', 'condition']].copy()
                records = []
                for sample in sorted(obs['sample'].astype(str).unique()):
                    mask = obs['sample'].astype(str).eq(sample).to_numpy()
                    a = assignment[mask]
                    nn = np.asarray(a.sum(axis=0)).ravel().astype(int)
                    sums = np.asarray(a.T @ selected[mask])
                    det = np.asarray(a.T @ (selected[mask] > 0).astype(float))
                    logs = np.asarray(a.T @ normalized[mask])
                    conditions = obs.loc[mask, 'condition'].astype(str).unique().tolist()
                    if len(conditions) != 1:
                        raise ValueError('Mouse has conflicting genotypes')
                    for k, atom in enumerate(index.itertuples(index=False)):
                        for j, (gid, sym) in enumerate(zip(all_ids, all_symbols)):
                            records.append(dict(tool=atom.tool, branch=branch, population=atom.population,
                                cluster_scope='original_frozen', cluster_id=atom.cluster_id,
                                sample=sample, condition=conditions[0], gene_id=gid, gene_symbol=sym,
                                n_nuclei=int(nn[k]), sum_counts=float(sums[k,j]),
                                n_detected=float(det[k,j]), sum_log1p_normalized=float(logs[k,j])))
                frame = pd.DataFrame(records)
                old_audit = []
                for chunk in pd.read_csv(prepared / 'gene_cluster_mouse_aggregates.tsv.gz', sep='\t',
                                         dtype={'cluster_id': str}, chunksize=150000):
                    old_audit.append(chunk.loc[chunk.gene_symbol.isin(AUDIT_GENES)])
                old_audit = pd.concat(old_audit, ignore_index=True)
                check = frame.loc[frame.gene_symbol.isin(AUDIT_GENES)].merge(old_audit,
                    on=KEYS + ['sample', 'gene_id'], suffixes=('_fresh', '_native'), validate='one_to_one')
                if len(check) != len(index) * 12 * len(AUDIT_GENES):
                    raise ValueError('Native audit lacks complete atom/mouse coverage')
                for field in ['n_nuclei', 'sum_counts', 'n_detected']:
                    if not np.array_equal(check[field + '_fresh'], check[field + '_native']):
                        raise ValueError('Native aggregate audit differs: ' + field)
                log_error = np.abs(check.sum_log1p_normalized_fresh - check.sum_log1p_normalized_native)
                if not np.allclose(check.sum_log1p_normalized_fresh, check.sum_log1p_normalized_native, rtol=2e-12, atol=1e-10):
                    raise ValueError('Native full-library log-sum audit failed')
                fresh = frame.loc[frame.gene_id.isin(ids)].reset_index(drop=True)
                expected = {(r.tool, r.population, r.cluster_id, r.gene_id, r.sample): int(r.sum_counts)
                            for r in fresh.itertuples(index=False)}
                checked = set()
                samples = fresh['sample'].unique().tolist()
                # Independent check of all new gene x atom x mouse sums against saved full-gene native arrays.
                for atom in index.itertuples(index=False):
                    counts_file = Path(atom.counts_file)
                    if sha_file(counts_file) != atom.counts_sha256:
                        raise ValueError('Frozen full-gene native atom count array changed')
                    hits = set()
                    with gzip.open(counts_file, 'rt') as handle:
                        for row in csv.DictReader(handle, delimiter='\t'):
                            if row['gene_id'] not in ids_set:
                                continue
                            hits.add(row['gene_id'])
                            for sample in samples:
                                key = (atom.tool, atom.population, atom.cluster_id, row['gene_id'], sample)
                                if key in checked or int(row[sample]) != expected[key]:
                                    raise ValueError('New marker atom count differs from native full-gene array')
                                checked.add(key)
                    if hits != ids_set:
                        raise ValueError('Native atom array lacks exact new gene IDs')
                if checked != set(expected):
                    raise ValueError('Independent new-gene count audit incomplete')
                mouse, summaries = summarize(fresh, cluster_counts, 5, .1)
                mouse.to_csv(mouse_handle, sep='\t', index=False, na_rep='NA', header=header)
                header = False
                pooled.append(summaries[0]); wt.append(summaries[1])
                branch_checks.append(dict(branch=branch, source_sha256=source['sha256'],
                    n_source_nuclei=n, n_source_genes=data.n_vars, n_new_genes=len(ids),
                    all_source_libraries_exact=True, original_source_axis_reused_exactly=True,
                    decoded_new_vectors_equal_native_counts=True,
                    n_old_gene_atom_mouse_checks=len(check), max_abs_native_log_sum_error=float(log_error.max()),
                    old_count_and_detection_agreement='EXACT',
                    n_independent_new_gene_atom_mouse_integer_checks=len(checked),
                    independent_all_new_gene_integer_sum_agreement='EXACT',
                    population_index_sha256=sha_file(prepared / 'POPULATION_INDEX.tsv'),
                    population_memberships_sha256=sha_file(prepared / 'POPULATION_MEMBERSHIPS.tsv.gz'),
                    source_prepare_receipt_sha256=sha_file(receipt_file)))
                print(json.dumps(dict(branch=branch, status='PASS', new_genes=len(ids),
                    independent_integer_checks=len(checked))), flush=True)
            finally:
                data.file.close()

    pooled = add_reference_scaling(add_mask_and_contrasts(pd.concat(pooled, ignore_index=True), annotations))
    wt = add_reference_scaling(add_mask_and_contrasts(pd.concat(wt, ignore_index=True), annotations))
    pooled.to_csv(out / 'native/new_gene_cluster_summary.tsv', sep='\t', index=False, na_rep='NA')
    wt.to_csv(out / 'native/new_gene_cluster_summary_WT.tsv', sep='\t', index=False, na_rep='NA')
    for ds in manifest['datasets']:
        select = lambda x: x.loc[x.tool.eq(ds['tool']) & x.branch.eq(ds['branch']) & x.population.eq(ds['population'])]
        p, w = select(pooled), select(wt)
        w_lookup = {(r.cluster_id, r.gene_id): r._asdict() for r in w.itertuples(index=False)}
        payload = dict(schema='satb2-user-marker-extension-1', dataset_id=ds['id'],
            tool=ds['tool'], branch=ds['branch'], population=ds['population'],
            cluster_scope='original_frozen', gene_ids=ids,
            base_marker_shard=ds['markers_shard'],
            base_marker_payload_sha256=manifest['shards'][ds['markers_shard']]['payload_sha256'],
            rows=[], WT_rows=[], cKO_rows=[], scientific_model_fitting_performed=False)
        for _, row in p.iterrows():
            ar, wr = marker_row(row), marker_row(w_lookup[(row.cluster_id, row.gene_id)])
            payload['rows'].append(ar); payload['WT_rows'].append(wr); payload['cKO_rows'].append(cko_row(ar, wr))
        for scope in ['rows', 'WT_rows', 'cKO_rows']:
            keys = {(r['cluster_id'], r['gene_id']) for r in payload[scope]}
            if len(keys) != len(payload[scope]) or len(payload[scope]) != p.cluster_id.nunique() * len(ids):
                raise ValueError('Incomplete or duplicated new marker atom/gene scope')
        ident = 'markers_v34__' + ds['id']
        catalogs['shards'][ident] = emit(marker_dir, ident, payload)
        catalogs['datasets'].append(dict(id=ds['id'], tool=ds['tool'], branch=ds['branch'],
            population=ds['population'], cluster_scope='original_frozen', markers_shard=ident,
            base_markers_shard=ds['markers_shard']))
    if len(catalogs['datasets']) != 24 or len(vectors['shards']) != 3 * len(ids):
        raise ValueError('Incomplete new additive catalogs')
    for relative, original_sha in source_invariants.items():
        if sha_file(repo / relative) != original_sha:
            raise ValueError('Immutable existing catalog was modified during descriptive export')
    catalog_files = write_catalog(repo / 'site', 'MARKER_EXTENSION_CATALOG', 'SATB2_MARKER_EXTENSION_CATALOG', catalogs)
    catalog_files += write_catalog(repo / 'site', 'GENE_EXPRESSION_EXTENSION_CATALOG', 'SATB2_GENE_EXPRESSION_EXTENSION_CATALOG', vectors)
    receipt = dict(schema='satb2-user-marker-extension-export-1', status='COMPLETE', release_tag=TAG,
        scientific_model_fitting_performed=False, previous_scientific_results_modified=False,
        original_889_gene_marker_shards_modified=False, original_2670_expression_files_modified=False,
        original_source_axis_reused=True, new_descriptive_native_count_aggregation_performed=True,
        normalization=vectors['normalization'], zero_policy=vectors['zero_policy'],
        cKO_z_policy=catalogs['cKO_z_policy'], n_new_genes=len(ids), n_datasets=24,
        n_new_gene_vector_shards=len(vectors['shards']), n_new_marker_shards=len(catalogs['shards']),
        existing_catalog_sha256=source_invariants,
        producer_sha256=sha_file(__file__), new_genes_input_sha256=sha_file(args.new_genes),
        source_objects_manifest_sha256=sha_file(source_manifest_file),
        native_summary_code_sha256=sha_file(repo / 'analysis/v2/executed_workflow_snapshot/plots/summarize_aggregates.py'),
        initial_annotation_sha256=sha_file(annotations_file), branch_validation=branch_checks,
        new_catalog_sha256={str(p.relative_to(repo)):sha_file(p) for p in catalog_files},
        new_marker_shard_sha256={Path(d['path']).name:d['sha256'] for d in catalogs['shards'].values()},
        new_gene_vector_sha256={Path(d['path']).name:d['sha256'] for d in vectors['shards'].values()},
        native_output_sha256={str(p.relative_to(out)):sha_file(p) for p in sorted((out / 'native').glob('*'))})
    (out / 'COMPLETE.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k:receipt[k] for k in ['status', 'n_new_genes', 'n_datasets', 'n_new_gene_vector_shards']}), flush=True)


if __name__ == '__main__':
    main()
