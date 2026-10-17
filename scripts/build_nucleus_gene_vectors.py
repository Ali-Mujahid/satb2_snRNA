#!/usr/bin/env python3
"""Export exact sparse per-nucleus display counts from immutable native objects.

No filtering, normalization model, embedding, annotation or DE fitting occurs.
The client evaluates log1p(10000 * count / complete source library), or displays
the exact source count. Unknown source IDs remain unavailable, never zero.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

import anndata as ad
import numpy as np
from scipy import sparse

BRANCHES = ['original_full', 'paired_original', 'paired_cellbender']
RAW_PREFIX = 'https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.1.0/gene_expression/'


def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''): h.update(b)
    return h.hexdigest()


def compact(obj):
    return json.dumps(obj, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()


def emit(out, ident, payload):
    raw = compact(payload)
    packed = gzip.compress(raw, compresslevel=6, mtime=0)
    file = out / (ident + '.json.gz')
    file.write_bytes(packed)
    return dict(path=RAW_PREFIX + file.name, encoding='gzip-bytes',
                sha256=hashlib.sha256(packed).hexdigest(),
                payload_sha256=hashlib.sha256(raw).hexdigest(),
                gzip_bytes=len(packed), uncompressed_bytes=len(raw))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', type=Path, required=True)
    ap.add_argument('--object-manifest', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    repo, out = args.repo.resolve(), args.out.resolve()
    if out.exists() and any(out.iterdir()): raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = repo / 'site/data/manifest.json'
    manifest = json.loads(manifest_path.read_text())
    registry = manifest['genes']
    ids = [g['gene_id'] for g in registry]
    if len(ids) != 889 or len(set(ids)) != 889: raise ValueError('Expected 889 unique registry IDs')
    inputs = json.loads(args.object_manifest.read_text())['files']
    sources = {x['branch']: x for x in inputs if x['kind'] == 'shared_counts'}
    catalog = dict(schema='satb2-nucleus-expression-catalog-1',
        release_tag='atlas-v3.1.0', genes=registry,
        normalization=dict(target_sum=10000, transform='log1p', denominator='complete source library',
            expression='log1p(10000 * exact_source_gene_count / full_source_library_count)'),
        zero_policy='A present source nucleus with no sparse entry has exact count zero; an absent source nucleus is unavailable.',
        source_scope='Full source count cohorts; atlas rendering retains its unchanged frozen coordinate and exclusion subsets.',
        scientific_model_fitting_performed=False,
        branches={}, shards={})
    checks = []
    for branch in BRANCHES:
        source = sources[branch]
        path = Path(source['local_path'])
        if sha_file(path) != source['sha256']: raise ValueError('Native source checksum mismatch: ' + branch)
        data = ad.read_h5ad(path, backed='r')
        try:
            if data.uns['count_branch'] != branch: raise ValueError('Native branch metadata mismatch')
            if not data.obs_names.is_unique or not data.var_names.is_unique: raise ValueError('Native source IDs not unique')
            columns = data.var_names.get_indexer(ids)
            if (columns < 0).any(): raise ValueError('Registry gene absent from source')
            if any(str(data.var.iloc[j].gene_symbol) != g['gene_symbol'] for j, g in zip(columns, registry)):
                raise ValueError('Native registry symbol mismatch')
            n = data.n_obs
            library = np.empty(n, dtype=np.int64)
            selected_parts = []
            for start in range(0, n, 8192):
                stop = min(n, start + 8192)
                full = sparse.csr_matrix(data.X[start:stop, :])
                if not np.isfinite(full.data).all() or (full.data < 0).any() or not np.equal(full.data, np.floor(full.data)).all():
                    raise ValueError('Invalid or noninteger native counts')
                sums = np.asarray(full.sum(axis=1)).ravel()
                if (sums <= 0).any() or (sums > 2**53 - 1).any(): raise ValueError('Unsafe full library counts')
                library[start:stop] = sums.astype(np.int64)
                selected_parts.append(full[:, columns].astype(np.int64))
            if not np.array_equal(library, data.obs.branch_UMIs_verified.to_numpy(dtype=np.int64)):
                raise ValueError('Full 32,285-gene library sums disagree with verified native source')
            selected = sparse.vstack(selected_parts, format='csc')
            selected.eliminate_zeros()
            selected.sort_indices()
            axis_id = 'axis__' + branch
            axis = dict(schema='satb2-nucleus-expression-axis-1', branch=branch,
                        nucleus_ids=data.obs_names.tolist(), library_counts=library.tolist())
            axis_desc = emit(out, axis_id, axis)
            catalog['shards'][axis_id] = axis_desc
            gene_shards = {}
            # The seed is fixed so validation is reproducible, independent of sparse nonzero layout.
            rng = np.random.default_rng(20261017)
            random_rows = np.sort(rng.choice(n, size=min(257, n), replace=False))
            random_cols = np.sort(rng.choice(len(ids), size=61, replace=False))
            reconstructed = np.zeros((len(random_rows), len(random_cols)), dtype=np.int64)
            source_sums = np.asarray(selected.sum(axis=0)).ravel().astype(np.int64)
            source_detected = np.diff(selected.indptr)
            emitted_sums, emitted_detected = [], []
            for j, gene in enumerate(registry):
                begin, end = selected.indptr[j:j+2]
                ix = selected.indices[begin:end].astype(np.int64)
                vals = selected.data[begin:end].astype(np.int64)
                if (vals <= 0).any() or (vals > 2**53-1).any() or (np.diff(ix) <= 0).any():
                    raise ValueError('Unsafe sparse gene count vector')
                ident = 'gene__' + branch + '__' + gene['gene_id']
                body = dict(schema='satb2-nucleus-gene-counts-1', branch=branch,
                    gene_id=gene['gene_id'], gene_symbol=gene['gene_symbol'], n_source_nuclei=n,
                    axis_payload_sha256=axis_desc['payload_sha256'], indices=ix.tolist(), counts=vals.tolist())
                desc = emit(out, ident, body)
                catalog['shards'][ident] = desc
                gene_shards[gene['gene_id']] = ident
                # Decode the written bytes, rather than checking in-memory arrays alone.
                decoded = json.loads(gzip.decompress((out / (ident + '.json.gz')).read_bytes()))
                emitted_sums.append(sum(decoded['counts']))
                emitted_detected.append(len(decoded['indices']))
                if j in random_cols:
                    lookup = dict(zip(decoded['indices'], decoded['counts']))
                    reconstructed[:, np.flatnonzero(random_cols == j)[0]] = [lookup.get(int(i), 0) for i in random_rows]
            if not np.array_equal(source_sums, emitted_sums) or not np.array_equal(source_detected, emitted_detected):
                raise ValueError('Whole-branch per-gene sums/detection changed after export')
            # Independent backed source reads of random rows and gene columns.
            direct = np.vstack([sparse.csr_matrix(data.X[int(i):int(i)+1, :])[:, columns[random_cols]].toarray()
                                for i in random_rows]).astype(np.int64)
            if not np.array_equal(direct, reconstructed): raise ValueError('Independent native nucleus/gene spotcheck failed')
            source_axis = {x: i for i, x in enumerate(data.obs_names)}
            coordinate_checks = []
            aggregate_checks = []
            # Every native coordinate ID must map exactly; populations and masks stay unchanged.
            for dataset in [d for d in manifest['datasets'] if d['branch'] == branch]:
                ident = dataset['coordinates_shard']
                desc = manifest['shards'][ident]
                f = repo / 'site' / desc['path']
                if sha_file(f) != desc['sha256']: raise ValueError('Existing coordinate payload checksum mismatch')
                coords = json.loads(gzip.decompress(f.read_bytes()))
                nuclei = coords['nuclei']
                cids = [r['nucleus_id'] for r in nuclei]
                if len(set(cids)) != len(cids): raise ValueError('Existing coordinate IDs not unique')
                if any(x not in source_axis for x in cids): raise ValueError('Coordinate source-axis ID not found')
                coordinate_checks.append(dict(dataset=dataset['id'], n_nuclei=len(cids),
                    all_coordinate_ids_present_exactly_once=True, coordinate_payload_sha256=desc.get('payload_sha256'),
                    coordinate_gzip_sha256=desc['sha256']))
                # A separate three-gene sum/mean audit against the existing atom summaries.
                marker_desc = manifest['shards'][dataset['markers_shard']]
                marker_file = repo / 'site' / marker_desc['path']
                if sha_file(marker_file) != marker_desc['sha256']: raise ValueError('Existing marker payload checksum mismatch')
                markers = json.loads(gzip.decompress(marker_file.read_bytes()))
                audits = [g for g in registry if g['gene_symbol'] in ['Ccnd2', 'Pvalb', 'Sst']]
                audit_js = [ids.index(g['gene_id']) for g in audits]
                rows = {}
                for r in nuclei: rows.setdefault(str(r['cluster_id']), []).append(source_axis[r['nucleus_id']])
                wanted = {(str(r['cluster_id']), r['gene_id']): r for r in markers['rows']}
                max_error = 0.0
                for cid, cell_indices in rows.items():
                    counts = selected[cell_indices, :][:, audit_js].toarray()
                    sums = counts.sum(axis=0)
                    detected = (counts > 0).sum(axis=0)
                    means = np.log1p(counts * (10000.0/library[cell_indices, None])).mean(axis=0)
                    for z, g in enumerate(audits):
                        old = wanted[(cid, g['gene_id'])]
                        if int(old['n_nuclei']) != len(cell_indices) or old['sum_counts'] != sums[z] or old['n_detected'] != detected[z]:
                            raise ValueError('Coordinate expression sums disagree with frozen marker summary')
                        err = abs(old['mean_log1p_normalized'] - means[z])
                        max_error = max(max_error, err)
                        if not np.isclose(old['mean_log1p_normalized'], means[z], rtol=2e-12, atol=1e-10):
                            raise ValueError('Per-nucleus normalization disagrees with native cluster marker mean')
                aggregate_checks.append(dict(dataset=dataset['id'], genes=[g['gene_symbol'] for g in audits],
                    n_cluster_gene_checks=len(rows)*len(audits), counts_and_detection='EXACT',
                    max_abs_mean_log1p_error=max_error, marker_gzip_sha256=marker_desc['sha256']))
            catalog['branches'][branch] = dict(axis_shard=axis_id, genes=gene_shards,
                n_source_nuclei=n, source_count_object_sha256=source['sha256'],
                source_count_object_release_asset=source['asset_name'])
            checks.append(dict(branch=branch, source_sha256=source['sha256'],
                source_release_asset=source['asset_name'], n_source_nuclei=n, n_source_genes=data.n_vars,
                exported_genes=len(ids), all_registry_ids_and_symbols_unique_and_exact=True,
                all_full_libraries_equal_native_verified_libraries=True,
                source_889_gene_sums_and_detection_equal_decoded_vectors=True,
                independent_backed_native_cell_gene_checks=int(direct.size),
                independent_backed_native_cell_gene_agreement='EXACT',
                coordinate_checks=coordinate_checks, native_summary_checks=aggregate_checks,
                compressed_gene_vector_bytes=sum(catalog['shards'][v]['gzip_bytes'] for v in gene_shards.values())))
            print(json.dumps({'branch': branch, 'status': 'PASS',
                'gene_gzip_bytes': checks[-1]['compressed_gene_vector_bytes'],
                'axis_gzip_bytes': axis_desc['gzip_bytes']}), flush=True)
        finally: data.file.close()
    if len(catalog['shards']) != 3*890: raise ValueError('Incomplete exported catalog')
    catalog_path = repo / 'site/data/GENE_EXPRESSION_CATALOG.json'
    if catalog_path.exists(): raise FileExistsError(catalog_path)
    catalog_path.write_bytes(compact(catalog) + b'\n')
    (repo / 'site/data/gene_expression_catalog.js').write_text('window.SATB2_GENE_EXPRESSION_CATALOG=' + compact(catalog).decode() + ';\n')
    receipt = dict(schema='satb2-nucleus-gene-vector-export-1', status='COMPLETE',
        source_manifest_sha256=sha_file(manifest_path), producer_sha256=sha_file(Path(__file__)),
        scientific_model_fitting_performed=False, source_count_objects_modified=False,
        coordinate_membership_annotation_QC_DESeq2_modified=False,
        descriptive_per_nucleus_count_export_performed=True,
        normalization=catalog['normalization'], zero_policy=catalog['zero_policy'],
        n_shards=len(catalog['shards']), n_genes=889,
        total_compressed_bytes=sum(s['gzip_bytes'] for s in catalog['shards'].values()),
        release_tag=catalog['release_tag'], branches=checks,
        catalog_sha256=sha_file(catalog_path),
        output_sha256={Path(s['path']).name:s['sha256'] for s in catalog['shards'].values()})
    (out / 'PROVENANCE.json').write_text(json.dumps(receipt, indent=2) + '\n')
    (repo / 'site/data/GENE_EXPRESSION_PROVENANCE.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status':'COMPLETE','n_shards':len(catalog['shards']),
        'total_gzip_bytes':receipt['total_compressed_bytes']}), flush=True)


if __name__ == '__main__': main()
