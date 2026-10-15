#!/usr/bin/env python3
"""Add eleven display markers to frozen SATB2 native atomic populations.

This explicitly authorized extension reads existing integer count AnnData files,
audits their complete libraries, and computes descriptive marker sums. It never
fits QC, embeddings, clusters, annotations or differential expression. Existing
878-gene row objects are preserved exactly and eleven genes are appended.
"""
from __future__ import annotations
import argparse
import base64
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import struct
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ADDITIONS = 'Unc5c Cpne5 Nrxn1 Tcf4 Fgf12 Ptpre Gria4 Gria1 Ctnnd2 Mgat4c Ptprd'.split()
AUDIT_GENES = ['Ccnd2', 'Pvalb', 'Sst']
BRANCHES = ['original_full', 'paired_original', 'paired_cellbender']
SCOPES = ['rows', 'WT_rows', 'cKO_rows']
FIELDS = ['cluster_id','gene_id','gene_symbol','n_nuclei','sum_counts','n_detected',
          'sum_log1p_normalized','mean_log1p_normalized','fraction_detected',
          'mean_source_counts','scaled_expression_z','n_mice_with_nuclei',
          'n_mice_eligible_for_expression_review','n_mice_expression_supported',
          'pct_difference_vs_rest','mean_log1p_difference_vs_rest','keep_for_plot',
          'measurement_status','scale_reference_branch','scale_reference_mean',
          'scale_reference_sd','scale_definition']
STRING_FIELDS = {'cluster_id','gene_id','gene_symbol','measurement_status',
                 'scale_reference_branch','scale_definition'}

def sha_file(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def compact(x):
    return json.dumps(x, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()

def json_value(v):
    if pd.isna(v): return None
    if isinstance(v, np.generic): return v.item()
    return v

def marker_row(r):
    return {k: str(r[k]) if k in STRING_FIELDS else bool(r[k]) if k == 'keep_for_plot'
            else float(r[k]) if pd.notna(r[k]) else None for k in FIELDS}

def cko_row(a, w):
    r = {k:a[k] for k in ['cluster_id','gene_id','gene_symbol','keep_for_plot']}
    for k in ['n_nuclei','sum_counts','n_detected','sum_log1p_normalized',
              'n_mice_with_nuclei','n_mice_eligible_for_expression_review','n_mice_expression_supported']:
        r[k] = a[k] - w[k]
        if r[k] < -1e-7: raise ValueError('Negative cKO numerator: ' + k)
        if r[k] < 0: r[k] = 0.0
    n = r['n_nuclei']
    r.update(mean_log1p_normalized=r['sum_log1p_normalized']/n if n else None,
             fraction_detected=r['n_detected']/n if n else None,
             mean_source_counts=r['sum_counts']/n if n else None,
             scaled_expression_z=None, pct_difference_vs_rest=None,
             mean_log1p_difference_vs_rest=None,
             measurement_status='FROZEN_ATOMIC_MARKER_EXTENSION_POOLED_MINUS_WT;NO_NEW_MODEL',
             scale_reference_branch=None, scale_reference_mean=None,
             scale_reference_sd=None, scale_definition='UNAVAILABLE; cKO z not computed')
    return r

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', type=Path, required=True)
    ap.add_argument('--native-results', type=Path, required=True)
    ap.add_argument('--object-manifest', type=Path, required=True)
    ap.add_argument('--panel', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args(); repo = a.repo.resolve(); out = a.out.resolve()
    if out.exists() and any(out.iterdir()): raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    for d in ['hosted', 'offline', 'native', 'validation']: (out/d).mkdir()
    summary_dir = repo/'analysis/v2/executed_workflow_snapshot/plots'
    sys.path.insert(0, str(summary_dir))
    from summarize_aggregates import summarize, add_mask_and_contrasts, add_reference_scaling
    from plot_common import KEYS, read_table
    manifest = json.loads((repo/'site/data/manifest.json').read_text())
    base_registry = {g['gene_symbol']: g['gene_id'] for g in manifest['genes']}
    if any(g in base_registry for g in ADDITIONS): raise ValueError('Some extension genes already registered')
    panel_symbols = a.panel.read_text().split()
    if len(panel_symbols) != 75 or len(set(panel_symbols)) != 75: raise ValueError('Expected exact unique screenshot75')
    if set(panel_symbols).difference(base_registry) != set(ADDITIONS): raise ValueError('Panel missing gene list changed')
    object_inputs = json.loads(a.object_manifest.read_text())['files']
    sources = {x['branch']:x for x in object_inputs if x['kind'] == 'shared_counts'}
    annotations_path = repo/'analysis/data/annotation_initial/cluster_annotations.tsv'
    annotations = read_table(annotations_path)
    pooled, wt, mouse_parts, sources_audit, check_parts = [], [], [], [], []
    new_gene_ids = None
    audit_symbols = ADDITIONS + AUDIT_GENES
    for branch in BRANCHES:
        prepared = a.native_results/'prepared'/branch
        receipt = json.loads((prepared/'COMPLETE.json').read_text())
        for name in ['POPULATION_INDEX.tsv','POPULATION_MEMBERSHIPS.tsv.gz','cluster_mouse_counts.tsv']:
            if sha_file(prepared/name) != receipt['output_sha256'][name]: raise ValueError('Prepared input hash mismatch: '+name)
        source = sources[branch]; path = Path(source['local_path'])
        if sha_file(path) != source['sha256']: raise ValueError('Full count object hash mismatch')
        count_data = ad.read_h5ad(path, backed='r')
        try:
            if count_data.uns['count_branch'] != branch: raise ValueError('Native count branch mismatch')
            if not count_data.obs_names.is_unique or not count_data.var_names.is_unique: raise ValueError('Duplicate native IDs')
            symbol_col = count_data.var['gene_symbol'].astype(str)
            ids = []
            for g in audit_symbols:
                hits = count_data.var_names[symbol_col.eq(g)].tolist()
                if len(hits) != 1: raise ValueError('Gene symbol has nonunique/absent stable ID: ' + g)
                ids.append(hits[0])
            if new_gene_ids is None: new_gene_ids = dict(zip(ADDITIONS, ids[:len(ADDITIONS)]))
            elif new_gene_ids != dict(zip(ADDITIONS, ids[:len(ADDITIONS)])): raise ValueError('Gene axis differs across branches')
            gene_j = count_data.var_names.get_indexer(ids)
            n = count_data.n_obs
            selected = np.zeros((n, len(ids)), dtype=np.float64)
            library = np.empty(n, dtype=np.float64)
            for start in range(0, n, 8192):
                stop = min(start+8192, n); x = count_data.X[start:stop, :]
                if (x.data < 0).any() or not np.isfinite(x.data).all() or not np.equal(x.data, np.floor(x.data)).all():
                    raise ValueError('Native source contains noninteger/invalid values')
                library[start:stop] = np.asarray(x.sum(axis=1)).ravel()
                selected[start:stop, :] = x[:, gene_j].toarray()
            expected = count_data.obs['branch_UMIs_verified'].to_numpy(dtype=float)
            if not np.array_equal(library, expected) or (library <= 0).any(): raise ValueError('Full source library audit failed')
            obs = count_data.obs[['sample', 'condition']].copy()
            index = read_table(prepared/'POPULATION_INDEX.tsv')
            memberships = read_table(prepared/'POPULATION_MEMBERSHIPS.tsv.gz')
            counts = read_table(prepared/'cluster_mouse_counts.tsv')
            if memberships.duplicated(['population_key','nucleus_id']).any(): raise ValueError('Duplicate membership')
            oi = count_data.obs_names.get_indexer(memberships.nucleus_id)
            pi = pd.Index(index.population_key).get_indexer(memberships.population_key)
            if (oi < 0).any() or (pi < 0).any(): raise ValueError('Membership ID outside source')
            for k, r in index.iterrows():
                member_ids = sorted(memberships.loc[memberships.population_key.eq(r.population_key), 'nucleus_id'])
                digest = hashlib.sha256(('\n'.join(member_ids)+'\n' if member_ids else '').encode()).hexdigest()
                if len(member_ids) != r.n_nuclei or digest != r.membership_sha256: raise ValueError('Membership hash/count changed')
            assignment = sparse.csr_matrix((np.ones(len(oi)), (oi, pi)), shape=(n, len(index)))
            normalized = np.log1p(selected * (10000.0/library[:, None]))
            records = []
            # All twelve mice explicitly retained, including zero-nucleus atoms.
            for sample in sorted(obs['sample'].astype(str).unique()):
                mask = obs['sample'].astype(str).eq(sample).to_numpy()
                sample_assign = assignment[mask]
                nn = np.asarray(sample_assign.sum(axis=0)).ravel().astype(int)
                sums = np.asarray(sample_assign.T @ selected[mask])
                det = np.asarray(sample_assign.T @ (selected[mask] > 0).astype(float))
                logs = np.asarray(sample_assign.T @ normalized[mask])
                conditions = obs.loc[mask,'condition'].astype(str).unique().tolist()
                if len(conditions) != 1: raise ValueError('Conflicting source condition')
                for k, r in index.iterrows():
                    for j, (gid, sym) in enumerate(zip(ids, audit_symbols)):
                        records.append(dict(tool=r.tool, branch=branch, population=r.population,
                                            cluster_scope='original_frozen', cluster_id=r.cluster_id,
                                            sample=sample, condition=conditions[0], gene_id=gid,
                                            gene_symbol=sym, n_nuclei=int(nn[k]),
                                            sum_counts=float(sums[k,j]), n_detected=float(det[k,j]),
                                            sum_log1p_normalized=float(logs[k,j])))
            frame = pd.DataFrame(records)
            # Compare three already-packaged genes against original native R sums
            # for every atomic cluster and mouse (17,604 rows across branches).
            wanted = []
            for c in pd.read_csv(prepared/'gene_cluster_mouse_aggregates.tsv.gz', sep='\t',
                                 dtype={'cluster_id':str}, chunksize=150000):
                wanted.append(c.loc[c.gene_symbol.isin(AUDIT_GENES)])
            old_check = pd.concat(wanted, ignore_index=True)
            keys = KEYS + ['sample', 'gene_id']
            check = frame.loc[frame.gene_symbol.isin(AUDIT_GENES)].merge(old_check,
                on=keys, suffixes=('_fresh','_native'), validate='one_to_one')
            if len(check) != len(index)*12*len(AUDIT_GENES): raise ValueError('Native audit coverage missing')
            for field in ['n_nuclei', 'sum_counts', 'n_detected']:
                if not np.array_equal(check[field+'_fresh'], check[field+'_native']): raise ValueError('Exact native audit failed: '+field)
            errors = np.abs(check.sum_log1p_normalized_fresh - check.sum_log1p_normalized_native)
            if not np.allclose(check.sum_log1p_normalized_fresh, check.sum_log1p_normalized_native, rtol=2e-12, atol=1e-10):
                raise ValueError('Native full-library normalization audit failed')
            check_parts.append(dict(branch=branch, n_rows=len(check), integer_sums='EXACT',
                                    max_abs_log_sum_error=float(errors.max()),
                                    relative_tolerance=2e-12, absolute_tolerance=1e-10))
            frame = frame.loc[frame.gene_symbol.isin(ADDITIONS)].reset_index(drop=True)
            mouse, summaries = summarize(frame, counts, 5, .1)
            mouse_parts.append(mouse); pooled.append(summaries[0]); wt.append(summaries[1])
            sources_audit.append(dict(branch=branch, count_object_sha256=source['sha256'],
                count_object_release_asset=source['asset_name'], source_nuclei=n, source_genes=count_data.n_vars,
                all_source_library_sums_exact=True, atomic_memberships_sha256=sha_file(prepared/'POPULATION_MEMBERSHIPS.tsv.gz'),
                atomic_index_sha256=sha_file(prepared/'POPULATION_INDEX.tsv'),
                source_prepare_receipt_sha256=sha_file(prepared/'COMPLETE.json')))
            print(json.dumps({'branch':branch,'n_mouse_gene_rows':len(mouse), 'native_spotcheck':'PASS',
                              'max_abs_log_sum_error':float(errors.max())}), flush=True)
        finally: count_data.file.close()
    pooled = add_reference_scaling(add_mask_and_contrasts(pd.concat(pooled, ignore_index=True), annotations))
    wt = add_reference_scaling(add_mask_and_contrasts(pd.concat(wt, ignore_index=True), annotations))
    pd.concat(mouse_parts, ignore_index=True).to_csv(out/'native/extension_gene_cluster_mouse_summary.tsv.gz', sep='\t',index=False,na_rep='NA')
    pooled.to_csv(out/'native/extension_gene_cluster_summary.tsv',sep='\t',index=False,na_rep='NA')
    wt.to_csv(out/'native/extension_gene_cluster_summary_WT.tsv',sep='\t',index=False,na_rep='NA')
    patch = {'genes':[{'gene_id':new_gene_ids[s], 'gene_symbol':s} for s in ADDITIONS],
             'gene_sets': {'Screenshot75_20261015':[(base_registry|new_gene_ids)[s] for s in panel_symbols]},
             'screenshot_symbols_in_order': panel_symbols,
             'screenshot_blocks':[{'label':f'Screenshot block {i+1}',
                 'gene_symbols':panel_symbols[i*5:(i+1)*5],
                 'gene_ids':[(base_registry|new_gene_ids)[s] for s in panel_symbols[i*5:(i+1)*5]]}
                 for i in range(15)], 'shards':{}}
    shard_checks = []
    for ds in manifest['datasets']:
        ident = ds['markers_shard']; original = repo/'site/data/shards'/f'{ident}.json.gz'
        payload = json.loads(gzip.decompress(original.read_bytes()))
        base_bytes = {k:compact(payload[k]) for k in SCOPES}
        orig_lengths = {k:len(payload[k]) for k in SCOPES}
        p = pooled.loc[pooled.tool.eq(ds['tool']) & pooled.branch.eq(ds['branch']) & pooled.population.eq(ds['population'])]
        w = wt.loc[wt.tool.eq(ds['tool']) & wt.branch.eq(ds['branch']) & wt.population.eq(ds['population'])]
        w_lookup = {(r.cluster_id,r.gene_id): r for r in w.itertuples(index=False)}
        for _, r in p.iterrows():
            wr = w_lookup[(r.cluster_id,r.gene_id)]._asdict()
            ar, wrow = marker_row(r), marker_row(wr)
            payload['rows'].append(ar); payload['WT_rows'].append(wrow); payload['cKO_rows'].append(cko_row(ar,wrow))
        for k in SCOPES:
            if compact(payload[k][:orig_lengths[k]]) != base_bytes[k]: raise ValueError('Existing marker row text/numbers changed')
            if len({(r['cluster_id'],r['gene_id']) for r in payload[k]}) != len(payload[k]): raise ValueError('Repeated extended marker key')
            clusters = {r['cluster_id'] for r in payload[k]}
            if len(payload[k]) != len(clusters)*889: raise ValueError('Expected every atomic cluster x889 genes')
        payload['marker_extension'] = dict(schema='satb2-screenshot-markers-1', appended_symbols=ADDITIONS,
            source='Existing full-gene integer counts, exact frozen atomic memberships and complete source library',
            new_model_fitting_performed=False, new_descriptive_count_aggregation_performed=True,
            existing_878_marker_rows_numeric_and_text_unchanged=True,
            cKO_z_policy='Unavailable exactly as original atlas; use mean expression for cKO')
        raw = compact(payload); packed = gzip.compress(raw, compresslevel=6, mtime=0)
        gz_path = out/'hosted'/f'{ident}.json.gz'; gz_path.write_bytes(packed)
        script = ('window.SATB2_SHARDS=window.SATB2_SHARDS||{};window.SATB2_SHARDS['+json.dumps(ident)+']={encoding:"gzip-base64",data:'+json.dumps(base64.b64encode(packed).decode())+'};\n').encode()
        js_path = out/'offline'/f'{ident}.js'; js_path.write_bytes(script)
        common = dict(payload_sha256=hashlib.sha256(raw).hexdigest(), gzip_bytes=len(packed),uncompressed_bytes=len(raw))
        patch['shards'][ident] = dict(hosted=dict(path=f'data/shards/{ident}.json.gz', encoding='gzip-bytes',sha256=sha_file(gz_path),**common),
                                     offline=dict(path=f'data/shards/{ident}.js',encoding='gzip-base64',sha256=sha_file(js_path),**common))
        shard_checks.append(dict(id=ident, original_gene_count=878, extended_gene_count=889,
                                 old_rows_per_scope=orig_lengths, added_rows_per_scope=len(p),
                                 existing_rows_canonical_json_sha256={k:hashlib.sha256(v).hexdigest() for k,v in base_bytes.items()},
                                 existing_rows_bitwise_and_text_preservation='PASS',
                                 source_gzip_sha256=sha_file(original), extended_gzip_sha256=sha_file(gz_path)))
    (out/'manifest_patch.json').write_text(json.dumps(patch, indent=2)+'\n')
    (out/'validation/native_aggregate_comparison.json').write_text(json.dumps(check_parts,indent=2)+'\n')
    (out/'validation/existing_marker_preservation.json').write_text(json.dumps(shard_checks,indent=2)+'\n')
    receipt = dict(status='COMPLETE', schema='satb2-screenshot-marker-extension-1',
        new_descriptive_marker_count_aggregation_performed=True, scientific_model_fitting_performed=False,
        QC_clustering_embeddings_annotation_and_DESeq2_unchanged=True, source_run_modified=False,
        n_datasets=len(shard_checks), n_original_genes=878,n_extended_genes=889,
        n_screenshot_genes=75, addition_symbols=ADDITIONS, native_comparison=check_parts,
        initial_annotation_sha256=sha_file(annotations_path), native_source_inputs=sources_audit,
        producer_sha256=sha_file(__file__), native_summary_code_sha256=sha_file(summary_dir/'summarize_aggregates.py'),
        output_sha256={str(p.relative_to(out)):sha_file(p) for p in sorted(out.rglob('*')) if p.is_file()})
    (out/'COMPLETE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','n_datasets','n_original_genes','n_extended_genes','n_screenshot_genes']}),flush=True)

if __name__ == '__main__': main()
