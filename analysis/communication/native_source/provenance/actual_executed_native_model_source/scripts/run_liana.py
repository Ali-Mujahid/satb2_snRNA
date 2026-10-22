"""Native LIANA scoring and endpoint-specific evidence joins; no new DE tests."""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import json
import shutil
import time
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd

KEY = ['source', 'target', 'ligand_complex', 'receptor_complex']
SEVEN = KEY + ['lr_logfc', 'lr_gmeans', 'magnitude_rank']


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def write_tsv(frame, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False, na_rep='NA')


def pairs_for(groups):
    ca1 = [x for x in groups if x.startswith('CA1_')]
    inhibitory = [x for x in groups if x.startswith('IN_')]
    return pd.DataFrame([(x, y) for x in ca1 for y in inhibitory]
                        + [(y, x) for x in ca1 for y in inhibitory],
                        columns=['source', 'target'])


def endpoint_evidence(config, group_meta, features):
    """Keep individual stable-ID evidence; ambiguous symbols never gain primary support."""
    assert 'analysis_symbol' in features, 'Prepared feature mapping contract mismatch'
    feature_counts = features.groupby('analysis_symbol').size().to_dict()
    endpoint_by_key = {x['population_key']: x for x in config['endpoints']}
    evidence, exports, snapshots = {}, [], []
    for group, meta in group_meta.items():
        key = meta['endpoint_key']
        endpoint = endpoint_by_key[key]
        tested = str(endpoint['status']).startswith('TESTED')
        mapping = {}
        if tested:
            audit = pd.read_csv(endpoint['stable_audit_path'], sep='\t')
            audit['symbol'] = audit['source_gene_symbol'].fillna('').astype(str)
            for symbol, rows in audit.groupby('symbol', sort=False):
                if not symbol:
                    continue
                details = []
                for row in rows.to_dict('records'):
                    details.append({
                        'gene_id': str(row['gene_id']),
                        'log2FoldChange': None if pd.isna(row['log2FoldChange']) else float(row['log2FoldChange']),
                        'padj': None if pd.isna(row['padj']) else float(row['padj']),
                    })
                ambiguous = feature_counts.get(symbol, 0) != 1 or len(rows) != 1
                significant = bool((rows['padj'] < config['fdr_cutoff']).fillna(False).any())
                mapping[symbol] = {'symbol': symbol, 'tested_endpoint': True,
                                   'ambiguous_symbol': bool(ambiguous),
                                   'significant_any_stable_id': significant,
                                   'primary_deg_support': significant and not ambiguous,
                                   'stable_id_results': details}
        evidence[group] = {'endpoint_key': key, 'endpoint_tested': tested, 'genes': mapping}
        for field in ['all_results_path', 'deg_path', 'stable_audit_path']:
            path = Path(endpoint[field])
            if path.exists():
                expected = endpoint[field + '_sha256']
                assert sha(path) == expected, f'Changed DESeq2 source {path}'
                exports.append((group, field, path))
        snapshots.append({'group': group, 'endpoint_key': key, 'status': endpoint['status'],
                          'reason': endpoint.get('reason', ''), 'native_DEG_FDR05': endpoint['DEG_FDR05']})
    return evidence, exports, pd.DataFrame(snapshots)


def complex_evidence(group, complex_name, evidence, feature_counts):
    ep = evidence[group]
    details = []
    for symbol in str(complex_name).split('_'):
        entry = ep['genes'].get(symbol)
        if entry is None:
            entry = {'symbol': symbol, 'tested_endpoint': ep['endpoint_tested'],
                     'ambiguous_symbol': feature_counts.get(symbol, 0) != 1,
                     'significant_any_stable_id': False, 'primary_deg_support': False,
                     'stable_id_results': [], 'gene_result_available': False}
        details.append(entry)
    any_de = any(x['primary_deg_support'] for x in details)
    ambiguous = any(x['ambiguous_symbol'] for x in details)
    return any_de, ambiguous, json.dumps(details, separators=(',', ':'), allow_nan=False)


def annotate(frame, scope, condition, group_meta, evidence, feature_counts, detections):
    frame = frame.copy(deep=True)
    frame['scope'] = scope
    frame['condition'] = condition
    frame['direction'] = np.where(frame['source'].str.startswith('CA1_'), 'CA1_to_IN', 'IN_to_CA1')
    for side in ['source', 'target']:
        frame[side + '_label'] = frame[side].map(lambda x: group_meta[x]['label'])
        frame[side + '_review'] = frame[side].map(lambda x: group_meta[x]['identity_review'])
        frame[side + '_de_endpoint'] = frame[side].map(lambda x: evidence[x]['endpoint_key'])
        frame[side + '_de_tested'] = frame[side].map(lambda x: evidence[x]['endpoint_tested'])
    cache = {}
    for side, component in [('source', 'ligand'), ('target', 'receptor')]:
        keys = list(zip(frame[side], frame[component + '_complex']))
        for k in set(keys):
            if k not in cache:
                cache[k] = complex_evidence(*k, evidence, feature_counts)
        frame[component + '_deg_any'] = [cache[k][0] for k in keys]
        frame[component + '_de_symbol_ambiguous'] = [cache[k][1] for k in keys]
        frame[component + '_deseq2_components_json'] = [cache[k][2] for k in keys]
        frame[component + '_all_subunits_deg'] = [
            all(x['primary_deg_support'] for x in json.loads(cache[k][2])) for k in keys]
        frame[component + '_all_subunit_min_prop'] = [
            min(float(detections[group].get(gene, 0)) for gene in str(complex_name).split('_'))
            for group, complex_name in keys]
    frame['deg_supported'] = frame['ligand_deg_any'] | frame['receptor_deg_any']
    frame['both_sides_deg_supported'] = frame['ligand_deg_any'] & frame['receptor_deg_any']
    frame['lr_prop'] = frame[['ligand_all_subunit_min_prop', 'receptor_all_subunit_min_prop']].min(axis=1)
    return frame


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--prepared', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    config = json.loads(Path(args.config).read_text())
    prepared, out = Path(args.prepared), Path(args.outdir)
    assert not out.exists() or not any(out.iterdir()), 'Refusing existing nonempty output'
    out.mkdir(parents=True, exist_ok=True)
    assert importlib.metadata.version('liana') == config['liana_version']
    assert sha(config['resource_file']) == config['resource_sha256']
    resource = pd.read_csv(config['resource_file'], sep='\t')
    assert len(resource) == 3989 and not resource.duplicated().any()
    assert Path(prepared / 'COMPLETE.json').exists()
    expression = ad.read_h5ad(prepared / 'lr_expression.h5ad')
    assert expression.n_obs == config['expected_nuclei']
    assert np.all(np.isfinite(expression.X.data)) and np.all(expression.X.data >= 0)
    meta_frame = pd.read_csv(prepared / 'group_metadata.tsv', sep='\t').fillna('')
    group_meta = meta_frame.set_index('group').to_dict('index')
    features = pd.read_csv(prepared / 'feature_mapping.tsv', sep='\t', keep_default_na=False)
    evidence, exports, endpoint_summary = endpoint_evidence(config, group_meta, features)
    feature_counts = features.groupby('analysis_symbol').size().to_dict()
    for group, kind, path in exports:
        dest = out / 'deseq2_source_exports' / group
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest / path.name)
    write_tsv(endpoint_summary, out / 'tables/DESEQ2_ENDPOINTS.tsv')
    (out / 'tables/deseq2_component_evidence.json').write_text(json.dumps(evidence, indent=2, allow_nan=False))
    native_dir = out / 'native_scores'
    native_dir.mkdir()
    all_frames, mouse_frames, elig_frames, condition_pairs, timings = [], [], [], [], []
    for scope in ['IN_clusters', 'IN_total']:
        data = expression.copy()
        data.obs['ccc_group'] = data.obs['ccc_group'].astype(str)
        if scope == 'IN_total':
            data.obs.loc[data.obs['ccc_group'].str.startswith('IN_'), 'ccc_group'] = 'IN_Total'
        data.obs['ccc_group'] = data.obs['ccc_group'].astype('category')
        for condition in ['pooled', 'WT', 'cKO']:
            subset = data.copy() if condition == 'pooled' else data[data.obs['condition'].astype(str) == condition].copy()
            subset.obs['ccc_group'] = subset.obs['ccc_group'].cat.remove_unused_categories()
            pairs = pairs_for(list(subset.obs['ccc_group'].cat.categories))
            n_per_group = subset.obs['ccc_group'].astype(str).value_counts()
            for row in pairs.to_dict('records'):
                condition_pairs.append(dict(scope=scope, condition=condition,
                                            source=row['source'], target=row['target'],
                                            source_nuclei=int(n_per_group[row['source']]),
                                            target_nuclei=int(n_per_group[row['target']])))
            options = dict(groupby='ccc_group', resource=resource, resource_name='mouseconsensus',
                           expr_prop=config['expr_prop'], min_cells=config['min_cells'],
                           groupby_pairs=pairs, base=np.e, use_raw=False, layer=None,
                           seed=config['seed'], n_jobs=config['threads'],
                           return_all_lrs=False, inplace=False, verbose=True)
            start = time.time()
            print('START_CONSENSUS', scope, condition, subset.n_obs, flush=True)
            native = li.mt.rank_aggregate(subset, **options, n_perms=config['n_perms'],
                                          supp_columns=['ligand_logfc', 'receptor_logfc'])
            assert isinstance(native, pd.DataFrame) and not native.empty
            native = native.copy(deep=True)
            assert not native.duplicated(KEY).any()
            write_tsv(native, native_dir / f'{scope}_{condition}_rank_aggregate.tsv.gz')
            gm = li.mt.geometric_mean(subset, **options, n_perms=None)
            assert not gm.duplicated(KEY).any()
            write_tsv(gm, native_dir / f'{scope}_{condition}_geometric_mean.tsv.gz')
            merged = native.merge(gm[KEY + ['lr_gmeans']], on=KEY, how='left', validate='one_to_one')
            assert len(merged) == len(native) and merged['lr_gmeans'].notna().all()
            assert merged['magnitude_rank'].between(0, 1).all()
            detections = {}
            for group in subset.obs['ccc_group'].cat.categories:
                mask = subset.obs['ccc_group'].astype(str).to_numpy() == group
                detections[group] = pd.Series(np.asarray((subset.X[mask] > 0).mean(axis=0)).ravel(),
                                              index=subset.var_names)
            annotated = annotate(merged, scope, condition, group_meta, evidence, feature_counts, detections)
            assert annotated['lr_prop'].ge(config['expr_prop'] - 1e-7).all()
            all_frames.append(annotated)
            timings.append({'scope': scope, 'condition': condition, 'n_nuclei': subset.n_obs,
                            'n_interactions': len(native), 'seconds': time.time() - start})
            print('COMPLETE_CONSENSUS', scope, condition, len(native), round(time.time() - start, 1), flush=True)
            del subset, native, gm, merged, annotated, detections
            gc.collect()
        # Expression/detection support is estimated independently within each mouse.
        for sample in sorted(data.obs['sample'].astype(str).unique()):
            subset = data[data.obs['sample'].astype(str) == sample].copy()
            counts = subset.obs['ccc_group'].astype(str).value_counts()
            eligible = sorted(counts[counts >= config['min_cells']].index)
            pairs = pairs_for(eligible)
            conds = subset.obs['condition'].astype(str).unique()
            assert len(conds) == 1
            condition = conds[0]
            for row in pairs.to_dict('records'):
                elig_frames.append(dict(scope=scope, sample=sample, condition=condition,
                                        source=row['source'], target=row['target'],
                                        source_nuclei=int(counts[row['source']]), target_nuclei=int(counts[row['target']])))
            subset.obs['ccc_group'] = subset.obs['ccc_group'].cat.remove_unused_categories()
            print('MOUSE_SUPPORT', scope, sample, len(pairs), flush=True)
            if not pairs.empty:
                gm = li.mt.geometric_mean(subset, groupby='ccc_group', resource=resource,
                                         expr_prop=config['expr_prop'], min_cells=config['min_cells'],
                                         groupby_pairs=pairs, base=np.e, use_raw=False,
                                         n_perms=None, n_jobs=config['threads'], inplace=False,
                                         return_all_lrs=False, verbose=False)
                gm['scope'], gm['sample'], gm['condition'] = scope, sample, condition
                mouse_frames.append(gm)
            del subset
            gc.collect()
        del data
        gc.collect()
    combined = pd.concat(all_frames, ignore_index=True)
    mouse = pd.concat(mouse_frames, ignore_index=True)
    eligibility = pd.DataFrame(elig_frames)
    support = mouse.groupby(['scope', 'condition'] + KEY, observed=True)['sample'].nunique().rename('support').reset_index()
    eligible = eligibility.groupby(['scope', 'condition', 'source', 'target'])['sample'].nunique().rename('eligible').reset_index()
    for condition in ['WT', 'cKO']:
        s = support[support['condition'] == condition].drop(columns='condition').rename(columns={'support': 'sample_support_' + condition})
        e = eligible[eligible['condition'] == condition].drop(columns='condition').rename(columns={'eligible': 'sample_eligible_' + condition})
        combined = combined.merge(s, on=['scope'] + KEY, how='left', validate='many_to_one')
        combined = combined.merge(e, on=['scope', 'source', 'target'], how='left', validate='many_to_one')
        combined['sample_support_' + condition] = combined['sample_support_' + condition].fillna(0).astype(int)
        combined['sample_eligible_' + condition] = combined['sample_eligible_' + condition].fillna(0).astype(int)
        assert (combined['sample_support_' + condition] <= combined['sample_eligible_' + condition]).all()
    combined['support_ge3_either_genotype'] = (combined['sample_support_WT'] >= 3) | (combined['sample_support_cKO'] >= 3)
    combined['support_ge3_both_genotypes'] = (combined['sample_support_WT'] >= 3) & (combined['sample_support_cKO'] >= 3)
    assert not combined.duplicated(['scope', 'condition'] + KEY).any()
    # All three views are projections of one directed scoring universe.
    for scope in ['IN_clusters', 'IN_total']:
        for condition in ['pooled', 'WT', 'cKO']:
            f = combined[(combined.scope == scope) & (combined.condition == condition)]
            for direction in ['CA1_to_IN', 'IN_to_CA1', 'bidirectional']:
                d = f if direction == 'bidirectional' else f[f.direction == direction]
                d = d.sort_values(['magnitude_rank'] + KEY)
                prefix = out / 'tables' / scope / condition / direction
                write_tsv(d[SEVEN], str(prefix) + '_rank_aggregate.tsv')
                write_tsv(d.loc[d.deg_supported, SEVEN], str(prefix) + '_DESeq2_supported_rank_aggregate.tsv')
    write_tsv(combined, out / 'tables/interaction_evidence.tsv.gz')
    write_tsv(mouse, out / 'tables/mouse_expression_support.tsv.gz')
    write_tsv(eligibility, out / 'tables/mouse_pair_eligibility.tsv')
    write_tsv(pd.DataFrame(condition_pairs), out / 'tables/condition_pair_eligibility.tsv')
    write_tsv(pd.DataFrame(timings), out / 'tables/scoring_timings.tsv')
    # Comparison is descriptive and aligns exact LR keys; absence remains NA.
    rows = []
    compare_cols = ['lr_gmeans', 'magnitude_rank', 'specificity_rank', 'cellphone_pvals']
    for scope in ['IN_clusters', 'IN_total']:
        frames = []
        for condition in ['WT', 'cKO']:
            f = combined[(combined.scope == scope) & (combined.condition == condition)]
            frames.append(f[KEY + [x for x in compare_cols if x in f]].rename(columns={x: condition + '_' + x for x in compare_cols}))
        comparison = frames[0].merge(frames[1], on=KEY, how='outer', validate='one_to_one')
        comparison['scope'] = scope
        comparison['descriptive_gmean_difference_cKO_minus_WT'] = comparison['cKO_lr_gmeans'] - comparison['WT_lr_gmeans']
        comparison['genotype_interaction_pvalue'] = np.nan
        comparison['genotype_interaction_test_status'] = 'NOT_TESTED_DESCRIPTIVE_COMPARISON_ONLY'
        rows.append(comparison)
    write_tsv(pd.concat(rows, ignore_index=True), out / 'tables/WT_cKO_descriptive_score_comparison.tsv.gz')
    settings = {k: config[k] for k in ['liana_version', 'resource_name', 'resource_sha256', 'expr_prop', 'min_cells', 'n_perms', 'seed', 'threads', 'rank_cutoff', 'fdr_cutoff']}
    (out / 'SCORING_SETTINGS.json').write_text(json.dumps(settings, indent=2))
    counts = combined.assign(priority=combined.magnitude_rank <= config['rank_cutoff']).groupby(['scope', 'condition', 'direction']).agg(
        expression_qualified=('magnitude_rank', 'size'), rank_priority=('priority', 'sum'),
        DEG_supported=('deg_supported', 'sum'), both_sides_DEG=('both_sides_deg_supported', 'sum')).reset_index()
    priority = combined[(combined.magnitude_rank <= config['rank_cutoff']) & combined.deg_supported]
    pc = priority.groupby(['scope', 'condition', 'direction']).size().rename('DEG_supported_rank_priority').reset_index()
    counts = counts.merge(pc, on=['scope', 'condition', 'direction'], how='left').fillna({'DEG_supported_rank_priority': 0})
    write_tsv(counts, out / 'tables/RESULT_COUNTS.tsv')
    (out / 'COMPLETE.json').write_text(json.dumps({'status': 'SCORED_REVIEW_REQUIRED', 'rows': len(combined),
                                                'mouse_score_rows': len(mouse), 'settings': settings}, indent=2))
    print('COMPLETE', len(combined), 'mouse_support_rows', len(mouse), flush=True)


if __name__ == '__main__':
    main()
