#!/usr/bin/env python3
"""Integrate the frozen 11-gene descriptive extension into owned atlas derivatives.

Never modifies the original analyses or fitted results. Input outputs are SHA-bound.
"""
import argparse, copy, hashlib, json, shutil
from pathlib import Path

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, separators=(',', ':'), ensure_ascii=False)+'\n')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--extension', type=Path, required=True)
    ap.add_argument('--offline', type=Path, required=True)
    ap.add_argument('--baseline', type=Path, required=True)
    a=ap.parse_args(); root=Path(__file__).resolve().parents[1]; ext=a.extension
    receipt=json.loads((ext/'COMPLETE.json').read_text())
    assert receipt['status']=='COMPLETE' and not receipt['scientific_model_fitting_performed']
    for rel, digest in receipt['output_sha256'].items():
        assert sha(ext/rel)==digest, rel
    patch=json.loads((ext/'manifest_patch.json').read_text())
    assert len(patch['genes'])==11 and len(patch['shards'])==24
    baseline=json.loads(a.baseline.read_text())
    assert len(baseline['genes'])==878
    preset={'id':'screenshot_20261015_75', 'label':'Screenshot · all 75 genes in original order',
            'gene_ids':patch['gene_sets']['Screenshot75_20261015'],
            'gene_blocks':patch['screenshot_blocks'],
            'source':'User screenshot 2026-10-15 14:53:53; B1–B15 are image gene blocks, not atlas cluster identities.'}
    assert len(set(preset['gene_ids']))==75
    genes=baseline['genes']+patch['genes']; symbols={x['gene_id']:x['gene_symbol'] for x in genes}
    assert len(symbols)==889 and [symbols[x] for x in preset['gene_ids']]==patch['screenshot_symbols_in_order']
    metadata={'version':'atlas-v3.0.0','receipt':'data/MARKER_EXTENSION_PROVENANCE.json',
              'new_descriptive_marker_count_aggregation_performed':True,
              'new_scientific_model_fitting_performed':False,
              'original_marker_rows_preserved':878,'added_genes':patch['genes'],
              'preset_id':preset['id'],'source_fingerprint_unchanged':True,
              'QC_clustering_embeddings_annotations_and_DESeq2_unchanged':True}
    changes=[]
    for target, transport in [(root/'site','hosted'),(a.offline,'offline')]:
        current=json.loads((target/'data/manifest.json').read_text())
        # Baseline metadata may differ only in web versus classic transport and paths.
        assert len(current['genes']) in (878,889)
        current['genes']=copy.deepcopy(genes)
        current['gene_sets']=[preset]+[p for p in current['gene_sets'] if p['id']!=preset['id']]
        current['marker_extension']=metadata
        contract=current.setdefault('measurement_contract',{})
        contract['new_RNA_aggregation_performed']=True
        contract['aggregation_scope']='11 added descriptive marker summaries from unchanged full-gene source counts; native full-library normalization.'
        contract['new_scientific_fits_performed']=False
        for key, descriptors in patch['shards'].items():
            previous=current['shards'][key]
            prior=previous.get('v2_baseline', {k:previous[k] for k in ['sha256','payload_sha256','source_script_sha256'] if k in previous})
            desc={**previous, **descriptors[transport], 'v2_baseline':prior,
                  'marker_extension':'eleven appended genes; existing rows preserved'}
            if transport=='hosted': desc['source_script_sha256']=descriptors['offline']['sha256']
            current['shards'][key]=desc
            source=ext/transport/(key+('.json.gz' if transport=='hosted' else '.js'))
            destination=target/desc['path']; shutil.copyfile(source,destination)
            assert sha(destination)==desc['sha256']
        if transport=='hosted':
            old=current.get('web_transport',{})
            current['web_transport']={**old,'version':'atlas-web-gzip-3','release_tag':'atlas-v3.0.0',
                'scientific_payloads_changed':True,'native_gzip_bytes_preserved':False,
                'change_scope':'24 marker shards extended by 11 descriptive genes; all other shard bytes unchanged',
                'unchanged_exact_downloads_release_tag':'atlas-v2.1.0'}
        write_json(target/'data/manifest.json',current)
        (target/'data/manifest.js').write_text('window.SATB2_DATA='+json.dumps(current,separators=(',',':'),ensure_ascii=False)+';\n')
        write_json(target/'data/MARKER_EXTENSION_PROVENANCE.json',receipt)
        for rel in ['validation/existing_marker_preservation.json','validation/native_aggregate_comparison.json']:
            dest=target/'data/marker_extension'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ext/rel,dest)
        if transport=='offline':
            for source in (ext/'native').iterdir():
                dest=target/'data/marker_extension/native'/source.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest)
        changes.append({'transport':transport,'marker_shards':24,'genes':889,'preset_genes':75})
    history=root/'site/data/WEB_DATA_PROVENANCE.json'
    prior=json.loads(history.read_text())
    if 'historical_v2_transport' not in prior:
        write_json(history,{'status':'PASS','version':'atlas-v3.0.0','historical_v2_transport':prior,
            'current_extension':metadata,'n_total_shards':3140,'n_extended_marker_shards':24,
            'all_other_shards_byte_identical_to_v2':True,
            'new_scientific_fitting_performed':False,'RNA_descriptive_aggregation_performed':True})
    print(json.dumps({'status':'PASS','integrated':changes}))

if __name__=='__main__':main()
