#!/usr/bin/env python3
"""Read-only exact eleven-gene audit against every native atomic count array."""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json
from pathlib import Path
import pandas as pd

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--extension',type=Path,required=True)
    p.add_argument('--native-results',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    frame=pd.read_csv(a.extension/'native/extension_gene_cluster_mouse_summary.tsv.gz',sep='\t',dtype={'cluster_id':str})
    genes=set(frame.gene_id)
    expected={(r.branch,r.tool,r.population,r.cluster_id,r.gene_id,r.sample):r.sum_counts for r in frame.itertuples()}
    checked=set(); rows=[]; cells=0
    for branch in ['original_full','paired_original','paired_cellbender']:
        folder=a.native_results/'prepared'/branch
        index=pd.read_csv(folder/'POPULATION_INDEX.tsv',sep='\t',dtype={'cluster_id':str})
        for atom in index.itertuples():
            source=Path(atom.counts_file)
            if sha(source)!=atom.counts_sha256:raise ValueError('Native source count checksum changed')
            matched=0
            with gzip.open(source,'rt') as h:
                for gene in csv.DictReader(h,delimiter='\t'):
                    if gene['gene_id'] not in genes:continue
                    matched+=1
                    for sample in frame['sample'].unique():
                        k=(branch,atom.tool,atom.population,atom.cluster_id,gene['gene_id'],sample)
                        if k in checked:raise ValueError('Duplicate count audit key')
                        value=int(gene[sample])
                        if value!=expected[k]:raise ValueError('Native atomic count differs: '+str(k))
                        checked.add(k);cells+=1
            if matched!=11:raise ValueError('Native count array did not contain all eleven genes')
        rows.append(dict(branch=branch,n_atomic_arrays=len(index),n_genes=11,
                         n_mice=12,new_gene_integer_count_comparison='EXACT_ALL'))
        print(json.dumps(rows[-1]),flush=True)
    if checked!=set(expected):raise ValueError('New gene audit coverage incomplete')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(dict(status='PASS',n_checked_count_values=cells,branches=rows,
        native_fullgene_count_arrays_read_only=True,new_scientific_execution=False,
        extension_complete_sha256=sha(a.extension/'COMPLETE.json'),producer_sha256=sha(Path(__file__))),indent=2)+'\n')

if __name__=='__main__':main()
