"""Independent Python check of every native R atom and final union, all32k genes×12 mice.

Scientific count aggregation executes only as a registered native workflow rule.
Raw shared .X integer counts are streamed; normalized expression is never used.
"""
from pathlib import Path
import argparse, gc, hashlib, json
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

def read(p):return pd.read_csv(p,sep='\t',dtype={'gene_id':str,'nucleus_id':str,'population_key':str},keep_default_na=False,na_values=['NA'])
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def arrays(path,genes,samples):
    z=read(path)
    if not np.array_equal(z.gene_id.to_numpy(),genes):raise ValueError('gene order')
    x=z[samples].to_numpy(dtype=np.float64)
    if not np.isfinite(x).all() or (x<0).any() or (x!=np.floor(x)).any() or (x>np.iinfo(np.int32).max).any():raise ValueError('invalid/noninteger count')
    return x.astype(np.int64)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--outdir',required=True);a=ap.parse_args()
    cfg=json.loads(Path(a.config).read_text());out=Path(a.outdir)
    if out.exists() and any(out.iterdir()):raise ValueError('No overwrite')
    out.mkdir(parents=True,exist_ok=True);sc=cfg['deseq2'];receipts=[];total=0;validated_sources={}
    for branch in cfg['branches']:
        prepared=Path(sc['prepared_root_by_branch'][branch]);finalized=Path(sc['finalized_root_by_branch'][branch])
        ix=read(prepared/'POPULATION_INDEX.tsv');mm=read(prepared/'POPULATION_MEMBERSHIPS.tsv.gz')
        for p in [prepared/'POPULATION_INDEX.tsv',prepared/'POPULATION_MEMBERSHIPS.tsv.gz',finalized/'POPULATION_INDEX.tsv',finalized/'POPULATION_ATOMS.tsv']:
            validated_sources[str(p.resolve())]=sha(p)
        samples=read(sc['sample_manifest_tsv_by_branch'][branch])['sample'].tolist()
        h5=Path(cfg['previous_run_root'])/'results/shared'/branch/'branch_counts.h5ad'
        validated_sources[str(h5.resolve())]=sha(h5)
        h=ad.read_h5ad(h5,backed='r');genes=h.var['gene_id'].astype(str).to_numpy() if 'gene_id' in h.var else h.var_names.astype(str).to_numpy()
        ids=h.obs_names.astype(str);positions=pd.Series(np.arange(h.n_obs),index=ids)
        if positions.index.has_duplicates:raise ValueError('Duplicate shared nucleus')
        keys=ix.population_key.tolist();ki={k:i for i,k in enumerate(keys)};si={s:i for i,s in enumerate(samples)}
        rr=positions.loc[mm.nucleus_id].to_numpy();cc=mm.population_key.map(ki).to_numpy()
        sm=h.obs['sample'].astype(str).map(si).to_numpy()
        assignment=sparse.csr_matrix((np.ones(len(rr),dtype=np.int64),(rr,cc*len(samples)+sm[rr])),shape=(h.n_obs,len(keys)*len(samples)))
        accum=np.zeros((h.n_vars,len(keys)*len(samples)),dtype=np.int64)
        for start in range(0,h.n_obs,2048):
            end=min(start+2048,h.n_obs);x=h.X[start:end,:].tocsr()
            if np.any(x.data<0) or np.any(x.data!=np.floor(x.data)):raise ValueError('Shared source noninteger')
            accum+=(x.T.astype(np.int64)@assignment[start:end,:]).toarray()
        h.file.close();del h,assignment;gc.collect()
        expected={}
        for i,row in enumerate(ix.itertuples(index=False)):
            if sha(row.counts_file)!=row.counts_sha256:raise ValueError('Count source changed')
            x=accum[:,i*len(samples):(i+1)*len(samples)];actual=arrays(row.counts_file,genes,samples)
            bad=int(np.count_nonzero(x!=actual));entries=int(x.size);total+=entries
            receipts.append(dict(branch=branch,scope='atomic',population_key=row.population_key,entries=entries,differing_entries=bad))
            if bad:raise ValueError(f'Independent atomic mismatch {branch}/{row.population_key}:{bad}')
            expected[row.population_key]=x
        fi=read(finalized/'POPULATION_INDEX.tsv');fa=read(finalized/'POPULATION_ATOMS.tsv')
        for row in fi.itertuples(index=False):
            if sha(row.counts_file)!=row.counts_sha256:raise ValueError('Final union source changed')
            atoms=fa.loc[fa.population_key.eq(row.population_key),'atomic_population_key'].tolist()
            if len(atoms)!=len(set(atoms)):raise ValueError('duplicate final atom')
            x=sum((expected[k] for k in atoms),np.zeros((len(genes),len(samples)),dtype=np.int64))
            actual=arrays(row.counts_file,genes,samples)
            bad=int(np.count_nonzero(x!=actual));entries=int(x.size);total+=entries
            receipts.append(dict(branch=branch,scope='final_union',population_key=row.population_key,entries=entries,differing_entries=bad))
            if bad:raise ValueError(f'Independent union mismatch {branch}/{row.population_key}:{bad}')
        del accum,expected;gc.collect();print('verified',branch,flush=True)
    pd.DataFrame(receipts).to_csv(out/'EXACT_COUNT_CHECKS.tsv',sep='\t',index=False)
    (out/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',exact_integer_counts_PASS=True,config_sha256=sha(a.config),validated_sources=validated_sources,n_population_arrays=len(receipts),n_entries_checked=total,n_differing_entries=0,check='independent Python raw shared .X streaming; all genes and all12 mice; final disjoint unions'),indent=2)+'\n')

if __name__=='__main__':main()
