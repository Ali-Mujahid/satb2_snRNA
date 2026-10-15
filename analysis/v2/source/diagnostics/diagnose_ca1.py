"""Native fixed-input CA1/ProS graph, PCA, Harmony and Leiden controls.

Controls never replace primary clustering or DE membership. Rotations/signs of
PCA are removed before comparing coordinates. Same seed integers across tools
do not guarantee identical random-number streams.
"""
from pathlib import Path
import argparse,json,gc,hashlib
import anndata as ad
import numpy as np,pandas as pd
import scanpy as sc,harmonypy
from scipy import sparse,linalg,io
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score

def read(p):return pd.read_csv(p,sep='\t',dtype={'nucleus_id':str,'cluster':str,'cluster_harmony':str,'cluster_unintegrated':str})
def ari(x,y):return float(adjusted_rand_score(x,y))
def knn(x):return NearestNeighbors(n_neighbors=20,n_jobs=4).fit(x).kneighbors(return_distance=False)
def near_overlap(x,y):
    a,b=knn(x),knn(y)
    return float(np.mean([len(set(i)&set(j))/20 for i,j in zip(a,b)]))
def latent_metrics(x,y):
    x=x-x.mean(0);y=y-y.mean(0);r,_=linalg.orthogonal_procrustes(x,y)
    return {'principal_angle_max_degrees':float(np.degrees(linalg.subspace_angles(x,y)).max()),
            'principal_angle_median_degrees':float(np.median(np.degrees(linalg.subspace_angles(x,y)))),
            'rotation_aligned_relative_error':float(linalg.norm(x@r-y)/linalg.norm(y)),
            'euclidean_knn20_mean_overlap':near_overlap(x,y)}
def graph_metrics(g,h):
    a=g.copy();a.setdiag(0);a.eliminate_zeros();a.data=np.ones_like(a.data)
    b=h.copy();b.setdiag(0);b.eliminate_zeros();b.data=np.ones_like(b.data)
    n=float(a.multiply(b).sum());u=float((a+b).astype(bool).sum())
    return {'directed_stored_edge_jaccard':n/u if u else 0,
            'scanpy_stored_self_loops':int(np.count_nonzero(g.diagonal())),
            'seurat_stored_self_loops':int(np.count_nonzero(h.diagonal())),
            'scanpy_mean_degree':float(np.diff(a.indptr).mean()),'seurat_mean_degree':float(np.diff(b.indptr).mean()),
            'mean_nucleus_neighbor_jaccard':float(np.mean([len(set(a.indices[a.indptr[i]:a.indptr[i+1]])&set(b.indices[b.indptr[i]:b.indptr[i+1]]))/max(1,len(set(a.indices[a.indptr[i]:a.indptr[i+1]])|set(b.indices[b.indptr[i]:b.indptr[i+1]]))) for i in range(a.shape[0])]))}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--graphs',required=True);ap.add_argument('--outdir',required=True);a=ap.parse_args()
    cfg=json.loads(Path(a.config).read_text());old=Path(cfg['previous_run_root'])/'results';graphs=Path(a.graphs);out=Path(a.outdir)
    if out.exists() and any(out.iterdir()):raise ValueError('No overwrite')
    out.mkdir(parents=True,exist_ok=True);latent=[];metrics=[];iterations=[];history=[]
    for branch in cfg['branches']:
        pyroot=old/'scanpy'/branch/'CA1_ProS';rroot=old/'seurat'/branch/'annotate_subsets/CA1_PROS_CANDIDATE'
        h=ad.read_h5ad(pyroot/'scanpy.h5ad',backed='r');ids=h.obs_names.astype(str);obs=h.obs.copy()
        py_md=read(pyroot/'nucleus_metadata.tsv').set_index('nucleus_id').loc[ids]
        r_md=read(rroot/'nucleus_metadata.tsv').set_index('nucleus_id').loc[ids]
        labels={}
        for rep,key in [('unintegrated','X_pca'),('harmony','X_pca_harmony')]:
            x=np.asarray(h.obsm[key],dtype=float);r=read(rroot/(('pca' if rep=='unintegrated' else 'harmony')+'.tsv')).set_index('nucleus_id').loc[ids].to_numpy(dtype=float)
            latent.append(dict(branch=branch,representation=rep,**latent_metrics(x,r)))
            g=h.obsp[f'neighbors_{rep}_connectivities'].tocsr()
            rid=pd.Index((graphs/f'{branch}_{rep}_ids.txt').read_text().splitlines());pos=rid.get_indexer(ids)
            if (pos<0).any() or len(rid)!=len(ids):raise ValueError('Graph vertex mismatch')
            rg=io.mmread(graphs/f'{branch}_{rep}.mtx').tocsr()[pos,:][:,pos].tocsr()
            metrics.append(dict(branch=branch,representation=rep,**graph_metrics(g,rg)))
            pc=f'cluster_{rep}';rc=f'CA1_PROS_CANDIDATE_{rep}_r0p4'
            for graph_name,graph in [('scanpy_connectivity',g),('seurat_SNN',rg)]:
                control=ad.AnnData(sparse.csr_matrix((len(ids),0)),obs=pd.DataFrame(index=ids))
                for ni in [2,10,30]:
                    sc.tl.leiden(control,adjacency=graph,resolution=.4,random_state=20261013,n_iterations=ni,flavor='igraph',directed=False,key_added='control')
                    z=control.obs.control.astype(str).to_numpy();labels[(rep,graph_name,ni)]=z
                    iterations.append(dict(branch=branch,representation=rep,graph=graph_name,runner='scanpy_igraph',n_iterations=ni,
                                           graph_conversion='all stored nonzero adjacency entries to undirected edges; symmetric entries create parallel pairs; stored diagonal loops retained',
                                           n_clusters=len(set(z)),ARI_vs_original_Scanpy=ari(z,py_md[pc]),ARI_vs_original_Seurat=ari(z,r_md[rc]),
                                           ARI_vs_same_graph_iter2=ari(z,labels[(rep,graph_name,2)])))
                    pd.DataFrame({'nucleus_id':ids,'cluster':z}).to_csv(out/f'{branch}_{rep}_{graph_name}_iter{ni}.tsv.gz',sep='\t',index=False)
                del control
            # Standardize conversion as single lower off-diagonal weighted
            # edges for both source graphs under the same Python runner.
            # Remaining source-graph differences also include upstream PCs.
            for graph_name,graph in [('canonical_scanpy_connectivity',g),('canonical_seurat_SNN',rg)]:
                control=ad.AnnData(sparse.csr_matrix((len(ids),0)),obs=pd.DataFrame(index=ids))
                canonical=sparse.tril(graph,k=-1,format='csr')
                sc.tl.leiden(control,adjacency=canonical,resolution=.4,random_state=20261013,n_iterations=10,flavor='igraph',directed=False,key_added='control')
                z=control.obs.control.astype(str).to_numpy();labels[(rep,graph_name,10)]=z
                iterations.append(dict(branch=branch,representation=rep,graph=graph_name,runner='scanpy_igraph_canonical_lower_no_loops',n_iterations=10,
                                       graph_conversion='lower off-diagonal once; one undirected edge per pair; loops removed; all vertices preserved',
                                       n_clusters=len(set(z)),ARI_vs_original_Scanpy=ari(z,py_md[pc]),ARI_vs_original_Seurat=ari(z,r_md[rc]),
                                       ARI_vs_same_graph_iter2=np.nan))
                pd.DataFrame({'nucleus_id':ids,'cluster':z}).to_csv(out/f'{branch}_{rep}_{graph_name}_iter10.tsv.gz',sep='\t',index=False)
                del control,canonical
            for ni in [2,10,30]:
                z=read(graphs/f'{branch}_{rep}_leiden{ni}.tsv').set_index('nucleus_id').loc[ids]['cluster'].astype(str).to_numpy()
                labels[(rep,'native_seurat',ni)]=z
                iterations.append(dict(branch=branch,representation=rep,graph='seurat_SNN',runner='native_seurat_igraph_GroupSingletons',n_iterations=ni,
                                       graph_conversion='lower off-diagonal once; one undirected edge per pair; loops removed; native singleton grouping',
                                       n_clusters=len(set(z)),ARI_vs_original_Scanpy=ari(z,py_md[pc]),ARI_vs_original_Seurat=ari(z,r_md[rc]),
                                       ARI_vs_same_graph_iter2=ari(z,labels[(rep,'native_seurat',2)])))
        # Repeat the exact earlier Harmony call from saved PCs, extending only
        # its iteration cap. This recovers missing Python objective history.
        input_pc=np.asarray(h.obsm['X_pca']);saved_pc=np.asarray(h.obsm['X_pca_harmony'])
        px=np.asarray(input_pc,dtype=float);saved=np.asarray(saved_pc,dtype=float);reruns={}
        for cap in [10,50]:
            hr=harmonypy.run_harmony(px,obs,'sample',max_iter_harmony=cap,theta=2,sigma=.1,lamb=1,nclust=100,max_iter_kmeans=4,epsilon_harmony=.01,ncores=4,random_state=20261013,verbose=False)
            corrected=np.asarray(hr.Z_corr);reruns[cap]=corrected
            objective=[float(v) for v in hr.objective_harmony]
            history.append(dict(branch=branch,max_iter_harmony=cap,completed_harmony_iterations=len(objective)-1,
                                input_PC_storage_dtype=str(input_pc.dtype),saved_Harmony_storage_dtype=str(saved_pc.dtype),
                                objective_history=json.dumps(objective),kmeans_rounds=json.dumps([int(v) for v in hr.kmeans_rounds]),
                                relative_error_vs_saved=float(linalg.norm(corrected-saved)/linalg.norm(saved)),
                                max_abs_error_after_matching_saved_dtype=float(np.max(np.abs(corrected.astype(saved_pc.dtype)-saved_pc))),
                                final_objective=objective[-1]))
        latent.append(dict(branch=branch,representation='Python_Harmony_cap10_vs50',**latent_metrics(reruns[10],reruns[50])))
        h.file.close();del h,labels,reruns;gc.collect();print('CA1 controls complete',branch,flush=True)
    pd.DataFrame(latent).to_csv(out/'LATENT_SPACE_COMPARISON.tsv',sep='\t',index=False)
    pd.DataFrame(metrics).to_csv(out/'NEIGHBOR_GRAPH_COMPARISON.tsv',sep='\t',index=False)
    pd.DataFrame(iterations).to_csv(out/'LEIDEN_FIXED_GRAPH_CONTROLS.tsv',sep='\t',index=False)
    pd.DataFrame(history).to_csv(out/'PYTHON_HARMONY_CONVERGENCE.tsv',sep='\t',index=False)
    (out/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',scope='CA1_ProS_fixed_saved_PC_and_graph_controls',primary_memberships_unchanged=True,
        interpretation='Iteration controls isolate iteration effects within each fixed graph and runner. Original Python adjacency conversion keeps symmetric parallel pairs and stored loops; native Seurat uses lower off-diagonal edges once and singleton grouping. Additional canonical-lower controls standardize graph conversion under the Python runner. Source graph differences still include upstream latent-space differences; graph metrics alone do not establish causality. Molecular atlas/marker evidence is assessed independently.'),indent=2)+'\n')

if __name__=='__main__':main()
