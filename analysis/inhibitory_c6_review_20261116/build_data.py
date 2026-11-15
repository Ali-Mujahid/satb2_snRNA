#!/usr/bin/env python3
"""Derive C6 comparison data from existing frozen original_full Seurat outputs."""
import csv,json,collections,math
from pathlib import Path
meta=Path('/Users/mujahid357/Documents/codex_new/SATB2_SEURAT_INHIBITORY_DESEQ2_PACKAGE_20261109/01_nucleus_metadata/rendered_nucleus_metadata.tsv')
marker=Path('/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0/results/marker_summaries/gene_cluster_summary.tsv')
out=Path(__file__).with_name('data.json')
ids={'2','6','10','17'}; cluster=collections.Counter(); geno=collections.Counter(); sup=collections.Counter(); mice=collections.Counter(); mice_228=collections.Counter(); subclasses=collections.Counter(); probs=collections.defaultdict(list)
with meta.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  c=r['cluster_id'];
  if c not in ids: continue
  assert r['tool']=='seurat' and r['population']=='Inhibitory' and r['branch']=='original_full' and r['cluster_scope']=='original_frozen'
  cluster[c]+=1; geno[(c,r['condition'])]+=1; sup[(c,r['reference_supertype_name'])]+=1; subclasses[(c,r['reference_subclass_name'])]+=1
  key=(c,r['sample'],r['condition']); mice[key]+=1
  if r['reference_supertype_name'].startswith('0228 '):mice_228[key]+=1
  try: probs[c].append(float(r['reference_supertype_bootstrapping_probability']))
  except ValueError: pass
expected={'2':1307,'6':501,'10':331,'17':116}
assert dict(cluster)==expected,(cluster,expected)
assert (geno[('6','WT')],geno[('6','cKO')])==(214,287),geno
assert sup[('6','0228 Sst Gaba_15')]==237 and sup[('6','0229 Sst Gaba_16')]==90 and sup[('6','0232 Sst Gaba_19')]==62
assert subclasses[('6','053 Sst Gaba')]==488

def number(s):
 try:
  x=float(s)
  return x if math.isfinite(x) else None
 except (ValueError,TypeError):return None
markers=[]
with marker.open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['tool']!='seurat' or r['branch']!='original_full' or r['population']!='Inhibitory' or r['cluster_scope']!='original_frozen' or r['condition_scope']!='ALL_GENOTYPES' or r['cluster_id'] not in ids: continue
  c=r['cluster_id']; n=int(r['n_nuclei']); d=int(r['n_detected']); assert n==expected[c] and 0<=d<=n
  frac=number(r['fraction_detected']); assert frac is not None and abs(frac-d/n)<1e-9
  markers.append(dict(cluster='C'+c,gene=r['gene_symbol'],gene_id=r['gene_id'],n_nuclei=n,n_detected=d,fraction_detected=frac,mean_log1p_normalized=number(r['mean_log1p_normalized']),scaled_expression_z=number(r['scaled_expression_z'])))
markers.sort(key=lambda x:(x['gene'].lower(),x['gene_id'],int(x['cluster'][1:])))
assert len(markers)>3000,len(markers)
assert len({(r['cluster'],r['gene_id']) for r in markers})==len(markers)
result=dict(clusters=[dict(cluster='C'+c,n=cluster[c],wt=geno[(c,'WT')],cko=geno[(c,'cKO')]) for c in ['2','6','10','17']],
 supertypes=[dict(cluster='C'+c,supertype=s,n=n) for (c,s),n in sorted(sup.items(),key=lambda p:(int(p[0][0]),-p[1],p[0][1]))],
 mice=[dict(cluster='C'+c,sample=s,condition=condition,n=n,supertype_0228=mice_228[(c,s,condition)]) for (c,s,condition),n in sorted(mice.items(),key=lambda p:(int(p[0][0]),p[0][2],p[0][1]))],
 markers=markers,
 subclasses=[dict(cluster='C'+c,subclass=s,n=n) for (c,s),n in sorted(subclasses.items(),key=lambda p:(int(p[0][0]),-p[1],p[0][1]))],
 provenance=dict(tool='Seurat v5',branch='original_full',population='Inhibitory',cluster_scope='original_frozen',condition_scope='ALL_GENOTYPES',metadata_source=str(meta),marker_summary_source=str(marker),nuclei_unit='nucleus',biological_replicate='mouse',mean_expression='mean log1p normalized source counts per 10,000 full-library UMIs, saved original full-cohort summary',z_score='gene-scaled z of original full cohort per-cluster mean log1p-normalized expression; equal-cluster weight; population standard deviation; precomputed over all inhibitory clusters, not recalculated for four displayed clusters',fraction='n_detected / n_nuclei; pooled genotypes; detection is not intensity',reference='Saved existing MapMyCells labels and bootstrap support; reference-mapped labels are predictions, not anatomical validation',missing='Markers absent from this saved summary are unavailable, not necessarily zero; no imputation.'))
out.write_text(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n')
print('clusters',result['clusters']);print('supertype rows',len(result['supertypes']),'mouse rows',len(result['mice']),'markers',len(markers),'file bytes',out.stat().st_size)
for g in ['Sst','Cdh7','Pcdh11x','Grin3a','Trpc6','Reln','Npy','Nos1','Elfn1','Chodl','Chrna2']:
 a=[x for x in markers if x['gene']==g];print(g,[(x['cluster'],round(x['fraction_detected'],4)) for x in a])
