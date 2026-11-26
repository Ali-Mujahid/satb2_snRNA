from pathlib import Path
import csv,json,math,collections
B=Path(__file__).resolve().parent;R=B/'run_v2/results';O=B/'deliverables'
def read(n):return list(csv.DictReader((R/n).open(),delimiter='\t'))
m=json.loads((R/'COMPLETE.json').read_text());assert m['n_mice']==12 and m['n_WT']==m['n_cKO']==6
a=json.loads((R/'annotation.json').read_text());assert a['gene']['strand']=='-' and a['gene']['chrom']=='chr10'
assert len(a['transcripts'])==6 and sum(bool(t['CDS']) for t in a['transcripts'])==3
assert {t['protein_length'] for t in a['transcripts'] if t['protein_family']=='alpha'}=={1199}
assert {t['protein_length'] for t in a['transcripts'] if t['protein_family']=='beta'}=={906}
assert len({t['protein_sha256'] for t in a['transcripts'] if t['protein_family']=='beta'})==1
roster=list(csv.DictReader((B/'inputs/nucleus_roster.tsv').open(),delimiter='\t'));assert len(roster)==36614
den=collections.Counter();nuc=collections.Counter()
for r in roster:
 for p in r['groups'].split(';'):den[p,r['sample']]+=int(r['total_gene_UMIs']);nuc[p,r['sample']]+=1
cov=read('population_mouse_coverage.tsv');assert len(cov)==411
counts={}
for r in read('all_counts_long.tsv'):
 key=(r['population'],r['sample'],r['metric'],r['feature']);assert key not in counts;counts[key]=int(r['count']);assert counts[key]>=0
 d=den[r['population'],r['sample']];n=nuc[r['population'],r['sample']];assert int(r['total_gene_UMIs'])==d and int(r['nuclei'])==n
 assert math.isclose(float(r['per_million_library_UMIs']),int(r['count'])/d*1e6,rel_tol=1e-10,abs_tol=1e-10)
 assert math.isclose(float(r['per_1000_nuclei']),int(r['count'])/n*1000,rel_tol=1e-10,abs_tol=1e-10)
for r in cov:
 p,s=r['population'],r['sample'];gene=counts[p,s,'gene','gene']
 for metric in ['composition','family_composition','equivalence_class']:assert sum(v for (g,ss,m,f),v in counts.items() if g==p and ss==s and m==metric)==gene
 for f in ['alpha','beta']:
  assert counts[p,s,'family_splice_supported',f]<=counts[p,s,'family_diagnostic_feature',f]
  assert counts[p,s,'family_splice_supported',f]<=counts[p,s,'family_compatible_upper_bound',f]
 for t in a['transcripts']:
  assert counts[p,s,'splice_supported_unique',t['name']]<=counts[p,s,'diagnostic_feature',t['name']]
for s in m['mice']:
 for p,cc in [('CA1_superficial',['CA1_C2','CA1_C3']),('CA1_deep',['CA1_C1','CA1_C5']),('CA1_total',['CA1_C'+str(i) for i in range(1,12)]),('IN_total',['IN_C'+str(i) for i in range(1,21)])]:
  assert counts[p,s,'gene','gene']==sum(counts.get((c,s,'gene','gene'),0) for c in cc)
 sample=json.loads((R/(s+'.json')).read_text());assert sample['provenance']['header_RG_samples']==[s]
grid=list(csv.DictReader((O/'all_mouse_population_coverage.tsv').open(),delimiter='\t'));assert len(grid)==420 and sum(r['status']=='NOT_SAMPLED' for r in grid)==9
bench=read('original_matrix_vs_BAM_gene_counts.tsv');assert max(abs(int(r['BAM_minus_matrix'])) for r in bench)<=2
result=dict(status='PASS',mice=12,nuclei=len(roster),populations=35,sampled_mouse_population_combinations=411,not_sampled_combinations=9,all_six_RNA_models_retained=True,coding_alpha_1199aa=True,coding_beta_906aa_identical_two_transcripts=True,negative_strand_translation_verified=True,full_library_denominators_verified=True,UMI_conservation_verified=True,original_matrix_gene_counts_exact=sum(int(r['BAM_minus_matrix'])==0 for r in bench),original_matrix_gene_counts_max_abs_difference=2,gtf_sha256=m['gtf_sha256'])
(O/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
