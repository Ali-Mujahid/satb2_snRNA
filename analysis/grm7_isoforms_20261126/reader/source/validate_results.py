from pathlib import Path
import csv,json,math,collections
B=Path(__file__).resolve().parent;R=B/'run_v2/results';O=B/'deliverables'
manifest=json.loads((R/'COMPLETE.json').read_text());assert manifest['n_mice']==12 and manifest['n_WT']==manifest['n_cKO']==6
ann=json.loads((R/'annotation.json').read_text());assert len(ann['transcripts'])==5
coverage=list(csv.DictReader((R/'population_mouse_coverage.tsv').open(),delimiter='\t'));assert len(coverage)==411
by={(r['population'],r['sample']):r for r in coverage}
allcounts={}
for r in csv.DictReader((R/'all_counts_long.tsv').open(),delimiter='\t'):
 k=(r['population'],r['sample'],r['metric'],r['feature']);assert k not in allcounts;allcounts[k]=int(r['count']);assert int(r['count'])>=0
 assert math.isclose(float(r['per_1000_nuclei']),int(r['count'])/int(r['nuclei'])*1000,rel_tol=1e-10,abs_tol=1e-10)
for p,s in by:
 gene=allcounts[p,s,'gene','gene'];assert gene==int(by[p,s]['Grm7_locus_UMIs'])
 for metric in ['composition','equivalence_class']:assert sum(n for (g,ss,m,f),n in allcounts.items() if g==p and ss==s and m==metric)==gene
 for t in ann['transcripts']:
  strict=allcounts[p,s,'splice_supported_unique',t['name']];assert strict<=allcounts[p,s,'diagnostic_feature',t['name']] and strict<=allcounts[p,s,'compatible_upper_bound',t['name']]
 for (g,ss,m,f),n in allcounts.items():
  if g==p and ss==s and m=='feature_UMIs':assert n<=gene
for s in manifest['mice']:
 for p,clusters in [('CA1_superficial',['CA1_C2','CA1_C3']),('CA1_deep',['CA1_C1','CA1_C5']),('CA1_total',['CA1_C'+str(i) for i in range(1,12)]),('IN_total',['IN_C'+str(i) for i in range(1,21)])]:
  assert allcounts[p,s,'gene','gene']==sum(allcounts.get((g,s,'gene','gene'),0) for g in clusters)
 a=json.loads((B/'run/results'/(s+'.json')).read_text());b=json.loads((R/(s+'.json')).read_text())
 key=lambda r:{(x['group'],x['metric'],x['feature']):x['count'] for x in r['counts'] if x['metric'] in ['gene','feature_UMIs','feature_records','records','diagnostic_feature']}
 assert key(a)==key(b),'Gene or feature counts changed unexpectedly'
 assert a['provenance']['index_sha256']==b['provenance']['index_sha256']
 assert a['provenance']['header_sha256']==b['provenance']['header_sha256']
assert sum(r['status']=='NOT_SAMPLED' for r in csv.DictReader((O/'all_mouse_population_coverage.tsv').open(),delimiter='\t'))==9
result=dict(status='PASS',mice=12,nuclei=36614,populations=35,sampled_mouse_population_combinations=411,not_sampled_combinations=9,transcripts=5,original_feature_counts_preserved=True,normalizations_verified=True,UMI_conservation_verified=True,reference_sha256=manifest['gtf_sha256'])
(O/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
