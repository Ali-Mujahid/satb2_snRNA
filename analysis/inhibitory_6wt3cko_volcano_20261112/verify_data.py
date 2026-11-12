from pathlib import Path
import json,pandas as pd,zipfile
p=Path('analysis/inhibitory_6wt3cko_volcano_20261112');d=json.loads((p/'summary.json').read_text());ss=d['summary'];assert len(ss)==21 and len(d['mouse'])==252
m=pd.DataFrame(d['mouse']);total=0
for x in ss:
 pop=x['population'];s=m[(m.population==pop)&m.selected];assert len(s)==9 and int(s[s.condition=='WT'].n_nuclei.sum())==x['WT_nuclei'] and int(s[s.condition=='cKO'].n_nuclei.sum())==x['cKO_nuclei']
 if x['status']=='FITTED':
  a=json.loads((p/'data'/f'{pop}.json').read_text());assert a['population']==pop and len(a['rows'])==len({z[0] for z in a['rows']});hits=[z for z in a['rows'] if z[7] is not None and z[7]<.05];assert len(hits)==x['n_deg'];assert sum(z[3]>0 for z in hits)==x['n_up'];assert sum(z[3]<0 for z in hits)==x['n_down']; both=[z for z in hits if z[9] is not None and z[9]<.05];assert len(both)==x['overlap_total'];assert sum(z[3]>0 and z[8]>0 for z in both)==x['overlap_up'];assert sum(z[3]<0 and z[8]<0 for z in both)==x['overlap_down'];assert sum((z[3]*z[8])<0 for z in both)==x['overlap_opposite'];assert len(hits)-len(both)==x['selected_only'];
  if pop!='TOTAL':total+=len(hits)
 else:assert pop in ('C18','C20') and not (p/'data'/f'{pop}.json').exists()
assert total==198 and next(x for x in ss if x['population']=='TOTAL')['n_deg']==60
with zipfile.ZipFile(p/'all_cluster_DESeq2_results.zip') as z:assert z.testzip() is None and len(z.namelist())==80
print('PASS: 21 populations, 252 mouse rows, exact 198 cluster calls, all overlaps/directions, 60 total and ZIP integrity')
