from pathlib import Path
import csv,math
D=Path(__file__).resolve().parents[1]
A={(r['population'],r['sample']):int(r['full_source_library_UMIs']) for r in csv.DictReader((D/'library_denominator_audit.tsv').open(),delimiter='\t')}
for r in csv.DictReader((D/'all_counts_long.tsv').open(),delimiter='\t'):
    assert int(r['total_gene_UMIs'])==A[r['population'],r['sample']]
    assert math.isclose(float(r['per_million_library_UMIs']),int(r['count'])/int(r['total_gene_UMIs'])*1e6,rel_tol=1e-10,abs_tol=1e-10)
print('Full-library CPM arithmetic: PASS')
