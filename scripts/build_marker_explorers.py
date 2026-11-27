"""Build lightweight readers over unchanged, existing Seurat marker shards."""
from pathlib import Path
import csv
import gzip
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / 'site'
AUDIT = ROOT / 'analysis/marker_explorers_20261127'
AUDIT.mkdir(parents=True, exist_ok=True)
CAT = json.loads((SITE/'data/marker_extension_catalog.js').read_text().split('=', 1)[1].rstrip(';\n'))
SOURCE = ROOT.parent/'SATB2_scanpy_seurat_mapmycells_20261014_revision/runs/20261014_v2_0_0/results/marker_summaries/representative_markers.tsv'
with SOURCE.open() as f:
    REPRESENTATIVES = list(csv.DictReader(f, delimiter='\t'))

CA1 = 'Ndst4 Col11a1 Zfp804a Hdac9 Kcnip4 Prkca Pde10a Sorcs3 Epha6 Ntrk2 Cdh10 Stac Man1a Hunk Cacna2d3 Nr3c2 Brinp1 R3hdm1 Spock3 Plekhg1 Srgap1 Slc24a2 Enox1 Garem1 Sik3 Ap2b1 Dock4 Stxbp5l Pip5k1b Nr4a1 Fam189a1 Zbtb16 Adgrb3 Ptprn Nr4a2 Malat1 Epha3 Kcnab1 Robo2 March1 Ppp3ca Arl15 Adcy2 Phactr1 Thsd7b Ano3 Grm1 Kcnc2 Tafa1 Nlgn1 Ralyl Mdga1 Gpc5 Galntl6 Ntm Fras1 Gabrg3 Opcml Pde4d Pde1a Ccser1 Plekha7 Lingo2 Pcdh7 Ptprk Asic2 Cntnap5a Alcam Dcc Ptprd Fstl4 Sorbs2 Fgf14 Erc2 Tenm4 Ryr2 Rora Nrg2 Neto1 Ak5 Rab3c Agbl4 Nav2 Cntn4 Arhgap26 Zfpm2 Zmat4 Smad3'.split()
IN = 'Slit2 Sox5 Sox6 Cntnap4 Fgf12 Mgat4c Sst Lypd6b Lypd6 Vip Sorcs3 Synpr Nrxn1 Sema5b Cacna2d1 Unc5c Ptprd Vcan Hapln1 Sema3c Tcf4 Gria4 Pvalb Grin3a Cdh7 Ctnnd2 Kit Egfr Ndnf Satb1 Glis3 Npas3 Rgs12 Syt6 Eps8 Cpne5 Npas1 Cxcl14 Ptprt Etl4 Gria1 Cux1 Cux2 Calb1 Mctp2 Mef2c Ptpre Unc5d Nectin3 Cnr1 Frem1 Itga8 Ntng1 Chrm2 Cpne4 Pde3a Stk32b Zfp57 Zbtb20 Ccbe1 Meis2 Ano2 Ano1 Cep112 Zic1 Zic4 Ephb1 Raver2 Nrg1 Nos1 Npy Rmst Rfx3'.split()
IN_LABELS = ['PV basket-like','SST O-LM-like','VIP IS-3-like','Ivy CGE-like','PV bistratified-like','SST heterogeneous','Ivy MGE-like','Pvalb/Sst mixed','CCK R–LM-like','SST hippocamposeptal-like','PV chandelier-like','Sncg/Ndnf mixed','RHP–COA Ndnf-like','CCK / IS-3 mixed','Trilaminar-like','IMN/Sst mixed','SST / Zbtb20 mixed','LSX-associated mixed','Sst Chodl/Nfib mixed','GABA/myeloid mixed']
IN_COLORS = ['#D84B08','#FA777D','#F5A544','#D8BC18','#A6A848','#00A884','#36ABAB','#75A6D4','#7878DB','#A552E2','#962DC8','#722081','#D585DF','#B985B2','#DC3EA1','#BF2A9D','#FF507C','#B951C7','#D523A7','#7E1847']


def sha(b): return hashlib.sha256(b).hexdigest()


def build(pop, slug, output, symbols):
    base = SITE/f'data/shards/markers__seurat__original_full__{pop}.json.gz'
    base_bytes = base.read_bytes()
    base_raw = gzip.decompress(base_bytes)
    d = json.loads(base_raw)
    ext_key = f'markers_v34__seurat__original_full__{pop}'
    desc = CAT['shards'][ext_key]
    ext_path = ROOT/'marker_extensions_v34'/f'{ext_key}.json.gz'
    ext_bytes = ext_path.read_bytes()
    ext_raw = gzip.decompress(ext_bytes)
    assert sha(ext_bytes) == desc['sha256'] and sha(ext_raw) == desc['payload_sha256']
    e = json.loads(ext_raw)
    sources = [dict(path=str(base.relative_to(SITE)), sha256=sha(base_bytes), payload_sha256=sha(base_raw)),
               dict(path=desc['path'], sha256=sha(ext_bytes), payload_sha256=sha(ext_raw))]
    rows = d['rows'] + e['rows']
    genes = {r['gene_id']: r['gene_symbol'] for r in rows}
    by_symbol = {s: g for g, s in genes.items()}
    assert len(by_symbol) == len(genes)
    assert len(symbols) == len(set(symbols))
    assert set(symbols) <= set(by_symbol), set(symbols) - set(by_symbol)
    ids = [by_symbol[s] for s in symbols]
    counts = {r['cluster_id']: int(r['n_nuclei']) for r in rows}
    clusters = []
    for c in sorted(counts, key=int):
        if pop == 'CA1_ProS':
            label = 'Superficial-associated' if c in ['2','3'] else 'Deep-associated' if c in ['1','5'] else 'CA1-prosubiculum' if c == '8' else 'Identity-conflicting' if c == '11' else 'Mixed / weak'
            color = '#d53c57' if c in ['2','3'] else '#3789c2' if c in ['1','5'] else '#9065b6' if c == '8' else '#a5b3bc' if c == '11' else '#f4a35f'
        else:
            label, color = IN_LABELS[int(c)-1], IN_COLORS[int(c)-1]
        clusters.append(dict(id=c, label=label, color=color, n=counts[c]))
    reps = [r for r in REPRESENTATIVES if r['tool']=='seurat' and r['branch']=='original_full' and r['population']==pop and int(r['representative_rank']) <= 2 and (pop!='CA1_ProS' or r['cluster_id']!='7')]
    reps.sort(key=lambda r: (int(r['cluster_id']), int(r['representative_rank'])))
    top2 = list(dict.fromkeys(r['gene_id'] for r in reps))
    presets = [dict(id='screenshot', label=f'All {len(ids)} screenshot markers', genes=ids),
               dict(id='top2', label='Saved top-two representatives', genes=top2)]
    mapping = [dict(screenshot_position=i+1, gene_symbol=s, gene_id=by_symbol[s], mapping='exact_unique_source_symbol') for i,s in enumerate(symbols)]
    config = dict(schema='satb2-marker-explorer-3', id=slug, population=pop,
        title='CA1 / ProS' if pop=='CA1_ProS' else 'Inhibitory neurons', provenance_path=output+'/provenance.json',
        sources=sources, catalog=sorted([[g,s] for g,s in genes.items()], key=lambda a:a[1].lower()),
        clusters=clusters, presets=presets, default_preset='screenshot',
        default_clusters=[c['id'] for c in clusters if pop!='CA1_ProS' or c['id']!='7'],
        default_row_label='id' if pop=='CA1_ProS' else 'both',
        screenshot_markers=symbols, n_mice=12, genotypes='6 WT + 6 cKO',
        no_models_refitted=True, reference='Existing pooled full-cluster reference; hiding rows or genes never recomputes it.',
        marker_selection='Exact screenshot panel plus the complete existing cached marker catalog. Selection changes presentation only.',
        screenshot_file='Screenshot 2026-11-27 at 10.03.41.png' if pop=='CA1_ProS' else 'Screenshot 2026-11-27 at 10.03.55.png')
    out = SITE/output
    out.mkdir(parents=True, exist_ok=True)
    (out/'data.js').write_text('window.SATB2_MARKER_EXPLORER_DATA='+json.dumps(config,separators=(',',':'),ensure_ascii=False)+';\n')
    (out/'provenance.json').write_text(json.dumps(config,indent=2,ensure_ascii=False)+'\n')
    with (out/'marker_selection.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(mapping[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(mapping)
    (AUDIT/f'{slug}_screenshot_mapping.json').write_text(json.dumps(mapping,indent=2)+'\n')
    for condition in ['rows','WT_rows','cKO_rows']:
        combined=d[condition]+e[condition]
        keys={(r['cluster_id'],r['gene_id']) for r in combined}
        assert len(keys)==len(combined)
        assert all((c['id'],g) in keys for c in clusters for g in ids)
    print(slug, len(ids), 'screenshot markers;', len(genes), 'cached genes;', len(clusters), 'available clusters')
    return config


if __name__ == '__main__':
    a=build('CA1_ProS','ca1_markers','ca1_sublayers/top_markers',CA1)
    b=build('Inhibitory','inhibitory_markers','inhibitory_markers',IN)
    (AUDIT/'BUILD_AUDIT.json').write_text(json.dumps(dict(seurat_only=True,no_models_refitted=True,
        datasets=[dict(id=x['id'],screenshot_markers=len(x['screenshot_markers']),catalog_genes=len(x['catalog']),clusters=len(x['clusters']),sources=x['sources']) for x in [a,b]]),indent=2)+'\n')
