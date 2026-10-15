#!/usr/bin/env python3
"""Index/copy existing rendered pages, and audit produced mouse display inputs.
No science fits, source/count changes or raw nucleus aggregation.
"""
from pathlib import Path
import argparse,collections,html,importlib.util,json,shutil,time
import csv

def main(args):
 study=Path(args.study_root).resolve();saved=study/'reports/saved_results_v2';bundle=study/'reports/offline_dashboard_v2';native=study/'runs/20261014_v2_0_0/results'
 spec=importlib.util.spec_from_file_location('saved_export',study/'supplementary_code/build_saved_results_v2.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);e=m.Evidence()
 original_receipt=e.bind(saved/'COMPLETE.json');prior=json.load(original_receipt.open());cat_path=e.bind(bundle/'data/MOUSE_CATALOG.json');cat=json.load(cat_path.open())
 if prior['status']!='COMPLETE' or cat['status']!='COMPLETE':raise ValueError('Base export/mouse catalog not complete')
 if (saved/'PAGE_FIGURE_INDEX.tsv').exists():raise ValueError('Refusing repeated figure augmentation')
 # Independent all-endpoint factor/hash and full-library count checks.
 producer={b:json.load(e.bind(native/'deseq2'/b/'COMPLETE.json').open()) for b in ['original_full','paired_original','paired_cellbender']};qa=[]
 for fid,f in cat['fits'].items():
  audit=e.bind(f['source_endpoint_audit_path'],f['source_endpoint_audit_sha256']);endpoint=audit.parent.parent
  cov=m.read_tsv(e.bind(endpoint/'AUDIT/MICE_ALL12_AUDIT.tsv'));lookup={c['sample']:c for c in cov};totals=cat['shards'][f['counts_shard']]['all12_integer_column_totals']
  if [int(lookup[s]['total_UMIs']) for s in f['sample_order']]!=totals:raise ValueError('Packed fullgene column totals disagree with native population UMI coverage: '+fid)
  sf_path=f.get('source_size_factors_path');native_sf_sha=None
  if sf_path:
   p=e.bind(sf_path,f['source_size_factors_sha256']);relative=str(p.relative_to(native/'deseq2'/f['branch']));native_sf_sha=producer[f['branch']]['output_sha256'].get(relative)
   if native_sf_sha!=f['source_size_factors_sha256']:raise ValueError('Native immutable producer does not bind the exact source size-factor file')
  qa.append({'fit_id':fid,'count_source_column_totals_match_native_coverage':True,'recorded_size_factor_native_producer_SHA_matches':bool(native_sf_sha) if sf_path else 'NOT_FITTED_NO_FACTORS','status':f['status']})
 m.write_tsv(saved/'MOUSE_COUNT_SIZE_FACTOR_QA.tsv',qa)
 # Give every main book a useful tool/branch/population title.
 books=m.read_tsv(e.bind(native/'plots/PLOT_BOOK_INDEX.tsv'));bookmap={str(Path(x['pdf']).resolve()):x for x in books};figures=m.read_tsv(saved/'FIGURE_INDEX.tsv')
 for r in figures:
  if r['source_pdf'] in bookmap:
   b=bookmap[r['source_pdf']];r['title']=b['tool'].title()+' · '+m.POP_NAMES[b['population']]+' · '+b['branch']+' · '+b['cluster_scope']+' UMAP and marker book ('+b['pages']+' pages)'
 # Distinct native UMAP, representative-marker and supplied/canonical marker pages.
 pdflookup={r['source_pdf']:r['saved_pdf'] for r in figures};pages=[]
 for b in books:
  if b['cluster_scope']!='original_frozen':continue
  pi=Path(b['pdf']).parent/'PAGE_INDEX.tsv';e.bind(pi)
  for p in m.read_tsv(pi):
   name=p['name'];kind='Representative marker dotplot' if name.startswith('representative_markers') else 'Provided-sheet + canonical marker dotplot' if name.startswith('supplied_panel') else 'UMAP'
   outdir=saved/'figures/pages'/b['tool']/b['branch']/b['population'];png=e.clone(p['png'],outdir/Path(p['png']).name);svg=e.clone(p['svg'],outdir/Path(p['svg']).name)
   pages.append({'tool':b['tool'],'branch':b['branch'],'population':b['population'],'population_display_name':m.POP_NAMES[b['population']],'cluster_scope':'original_frozen','figure_kind':kind,'native_page_name':name,'pdf_page':int(p['page']),'saved_pdf':pdflookup[str(Path(b['pdf']).resolve())]+'#page='+p['page'],'saved_png':str(png.relative_to(saved)),'saved_svg':str(svg.relative_to(saved)),'source_png':p['png'],'source_svg':p['svg'],'png_sha256':m.sha(png),'svg_sha256':m.sha(svg),'displayed_nuclei':b['displayed_nuclei'],'physical_layer_assignment_validated':False})
 # Existing explicitly program-colored UMAPs: all tools/branches, CA1-ProS and IN.
 program_dir=study/'reports/visual_guide_v4/programs'
 for tool in ['scanpy','seurat']:
  for branch in ['original_full','paired_original','paired_cellbender']:
   for pop in ['CA1_ProS','Inhibitory']:
    stem=tool+'__'+branch+'__'+pop+'__expression_program';targets={}
    for ext in ['pdf','png','svg']:
     source=program_dir/(stem+'.'+ext)
     if not source.exists():raise ValueError('Completed program-render source absent: '+str(source))
     target=e.clone(source,saved/'figures/expression_programs'/source.name);targets[ext]=str(target.relative_to(saved))
    label='CA1 deep/superficial-associated expression programs' if pop=='CA1_ProS' else 'Inhibitory parent-associated expression programs'
    pages.append({'tool':tool,'branch':branch,'population':pop,'population_display_name':m.POP_NAMES[pop],'cluster_scope':'original_frozen','figure_kind':'Expression-program UMAP','native_page_name':label,'pdf_page':1,'saved_pdf':targets['pdf'],'saved_png':targets['png'],'saved_svg':targets['svg'],'source_png':str(program_dir/(stem+'.png')),'source_svg':str(program_dir/(stem+'.svg')),'png_sha256':m.sha(program_dir/(stem+'.png')),'svg_sha256':m.sha(program_dir/(stem+'.svg')),'displayed_nuclei':'SEE_SOURCE_PROGRAM_METADATA','physical_layer_assignment_validated':False})
 # Existing completed effect SVG/PNG output, alongside already saved PDFs.
 for f in m.read_tsv(study/'reports/effect_figures/FIGURE_INDEX.tsv'):
  for ext in ['png','svg']:
   src=f[ext+'_path'];e.clone(src,saved/'figures/effects'/Path(src).name)
 m.write_tsv(saved/'PAGE_FIGURE_INDEX.tsv',pages);m.write_tsv(saved/'FIGURE_INDEX.tsv',figures)
 # Rebuild own HTML with named books, then add a separate discoverable page table.
 primary=m.read_tsv(saved/'PRIMARY_RESULTS_INDEX.tsv');sensitivity=m.read_tsv(saved/'SENSITIVITY_RESULTS_INDEX.tsv');m.build_html(saved,primary,sensitivity,figures)
 rows=[]
 for p in pages:
  esc=html.escape
  key=' '.join(str(p[k]) for k in ['tool','branch','population_display_name','figure_kind','native_page_name'])
  links=' · '.join('<a href="'+esc(p['saved_'+ext])+'">'+ext.upper()+'</a>' for ext in ['pdf','png','svg'])
  rows.append('<tr data-figure-search="'+esc(key.lower())+'"><td>'+esc(p['tool'])+'</td><td>'+esc(p['branch'])+'</td><td>'+esc(p['population_display_name'])+'</td><td>'+esc(p['figure_kind'])+'</td><td>'+esc(p['native_page_name'])+'</td><td>'+links+'</td></tr>')
 section='<section id="figure-pages"><h2>UMAP and marker pages: direct files</h2><p>Every original-frozen dataset has separate UMAP, representative-marker, and provided-sheet/canonical-panel pages below. CA1 and inhibitory program-colored UMAPs are separate entries. PDF links open the exact book page where supported; PNG and SVG are standalone saved pages.</p><p>The reviewed marker panel includes123 mapped workbook rows (117 unique stable genes), with canonical additions giving145 mapped rows (139 unique genes). Two workbook AC129186.1 rows were unmapped. Program names describe expression associations rather than validated physical sublayers.</p><p><a href="PAGE_FIGURE_INDEX.tsv">Full page figure index</a> · <input id="figure-search" placeholder="Filter figures: original_full CA1 representative, supplied, superficial…"></p><div class="tablebox"><table><thead><tr><th>Tool</th><th>Count branch</th><th>Population</th><th>Figure type</th><th>Page / program</th><th>Saved formats</th></tr></thead><tbody>'+''.join(rows)+'</tbody></table></div></section>'
 hp=saved/'index.html';text=hp.read_text();text=text.replace('<h2>Per-total, per-cluster and reviewed category DEG lists</h2>',section+'<h2>Per-total, per-cluster and reviewed category DEG lists</h2>');text=text.replace("const q=document.querySelector('#search'),rs=[...document.querySelectorAll('tbody tr')];", "const q=document.querySelector('#search'),rs=[...document.querySelectorAll('tr[data-search]')];")
 text=text.replace('</body>',"<script>const fq=document.querySelector('#figure-search'),frs=[...document.querySelectorAll('tr[data-figure-search]')];fq.addEventListener('input',()=>{const t=fq.value.toLowerCase();for(const r of frs)r.hidden=!r.dataset.figureSearch.includes(t)});</script></body>")
 hp.write_text(text,encoding='utf-8')
 # Keep the preliminary engineering receipt, and bind final navigation additions.
 history=study/'reports/saved_results_v2_build_history'/('pre_figure_page_augmentation_'+str(int(time.time())));history.mkdir(parents=True);shutil.copy2(original_receipt,history/'COMPLETE.json')
 e.inputs[str(history/'COMPLETE.json')]=e.inputs.pop(str(original_receipt));exports=m.read_tsv(saved/'SOURCE_FILE_EXPORT_INDEX.tsv');exports+=e.exports;m.write_tsv(saved/'SOURCE_FILE_EXPORT_INDEX.tsv',exports)
 prior['source_figure_augmentation_code_sha256']=m.sha(__file__);prior['input_sha256'].update(e.inputs);prior['n_standalone_figure_pages']=len(pages);prior['n_main_book_datasets']=sum(b['cluster_scope']=='original_frozen' for b in books);prior['standalone_figure_kinds']=dict(collections.Counter(p['figure_kind'] for p in pages));prior['all1060_native_size_factor_SHA_and_fullgene_column_coverage_QA']='PASS';prior['pre_navigation_receipt_archive']=str(history/'COMPLETE.json');prior['source_isolated_copy_methods']=dict(collections.Counter(x['copy_method'] for x in exports));prior['output_sha256']={str(p.relative_to(saved)):m.sha(p) for p in saved.rglob('*') if p.is_file() and p.name!='COMPLETE.json'};m.dump(saved/'COMPLETE.json',prior)
 m.dump(saved/'PAGE_FIGURE_AUGMENTATION_COMPLETE.json',{'status':'COMPLETE','source_code_sha256':m.sha(__file__),'n_standalone_figure_pages':len(pages),'all1060_model_source_QA_PASS':True,'output_sha256':{str(p.relative_to(saved)):m.sha(p) for p in [saved/'PAGE_FIGURE_INDEX.tsv',saved/'FIGURE_INDEX.tsv',saved/'MOUSE_COUNT_SIZE_FACTOR_QA.tsv',saved/'index.html']},'input_sha256':e.inputs,'native_outputs_modified':False,'new_fitting_performed':False})
 print(json.dumps({k:prior[k] for k in ['n_standalone_figure_pages','n_main_book_datasets','standalone_figure_kinds','all1060_native_size_factor_SHA_and_fullgene_column_coverage_QA']},indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--study-root',required=True);main(ap.parse_args())
