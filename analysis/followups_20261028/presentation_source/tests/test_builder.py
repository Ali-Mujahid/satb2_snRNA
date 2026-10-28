import contextlib, gzip, hashlib, importlib.util, io, json, pathlib, re, tempfile, unittest
from unittest import mock

SOURCE=pathlib.Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('reader_builder',SOURCE/'build_reader.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def rec(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

class ReaderTests(unittest.TestCase):
 def add_layout(self,base,cfg):
  native=pathlib.Path(cfg['tasks']['cellchat']['native_root']);root=base/'layout';root.mkdir()
  scientific=json.loads((native/'FINAL_VALIDATION.json').read_text());mapping=[];outputs=[]
  for record in scientific['files']:
   original=native/record['path']
   if original.suffix not in ('.png','.pdf'):continue
   corrected=root/record['path'];corrected.parent.mkdir(parents=True,exist_ok=True);corrected.write_bytes(b'READABLE_LAYOUT_TEST_ONLY_'+original.read_bytes())
   mapping.append({'native_original':rec(original),'corrected_presentation':rec(corrected)})
   outputs.append(rec(corrected))
  mapping_path=root/'NATIVE_TO_PRESENTATION_PLOT_MAPPING.json';mapping_path.write_text(json.dumps(mapping));outputs.append(rec(mapping_path))
  seal=root/'FINAL_VALIDATION.json';seal.write_text(json.dumps({'schema':'satb2-cellchat-layout-only-final/v1','status':'COMPLETE_PENDING_VISUAL_REVIEW',
   'scientific_fit_or_rescore':False,'scientific_aggregation_formulas_changed':False,'colors_and_numerical_scales_changed':False,
   'original_native_outputs_modified':False,'new_statistical_results':False,'native_final_validation':cfg['tasks']['cellchat']['terminal_seal'],'outputs':outputs}))
  review=base/'layout_visual_review.md';review.write_text('TEST ONLY: independent artificial layout review, no biological study result')
  cfg['tasks']['cellchat']['presentation_override']={'approved_for_display':True,'presentation_root':str(root),'terminal_seal':rec(seal),'visual_review_receipts':[rec(review)]}
  return mapping_path,seal

 def index(self,out):
  return json.loads(re.search(r'const INDEX=(.*?);\nconst el=',(out/'index.html').read_text()).group(1))

 def fixture(self,base):
  tasks={}
  for kind,(schema,states,key) in module.SEALS.items():
   root=base/kind;root.mkdir()
   relative={'composition':'fits/seurat__inhibitory_clusters/all12/raw_mouse_proportions.png','cellchat':'figures/triMean/CA1_to_IN_dotplot_page1.png','monocle3':'monocle3/CA1/models/CA1_DL_SL/figures/unrooted_condition.png'}[kind]
   png=root/relative;png.parent.mkdir(parents=True);png.write_bytes(b'TEST_ONLY_PNG_BYTES')
   pdf=png.with_suffix('.pdf');pdf.write_bytes(b'TEST_ONLY_PDF_BYTES')
   table=png.parent/'saved_native_table.tsv';table.write_text('gene\tvalue\nGeneA\tNA\nGeneB\t0.003\n')
   model=root/'large_native_model.rds';model.write_bytes(b'TEST_ONLY_MODEL_BYTES')
   rows=[]
   for path in [png,pdf,table,model]:
    r=rec(path);r['path']=str(path.relative_to(root));rows.append(r)
   state='COMPLETE_WITH_INCONCLUSIVE_OR_UNTESTABLE_SCOPES' if kind=='composition' else 'COMPLETE'
   seal=root/'FINAL_VALIDATION.json';seal.write_text(json.dumps({'schema':schema,'status':state,key:rows}))
   review=base/(kind+'_review.md');review.write_text('Test fixture review; not a real study result.\n')
   tasks[kind]={'native_root':str(root),'terminal_seal':rec(seal),'headline':'TEST ONLY','findings':['Test-only saved result fixture'],'limitations':['No biological interpretation'],'review_receipts':[rec(review)]}
  assessment=base/'blocked.json';assessment.write_text(json.dumps({'schema':'satb2-scvelo-starting-point-assessment/v1','task_status':'INPUT_BLOCKED_NOT_EXECUTED','scientific_jobs_started':False}))
  review=base/'velocity_review.md';review.write_text('Test only blocked input review')
  tasks['velocity']={'status':'INPUT_BLOCKED','assessment_receipt':rec(assessment),'headline':'Input blocked','findings':['No analysis run'],'limitations':['Missing measured velocity layers'],'review_receipts':[rec(review)]}
  cfg={'schema':'satb2-followups-portable-presentation/v1','approved_for_portable_presentation':True,'maximum_package_bytes':1_000_000,'tasks':tasks}
  p=base/'approved.json';p.write_text(json.dumps(cfg));return p,cfg

 def test_native_transport_all_four_sections_and_local_pdf_links(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);out=base/'portable'
   with contextlib.redirect_stdout(io.StringIO()):result=module.build(p,out)
   self.assertTrue(result['no_fit_rescore_or_native_mutation']);self.assertEqual(len(result['native_seals']),3)
   content=(out/'index.html').read_text();self.assertIn('INPUT_BLOCKED',content);self.assertIn('Original PDF',content)
   self.assertNotIn('/*__REVIEWED_NATIVE_NAVIGATION__*/',content)
   original=base/'composition/fits/seurat__inhibitory_clusters/all12/raw_mouse_proportions.pdf'
   copied=out/'native/composition/fits/seurat__inhibitory_clusters/all12/raw_mouse_proportions.pdf'
   self.assertEqual(original.read_bytes(),copied.read_bytes())
   self.assertIn('GeneA',content);self.assertIn('COMPLETE_WITH_INCONCLUSIVE',content)
   self.assertEqual(len(json.loads((out/'NATIVE_ARCHIVE_ONLY_FILES.json').read_text())),3)
   task=next(t for t in self.index(out)['tasks'] if t['id']=='cellchat')
   self.assertNotIn('default_rendering',task)
   self.assertIn('vertical jitter',task['original_layout_warning'])
   self.assertIn('does not equal the saved raw score',task['original_layout_warning'])
   self.assertEqual({a['rendering'] for a in task['assets'] if a['extension']=='.png'},{'Archival native (known display issues)'})

 def test_verified_layout_default_preserves_original_and_relative_corrected_pdf(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);self.add_layout(base,cfg);p.write_text(json.dumps(cfg));out=base/'reader'
   with contextlib.redirect_stdout(io.StringIO()):module.build(p,out)
   task=next(t for t in self.index(out)['tasks'] if t['id']=='cellchat')
   self.assertEqual(task['default_rendering'],'Readable layout')
   self.assertEqual({a['rendering'] for a in task['assets'] if a['extension']=='.png'},{'Archival native (known display issues)','Readable layout'})
   pdf=next(a for a in task['assets'] if a['extension']=='.pdf' and a['rendering']=='Readable layout')
   self.assertTrue(pdf['href'].startswith('presentation/cellchat/figures/triMean/'))
   self.assertFalse(pathlib.Path(pdf['href']).is_absolute());self.assertTrue((out/pdf['href']).exists())
   self.assertTrue((out/pdf['original_href']).exists())
   self.assertNotEqual((out/pdf['href']).read_bytes(),(out/pdf['original_href']).read_bytes())
   self.assertTrue((out/task['layout_provenance']['mapping_href']).exists())
   self.assertIn('Archival PDF (known display issues; preserved)',(out/'index.html').read_text())

 def test_layout_missing_visual_review_or_mapping_mismatch_fails_closed(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);mapping_path,seal=self.add_layout(base,cfg)
   saved=cfg['tasks']['cellchat']['presentation_override']['visual_review_receipts']
   cfg['tasks']['cellchat']['presentation_override']['visual_review_receipts']=[];p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):module.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())
   cfg['tasks']['cellchat']['presentation_override']['visual_review_receipts']=saved
   mapping=json.loads(mapping_path.read_text());mapping[0]['native_original']['sha256']='0'*64;mapping_path.write_text(json.dumps(mapping))
   data=json.loads(seal.read_text());data['outputs'][-1]=rec(mapping_path);seal.write_text(json.dumps(data));cfg['tasks']['cellchat']['presentation_override']['terminal_seal']=rec(seal);p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):module.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())

 def test_unapproved_or_pending_seal_fails_before_destination_created(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);cfg['approved_for_portable_presentation']=False;p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):module.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())
   cfg['approved_for_portable_presentation']=True;seal=pathlib.Path(cfg['tasks']['cellchat']['terminal_seal']['path']);x=json.loads(seal.read_text());x['status']='RUNNING';seal.write_text(json.dumps(x));cfg['tasks']['cellchat']['terminal_seal']=rec(seal);p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):module.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())

 def test_modified_native_asset_and_escaping_receipt_rejected(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);native=base/'cellchat/figures/triMean/CA1_to_IN_dotplot_page1.png';native.write_bytes(b'CORRUPTED')
   with self.assertRaises(ValueError):module.build(p,base/'bad')
   external=base/'outside.tsv';external.write_text('TEST_ONLY')
   seal=base/'composition/FINAL_VALIDATION.json';x=json.loads(seal.read_text());x['outputs'][0]=rec(external);seal.write_text(json.dumps(x));cfg['tasks']['composition']['terminal_seal']=rec(seal);p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):module.build(p,base/'bad')

 def test_velocity_cannot_be_silently_unblocked(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);cfg['tasks']['velocity']['status']='COMPLETE';p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):module.build(p,base/'bad')

 def test_archive_javascript_or_credential_urls_are_rejected(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base)
   for url in ['javascript:alert(1)','https://user:password@example.org/archive.zip']:
    cfg['tasks']['composition']['archive_url']=url;p.write_text(json.dumps(cfg))
    with self.assertRaises(ValueError):module.build(p,base/'bad')

 def test_budget_and_existing_destination_no_silent_drops(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);cfg['maximum_package_bytes']=10;p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):module.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())
   cfg['maximum_package_bytes']=1_000_000;p.write_text(json.dumps(cfg));(base/'bad').mkdir()
   with self.assertRaises(ValueError):module.build(p,base/'bad')

 def test_synthetic_study_roots_fail_and_script_content_is_escaped(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base)
   cfg['tasks']['cellchat']['headline']='</script><script>alert(1)</script>';p.write_text(json.dumps(cfg))
   with contextlib.redirect_stdout(io.StringIO()):module.build(p,base/'escaped')
   self.assertIn('\\u003c/script>',(base/'escaped/index.html').read_text())
   path=base/'synthetic_study';(base/'cellchat').rename(path);cfg['tasks']['cellchat']['native_root']=str(path);cfg['tasks']['cellchat']['terminal_seal']=rec(path/'FINAL_VALIDATION.json');p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):module.build(p,base/'bad')

 def test_complete_large_gene_table_gzip_is_lossless_deterministic_and_bounded(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base)
   root=pathlib.Path(cfg['tasks']['monocle3']['native_root'])
   source=root/'monocle3/CA1/models/CA1_DL_SL/graph_test_all.tsv'
   original=('gene_symbol\tq_value\tstatus\n'+'GeneA\t0\tOK\n'+'GeneB\tNA\tFAIL\n'*100).encode()
   source.write_bytes(original)
   axes=root/'input_cells.tsv';axes.write_bytes(b'nucleus_id\n'+b'Mouse_Nucleus\n'*100)
   seal=root/'FINAL_VALIDATION.json';data=json.loads(seal.read_text())
   for path in (source,axes):
    row=rec(path);row['path']=path.relative_to(root).as_posix();data['outputs'].append(row)
   seal.write_text(json.dumps(data));cfg['tasks']['monocle3']['terminal_seal']=rec(seal)
   cfg['maximum_table_bytes']=64;p.write_text(json.dumps(cfg));out=base/'portable'
   with contextlib.redirect_stdout(io.StringIO()):result=module.build(p,out)
   encoded=result['lossless_compressed_native_tables'];self.assertEqual(len(encoded),1)
   record=encoded[0];target=out/record['transport']['path']
   self.assertEqual(gzip.decompress(target.read_bytes()),original)
   self.assertEqual(record['decoded_sha256'],hashlib.sha256(original).hexdigest())
   self.assertEqual(source.read_bytes(),original)
   second=base/'again.tsv.gz';module.lossless_gzip(source,second)
   self.assertEqual(target.read_bytes(),second.read_bytes())
   task=next(t for t in self.index(out)['tasks'] if t['id']=='monocle3')
   asset=next(a for a in task['assets'] if a['extension']=='.gz')
   self.assertTrue((out/asset['href']).is_file());self.assertEqual(asset['original_native_receipt']['sha256'],record['decoded_sha256'])
   self.assertLess(target.stat().st_size,source.stat().st_size)
   self.assertTrue(any(a['native_path']=='input_cells.tsv' for a in json.loads((out/'NATIVE_ARCHIVE_ONLY_FILES.json').read_text())))
   self.assertFalse((out/'native/monocle3/input_cells.tsv').exists())

 def test_fresh_destination_inside_scientific_or_layout_run_is_rejected_before_writes(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);self.add_layout(base,cfg)
   run=base/'existing_run';run.mkdir();results=run/'results';(base/'composition').rename(results)
   cfg['tasks']['composition']['native_root']=str(results);cfg['tasks']['composition']['terminal_seal']=rec(results/'FINAL_VALIDATION.json');p.write_text(json.dumps(cfg))
   targets=[results/'unsealed_child',run/'workflow/new_child',base/'layout/unsealed_child']
   for target in targets:
    with self.assertRaisesRegex(ValueError,'protected native/layout'):module.build(p,target)
    self.assertFalse(target.exists())

 def test_compression_scratch_cannot_be_created_inside_a_native_tree(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);root=base/'composition'
   with mock.patch.object(module.tempfile,'gettempdir',return_value=str(root)):
    with self.assertRaisesRegex(ValueError,'protected native/layout'):module.build(p,base/'reader')
   self.assertFalse((base/'reader').exists());self.assertFalse(list(root.glob('satb2_reader_transport_*')))

 def test_cellchat_boxplot_direction_comes_from_the_sealed_plot_manifest(self):
  with tempfile.TemporaryDirectory(prefix='reader_unit_') as t:
   base=pathlib.Path(t);p,cfg=self.fixture(base);root=base/'cellchat';folder=root/'figures/triMean'
   png=folder/'pathway_mouse_boxplots.png';png.write_bytes(b'TEST_ONLY_BOX_PNG');pdf=png.with_suffix('.pdf');pdf.write_bytes(b'TEST_ONLY_BOX_PDF')
   plot_manifest=folder/'PLOT_MANIFEST.tsv';plot_manifest.write_text('base\ttype\tdirection\tcondition\npathway_mouse_boxplots\tboxplot\tbidirectional\tboth\n')
   seal=root/'FINAL_VALIDATION.json';data=json.loads(seal.read_text())
   for path in (png,pdf,plot_manifest):
    row=rec(path);row['path']=path.relative_to(root).as_posix();data['files'].append(row)
   seal.write_text(json.dumps(data));cfg['tasks']['cellchat']['terminal_seal']=rec(seal);p.write_text(json.dumps(cfg))
   with contextlib.redirect_stdout(io.StringIO()):module.build(p,base/'reader')
   task=next(t for t in self.index(base/'reader')['tasks'] if t['id']=='cellchat')
   boxplots=[a for a in task['assets'] if a['relative_native_path'].endswith(('pathway_mouse_boxplots.png','pathway_mouse_boxplots.pdf'))]
   self.assertEqual(len(boxplots),2);self.assertEqual({a['direction'] for a in boxplots},{'bidirectional'})

if __name__=='__main__':unittest.main()
