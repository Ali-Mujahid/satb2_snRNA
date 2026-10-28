import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

import test_builder as portable_tests

SOURCE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SOURCE))
spec=importlib.util.spec_from_file_location('online_builder',SOURCE/'build_online_reader.py')
online=importlib.util.module_from_spec(spec);spec.loader.exec_module(online)


class OnlineTests(unittest.TestCase):
 def fixture(self,base):
  fixture=portable_tests.ReaderTests();config,_=fixture.fixture(base)
  portable=base/'portable'
  with contextlib.redirect_stdout(io.StringIO()):fixture_module=portable_tests.module.build(config,portable)
  proposal=online.plan_release_assets(portable,'followups-TEST-ONLY-v1')
  rows=[]
  for row in proposal['assets']:
   item={k:v for k,v in row.items() if k not in ('local_path','asset_name')}
   item.update(remote_verified_bytes=row['bytes'],remote_verified_sha256=row['sha256'],verification_method='github_server_digest')
   rows.append(item)
  archive=base/'portable-TEST-ONLY.zip'
  with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
   for p in portable.rglob('*'):
    if p.is_file():z.write(p,p.relative_to(portable).as_posix())
  review=base/'portable_visual_review.md';review.write_text('TEST ONLY reviewed fake fixture, no actual study build')
  remote=base/'release_verification.json';remote.write_text(json.dumps({'test_only':True,'mock_verified_assets':rows}))
  site=base/'existing-site';site.mkdir();(site/'old.html').write_text('Preserved existing site bytes')
  zr=online.receipt(archive)
  cfg={'schema':'satb2-followups-online-transport/v1','approved_for_online_presentation':True,
   'portable_root':str(portable),'portable_manifest':online.receipt(portable/'READER_MANIFEST.json'),
   'portable_visual_review_receipts':[online.receipt(review)],'remote_verification_receipts':[online.receipt(remote)],
   'verified_release_assets':rows,'existing_site_root':str(site),'maximum_online_bytes':1_000_000,'maximum_site_bytes':2_000_000,
   'portable_zip':{'url':'https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/followups-TEST-ONLY-v1/portable.zip',
    'bytes':zr['bytes'],'sha256':zr['sha256'],'local_receipt':zr,'remote_verified_bytes':zr['bytes'],
    'remote_verified_sha256':zr['sha256'],'verification_method':'download_sha256'}}
  p=base/'online-approved.json';p.write_text(json.dumps(cfg));return p,cfg,portable,site

 def test_small_online_transport_preserves_status_values_and_verified_remote_bytes(self):
  with tempfile.TemporaryDirectory(prefix='reader_online_unit_') as t:
   base=Path(t);p,cfg,portable,site=self.fixture(base);before=(site/'old.html').read_bytes()
   out=base/'online';result=online.build(p,out)
   self.assertFalse(result['scientific_values_or_native_outputs_changed']);self.assertEqual((site/'old.html').read_bytes(),before)
   self.assertIn('INPUT_BLOCKED',(out/'index.html').read_text());self.assertIn('Online preview and complete offline download',(out/'index.html').read_text())
   nav=online.read_navigation(out);self.assertEqual(nav['delivery']['mode'],'online')
   for task in nav['tasks']:
    for a in task['assets']:
     if a['extension'] in ('.png','.pdf'):
      self.assertTrue(a['href'].startswith('https://github.com/Ali-Mujahid/satb2_snRNA/releases/download/'))
      self.assertTrue((portable/a['portable_href']).exists());self.assertFalse((out/a['portable_href']).exists())
   self.assertTrue((out/'READER_MANIFEST.json').exists())
   self.assertTrue((out/'NATIVE_ARCHIVE_ONLY_FILES.json').exists())
   self.assertIn('Online transport receipts',(out/'index.html').read_text())
   self.assertLess(sum(x.stat().st_size for x in out.rglob('*') if x.is_file()),cfg['maximum_online_bytes'])

 def test_missing_approval_review_or_mismatched_remote_verification_fails_before_writes(self):
  with tempfile.TemporaryDirectory(prefix='reader_online_unit_') as t:
   base=Path(t);p,cfg,portable,site=self.fixture(base)
   cfg['approved_for_online_presentation']=False;p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):online.build(p,base/'bad')
   cfg['approved_for_online_presentation']=True;cfg['portable_visual_review_receipts']=[];p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):online.build(p,base/'bad')
   cfg['portable_visual_review_receipts']=[online.receipt(base/'portable_visual_review.md')]
   cfg['verified_release_assets'][0]['remote_verified_sha256']='0'*64;p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):online.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())

 def test_unsafe_urls_paths_and_unmanifested_assets_fail(self):
  for value in ['javascript:alert(1)','https://user:secret@github.com/Ali-Mujahid/satb2_snRNA/releases/download/t/a.png','https://example.org/a.png','https://github.com/Someone/else/releases/download/t/a.png']:
   with self.assertRaises(ValueError):online.safe_release_url(value)
  for value in ['../escape.png','%2E%2E/escape.png','https://example.org/a.png','/absolute','folder\\file']:
   with self.assertRaises(ValueError):online.relative_file(value)
  with tempfile.TemporaryDirectory(prefix='reader_online_unit_') as t:
   base=Path(t);p,cfg,portable,site=self.fixture(base)
   cfg['verified_release_assets'][0]['portable_path']='not_a_portable_file.png';p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):online.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())

 def test_incomplete_or_wrong_offline_zip_fails_even_if_zip_transport_matches(self):
  with tempfile.TemporaryDirectory(prefix='reader_online_unit_') as t:
   base=Path(t);p,cfg,portable,site=self.fixture(base);bad=base/'wrong.zip'
   with zipfile.ZipFile(bad,'w') as z:z.writestr('index.html','Wrong/incomplete package')
   rec=online.receipt(bad);cfg['portable_zip'].update(local_receipt=rec,bytes=rec['bytes'],sha256=rec['sha256'],remote_verified_bytes=rec['bytes'],remote_verified_sha256=rec['sha256'])
   p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):online.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())

 def test_site_budget_existing_destination_and_mutated_portable_files_are_rejected(self):
  with tempfile.TemporaryDirectory(prefix='reader_online_unit_') as t:
   base=Path(t);p,cfg,portable,site=self.fixture(base);cfg['maximum_site_bytes']=10;p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):online.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())
   cfg['maximum_site_bytes']=2_000_000;p.write_text(json.dumps(cfg))
   with self.assertRaises(ValueError):online.build(p,site/'bad')
   with self.assertRaisesRegex(ValueError,'protected native/layout'):online.build(p,base/'composition/unsealed_online_child')
   self.assertFalse((base/'composition/unsealed_online_child').exists())
   (base/'bad').mkdir()
   with self.assertRaises(ValueError):online.build(p,base/'bad')
   a=next(a for a in cfg['verified_release_assets'] if a['portable_path'].endswith('.png'))
   (portable/a['portable_path']).write_bytes(b'CHANGED')
   with self.assertRaises(ValueError):online.build(p,base/'another')
   self.assertFalse((base/'another').exists())

 def test_portable_zip_symlink_metadata_is_rejected_even_with_identical_file_bytes(self):
  with tempfile.TemporaryDirectory(prefix='reader_online_unit_') as t:
   base=Path(t);p,cfg,portable,site=self.fixture(base);bad=base/'symlink.zip'
   with zipfile.ZipFile(bad,'w') as z:
    for f in portable.rglob('*'):
     if not f.is_file():continue
     if f.name=='index.html':
      info=zipfile.ZipInfo(f.relative_to(portable).as_posix());info.create_system=3;info.external_attr=(0o120777<<16)
      z.writestr(info,f.read_bytes())
     else:z.write(f,f.relative_to(portable).as_posix())
   rec=online.receipt(bad);cfg['portable_zip'].update(local_receipt=rec,bytes=rec['bytes'],sha256=rec['sha256'],remote_verified_bytes=rec['bytes'],remote_verified_sha256=rec['sha256'])
   p.write_text(json.dumps(cfg))
   with self.assertRaisesRegex(ValueError,'symbolic-link'):online.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())

 def test_declared_zip_prefix_is_supported_but_unsafe_directory_entries_are_rejected(self):
  with tempfile.TemporaryDirectory(prefix='reader_online_unit_') as t:
   base=Path(t);p,cfg,portable,site=self.fixture(base);archive=base/'prefixed.zip'
   with zipfile.ZipFile(archive,'w') as z:
    z.writestr('offline_reader/','')
    for f in portable.rglob('*'):
     if f.is_file():z.write(f,'offline_reader/'+f.relative_to(portable).as_posix())
   rec=online.receipt(archive);cfg['portable_zip'].update(local_receipt=rec,archive_prefix='offline_reader',bytes=rec['bytes'],sha256=rec['sha256'],remote_verified_bytes=rec['bytes'],remote_verified_sha256=rec['sha256']);p.write_text(json.dumps(cfg))
   online.build(p,base/'valid-online')
   with zipfile.ZipFile(archive,'a') as z:z.writestr('../../outside/','')
   rec=online.receipt(archive);cfg['portable_zip'].update(local_receipt=rec,bytes=rec['bytes'],sha256=rec['sha256'],remote_verified_bytes=rec['bytes'],remote_verified_sha256=rec['sha256']);p.write_text(json.dumps(cfg))
   with self.assertRaisesRegex(ValueError,'Unsafe portable relative path'):online.build(p,base/'bad')
   self.assertFalse((base/'bad').exists())

 def test_asset_plan_is_explicitly_unverified_and_has_unique_stable_names(self):
  with tempfile.TemporaryDirectory(prefix='reader_online_unit_') as t:
   base=Path(t);p,cfg,portable,site=self.fixture(base)
   a=online.plan_release_assets(portable,'followups-TEST-ONLY-v1');b=online.plan_release_assets(portable,'followups-TEST-ONLY-v1')
   self.assertEqual(a,b);self.assertIn('NOT_VERIFIED',a['status'])
   self.assertEqual(len(a['assets']),len({r['asset_name'] for r in a['assets']}))
   self.assertTrue(all(r['remote_verified_sha256'] is None for r in a['assets']))


if __name__=='__main__':unittest.main()
