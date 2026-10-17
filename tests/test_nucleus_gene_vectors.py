"""Focused guards for the exact sparse per-nucleus expression export."""
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

MODULE_PATH = Path(__file__).resolve().parents[1] / 'scripts/build_nucleus_gene_vectors.py'
SPEC = importlib.util.spec_from_file_location('build_nucleus_gene_vectors', MODULE_PATH)
PRODUCER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PRODUCER)


class SparseVectorGuards(unittest.TestCase):
    def test_serialized_counts_and_both_integrity_hashes(self):
        import hashlib
        body = {'indices': [0, 2], 'counts': [2, 100001], 'branch': 'original_full'}
        with tempfile.TemporaryDirectory() as folder:
            desc = PRODUCER.emit(Path(folder), 'check', body)
            packed = (Path(folder) / 'check.json.gz').read_bytes()
            raw = gzip.decompress(packed)
            self.assertEqual(json.loads(raw), body)
            self.assertEqual(desc['sha256'], hashlib.sha256(packed).hexdigest())
            self.assertEqual(desc['payload_sha256'], hashlib.sha256(raw).hexdigest())
            self.assertEqual(desc['gzip_bytes'], len(packed))
            self.assertEqual(desc['uncompressed_bytes'], len(raw))

    def test_actual_zero_and_unknown_nucleus_are_distinct(self):
        axis = {'a': 0, 'b': 1, 'c': 2}
        libraries = [100, 200, 400]
        nonzero = {0: 2, 2: 1}

        def expression(nucleus_id):
            if nucleus_id not in axis:
                return None
            i = axis[nucleus_id]
            return math.log1p(10000 * nonzero.get(i, 0) / libraries[i])

        self.assertAlmostEqual(expression('a'), math.log1p(200))
        self.assertEqual(expression('b'), 0)
        self.assertAlmostEqual(expression('c'), math.log1p(25))
        self.assertIsNone(expression('missing'))

    def test_all_889_source_genes_have_three_integrity_checked_vectors(self):
        root = MODULE_PATH.parents[1]
        path = root / 'site/data/GENE_EXPRESSION_CATALOG.json'
        if not path.exists(): self.skipTest('Full native export has not completed yet')
        catalog = json.loads(path.read_text())
        self.assertEqual(len(catalog['genes']), 889)
        self.assertEqual(len(catalog['shards']), 2670)
        ids = {g['gene_id'] for g in catalog['genes']}
        self.assertEqual(len(ids), 889)
        for branch, item in catalog['branches'].items():
            self.assertEqual(set(item['genes']), ids)
            axis_desc = catalog['shards'][item['axis_shard']]
            axis_file = root / 'gene_expression' / Path(axis_desc['path']).name
            self.assertEqual(PRODUCER.sha_file(axis_file), axis_desc['sha256'])
            axis_raw = gzip.decompress(axis_file.read_bytes())
            self.assertEqual(hashlib.sha256(axis_raw).hexdigest(), axis_desc['payload_sha256'])
            self.assertEqual(len(axis_raw), axis_desc['uncompressed_bytes'])
            axis = json.loads(axis_raw)
            self.assertEqual(axis['branch'], branch)
            self.assertEqual(len(axis['nucleus_ids']), item['n_source_nuclei'])
            self.assertEqual(len(set(axis['nucleus_ids'])), item['n_source_nuclei'])
            self.assertTrue(all(isinstance(n, int) and n > 0 for n in axis['library_counts']))
            for gene, shard_id in item['genes'].items():
                desc = catalog['shards'][shard_id]
                file = root / 'gene_expression' / Path(desc['path']).name
                self.assertEqual(PRODUCER.sha_file(file), desc['sha256'])
                raw = gzip.decompress(file.read_bytes())
                self.assertEqual(hashlib.sha256(raw).hexdigest(), desc['payload_sha256'])
                self.assertEqual(len(raw), desc['uncompressed_bytes'])
                vector = json.loads(raw)
                self.assertEqual(vector['gene_id'], gene)
                self.assertEqual(vector['branch'], branch)
                self.assertEqual(vector['axis_payload_sha256'], axis_desc['payload_sha256'])
                self.assertEqual(vector['n_source_nuclei'], len(axis['nucleus_ids']))
                self.assertEqual(len(vector['indices']), len(vector['counts']))
                self.assertTrue(all(isinstance(c, int) and 0 < c <= 2**53-1 for c in vector['counts']))
                self.assertTrue(all(a < b for a, b in zip(vector['indices'], vector['indices'][1:])))
                self.assertTrue(all(0 <= i < len(axis['nucleus_ids']) for i in vector['indices']))


if __name__ == '__main__': unittest.main()
