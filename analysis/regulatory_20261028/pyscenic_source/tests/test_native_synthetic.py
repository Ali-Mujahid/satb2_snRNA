#!/usr/bin/env python3
"""Synthetic-only native API/IO tests; never loads SATB2 expression or references."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from ctxcore.genesig import Regulon
from ctxcore.rnkdb import RankingDatabase

SOURCE = Path(__file__).resolve().parents[1] / "scripts/run_pyscenic_native.py"
spec = importlib.util.spec_from_file_location("scenic_producer", SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ToyRankingDatabase(RankingDatabase):
    """A synthetic motif ranking object, intentionally unrelated to real refs."""
    def __init__(self):
        super().__init__("SYNTHETIC_ONLY_TOY_MOTIFS")
        genes = ["G" + str(i) for i in range(2000)]
        rng = np.random.RandomState(71)
        rankings = np.vstack([np.arange(2000)] + [rng.permutation(2000) for _ in range(99)])
        self.df = pd.DataFrame(rankings.astype(np.int16), index=["motif" + str(i) for i in range(100)], columns=genes)

    @property
    def total_genes(self):
        return 2000

    @property
    def genes(self):
        return tuple(self.df.columns)

    def load_full(self):
        return self.df

    def load(self, signature):
        return self.df[[g for g in self.df.columns if g in signature.genes]]


class NativeSynthetic(unittest.TestCase):
    def test_unambiguous_genotype_blind_full_library(self):
        counts = sparse.csc_matrix(np.array([[1, 0, 2, 1], [0, 4, 0, 1], [1, 1, 1, 1], [1, 1, 0, 3], [1, 1, 0, 2]], float))
        genes = pd.DataFrame({"gene_symbol": ["TF", "A", "B", "DUP", "DUP"], "duplicate_symbol": [False, False, False, True, True], "symbol_missing": False})
        cells = pd.DataFrame({"endpoint": ["E1", "E1", "E2", "E2"], "condition": ["WT", "cKO", "WT", "cKO"], "library_umi_full_genes": np.asarray(counts.sum(axis=0)).ravel()})
        x, names, audit = m.prepare_matrix(counts, genes, cells)
        x2, names2, _ = m.prepare_matrix(counts, genes, cells.assign(condition=["changed"] * 4))
        np.testing.assert_array_equal(x.toarray(), x2.toarray())
        self.assertEqual(names, names2)
        self.assertNotIn("DUP", names)
        np.testing.assert_allclose(x.toarray()[0, names.index("TF")], np.log1p(10000 / 4))
        self.assertFalse(audit.iloc[3].eligible_full_RNA_background)

    def test_sparse_native_grnboost2_seed(self):
        x = sparse.csc_matrix(np.random.RandomState(21).poisson(2, (40, 24)).astype(float))
        genes = ["G" + str(i) for i in range(24)]
        a = m.native_grn(x, genes, ["G0", "G1", "G2"], 1337, 2)
        b = m.native_grn(x, genes, ["G0", "G1", "G2"], 1337, 2)
        pd.testing.assert_frame_equal(a, b, check_exact=True)
        self.assertEqual(set(a.target), set(genes))

    def test_native_positive_modules(self):
        from pyscenic.utils import add_correlation, modules_from_adjacencies
        rng = np.random.RandomState(31)
        latent = np.arange(100, dtype=float)
        df = pd.DataFrame({"TF": latent, **{"G" + str(i): latent + rng.normal(0, .1, 100) for i in range(30)}, "NEG": -latent})
        adj = pd.DataFrame({"TF": ["TF"] * 31, "target": ["G" + str(i) for i in range(30)] + ["NEG"], "importance": np.arange(1, 32, dtype=float)})
        corr = add_correlation(adj, df, rho_threshold=.03, mask_dropouts=False)
        regs = modules_from_adjacencies(corr, df, min_genes=20, rho_dichotomize=True, keep_only_activating=True)
        self.assertTrue(regs)
        self.assertTrue(all("activating" in r.context and "NEG" not in r.genes for r in regs))
        self.assertTrue(all("TF" in r.genes for r in regs))

    def test_native_ctx_pruning_and_regulon_conversion(self):
        from pyscenic.prune import prune2df
        from pyscenic.transform import df2regulons
        module = Regulon(name="G0_module", transcription_factor="G0", gene2weight={"G" + str(i): 1.0 for i in range(30)}, gene2occurrence=[], context=frozenset(["activating"]))
        with tempfile.TemporaryDirectory() as tmp:
            ann = Path(tmp) / "SYNTHETIC_motif_annotations.tsv"
            pd.DataFrame({"#motif_id": ["motif0"], "gene_name": ["G0"], "motif_similarity_qvalue": [0.0], "orthologous_identity": [1.0], "description": ["SYNTHETIC_ONLY_DIRECT_ANNOTATION"]}).to_csv(ann, sep="\t", index=False)
            cluster, client = m.local_client(2)
            try:
                motifs = prune2df([ToyRankingDatabase()], [module], str(ann), rank_threshold=1500, auc_threshold=.05, nes_threshold=3.0,
                                  motif_similarity_fdr=.001, orthologuous_identity_threshold=0.0, client_or_address=client, num_workers=2, module_chunksize=1)
            finally:
                client.close()
                cluster.close()
            self.assertFalse(motifs.empty)
            regs = df2regulons(motifs)
            self.assertEqual(regs[0].name, "G0(+)")
            plain = m.motif_table_plain(motifs)
            self.assertIn("TargetGenes", plain.columns)
            p = Path(tmp) / "native.pkl.gz"
            m.pickle_write(p, regs)
            self.assertEqual(m.pickle_read(p)[0].gene2weight, regs[0].gene2weight)

    def test_batched_native_auc_exact_with_zeros_and_ties(self):
        rng = np.random.RandomState(13)
        values = rng.poisson(.5, (37, 100)).astype(float)
        values[:, 1] = values[:, 2]
        values[0, :] = 0
        genes = ["G" + str(i) for i in range(100)]
        ids = ["SYNTHETIC_NUCLEUS_" + str(i) for i in range(37)]
        regs = [Regulon(name="G0(+)", transcription_factor="G0", gene2weight={"G" + str(i): i + 1 for i in range(15)}, gene2occurrence=[], context=frozenset(["activating"])),
                Regulon(name="G20(+)", transcription_factor="G20", gene2weight={"G" + str(i): 1 for i in range(20, 40)}, gene2occurrence=[], context=frozenset(["activating"]))]
        a = m.native_auc_batched(sparse.csc_matrix(values), genes, ids, regs, batch_size=1000)
        b = m.native_auc_batched(sparse.csc_matrix(values), genes, ids, regs, batch_size=7)
        pd.testing.assert_frame_equal(a, b, check_exact=True)
        self.assertTrue((a.max() > a.min()).all())
        self.assertTrue(np.isfinite(a.to_numpy()).all())
        without_self = regs[0].copy(gene2weight={g: w for g, w in regs[0].gene2weight.items() if g != "G0"})
        self.assertEqual(len(without_self), len(regs[0]) - 1)
        self.assertIn("G0", regs[0].genes)

    def test_output_and_receipt_guards(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = m.fresh_dir(root / "out")
            m.table_write(pd.DataFrame({"x": [1]}), out / "x.tsv")
            r = m.receipt(out / "x.tsv")
            self.assertEqual(r, m.verify(r))
            with self.assertRaises(ValueError):
                m.fresh_dir(out)
            with self.assertRaises(ValueError):
                m.table_write(pd.DataFrame({"x": [2]}), out / "x.tsv")
            (out / "x.tsv").write_text("changed")
            with self.assertRaises(ValueError):
                m.verify(r)
            (root / "link").symlink_to(out / "x.tsv")
            with self.assertRaises(ValueError):
                m.receipt(root / "link")


if __name__ == "__main__":
    unittest.main(verbosity=2)
