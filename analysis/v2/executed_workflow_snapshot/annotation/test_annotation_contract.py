#!/usr/bin/env python3
"""Synthetic contract checks only: no user counts, fits, or biological results."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

import annotate_clusters as a


def fixture_metadata(tool, branch, population):
    ids = range(420) if branch == "original_full" else range(4, 420)
    rows = []
    for i in ids:
        sample = "m" + str(i % 12)
        cluster = ("0" if i < 100 else "1" if i < 220 else "2") if tool == "scanpy" else ("1" if i < 99 else "2" if i < 220 else "3")
        broad = "Inhibitory" if population == "Inhibitory" else "Excitatory"
        rows.append({"nucleus_id": f"{sample}_b{i}", "sample": sample, "condition": "WT" if i % 12 < 6 else "cKO",
                     "current_cluster_id": cluster, "source_cluster_id": cluster,
                     "reference_broad_prediction": broad, "reference_subclass_name": "016 CA1-ProS Glut" if broad == "Excitatory" else "053 Sst Gaba",
                     "reference_supertype_name": "toy_supertype", "reference_broad_high_support": True,
                     "reference_subclass_high_support": True})
    return pd.DataFrame(rows)


class AnnotationContract(unittest.TestCase):
    def test_iterative_shared_mask_and_rare_population_retention(self):
        ds = {(t, b, p): fixture_metadata(t, b, p) for t in a.TOOLS for b in a.BRANCHES for p in a.POPULATIONS}
        masks, removed, iterations = a.shared_cluster_masks(ds, {"cluster_filter_mode": "shared_iterative_intersection"})
        self.assertEqual(len(masks[("original_full", "allcells")]), 320)
        self.assertEqual(iterations.query('branch=="original_full" and population=="allcells"').n_removed.tolist(), [99, 1])
        self.assertEqual(masks[("paired_original", "allcells")], masks[("paired_cellbender", "allcells")])
        self.assertEqual(len(masks[("original_full", "Inhibitory")]), 420)
        self.assertEqual(len(masks[("original_full", "CA1_ProS")]), 420)
        self.assertTrue(removed.population.isin(["allcells", "Excitatory"]).all())

    def test_Satb2_excluded_from_layer_program(self):
        panel = pd.DataFrame({"cluster_id": ["0", "1", "0", "1"], "gene_id": ["g1", "g1", "g2", "g2"],
                              "gene_symbol": ["Satb2", "Satb2", "Asic2", "Asic2"], "fraction_detected": [1., 0., .8, .2]})
        defs = pd.DataFrame({"set": [a.LAYER_SETS[0], a.LAYER_SETS[0]], "gene_symbol": ["Satb2", "Asic2"]})
        scores = a.program_scores(panel, defs)
        self.assertTrue(scores.mapped_program_genes.eq(1).all())
        self.assertEqual(scores.set_index("cluster_id").loc["0", "program_score"], 1.)

    def test_source_change_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "input.tsv"
            p.write_text("x\n1\n")
            reader = a.Inputs()
            reader.table(p)
            p.write_text("x\n2\n")
            with self.assertRaises(ValueError): reader.validate()

    def test_initial_to_final_category_union_and_paired_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            previous = root / "previous"
            shared = previous / "results/shared/original_full"
            shared.mkdir(parents=True)
            inputs_dir = previous / "results/inputs"
            inputs_dir.mkdir()
            names = sorted(set(sum(a.COARSE.values(), []) + sum(a.FINE_IN.values(), []) +
                               ["Asic2", "Epha3", "Fstl4", "Ndst4", "Col11a1", "Fras1", "Satb2"]))
            genes = pd.DataFrame({"gene_id": [f"toy_g{i}" for i in range(len(names))], "gene_symbol": names})
            genes.to_csv(shared / "gene_features.tsv", sep="\t", index=False)
            gene_lookup = genes.set_index("gene_symbol").gene_id.to_dict()
            rows = []
            for group, symbols in {**a.COARSE, a.LAYER_SETS[0]: ["Asic2", "Epha3", "Fstl4"],
                                   a.LAYER_SETS[1]: ["Ndst4", "Col11a1", "Fras1", "Satb2"]}.items():
                for symbol in symbols: rows.append({"group": group, "set": group, "gene_symbol": symbol,
                                                  "source": "synthetic", "selection": "synthetic", "eligible": True})
            supplied = pd.DataFrame(rows)
            supplied.to_csv(inputs_dir / "marker_review_panels.tsv", sep="\t", index=False)
            cfg = {"previous_run_root": str(previous), "initial_annotation_dir": str(root / "initial"),
                   "marker_summary_dir": str(root / "summaries"), "cluster_filter_mode": "shared_iterative_intersection"}
            initial_dir, final_dir, summary_dir = root / "initial", root / "final", root / "summaries"
            for p in [initial_dir, final_dir, summary_dir]: p.mkdir()

            def load(_reader, _previous, tool, branch, population):
                return fixture_metadata(tool, branch, population), root

            def markers(_reader, _previous, tool, branch, population):
                marker, panel = [], []
                groups = fixture_metadata(tool, branch, population).current_cluster_id.unique()
                for group in groups:
                    for symbol in ["Mertk", "Inpp5d", "Camk2a", "Gad1", "Gad2"]:
                        marker.append({"cluster_id": group, "gene_id": gene_lookup[symbol], "gene_symbol": symbol,
                                       "fraction_in": .8, "fraction_rest": .1, "marker_log2FC": 1., "pct_difference": .7})
                    for r in supplied.itertuples(index=False):
                        panel.append({"cluster_id": group, "gene_id": gene_lookup[r.gene_symbol], "gene_symbol": r.gene_symbol,
                                      "set": r.set, "fraction_detected": .5, "source_mean": 1.})
                return pd.DataFrame(marker), pd.DataFrame(panel), "SYNTHETIC_ONLY"

            reader = a.Inputs()
            with patch.object(a, "source_dataset", load), patch.object(a, "marker_source", markers):
                result = a.initial(reader, cfg, initial_dir)
            self.assertEqual(result["source_dataset_count"], 24)
            receipt = {"status": "ANNOTATION_INITIAL_COMPLETE_REVIEW_REQUIRED", "stage": "initial",
                       "output_sha256": {p.name: a.digest(p) for p in initial_dir.iterdir()}}
            (initial_dir / "COMPLETE.json").write_text(json.dumps(receipt))
            ct = pd.read_csv(initial_dir / "cluster_annotations.tsv", sep="\t", dtype={k: str for k in a.KEYS})
            atomic_members = pd.read_csv(initial_dir / "atomic_population_memberships.tsv.gz", sep="\t", dtype=str)
            union = pd.read_csv(initial_dir / "gene_panel.tsv", sep="\t", dtype=str)
            sr, wr = [], []
            for row in ct[ct.cluster_scope.eq("original_frozen")].itertuples(index=False):
                m = atomic_members[(atomic_members.tool == row.tool) & (atomic_members.branch == row.branch) &
                                   (atomic_members.population == row.population) & (atomic_members.cluster_id == row.cluster_id)]
                wt = m[m.condition == "WT"]
                for gene in union.itertuples(index=False):
                    record = {k: getattr(row, k) for k in a.KEYS}
                    value = 2. if gene.gene_symbol in ["Mertk", "Inpp5d"] else 1.
                    sr.append({**record, "gene_id": gene.gene_id, "gene_symbol": gene.gene_symbol, "n_nuclei": row.n_nuclei,
                               "mean_log1p_normalized": value if row.n_nuclei else 0., "fraction_detected": .5 if row.n_nuclei else 0.,
                               "pct_difference_vs_rest": .2 if gene.gene_symbol in ["Mertk", "Inpp5d"] else 0.})
                    wr.append({**record, "gene_id": gene.gene_id, "gene_symbol": gene.gene_symbol, "n_WT_nuclei": len(wt),
                               "n_WT_mice": wt['sample'].nunique(), "mean_log1p_normalized": value if len(wt) else np.nan,
                               "fraction_detected": .5 if len(wt) else np.nan})
            pd.DataFrame(sr).to_csv(summary_dir / "gene_cluster_summary.tsv", sep="\t", index=False)
            pd.DataFrame(wr).to_csv(summary_dir / "gene_cluster_summary_WT.tsv", sep="\t", index=False)
            (summary_dir / "COMPLETE.json").write_text(json.dumps({"status": "COMPLETE", "stage": "native_aggregate_summary",
                                                                  "scientific_fitting_performed": False,
                                                                  "raw_nucleus_count_aggregation_performed": False}))
            result = a.final(a.Inputs(), cfg, final_dir)
            self.assertFalse(result["atomic_cluster_IDs_and_initial_inference_masks_changed"])
            categories = pd.read_csv(final_dir / "category_definition.tsv", sep="\t")
            checks = pd.read_csv(final_dir / "paired_category_membership_checks.tsv", sep="\t")
            self.assertTrue(checks.same_nucleus_membership.all())
            self.assertTrue(categories.physical_layer_assignment.eq("UNDETERMINED").all())
            members = pd.read_csv(final_dir / "category_memberships.tsv.gz", sep="\t", dtype=str)
            for key, z in members.groupby(["tool", "branch", "population", "category_family"]):
                self.assertFalse(z.nucleus_id.duplicated().any())
                expected = atomic_members[(atomic_members.tool == key[0]) & (atomic_members.branch == key[1]) &
                                          (atomic_members.population == key[2])]
                self.assertEqual(set(z.nucleus_id), set(expected.nucleus_id))


if __name__ == "__main__":
    unittest.main()
