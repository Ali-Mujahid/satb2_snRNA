#!/usr/bin/env python3
"""Review existing clusters using source-bound reference and expression evidence.

This module never fits a model, aggregates nucleus RNA counts, or changes clusters.
Initial assignments use already produced marker summaries. Final assignments keep
atomic cluster IDs/masks and refine labels from full-cohort and WT-only evidence. Labels describe
evidence; confidence categories are not calibrated annotation probabilities.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

BRANCHES = ("original_full", "paired_original", "paired_cellbender")
POPULATIONS = ("allcells", "Excitatory", "Inhibitory", "CA1_ProS")
TOOLS = ("scanpy", "seurat")
KEYS = ["tool", "branch", "population", "cluster_scope", "cluster_id"]
NUCLEUS_LABEL_FIELDS = ["broad_label", "reviewed_label", "identity_label", "label_confidence", "layer_program_label",
                        "expression_program_label", "program_evidence_status", "physical_layer_assignment",
                        "identity_review_status", "uncertainty_reason", "annotation_evidence_status", "annotation_method"]
ALIAS = {"allcells": "log_global", "Excitatory": "annotate_subsets/EX",
         "Inhibitory": "annotate_subsets/IN", "CA1_ProS": "annotate_subsets/CA1_PROS_CANDIDATE"}
MARKER_FIELD = {"allcells": "log_harmony_r0p4", "Excitatory": "EX_harmony_r0p4",
                "Inhibitory": "IN_harmony_r0p4", "CA1_ProS": "CA1_PROS_CANDIDATE_harmony_r0p4"}
COARSE = {
    "Neuron": ["Snap25", "Rbfox3", "Syt1"],
    "Excitatory": ["Slc17a7", "Slc17a6", "Camk2a"],
    "Inhibitory": ["Gad1", "Gad2", "Slc32a1"],
    "Astrocyte": ["Aqp4", "Slc1a3", "Aldoc", "Atp1a2", "Ndrg2"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog"],
    "OPC": ["Pdgfra", "Cspg4"],
    "Microglia": ["C1qa", "C1qb", "P2ry12", "Mertk", "Inpp5d", "Hexb", "Lrmda"],
    "Vascular": ["Cldn5", "Pecam1", "Rgs5", "Pdgfrb", "Vtn", "Kcnj8"],
    "Ependymal": ["Foxj1", "Dynlrb2", "Ttr", "Kl"],
    "Other_immune": ["Ptprc", "Cd74", "Ctss"],
}
FINE_IN = {"Pvalb": ["Pvalb", "Gad1", "Gad2"], "Sst": ["Sst", "Gad1", "Gad2"],
           "Vip": ["Vip", "Gad1", "Gad2"], "Sncg": ["Sncg", "Gad1", "Gad2"],
           "Ndnf": ["Ndnf", "Gad1", "Gad2"], "Lamp5 Lhx6": ["Lamp5", "Lhx6", "Gad1"],
           "Lamp5": ["Lamp5", "Gad1", "Gad2"]}
LAYER_SETS = ("Superficial_associated_VEH_program", "Deep_associated_VEH_program")
NON_IDENTITY_SETS = {"Neuron", "Excitatory", "Inhibitory", "Astrocyte", "Oligodendrocyte",
                     "Microglia", "OPC", "Endothelial", "CA1_population_candidates", *LAYER_SETS}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def id_digest(ids):
    # Deliberately defines the empty-set convention explicitly.
    return hashlib.sha256("".join(x + "\n" for x in sorted(map(str, ids))).encode()).hexdigest()


class Inputs:
    def __init__(self):
        self.hashes = {}

    def bind(self, path):
        p = Path(path).resolve(strict=True)
        if not p.is_file():
            raise ValueError(f"Expected regular input file: {p}")
        before = digest(p)
        old = self.hashes.setdefault(str(p), before)
        if old != before:
            raise ValueError(f"Input changed since first read: {p}")
        return p

    def table(self, path, **kwargs):
        p = self.bind(path)
        d = pd.read_csv(p, sep="\t", **kwargs)
        if digest(p) != self.hashes[str(p)]:
            raise ValueError(f"Input changed during read: {p}")
        return d

    def obj(self, path):
        p = self.bind(path)
        d = json.loads(p.read_text())
        if digest(p) != self.hashes[str(p)]:
            raise ValueError(f"Input changed during read: {p}")
        return d

    def validate(self):
        for p, expected in self.hashes.items():
            if digest(p) != expected:
                raise ValueError(f"Input changed before completion: {p}")


def truth(s):
    x = s.astype(str).str.lower()
    if not x.isin(["true", "false", "1", "0"]).all():
        raise ValueError("Invalid boolean values")
    return x.isin(["true", "1"])


def sort_key(x):
    return (0, int(x)) if str(x).isdigit() else (1, str(x))


def safe(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))


def weighted_mean(x, weights):
    z = pd.to_numeric(x, errors="coerce")
    ok = np.isfinite(z)
    return float(np.average(z[ok], weights=np.asarray(weights)[ok])) if ok.any() else np.nan


def atlas_profile(md, field):
    counts = md[field].fillna("UNKNOWN").astype(str).value_counts()
    if not len(counts):
        return {"UNKNOWN": 1.0}
    return {k: float(v / len(md)) for k, v in counts.items()}


def entropy(profile):
    p = np.asarray(list(profile.values()), dtype=float)
    return float(-(p[p > 0] * np.log2(p[p > 0])).sum())


def source_dataset(inputs, previous, tool, branch, population):
    root = previous / "results"
    if tool == "scanpy":
        folder = root / tool / branch / population
        receipt = inputs.obj(folder / "COMPLETE.json")
        if receipt.get("status") != "COMPLETE" or receipt.get("count_branch") != branch or receipt.get("population") != population:
            raise ValueError(f"Incomplete/mismatched Scanpy producer: {folder}")
        field = "cluster_harmony"
    else:
        folder = root / tool / branch / ALIAS[population]
        receipt_folder = folder if population == "allcells" else folder.parent
        receipt = inputs.obj(receipt_folder / "COMPLETE.json")
        stage = "log_global" if population == "allcells" else "annotate_subsets"
        if receipt.get("status") != "EXECUTED_REVIEW_REQUIRED" or receipt.get("count_source") != branch or receipt.get("stage") != stage:
            raise ValueError(f"Incomplete/mismatched Seurat producer: {folder}")
        field = "log_harmony_r0p4"
    d = inputs.table(folder / "nucleus_metadata.tsv", dtype={"nucleus_id": str, "sample": str, "condition": str, field: str},
                     usecols=lambda c: c in {"nucleus_id", "sample", "condition", field} or c.startswith("reference_"))
    required = {"nucleus_id", "sample", "condition", field, "reference_broad_prediction", "reference_subclass_name",
                "reference_subclass_high_support", "reference_broad_high_support", "reference_supertype_name"}
    if not required.issubset(d) or d.nucleus_id.duplicated().any() or d[field].isna().any():
        raise ValueError(f"Invalid metadata schema/IDs/clusters: {folder}")
    shared = root / "shared" / branch
    receipt = inputs.obj(shared / "COMPLETE.json")
    if receipt.get("status") != "COMPLETE" or receipt.get("count_branch") != branch:
        raise ValueError(f"Incomplete shared producer: {shared}")
    membership = inputs.table(shared / f"membership_{population}.tsv", dtype=str)
    if membership.nucleus_id.duplicated().any() or set(d.nucleus_id) != set(membership.nucleus_id):
        raise ValueError(f"Metadata differs from authoritative {branch}/{population} membership")
    bind = d.set_index("nucleus_id").loc[membership.nucleus_id]
    for key in ("sample", "condition"):
        if not np.array_equal(bind[key].to_numpy(), membership[key].to_numpy()):
            raise ValueError(f"Conflicting {key} binding in {folder}")
    if d['sample'].nunique() != 12 or set(d.condition) != {"WT", "cKO"}:
        raise ValueError(f"Unexpected sample/condition coverage: {folder}")
    d["current_cluster_id"] = d[field].astype(str)
    return d, folder


def marker_source(inputs, previous, tool, branch, population):
    root = previous / "results"
    if tool == "scanpy":
        folder = root / tool / branch / population
        p = folder / "cluster_markers_harmony.tsv"
        d = inputs.table(p, dtype={"cluster": str})
        d = d.rename(columns={"cluster": "cluster_id", "pct_nz_group": "fraction_in", "pct_nz_reference": "fraction_rest",
                              "logfoldchanges": "marker_log2FC", "pvals_adj": "marker_adjusted_p"})
        panel = inputs.table(folder / "marker_panel_cluster_expression.tsv", dtype={"cluster": str})
        panel = panel.rename(columns={"cluster": "cluster_id", "marker_group": "set", "mean_log_normalized_expression": "source_mean"})
        scope = "full_cohort_panel_detection; sampled_mouse_balanced_rank_markers"
    else:
        marker_receipt = inputs.obj(root / tool / branch / "markers/COMPLETE.json")
        if marker_receipt.get("status") != "EXECUTED_REVIEW_REQUIRED" or marker_receipt.get("stage") != "markers" or marker_receipt.get("count_source") != branch:
            raise ValueError("Incomplete source Seurat marker producer")
        folder = root / tool / branch / "markers" / MARKER_FIELD[population]
        d = inputs.table(folder / "ALL_EXPLORATORY_MARKERS.tsv", dtype={"cluster": str})
        d = d.rename(columns={"cluster": "cluster_id", "pct.1": "fraction_in", "pct.2": "fraction_rest",
                              "avg_log2FC": "marker_log2FC", "p_val_adj": "marker_adjusted_p"})
        rows = []
        genes = inputs.table(root / "shared/original_full/gene_features.tsv", dtype=str).set_index("gene_id").gene_symbol.to_dict()
        for path in sorted(folder.glob("curated_*_NUMERICAL_VALUES.tsv")):
            z = inputs.table(path, dtype={"id": str, "features.plot": str})
            for _, row in z.iterrows():
                rows.append({"cluster_id": str(row["id"]), "set": path.name[len("curated_"):-len("_NUMERICAL_VALUES.tsv")],
                             "gene_id": row["features.plot"], "gene_symbol": genes[row["features.plot"]],
                             "fraction_detected": float(row["pct.exp"]) / 100., "source_mean": float(row["avg.exp"])})
        panel = pd.DataFrame(rows)
        scope = "sampled_mouse_balanced_Seurat_DotPlot_detection_and_rank_markers"
    if not {"gene_id", "gene_symbol", "cluster_id", "fraction_in", "fraction_rest", "marker_log2FC"}.issubset(d):
        raise ValueError("Source marker schema mismatch")
    d["pct_difference"] = pd.to_numeric(d.fraction_in) - pd.to_numeric(d.fraction_rest)
    d = d.sort_values(["cluster_id", "pct_difference", "marker_log2FC", "gene_id"], ascending=[True, False, False, True])
    # An exported top rank is evidence only; absent rows do not mean zero expression.
    return d, panel, scope


def program_scores(panel, panels, value="fraction_detected"):
    if panel.empty:
        return pd.DataFrame(columns=["cluster_id", "set", "program_score", "mapped_program_genes", "detectable_program_genes"])
    fields = list(dict.fromkeys(["cluster_id", "gene_id", "gene_symbol", "fraction_detected", value]))
    z = panel[fields].drop_duplicates(["cluster_id", "gene_id"])
    # Standardize each gene's detection fraction across this tool/population's clusters.
    # This describes relative expression programs, never anatomical position.
    z["gene_detection_z"] = z.groupby("gene_id")[value].transform(
        lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 1e-12 else s * 0)
    definitions = panels[["set", "gene_symbol"]].drop_duplicates()
    # Satb2 is the experimental target: it cannot define layer-associated memberships.
    definitions = definitions[~(definitions['set'].isin(LAYER_SETS) & definitions.gene_symbol.eq("Satb2"))]
    out = z.merge(definitions, on="gene_symbol", how="inner", validate="many_to_many")
    out["detectable"] = out.fraction_detected >= .10
    return out.groupby(["cluster_id", "set"], observed=True).agg(
        program_score=("gene_detection_z", "mean"), mapped_program_genes=("gene_id", "nunique"),
        detectable_program_genes=("detectable", "sum")).reset_index()


def program_assignment(scores, cluster, population):
    z = scores[scores.cluster_id.eq(str(cluster))].copy()
    if population == "CA1_ProS":
        z = z[z['set'].isin(LAYER_SETS)].set_index("set")
        if len(z) != 2 or (z.mapped_program_genes < 3).any():
            return "Unresolved layer-associated expression program", "INSUFFICIENT_PANEL_COVERAGE", np.nan, np.nan
        superficial = float(z.loc[LAYER_SETS[0], "program_score"])
        deep = float(z.loc[LAYER_SETS[1], "program_score"])
        delta = superficial - deep
        if abs(delta) >= .50 and max(superficial, deep) >= .30:
            winner = LAYER_SETS[0] if delta > 0 else LAYER_SETS[1]
            if int(z.loc[winner, "detectable_program_genes"]) >= 3:
                label = "Superficial-associated expression program enriched" if delta > 0 else "Deep-associated expression program enriched"
                return label, "RELATIVE_PROGRAM_SUPPORT_ONLY", superficial, deep
        return "Mixed or weak layer-associated expression programs", "NO_DOMINANT_PROGRAM", superficial, deep
    if population == "Inhibitory":
        z = z[~z['set'].isin(NON_IDENTITY_SETS)].sort_values(["program_score", "set"], ascending=[False, True])
        if len(z) and z.iloc[0].program_score >= .30 and z.iloc[0].detectable_program_genes >= 2:
            winner = z.iloc[0]['set']
            runner = str(z.iloc[1]['set']) if len(z) > 1 and z.iloc[0].program_score - z.iloc[1].program_score < .30 else ""
            label = f"{winner} expression program" + (f" / {runner} mixed" if runner else " enriched")
            return label, "WORKBOOK_PROGRAM_CANDIDATE_ONLY", np.nan, np.nan
        return "Unresolved inhibitory expression program", "NO_DOMINANT_PROGRAM", np.nan, np.nan
    return "Not assessed in this population", "NOT_APPLICABLE", np.nan, np.nan


def supported_genes(panel, cluster, symbols, threshold=.20):
    z = panel[panel.cluster_id.eq(str(cluster)) & panel.gene_symbol.isin(symbols)]
    return sorted(set(z.loc[z.fraction_detected >= threshold, "gene_symbol"]))


def ranked_signature(marker, cluster, symbols):
    z = marker[marker.cluster_id.eq(str(cluster))]
    z = z[(z.pct_difference >= .10) & (z.marker_log2FC >= .50)].head(20)
    return sorted(set(z.loc[z.gene_symbol.isin(symbols), "gene_symbol"]))


def review_cluster(md, cluster, marker, panel, scores, population, scope):
    b = atlas_profile(md, "reference_broad_prediction")
    s = atlas_profile(md, "reference_subclass_name")
    u = atlas_profile(md, "reference_supertype_name")
    broad, broad_frac = next(iter(b.items()))
    subclass, sub_frac = next(iter(s.items()))
    supertype, super_frac = next(iter(u.items()))
    broad_hs = float(truth(md.reference_broad_high_support).mean())
    sub_hs = float(truth(md.reference_subclass_high_support).mean())
    canonical = supported_genes(panel, cluster, COARSE.get(broad, []))
    ranked = {name: ranked_signature(marker, cluster, genes) for name, genes in COARSE.items() if name != "Neuron"}
    conflicting = [name for name, genes in ranked.items() if len(genes) >= 2 and name != broad and
                   not (broad == "OPC" and name == "Oligodendrocyte") and not (broad == "Oligodendrocyte" and name == "OPC")]
    dominant = sub_frac >= .80 and sub_hs >= .80
    if broad_frac >= .80 and broad_hs >= .80:
        identity = subclass if dominant else f"{broad}; mixed atlas subclasses"
        status = "REFERENCE_SUPPORTED" if dominant else "REFERENCE_MIXED_SUBCLASSES"
        confidence = "REFERENCE_AND_MARKERS_CORROBORATE_COARSE" if len(canonical) >= 2 else "REFERENCE_SUPPORTED_MARKER_REVIEW_REQUIRED"
    else:
        identity = "Mixed/uncertain " + " / ".join(list(b)[:2])
        status = "REFERENCE_MIXED_OR_LOW_SUPPORT"
        confidence = "UNCERTAIN"
    if conflicting:
        identity += "; " + "/".join(conflicting) + "-associated signature review"
        status = "CONFLICTING_MARKER_SIGNATURE"
        confidence = "CONFLICTING_EVIDENCE"
    fine_support = []
    if population == "Inhibitory":
        for name in sorted(FINE_IN, key=len, reverse=True):
            if name in subclass:
                fine_support = supported_genes(panel, cluster, FINE_IN[name], .10)
                break
    expression, program_status, superficial, deep = program_assignment(scores, cluster, population)
    if population == "CA1_ProS" and conflicting:
        expression = "Layer program not assigned: conflicting identity signatures"
        program_status = "CONFLICTING_IDENTITY_EVIDENCE"
    counts = md.groupby(["sample", "condition"], observed=True).size()
    eligible = counts[counts >= 20].reset_index(name="n_nuclei")
    n_wt = int(eligible.condition.eq("WT").sum()) if len(eligible) else 0
    n_cko = int(eligible.condition.eq("cKO").sum()) if len(eligible) else 0
    present = counts[counts >= 1].reset_index(name="n_nuclei")
    n_wt_nonempty = int(present.condition.eq("WT").sum())
    n_cko_nonempty = int(present.condition.eq("cKO").sum())
    uncertainties = []
    if broad_frac < .80: uncertainties.append("mixed atlas broad labels")
    if sub_frac < .80: uncertainties.append("mixed atlas subclasses")
    if len(canonical) < 2: uncertainties.append("fewer than two coarse canonical panel genes above descriptive detection flag")
    if conflicting: uncertainties.append("ranked nonmatching lineage signatures; ambient/multiplet/state causes unresolved")
    if len(md) < 100: uncertainties.append("small cluster (<100 nuclei); instability/replicate support needs review")
    if population == "Inhibitory":
        if n_wt_nonempty < 3 or n_cko_nonempty < 3: uncertainties.append("fewer than three mice/genotype with any nuclei; model eligibility determined separately")
        elif n_wt < 3 or n_cko < 3: uncertainties.append("rare-group marker coverage is sparse; all-nucleus IN DE eligibility is assessed separately")
    elif n_wt < 3 or n_cko < 3: uncertainties.append("fewer than three mice/genotype with >=20 nuclei; DE eligibility determined separately")
    if population == "CA1_ProS": uncertainties.append("CA1-ProS reference candidate mask includes prosubiculum; physical layer undetermined")
    top = marker[marker.cluster_id.eq(str(cluster))].head(7)
    return {"broad_label": broad if broad_frac >= .80 else "Mixed/uncertain",
            "reviewed_label": identity, "identity_label": identity, "label_confidence": confidence,
            "layer_program_label": expression if population == "CA1_ProS" else "Not assigned",
            "expression_program_label": expression, "program_evidence_status": program_status,
            "superficial_program_detection_z": superficial, "deep_program_detection_z": deep,
            "physical_layer_assignment": "UNDETERMINED", "anatomy_caveat": "Expression programs do not establish anatomical layer or morphology",
            "atlas_majority": subclass, "atlas_majority_fraction": sub_frac,
            "atlas_broad_majority_fraction": broad_frac, "atlas_supertype_majority": supertype,
            "atlas_supertype_majority_fraction": super_frac, "atlas_broad_high_support_fraction": broad_hs,
            "atlas_subclass_high_support_fraction": sub_hs, "atlas_subclass_entropy_bits": entropy(s),
            "atlas_supertype_entropy_bits": entropy(u), "atlas_fractions": json.dumps(s, sort_keys=True),
            "atlas_broad_fractions": json.dumps(b, sort_keys=True), "atlas_supertype_fractions": json.dumps(u, sort_keys=True),
            "marker_evidence": json.dumps({"coarse_detected": canonical, "fine_IN_detected": fine_support,
                                           "ranked_conflicting_signatures": {k: ranked[k] for k in conflicting},
                                           "top7_source_markers": top.gene_symbol.tolist()}, sort_keys=True),
            "identity_review_status": status, "uncertainty_reason": "; ".join(uncertainties) or "reference-marker expression evidence remains provisional",
            "annotation_evidence_status": scope, "n_mice_represented": md['sample'].nunique(),
            "n_WT_mice_ge20": n_wt, "n_cKO_mice_ge20": n_cko,
            "n_WT_mice_nonempty": n_wt_nonempty, "n_cKO_mice_nonempty": n_cko_nonempty,
            "membership_sha256": id_digest(md.nucleus_id), "annotation_method": "hierarchical_reference_plus_expression_evidence_review_v1"}


def build_panel(inputs, previous, marker_data, supplied):
    genes = inputs.table(previous / "results/shared/original_full/gene_features.tsv", dtype=str)
    if genes.gene_id.duplicated().any():
        raise ValueError("Source gene IDs are not unique")
    rows = []
    for row in supplied.itertuples(index=False):
        for gid in genes.loc[genes.gene_symbol.eq(row.gene_symbol), "gene_id"]:
            rows.append({"gene_id": gid, "gene_symbol": row.gene_symbol, "panel_role": "supplied_or_prior_canonical_panel",
                         "panel_set": row.set, "source": row.source, "selection_evidence": row.selection})
    additional = set(sum(COARSE.values(), []) + sum(FINE_IN.values(), []))
    for symbol in sorted(additional):
        for gid in genes.loc[genes.gene_symbol.eq(symbol), "gene_id"]:
            rows.append({"gene_id": gid, "gene_symbol": symbol, "panel_role": "additional_canonical_review",
                         "panel_set": "canonical_identity_review", "source": "explicit review flag panel; not identity truth",
                         "selection_evidence": "needed for broad/fine lineage conflict checks; no physical identity gate"})
    for (tool, population), d in marker_data.items():
        # Keep up to 20 previously ranked candidates/cluster to allow full-cohort support filtering for plots.
        for cluster, z in d.groupby("cluster_id", observed=True):
            for row in z.head(20).itertuples(index=False):
                rows.append({"gene_id": row.gene_id, "gene_symbol": row.gene_symbol,
                             "panel_role": "previous_ranked_representative_candidate", "panel_set": f"{tool}/{population}/{cluster}",
                             "source": "original_full existing descriptive cluster-marker rank", "selection_evidence":
                             f"pct_difference={row.pct_difference:.8g};log2FC={row.marker_log2FC:.8g};sampled_nucleus_exploration_not_genotype_DE"})
    audit = pd.DataFrame(rows).drop_duplicates()
    if not set(audit.gene_id).issubset(set(genes.gene_id)):
        raise ValueError("Annotation panel includes unknown source gene ID")
    union = audit[["gene_id", "gene_symbol"]].drop_duplicates("gene_id").sort_values("gene_id")
    return union, audit


def shared_cluster_masks(datasets, cfg):
    if cfg.get("cluster_filter_mode", "shared_iterative_intersection") != "shared_iterative_intersection":
        raise ValueError("Main All/EX filtering requires shared_iterative_intersection")
    masks, exclusions, iterations = {}, [], []
    for population in POPULATIONS:
        for branch in ["original_full", "paired_original"]:
            left, right = datasets[("scanpy", branch, population)], datasets[("seurat", branch, population)]
            if set(left.nucleus_id) != set(right.nucleus_id):
                raise ValueError("Cross-tool source population nucleus IDs differ before filtering")
            remaining = set(left.nucleus_id)
            iteration = 0
            if population in ["allcells", "Excitatory"]:
                while True:
                    bad, reasons = set(), {}
                    for tool, d in [("scanpy", left), ("seurat", right)]:
                        z = d[d.nucleus_id.isin(remaining)]
                        sizes = z.source_cluster_id.value_counts()
                        small = set(sizes[sizes < 100].index)
                        for row in z[z.source_cluster_id.isin(small)].itertuples(index=False):
                            bad.add(row.nucleus_id)
                            reasons.setdefault(row.nucleus_id, []).append({"tool": tool, "original_cluster_id": row.source_cluster_id,
                                                                        "n_nuclei_before_iteration": int(sizes[row.source_cluster_id])})
                    if not bad: break
                    iteration += 1
                    for nucleus_id in sorted(bad):
                        exclusions.append({"branch": branch, "population": population, "nucleus_id": nucleus_id,
                                           "exclusion_iteration": iteration, "exclusion_reason": "REQUESTED_SHARED_ALL_OR_EX_CLUSTER_LT100",
                                           "offending_tool_cluster_counts": json.dumps(reasons[nucleus_id], sort_keys=True)})
                    iterations.append({"branch": branch, "population": population, "iteration": iteration,
                                       "n_nuclei_before": len(remaining), "n_removed": len(bad), "n_nuclei_after": len(remaining - bad)})
                    remaining -= bad
                    if iteration > left.source_cluster_id.nunique() + right.source_cluster_id.nunique():
                        raise RuntimeError("Cluster-filter loop failed to converge")
            masks[(branch, population)] = remaining
            if branch == "paired_original":
                masks[("paired_cellbender", population)] = set(remaining)
                for e in [x.copy() for x in exclusions if x['branch'] == branch and x['population'] == population]:
                    e["branch"] = "paired_cellbender"
                    exclusions.append(e)
    for tool in TOOLS:
        for branch in BRANCHES:
            for population in POPULATIONS:
                d = datasets[(tool, branch, population)]
                z = d[d.nucleus_id.isin(masks[(branch, population)])]
                if population in ["allcells", "Excitatory"] and (z.source_cluster_id.value_counts() < 100).any():
                    raise ValueError("Retained atomic cluster violates shared minimum size")
    return masks, pd.DataFrame(exclusions, columns=["branch", "population", "nucleus_id", "exclusion_iteration", "exclusion_reason", "offending_tool_cluster_counts"]), pd.DataFrame(iterations, columns=["branch", "population", "iteration", "n_nuclei_before", "n_removed", "n_nuclei_after"])


def initial(inputs, cfg, out):
    previous = Path(cfg["previous_run_root"]).resolve(strict=True)
    # These fixed scopes reflect the requested complete matched-tool revision.
    for name, expected in (("branches", BRANCHES), ("tools", TOOLS), ("populations", POPULATIONS)):
        value = cfg.get(name, expected)
        if isinstance(value, dict): value = value.keys()
        if set(value) != set(expected): raise ValueError(f"Incomplete/unexpected {name}")
    if float(cfg.get("main_resolution", cfg.get("resolution", .4))) != .4 or cfg.get("representation", "harmony").lower() != "harmony":
        raise ValueError("Annotation contract is original Harmony resolution .4")
    supplied = inputs.table(previous / "results/inputs/marker_review_panels.tsv", dtype=str)
    if "eligible" in supplied: supplied = supplied[truth(supplied.eligible)].copy()
    clusters, nuclei, scores_all, source_marker_data, agreements, kept = [], [], [], {}, [], []
    datasets = {}
    for tool in TOOLS:
        for population in POPULATIONS:
            original, _ = source_dataset(inputs, previous, tool, "original_full", population)
            for branch in BRANCHES:
                d = original.copy(deep=False) if branch == "original_full" else source_dataset(inputs, previous, tool, branch, population)[0]
                source = original.set_index("nucleus_id").loc[d.nucleus_id]
                for key in ["sample", "condition"] + [c for c in d if c.startswith("reference_")]:
                    left, right = source[key].fillna("NA").astype(str).to_numpy(), d[key].fillna("NA").astype(str).to_numpy()
                    if not np.array_equal(left, right): raise ValueError(f"Source reference/sample changed: {tool}/{population}/{branch}/{key}")
                d = d.copy(deep=False)
                d["source_cluster_id"] = source.current_cluster_id.to_numpy()
                datasets[(tool, branch, population)] = d
    masks, excluded, filter_iterations = shared_cluster_masks(datasets, cfg)
    for tool in TOOLS:
        for population in POPULATIONS:
            original = datasets[(tool, "original_full", population)]
            m, panel, evidence_scope = marker_source(inputs, previous, tool, "original_full", population)
            source_marker_data[(tool, population)] = m
            score = program_scores(panel, supplied)
            score["tool"], score["branch"], score["population"], score["cluster_scope"] = tool, "original_full", population, "original_frozen"
            scores_all.append(score)
            original_map = {}
            for cluster in sorted(original.current_cluster_id.unique(), key=sort_key):
                z = original[original.current_cluster_id.eq(cluster)]
                original_map[cluster] = review_cluster(z, cluster, m, panel, score, population, evidence_scope)
            for branch in BRANCHES:
                d = datasets[(tool, branch, population)]
                shared_keep = d.nucleus_id.isin(masks[(branch, population)])
                frozen_keep = {}
                for cluster in sorted(original_map, key=sort_key):
                    all_z = d[d.source_cluster_id.eq(cluster)]
                    z = all_z[all_z.nucleus_id.isin(masks[(branch, population)])]
                    keep = bool(len(z))
                    frozen_keep[cluster] = keep
                    rec = {"tool": tool, "branch": branch, "population": population, "cluster_scope": "original_frozen",
                           "cluster_id": cluster, "source_cluster_id": cluster, "source_count_branch": "original_full",
                           "n_nuclei": len(z), "source_cluster_n_nuclei_before_filter": len(all_z),
                           "original_cluster_n_nuclei": int((original.current_cluster_id == cluster).sum()),
                           **original_map[cluster], "membership_sha256": id_digest(z.nucleus_id), "keep_for_analysis": keep,
                           "keep_for_plot": keep, "exclusion_reason": "" if keep else "EMPTY_AFTER_REQUESTED_SHARED_FILTER_OR_PAIRING",
                           "assignment_carried_from_original": branch != "original_full"}
                    clusters.append(rec)
                labels = pd.DataFrame([{"cluster_id": k, **v} for k, v in original_map.items()]).set_index("cluster_id")
                labels = labels.loc[d.source_cluster_id].reset_index(drop=True)
                cols = ["nucleus_id", "sample", "condition", "source_cluster_id", "current_cluster_id"] + [c for c in d if c.startswith("reference_")]
                nd = d[cols].reset_index(drop=True).copy()
                for key in NUCLEUS_LABEL_FIELDS: nd[key] = labels[key].to_numpy()
                nd["tool"], nd["branch"], nd["population"], nd["cluster_scope"] = tool, branch, population, "original_frozen"
                nd["cluster_id"] = d.source_cluster_id.to_numpy()
                nd["keep_for_analysis"] = shared_keep.to_numpy()
                nd["keep_for_plot"] = nd.keep_for_analysis
                nd["exclusion_reason"] = np.where(nd.keep_for_analysis, "", "REQUESTED_SHARED_ALL_OR_EX_CLUSTER_LT100")
                nd["source_count_branch"] = "original_full"
                nuclei.append(nd)
                kept.append({"tool": tool, "branch": branch, "population": population, "cluster_scope": "original_frozen",
                             "n_nuclei": len(d), "kept_nuclei": int(nd.keep_for_analysis.sum()), "kept_clusters": sum(frozen_keep.values()),
                             "excluded_nuclei": int((~nd.keep_for_analysis).sum()), "excluded_clusters": len(frozen_keep) - sum(frozen_keep.values()),
                             "kept_membership_sha256": id_digest(nd.loc[nd.keep_for_analysis, "nucleus_id"])})
                if branch == "original_full": continue
                if tool == "seurat":
                    # The accepted paired Seurat marker stage explicitly used FROZEN
                    # ORIGINAL labels, not its fresh branch cluster numbers. Never
                    # attach those marker rows to unrelated fresh-numbered groups.
                    fm = pd.DataFrame(columns=["cluster_id", "gene_id", "gene_symbol", "pct_difference", "marker_log2FC"])
                    fp = pd.DataFrame(columns=["cluster_id", "gene_id", "gene_symbol", "fraction_detected"])
                    fs = "REFERENCE_METADATA_ONLY_FRESH_PAIRED_SEURAT; existing_paired_marker_artifacts_use_frozen_original_labels"
                else:
                    fm, fp, fs = marker_source(inputs, previous, tool, branch, population)
                sp = program_scores(fp, supplied)
                sp["tool"], sp["branch"], sp["population"], sp["cluster_scope"] = tool, branch, population, "fresh_branch"
                if not sp.empty: scores_all.append(sp)
                fresh_map = {}
                for cluster in sorted(d.current_cluster_id.unique(), key=sort_key):
                    z = d[d.current_cluster_id.eq(cluster)]
                    fresh_map[cluster] = review_cluster(z, cluster, fm, fp, sp, population, fs)
                    display = not (population in ["allcells", "Excitatory"] and len(z) < 100)
                    mix = z.source_cluster_id.value_counts(normalize=True).to_dict()
                    clusters.append({"tool": tool, "branch": branch, "population": population, "cluster_scope": "fresh_branch",
                                     "cluster_id": cluster, "source_cluster_id": "MULTIPLE_ORIGINAL_CLUSTERS",
                                     "source_count_branch": branch, "n_nuclei": len(z), **fresh_map[cluster],
                                     "keep_for_analysis": False, "keep_for_plot": display,
                                     "exclusion_reason": "DISPLAY_ONLY_FRESH_BRANCH_PARTITION",
                                     "original_cluster_membership_fractions": json.dumps(mix, sort_keys=True),
                                     "assignment_carried_from_original": False})
                flabels = pd.DataFrame([{"cluster_id": k, **v} for k, v in fresh_map.items()]).set_index("cluster_id").loc[d.current_cluster_id].reset_index(drop=True)
                fn = d[cols].reset_index(drop=True).copy()
                for key in NUCLEUS_LABEL_FIELDS: fn[key] = flabels[key].to_numpy()
                fn["tool"], fn["branch"], fn["population"], fn["cluster_scope"] = tool, branch, population, "fresh_branch"
                fn["cluster_id"] = d.current_cluster_id.to_numpy()
                fn["keep_for_analysis"] = False
                fn["keep_for_plot"] = fn.cluster_id.map(d.current_cluster_id.value_counts()).ge(100).to_numpy() if population in ["allcells", "Excitatory"] else True
                fn["exclusion_reason"] = "DISPLAY_ONLY_FRESH_BRANCH_PARTITION"
                fn["source_count_branch"] = branch
                nuclei.append(fn)
            pair_a = datasets[(tool, "paired_original", population)].set_index("nucleus_id")
            pair_b = datasets[(tool, "paired_cellbender", population)].set_index("nucleus_id")
            if set(pair_a.index) != set(pair_b.index) or not np.array_equal(pair_a.source_cluster_id, pair_b.loc[pair_a.index].source_cluster_id):
                raise ValueError("Paired raw/corrected frozen memberships are not identical")
            agreements.append({"tool": tool, "population": population, "paired_original_nuclei": len(pair_a),
                               "paired_cellbender_nuclei": len(pair_b), "same_frozen_IDs": True,
                               "same_frozen_labels": True, "same_analysis_exclusion_mask": True})
    ct, nt = pd.DataFrame(clusters), pd.concat(nuclei, ignore_index=True)
    if ct.duplicated(KEYS).any() or nt.duplicated(KEYS + ["nucleus_id"]).any(): raise ValueError("Duplicate annotation keys")
    # Annotations preserve excluded cells, and exactly enumerate every produced cluster.
    for tool in TOOLS:
        for population in POPULATIONS:
            for branch in BRANCHES:
                q = nt[(nt.tool == tool) & (nt.population == population) & (nt.branch == branch) & (nt.cluster_scope == "original_frozen")]
                if set(q.nucleus_id) != set(datasets[(tool, branch, population)].nucleus_id): raise ValueError("Incomplete annotation nucleus roster")
    union, audit = build_panel(inputs, previous, source_marker_data, supplied)
    write_table(out / "cluster_annotations.tsv", ct)
    write_table(out / "nucleus_annotations.tsv.gz", nt)
    write_table(out / "annotation_program_scores.tsv", pd.concat(scores_all, ignore_index=True))
    write_table(out / "annotation_retention.tsv", pd.DataFrame(kept))
    write_table(out / "paired_frozen_membership_checks.tsv", pd.DataFrame(agreements))
    write_table(out / "cluster_filter_exclusions.tsv.gz", excluded)
    write_table(out / "cluster_filter_iterations.tsv", filter_iterations)
    write_table(out / "gene_panel.tsv", union)
    write_table(out / "plot_gene_panel.tsv", union)
    write_table(out / "gene_panel_selection_audit.tsv", audit)
    source_genes = inputs.table(previous / "results/shared/original_full/gene_features.tsv", dtype=str)
    panel_membership = supplied.merge(source_genes, on="gene_symbol", how="left", validate="many_to_many")
    panel_membership["mapping_status"] = np.where(panel_membership.gene_id.isna(), "MISSING_EXACT_SOURCE_SYMBOL", "MAPPED_STABLE_ID")
    write_table(out / "panel_membership.tsv", panel_membership)
    atomic = ct[ct.cluster_scope.eq("original_frozen")][["tool", "branch", "population", "cluster_id", "n_nuclei", "membership_sha256", "keep_for_analysis"]].copy()
    atomic["population_key"] = [f"{r.tool}__{r.population}__H0p4__{r.cluster_id}" for r in atomic.itertuples(index=False)]
    atomic["model_enabled"] = atomic.keep_for_analysis
    write_table(out / "atomic_population_index.tsv", atomic)
    members = nt[nt.cluster_scope.eq("original_frozen") & nt.keep_for_analysis][["tool", "branch", "population", "cluster_id", "nucleus_id", "sample", "condition"]].copy()
    members["population_key"] = [f"{r.tool}__{r.population}__H0p4__{r.cluster_id}" for r in members.itertuples(index=False)]
    write_table(out / "atomic_population_memberships.tsv.gz", members)
    # All marker sources are retained separately for audit; no unsupported zero padding.
    missing = supplied[~supplied.gene_symbol.isin(union.gene_symbol)]
    write_table(out / "supplied_panel_missing_genes.tsv", missing)
    return {"n_cluster_rows": len(ct), "n_nucleus_annotation_rows": len(nt),
            "original_full_cluster_rows": len(ct[(ct.branch == "original_full") & (ct.cluster_scope == "original_frozen")]),
            "source_gene_panel_size": len(union), "source_dataset_count": len(datasets), "paired_frozen_binding_checks": 8,
            "assignment_definition": "original atomic clusters frozen; shared iterative All/EX intersection; paired-original/CellBender same restricted IDs; fresh branch display evidence separate",
            "full_cohort_corroboration_pending": True}


def final(inputs, cfg, out):
    initial_dir = Path(cfg["initial_annotation_dir"]).resolve(strict=True)
    previous = inputs.obj(initial_dir / "COMPLETE.json")
    if previous.get("status") != "ANNOTATION_INITIAL_COMPLETE_REVIEW_REQUIRED" or previous.get("stage") != "initial":
        raise ValueError("Final annotation requires completed initial fixed map")
    for name, expected in previous.get("output_sha256", {}).items():
        if digest(inputs.bind(initial_dir / name)) != expected: raise ValueError(f"Initial annotation output changed: {name}")
    summary_dir = Path(cfg["marker_summary_dir"]).resolve(strict=True)
    summary_receipt = inputs.obj(summary_dir / "COMPLETE.json")
    if summary_receipt.get("status") != "COMPLETE" or summary_receipt.get("stage") != "native_aggregate_summary":
        raise ValueError("Full-cohort marker summary producer is incomplete")
    if summary_receipt.get("scientific_fitting_performed") is not False or summary_receipt.get("raw_nucleus_count_aggregation_performed") is not False:
        raise ValueError("Marker summary receipt lacks the declared arithmetic-only evidence scope")
    for name, expected in summary_receipt.get("output_sha256", {}).items():
        if digest(inputs.bind(summary_dir / name)) != expected: raise ValueError(f"Full-cohort marker summary output changed: {name}")
    summary = inputs.table(summary_dir / "gene_cluster_summary.tsv", dtype={k: str for k in KEYS + ["gene_id", "gene_symbol"]})
    wt = inputs.table(summary_dir / "gene_cluster_summary_WT.tsv", dtype={k: str for k in KEYS + ["gene_id", "gene_symbol"]})
    if not set(KEYS + ["gene_id", "gene_symbol", "n_nuclei", "mean_log1p_normalized", "fraction_detected"]).issubset(summary):
        raise ValueError("Full-cohort marker summary schema mismatch")
    if not set(KEYS + ["gene_id", "gene_symbol", "n_WT_nuclei", "n_WT_mice", "mean_log1p_normalized", "fraction_detected"]).issubset(wt):
        raise ValueError("WT-only marker summary schema mismatch")
    if summary.duplicated(KEYS + ["gene_id"]).any() or wt.duplicated(KEYS + ["gene_id"]).any(): raise ValueError("Duplicate full-cohort/WT gene-summary keys")
    ct = inputs.table(initial_dir / "cluster_annotations.tsv", dtype={k: str for k in KEYS})
    nt = inputs.table(initial_dir / "nucleus_annotations.tsv.gz", dtype={k: str for k in KEYS + ["nucleus_id"]})
    supplied = inputs.table(initial_dir / "panel_membership.tsv", dtype=str)
    panels = supplied[["set", "gene_symbol"]].drop_duplicates()
    union = inputs.table(initial_dir / "gene_panel.tsv", dtype=str)
    expected_genes = set(union.gene_id)
    allowed_atomic = ct[ct.cluster_scope.eq("original_frozen")][KEYS]
    if len(summary[KEYS].drop_duplicates().merge(allowed_atomic, on=KEYS, how="left", indicator=True).query('_merge == "left_only"')):
        raise ValueError("Native full-cohort summaries include unknown atomic group")
    by_summary = {key: z for key, z in summary.groupby(KEYS, observed=True)}
    by_wt = {key: z for key, z in wt.groupby(KEYS, observed=True)}
    full_audit, program_rows, original_final, changes = [], [], {}, []
    for tool in TOOLS:
        for population in POPULATIONS:
            base = ct[(ct.tool == tool) & (ct.population == population) & ct.branch.eq("original_full") & ct.cluster_scope.eq("original_frozen")]
            wp = wt[(wt.tool == tool) & (wt.population == population) & wt.branch.eq("original_full") & wt.cluster_scope.eq("original_frozen")].copy()
            scores = program_scores(wp, panels, value="mean_log1p_normalized")
            scores["tool"], scores["population"], scores["branch"], scores["cluster_scope"] = tool, population, "original_full", "original_frozen"
            scores["score_definition"] = "equal-gene mean of WT-only mean-log1p expression z across original atomic clusters"
            scores["physical_layer_assignment"] = "UNDETERMINED"
            program_rows.append(scores)
            for row in base.itertuples(index=False):
                key = tuple(getattr(row, k) for k in KEYS)
                z = by_summary.get(key, pd.DataFrame())
                wz = by_wt.get(key, pd.DataFrame())
                if int(row.n_nuclei) > 0 and (z.empty or wz.empty): raise ValueError("Missing original full/WT atomic marker evidence")
                if len(z) and set(z.gene_id) != expected_genes: raise ValueError("Original native marker summary gene-panel coverage differs")
                if len(wz) and set(wz.gene_id) != expected_genes: raise ValueError("WT native marker summary gene-panel coverage differs")
                native_canonical = sorted(z.loc[z.gene_symbol.isin(COARSE.get(row.broad_label, [])) & z.fraction_detected.ge(.20), "gene_symbol"].unique()) if len(z) else []
                new = {k: getattr(row, k) for k in NUCLEUS_LABEL_FIELDS}
                new["initial_reviewed_label"] = row.reviewed_label
                new["initial_layer_program_label"] = row.layer_program_label
                conflicts = {}
                if len(z) and "pct_difference_vs_rest" in z:
                    for label, genes in COARSE.items():
                        if label in ["Neuron", row.broad_label]: continue
                        flag = z[z.gene_symbol.isin(genes) & z.fraction_detected.ge(.20) & z.pct_difference_vs_rest.ge(.10)]
                        if len(flag.gene_id.unique()) >= 2: conflicts[label] = sorted(flag.gene_symbol.unique())
                if conflicts:
                    new["reviewed_label"] = str(row.atlas_majority) + "; " + "/".join(sorted(conflicts)) + "-associated signature review"
                    new["identity_label"] = new["reviewed_label"]
                    new["identity_review_status"] = "FULL_COHORT_CONFLICTING_MARKER_SIGNATURE"
                    new["label_confidence"] = "CONFLICTING_EVIDENCE"
                    new["uncertainty_reason"] = "Full-cohort nonmatching lineage expression; ambient/multiplet/state causes unresolved. " + str(row.uncertainty_reason)
                elif len(native_canonical) >= 2 and row.identity_review_status != "CONFLICTING_MARKER_SIGNATURE":
                    new["label_confidence"] = "REFERENCE_AND_FULL_COHORT_COARSE_MARKERS_CORROBORATE"
                    new["uncertainty_reason"] = "Source atlas finer subclasses, local state, ambient and multiplet explanations remain provisional; full-cohort coarse marker support is descriptive."
                if int(row.n_nuclei) < 100 and row.population in ["Inhibitory", "CA1_ProS"]:
                    if new["label_confidence"] != "CONFLICTING_EVIDENCE": new["label_confidence"] = "RARE_REFERENCE_MARKER_EVIDENCE_REVIEW_REQUIRED"
                if len(wz):
                    wt_n = int(wz.n_WT_nuclei.iloc[0])
                    wt_mice = int(wz.n_WT_mice.iloc[0])
                    if not wz.n_WT_nuclei.eq(wt_n).all() or not wz.n_WT_mice.eq(wt_mice).all(): raise ValueError("Inconsistent WT atomic denominators")
                else: wt_n, wt_mice = 0, 0
                expression, ps, superficial, deep = program_assignment(scores, row.cluster_id, population)
                if population in ["Inhibitory", "CA1_ProS"] and (wt_n < 20 or wt_mice < 3):
                    expression, ps = "Unresolved expression program: sparse WT marker evidence", "SPARSE_WT_MARKER_EVIDENCE"
                if population == "CA1_ProS" and (conflicts or row.identity_review_status == "CONFLICTING_MARKER_SIGNATURE"):
                    expression, ps = "Layer program not assigned: conflicting identity signatures", "CONFLICTING_IDENTITY_EVIDENCE"
                new["expression_program_label"] = expression
                new["layer_program_label"] = expression if population == "CA1_ProS" else "Not assigned"
                new["program_evidence_status"] = ps
                new["annotation_evidence_status"] = "original_full_reference_and_native_full_cohort_markers; program_definition_from_original_full_WT_only"
                new["program_annotation_reference_branch"] = "original_full"
                new["program_assignment_condition"] = "WT"
                new["n_WT_program_evidence_nuclei"] = wt_n
                new["n_WT_program_evidence_mice"] = wt_mice
                new["superficial_program_WT_expression_z"] = superficial
                new["deep_program_WT_expression_z"] = deep
                new["final_identity_conflicting_genes"] = json.dumps(conflicts, sort_keys=True)
                new["final_coarse_canonical_genes"] = ";".join(native_canonical)
                original_final[(tool, population, row.cluster_id)] = new
                changes.append({"tool": tool, "population": population, "cluster_id": row.cluster_id,
                                "initial_reviewed_label": row.reviewed_label, "final_reviewed_label": new["reviewed_label"],
                                "initial_program": row.expression_program_label, "final_program": expression,
                                "atomic_cluster_ID_changed": False, "inference_membership_changed_after_initial": False})
    for row in ct.itertuples(index=False):
        key = tuple(getattr(row, k) for k in KEYS)
        z = by_summary.get(key, pd.DataFrame())
        if row.cluster_scope != "original_frozen":
            # Fresh branch summaries may be deferred; its source-specific marker review remains visible.
            full_audit.append({**{k: getattr(row, k) for k in KEYS}, "full_cohort_marker_status": "SOURCE_MARKER_REVIEW_ONLY_FRESH_BRANCH",
                               "full_cohort_canonical_genes": "", "full_cohort_conflicting_genes": "", "atomic_cluster_ID_changed": False})
            continue
        if int(row.n_nuclei) == 0:
            status, genes = "EMPTY_AFTER_SHARED_FILTER_OR_PAIRING_RETAINED_IN_REGISTRY", []
        else:
            if z.empty or not z.n_nuclei.eq(int(row.n_nuclei)).all(): raise ValueError("Full-cohort marker summary missing/mismatched cluster denominator")
            if set(z.gene_id) != expected_genes: raise ValueError("Native marker summary branch gene-panel coverage differs")
            if not z.fraction_detected.between(0, 1).all() or not np.isfinite(z.mean_log1p_normalized).all(): raise ValueError("Invalid native summary metrics")
            genes = sorted(z.loc[z.gene_symbol.isin(COARSE.get(row.broad_label, [])) & z.fraction_detected.ge(.20), "gene_symbol"].unique())
            status = "FULL_COHORT_COARSE_MARKERS_CORROBORATE" if len(genes) >= 2 else "FULL_COHORT_LOW_COARSE_CANONICAL_SUPPORT"
        conflicting = []
        if not z.empty and "pct_difference_vs_rest" in z:
            for label, panel in COARSE.items():
                if label in ["Neuron", row.broad_label]: continue
                flag = z[z.gene_symbol.isin(panel) & z.fraction_detected.ge(.20) & z.pct_difference_vs_rest.ge(.10)]
                if len(flag.gene_id.unique()) >= 2: conflicting.extend(label + ":" + str(x) for x in flag.gene_symbol)
            if conflicting: status = "FULL_COHORT_NONMATCHING_SIGNATURE_REVIEW_REQUIRED"
        full_audit.append({**{k: getattr(row, k) for k in KEYS}, "full_cohort_marker_status": status,
                           "full_cohort_canonical_genes": ";".join(genes), "full_cohort_conflicting_genes": ";".join(conflicting),
                           "atomic_cluster_ID_changed": False})
    audit = pd.DataFrame(full_audit)
    ct = ct.merge(audit, on=KEYS, validate="one_to_one")
    for i, row in ct.iterrows():
        if row.cluster_scope != "original_frozen": continue
        fields = original_final[(row.tool, row.population, row.cluster_id)]
        for k, value in fields.items(): ct.loc[i, k] = value
    # Final map is defined once from original-full WT evidence, then carried to paired cohorts.
    fresh_keys = ct.cluster_scope.eq("fresh_branch")
    ct.loc[fresh_keys, "annotation_evidence_status"] = ct.loc[fresh_keys, "annotation_evidence_status"].astype(str) + "; secondary_fresh_partition_source_review_not_inference_category"
    drop = [x for x in NUCLEUS_LABEL_FIELDS if x in nt]
    nt = nt.drop(columns=drop).merge(ct[KEYS + NUCLEUS_LABEL_FIELDS + ["full_cohort_marker_status"]], on=KEYS, validate="many_to_one")
    write_table(out / "cluster_annotations.tsv", ct)
    write_table(out / "nucleus_annotations.tsv.gz", nt)
    write_table(out / "full_cohort_marker_corroboration.tsv", audit)
    write_table(out / "annotation_program_scores_WT.tsv", pd.concat(program_rows, ignore_index=True))
    write_table(out / "label_refinement_audit.tsv", pd.DataFrame(changes))
    definitions, category_members = [], []
    for (tool, branch, population), block in ct[ct.cluster_scope.eq("original_frozen")].groupby(["tool", "branch", "population"], observed=True):
        block = block[truth(block.keep_for_analysis)]
        family_fields = [("REVIEWED_COARSE", "broad_label"), ("REVIEWED_IDENTITY", "reviewed_label")]
        if population == "CA1_ProS": family_fields.append(("CA1_LAYER_EXPRESSION_PROGRAM", "layer_program_label"))
        if population == "Inhibitory": family_fields.append(("IN_EXPRESSION_PROGRAM", "expression_program_label"))
        cells = nt[(nt.tool == tool) & (nt.branch == branch) & (nt.population == population) & nt.cluster_scope.eq("original_frozen") & truth(nt.keep_for_analysis)]
        for family, field in family_fields:
            seen = set()
            for label, group in block.groupby(field, observed=True, dropna=False):
                ids = sorted(group.cluster_id.astype(str), key=sort_key)
                if seen.intersection(ids): raise ValueError("Overlapping atomic definitions within annotation category family")
                seen.update(ids)
                members = cells[cells.cluster_id.isin(ids)][["nucleus_id", "sample", "condition", "source_cluster_id"]].copy()
                for k, v in [("tool", tool), ("branch", branch), ("population", population), ("category_family", family), ("label", str(label))]: members[k] = v
                definitions.append({"tool": tool, "branch": branch, "population": population, "category_family": family, "label": str(label),
                                    "member_cluster_ids": json.dumps(ids), "cluster_scope": "original_frozen", "source_branch": "original_full",
                                    "assignment_basis": "hierarchical_reference+original_full_WT_program+full_cohort_marker_review",
                                    "program_assignment_count_branch": "original_full", "program_assignment_condition": "WT",
                                    "n_nuclei": len(members), "membership_sha256": id_digest(members.nucleus_id),
                                    "physical_layer_assignment": "UNDETERMINED"})
                category_members.append(members)
            if seen != set(block.cluster_id): raise ValueError("Not all retained atoms reviewed within category family")
    define = pd.DataFrame(definitions)
    cm = pd.concat(category_members, ignore_index=True)
    cat_keys = ["tool", "branch", "population", "category_family", "label"]
    if define.duplicated(cat_keys).any() or cm.duplicated(cat_keys + ["nucleus_id"]).any(): raise ValueError("Duplicate final category definitions/members")
    for family, members in cm.groupby(["tool", "branch", "population", "category_family"], observed=True):
        if members.nucleus_id.duplicated().any(): raise ValueError("Overlapping nucleus assignments within category family")
    paired_checks = []
    for (tool, population, family), group in define.groupby(["tool", "population", "category_family"], observed=True):
        left = group[group.branch.eq("paired_original")].set_index("label")
        right = group[group.branch.eq("paired_cellbender")].set_index("label")
        if set(left.index) != set(right.index): raise ValueError("Paired count branches have different final categories")
        for label, row in left.iterrows():
            other = right.loc[label]
            for col in ["member_cluster_ids", "n_nuclei", "membership_sha256", "assignment_basis"]:
                if row[col] != other[col]: raise ValueError("Paired final category definitions/memberships differ")
            paired_checks.append({"tool": tool, "population": population, "category_family": family, "label": label,
                                  "same_atomic_IDs": True, "same_nucleus_membership": True,
                                  "n_nuclei": int(row.n_nuclei), "membership_sha256": row.membership_sha256,
                                  "program_assignment_count_branch": "original_full", "program_assignment_condition": "WT"})
    write_table(out / "category_definition.tsv", define)
    write_table(out / "category_memberships.tsv.gz", cm)
    write_table(out / "paired_category_membership_checks.tsv", pd.DataFrame(paired_checks))
    for name in ["annotation_program_scores.tsv", "annotation_retention.tsv", "paired_frozen_membership_checks.tsv",
                 "cluster_filter_exclusions.tsv.gz", "cluster_filter_iterations.tsv", "atomic_population_index.tsv",
                 "atomic_population_memberships.tsv.gz", "gene_panel.tsv", "plot_gene_panel.tsv", "panel_membership.tsv",
                 "gene_panel_selection_audit.tsv", "supplied_panel_missing_genes.tsv"]:
        p = inputs.bind(initial_dir / name)
        (out / name).write_bytes(p.read_bytes())
    return {"n_cluster_rows": len(ct), "n_nucleus_annotation_rows": len(nt), "atomic_cluster_IDs_and_initial_inference_masks_changed": False,
            "reviewed_label_map_refined_before_category_DE": True, "n_category_definitions": len(define),
            "native_full_cohort_marker_audit_rows": len(audit), "full_cohort_corroboration_pending": False,
            "scientific_acceptance": "REVIEW_REQUIRED; original reference and expression evidence are provisional"}


def write_table(path, data):
    data.to_csv(path, sep="\t", index=False, compression="gzip" if str(path).endswith(".gz") else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--stage", choices=["initial", "final"], default="initial")
    args = ap.parse_args()
    out = args.outdir.resolve()
    if out.exists() and any(out.iterdir()): raise ValueError("Refusing nonempty output directory; choose a new native stage directory")
    out.mkdir(parents=True, exist_ok=True)
    inputs = Inputs()
    cfg = inputs.obj(args.config)
    result = initial(inputs, cfg, out) if args.stage == "initial" else final(inputs, cfg, out)
    inputs.validate()
    outputs = {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file()}
    receipt = {"status": "ANNOTATION_INITIAL_COMPLETE_REVIEW_REQUIRED" if args.stage == "initial" else "ANNOTATION_FINAL_COMPLETE_REVIEW_REQUIRED",
               "stage": args.stage, "source_code_sha256": digest(Path(__file__)), "config_sha256": digest(args.config),
               "input_sha256": inputs.hashes, "output_sha256": outputs,
               "scientific_model_fit_performed": False, "raw_nucleus_RNA_aggregation_performed": False,
               "source_cluster_IDs_and_MapMyCells_predictions_modified": False, "calibrated_annotation_confidence": False,
               "physical_layer_assignment_performed": False, **result}
    (out / "COMPLETE.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": receipt["status"], "outdir": str(out), **result}, allow_nan=False))


if __name__ == "__main__":
    main()
