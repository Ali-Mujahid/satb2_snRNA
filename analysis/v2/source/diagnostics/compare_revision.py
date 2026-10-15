#!/usr/bin/env python3
"""Native final review of produced annotation/count/DESeq2/diagnostic artifacts.

No scientific fitting, graph construction, or raw-nucleus aggregation occurs here.
Cross-tool DE uses the same DESeq2 implementation; matched fine populations can
differ. Availability, low coverage, mixed labels, and untested endpoints stay visible.
"""
from __future__ import annotations
import argparse
from collections import OrderedDict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr, pearsonr

BRANCHES = ["original_full", "paired_original", "paired_cellbender"]
TOOLS = ["scanpy", "seurat"]
POPS = ["allcells", "Excitatory", "Inhibitory", "CA1_ProS"]
SEMANTIC = ["tool", "population", "analysis_kind", "category_family", "analysis_label"]
STATS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
SEVEN = ["gene", *STATS]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""): h.update(b)
    return h.hexdigest()


def member_sha(ids):
    return hashlib.sha256("".join(str(i) + "\n" for i in sorted(ids)).encode()).hexdigest()


def truth(x):
    s = x.astype(str).str.lower()
    if not s.isin(["true", "false", "1", "0"]).all(): raise ValueError("Invalid boolean flag")
    return s.isin(["true", "1"])


def identity(row):
    return tuple(str(row[k]) if pd.notna(row[k]) else "" for k in SEMANTIC)


class Evidence:
    def __init__(self): self.hashes, self.expected = {}, {}
    def bind(self, path):
        p = Path(path).resolve(strict=True)
        value = sha(p)
        if str(p) in self.expected and self.expected[str(p)] != value: raise ValueError(f"Producer output changed: {p}")
        if str(p) in self.hashes and self.hashes[str(p)] != value: raise ValueError(f"Input changed during review: {p}")
        self.hashes[str(p)] = value
        return p
    def table(self, path, **kwargs):
        p = self.bind(path)
        z = pd.read_csv(p, sep="\t", keep_default_na=False, na_values=["NA"], **kwargs)
        if sha(p) != self.hashes[str(p)]: raise ValueError(f"Input changed while parsing: {p}")
        return z
    def json(self, path):
        p = self.bind(path); z = json.loads(p.read_text())
        if sha(p) != self.hashes[str(p)]: raise ValueError(f"Input changed while parsing: {p}")
        return z
    def stage(self, root, status, stage=None, branch=None):
        root = Path(root).resolve(strict=True)
        r = self.json(root / "COMPLETE.json")
        if r.get("status") != status or (stage is not None and r.get("stage") != stage) or (branch is not None and r.get("branch") != branch):
            raise ValueError(f"Incomplete/mismatched required producer: {root}")
        for name, expected in r.get("output_sha256", {}).items():
            p = (root / name).resolve(strict=True)
            if not p.is_relative_to(root): raise ValueError("Producer output escapes stage root")
            self.expected[str(p)] = expected
        return r
    def finish(self):
        for path, expected in self.hashes.items():
            if sha(path) != expected: raise ValueError(f"Input changed before completion: {path}")


def effect_metrics(left, right):
    z = left.set_index("gene_id").join(right.set_index("gene_id"), how="inner", lsuffix="_left", rsuffix="_right")
    a, b = z.log2FoldChange_left.to_numpy(float), z.log2FoldChange_right.to_numpy(float)
    finite = np.isfinite(a) & np.isfinite(b)
    sign = np.sign(a) == np.sign(b)
    p = np.isfinite(z.pvalue_left) & np.isfinite(z.pvalue_right)
    fdr = np.isfinite(z.padj_left) & np.isfinite(z.padj_right)
    significant_a = z.padj_left.lt(.05) & np.isfinite(z.padj_left)
    significant_b = z.padj_right.lt(.05) & np.isfinite(z.padj_right)
    both = significant_a & significant_b
    union = (significant_a | significant_b) & fdr
    variable = finite.sum() >= 3 and np.std(a[finite]) > 0 and np.std(b[finite]) > 0
    return {"left_gene_rows": len(left), "right_gene_rows": len(right), "common_gene_rows": len(z),
            "common_finite_effects": int(finite.sum()), "common_finite_pvalues": int(p.sum()), "common_finite_FDR": int(fdr.sum()),
            "effect_Spearman": float(spearmanr(a[finite], b[finite]).statistic) if variable else np.nan,
            "effect_Pearson": float(pearsonr(a[finite], b[finite]).statistic) if variable else np.nan,
            "effect_same_sign_fraction": float(sign[finite].mean()) if finite.any() else np.nan,
            "median_abs_effect_difference": float(np.median(np.abs(a[finite] - b[finite]))) if finite.any() else np.nan,
            "FDR05_left_on_common": int(significant_a.sum()), "FDR05_right_on_common": int(significant_b.sum()),
            "FDR05_left_with_both_FDR_available": int((significant_a & fdr).sum()), "FDR05_right_with_both_FDR_available": int((significant_b & fdr).sum()),
            "FDR05_both": int(both.sum()), "FDR05_both_same_sign": int((both & finite & sign).sum()),
            "FDR05_Jaccard_on_common_available_rows": float(both.sum()/union.sum()) if union.any() else np.nan,
            "unavailable_FDR_retained_as_missing": True, "statistical_implementation_independence": False}


def same_sign_finite(values, reference):
    values = np.asarray(values, float)
    return np.isfinite(values) & np.isfinite(reference) & (np.sign(values) == np.sign(reference)) & (np.sign(reference) != 0)


class Results:
    def __init__(self, e, cfg):
        self.e, self.cfg = e, cfg
        self.root = Path(cfg["results_root"]).resolve(strict=True)
        self.indices, self.members, self.rosters, self.model_roots, self.families, self.atoms = {}, {}, {}, {}, {}, {}
        self.cache, self.endpoint_audits, self.contracts = OrderedDict(), [], []
        for branch in BRANCHES:
            finalized = Path(cfg["deseq2"]["finalized_root_by_branch"][branch])
            e.stage(finalized, "EXECUTED_REVIEW_REQUIRED", "finalize_populations", branch)
            model = self.root / "deseq2" / branch
            e.stage(model, "EXECUTED_REVIEW_REQUIRED", "deseq2", branch)
            self.model_roots[branch] = model
            ix = e.table(finalized / "POPULATION_INDEX.tsv", dtype={k: str for k in SEMANTIC + ["population_key", "membership_sha256"]})
            for c in ["category_family", "analysis_label"]: ix[c] = ix[c].fillna("").astype(str)
            if ix.population_key.duplicated().any() or len({identity(r) for _, r in ix.iterrows()}) != len(ix): raise ValueError("Duplicate population identity")
            mi = e.table(finalized / "POPULATION_MEMBERSHIPS.tsv.gz", dtype=str)
            if mi.duplicated(["population_key", "nucleus_id"]).any() or not set(mi.population_key).issubset(set(ix.population_key)): raise ValueError("Invalid final population memberships")
            member = mi.groupby("population_key", observed=True).nucleus_id.apply(set).to_dict()
            for _, row in ix.iterrows():
                ids = member.setdefault(row.population_key, set())
                if len(ids) != int(row.n_nuclei) or member_sha(ids) != row.membership_sha256: raise ValueError("Population membership/count hash mismatch")
            model_ix = e.table(model / "POPULATION_INDEX.tsv", dtype={k: str for k in SEMANTIC + ["population_key"]})
            for c in ["category_family", "analysis_label"]: model_ix[c] = model_ix[c].fillna("").astype(str)
            if set(map(identity, (r for _, r in ix.iterrows()))) != set(map(identity, (r for _, r in model_ix.iterrows()))): raise ValueError("Model/finalized population identities differ")
            audit = e.table(model / "PSEUDOBULK_EXECUTION_AUDIT.tsv", dtype={"population_key": str})
            if audit.status.eq("ERROR").any(): raise ValueError("Scientific model ERROR present")
            self.indices[branch], self.members[branch] = ix, member
            atoms = e.table(finalized / "POPULATION_ATOMS.tsv", dtype=str)
            if atoms.duplicated(["population_key", "atomic_population_key"]).any(): raise ValueError("Duplicate final population atom")
            self.atoms[branch] = atoms.groupby('population_key').atomic_population_key.apply(list).to_dict()
            self.families[branch] = {identity(r): r for _, r in ix.iterrows()}
            manifest = e.table(cfg["deseq2"]["sample_manifest_tsv_by_branch"][branch], dtype={"sample": str, "condition": str})
            if len(manifest) != 12 or manifest['sample'].duplicated().any() or manifest.condition.value_counts().to_dict() != {"WT": 6, "cKO": 6}: raise ValueError("Incorrect all12 biological sample roster")
            self.rosters[branch] = manifest.set_index("sample").condition.to_dict()
        acceptance = e.stage(self.root / "count_validation", "COMPLETE")
        if acceptance.get("exact_integer_counts_PASS") is not True or acceptance.get("n_differing_entries") != 0: raise ValueError("Independent exact count acceptance is absent/failed")
        checks = e.table(self.root / "count_validation/EXACT_COUNT_CHECKS.tsv", dtype={"population_key": str})
        for branch, ix in self.indices.items():
            actual = checks[checks.branch.eq(branch) & checks.scope.eq("final_union")]
            if set(actual.population_key) != set(ix.population_key) or actual.differing_entries.ne(0).any() or actual.population_key.duplicated().any(): raise ValueError("Independent count acceptance lacks exact final-population coverage")
        self.count_acceptance = acceptance

    def endpoint(self, branch, row, label="PRIMARY"):
        key = (branch, row.tool, row.population_key, label)
        if key in self.cache: self.cache.move_to_end(key); return self.cache[key]
        folder = self.model_roots[branch] / row.tool / row.population_key / label
        audit_path = folder / "AUDIT/EXECUTION_AUDIT.tsv"
        if not audit_path.exists():
            required = label == "PRIMARY" or (branch == 'original_full' and str(row.core_total).lower() in ['true','1'])
            if required: raise ValueError(f"Required completed endpoint audit missing: {audit_path}")
            return None, {"status": "NOT_REQUESTED_OR_MISSING_ENDPOINT", "reason": label, "warning_count": 0}
        ad = self.e.table(audit_path)
        if len(ad) != 1: raise ValueError("Endpoint audit is not one row")
        audit = ad.iloc[0].to_dict()
        self.endpoint_audits.append({"branch": branch, "tool": row.tool, "population": row.population, "population_key": row.population_key, "analysis_kind": row.analysis_kind, "analysis_label": row.analysis_label, "endpoint_folder": label, **audit})
        if audit["status"] != "TESTED_REVIEW_REQUIRED":
            if audit["status"] != "NOT_TESTED": raise ValueError("Unexpected endpoint model status")
            self.cache[key] = (None, audit); return None, audit
        settings = self.e.table(folder / "AUDIT/MODEL_SETTINGS.tsv", dtype=str).set_index("setting").value.to_dict()
        minimum = 1 if row.population == "Inhibitory" else 20
        expect = {"method": "DESeq2", "numerator": "cKO", "denominator": "WT", "biological_unit": "independent_mouse",
                  "normalization": "poscounts", "statistic": "signed_Wald", "min_nuclei_per_mouse": str(minimum),
                  "min_mice_per_genotype": "3", "min_residual_df": "2", "automatic_count_replacement": "disabled",
                  "primary_gene_count_cutoff": "5" if minimum == 1 else "10",
                  "export_fields": ",".join(SEVEN)}
        if any(settings.get(k) != v for k, v in expect.items()): raise ValueError("DESeq2 model policy/settings differ from declared contract")
        expected_counts = "CellBender_integer_count_sensitivity_fixed_original_cluster_union" if branch == 'paired_cellbender' else 'original_CellRanger_integer_UMI_' + branch
        if settings.get('counts') != expected_counts: raise ValueError("DESeq2 count-branch policy binding differs")
        if settings.get("independentFiltering") != ("FALSE" if label.startswith("LOO_omit_") else "TRUE"): raise ValueError("Endpoint filtering policy mismatch")
        if audit["n_WT"] < 3 or audit["n_cKO"] < 3 or audit["residual_df"] < 2: raise ValueError("Tested endpoint violates replicate/design guard")
        if label == "PRIMARY" and str(audit["design"]).replace(" ", "") != "~condition": raise ValueError("Primary design is not condition-only")
        mice = self.e.table(folder / "AUDIT/MICE_ALL12_AUDIT.tsv", dtype={"sample": str, "condition": str})
        if set(mice['sample']) != set(self.rosters[branch]) or mice['sample'].duplicated().any(): raise ValueError("Endpoint all12 mouse audit mismatch")
        if any(mice.condition != mice['sample'].map(self.rosters[branch])): raise ValueError("Endpoint genotype map changed")
        selected = mice[truth(mice.included_in_model)]
        if (selected.n_nuclei < minimum).any() or (selected.total_UMIs <= 0).any(): raise ValueError("Tested endpoint violates mouse count/library policy")
        if selected.condition.eq("WT").sum() != audit["n_WT"] or selected.condition.eq("cKO").sum() != audit["n_cKO"]: raise ValueError("Endpoint animal count audit mismatch")
        z = self.e.table(folder / "AUDIT/ALL_RESULTS_STABLE_ID_AUDIT.tsv", dtype={"gene_id": str})
        display = self.e.table(folder / "ALL_RESULTS.txt")
        if list(display) != SEVEN or not set(["gene_id", *STATS]).issubset(z) or z.gene_id.duplicated().any() or len(z) != len(display): raise ValueError("DESeq2 seven-field/stable-ID export contract mismatch")
        if "export_row" in z and not np.array_equal(z.export_row, np.arange(1, len(z)+1)): raise ValueError("Stable-ID export row binding mismatch")
        if len(z) != int(audit['n_endpoint_genes']): raise ValueError("Endpoint modeled gene count differs from stable-ID results")
        for column in STATS:
            if not np.allclose(z[column], display[column], rtol=0, atol=0, equal_nan=True): raise ValueError("Stable-ID/display numeric export differs")
        if "display_symbol" in z:
            if not np.array_equal(z.display_symbol.fillna("MISSING"), display.gene.fillna("MISSING")): raise ValueError("Symbol/stable-ID export binding mismatch")
        for count, column in [("n_pvalue", "pvalue"), ("n_padj", "padj")]:
            if int(z[column].notna().sum()) != int(audit[count]): raise ValueError("Endpoint finite result availability count mismatch")
        if int(z.padj.lt(.05).sum()) != int(audit["DEG_FDR05"]): raise ValueError("Endpoint FDR05 count differs from audit")
        self.contracts.append({"branch": branch, "tool": row.tool, "population_key": row.population_key, "endpoint": label,
                               "seven_export_fields_exact": True, "stable_ID_binding_exact": True, "n_result_rows": len(z),
                               "n_WT": int(audit["n_WT"]), "n_cKO": int(audit["n_cKO"]), "min_nuclei_policy": minimum,
                               "coverage_flag": audit.get("coverage_flag", ""), "warning_count": int(audit["warning_count"]),
                               "execution_origin": audit.get("execution_origin", "UNKNOWN")})
        self.cache[key] = (z, audit)
        if len(self.cache) > 64: self.cache.popitem(last=False)
        return z, audit

    def universe(self, branch, row):
        p = self.model_roots[branch] / row.tool / row.population_key / "AUDIT/FROZEN_PRIMARY_GENE_PREFILTER_AUDIT.tsv"
        z = self.e.table(p, dtype={"gene_id": str})
        return set(z.loc[truth(z.passes_primary_prefilter), "gene_id"])


def label_reliability(annotations, registry):
    keys = ["tool", "population", "cluster_scope", "cluster_id"]
    reference = registry[registry.branch.eq('original_full') & registry.cluster_scope.eq('original_frozen')]
    z = annotations.merge(reference[keys + ["n_supported_representative_genes", "marker_review_status"]], on=keys, how="left", validate="many_to_one")
    z['representative_support_reference_branch'] = 'original_full'
    z["n_supported_representative_genes"] = z.n_supported_representative_genes.fillna(0).astype(int)
    z["coarse_label_evidence_supported"] = z.atlas_broad_majority_fraction.ge(.8) & z.atlas_broad_high_support_fraction.ge(.8) & ~z.identity_review_status.str.contains("CONFLICT|MIXED_OR_LOW", na=False)
    z["fine_label_evidence_supported"] = (z.coarse_label_evidence_supported & z.atlas_majority_fraction.ge(.8) & z.atlas_subclass_high_support_fraction.ge(.8)
                                          & z.n_supported_representative_genes.ge(2) & z.n_nuclei.ge(100)
                                          & ~z.identity_review_status.str.contains("CONFLICT|MIXED", na=False)
                                          & z.cluster_scope.eq("original_frozen"))
    z["physical_layer_validated"] = False
    z["fine_label_scope"] = np.where(z.fine_label_evidence_supported, "CONSERVATIVE_REFERENCE_SUBCLASS_WITH_MARKER_SUPPORT;ANATOMY_NOT_VALIDATED", "UNCERTAIN_MIXED_RARE_OR_WEAK_MARKER_SUPPORT")
    return z


def population_label_evidence(row, branch, data, label_rows, atomic_lookup):
    if row.analysis_kind == 'total': return {'label_evidence_scope':'DEFINED_TOTAL_COHORT;NO_FINE_IDENTITY_CLAIM','fine_label_evidence_supported':False,'label_uncertainty_present':False}
    atoms = data.atoms[branch].get(row.population_key, [])
    records = [label_rows.get((branch, atomic_lookup[x]['tool'], atomic_lookup[x]['population'], str(atomic_lookup[x]['cluster_id']))) for x in atoms if x in atomic_lookup]
    good = len(records) == len(atoms) and bool(atoms) and all(r is not None and bool(r['fine_label_evidence_supported']) for r in records)
    if row.analysis_kind in ['layer_program','expression_program'] and any(x in str(row.analysis_label).lower() for x in ['mixed','weak','unresolved','not assigned','conflicting']): good = False
    uncertainty = not good or any(r is not None and ('CONFLICT' in str(r['identity_review_status']) or 'MIXED' in str(r['identity_review_status'])) for r in records)
    return {'label_evidence_scope':'ALL_DECLARED_MEMBER_ATOMS_REFERENCE_AND_MARKER_REVIEW;ANATOMY_NOT_VALIDATED',
            'fine_label_evidence_supported':good,'label_uncertainty_present':uncertainty,'n_atomic_clusters':len(atoms)}


def crosswalk(left, right, memberships):
    overlaps = np.asarray([[len(memberships[x] & memberships[y]) for y in right.population_key] for x in left.population_key], dtype=np.int64)
    if not overlaps.size: return []
    rows, cols = linear_sum_assignment(overlaps, maximize=True)
    out = []
    for i, j in zip(rows, cols):
        a, b = left.iloc[i], right.iloc[j]; n = int(overlaps[i, j]); ln, rn = len(memberships[a.population_key]), len(memberships[b.population_key])
        out.append({"left": a, "right": b, "shared_nuclei": n, "left_fraction": n/ln if ln else 0., "right_fraction": n/rn if rn else 0.,
                    "jaccard": n/(ln+rn-n) if ln+rn-n else np.nan, "same_membership": memberships[a.population_key] == memberships[b.population_key],
                    "match_pass": n >= 20 and ln > 0 and rn > 0 and n/ln >= .5 and n/rn >= .5})
    return out


def core_reliability(data, row):
    original, original_audit = data.endpoint("original_full", row)
    if original is None: return pd.DataFrame(), {"tool": row.tool, "population": row.population, "status": original_audit['status'], "reason": original_audit.get('reason', '')}
    base = original[["gene_id", *STATS]].rename(columns={c: "original_full_" + c for c in STATS}).set_index("gene_id")
    audits = {"original_full": original_audit}
    for branch in ["paired_original", "paired_cellbender"]:
        match = data.families[branch].get(identity(row))
        z, ad = data.endpoint(branch, match) if match is not None else (None, {"status": "MISSING_POPULATION"})
        audits[branch] = ad
        if z is not None: base = base.join(z[["gene_id", *STATS]].set_index("gene_id").add_prefix(branch + "_"), how="left")
        else:
            for c in STATS: base[branch + "_" + c] = np.nan
    primary_sign = base.original_full_log2FoldChange.to_numpy(float)
    primary3 = np.ones(len(base), dtype=bool)
    for branch in BRANCHES:
        primary3 &= (base[branch + "_padj"].to_numpy(float) < .05) & same_sign_finite(base[branch + "_log2FoldChange"], primary_sign)
        base[branch + "_endpoint_status"] = audits[branch]["status"]
        base[branch + "_warning_count"] = audits[branch].get("warning_count", np.nan)
    partition_pass = np.ones(len(base), dtype=bool)
    for endpoint in ["KNOWN_PARTITION_ADJUSTED", "MAIN_PARTITION_ONLY"]:
        z, ad = data.endpoint("original_full", row, endpoint)
        prefix = endpoint.lower()
        columns = z.set_index("gene_id") if z is not None else pd.DataFrame(index=base.index)
        for c in ["log2FoldChange", "pvalue", "padj"]: base[prefix + "_" + c] = columns[c].reindex(base.index) if c in columns else np.nan
        base[prefix + "_endpoint_status"] = ad['status']
        base[prefix + "_warning_count"] = ad.get("warning_count", np.nan)
        partition_pass &= same_sign_finite(base[prefix + "_log2FoldChange"], primary_sign) & (ad['status'] == "TESTED_REVIEW_REQUIRED")
    effects, p_available, fdr_available, loo_statuses = [], [], [], []
    for sample in data.rosters["original_full"]:
        z, ad = data.endpoint("original_full", row, "LOO_omit_" + sample)
        z = z.set_index("gene_id").reindex(base.index) if z is not None else pd.DataFrame(index=base.index)
        effects.append(z.log2FoldChange.to_numpy(float) if "log2FoldChange" in z else np.full(len(base), np.nan))
        p_available.append(z.pvalue.notna().to_numpy() if "pvalue" in z else np.zeros(len(base), bool))
        fdr_available.append(z.padj.notna().to_numpy() if "padj" in z else np.zeros(len(base), bool))
        loo_statuses.append(ad['status'])
    matrix = np.vstack(effects).T
    finite = np.isfinite(matrix)
    same = finite & (np.sign(matrix) == np.sign(primary_sign[:, None])) & (np.sign(primary_sign[:, None]) != 0)
    base["n_LOO_requested"] = 12
    base["n_LOO_endpoint_TESTED"] = sum(s == "TESTED_REVIEW_REQUIRED" for s in loo_statuses)
    base["n_LOO_effect_estimable"] = finite.sum(axis=1)
    base["n_LOO_same_sign"] = same.sum(axis=1)
    base["n_LOO_pvalue_available"] = np.vstack(p_available).sum(axis=0)
    base["n_LOO_FDR_available"] = np.vstack(fdr_available).sum(axis=0)
    base["min_LOO_log2FoldChange"] = np.where(finite.any(axis=1), np.min(np.where(finite, matrix, np.inf), axis=1), np.nan)
    base["max_LOO_log2FoldChange"] = np.where(finite.any(axis=1), np.max(np.where(finite, matrix, -np.inf), axis=1), np.nan)
    base["PRIMARY_FDR05_all_three_count_branches_same_sign"] = primary3
    base["both_partition_sensitivities_estimable_same_sign"] = partition_pass
    base["all12_LOO_estimable_same_sign"] = same.sum(axis=1) == 12
    base["stringent_reliability_flag"] = primary3 & partition_pass & (same.sum(axis=1) == 12) & (base.n_LOO_endpoint_TESTED == 12)
    base["reliability_scope"] = "CORE_TOTAL_ONLY;WITHIN_ENDPOINT_BH;SAME_DESEQ2_ENGINE;NO_GLOBAL_DISCOVERY_FDR"
    base["tool"], base["population"], base["population_key"] = row.tool, row.population, row.population_key
    return base.reset_index(), {"tool": row.tool, "population": row.population, "status": "REVIEWED", "n_genes": len(base), "n_stringent": int(base.stringent_reliability_flag.sum()), "n_PRIMARY_FDR05_all3": int(primary3.sum()), "LOO_endpoint_TESTED": sum(s == "TESTED_REVIEW_REQUIRED" for s in loo_statuses)}


def write_table(out, name, frame):
    frame.to_csv(out / name, sep="\t", index=False, na_rep="NA", compression="gzip" if name.endswith(".gz") else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument("--config", required=True, type=Path); ap.add_argument("--outdir", required=True, type=Path); args = ap.parse_args()
    out = args.outdir.resolve()
    if out.exists() and any(out.iterdir()): raise ValueError("Refusing nonempty comparison destination")
    out.mkdir(parents=True, exist_ok=True)
    e = Evidence(); cfg = e.json(args.config)
    if cfg.get("branches") != BRANCHES or cfg.get("tools") != TOOLS or set(cfg.get("populations", [])) != set(POPS): raise ValueError("Incomplete matched comparison configuration")
    data = Results(e, cfg); root = data.root
    annotation_dir = Path(cfg["final_annotation_dir"])
    e.stage(annotation_dir, "ANNOTATION_FINAL_COMPLETE_REVIEW_REQUIRED", "final")
    annotations = e.table(annotation_dir / "cluster_annotations.tsv", dtype={"cluster_id": str})
    summary_dir = Path(cfg["marker_summary_dir"]); e.stage(summary_dir, "COMPLETE", "native_aggregate_summary")
    registry = e.table(summary_dir / "representative_cluster_registry.tsv", dtype={"cluster_id": str})
    reliability = label_reliability(annotations, registry)
    write_table(out, "LABEL_RELIABILITY_REVIEW.tsv", reliability)
    label_rows = {(r.branch, r.tool, r.population, str(r.cluster_id)): r._asdict() for r in reliability[reliability.cluster_scope.eq('original_frozen')].itertuples(index=False)}
    atomic_index = e.table(Path(cfg['initial_annotation_dir'])/'atomic_population_index.tsv', dtype={'population_key':str,'cluster_id':str})
    atomic_lookup = {r.population_key:r._asdict() for r in atomic_index.itertuples(index=False)}
    diagnostic_dir = Path(cfg["ca1_diagnostic_dir"]); diagnostic_receipt = e.stage(diagnostic_dir, "COMPLETE")
    if diagnostic_receipt.get("primary_memberships_unchanged") is not True: raise ValueError("CA1 controls changed primary memberships")
    for name in ["LATENT_SPACE_COMPARISON.tsv", "NEIGHBOR_GRAPH_COMPARISON.tsv", "LEIDEN_FIXED_GRAPH_CONTROLS.tsv", "PYTHON_HARMONY_CONVERGENCE.tsv"]:
        write_table(out, name, e.table(diagnostic_dir / name))
    plot_dir = root / "plots"; e.stage(plot_dir, "COMPLETE", "reviewed_annotation_and_native_expression_rendering")
    plots = e.table(plot_dir / "PLOT_BOOK_INDEX.tsv"); write_table(out, "PLOT_BOOK_INDEX.tsv", plots)
    pairs, matches, skipped, exact = [], [], [], []
    def pair(scope, lb, lr, rb, rr, overlap=None):
        a, aa = data.endpoint(lb, lr); b, ba = data.endpoint(rb, rr)
        base = {"comparison_scope": scope, "population": lr.population, "left_branch": lb, "right_branch": rb,
                "left_tool": lr.tool, "right_tool": rr.tool, "left_population_key": lr.population_key, "right_population_key": rr.population_key,
                "left_analysis_label": lr.analysis_label, "right_analysis_label": rr.analysis_label,
                "left_status": aa['status'], "right_status": ba['status'], "left_warning_count": aa.get('warning_count', 0), "right_warning_count": ba.get('warning_count', 0),
                "left_coverage_flag": aa.get('coverage_flag', ''), "right_coverage_flag": ba.get('coverage_flag', '')}
        for prefix, row, branch in [('left',lr,lb),('right',rr,rb)]:
            base.update({prefix+'_'+k:v for k,v in population_label_evidence(row,branch,data,label_rows,atomic_lookup).items()})
        if overlap: base.update({k: v for k, v in overlap.items() if k not in ["left", "right"]})
        if a is None or b is None: skipped.append({**base, "reason": "MODEL_NOT_TESTED_OR_MISSING"}); return
        base['fine_population_interpretation_guard_pass'] = bool(base['left_fine_label_evidence_supported'] and base['right_fine_label_evidence_supported'] and
                                                                (overlap is None or overlap.get('jaccard',1.) >= .5) and
                                                                not ('LOW_NUCLEUS' in str(base['left_coverage_flag']) or 'LOW_NUCLEUS' in str(base['right_coverage_flag'])))
        pairs.append({**base, **effect_metrics(a, b)})
    for branch in BRANCHES:
        ix, mm = data.indices[branch], data.members[branch]
        for pop in POPS:
            left = ix[(ix.tool == 'scanpy') & (ix.population == pop) & (ix.analysis_kind == 'total')]
            right = ix[(ix.tool == 'seurat') & (ix.population == pop) & (ix.analysis_kind == 'total')]
            if len(left) != 1 or len(right) != 1: raise ValueError("Expected one total/tool/population/branch")
            l, r = left.iloc[0], right.iloc[0]
            if mm[l.population_key] != mm[r.population_key]: raise ValueError("Shared total cross-tool nucleus IDs differ")
            arrays = []
            for row in [l, r]:
                if sha(e.bind(row.counts_file)) != row.counts_sha256: raise ValueError("Exact total count source hash changed")
                z = e.table(row.counts_file, dtype={"gene_id": str})
                samples = list(data.rosters[branch])
                x = z[samples].to_numpy()
                if (x < 0).any() or not np.all(np.isfinite(x)) or not np.all(x == np.floor(x)): raise ValueError("Total array is not nonnegative integer")
                arrays.append((z.gene_id.to_numpy(), x.astype(np.int64)))
            if not np.array_equal(arrays[0][0], arrays[1][0]) or not np.array_equal(arrays[0][1], arrays[1][1]): raise ValueError("Same shared total RNA counts differ across tools")
            exact.append({"branch": branch, "population": pop, "same_nuclei": True, "same_gene_ID_axis": True, "integer_entries_compared": int(arrays[0][1].size), "differing_integer_entries": 0})
            pair("SAME_SHARED_TOTAL_AND_SAME_DESEQ2_ENGINE", branch, l, branch, r)
            for kind, family in ix[(ix.population == pop) & ix.analysis_kind.ne('total')][['analysis_kind','category_family']].drop_duplicates().itertuples(index=False):
                lset = ix[(ix.tool == 'scanpy') & (ix.population == pop) & (ix.analysis_kind == kind) & (ix.category_family == family) & truth(ix.model_enabled)]
                rset = ix[(ix.tool == 'seurat') & (ix.population == pop) & (ix.analysis_kind == kind) & (ix.category_family == family) & truth(ix.model_enabled)]
                assignments = crosswalk(lset, rset, mm)
                chosen_l, chosen_r = set(), set()
                for assignment in assignments:
                    lrow, rrow = assignment['left'], assignment['right']; chosen_l.add(lrow.population_key); chosen_r.add(rrow.population_key)
                    desc = {"branch": branch, "population": pop, "analysis_kind": kind, "category_family": family,
                            "left_population_key": lrow.population_key, "right_population_key": rrow.population_key,
                            "left_label": lrow.analysis_label, "right_label": rrow.analysis_label,
                            **{k:v for k,v in assignment.items() if k not in ['left','right']}}
                    matches.append(desc)
                    if assignment['match_pass']: pair("FROZEN_POPULATION_DEFINITION_COMPARISON_SAME_DESEQ2_ENGINE", branch, lrow, branch, rrow, assignment)
                    else: skipped.append({**desc, "reason": "RECIPROCAL_OVERLAP_LT0.5_OR_COMMON_NUCLEI_LT20"})
                for side, z, chosen in [('left', lset, chosen_l), ('right', rset, chosen_r)]:
                    for _, row in z[~z.population_key.isin(chosen)].iterrows(): skipped.append({"branch": branch, "population": pop, "analysis_kind": kind, "category_family": family, "side": side, "population_key": row.population_key, "reason": "UNMATCHED_FINE_POPULATION"})
    if set(data.families['paired_original']) != set(data.families['paired_cellbender']): raise ValueError("Paired semantic population registry differs")
    for key, l in data.families['paired_original'].items():
        r = data.families['paired_cellbender'][key]
        if data.members['paired_original'][l.population_key] != data.members['paired_cellbender'][r.population_key]: raise ValueError("Paired raw/corrected fixed population membership changed")
        if data.universe('paired_original', l) != data.universe('paired_cellbender', r): raise ValueError("Paired raw/corrected frozen gene universe changed")
        pair("PAIRED_FIXED_MEMBERSHIP_COUNT_CORRECTION_SENSITIVITY", 'paired_original', l, 'paired_cellbender', r, {"same_membership": True, "same_frozen_gene_universe": True})
    core, scopes = [], []
    for _, row in data.indices['original_full'][truth(data.indices['original_full'].core_total)].iterrows():
        z, summary = core_reliability(data, row); core.append(z); scopes.append(summary)
    core_all = pd.concat(core, ignore_index=True)
    write_table(out, 'CORE_WT_cKO_RELIABILITY_ALL_GENES.tsv.gz', core_all)
    write_table(out, 'CORE_STRINGENT_RELIABILITY_GENES.tsv', core_all[core_all.stringent_reliability_flag] if 'stringent_reliability_flag' in core_all else core_all)
    write_table(out, 'CORE_RELIABILITY_SCOPE.tsv', pd.DataFrame(scopes))
    write_table(out, 'EFFECT_COMPARISONS.tsv', pd.DataFrame(pairs)); write_table(out, 'FINE_POPULATION_CROSSWALK.tsv', pd.DataFrame(matches)); write_table(out, 'SKIPPED_COMPARISONS.tsv', pd.DataFrame(skipped))
    write_table(out, 'SHARED_TOTAL_EXACT_COUNTS.tsv', pd.DataFrame(exact)); write_table(out, 'ENDPOINT_EXPORT_CONTRACTS.tsv', pd.DataFrame(data.contracts).drop_duplicates()); write_table(out, 'ENDPOINT_AUDITS_REVIEWED.tsv', pd.DataFrame(data.endpoint_audits).drop_duplicates())
    report = ["Takeaway: Matched snRNA tool comparisons are complete; fine labels and count-sensitive effects require qualified review.", "",
              "This native comparison reviewed the completed 12-mouse WT/cKO revision using fixed original Harmony clusters and shared main analysis cohorts. Exact shared totals agree across tools; both population definitions use the same DESeq2 statistical implementation. Fine population comparisons require reciprocal membership support and retain unmatched, mixed, rare or weak-marker groups. Core reliability requires FDR<0.05 with the same effect direction in all three count branches, both original-full partition sensitivities estimable with the same sign, and all twelve leave-one-mouse-out effects estimable with the same sign. Missing FDR/untested endpoints remain unavailable, and overlapping population tests do not receive a global discovery correction here. Expression-associated CA1/IN program labels remain provisional and do not validate anatomy or morphology.", "",
              f"Reviewed {len(exact)} shared-total count pairs, {len(pairs)} effect comparisons and {len(matches)} fine membership assignments; {len(skipped)} comparisons were skipped explicitly.", "",
              "DESeq2 exports exactly gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj. Display gene symbols retain duplicate/missing symbols with stable-ID row correspondence; comparisons use stable IDs. IN models use every available nucleus in nonempty positive-library mouse aggregates (minimum1), counts>=5 in>=3 eligible mice, >=3 mice/genotype and >=2 residual df; contributions below5 nuclei remain exploratory coverage flags. Other populations use minimum20 nuclei/mouse. Marker identity confidence is separate from these model eligibility rules.", "",
              "Core totals receive original-full PRIMARY, known-partition, main-partition and all12 LOO checks; fine cluster/category models receive PRIMARY only. Stringent core flags are not extended to fine groups without those sensitivity endpoints. The label review separately requires reference majority/high support, representative marker evidence and an absence of mixed/conflicting signatures; coarse marker corroboration alone does not certify a fine subtype.", "",
              "Results: CORE_WT_cKO_RELIABILITY_ALL_GENES.tsv.gz, CORE_RELIABILITY_SCOPE.tsv, EFFECT_COMPARISONS.tsv, FINE_POPULATION_CROSSWALK.tsv and SKIPPED_COMPARISONS.tsv. QC/provenance: LABEL_RELIABILITY_REVIEW.tsv, SHARED_TOTAL_EXACT_COUNTS.tsv, ENDPOINT_EXPORT_CONTRACTS.tsv and ENDPOINT_AUDITS_REVIEWED.tsv. CA1 saved-input technical control tables are retained with their actual metrics; graph/iteration controls do not alter primary memberships, and graph differences alone are not causal attribution.", "",
              "Figures are indexed by PLOT_BOOK_INDEX.tsv with exact native PDF paths. These main frozen-population annotated UMAP and marker books use current shared masks and full-cohort expression sums; secondary fresh views cannot borrow frozen-group dots. WT-only expression programs define memberships before paired count comparisons; Satb2 is excluded from layer definitions. Full source plots/annotation evidence and independent count validation remain separate from biological acceptance."]
    (out / 'comparison_review.md').write_text('\n'.join(report) + '\n')
    e.finish()
    output_hashes = {p.name: sha(p) for p in out.iterdir() if p.is_file()}
    receipt = {"status": "EXECUTED_REVIEW_REQUIRED", "stage": "final_produced_artifact_comparison", "config_sha256": sha(args.config), "source_code_sha256": sha(__file__),
               "input_sha256": e.hashes, "output_sha256": output_hashes, "exact_integer_counts_PASS": True, "same_shared_total_counts_PASS": True,
               "all_required_producer_stages_completed": True, "n_effect_comparisons": len(pairs), "n_fine_assignments": len(matches), "n_skipped_comparisons": len(skipped),
               "scientific_fitting_performed": False, "raw_nucleus_count_aggregation_performed": False, "statistical_implementation_independence": False,
               "physical_layer_annotation_validated": False, "automatic_biological_consistency_verdict": None}
    (out / 'COMPLETE.json').write_text(json.dumps(receipt, indent=2, allow_nan=False) + '\n')
    print(json.dumps({"status": receipt['status'], "n_effect_comparisons": len(pairs), "n_skipped_comparisons": len(skipped)}))


if __name__ == '__main__': main()
