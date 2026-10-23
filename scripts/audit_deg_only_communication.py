#!/usr/bin/env python3
"""Independently audit a strict, endpoint-matched DEG-only LIANA projection.

Reads completed native scores and original DESeq2 exports only. No analysis fit,
expression transformation, native-file mutation, or rank recomputation occurs.
Every ligand AND receptor complex subunit must be an unambiguous significant
gene in its own source/target population. The five endpoints are intentionally
fixed to the user's October 23 correction. Standard library only.
"""
from __future__ import annotations

import argparse
import base64
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil


ENDPOINTS = {
    "CA1_Deep": ("P0244_seurat_CA1_ProS_LAYER_PROGRAM", "Deep-associated expression program enriched"),
    "CA1_Superficial": ("P0247_seurat_CA1_ProS_LAYER_PROGRAM", "Superficial-associated expression program enriched"),
    "IN_C1": ("P0140_seurat_Inhibitory_CLUSTER", "1"),
    "IN_C6": ("P0145_seurat_Inhibitory_CLUSTER", "6"),
    "IN_C8": ("P0147_seurat_Inhibitory_CLUSTER", "8"),
}
CONDITIONS = ("pooled", "WT", "cKO")
DIRECTIONS = ("CA1_to_IN", "IN_to_CA1")
SEVEN = ("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), reader.fieldnames


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def significant(row):
    padj = number(row.get("padj"))
    return padj is not None and padj < .05


def sign(value):
    value = number(value)
    if value is None:
        return "unavailable"
    return "up_in_cKO" if value > 0 else "down_in_cKO" if value < 0 else "zero_effect"


def complex_direction(components):
    signs = {row["direction"] for row in components}
    if len(signs) == 1:
        return next(iter(signs))
    return "mixed_subunit_directions"


def build_projection(repo, project, native_rows, audit, audit_path):
    """Add a small strict transport and downloads, leaving the v3.5 payload intact."""
    target = repo / "site/communication/data"
    old_payload = target / "data.js"
    old_manifest = target / "manifest.json"
    baseline_hashes = {str(path.relative_to(repo)): digest(path) for path in (old_payload, old_manifest)}
    baseline = json.loads(old_manifest.read_text())
    columns = baseline["columns"] + [
        "ligand_deg_direction", "receptor_deg_direction", "ligand_deseq2_log2fc", "receptor_deseq2_log2fc",
        "ligand_deseq2_padj", "receptor_deseq2_padj", "deg_direction_pair",
    ]
    column_types = dict(baseline["column_types"])
    column_types.update({col: "string" if col.endswith("direction") or col == "deg_direction_pair" else "number" for col in columns if col not in column_types})
    def typed(col, value):
        if value in (None, "", "NA", "NaN", "nan", "None", "null"):
            return None
        kind = column_types[col]
        if kind == "string":
            return str(value)
        if kind == "boolean":
            assert value in (True, False, "True", "False")
            return value in (True, "True")
        return int(value) if kind == "integer" else float(value)
    def key(row):
        return tuple(row[col] for col in ("scope", "condition", "source", "target", "ligand_complex", "receptor_complex"))
    strict = {key(row): row for row in audit["strict_saved_records"]}
    assert len(strict) == 47
    selected, components = [], {group: {} for group in ENDPOINTS}
    for native in native_rows:
        evidence = strict.get(key(native))
        if evidence is None:
            continue
        row = {col: typed(col, native[col]) for col in baseline["columns"]}
        row.update(scope=native["scope"], condition=native["condition"])
        signs = []
        for side, role in (("source", "ligand"), ("target", "receptor")):
            parts = evidence[role + "_components"]
            original_components = json.loads(native[role + "_deseq2_components_json"])
            assert all(part["primary_deg_support"] and not part["ambiguous_symbol"] for part in original_components)
            # The historical native component JSON rounded floating-point fields
            # during serialization. Read this correction's signed gene values
            # directly from the byte-verified original stable-ID DESeq2 audits.
            exact_parts = {part["gene_symbol"]: part for part in parts}
            for part in original_components:
                exact = exact_parts[part["symbol"]]
                assert len(part["stable_id_results"]) == 1
                assert part["stable_id_results"][0]["gene_id"] == exact["gene_id"]
                for field in ("log2FoldChange", "padj"):
                    assert math.isclose(part["stable_id_results"][0][field], exact[field], rel_tol=1e-10, abs_tol=1e-13)
                part["stable_id_results"] = [{"gene_id": exact["gene_id"], "log2FoldChange": exact["log2FoldChange"], "padj": exact["padj"]}]
            population, complex_name = native[side], native[role + "_complex"]
            previous = components[population].setdefault(complex_name, original_components)
            assert previous == original_components
            effect = complex_direction(parts)
            compact = {"up_in_cKO": "up", "down_in_cKO": "down", "mixed_subunit_directions": "mixed"}.get(effect, "unavailable")
            signs.append(compact)
            row[role + "_deg_direction"] = compact
            row[role + "_deseq2_log2fc"] = parts[0]["log2FoldChange"] if len(parts) == 1 else None
            row[role + "_deseq2_padj"] = parts[0]["padj"] if len(parts) == 1 else None
        row["deg_direction_pair"] = "→".join(signs)
        assert row["ligand_deseq2_log2fc"] is not None and row["receptor_deseq2_log2fc"] is not None
        assert row["deg_supported"] and row["both_sides_deg_supported"] and row["ligand_all_subunits_deg"] and row["receptor_all_subunits_deg"]
        selected.append(row)
    # The C8 explanation is embedded even though it has no resource-supported edge.
    c8 = audit["endpoint_summaries"]["IN_C8"]["DEGs"]
    assert len(c8) == 1 and c8[0]["gene_symbol"] == "Calb1"
    components["IN_C8"]["Calb1"] = [{"symbol": "Calb1", "tested_endpoint": True, "ambiguous_symbol": False, "significant_any_stable_id": True, "primary_deg_support": True, "in_mouse_resource": False, "stable_id_results": [{"gene_id": c8[0]["gene_id"], "log2FoldChange": c8[0]["log2FoldChange"], "padj": c8[0]["padj"]}]}]
    raw_pairs, _ = read_tsv(repo / "communication_data/full_result_tables/condition_pair_eligibility.tsv")
    datasets = {}
    for condition in CONDITIONS:
        rows = [row for row in selected if row["condition"] == condition]
        pairs = [{"source": row["source"], "target": row["target"], "source_nuclei": int(row["source_nuclei"]), "target_nuclei": int(row["target_nuclei"])} for row in raw_pairs if row["scope"] == "IN_clusters" and row["condition"] == condition and row["source"] in ENDPOINTS and row["target"] in ENDPOINTS]
        assert len(pairs) == 12
        datasets["IN_clusters|" + condition] = {"scope": "IN_clusters", "condition": condition, "rows": [[row[col] for col in columns] for row in rows], "eligible_pairs": pairs}
    by_node = {node["id"]: node for node in baseline["nodes"]}
    nodes = []
    for group in ENDPOINTS:
        node = dict(by_node[group])
        summary = audit["endpoint_summaries"][group]
        node.update(deg_up_in_cKO=summary["up_in_cKO"], deg_down_in_cKO=summary["down_in_cKO"], deg_resource_genes=len(summary["LR_resource_DEG_symbols"]), strict_only=True)
        if group == "IN_C8":
            node["zero_strict_interactions_reason"] = audit["c8_result"]
        nodes.append(node)
    provenance = json.loads(json.dumps(baseline["provenance"]))
    provenance.update(
        correction_date="2026-10-23", atlas_version="3.6.0",
        selection_policy="DEG-only: all ligand and receptor complex subunits must be unambiguous padj<0.05 genes in their exact matching sender/receiver DESeq2 endpoints; only CA1 Deep/Superficial and IN C1/C6/C8",
        source_memberships="Frozen Seurat original_full molecular CA1 deep/superficial programs and inhibitory C1/C6/C8; IN_Total and every other inhibitory endpoint are excluded",
        derivation="Corrected strict projection of completed native pooled/WT/cKO IN_clusters fits. Full-gene expression normalization, scoring background, ranks and all native scientific source files remain unchanged. No refits or expression recomputation.",
        endpoint_scope=list(ENDPOINTS), no_new_scientific_fits=True, genotype_interaction_test=False,
    )
    provenance["score_semantics"].update(
        deg_supported="Every retained ligand and receptor subunit is a matching endpoint DEG (not the legacy OR policy)",
        ligand_deg_direction="Original sender endpoint gene DESeq2 sign: positive log2(cKO/WT)=up; negative=down",
        receptor_deg_direction="Original receiver endpoint gene DESeq2 sign: positive log2(cKO/WT)=up; negative=down",
        deg_direction_pair="Ordered ligand→receptor gene effect directions. This is not a signed interaction activity change.",
    )
    provenance["interpretation_caveats"] = [
        "Genotype effects belong to component genes in the original mouse-level DESeq2 fit, including when viewing WT-only or cKO-only native scores; no genotype interaction test was fitted.",
        "Native lr_logfc remains group-versus-rest specificity, not genotype fold change. Scores and ranks were not refitted or reranked using a DEG-only background.",
        "Deep/superficial are molecular expression programs, not measured physical layers. C8 has an identity-review flag.",
        "Calb1 is C8's sole DEG and absent from mouseconsensus. Zero retained C8 interactions is not proof that C8 does not communicate.",
        "RNA/prior-compatible LR candidates do not establish contact, functional synapses, receptor activation, inhibition or net activity. Reversed Nrxn/Nlgn resource orientations are not independent physiological routes.",
    ]
    payload = {"schema": baseline["schema"], "columns": columns, "datasets": datasets, "nodes": nodes, "components": components, "provenance": provenance}
    raw = json.dumps(payload, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()
    compressed = gzip.compress(raw, compresslevel=9, mtime=0)
    script = ("/* Strict matching-endpoint DEG-only projection of completed native results; no refits. */\nwindow.CCC_DATA_GZIP=" + json.dumps(base64.b64encode(compressed).decode()) + ";\n").encode()
    data_path = target / "deg_only_data.js"
    data_path.write_bytes(script)
    export_root = repo / "communication_data/deg_only_v36"
    export_root.mkdir(parents=True, exist_ok=True)
    files = []
    def write_tsv(path, rows, headers):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({col: "NA" if row.get(col) is None else row.get(col) for col in headers})
        files.append(path)
    seven_lr = provenance["native_seven_column_projection"]
    extended = ["scope", "condition", *columns]
    for condition in CONDITIONS:
        for direction in (*DIRECTIONS, "bidirectional"):
            rows = [row for row in selected if row["condition"] == condition and (direction == "bidirectional" or row["direction"] == direction)]
            rows.sort(key=lambda row: (row["magnitude_rank"], row["source"], row["target"], row["ligand_complex"], row["receptor_complex"]))
            write_tsv(export_root / condition / (direction + "_rank_aggregate.tsv"), rows, seven_lr)
            write_tsv(export_root / condition / (direction + "_signed_extended.tsv"), rows, extended)
    write_tsv(export_root / "interaction_evidence_signed.tsv", selected, extended)
    write_tsv(export_root / "CONDITION_DIRECTION_COUNTS.tsv", [row for row in audit["counts"] if row["scope"] == "IN_clusters"], list(audit["counts"][0]))
    signs = []
    for group, summary in audit["endpoint_summaries"].items():
        origin = repo / "communication_data/existing_DESeq2_exports" / group
        destination = export_root / "DESeq2" / group
        destination.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(origin / "DEG_FDR05.txt", destination / "DEG_FDR05.txt")
        files.append(destination / "DEG_FDR05.txt")
        assert digest(origin / "DEG_FDR05.txt") == digest(destination / "DEG_FDR05.txt")
        degs, _ = read_tsv(origin / "DEG_FDR05.txt")
        stable, stable_columns = read_tsv(origin / "ALL_RESULTS_STABLE_ID_AUDIT.tsv")
        write_tsv(destination / "UP_in_cKO.txt", [row for row in degs if sign(row["log2FoldChange"]) == "up_in_cKO"], list(SEVEN))
        write_tsv(destination / "DOWN_in_cKO.txt", [row for row in degs if sign(row["log2FoldChange"]) == "down_in_cKO"], list(SEVEN))
        write_tsv(destination / "DEG_FDR05_STABLE_ID_AUDIT.tsv", [row for row in stable if significant(row)], stable_columns)
        for row in summary["DEGs"]:
            signs.append({"population": group, "population_key": summary["population_key"], "analysis_label": summary["analysis_label"], **row})
    write_tsv(export_root / "ENDPOINT_GENE_DIRECTIONS.tsv", signs, list(signs[0]))
    pair_counts = []
    for condition in CONDITIONS:
        for direction in (*DIRECTIONS, "bidirectional"):
            rows = [row for row in selected if row["condition"] == condition and (direction == "bidirectional" or row["direction"] == direction)]
            for rank_label, rank_limit in (("all_saved_ranks", 1), ("magnitude_rank_le_0.05", .05)):
                counter = Counter(row["deg_direction_pair"] for row in rows if row["magnitude_rank"] <= rank_limit)
                for direction_pair in ("up→up", "up→down", "down→up", "down→down"):
                    pair_counts.append({"condition": condition, "direction": direction, "rank_selection": rank_label, "deg_direction_pair": direction_pair, "records": counter[direction_pair]})
    write_tsv(export_root / "DIRECTION_PAIR_COUNTS.tsv", pair_counts, list(pair_counts[0]))
    components_path = export_root / "strict_component_evidence.json"
    components_path.write_text(json.dumps(components, indent=2, allow_nan=False) + "\n")
    files.append(components_path)
    base_url = "https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.6.0/communication_data/deg_only_v36/"
    manifest = {
        "schema": "satb2-liana-strict-deg-atlas/v1", "atlas_version": "3.6.0", "policy": audit["policy"], "columns": columns, "column_types": column_types,
        "datasets": {key: {"rows": len(value["rows"]), "eligible_pairs": len(value["eligible_pairs"])} for key, value in datasets.items()},
        "nodes": nodes, "provenance": provenance,
        "endpoint_summaries": {group: {key: value for key, value in summary.items() if key != "DEGs"} for group, summary in audit["endpoint_summaries"].items()},
        "counts": audit["counts"], "direction_pair_counts": pair_counts, "c8_result": audit["c8_result"],
        "components_layout": "components[population_id][complex_name] gives significant matching endpoint subunits only. Gene log2FoldChange/padj are read directly from byte-verified original stable-ID DESeq2 audits, avoiding rounded historical JSON transport values.",
        "transport": {"file": "deg_only_data.js", "encoding": "window.CCC_DATA_GZIP: base64(gzip(UTF-8 JSON))", "json_bytes": len(raw), "gzip_bytes": len(compressed), "script_bytes": len(script), "json_sha256": hashlib.sha256(raw).hexdigest(), "gzip_sha256": hashlib.sha256(compressed).hexdigest(), "script_sha256": hashlib.sha256(script).hexdigest()},
        "validation": {"native_source_rows": len(native_rows), "strict_saved_records": len(selected), "unique_population_LR_keys": 16, "allowed_endpoints_only": True, "all_ligand_AND_receptor_subunits_DEG": True, "C8_explicit_eligible_zero_edge_population": True, "native_scores_and_ranks_unchanged": True, "legacy_data_payload_unchanged": True, "no_model_refits": True},
        "downloads": {"native_downloads_base": base_url, "native_tables": base_url + "index.html", "original_DESeq2": base_url + "index.html#gene-deg", "native_evidence": base_url + "interaction_evidence_signed.tsv", "component_evidence": base_url + "strict_component_evidence.json", "readable_manifest": "DEG_ONLY_MANIFEST.json", "audit": base_url + "DEG_ONLY_ENDPOINT_AUDIT.json"},
        "source_hashes": audit["source_hashes"], "legacy_preserved_hashes": baseline_hashes,
        "export_files": [{"path": str(path.relative_to(export_root)), "bytes": path.stat().st_size, "sha256": digest(path), "url": base_url + str(path.relative_to(export_root))} for path in files],
    }
    manifest_path = target / "DEG_ONLY_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    shutil.copyfile(manifest_path, export_root / "MANIFEST.json")
    shutil.copyfile(audit_path, export_root / "DEG_ONLY_ENDPOINT_AUDIT.json")
    write_results(repo, export_root, audit, selected)
    # The whole download directory is immutable at the tagged Git commit.
    links = []
    for path in sorted(export_root.rglob("*")):
        if path.is_file() and path.name != "index.html":
            relative = str(path.relative_to(export_root))
            links.append('<li><a href="' + base_url + relative + '">' + relative + '</a></li>')
    (export_root / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>Strict DEG-only interaction downloads</title><main style="max-width:75rem;margin:3rem auto;font:17px system-ui;line-height:1.6"><h1>Strict DEG-only interaction downloads · v3.6.0</h1><p>All files contain only matching-endpoint significant ligands and receptors from CA1 Deep/Superficial and IN C1/C6/C8. Native full-universe scores/ranks are preserved, with no refits. Original DEG_FDR05.txt copies retain exactly seven columns and source bytes.</p><p>Interaction tables include all saved native ranks. Apply magnitude_rank≤0.05 explicitly to obtain pooled 8 / WT 9 / cKO 7 records.</p><p id="gene-deg">Positive gene log2(cKO/WT) is up in cKO; negative is down. The ligand→receptor sign pair is not an interaction-activity direction.</p><p><a href="https://ali-mujahid.github.io/satb2_snRNA/communication/">Return to atlas</a></p><ul>' + '\n'.join(links) + '</ul></main>')
    for path, prior_hash in baseline_hashes.items():
        assert digest(repo / path) == prior_hash, (path, "legacy payload changed")
    print(json.dumps({"strict_transport": str(data_path), "transport_rows": len(selected), "transport_bytes": len(script), "manifest": str(manifest_path), "exports": str(export_root), "legacy_payload_unchanged": True}, indent=2))


def write_results(repo, export_root, audit, rows):
    text = """# Corrected DEG-only CA1–inhibitory communication projection

The previous default used a ligand-DEG OR receptor-DEG rule, which retained non-DEG partners and other inhibitory populations. This corrected dataset requires **every ligand AND receptor subunit to be a significant, unambiguous DEG at its own exact endpoint**. Only CA1 Deep, CA1 Superficial, inhibitory C1, C6 and C8 are eligible. Other inhibitory clusters and inhibitory total are excluded.

The corrected atlas uses only the saved DESeq2 `padj < 0.05` genes from the five requested endpoints. It projects completed native LIANA results; it does not refit expression, DESeq2 or LIANA, change expression normalization, or recalculate ranks using a restricted background. Native LIANA scores and ranks retain their original full-gene/background interpretation.

## Exact endpoint DEG counts

| Endpoint | Original DEGs | Up in cKO | Down in cKO | Unambiguous resource DEG genes |
|---|---:|---:|---:|---:|
"""
    for group, summary in audit["endpoint_summaries"].items():
        text += f'| {group} | {summary["native_DEG_FDR05"]} | {summary["up_in_cKO"]} | {summary["down_in_cKO"]} | {len(summary["LR_resource_DEG_symbols"])} |\n'
    text += """
CA1 Deep uses `P0244_seurat_CA1_ProS_LAYER_PROGRAM`, label `Deep-associated expression program enriched`; CA1 Superficial uses `P0247_seurat_CA1_ProS_LAYER_PROGRAM`, label `Superficial-associated expression program enriched`. IN C1/C6/C8 use `P0140`, `P0145`, `P0147` Seurat inhibitory cluster endpoints, respectively. Every copied DEG file is byte-identical to the original. Superficial significant symbols Pakap and Gm35438 are ambiguous across stable IDs and excluded from primary interaction support; neither occurs in the LR resource.

## Strict population-specific interaction counts

| Native scoring condition | CA1→IN, all saved ranks | IN→CA1, all saved ranks | Combined | Combined at magnitude rank≤0.05 |
|---|---:|---:|---:|---:|
| Pooled | 6 | 10 | 16 | 8 |
| WT | 5 | 10 | 15 | 9 |
| cKO | 6 | 10 | 16 | 7 |

These are 47 saved records across the three conditions and 16 distinct population × LR keys. WT/cKO counts reflect native expression qualification and ranking; their difference is descriptive and is not evidence of differential interaction activity. At magnitude rank≤0.05, the directional counts are pooled 3/5, WT 3/6 and cKO 3/4 (CA1→IN / IN→CA1). The stricter replicated/context-specific criteria retain the two C6 Cntn3→CA1 Ptprg candidates in each condition.

### Why C8 has zero strict interactions

C8's sole DEG is **Calb1**, up in cKO (log2FC +1.93206760368337; padj 0.0154178430986269). Calb1 is absent from the mouseconsensus ligand/receptor resource. C8 remains an explicit eligible population with tested directions and zero retained edges. This is not an omission and does not establish an absence of biological communication. Its original mixed-identity review flag remains attached.

## Gene up/down distinction

Positive original gene DESeq2 log2(cKO/WT) means **up in cKO**; negative means **down in cKO**. Both endpoint effects and adjusted P values accompany each interaction. `down→up`, for example, means the sender's ligand DEG is down in cKO and the receiver's receptor DEG is up. It does not assign an up/down activity change to the entire interaction. For a multi-subunit complex with opposing component signs, report mixed subunit directions rather than inventing a single complex effect. Every retained complex here is a single-subunit symbol.

The strongest strict candidates are IN C6 **Cntn3→Ptprg**, directed to both CA1 programs. Cntn3 is down in C6 (log2FC −0.5383457365; padj 0.0037005639), whereas Ptprg is up in Deep (+0.3966628983; padj 0.0164890846) and Superficial (+0.4904369824; padj 2.9270547465×10⁻8). Both endpoints have significant gene effects and all six WT and six cKO mice support the expression gate. The opposite signs still do not establish net receptor activity.

All pooled strict candidates are listed below. Duplicate reverse Nrxn/Nlgn resource orientations do not represent independent physiological pathways.

| Source | Target | Ligand | Receptor | Gene effect pair | Ligand log2(cKO/WT) | Receptor log2(cKO/WT) | Native magnitude rank |
|---|---|---|---|---|---:|---:|---:|
"""
    for row in sorted((row for row in rows if row["condition"] == "pooled"), key=lambda row: row["magnitude_rank"]):
        text += f'| {row["source"]} | {row["target"]} | {row["ligand_complex"]} | {row["receptor_complex"]} | {row["deg_direction_pair"]} | {row["ligand_deseq2_log2fc"]:+.5f} | {row["receptor_deseq2_log2fc"]:+.5f} | {row["magnitude_rank"]:.6g} |\n'
    text += """
## Statistical meanings and limits

- Gene DESeq2 `log2FoldChange` is the signed cKO/WT component effect from six independent mice per genotype.
- Native LIANA `lr_logfc` is a source/target group-versus-rest specificity score, **not** a genotype fold change.
- Native `magnitude_rank` and `specificity_rank` are context priority ranks, not genotype P values or FDRs. Smaller ranks are higher priority.
- The LIANA CellPhoneDB scorer estimates cell-label specificity inside a scoring condition; it is not a standalone CellPhoneDB differential-genotype test.
- No genotype interaction-activity P value was fitted. Transcript expression and curated priors do not establish protein localization, contacts, activation/inhibition, functional synapses or anatomical layers.
- Every original biological source file, the original v3.5 payload, and all prior plots remain preserved. Version 3.6 is an explicitly labeled corrected selection.

## Downloads and reproducibility

The condition/direction `rank_aggregate.tsv` files keep exactly the requested seven native score columns: `source, target, ligand_complex, receptor_complex, lr_logfc, lr_gmeans, magnitude_rank`. Separate signed extended tables add ligand and receptor gene log2FC, adjusted P, up/down and ordered effect pairs. `DESeq2/{population}/DEG_FDR05.txt` keeps the original seven gene columns and bytes; `UP_in_cKO.txt` and `DOWN_in_cKO.txt` are endpoint-matched subsets. Significant stable-ID audit rows are also supplied.

`DEG_ONLY_ENDPOINT_AUDIT.json` contains the independent source-hash and endpoint reconstruction checks. `MANIFEST.json` gives compact transport hashes, source hashes, sign counts, rules and download inventory. The read-only script `scripts/audit_deg_only_communication.py` reproduces this projection using completed native outputs; it performs no scientific fits.
"""
    (repo / "docs/DEG_ONLY_RESULTS.md").write_text(text)
    (export_root / "RESULTS_AND_INTERPRETATION.md").write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--native-project", type=Path, default=Path("/Users/mujahid357/Documents/codex_new/SATB2_CCC_20261022"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    repo, project = args.repo.resolve(), args.native_project.resolve()
    output = args.output or repo / "docs/DEG_ONLY_ENDPOINT_AUDIT.json"
    config = json.loads((project / "config/analysis.json").read_text())
    prepared = project / "runs/20261022_liana_v1/results/prepared"
    feature_rows, _ = read_tsv(prepared / "feature_mapping.tsv")
    symbol_ids = defaultdict(set)
    for row in feature_rows:
        symbol_ids[row["source_gene_symbol"]].add(row["gene_id"])
    ambiguous = {symbol for symbol, ids in symbol_ids.items() if len(ids) > 1}
    ambiguity_rows, _ = read_tsv(prepared / "symbol_ambiguities.tsv")
    assert ambiguous == {row["source_gene_symbol"] for row in ambiguity_rows}
    resource_path = project / "references/mouseconsensus.tsv"
    resource_rows, _ = read_tsv(resource_path)
    resource_symbols = {subunit for row in resource_rows for side in ("ligand", "receptor") for subunit in row[side].split("_")}
    resource_keys = {(row["ligand"], row["receptor"]) for row in resource_rows}
    source_hashes, summaries, deg_maps = [], {}, {}
    for group, (population_key, label) in ENDPOINTS.items():
        endpoint = next(row for row in config["endpoints"] if row["population_key"] == population_key and row["analysis_label"] == label)
        folder = repo / "communication_data/existing_DESeq2_exports" / group
        for filename, config_key in (("ALL_RESULTS.txt", "all_results_path"), ("DEG_FDR05.txt", "deg_path"), ("ALL_RESULTS_STABLE_ID_AUDIT.tsv", "stable_audit_path")):
            local = folder / filename
            native = Path(endpoint[config_key])
            saved_hash = endpoint[config_key + "_sha256"]
            assert digest(local) == saved_hash == digest(native), (group, filename, "source bytes differ")
            source_hashes.append({"endpoint": group, "repo_path": str(local.relative_to(repo)), "native_path": str(native), "sha256": saved_hash})
        native_deg, fields = read_tsv(folder / "DEG_FDR05.txt")
        assert fields == list(SEVEN), (group, fields)
        assert all(significant(row) for row in native_deg)
        stable, _ = read_tsv(folder / "ALL_RESULTS_STABLE_ID_AUDIT.tsv")
        stable_deg = [row for row in stable if significant(row)]
        # Match all native significant rows, including any ambiguous symbols.
        assert len(native_deg) == len(stable_deg) == int(endpoint["DEG_FDR05"])
        projected = Counter(tuple(row["display_symbol"] if col == "gene" else row[col] for col in SEVEN) for row in stable_deg)
        assert projected == Counter(tuple(row[col] for col in SEVEN) for row in native_deg)
        eligible = [row for row in stable_deg if row["source_gene_symbol"] not in ambiguous and row["source_gene_symbol"] and row["source_gene_symbol"] != "NA"]
        assert len({row["source_gene_symbol"] for row in eligible}) == len(eligible)
        deg_maps[group] = {row["source_gene_symbol"]: row for row in eligible}
        rows = []
        for row in stable_deg:
            symbol = row["source_gene_symbol"]
            rows.append({"gene_symbol": symbol, "gene_id": row["gene_id"], "log2FoldChange": number(row["log2FoldChange"]), "padj": number(row["padj"]), "direction": sign(row["log2FoldChange"]), "symbol_ambiguous": symbol in ambiguous, "in_mouse_resource": symbol in resource_symbols})
        summaries[group] = {
            "population_key": population_key, "analysis_label": label,
            "native_DEG_FDR05": len(stable_deg), "up_in_cKO": sum(row["direction"] == "up_in_cKO" for row in rows),
            "down_in_cKO": sum(row["direction"] == "down_in_cKO" for row in rows),
            "unambiguous_DEGs": len(eligible), "ambiguous_significant_symbols": sorted({row["source_gene_symbol"] for row in stable_deg if row["source_gene_symbol"] in ambiguous}),
            "LR_resource_DEG_symbols": sorted(deg_maps[group].keys() & resource_symbols), "DEGs": rows,
        }
    table = repo / "communication_data/full_result_tables/interaction_evidence.tsv.gz"
    native_table = project / "runs/20261022_liana_v1/results/communication/tables/interaction_evidence.tsv.gz"
    assert digest(table) == digest(native_table)
    prior_digest = digest(table)
    native_rows, _ = read_tsv(table)
    assert len(native_rows) == 64936
    kept, legacy_counts, new_counts, endpoint_pair_counts = [], Counter(), Counter(), Counter()
    def components(group, complex_name):
        result = []
        for subunit in complex_name.split("_"):
            row = deg_maps[group].get(subunit)
            if row is None:
                return None
            result.append({"gene_symbol": subunit, "gene_id": row["gene_id"], "log2FoldChange": number(row["log2FoldChange"]), "padj": number(row["padj"]), "direction": sign(row["log2FoldChange"])})
        return result
    for row in native_rows:
        key = (row["scope"], row["condition"], row["direction"])
        if row["deg_supported"] == "True" and float(row["magnitude_rank"]) <= .05:
            legacy_counts[key] += 1
        if row["scope"] != "IN_clusters" or row["source"] not in deg_maps or row["target"] not in deg_maps:
            continue
        ligand = components(row["source"], row["ligand_complex"])
        receptor = components(row["target"], row["receptor_complex"])
        if not ligand or not receptor:
            continue
        assert (row["ligand_complex"], row["receptor_complex"]) in resource_keys
        assert row["ligand_all_subunits_deg"] == row["receptor_all_subunits_deg"] == "True"
        new_counts[key + ("all_saved_ranks",)] += 1
        endpoint_pair_counts[row["condition"], row["source"], row["target"]] += 1
        if float(row["magnitude_rank"]) <= .05:
            new_counts[key + ("magnitude_rank_le_0.05",)] += 1
        strong = (float(row["magnitude_rank"]) <= .05 and float(row["specificity_rank"]) <= .05 and float(row["cellphone_pvals"]) <= .05 and int(row["sample_support_WT"]) >= 3 and int(row["sample_support_cKO"]) >= 3 and all(row[side + "_review"] == "REFERENCE_SUPPORTED" for side in ("source", "target") if row[side].startswith("IN_")))
        if strong:
            new_counts[key + ("replicated_specific",)] += 1
        ligand_sign, receptor_sign = complex_direction(ligand), complex_direction(receptor)
        pair_sign = "both_up_in_cKO" if ligand_sign == receptor_sign == "up_in_cKO" else "both_down_in_cKO" if ligand_sign == receptor_sign == "down_in_cKO" else "opposite_directions" if {ligand_sign, receptor_sign} == {"up_in_cKO", "down_in_cKO"} else "mixed_subunit_directions"
        kept.append({
            "scope": row["scope"], "condition": row["condition"], "direction": row["direction"], "source": row["source"], "target": row["target"],
            "ligand_complex": row["ligand_complex"], "receptor_complex": row["receptor_complex"],
            "ligand_direction": ligand_sign, "receptor_direction": receptor_sign, "component_direction_pattern": pair_sign,
            "ligand_components": ligand, "receptor_components": receptor,
            **{col: number(row[col]) for col in ("lr_logfc", "lr_gmeans", "magnitude_rank", "specificity_rank", "cellphone_pvals", "lr_prop")},
            "passes_magnitude_rank_0.05": float(row["magnitude_rank"]) <= .05, "passes_replicated_specific": strong,
            **{col: int(row[col]) for col in ("sample_support_WT", "sample_support_cKO", "sample_eligible_WT", "sample_eligible_cKO")},
        })
    counts = []
    for scope in ("IN_clusters", "IN_total"):
        for condition in CONDITIONS:
            for direction in DIRECTIONS:
                key = (scope, condition, direction)
                counts.append({"scope": scope, "condition": condition, "direction": direction, "legacy_any_endpoint_DEG_rank_0.05": legacy_counts[key], **{metric: new_counts[key + (metric,)] for metric in ("all_saved_ranks", "magnitude_rank_le_0.05", "replicated_specific")}, "endpoint_scope_status": "eligible_requested_population_pairs" if scope == "IN_clusters" else "excluded_not_one_of_requested_five_endpoints"})
    assert len(kept) == 47
    expected = {"pooled": (6, 10, 3, 5), "WT": (5, 10, 3, 6), "cKO": (6, 10, 3, 4)}
    for condition, values in expected.items():
        actual = tuple(new_counts["IN_clusters", condition, direction, metric] for metric in ("all_saved_ranks", "magnitude_rank_le_0.05") for direction in DIRECTIONS)
        assert actual == values, (condition, actual, values)
    assert summaries["IN_C8"]["LR_resource_DEG_symbols"] == []
    assert digest(table) == prior_digest == digest(native_table)
    result = {
        "schema": "satb2-strict-deg-only-independent-audit/v1", "status": "PASS", "created_date": "2026-10-23",
        "policy": {"DESeq2_FDR": "padj < 0.05", "effect_direction": "signed original gene DESeq2 log2(cKO / WT)", "endpoint_scope": list(ENDPOINTS), "complex_rule": "Every ligand AND every receptor complex subunit is an unambiguous DEG at its exact matching sender or receiver endpoint", "native_rank_rule": "Original complete-universe LIANA scores and ranks are retained; no rank recomputation after DEG selection", "total_inhibitory_fit": "Excluded because IN_Total is not one of the five requested endpoints", "pseudobulk_or_LIANA_refits": 0},
        "original_issue": "The previous default used a logical OR: any matching ligand DEG OR receptor DEG. It admitted interactions with a non-DEG partner and populations outside C1/C6/C8. The new strict projection requires a logical AND across all complex subunits at both exact endpoints.",
        "sign_interpretation": {"up_in_cKO": "positive gene DESeq2 log2(cKO/WT)", "down_in_cKO": "negative gene DESeq2 log2(cKO/WT)", "complex_mixed": "Subunits have different signs; retain component-level signs without assigning one complex direction", "opposite_directions": "Ligand and receptor gene effects differ in sign; this is not a signed interaction-activity change", "native_lr_logfc": "Native source/target group-versus-rest specificity; NOT genotype fold change", "genotype_interaction_test": "Not fitted. Genotype interaction significance and up/down activity are not inferred."},
        "endpoint_summaries": summaries, "counts": counts,
        "all_rank_unique_population_interaction_keys": len({(row["source"], row["target"], row["ligand_complex"], row["receptor_complex"]) for row in kept}),
        "strict_saved_records": kept,
        "c8_result": "Calb1 is the sole saved C8 DEG (up in cKO); it is absent from mouseconsensus ligand/receptor subunits, so C8 has zero strict resource-supported interactions. This is not an omitted analysis or proof of no communication.",
        "source_hashes": source_hashes + [{"repo_path": str(table.relative_to(repo)), "native_path": str(native_table), "sha256": prior_digest}, {"native_path": str(resource_path), "sha256": digest(resource_path)}, {"native_path": str(prepared / "feature_mapping.tsv"), "sha256": digest(prepared / "feature_mapping.tsv")}],
        "checks": {"original_DESeq2_exports_byte_identical": True, "standard_seven_column_DEGs_match_stable_ID_audit": True, "resource_keys_match": True, "every_retained_complex_subunit_significant_matching_endpoint": True, "native_saved_both_all_subunits_flags_match_independent_audit": True, "native_scores_untouched": True, "no_input_mutations": True},
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    build_projection(repo, project, native_rows, result, output)
    print(json.dumps({"status": "PASS", "report": str(output), "strict_records": len(kept), "unique_keys": result["all_rank_unique_population_interaction_keys"], "counts": counts}, indent=2))


if __name__ == "__main__":
    main()
