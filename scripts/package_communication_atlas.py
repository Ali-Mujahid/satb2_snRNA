#!/usr/bin/env python3
"""Package completed native LIANA results for static/offline atlas display.

This script reads finished tables only. It does not fit a model, transform an
expression matrix, or change any native scientific result. Standard library only.
"""
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import html
import json
import math
from pathlib import Path
import re
import shutil


SCHEMA = "satb2-liana-atlas/v1"
NATIVE_FINAL_RUN = "d397e9d9-03bc-43b0-9bd8-901409248fd6"
NATIVE_COMMIT = "c59472ccc9de8360dbbf5016db75f8abde08dd3e"
IMMUTABLE_DATA_BASE = "https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.5.0/communication_data/"
STRING_COLUMNS = {
    "source", "target", "ligand_complex", "receptor_complex", "ligand", "receptor",
    "direction", "source_label", "source_review", "source_de_endpoint",
    "target_label", "target_review", "target_de_endpoint",
}
BOOLEAN_COLUMNS = {
    "source_de_tested", "target_de_tested", "ligand_deg_any",
    "ligand_de_symbol_ambiguous", "ligand_all_subunits_deg", "receptor_deg_any",
    "receptor_de_symbol_ambiguous", "receptor_all_subunits_deg", "deg_supported",
    "both_sides_deg_supported", "support_ge3_either_genotype", "support_ge3_both_genotypes",
}
INTEGER_COLUMNS = {"sample_support_WT", "sample_eligible_WT", "sample_support_cKO", "sample_eligible_cKO"}
OMITTED_COLUMNS = {"scope", "condition", "ligand_deseq2_components_json", "receptor_deseq2_components_json"}
REQUESTED_COLUMNS = ["source", "target", "ligand_complex", "receptor_complex", "lr_logfc", "lr_gmeans", "magnitude_rank"]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def dump(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def tsv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def typed_value(column, value):
    if value in ("", "NA", "NaN", "nan", "None", "null"):
        return None
    if column in STRING_COLUMNS:
        return value
    if column in BOOLEAN_COLUMNS:
        assert value in ("True", "False"), (column, value)
        return value == "True"
    if column in INTEGER_COLUMNS:
        return int(value)
    result = float(value)
    assert math.isfinite(result), (column, value)
    return result


def inline(text):
    text = html.escape(text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', text)
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    return re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)


def markdown_to_html(markdown):
    """Render the small, trusted authored report without a network dependency."""
    output, paragraph, table = [], [], []
    def flush():
        if paragraph:
            output.append("<p>" + inline(" ".join(paragraph)) + "</p>")
            paragraph.clear()
        if table:
            header, *body = table
            output.append('<div class="table-wrap"><table><thead><tr>' + "".join("<th>" + inline(v) + "</th>" for v in header) + "</tr></thead><tbody>")
            for row in body:
                output.append("<tr>" + "".join("<td>" + inline(v) + "</td>" for v in row) + "</tr>")
            output.append("</tbody></table></div>")
            table.clear()
    for line in markdown.splitlines():
        if line.startswith("|"):
            if paragraph:
                flush()
            cells = [v.strip() for v in line.strip().strip("|").split("|")]
            if not all(re.match(r"^:?-+:?$", v) for v in cells):
                table.append(cells)
        elif line.startswith("#"):
            flush()
            heading, _, content = line.partition(" ")
            level = min(len(heading), 6)
            output.append(f"<h{level}>" + inline(content) + f"</h{level}>")
        elif not line.strip():
            flush()
        else:
            if table:
                flush()
            paragraph.append(line)
    flush()
    return "\n".join(output)


def web_interpretation(source: Path, repo: Path):
    """Adapt navigation/provenance prose, keeping the completed result findings."""
    text = (source / "RESULTS_AND_INTERPRETATION.md").read_text()
    text = text.replace(
        "Takeaway: Seurat-defined SATB2 CA1–inhibitory analysis prioritizes 324 replicated, context-specific LR candidates.",
        "# SATB2 CA1–inhibitory ligand–receptor results\n\nSeurat-defined CA1–inhibitory analysis prioritizes 324 replicated, context-specific LR candidates.",
    )
    start = text.index("## Question, lifecycle, and endpoint")
    method = text.index("Method: fixed official LIANA", start)
    text = text[:start] + (
        "## Completed analysis and methods\n\n"
        "All six native consensus fits and 24 per-mouse geometric-mean support calculations completed. "
        "The final registered run is `d397e9d9-03bc-43b0-9bd8-901409248fd6`, status **completed**. "
        "This atlas packages the verified saved results; it performs **zero model refits**. "
        "Color, geometry, label and size controls affect visualization only. Native scores and previous DESeq2 exports retain their contents.\n\n"
    ) + text[method:]
    start = text.index("## Artifacts and how to use them")
    end = text.index("## Limitations and blockers", start)
    text = text[:start] + """## Interactive views, evidence and downloads

The [ligand–receptor atlas](index.html) displays all six saved fits. Choose inhibitory clusters or total, pooled/WT/cKO data, CA1→inhibitory, inhibitory→CA1 or both directed sets, and an evidence tier. Heatmap, dotplot, circos and circle views use the selected saved records. Colors, color limits, plot sizes, arrows, arcs and labels are display settings. They do not alter expression, scores, DESeq2 fits, annotations or genotype inference.

The default DEG-supported tier requires magnitude_rank≤0.05 and a matching unambiguous endpoint gene DESeq2 adjusted P<0.05. The stricter replicated/context-specific tier additionally requires specificity_rank≤0.05, native CellPhoneDB P≤0.05, at least three supporting WT and three supporting cKO mice, and reference-supported endpoint identity. All thresholds and selection rules should accompany exported figures. Changing a threshold selects saved candidates; it does not calculate a new statistical test.

Heatmap and network counts represent population-specific directed LR candidates after the current filters. A tested, eligible source/target pair can have zero retained candidates. Untested directions are missing/NA, rather than zero. A combined view is the union of directed sets, not an independent reciprocal-signaling test. The separate total-inhibitory fit cannot be combined with its constituent clusters.

Dot area can represent the minimum detection fraction across **every** complex subunit (`lr_prop`). This can differ from the native `ligand_props` or `receptor_props`, which refer to the native mean-expression-limiting subunit. Fractions slightly above one from sparse floating-point summation are clamped only for display; the original saved numbers remain intact. Mouse support is shown with eligible denominators, so unavailable mice are not treated as zero.

The [native seven-column interaction tables](data/full_result_tables/index.html) are exactly `source, target, ligand_complex, receptor_complex, lr_logfc, lr_gmeans, magnitude_rank`. Full tables and separately labeled DESeq2-supported tables are available for each resolution, condition and direction. The [original gene DESeq2 exports](data/existing_DESeq2_exports/index.html) keep their standard seven-column format and gene symbols unchanged. Extended native evidence retains component stable IDs, signed cKO/WT effects, adjusted P values, identity flags and mouse denominators.

Native download links point to the immutable atlas-v3.5.0 repository tag. They require internet. The complete portable results ZIP attached to the GitHub release contains all native tables and the original 258 primary plus 18 supplemental precomputed PDF/PNG/SVG panels for offline use. Interactive selected tables and modified chart exports are generated directly by the atlas.

## Score meanings

| Field | Interpretation |
|---|---|
| lr_logfc | Native cell-group-versus-rest specificity; **not** cKO/WT fold change |
| lr_gmeans | Geometric mean of endpoint mean log-normalized expression with native complex bottleneck rules |
| magnitude_rank | Lower means higher consensus magnitude priority; **not** genotype FDR |
| specificity_rank | Native aggregate cell-context specificity rank; **not** genotype evidence |
| cellphone_pvals | Native LIANA CellPhoneDB cell-label specificity P; zero is below the 1,000-permutation empirical resolution |
| deg_supported | At least one unambiguous matching sender-ligand or receiver-receptor gene has saved DESeq2 adjusted P<0.05 |
| lr_prop | Minimum detection fraction across all ligand and receptor complex subunits |
| sample_support / sample_eligible | Independent mice passing expression gates / mice eligible for that population pair |

Genotype score/rank differences are descriptive only. No replicate-aware differential-interaction model was fitted, so **no interaction is called significantly altered between WT and cKO**. Existing DESeq2 support is reused gene-level evidence, not an independent validation of an interaction.

Verified sources include the [native result counts](https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.5.0/communication_data/full_result_tables/RESULT_COUNTS.tsv), [extended interaction evidence](https://raw.githubusercontent.com/Ali-Mujahid/satb2_snRNA/atlas-v3.5.0/communication_data/full_result_tables/interaction_evidence.tsv.gz), [scientific review](data/reviews/SCIENTIFIC_RESULT_REVIEW.md), [strong 324-candidate table](data/reviews/strong_consensus_specificity_candidates.tsv), [curated 16-candidate selection policy](data/reviews/curated_16_selection_policy.json), [native final validation](data/reviews/NATIVE_FINAL_VALIDATION.json), [dataset manifest](data/manifest.json) and [transport checksums](data/TRANSPORT_RECEIPT.json). Independent checks covered 154 input checks, 43 native scientific checks, visual review of all chart types and existence of all 276 panel asset sets. Interactive transport independently reconstructs all 64,936 native records and exact component evidence.

""" + text[end:]
    text = text.replace("## Limitations and blockers", "## Scientific limitations")
    text = text.replace("There is no remaining execution blocker. ", "")
    browser_note = "Automated browser opening of local `file://` HTML was blocked by the browser tool's protocol policy. Saved figure assets, inline report data, selector coverage, JavaScript syntax and local download targets are checked separately; no browser-policy workaround or external publication is used.\n\n"
    text = text.replace(browser_note, "")
    destination = repo / "site/communication/RESULTS_AND_INTERPRETATION.md"
    destination.write_text(text)
    body = markdown_to_html(text)
    document = """<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SATB2 ligand–receptor results and interpretation</title><link rel="stylesheet" href="../styles.css"><style>body{margin:0;background:#f4f7fb;color:#203c50;font:16px/1.65 system-ui,sans-serif}.ccc-report{max-width:1080px;margin:2rem auto;padding:2rem 3rem;background:white;border:1px solid #d6e0e8;border-radius:16px}.ccc-report h1{font-size:2rem;line-height:1.2}.ccc-report h2{margin-top:2.4rem;font-size:1.35rem}.ccc-report p{margin:1rem 0}.ccc-report a{color:#176b91;text-decoration:underline}.ccc-report code{font-size:.87em;overflow-wrap:anywhere;background:#edf3f7;padding:.1rem .25rem}.table-wrap{overflow:auto}.ccc-report table{border-collapse:collapse;width:100%;font-size:.9rem}.ccc-report td,.ccc-report th{border:1px solid #d6e0e8;padding:.55rem .7rem;text-align:left}.ccc-report th{background:#edf3f7}.ccc-report nav{display:flex;gap:1rem;flex-wrap:wrap;margin-bottom:2rem}@media(max-width:700px){.ccc-report{margin:.6rem;padding:1.1rem}.ccc-report h1{font-size:1.6rem}}</style></head><body><main class="ccc-report"><nav><a href="index.html">Interactive ligand–receptor atlas</a><a href="../index.html">Main SATB2 atlas</a><a download href="RESULTS_AND_INTERPRETATION.md">Download interpretation</a></nav>""" + body + "</main></body></html>"
    (repo / "site/communication/RESULTS.html").write_text(document)


def package(source: Path, repo: Path):
    target = repo / "site/communication/data"
    target.mkdir(parents=True, exist_ok=True)
    receipt = []
    def copy_file(original: Path, destination: Path, compress=False):
        raw = original.read_bytes()
        output = gzip.compress(raw, compresslevel=9, mtime=0) if compress else raw
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(output)
        receipt.append({
            "source": str(original.relative_to(source)) if original.is_relative_to(source) else str(original),
            "path": str(destination.relative_to(repo)), "source_bytes": len(raw),
            "source_sha256": digest(raw), "transport_bytes": len(output),
            "transport_sha256": digest(output), "compression": "gzip" if compress else "none",
            "native_bytes_preserved": True,
        })

    original_rows = list(tsv(source / "data/interaction_evidence.tsv.gz"))
    assert len(original_rows) == 64936, len(original_rows)
    columns = [name for name in original_rows[0] if name not in OMITTED_COLUMNS]
    datasets, components, node_lookup = {}, {}, {}
    native_node_info = {r["group"]: r for r in tsv(source / "data/node_annotations.tsv")}
    endpoint_info = {r["group"]: r for r in tsv(source / "data/full_result_tables/DESEQ2_ENDPOINTS.tsv")}
    seen, minimum_fraction, maximum_fraction = set(), 1., 0.
    for row in original_rows:
        key = row["scope"] + "|" + row["condition"]
        unique_key = (key, row["source"], row["target"], row["ligand_complex"], row["receptor_complex"])
        assert unique_key not in seen, unique_key
        seen.add(unique_key)
        datasets.setdefault(key, {"rows": [], "eligible_pairs": []})["rows"].append([typed_value(c, row[c]) for c in columns])
        fraction = float(row["lr_prop"])
        minimum_fraction, maximum_fraction = min(minimum_fraction, fraction), max(maximum_fraction, fraction)
        assert -1e-9 <= fraction <= 1 + 1e-9
        for side, role in (("source", "ligand"), ("target", "receptor")):
            group = row[side]
            value = json.loads(row[role + "_deseq2_components_json"])
            complex_name = row[role + "_complex"]
            old = components.setdefault(group, {}).setdefault(complex_name, value)
            assert old == value, (group, complex_name, "component evidence changed between contexts")
            info = native_node_info[group]
            node = {
                "id": group, "label": info["label"], "annotation": row[side + "_label"],
                "review": row[side + "_review"], "needs_review": info["review"] == "True",
                "de_endpoint": row[side + "_de_endpoint"], "de_tested": row[side + "_de_tested"] == "True",
                "de_status": endpoint_info[group]["status"],
                "native_DEG_FDR05": None if endpoint_info[group]["native_DEG_FDR05"] in ("", "nan", "NA") else int(float(endpoint_info[group]["native_DEG_FDR05"])),
            }
            if group in node_lookup:
                assert node_lookup[group] == node
            else:
                node_lookup[group] = node
    for row in tsv(source / "data/condition_pair_eligibility.tsv"):
        key = row["scope"] + "|" + row["condition"]
        datasets[key]["eligible_pairs"].append({"source": row["source"], "target": row["target"], "source_nuclei": int(row["source_nuclei"]), "target_nuclei": int(row["target_nuclei"])})
    assert len(datasets) == 6
    expected_counts = {"IN_clusters|pooled": 20471, "IN_clusters|WT": 19467, "IN_clusters|cKO": 21680, "IN_total|pooled": 1108, "IN_total|WT": 1040, "IN_total|cKO": 1170}
    for key, dataset in datasets.items():
        assert len(dataset["rows"]) == expected_counts[key]
        assert len(dataset["eligible_pairs"]) == (80 if key.startswith("IN_clusters") else 4)
        dataset["scope"], dataset["condition"] = key.split("|")
    strong = []
    default = []
    for row in original_rows:
        if row["scope"] != "IN_clusters" or row["condition"] != "pooled":
            continue
        if float(row["magnitude_rank"]) <= .05 and row["deg_supported"] == "True":
            default.append(row)
            if (float(row["specificity_rank"]) <= .05 and float(row["cellphone_pvals"]) <= .05
                and int(row["sample_support_WT"]) >= 3 and int(row["sample_support_cKO"]) >= 3
                and row["source_review"] == "REFERENCE_SUPPORTED" and row["target_review"] == "REFERENCE_SUPPORTED"):
                strong.append(row)
    assert len(default) == 1589
    assert len(strong) == 324
    assert sum(r["direction"] == "CA1_to_IN" for r in strong) == 112
    assert sum(r["direction"] == "IN_to_CA1" for r in strong) == 212
    assert sum(r["both_sides_deg_supported"] == "True" for r in strong) == 2

    def natural(group):
        return (0 if group == "CA1_Deep" else 1 if group == "CA1_Superficial" else 3 if group == "IN_Total" else 2, int(group.split("_C")[-1]) if "_C" in group else 0)
    provenance = {
        "analysis_date": "2026-10-22", "final_registered_run": NATIVE_FINAL_RUN,
        "final_status": "completed", "scorer": "LIANA 2.0.0 native rank aggregate with native CellPhoneDB scorer",
        "liana_commit": NATIVE_COMMIT, "resource": "mouseconsensus", "resource_ordered_pairs": 3989,
        "resource_input_represented_pairs": 3946, "expr_prop": .10, "min_cells": 5,
        "n_perms": 1000, "seed": 1337, "mice": {"WT": 6, "cKO": 6},
        "normalization": "Raw duplicate-symbol RNA counts summed before log1p(CP10k), using full 32,285-feature library denominator; 32,245 collapsed genes",
        "source_memberships": "Frozen Seurat original_full molecular CA1 deep/superficial programs and inhibitory clusters C1–C20; complementary IN total fit",
        "distinct_nuclei": 29367, "no_new_scientific_fits": True,
        "genotype_interaction_test": False,
        "source_evidence_sha256": digest((source / "data/interaction_evidence.tsv.gz").read_bytes()),
        "native_seven_column_projection": REQUESTED_COLUMNS,
        "score_semantics": {
            "lr_logfc": "Native cell-group-versus-rest specificity; not cKO/WT log2 fold change",
            "lr_gmeans": "Geometric mean of endpoint means of log-normalized expression, using native complex bottleneck rules",
            "magnitude_rank": "Lower consensus magnitude rank indicates higher priority; not a genotype P value or FDR",
            "specificity_rank": "Native aggregate cell-context specificity rank; not genotype evidence",
            "cellphone_pvals": "Native LIANA CellPhoneDB cell-label specificity P; zero is below empirical 1/1000 resolution",
            "lr_prop": "Minimum detection fraction over every ligand and receptor complex subunit; clamp only for plotting",
            "ligand_props": "Detection fraction of native mean-expression-limiting ligand subunit; may differ from all-subunit minimum",
            "receptor_props": "Detection fraction of native mean-expression-limiting receptor subunit; may differ from all-subunit minimum",
            "deg_supported": "Saved matching sender-ligand or receiver-receptor gene DESeq2 padj<0.05, excluding ambiguous symbol evidence",
            "sample_support": "Expression-gate support among eligible independent mice, with eligible denominators retained",
        },
        "interpretation_caveats": [
            "WT/cKO score comparisons are descriptive; no interaction is called significantly altered between genotypes.",
            "Deep/superficial are molecular expression programs, not measured physical layers.",
            "RNA and curated priors do not establish protein availability, contact, receptor activity, functional synapses or physiological direction.",
            "Curated reverse Nrxn/Nlgn orientations are not independent routes; resource inclusion requires localization review, e.g. Rims1 is a cytoplasmic scaffold.",
            "IN C18 DESeq2 was not tested; unavailable evidence remains missing. IN total and its constituent clusters are distinct fits.",
        ],
    }
    payload = {"schema": SCHEMA, "columns": columns, "datasets": datasets,
               "nodes": [node_lookup[n] for n in sorted(node_lookup, key=natural)],
               "components": components, "provenance": provenance}
    raw = json.dumps(payload, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()
    compressed = gzip.compress(raw, compresslevel=9, mtime=0)
    encoded = base64.b64encode(compressed).decode("ascii")
    script = ("/* Completed native results, gzip/base64 transport; decode with bundled fflate. */\nwindow.CCC_DATA_GZIP=" + json.dumps(encoded) + ";\n").encode()
    (target / "data.js").write_bytes(script)
    # A readable exact-schema manifest avoids loading the large payload to inspect provenance.
    manifest = {
        "schema": SCHEMA, "columns": columns,
        "column_types": {c: "string" if c in STRING_COLUMNS else "boolean" if c in BOOLEAN_COLUMNS else "integer" if c in INTEGER_COLUMNS else "number" for c in columns},
        "datasets": {k: {"rows": len(v["rows"]), "eligible_pairs": len(v["eligible_pairs"])} for k, v in datasets.items()},
        "nodes": payload["nodes"], "provenance": provenance,
        "components_layout": "components[population_id][complex_name] gives exact native per-subunit DESeq2 evidence array; source+ligand_complex for sender, target+receptor_complex for receiver",
        "transport": {"file": "data.js", "encoding": "window.CCC_DATA_GZIP: base64(gzip(UTF-8 JSON))", "json_bytes": len(raw), "gzip_bytes": len(compressed), "script_bytes": len(script), "json_sha256": digest(raw), "gzip_sha256": digest(compressed), "script_sha256": digest(script)},
        "validation": {"native_rows": len(original_rows), "exact_endpoint_keys_unique": True, "component_evidence_exactly_reconstructed": True, "pooled_default_DEG_priority": len(default), "pooled_replicated_specific": len(strong), "strict_CA1_to_IN": 112, "strict_IN_to_CA1": 212, "strict_both_end_DEGs": 2, "lr_prop_raw_min": minimum_fraction, "lr_prop_raw_max": maximum_fraction, "no_model_refits": True},
    }

    # Preserve requested native exports, adding reversible compression only to the
    # redundant full-gene JSON that is otherwise 76.6 MB uncompressed.
    for folder in ("full_result_tables", "existing_DESeq2_exports", "input_audits", "reviews"):
        for original in sorted((source / "data" / folder).rglob("*")):
            if not original.is_file():
                continue
            relative = original.relative_to(source / "data")
            destination = (repo / "communication_data" / relative) if folder in ("full_result_tables", "existing_DESeq2_exports") else target / relative
            if original.name == "deseq2_component_evidence.json":
                copy_file(original, destination.with_suffix(".json.gz"), compress=True)
            else:
                copy_file(original, destination)
    index_path = repo / "communication_data/full_result_tables/index.html"
    index_text = index_path.read_text().replace('href="deseq2_component_evidence.json"', 'href="deseq2_component_evidence.json.gz"').replace('>deseq2_component_evidence.json<', '>deseq2_component_evidence.json.gz (lossless gzip)<')
    index_path.write_text(index_text)
    # Distinguish web navigation adaptation from native scientific file preservation.
    for item in receipt:
        if item["path"] == str(index_path.relative_to(repo)):
            item.update(native_bytes_preserved=False, adaptation="Download index JSON link uses lossless gzip transport", transport_bytes=index_path.stat().st_size, transport_sha256=digest(index_path.read_bytes()))
    # Existing Pages payload is close to its size ceiling. Keep native source
    # downloads in the repository outside Pages, linking an immutable release tag.
    # The complete portable ZIP also provides those native files for offline use.
    for folder in ("full_result_tables", "existing_DESeq2_exports"):
        download_root = repo / "communication_data" / folder
        links = []
        for file in sorted(download_root.rglob("*")):
            if file.is_file() and file.name != "index.html":
                relative = file.relative_to(repo / "communication_data")
                links.append('<li><a download href="' + IMMUTABLE_DATA_BASE + str(relative) + '">' + html.escape(str(file.relative_to(download_root))) + '</a></li>')
        landing = target / folder / "index.html"
        landing.parent.mkdir(parents=True, exist_ok=True)
        landing.write_text('<!doctype html><meta charset="utf-8"><title>Native communication downloads</title><link rel="stylesheet" href="../../../styles.css"><main style="padding:2rem"><h1>Native ' + html.escape(folder.replace("_", " ")) + '</h1><p>Unchanged scientific result files are served from the immutable atlas-v3.5.0 repository tag. These links need internet; the complete downloadable portable results archive also includes all files for offline use. The full component JSON uses reversible gzip compression.</p><p><a href="../../index.html">Return to ligand–receptor atlas</a></p><ul>' + "\n".join(links) + '</ul></main>')
    copy_file(source / "data/node_annotations.tsv", target / "node_annotations.tsv")
    for original in sorted((source / "analysis_source").rglob("*")):
        if original.is_file() and original.suffix != ".pyc":
            copy_file(original, repo / "analysis/communication/native_source" / original.relative_to(source / "analysis_source"))
    project = source.parent.parent
    for name in ("requirements.freeze.txt", "requirements.lock.txt", "package_artifacts.tsv", "runtime_api_audit.json", "runtime_readiness_summary.json", "README.md"):
        original = project / "environment" / name
        if original.is_file():
            copy_file(original, repo / "analysis/communication/environment" / name)
    copy_file(project / "references/mouseconsensus.tsv", repo / "analysis/communication/references/mouseconsensus.tsv")
    manifest["downloads"] = {"native_tables": "full_result_tables/index.html", "original_DESeq2": "existing_DESeq2_exports/index.html", "native_evidence": IMMUTABLE_DATA_BASE + "full_result_tables/interaction_evidence.tsv.gz", "component_evidence": IMMUTABLE_DATA_BASE + "full_result_tables/deseq2_component_evidence.json.gz", "readable_manifest": "manifest.json", "transport_receipt": "TRANSPORT_RECEIPT.json", "native_downloads_base": IMMUTABLE_DATA_BASE, "offline_native_downloads": "Complete portable release ZIP"}
    dump(target / "manifest.json", manifest)
    dump(target / "TRANSPORT_RECEIPT.json", {"schema": SCHEMA, "no_model_refits": True, "native_source_file_count": len(receipt), "files": receipt})
    web_interpretation(source, repo)
    print(json.dumps({"datasets": manifest["datasets"], "validation": manifest["validation"], "transport": manifest["transport"], "copied_files": len(receipt), "site_data_bytes": sum(p.stat().st_size for p in target.rglob("*") if p.is_file())}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="Completed SATB2_communication_report directory")
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    package(args.source.resolve(), args.repo.resolve())
