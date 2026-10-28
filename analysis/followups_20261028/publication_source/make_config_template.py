#!/usr/bin/env python3
"""Write a DISABLED transport template; never package data or execute a workflow."""
import argparse
import hashlib
import json
from pathlib import Path


def rec(path):
    p = Path(path).resolve(strict=False)
    if not p.is_file():
        return {"path": str(p), "bytes": None, "sha256": None}
    h = hashlib.sha256()
    with p.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return {"path": str(p), "bytes": p.stat().st_size, "sha256": h.hexdigest()}


def template(project):
    f = Path(project).resolve(strict=True)
    base = f.parent
    strict = Path("/Volumes/Neurowissenschaften/-Ali_on_server/snr_p60/SATB2_FOLLOWUPS_20261028_sccoda_strict_v1_1_0")
    runs = {
        "initial": f / "runs/20261028_sccoda_v1_0_2",
        "strict": strict,
        "cellchat": f / "runs/20261028_cellchat_v1",
        "layout": f / "runs/20261028_cellchat_layout_v1",
        "monocle3": f / "runs/monocle3_20261028_v1",
    }
    science = {k: p / "results" for k, p in runs.items()}
    science["cellchat"] /= "cellchat"
    science["layout"] /= "layout"
    prepared = base / "SATB2_REGULATORY_20261026/runs/20261026_hdwgcna_v1_0_2_export_recovery/results/prepared"
    reviews = {
        "initial": f / "reviews/20261028_composition_v1_0_2_completed_review.md",
        "strict": f / "reviews/20261028_composition_strict_completed_review.md",  # root supplies its actual final review
        "cellchat": f / "reviews/20261028_cellchat_v1_completed_review.md",
        "layout": f / "reviews/20261028_cellchat_layout_completed_review.md",
        "monocle3": f / "reviews/20261028_monocle3_completed_review.md",
    }
    layout_display = f / "reviews/CellChat_layout_final_audit_20261028/DISPLAY_ACCEPTANCE.json"

    def gate(key, kind=None):
        g = {"kind": kind or key, "root": str(science[key]), "archive_prefix": f"native/{key}", "seal": rec(science[key] / "FINAL_VALIDATION.json"), "review_receipts": [rec(reviews[key])]}
        if key == "layout":
            g["display_review_receipts"] = [rec(layout_display)]
        return g

    def native_component(ident, keys, role, text):
        return {"id": ident, "archive_name": f"SATB2_{ident}_20261028_v1.zip", "role": role, "readme": text + "\n\nNative byte hashes and statuses are in ARCHIVE_MANIFEST.json. This transport archive is not additional scientific evidence. Original DESeq2 fits and prior analyses remain unchanged.\n", "gates": [gate(k, "composition" if k in ("initial", "strict") else k) for k in keys],
                "trees": [{"path": str(runs[k] / name), "archive_prefix": f"executed/{k}/{name}"} for k in keys for name in ("workflow", "config", "logs")]}

    components = [
        native_component("composition_initial", ["initial"], "initial_native_composition_inconclusive", "# Initial scCODA composition\n\nAll sixteen models completed computationally, but every original frozen sampling gate failed. Three native selections were withheld and none was validated. Use raw mouse-level proportions descriptively; no-difference conclusions are unsupported. Keep this archive separate from the stricter same-model check."),
        native_component("composition_strict", ["strict"], "same_model_stricter_sampler_native", "# Separate same-model scCODA check\n\nThis archive contains the actual reviewed terminal outputs of the single globally specified target-acceptance0.99/warmup4000/retained6000 check. Data, references, priors, contrasts, seeds, depth and inference gates were preserved. Consult each RESULT_STATUS and the actual review; computational completion does not imply validated proportion differences. There is no further automatic tuning."),
        native_component("cellchat", ["cellchat", "layout"], "strict_DEG_native_cellchat_and_reviewed_layout", "# Strict-DEG CellChat with reviewed readable layouts\n\nEvery retained ligand/receptor subunit is an original DEG in its matching source/target endpoint. Signs, within-mouse native scores, separate exploratory mouse-label tests and both averaging methods are preserved. The combined view is not another independent test. Original plots have documented clipping/vertical-jitter artifacts; use the separate mapped readable layouts. RNA predictions do not establish physical signaling or pathway activation. Full-cohort primary leads were not primary common-partition discoveries."),
        native_component("monocle3", ["monocle3"], "native_exploratory_state_graphs", "# Monocle3 expression-state graphs\n\nCA1 DL/SL and inhibitory C1/C6/C8 preserve their distinct native models, all twelve mice, graphs, gene modules, raw/display aggregates and complete native object/index directories. Arbitrary-reference graph distance is not established developmental time. Native graph association is not independent-mouse genotype gene-DE. Preserve original DESeq2 for genotype gene inference; CA1 anatomical and C8 identity uncertainty remain."),
    ]
    comp_cfg = json.loads((f / "config/composition.v1_0_2.json").read_text())
    cc_cfg = json.loads((f / "config/cellchat.v1.json").read_text())
    endpoint_config = Path(cc_cfg["existing_endpoint_config"]["path"])
    endpoints = json.loads(endpoint_config.read_text())["existing_deseq2"]
    files = []
    external = []
    for name in ("base_metadata", "nucleus_annotations", "cluster_annotations"):
        r = comp_cfg["inputs"][name]
        files.append({"member": f"inputs/composition/{name}/{Path(r['path']).name}", "receipt": {**r, "bytes": Path(r["path"]).stat().st_size}})
        external.append(r["path"])
    files.append({"member": "inputs/endpoint_selection/hdwgcna.v1_0_2_export_recovery.json", "receipt": rec(endpoint_config)})
    external.append(str(endpoint_config))
    for endpoint, r in endpoints.items():
        files.append({"member": f"inputs/original_DESeq2/{endpoint}/ALL_RESULTS.txt", "receipt": {k: r[k] for k in ("path", "bytes", "sha256")}})
        external.append(r["path"])
    for p in sorted((f / "environment").iterdir()):
        if p.is_file():
            files.append({"member": "environment_provenance/" + p.name, "receipt": rec(p)})
    software = json.loads((f / "environment/SOFTWARE_DISTRIBUTION_RECEIPT.json").read_text())
    for name, r in software["native_sources"].items():
        candidates = list((f / "environment/R-distributions").glob(f"{name}_*.tar.gz"))
        if len(candidates) != 1:
            raise ValueError(f"Expected one recorded native source archive: {name}")
        files.append({"member": "official_native_source/" + candidates[0].name, "receipt": {"path": str(candidates[0]), "bytes": r["bytes"], "sha256": r["sha256"]}})
    for p in sorted((f / "receipts").glob("*.json")):
        files.append({"member": "known_native_history/" + p.name, "receipt": rec(p)})
    for name in ("20261028_composition_v1_failed_review.md", "20261028_composition_v1_0_1_failed_review.md", "20261028_scvelo_blocked_review.md"):
        files.append({"member": "reviews/" + name, "receipt": rec(f / "reviews" / name)})
    for d in ("CellChat_actual_readonly", "monocle3_readonly_audit", "CellChat_layout_final_audit_20261028"):
        for p in sorted((f / "reviews" / d).glob("*")):
            if p.is_file() and p.suffix.lower() in (".md", ".json", ".tsv", ".py", ".r", ".log"):
                files.append({"member": f"independent_audits/{d}/{p.name}", "receipt": rec(p)})
    for p in sorted((f / "reviews/composition_readonly_audit/20261028_terminal_v1_0_2").glob("*")):
        if p.is_file():
            files.append({"member": "independent_audits/composition_initial/" + p.name, "receipt": rec(p)})
    shared = {"id": "shared_inputs_source", "archive_name": "SATB2_shared_inputs_source_20261028_v1.zip", "role": "shared_frozen_inputs_and_source", "deduplicate_identical_files": True,
              "readme": (Path(__file__).with_name("README.md")).read_text(),
              "gates": [{"kind": "prepared", "root": str(prepared), "archive_prefix": "inputs/prepared_CA1_IN", "seal": rec(prepared / "COMPLETE.json"), "review_receipts": [rec(reviews["cellchat"]), rec(reviews["monocle3"])]}],
              "files": files,
              "trees": [{"path": str(f / path), "archive_prefix": "source/" + path} for path in ("composition_source_v1_0_2", "composition_sampler_strict_source_v1_1_0", "monocle3_source_release_v1", "presentation_source", "publication_source", "assessment", "velocity_readiness")],
              "scope_note": "Only already frozen source/input/history files; the running strict sampler's result tree is excluded. Its final source/config/logs/results are included in its own later native archive."}
    components.append(shared)
    reader = f / "deliverables/followups_offline_reader_v1"
    components.append({"id": "portable_reader", "archive_name": "SATB2_followups_portable_offline_20261028_v1.zip", "role": "reviewed_portable_offline_reader", "readme": "# Portable offline reader\n\nOpen offline_reader/index.html and keep the full directory together. All packaged plot/PDF/table links are relative; full native model/source archives are separate. No scientific analysis occurs in the browser. Four tasks stay separate, with velocity explicitly unavailable.\n",
                       "gates": [{"kind": "portable", "root": str(reader), "archive_prefix": "offline_reader", "seal": rec(reader / "READER_MANIFEST.json"), "review_receipts": [rec(f / "reviews/20261028_followups_reader_completed_review.md")], "display_review_receipts": [rec(f / "reviews/20261028_followups_reader_rendered_review.json")]}]})
    return {"schema": "satb2-followups-archive-delivery/v1", "approved_for_archives": False,
            "template_note": "DISABLED. Root must supply actual strict terminal/review and reader paths/receipts, verify/freeze component definitions and explicitly approve a copied configuration. This file is not publication approval.",
            "destination_dir": "/Volumes/Neurowissenschaften/-Ali_on_server/snr_p60/SATB2_FOLLOWUPS_20261028_delivery_v1",
            "allowed_output_roots": ["/Volumes/Neurowissenschaften/-Ali_on_server/snr_p60"],
            "allowed_source_roots": [str(f), str(strict), str(prepared), *external],
            "forbidden_source_roots": [str(f / "environment/sccoda"), str(f / "environment/R-library")],
            "maximum_archive_bytes": 1_900_000_000, "components": components}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True, help="New disabled JSON template path; must not exist")
    args = parser.parse_args()
    if args.destination.exists():
        raise ValueError("Never overwrite an existing template/configuration")
    args.destination.write_text(json.dumps(template(args.project_root), indent=2) + "\n")
    print(args.destination)
