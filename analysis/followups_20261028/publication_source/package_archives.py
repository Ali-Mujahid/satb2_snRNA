#!/usr/bin/env python3
"""Read-only transport of reviewed native files into bounded deterministic ZIPs.

This program does not perform a scientific analysis. It verifies existing native
seals and root-supplied review receipts, then copies their exact bytes. It never
executes archived code or creates/changes files in a source or scientific run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import stat
import zipfile
from pathlib import Path, PurePosixPath

CHUNK = 1024 * 1024
ZIP_DATE = (2026, 10, 28, 0, 0, 0)
CAP = 1_900_000_000
NATIVE = {
    "composition": ("satb2-composition-final-validation/v1", {"COMPLETE", "COMPLETE_WITH_INCONCLUSIVE_OR_UNTESTABLE_SCOPES"}, "outputs"),
    "cellchat": ("satb2-cellchat-final-validation/v1", {"COMPLETE"}, "files"),
    "layout": ("satb2-cellchat-layout-only-final/v1", {"COMPLETE_PENDING_VISUAL_REVIEW"}, "outputs"),
    "monocle3": ("satb2-monocle3-final-validation/v1", {"COMPLETE"}, "outputs"),
    "prepared": ("satb2-regulatory-preparation/v1", {"COMPLETE"}, "files"),
    "portable": ("satb2-followups-portable-reader-manifest/v1", {"BUILT_PENDING_RENDERED_REVIEW"}, "files"),
}
SKIP = {"__pycache__", ".git", ".snakemake", "node_modules", ".cache"}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def record(path: Path) -> dict:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)}


def safe_member(value: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value or ":" in value:
        raise ValueError(f"Unsafe archive member: {value!r}")
    if any(ord(c) < 32 for c in value) or value.startswith("/"):
        raise ValueError(f"Unsafe archive member: {value!r}")
    if any(part in ("", ".", "..") for part in value.split("/")):
        raise ValueError(f"Unsafe archive member: {value!r}")
    if any(part != part.rstrip(" .") for part in value.split("/")):
        raise ValueError(f"Nonportable archive member: {value!r}")
    if PurePosixPath(value).is_absolute():
        raise ValueError(f"Absolute archive member: {value!r}")
    return value


def canonical(path, allowed: list[Path], *, directory=False) -> Path:
    p = Path(path)
    if not p.is_absolute():
        raise ValueError(f"Source must be absolute: {p}")
    # Reject links even when they happen to point inside an allowed root.
    for q in [p, *p.parents]:
        if q.is_symlink():
            raise ValueError(f"Symlink source is not allowed: {p}")
    p = p.resolve(strict=True)
    if not any(p.is_relative_to(root) for root in allowed):
        raise ValueError(f"Source outside explicitly allowed roots: {p}")
    if directory and not p.is_dir() or not directory and not p.is_file():
        raise ValueError(f"Wrong source type: {p}")
    return p


def verify(rec: dict, allowed: list[Path], base: Path | None = None) -> Path:
    p = Path(rec["path"])
    if not p.is_absolute():
        if base is None:
            raise ValueError("Relative source receipt has no native root")
        safe_member(rec["path"])
        p = base / p
    p = canonical(p, allowed)
    if base is not None and not p.is_relative_to(base):
        raise ValueError(f"Native file escapes its sealed root: {p}")
    size = rec.get("bytes", rec.get("size_bytes"))
    if size is None or p.stat().st_size != size or sha(p) != rec["sha256"]:
        raise ValueError(f"File receipt mismatch: {p}")
    return p


def join_member(prefix: str, suffix: str) -> str:
    return safe_member(prefix + "/" + suffix)


def gate_files(gate: dict, allowed: list[Path]):
    kind = gate["kind"]
    if kind not in NATIVE:
        raise ValueError(f"Unknown native seal kind: {kind}")
    root = canonical(gate["root"], allowed, directory=True)
    if "synthetic" in str(root).lower():
        raise ValueError("Synthetic result trees cannot become study archives")
    seal_path = verify(gate["seal"], allowed, root)
    seal = json.loads(seal_path.read_text())
    schema, states, key = NATIVE[kind]
    if seal.get("schema") != schema or seal.get("status") not in states:
        raise ValueError(f"Native component is not terminal: {seal_path}")
    reviews = gate.get("review_receipts", [])
    if not reviews:
        raise ValueError(f"Actual root review receipt required: {seal_path}")
    if kind in ("layout", "portable") and not gate.get("display_review_receipts"):
        raise ValueError("Actual rendered/display review required for layout/reader")
    if kind == "layout":
        accepted = False
        for review in gate["display_review_receipts"]:
            p = verify(review, allowed)
            if p.suffix.lower() != ".json":
                continue
            proof = json.loads(p.read_text())
            accepted |= (proof.get("schema") == "satb2-cellchat-layout-display-acceptance/v1"
                         and proof.get("approved_for_display") is True
                         and proof.get("layout_final_validation", {}).get("sha256") == gate["seal"]["sha256"])
        if not accepted:
            raise ValueError("Layout display acceptance must bind this exact native repair seal")
    out = seal.get(key, [])
    if isinstance(out, dict):
        out = [{"path": k, **v} for k, v in out.items()]
    if not out:
        raise ValueError(f"Native seal has no file receipts: {seal_path}")
    paths = [(p, str(p.relative_to(root)), rec) for rec in out for p in [verify(rec, allowed, root)]]
    if len({rel for _, rel, _ in paths}) != len(paths):
        raise ValueError("Duplicate native output receipt")
    paths.append((seal_path, str(seal_path.relative_to(root)), record(seal_path)))
    return root, seal, paths


def make_plan(configuration: Path, component_ids=None, *, require_fresh_destination=True) -> dict:
    cfg = json.loads(configuration.read_text())
    if cfg.get("schema") != "satb2-followups-archive-delivery/v1" or cfg.get("approved_for_archives") is not True:
        raise ValueError("A root-approved archive-delivery configuration is required")
    allowed = [Path(p).resolve(strict=True) for p in cfg["allowed_source_roots"]]
    forbidden = [Path(p).resolve(strict=False) for p in cfg.get("forbidden_source_roots", [])]
    output_allowed = [Path(p).resolve(strict=True) for p in cfg["allowed_output_roots"]]
    dest = Path(cfg["destination_dir"])
    if not dest.is_absolute() or require_fresh_destination and dest.exists():
        raise ValueError("Archive destination must be a fresh absolute directory")
    parent = dest.parent.resolve(strict=True)
    if not any(parent.is_relative_to(p) for p in output_allowed):
        raise ValueError("Archive destination is outside the approved output roots")
    if not os.access(parent, os.W_OK):
        raise ValueError("Archive destination parent is not writable")
    for p in [dest, *dest.parents]:
        if p.is_symlink():
            raise ValueError("Archive destination must not use symlinks")
    dest = parent / dest.name
    cap = int(cfg.get("maximum_archive_bytes", CAP))
    if not 0 < cap <= CAP:
        raise ValueError("Archive limit must be positive and at most 1.9 GB")
    components = []
    archive_names = set()
    configured_ids = [c["id"] for c in cfg["components"]]
    if len(set(configured_ids)) != len(configured_ids) or any(not re.fullmatch(r"[A-Za-z0-9_-]+", c) for c in configured_ids):
        raise ValueError("Component IDs must be unique safe labels")
    selected = set(configured_ids if component_ids is None else component_ids)
    if not selected or not selected.issubset(set(configured_ids)):
        raise ValueError("Select one or more configured component IDs")
    for component in cfg["components"]:
        if component["id"] not in selected:
            continue
        name = component["archive_name"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*\.zip", name) or name in archive_names:
            raise ValueError("Archive names must be unique safe ZIP basenames")
        archive_names.add(name)
        if not component.get("readme") or not component.get("gates"):
            raise ValueError("Every archive requires its interpretation README and terminal gates")
        entries = {}
        gates = []
        protected = []

        def add(path, member, expected=None):
            p = canonical(path, allowed)
            if any(p.is_relative_to(x) for x in forbidden):
                raise ValueError(f"Forbidden environment/cache source: {p}")
            member = safe_member(member)
            if member in entries:
                # A gate and explicit inclusion may name the same canonical file.
                if entries[member]["source"] != str(p):
                    raise ValueError(f"Archive member collision: {member}")
                return
            rec = record(p) if expected is None else {"path": str(p), "bytes": expected.get("bytes", expected.get("size_bytes")), "sha256": expected["sha256"]}
            executable = p.suffix.lower() in (".sh", ".py", ".r") and bool(p.stat().st_mode & 0o111)
            entries[member] = {"member": member, "source": str(p), "bytes": rec["bytes"], "sha256": rec["sha256"], "portable_mode": 0o755 if executable else 0o644}

        for g in component["gates"]:
            native_root, seal, files = gate_files(g, allowed)
            protected.append(native_root)
            gates.append({"kind": g["kind"], "native_schema": seal["schema"], "native_status": seal["status"], "seal": g["seal"]})
            prefix = safe_member(g["archive_prefix"])
            if g.get("include_native_outputs", True):
                for p, relative, rec in files:
                    add(p, join_member(prefix, Path(relative).as_posix()), rec)
            else:
                seal_path = verify(g["seal"], allowed, native_root)
                add(seal_path, join_member(prefix, seal_path.name), g["seal"])
            reviews = [*g["review_receipts"], *g.get("display_review_receipts", [])]
            for i, rec in enumerate(reviews, 1):
                p = verify(rec, allowed)
                add(p, join_member(prefix, f"reviews/{i:02d}_{p.name}"), rec)
        for item in component.get("files", []):
            p = verify(item["receipt"], allowed)
            add(p, item["member"], item["receipt"])
        for tree in component.get("trees", []):
            root = canonical(tree["path"], allowed, directory=True)
            if any(root.is_relative_to(p) for p in forbidden):
                raise ValueError("Cannot archive a mutable software environment")
            protected.append(root)
            for p in sorted(root.rglob("*")):
                relative = p.relative_to(root)
                if any(part in SKIP for part in relative.parts) or p.suffix in (".pyc", ".pyo"):
                    continue
                if p.is_symlink():
                    raise ValueError(f"Tree has a symlink: {p}")
                if p.is_file():
                    add(p, join_member(tree["archive_prefix"], relative.as_posix()))
        if any(dest.is_relative_to(p) or p.is_relative_to(dest) for p in protected):
            raise ValueError("Archive destination overlaps a scientific/source tree")
        reserved = {"ARCHIVE_MANIFEST.json", "README.md"}
        if reserved.intersection(entries):
            raise ValueError("Native/source paths collide with transport metadata")
        rows = sorted(entries.values(), key=lambda r: r["member"])
        if not rows:
            raise ValueError("No files selected for archive")
        deduplicate = bool(component.get("deduplicate_identical_files", False))
        if deduplicate and component.get("role") != "shared_frozen_inputs_and_source":
            raise ValueError("Native result/model paths must remain intact; deduplication is shared-input-only")
        by_hash = {}
        for row in rows:
            identity = (row["sha256"], row["bytes"])
            row["stored_member"] = by_hash.setdefault(identity, row["member"]) if deduplicate else row["member"]
        stored = [r for r in rows if r["member"] == r["stored_member"]]
        raw = sum(r["bytes"] for r in stored)
        # Conservative raw-size preflight; all planned components are comfortably
        # below this ceiling, rather than relying on a compression-ratio guess.
        if raw + 1_000_000 > cap:
            raise ValueError(f"Uncompressed preflight exceeds component budget: {name}")
        component_sha = hashlib.sha256(json.dumps(component, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        components.append({"id": component["id"], "component_configuration_sha256": component_sha, "name": name, "role": component["role"], "readme": component["readme"], "gates": gates, "files": rows, "stored_files": stored, "stored_uncompressed_bytes": raw})
    if not components:
        raise ValueError("No archive components configured")
    if shutil.disk_usage(parent).free < sum(c["stored_uncompressed_bytes"] + 1_000_000 for c in components):
        raise ValueError("Insufficient destination free space for conservative archive staging")
    return {"configuration": record(configuration.resolve()), "destination": str(dest), "maximum_archive_bytes": cap, "components": components}


def zip_info(member: str, mode=0o644) -> zipfile.ZipInfo:
    z = zipfile.ZipInfo(safe_member(member), ZIP_DATE)
    z.compress_type = zipfile.ZIP_DEFLATED
    z.external_attr = (stat.S_IFREG | mode) << 16
    z.create_system = 3
    return z


def build(configuration: Path, component_ids=None) -> dict:
    plan = make_plan(configuration, component_ids)
    dest = Path(plan["destination"])
    # Nothing is written until ALL components, terminal seals, actual reviews,
    # input hashes, memberships and size/containment conditions have passed.
    dest.mkdir()
    archives = []
    for component in plan["components"]:
        target = dest / component["name"]
        manifest = {"schema": "satb2-followups-archive-transport/v1", "component_id": component["id"], "component_configuration_sha256": component["component_configuration_sha256"], "role": component["role"], "scientific_calculations_performed": False,
                    "native_gates": component["gates"], "files": component["files"], "deduplication": "Aliases may share stored_member only in the shared-input archive; extraction helper restores every requested relative member.",
                    "byte_semantics": "Native bytes are unchanged; SHA256 refers to each decompressed member."}
        with zipfile.ZipFile(target, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as archive:
            for row in component["stored_files"]:
                source = Path(row["source"])
                if source.stat().st_size != row["bytes"]:
                    raise ValueError("Source size changed after preflight")
                digest = hashlib.sha256()
                with source.open("rb") as inp, archive.open(zip_info(row["member"], row["portable_mode"]), "w", force_zip64=True) as out:
                    for block in iter(lambda: inp.read(CHUNK), b""):
                        digest.update(block)
                        out.write(block)
                if digest.hexdigest() != row["sha256"] or source.stat().st_size != row["bytes"]:
                    raise ValueError("Source changed during read-only archive transport")
                if target.stat().st_size > plan["maximum_archive_bytes"]:
                    raise ValueError("Archive budget exceeded; partial fresh output retained for review")
            archive.writestr(zip_info("README.md"), component["readme"].encode())
            archive.writestr(zip_info("ARCHIVE_MANIFEST.json"), (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode())
        if target.stat().st_size > plan["maximum_archive_bytes"]:
            raise ValueError("Final archive budget exceeded; no original files were changed")
        # Read back every stored member. ZIP CRC alone does not replace SHA256.
        with zipfile.ZipFile(target) as archive:
            if archive.testzip() is not None:
                raise ValueError("ZIP CRC verification failed")
            for row in component["stored_files"]:
                digest = hashlib.sha256()
                with archive.open(row["member"]) as stream:
                    for block in iter(lambda: stream.read(CHUNK), b""):
                        digest.update(block)
                if digest.hexdigest() != row["sha256"]:
                    raise ValueError("Archive decompressed SHA256 mismatch")
        archives.append({"component_id": component["id"], "component_configuration_sha256": component["component_configuration_sha256"], "name": target.name, "bytes": target.stat().st_size, "sha256": sha(target), "requested_members": len(component["files"]), "stored_members": len(component["stored_files"]), "stored_uncompressed_bytes": component["stored_uncompressed_bytes"], "role": component["role"]})
    receipt = {"schema": "satb2-followups-archive-delivery-receipt/v1", "status": "BUILT_AND_READBACK_VERIFIED_NOT_UPLOADED", "configuration": plan["configuration"], "maximum_archive_bytes": plan["maximum_archive_bytes"], "archives": archives, "scientific_source_files_modified": False, "scientific_calculations_performed": False}
    (dest / "ARCHIVE_DELIVERY_RECEIPT.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    shutil.copyfile(Path(__file__).with_name("extract_reviewed_archive.py"), dest / "extract_reviewed_archive.py")
    return receipt


def verify_deliveries(configuration: Path, receipts: list[Path]) -> dict:
    """Require every configured component across separately staged deliveries.

    This checks local delivery completeness/integrity, never upload/publication
    state or scientific inference. Per-component config receipts allow separate
    fresh invocations while an unrelated component is still in progress.
    """
    plan = make_plan(configuration, require_fresh_destination=False)
    cfg = json.loads(configuration.read_text())
    roots = [Path(p).resolve(strict=True) for p in cfg["allowed_output_roots"]]
    expected = {c["id"]: c for c in plan["components"]}
    found = {}
    for receipt_path in receipts:
        path = canonical(receipt_path, roots)
        delivery = json.loads(path.read_text())
        if delivery.get("schema") != "satb2-followups-archive-delivery-receipt/v1" or delivery.get("status") != "BUILT_AND_READBACK_VERIFIED_NOT_UPLOADED":
            raise ValueError("Invalid local archive delivery receipt")
        for rec in delivery["archives"]:
            ident = rec["component_id"]
            if ident not in expected or ident in found:
                raise ValueError("Missing, duplicate or unexpected release component")
            component = expected[ident]
            if rec["name"] != component["name"] or rec["component_configuration_sha256"] != component["component_configuration_sha256"]:
                raise ValueError("Archived component configuration changed")
            archive_path = verify({**rec, "path": str(path.parent / rec["name"])}, roots)
            if archive_path.stat().st_size > plan["maximum_archive_bytes"]:
                raise ValueError("Archived release asset exceeds budget")
            with zipfile.ZipFile(archive_path) as archive:
                manifest = json.loads(archive.read("ARCHIVE_MANIFEST.json"))
                if manifest.get("component_id") != ident or manifest.get("component_configuration_sha256") != component["component_configuration_sha256"] or manifest.get("files") != component["files"] or manifest.get("native_gates") != component["gates"]:
                    raise ValueError("Archive manifest no longer matches reviewed native inputs")
                wanted = {r["stored_member"] for r in component["files"]} | {"README.md", "ARCHIVE_MANIFEST.json"}
                if len(archive.namelist()) != len(wanted) or set(archive.namelist()) != wanted:
                    raise ValueError("Unexpected or duplicate archive members")
                for row in component["stored_files"]:
                    if archive.getinfo(row["member"]).file_size != row["bytes"]:
                        raise ValueError("Archived member size mismatch")
                    digest = hashlib.sha256()
                    with archive.open(row["member"]) as stream:
                        for block in iter(lambda: stream.read(CHUNK), b""):
                            digest.update(block)
                    if digest.hexdigest() != row["sha256"]:
                        raise ValueError("Archived member SHA256 mismatch")
            found[ident] = {"archive": str(archive_path), "bytes": rec["bytes"], "sha256": rec["sha256"], "native_gates": component["gates"]}
    if set(found) != set(expected):
        raise ValueError(f"Release components still missing: {sorted(set(expected) - set(found))}")
    return {"schema": "satb2-followups-complete-local-delivery-verification/v1", "status": "ALL_CONFIGURED_COMPONENTS_READBACK_VERIFIED_NOT_PUBLICATION_APPROVAL", "components": found}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configuration", required=True, type=Path)
    parser.add_argument("--component", action="append", help="Selected configured component ID; repeat to select several (default all)")
    parser.add_argument("--validate-only", action="store_true", help="Verify terminal reviews/receipts and budgets without writing archives")
    parser.add_argument("--verify-deliveries", nargs="+", type=Path, help="Require ALL configured components across these local delivery receipts")
    args = parser.parse_args()
    if args.verify_deliveries:
        if args.component or args.validate_only:
            parser.error("Full delivery verification requires all components and cannot combine with other modes")
        print(json.dumps(verify_deliveries(args.configuration, args.verify_deliveries), indent=2))
    elif args.validate_only:
        plan = make_plan(args.configuration, args.component)
        print(json.dumps({"status": "VERIFIED_NO_ARCHIVES_WRITTEN", "components": [{"name": c["name"], "stored_uncompressed_bytes": c["stored_uncompressed_bytes"], "requested_members": len(c["files"])} for c in plan["components"]]}, indent=2))
    else:
        print(json.dumps(build(args.configuration, args.component), indent=2))
