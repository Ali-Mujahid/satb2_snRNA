#!/usr/bin/env python3
"""Package a frozen offline marker-atlas derivative; never run scientific analysis.

The two release ZIPs share offline_dashboard_v3/. Existing source archives are
read-only. Unchanged members retain their exact original compressed payloads.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import struct
import subprocess
import sys
import tempfile
import uuid
import zipfile
import zlib

CHUNK = 8 * 1024 * 1024
PREFIX = "offline_dashboard_v3"
MANIFEST = "BUNDLE_MANIFEST.json"
ASSET_LIMIT = 2 * 1024**3
SOURCE_PREFIX = "offline_dashboard_v2"


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def check(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def file_record(path, relative):
    h = hashlib.sha256()
    crc = 0
    size = 0
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
            crc = zlib.crc32(block, crc)
            size += len(block)
    return {"path": relative, "bytes": size, "sha256": h.hexdigest(), "crc32": crc & 0xffffffff}


def safe_relative(name):
    p = PurePosixPath(name)
    check(name and not name.startswith("/") and ".." not in p.parts,
          f"Unsafe relative member: {name!r}")
    check(not any(ord(c) < 32 or c == "\\" for c in name),
          f"Unsupported control/backslash in member: {name!r}")
    return name


def folder_records(root, include_manifest=False):
    check(root.is_dir() and not root.is_symlink(), "Bundle must be a real directory")
    records = {}
    for path in sorted(root.rglob("*")):
        check(not path.is_symlink(), f"Symlinks are not portable bundle members: {path.name}")
        if path.is_dir():
            continue
        check(path.is_file(), f"Non-file bundle member: {path.name}")
        rel = safe_relative(path.relative_to(root).as_posix())
        if include_manifest or rel != MANIFEST:
            records[rel] = file_record(path, rel)
    return records


def extension_record(path):
    raw = path.read_bytes()
    r = json.loads(raw)
    check(r.get("status") == "COMPLETE", "Marker-extension receipt is not COMPLETE")
    check(r.get("new_descriptive_marker_count_aggregation_performed") is True,
          "Marker-extension expression-aggregation scope is missing")
    check(r.get("scientific_model_fitting_performed") is False,
          "Marker extension must not fit models")
    check(r.get("QC_clustering_embeddings_annotation_and_DESeq2_unchanged") is True,
          "Marker extension must preserve QC/clustering/embedding/annotation/DESeq2")
    check(r.get("source_run_modified") is False, "Marker extension modified the source run")
    return {"receipt_name": path.name, "receipt_sha256": hashlib.sha256(raw).hexdigest(),
            "schema": r.get("schema", r.get("schema_version")),
            "n_datasets": r.get("n_datasets"), "n_original_genes": r.get("n_original_genes"),
            "n_extended_genes": r.get("n_extended_genes"),
            "n_screenshot_genes": r.get("n_screenshot_genes"),
            "new_descriptive_marker_count_aggregation_performed": True,
            "scientific_model_fitting_performed": False,
            "QC_clustering_embeddings_annotation_and_DESeq2_unchanged": True,
            "source_run_modified": False}


def prepare_manifest(root, extension_path):
    check(root.name == PREFIX, f"Bundle directory must be named {PREFIX}")
    extension = extension_record(extension_path)
    rows = folder_records(root)
    check("index.html" in rows, "Offline entrypoint index.html is absent")
    previous = digest(root / MANIFEST) if (root / MANIFEST).exists() else None
    manifest = {"schema": "satb2-offline-bundle-3", "generated_at": now(),
                "entrypoint": "index.html", "file_count": len(rows),
                "bytes": sum(r["bytes"] for r in rows.values()),
                "scientific_fitting_performed": False,
                "DE_model_fitting_performed": False,
                "embedding_recomputation_performed": False,
                "raw_nucleus_count_aggregation_performed": True,
                "count_aggregation_scope": "Descriptive added-marker summaries from existing selected-source count matrices only; no DE or embedding recomputation",
                "new_descriptive_marker_count_aggregation_performed": True,
                "previous_analysis_files_modified": False,
                "marker_extension": extension,
                "previous_folder_manifest_sha256": previous,
                "files": [{k: r[k] for k in ("path", "bytes", "sha256")} for r in rows.values()]}
    target = root / MANIFEST
    check(not target.is_symlink(), "Refusing a symlinked manifest")
    temp = root / (MANIFEST + ".tmp-" + uuid.uuid4().hex)
    temp.write_text(json.dumps(manifest, indent=2) + "\n")
    os.replace(temp, target)
    check(extension_record(extension_path) == extension, "Extension receipt changed during snapshot")
    return {"status": "MANIFEST_PREPARED", "bundle": root.name,
            "manifest_sha256": digest(target), "files_excluding_manifest": len(rows),
            "bytes_excluding_manifest": manifest["bytes"], "extension_receipt_sha256": extension["receipt_sha256"]}


def load_frozen_manifest(root, extension_path):
    raw = (root / MANIFEST).read_bytes()
    m = json.loads(raw)
    check(m.get("schema") == "satb2-offline-bundle-3", "Prepare the v3 manifest before packaging")
    check(m.get("marker_extension") == extension_record(extension_path), "Extension receipt binding changed")
    rows = folder_records(root)
    expected = {x["path"]: x for x in m["files"]}
    check(len(expected) == len(m["files"]) == m["file_count"], "Duplicate or invalid manifest rows")
    check(set(rows) == set(expected), "Frozen folder file set changed after manifest snapshot")
    for rel, record in rows.items():
        check(all(record[k] == expected[rel][k] for k in ("bytes", "sha256")),
              f"Frozen folder bytes changed after manifest snapshot: {rel}")
    check(sum(r["bytes"] for r in rows.values()) == m["bytes"], "Manifest byte total mismatch")
    rows[MANIFEST] = file_record(root / MANIFEST, MANIFEST)
    return m, rows, hashlib.sha256(raw).hexdigest()


def archive_map(zf, prefix):
    result = {}
    for info in zf.infolist():
        name = safe_relative(info.filename.rstrip("/"))
        check(name == prefix or name.startswith(prefix + "/"), "Source ZIP contains a foreign root")
        if info.is_dir():
            continue
        rel = name[len(prefix) + 1:]
        safe_relative(rel)
        check(rel not in result, f"Duplicate ZIP member: {rel}")
        check(not (info.flag_bits & 1), "Encrypted ZIP members are unsupported")
        result[rel] = info
    return result


def clone_copy(source, target):
    if sys.platform == "darwin":
        p = subprocess.run(["/bin/cp", "-c", str(source), str(target)], capture_output=True)
        if p.returncode == 0:
            return "APFS_clone_cp_c"
        if target.exists():
            target.unlink()
    shutil.copy2(source, target)
    return "byte_copy"


def extra_has_unicode_path(extra):
    i = 0
    while i + 4 <= len(extra):
        tag, n = struct.unpack_from("<HH", extra, i)
        if tag == 0x7075:
            return True
        i += 4 + n
    return False


def rename_headers_same_length(path, old, new):
    """Fallback: change filename bytes only, preserving offsets and all payloads.

    Python's ZIP parser supplies ZIP64-aware local/central offsets. We refuse
    Unicode-path override extra fields instead of leaving stale filename data.
    """
    check(len(old.encode()) == len(new.encode()), "Header fallback requires equal-length prefixes")
    with zipfile.ZipFile(path) as zf:
        infos = zf.infolist()
        central = zf.start_dir
        archive_map(zf, old)
    ob, nb = old.encode("ascii"), new.encode("ascii")
    with path.open("r+b") as fh:
        for info in infos:
            fh.seek(info.header_offset)
            head = fh.read(30)
            fields = struct.unpack("<IHHHHHIIIHH", head)
            check(fields[0] == 0x04034b50, "Invalid local ZIP header")
            name, extra = fh.read(fields[9]), fh.read(fields[10])
            check(not extra_has_unicode_path(extra), "Unicode-path extra field requires zipnote")
            check(name == ob or name.startswith(ob + b"/"), "Unexpected local filename root")
            fh.seek(info.header_offset + 30)
            fh.write(nb + name[len(ob):])
        fh.seek(central)
        for _ in infos:
            start = fh.tell()
            fields = struct.unpack("<IHHHHHHIIIHHHHHII", fh.read(46))
            check(fields[0] == 0x02014b50, "Invalid central ZIP header")
            name, extra = fh.read(fields[10]), fh.read(fields[11])
            check(not extra_has_unicode_path(extra), "Unicode-path extra field requires zipnote")
            check(name == ob or name.startswith(ob + b"/"), "Unexpected central filename root")
            end = start + 46 + fields[10] + fields[11] + fields[12]
            fh.seek(start + 46)
            fh.write(nb + name[len(ob):])
            fh.seek(end)
    return "same_length_header_fallback"


def rename_archive_prefix(path, old, new, method="auto"):
    note = shutil.which("zipnote")
    if method == "header" or not note:
        return rename_headers_same_length(path, old, new)
    text = subprocess.run([note, str(path)], check=True, capture_output=True).stdout
    lines = text.splitlines(keepends=True)
    patched = []
    count = 0
    for line in lines:
        patched.append(line)
        if line.startswith(b"@ ") and not line.startswith(b"@ ("):
            name = line[2:].rstrip(b"\r\n")
            check(name == old.encode() or name.startswith(old.encode() + b"/"), "Foreign zipnote entry")
            patched.append(b"@=" + new.encode() + name[len(old):] + b"\n")
            count += 1
    check(count > 0, "No ZIP members to rename")
    before_sha = digest(path)
    result = subprocess.run([note, "-w", str(path)], input=b"".join(patched), capture_output=True)
    if result.returncode != 0:
        check(digest(path) == before_sha, "Failed zipnote modified its staged archive; refusing fallback")
        rename_headers_same_length(path, old, new)
        return "same_length_header_fallback_after_zipnote_failure"
    return "zipnote_filename_rename"


def payload_digest(fh, info):
    fh.seek(info.header_offset)
    fields = struct.unpack("<IHHHHHIIIHH", fh.read(30))
    check(fields[0] == 0x04034b50, "Invalid compressed-member local header")
    fh.seek(fields[9] + fields[10], 1)
    h = hashlib.sha256()
    remaining = info.compress_size
    while remaining:
        block = fh.read(min(CHUNK, remaining))
        check(bool(block), "Truncated compressed payload")
        h.update(block)
        remaining -= len(block)
    return h.hexdigest()


def zipped_sha(zf, info):
    h = hashlib.sha256()
    with zf.open(info) as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)  # Reading through EOF verifies the member CRC.
    return h.hexdigest()


def update_zip(path, root, changed, removed):
    exe = shutil.which("zip")
    check(exe is not None, "Info-ZIP 'zip' is required for changed runtime members")
    if removed:
        listing = "".join(PREFIX + "/" + r + "\n" for r in removed).encode()
        subprocess.run([exe, "-q", "-d", str(path), "-@"], input=listing, check=True, capture_output=True)
    if changed:
        listing = "".join(PREFIX + "/" + r + "\n" for r in changed).encode()
        # No -u: a changed member must replace the old bytes even with an older mtime.
        subprocess.run([exe, "-q", "-D", str(path), "-@"], cwd=root.parent,
                       input=listing, check=True, capture_output=True)


def build_one(source, target, root, expected, baseline_hashes, role, rename_method):
    with zipfile.ZipFile(source) as before:
        old = archive_map(before, SOURCE_PREFIX)
        desired = {p: r for p, r in expected.items() if (p.startswith("saved_results/") == (role == "saved_results"))}
        changed, unchanged = [], []
        for rel, r in desired.items():
            info = old.get(rel)
            baseline = baseline_hashes.get(rel)
            same = info is not None and info.file_size == r["bytes"] and info.CRC == r["crc32"]
            if same:
                same = (baseline is not None and baseline.get("sha256") == r["sha256"]
                        and baseline.get("bytes") == r["bytes"])
                if baseline is None:
                    same = zipped_sha(before, info) == r["sha256"]
            (unchanged if same else changed).append(rel)
        removed = sorted(set(old) - set(desired))
        added = sorted(set(desired) - set(old))
    if role == "saved_results":
        check(not changed and not removed and not added,
              "saved_results differs from its preserved source archive; refusing a rename-only companion")
    copied = clone_copy(source, target)
    renamed = rename_archive_prefix(target, SOURCE_PREFIX, PREFIX, rename_method)
    update_zip(target, root, sorted(changed), removed)
    check(target.stat().st_size < ASSET_LIMIT, f"Release asset must be strictly below 2 GiB: {target.name}")
    proofs = []
    unchanged_set = set(unchanged)
    with zipfile.ZipFile(source) as before, zipfile.ZipFile(target) as after, source.open("rb") as bf, target.open("rb") as af:
        old = archive_map(before, SOURCE_PREFIX)
        actual = archive_map(after, PREFIX)
        check(set(actual) == set(desired), f"Final {role} member set mismatch")
        for rel, row in desired.items():
            info = actual[rel]
            check(info.file_size == row["bytes"] and info.CRC == row["crc32"], f"Final size/CRC mismatch: {rel}")
            check(zipped_sha(after, info) == row["sha256"], f"Final decompressed SHA mismatch: {rel}")
            if rel in unchanged_set:
                prev = old[rel]
                check((prev.CRC, prev.file_size, prev.compress_size, prev.compress_type) ==
                      (info.CRC, info.file_size, info.compress_size, info.compress_type),
                      f"Unchanged member compression metadata changed: {rel}")
                psha, nsha = payload_digest(bf, prev), payload_digest(af, info)
                check(psha == nsha, f"Unchanged compressed payload changed: {rel}")
                proofs.append({"path": rel, "compressed_sha256": psha,
                               "compressed_bytes": info.compress_size, "compressed_identity_PASS": True})
    return {"asset_name": target.name, "role": role, "bytes": target.stat().st_size,
            "sha256": digest(target), "copy_method": copied, "rename_method": renamed,
            "member_count": len(desired), "all_member_CRC_and_SHA_PASS": True,
            "under_2GiB": True, "unchanged_compressed_members": len(unchanged),
            "unchanged_compressed_identity_PASS": True, "added_members": added,
            "removed_owned_baseline_members": removed,
            "changed_or_added_members": [{"path": p, "bytes": desired[p]["bytes"],
                                          "sha256": desired[p]["sha256"]} for p in sorted(changed)],
            "unchanged_compressed_payload_proofs": proofs}


def package(root, extension_path, runtime_source, saved_source, outdir, rename_method="auto"):
    check(root.name == PREFIX, f"Bundle directory must be named {PREFIX}")
    check(not outdir.is_relative_to(root), "Package output directory must be outside the bundle")
    outdir.mkdir(parents=True, exist_ok=True)
    targets = [outdir / "SATB2_atlas_v3_runtime.zip", outdir / "SATB2_atlas_v3_saved_results.zip"]
    check(all(not p.exists() for p in targets), "Refusing an existing v3 output asset")
    receipt_path = outdir / "OFFLINE_V3_PACKAGE_RECEIPT.json"
    check(not receipt_path.exists(), "Refusing an existing package receipt")
    check(not (outdir / "OFFLINE_V3_PUBLIC_PROVENANCE.json").exists(), "Refusing existing public package provenance")
    check(runtime_source.resolve() != saved_source.resolve(), "Source archives must be distinct")
    for src in [runtime_source, saved_source]:
        check(src.is_file() and src.resolve() not in [p.resolve() for p in targets], "Invalid source/output archive path")
    manifest, expected, manifest_sha = load_frozen_manifest(root, extension_path)
    sources = {p.name: {"path": str(p), "sha256": digest(p), "bytes": p.stat().st_size}
               for p in [runtime_source, saved_source]}
    with zipfile.ZipFile(runtime_source) as zf:
        previous = json.loads(zf.read(SOURCE_PREFIX + "/" + MANIFEST))
    baseline = {r["path"]: r for r in previous["files"]}
    work = outdir / (".offline_v3_build_" + uuid.uuid4().hex)
    work.mkdir()
    assets = []
    try:
        for source, target, role in [(runtime_source, targets[0], "runtime"), (saved_source, targets[1], "saved_results")]:
            staged = work / target.name
            print("Verifying and packaging " + target.name, flush=True)
            assets.append(build_one(source, staged, root, expected, baseline, role, rename_method))
        # Source folder and source archives must remain byte-identical during packaging.
        _, again, again_sha = load_frozen_manifest(root, extension_path)
        check(again == expected and again_sha == manifest_sha, "Bundle changed during packaging")
        for src in [runtime_source, saved_source]:
            check(digest(src) == sources[src.name]["sha256"], "Preserved source archive changed")
        for staged, target in zip([work / p.name for p in targets], targets):
            check(not target.exists(), "Output asset appeared concurrently")
            os.rename(staged, target)
        receipt = {"schema": "satb2-offline-v3-package-1", "status": "PACKAGED_AND_VERIFIED",
                   "generated_at": now(), "bundle": PREFIX, "bundle_manifest_sha256": manifest_sha,
                   "folder_files_including_manifest": len(expected), "source_archives": sources,
                   "marker_extension": manifest["marker_extension"], "assets": assets,
                   "source_archives_unchanged_PASS": True, "source_folder_snapshot_unchanged_PASS": True,
                   "scientific_fitting_performed_by_packager": False,
                   "RNA_count_aggregation_performed_by_packager": False,
                   "existing_saved_results_file_bytes_unchanged_PASS": True,
                   "extract_together": "Extract both ZIPs into one parent directory; both use offline_dashboard_v3/",
                   "private_local_receipt": True}
        receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
        work.rmdir()
        public = {k: receipt[k] for k in ["schema", "status", "generated_at", "bundle", "bundle_manifest_sha256",
                  "folder_files_including_manifest", "marker_extension", "source_archives_unchanged_PASS",
                  "source_folder_snapshot_unchanged_PASS", "scientific_fitting_performed_by_packager",
                  "RNA_count_aggregation_performed_by_packager", "existing_saved_results_file_bytes_unchanged_PASS", "extract_together"]}
        public["source_archives"] = [{"asset_name": p.name, "sha256": sources[p.name]["sha256"],
                                      "bytes": sources[p.name]["bytes"]} for p in [runtime_source, saved_source]]
        public["assets"] = [{k: v for k, v in x.items() if k != "unchanged_compressed_payload_proofs"} for x in assets]
        (outdir / "OFFLINE_V3_PUBLIC_PROVENANCE.json").write_text(json.dumps(public, indent=2) + "\n")
        print(json.dumps({"status": receipt["status"], "receipt": str(receipt_path),
                          "assets": [{k: x[k] for k in ["asset_name", "bytes", "sha256", "member_count",
                                                               "unchanged_compressed_members"]} for x in assets]}, indent=2))
        return receipt
    except Exception:
        (work / "FAILED_BUILD.json").write_text(json.dumps({"status": "FAILED_PRESERVED_FOR_REVIEW", "generated_at": now()}) + "\n")
        raise


def self_test():
    """Fabricated ZIPs only; never opens the study folder or study archives."""
    with tempfile.TemporaryDirectory(prefix="satb2_v3_packaging_fixture_") as tmp:
        base = Path(tmp)
        old = base / SOURCE_PREFIX
        (old / "saved_results").mkdir(parents=True)
        fixture = {"index.html": b"old app", "unchanged.bin": bytes(range(256)) * 400,
                   "deleted-owned.txt": b"orphan", "saved_results/genes.tsv": b"gene\tvalue\nG\t1\n"}
        for name, data in fixture.items():
            (old / name).write_bytes(data)
        oldmanifest = {"files": [{"path": p, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()} for p, b in fixture.items()]}
        (old / MANIFEST).write_text(json.dumps(oldmanifest))
        sources = [base / "old_runtime.zip", base / "old_saved.zip"]
        for source, saved in zip(sources, [False, True]):
            with zipfile.ZipFile(source, "w", zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as z:
                for p in sorted(old.rglob("*")):
                    if p.is_file() and (p.relative_to(old).as_posix().startswith("saved_results/") == saved):
                        z.write(p, SOURCE_PREFIX + "/" + p.relative_to(old).as_posix())
        original_shas = [digest(p) for p in sources]
        new = base / PREFIX
        shutil.copytree(old, new)
        (new / "index.html").write_text("new UI older timestamp")
        os.utime(new / "index.html", (946684800, 946684800))
        (new / "deleted-owned.txt").unlink()
        (new / "new_marker.js").write_text("fixture only: no RNA data")
        ext = base / "fixture_extension.json"
        ext.write_text(json.dumps({"status": "COMPLETE", "schema": "synthetic-only",
                                  "new_descriptive_marker_count_aggregation_performed": True,
                                  "scientific_model_fitting_performed": False,
                                  "QC_clustering_embeddings_annotation_and_DESeq2_unchanged": True,
                                  "source_run_modified": False, "n_original_genes": 2,
                                  "n_extended_genes": 3, "n_screenshot_genes": 3, "n_datasets": 1}))
        prepare_manifest(new, ext)
        a = package(new, ext, *sources, base / "zipnote_out")
        b = package(new, ext, *sources, base / "header_fallback_out", rename_method="header")
        check([digest(p) for p in sources] == original_shas, "Fixture source archives changed")
        check(all(x["unchanged_compressed_identity_PASS"] for r in [a, b] for x in r["assets"]), "Fixture compressed preservation failed")
        # A stale snapshot must be rejected before producing final archives.
        (new / "new_marker.js").write_text("changed after snapshot")
        try:
            load_frozen_manifest(new, ext)
        except ValueError:
            pass
        else:
            raise AssertionError("Stale fixture manifest was accepted")
        (new / "new_marker.js").write_text("fixture only: no RNA data")
        (new / "saved_results/genes.tsv").write_text("changed saved results")
        prepare_manifest(new, ext)
        _, e, _ = load_frozen_manifest(new, ext)
        try:
            build_one(sources[1], base / "must_not_exist.zip", new, e,
                      {r["path"]: r for r in oldmanifest["files"]}, "saved_results", "auto")
        except ValueError:
            check(not (base / "must_not_exist.zip").exists(), "Saved-change refusal wrote an asset")
        else:
            raise AssertionError("Changed saved-results fixture was accepted")
    print("PASS synthetic ZIP prefix rename, compressed-byte preservation (automatic rename and forced header fallback), changed-old-mtime replacement, deletion/addition, folder SHA/CRC, stale-snapshot rejection and saved-results immutability; no study input read.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare-manifest", action="store_true")
    mode.add_argument("--package", action="store_true")
    mode.add_argument("--self-test", action="store_true")
    ap.add_argument("--bundle-root", type=Path)
    ap.add_argument("--extension-receipt", type=Path)
    ap.add_argument("--runtime-source-zip", type=Path)
    ap.add_argument("--saved-source-zip", type=Path)
    ap.add_argument("--outdir", type=Path)
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return
    check(args.bundle_root is not None and args.extension_receipt is not None,
          "--bundle-root and --extension-receipt are required")
    root, ext = args.bundle_root.resolve(), args.extension_receipt.resolve()
    if args.prepare_manifest:
        print(json.dumps(prepare_manifest(root, ext), indent=2))
    else:
        check(all([args.runtime_source_zip, args.saved_source_zip, args.outdir]),
              "--package requires both source archives and --outdir")
        package(root, ext, args.runtime_source_zip.resolve(), args.saved_source_zip.resolve(), args.outdir.resolve())


if __name__ == "__main__":
    main()
