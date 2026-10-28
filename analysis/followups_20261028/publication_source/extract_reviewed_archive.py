#!/usr/bin/env python3
"""Safely restore exact archived member bytes and shared-input aliases; no code runs."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path, PurePosixPath


def safe(value):
    if not isinstance(value, str) or not value or value.startswith("/") or "\\" in value or ":" in value or any(ord(c) < 32 for c in value) or any(x in ("", ".", "..") or x != x.rstrip(" .") for x in value.split("/")) or PurePosixPath(value).is_absolute():
        raise ValueError("Unsafe archive member")
    return value


def extract(archive_path, destination, expected_sha256=None, maximum_restored_bytes=3_800_000_000):
    if expected_sha256 is not None:
        digest = hashlib.sha256()
        with Path(archive_path).open("rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != expected_sha256:
            raise ValueError("Archive SHA256 does not match the delivery receipt")
    destination = Path(destination)
    if destination.exists() or not destination.is_absolute():
        raise ValueError("Use a fresh absolute extraction directory")
    parent = destination.parent.resolve(strict=True)
    if any(p.is_symlink() for p in [destination, *destination.parents]):
        raise ValueError("Extraction directory must not contain symlinks")
    destination = parent / destination.name
    with zipfile.ZipFile(archive_path) as archive:
        if len(archive.namelist()) != len(set(archive.namelist())):
            raise ValueError("Duplicate ZIP members")
        for item in archive.infolist():
            safe(item.filename)
            if (item.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError("ZIP symlink rejected")
        manifest = json.loads(archive.read("ARCHIVE_MANIFEST.json"))
        if manifest.get("schema") != "satb2-followups-archive-transport/v1":
            raise ValueError("Unknown archive transport manifest")
        rows = manifest["files"]
        if len(rows) > 100_000 or sum(int(r["bytes"]) for r in rows) > maximum_restored_bytes:
            raise ValueError("Restored file budget exceeded")
        if len({r["member"] for r in rows}) != len(rows):
            raise ValueError("Duplicate requested archive paths")
        for row in rows:
            safe(row["member"]); safe(row["stored_member"])
            if int(row["bytes"]) < 0 or len(row["sha256"]) != 64:
                raise ValueError("Invalid member receipt")
            if row.get("portable_mode", 0o644) not in (0o644, 0o755):
                raise ValueError("Unexpected restored file permissions")
            if archive.getinfo(row["stored_member"]).file_size != row["bytes"]:
                raise ValueError("Archive/member size mismatch")
        if set(archive.namelist()) != {r["stored_member"] for r in rows} | {"README.md", "ARCHIVE_MANIFEST.json"}:
            raise ValueError("Unexpected archive files")
        destination.mkdir()
        for row in rows:
            output = destination / row["member"]
            output.parent.mkdir(parents=True, exist_ok=True)
            digest = hashlib.sha256()
            with archive.open(row["stored_member"]) as inp, output.open("xb") as out:
                for block in iter(lambda: inp.read(1024 * 1024), b""):
                    digest.update(block); out.write(block)
            if digest.hexdigest() != row["sha256"]:
                raise ValueError("Restored member SHA256 mismatch; fresh partial destination retained")
            output.chmod(row.get("portable_mode", 0o644))
        for member in ["README.md", "ARCHIVE_MANIFEST.json"]:
            (destination / member).write_bytes(archive.read(member))
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--expected-sha256", help="Expected full ZIP SHA256 from ARCHIVE_DELIVERY_RECEIPT.json")
    args = parser.parse_args()
    print(extract(args.archive, args.destination, args.expected_sha256))
