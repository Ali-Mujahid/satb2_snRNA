#!/usr/bin/env python3
"""Create a new path-remapped JSON config; never execute scientific software.

Original configurations and scientific payloads are never modified. Paths in
referenced TSVs, receipts and AnnData/R objects are deliberately outside scope.
"""
import argparse
import hashlib
import json
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--path-map", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--drop-v1-reuse", action="store_true",
                        help="Remove v1 recovery-cache keys from the new copy.")
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    receipt = output.with_suffix(output.suffix + ".path_remap.json")
    if output.exists() or receipt.exists():
        parser.error("Refusing an existing output or path-remap receipt.")
    original = json.loads(args.config.read_text())
    spec = json.loads(args.path_map.read_text())
    variables = spec.get("variables", {})
    if any("<" in str(v) or ">" in str(v) for v in variables.values()):
        parser.error("Fill all placeholder variables in a separate path-map copy first.")
    prefix_map = []
    for item in spec["prefix_map"]:
        source = item["from"].rstrip("/")
        target = item["to"].format_map(variables).rstrip("/")
        if not target.startswith("/") or "<" in target or ">" in target:
            parser.error("Every replacement must be a filled absolute path.")
        prefix_map.append((source, target))
    prefix_map.sort(key=lambda pair: -len(pair[0]))
    changed, unmapped = [], []

    def remap(value, key="$"):
        if isinstance(value, dict):
            return {k: remap(v, key + "." + k) for k, v in value.items()}
        if isinstance(value, list):
            return [remap(v, key + "[" + str(i) + "]") for i, v in enumerate(value)]
        if not isinstance(value, str) or not value.startswith("/"):
            return value
        for source, target in prefix_map:
            if value == source or value.startswith(source + "/"):
                replacement = target + value[len(source):]
                changed.append({"key": key, "from": value, "to": replacement})
                return replacement
        unmapped.append({"key": key, "path": value})
        return value

    result = remap(original)
    dropped = []
    if args.drop_v1_reuse:
        for key in ("reuse_prepared_shared", "reuse_completed_scanpy"):
            if key in result:
                del result[key]
                dropped.append(key)
        unmapped = [x for x in unmapped if not any(
            x["key"].startswith("$." + key) for key in dropped)]
    if unmapped:
        parser.error("Add explicit mappings for these absolute paths: " +
                     json.dumps(unmapped, indent=2))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    receipt.write_text(json.dumps({
        "schema": "satb2-config-path-remap-1",
        "original_config_sha256": sha256(args.config),
        "path_map_sha256": sha256(args.path_map),
        "new_config_sha256": sha256(output),
        "changed_paths": changed,
        "dropped_recovery_keys": dropped,
        "scientific_execution_performed": False,
        "referenced_input_files_validated": False,
        "caveat": "JSON paths only. Referenced manifests still need their own "
                  "documented path remapping and all source-bound data/receipts."
    }, indent=2) + "\n")
    print("Created new configuration and path-remap receipt; no analysis executed.")


if __name__ == "__main__":
    main()
