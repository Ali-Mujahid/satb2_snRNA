"""Write the display overlay after an exact join to the deployed Seurat inhibitory shard."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

root = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("assignments", type=Path, help="Reviewed C6-only nucleus_assignments.tsv")
args = parser.parse_args()
shard = root / "site/data/shards/coordinates__seurat__original_full__Inhibitory.json.gz"
with gzip.open(shard, "rt") as handle:
    nuclei = json.load(handle)["nuclei"]
with args.assignments.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
parent = {row["nucleus_id"]: row for row in nuclei if str(row["cluster_id"]) == "6"}
groups = {r["nucleus_id"]: r["subcluster"] for r in rows}
assert len(nuclei) == len({r["nucleus_id"] for r in nuclei}) == 8911
assert len(rows) == len(groups) == 501 and parent.keys() == groups.keys()
assert Counter(groups.values()) == {"0": 291, "1": 110, "2": 100}
for row in rows:
    source = parent[row["nucleus_id"]]
    assert all(source[key] == row[key] for key in ("sample", "condition"))
    for dimension, axis in enumerate(("x", "y"), 1):
        assert source[f"{axis}_harmony"] == float(row[f"UMAP_harmony_{dimension}"])

payload = {
    "source_sha256": hashlib.sha256(args.assignments.read_bytes()).hexdigest(),
    "coordinate_sha256": hashlib.sha256(shard.read_bytes()).hexdigest(),
    "groups": {f"6.{i}": sorted(n for n, g in groups.items() if g == str(i)) for i in range(3)},
}
out = root / "site/data/c6_subclusters.js"
out.write_text("window.SATB2_C6_OVERLAY = " + json.dumps(payload, separators=(",", ":")) + ";\n")
