"""Annotate directed candidates using the four previously frozen reference lookups."""
import csv
from itertools import combinations
from pathlib import Path

RESOURCES = ("SynDive", "Synaptome.db", "G2Cdb", "Connectome-seq")


def write_table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def compare_pairs(data, output: Path):
    roles = []
    for pair in data["pairs"]:
        ligand, receptor = pair["ligand"], pair["receptor"]
        ca_sends = pair["source"].startswith("CA1_")
        assert ca_sends != pair["target"].startswith("CA1_")
        ca_role, in_role = ("ligand", "receptor") if ca_sends else ("receptor", "ligand")
        criteria = {
            "SynDive": bool(pair["FASS_common_strata"]),
            "Synaptome.db": pair["mouse_hippocampus_both"],
            "G2Cdb": bool(pair["G2C_shared_lists"]),
            "Connectome-seq": all(g in data["markers"] for g in (ligand, receptor)),
        }
        assert all(isinstance(x, bool) for x in criteria.values())
        hits = [name for name in RESOURCES if criteria[name]]
        row = {
            "CA1_sublayer": pair["source"] if ca_sends else pair["target"],
            "CA1_gene": pair[ca_role],
            "CA1_role": ca_role,
            "CA1_change": "up" if pair[ca_role + "_log2FC"] > 0 else "down",
            "CA1_log2FC": pair[ca_role + "_log2FC"],
            "CA1_FDR": pair[ca_role + "_FDR"],
            "IN_cluster": pair["target"] if ca_sends else pair["source"],
            "IN_gene": pair[in_role],
            "IN_role": in_role,
            "IN_change": "up" if pair[in_role + "_log2FC"] > 0 else "down",
            "IN_log2FC": pair[in_role + "_log2FC"],
            "IN_FDR": pair[in_role + "_FDR"],
            "reference_count": len(hits),
            "exact_reference_set": " + ".join(hits) or "None",
            **criteria,
            **pair,
        }
        roles.append(row)
    def sortkey(r):
        return (r["CA1_role"] != "ligand", r["ligand"], r["receptor"], r["CA1_sublayer"], r["IN_cluster"])
    roles.sort(key=sortkey)
    signatures = []
    for n in range(len(RESOURCES) + 1):
        for subset in combinations(RESOURCES, n):
            signature = " + ".join(subset) or "None"
            selected = [r for r in roles if r["exact_reference_set"] == signature]
            identities = sorted({(r["ligand"], r["receptor"]) for r in selected})
            signatures.append({
                "exact_reference_set": signature, "resource_count": n,
                "directed_population_rows": len(selected), "ordered_gene_pairs": len(identities),
                "unordered_gene_pairs": len({tuple(sorted(p)) for p in identities}),
                "ligand_to_receptor": " | ".join(f"{a} → {b}" for a,b in identities),
            })
    write_table(output / "CA1_IN_directed_pair_roles.tsv", roles)
    write_table(output / "exact_reference_overlap.tsv", signatures)
    return roles
