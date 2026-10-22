# Configuration

Use an absolute JSON `config_file` with `run_root` (new, absent),
`python_executable`, `inputs` mapping of path/bytes/sha256, and `endpoints`
(the 23 native DESeq2 PRIMARY endpoints and source hashes). `resource_file` is
an absolute pre-frozen ligand/receptor TSV with `resource_sha256`.
`liana_version` must equal 2.0.0. Set expected_nuclei=29367, expr_prop=0.10,
min_cells=5, n_perms=1000, seed=1337, threads=4, rank_cutoff=0.05 and
fdr_cutoff=0.05. The exact expected group counts are bound in the same JSON.

The prepared task environment and source inputs must already exist; this
workflow has no installation or download hooks. Inputs are checked by content
hash before and after execution. Endpoint statuses and mouse eligibility are
kept separate from communication ranks. Existing source result files remain
seven columns and are copied byte-for-byte into the new result tree.
