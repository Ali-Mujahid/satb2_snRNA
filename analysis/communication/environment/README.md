# New-task LIANA runtime — 2026-10-22

This environment is isolated under `SATB2_CCC_20261022/environment/.venv`; it does not inherit system or previous project site-packages. The existing Scanpy interpreter only created the new virtual environment. Use `run_python.sh` to keep Matplotlib, numba and temporary files inside this task.

## Version selection

PyPI currently advertises LIANA 1.10.0; its 2.0.0 endpoint returned 404. The official GitHub repository (now `scverse/liana`, redirected from `saezlab/liana-py`) exposes uppercase tag `V2.0.0` at immutable commit `c59472ccc9de8360dbbf5016db75f8abde08dd3e`. The task chose this fixed official release after source review identified a duplicate-score rank aggregation issue in 1.10.0. Version 1.10.0 was installed only during setup and replaced before any biological analysis. No local source patch was made.

The official source archive is preserved in `reference_cache/`; its SHA256 is `576d77d30eae572bf2f654d7d08c1e5fed2415272a166cda344d6abc14d91269`. `liana_source_manifest.json` records the tag API response, commit and archive identity. Since a GitHub archive omits Git metadata, the documented build-time `SETUPTOOLS_SCM_PRETEND_VERSION=2.0.0` supplied only the verified tag version. The first per-package override attempt was rejected by hatch-vcs and changed no installed distribution; its log remains for provenance. The successful second build log records the wheel hash.

## Offline readiness

Native imports and API signature inspection passed with Python socket connections blocked. LIANA bundles `omni_resource.csv`, and `select_resource('mouseconsensus')` required no network access. A copy of the 3,989-pair mouse ligand–receptor resource is preserved in `reference_cache/mouseconsensus.tsv`; downstream analysis should read this file and pass the DataFrame explicitly. No new ligand–receptor inference, biological fitting, normalization, UMAP or DESeq2 analysis was executed during runtime setup.

`runtime_api_audit.json` gives exact signatures, installed versions and resource checksums. `requirements.freeze.txt`, `requirements.lock.txt`, both pip installation reports, `package_artifacts.tsv` and `installed_liana_files.sha256` record the dependency closure and provenance. `pip check` passed. The historical standalone `docrep` package remains from the 1.10.0 setup dependency closure but is not required by the final LIANA 2.0.0 metadata.

## Previous environment integrity

The previous project's `pip freeze --all` output before and after setup is byte-identical (SHA256 `fe3de54b64846241ba250853c51fd70b9034bf2b38dad6ee41c5b6f719168332`). No install command targeted that environment or system Python. Preliminary metadata-only pip dry runs created disposable build metadata and the standard pip download cache before the task-local cache was configured; they did not install any package into the previous environment.

## Reproduction pattern

```bash
EXISTING_PYTHON=/Users/mujahid357/Documents/codex_new/SATB2_scanpy_seurat_mapmycells_20261013/.venv_scanpy/bin/python3.12
CCC_ENV=/absolute/path/to/new/environment
"$EXISTING_PYTHON" -m venv "$CCC_ENV/.venv"
# Use the preserved archive and dependency lock; keep build/cache directories task-local.
SETUPTOOLS_SCM_PRETEND_VERSION=2.0.0 TMPDIR="$CCC_ENV/tmp" PIP_CACHE_DIR="$CCC_ENV/wheel_cache" \
  "$CCC_ENV/.venv/bin/python" -m pip install --report "$CCC_ENV/pip_install_report.json" \
  -r /path/to/requirements.lock.txt
```

The lock is the macOS arm64 / Python 3.12 artifact closure. Its wheel hashes are specific to these artifacts; another operating system needs its own resolver/report rather than silently reusing incompatible wheels. Runtime package artifacts were downloaded only during setup. Actual analysis must use this prepared environment and cached resource without network.
