# Isolated stable pySCENIC environment

Environment: `/Users/mujahid357/Documents/codex_new/SATB2_REGULATORY_20261026/environment/pyscenic39`

Base interpreter: `/usr/local/bin/python3.9`, Python3.9.1, macOS x86_64. This new virtual environment does not modify the R, existing CCC, controller or system environments. The direct pins are in `requirements.direct.txt`; the complete68-package install resolution is in `requirements.resolved.txt` and `requirements.locked-hashes.txt`. `requirements.freeze.txt` also includes bootstrap pip/setuptools versions.

## Validation scope

- Dependency resolution succeeded before installation, against official PyPI only. Every installed distribution's recorded archive URL is on `files.pythonhosted.org` or `pypi.org`.
- `pip check`: no broken requirements.
- All18 import probes passed, including numpy, pandas, scipy, sklearn, pyarrow, numba, dask/dataframe, distributed, arboreto, ctxcore, loompy, umap and the pySCENIC CLI/pruning/AUCell modules.
- Function signatures were inspected for Arboreto GRNBoost2/create_graph, Dask from_delayed, FeatherRankingDatabase, prune2df and aucell. No Dask cluster, graph computation, GRN, motif pruning, AUCell scoring, clustering or synthetic scientific fit was started.
- The stable0.12.1 executable **does not support `pyscenic --version`**: its tested exit status is2 and its complete output is retained in `pyscenic_version.txt`. `pyscenic --help` succeeds with exit0 and its description reports0.12.1. The independent package metadata probe also confirms0.12.1.
- The pinned NumPy1.23.5 retains the aliases used by this stable release; its upstream deprecation notices and loompy/numba deprecation notices are retained. No third-party package source was patched.
- Static source inspection confirms Dask2023.3.2's `from_delayed` supports an empty list when `meta` is supplied, matching Arboreto0.1.6's `create_graph` metadata-table construction. This is an interface compatibility check, not an executed GRN test.
- No anndata package was added. The planned workflow can read verified MTX+genes/cells and full-library denominator metadata directly, using DataFrames. Loom/HDF5 support is present through loompy/h5py.

Runtime import readiness is established. Scientific workflow readiness, mouse cisTarget resource validation, immutable planning and approval remain separate tasks handled by the parent workflow.

## Provenance

`resolution.report.json` records the successful pre-install pip dry run. `install.report.json` records actual package archives, metadata, URLs and SHA256. `resolved_downloads.tsv` is a compact inventory. `bootstrap_manifest.json` records base interpreter/platform and official pip24.3.1 wheel URL/hash. `pip_upgrade.log`, `resolution.log`, `install.log` and all validation outputs preserve the setup evidence. `import_probe.json` records exact versions and inspected signatures. `SHA256SUMS` hashes these setup artifacts.

## Reproduction on the same platform

Choose a new empty destination; do not install over another analysis environment. The installed archive hashes target this macOS x86_64/Python3.9 platform. For a different platform, resolve the direct pins separately and save that platform's new report rather than substituting unmatched wheels.

```sh
/usr/local/bin/python3.9 -m venv /absolute/new/pyscenic39
/absolute/new/pyscenic39/bin/python -m pip install --index-url https://pypi.org/simple --no-cache-dir pip==24.3.1
/absolute/new/pyscenic39/bin/python -m pip install --index-url https://pypi.org/simple --no-cache-dir --no-deps --require-hashes -r requirements.locked-hashes.txt
/absolute/new/pyscenic39/bin/python -m pip check
```

The environment occupies approximately773MiB after installation; free disk was above the6GiB hold threshold throughout setup. No large motif-ranking resource was downloaded by this setup task.

Official references: [pySCENIC0.12.1 release](https://pypi.org/project/pyscenic/0.12.1/), [stable installation and CLI documentation](https://pyscenic.readthedocs.io/en/stable/installation.html).
