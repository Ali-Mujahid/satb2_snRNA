# Isolated Scanpy runtime setup

Required for the user-authorized Scanpy/Seurat comparison. New project-only .venv_scanpy; existing Python and R environments remain untouched. CPython3.12.14 arm64. Selected wheel-only dependencies and hashes are frozen in scanpy_requirements_lock.txt from the successful pip dry resolution.

Direct requirements:

- scanpy==1.12.4
- anndata==0.13.4
- scikit-misc==0.5.3
- igraph==1.0.0
- leidenalg==0.12.0
- harmonypy==2.1.0
- pandas==2.3.3
- h5py==3.16.0
- psutil==7.2.2

Budget: approximately2GB including wheel cache; run outputs separately budgeted20GB. Package/import verification only; no scientific execution. Harmony2 Python/R generation is matched.
