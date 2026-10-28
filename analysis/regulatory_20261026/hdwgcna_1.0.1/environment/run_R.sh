#!/bin/bash
set -euo pipefail
SATB2_REG_ENV="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
export R_LIBS_USER="$SATB2_REG_ENV/R-library"
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4
exec /usr/local/bin/Rscript --vanilla "$@"
