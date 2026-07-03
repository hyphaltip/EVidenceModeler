#!/bin/bash
set -euo pipefail

# Build the Singularity/Apptainer images. Must run from the repository root
# because the .def %files sections reference paths relative to the build cwd.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

# Prefer apptainer, fall back to singularity.
BUILDER="$(command -v apptainer || command -v singularity || true)"
[[ -n "${BUILDER}" ]] || { echo "apptainer/singularity not found" >&2; exit 1; }

TARGET="${1:-mariadb}"   # "mariadb" (default) or "plain"

case "${TARGET}" in
  mariadb) DEF="Singularity/evm-mariadb.def"; SIF="evm-mariadb.sif" ;;
  plain)   DEF="Singularity/evm.def";         SIF="evm.sif" ;;
  *) echo "usage: $0 [mariadb|plain]" >&2; exit 1 ;;
esac

echo "Building ${SIF} from ${DEF} using ${BUILDER}"
"${BUILDER}" build "${SIF}" "${DEF}"
echo "Done: ${REPO_ROOT}/${SIF}"
