#!/bin/bash
set -euo pipefail

# Build the Rust EVidenceModeler Docker image.
# Run from anywhere; the build context is the repository root.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Version comes from the Rust workspace manifest (single source of truth).
VERSION="$(grep -m1 '^version' "${REPO_ROOT}/evm/Cargo.toml" | sed -E 's/.*"([^"]+)".*/\1/')"

IMAGE="${IMAGE:-evidencemodeler}"
ORG="${ORG:-stajichlab}"

echo "Building ${ORG}/${IMAGE}:${VERSION} (and :latest)"
docker build \
    -f "${SCRIPT_DIR}/Dockerfile" \
    --build-arg EVM_VERSION="${VERSION}" \
    -t "${ORG}/${IMAGE}:${VERSION}" \
    -t "${ORG}/${IMAGE}:latest" \
    "${REPO_ROOT}"
