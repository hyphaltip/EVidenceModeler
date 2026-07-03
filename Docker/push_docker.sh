#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

VERSION="$(grep -m1 '^version' "${REPO_ROOT}/evm/Cargo.toml" | sed -E 's/.*"([^"]+)".*/\1/')"

IMAGE="${IMAGE:-evidencemodeler}"
ORG="${ORG:-stajichlab}"

docker push "${ORG}/${IMAGE}:${VERSION}"
docker push "${ORG}/${IMAGE}:latest"
