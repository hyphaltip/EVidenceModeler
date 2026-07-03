#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# 1. Build the Rust workspace (release) and install the eight binaries.
# ---------------------------------------------------------------------------
export CARGO_HOME="${CARGO_HOME:-${SRC_DIR}/.cargo}"
export CARGO_NET_OFFLINE="${CARGO_NET_OFFLINE:-false}"

cargo build --release --locked --manifest-path evm/Cargo.toml

BINS=(
  EVidenceModeler
  evidence_modeler
  partition_evm_inputs
  recombine_evm_outputs
  convert_EVM_outputs_to_GFF3
  gff3_file_to_proteins
  gff3_gene_prediction_file_validator
  augustus_to_evm_gff3
)

install -d "${PREFIX}/bin"
for b in "${BINS[@]}"; do
  install -m 0755 "evm/target/release/${b}" "${PREFIX}/bin/${b}"
done

# ---------------------------------------------------------------------------
# 2. Install the legacy Perl/Python helper converters + Perl libraries.
#    These are not on PATH by default; EVM_HOME points callers at them.
# ---------------------------------------------------------------------------
SHARE="${PREFIX}/opt/evidencemodeler"
install -d "${SHARE}"
cp -a EvmUtils PerlLib "${SHARE}/"

# ---------------------------------------------------------------------------
# 3. Activation scripts: export EVM_HOME so the Perl helpers and any tool
#    that expects $EVM_HOME (e.g. funannotate Perl fallback) can find them.
# ---------------------------------------------------------------------------
for CHANGE in activate deactivate; do
  install -d "${PREFIX}/etc/conda/${CHANGE}.d"
done

cat > "${PREFIX}/etc/conda/activate.d/evidencemodeler.sh" <<'EOF'
export EVM_HOME="${CONDA_PREFIX}/opt/evidencemodeler"
EOF

cat > "${PREFIX}/etc/conda/deactivate.d/evidencemodeler.sh" <<'EOF'
unset EVM_HOME
EOF
