# syntax=docker/dockerfile:1
#
# Multi-stage build for the Rust EVidenceModeler.
#   Stage 1 (builder): compiles the eight release binaries.
#   Stage 2 (runtime): minimal image with only the runtime deps needed by the
#                      Perl/Python helper converters (the Rust bins need none).
#
# Build context must be the repository root:
#   docker build -f Docker/Dockerfile -t evidencemodeler:2.5.0 .

ARG RUST_VERSION=1.83
ARG EVM_VERSION=2.5.0

# --------------------------------------------------------------------------- #
# Stage 1: builder
# --------------------------------------------------------------------------- #
FROM rust:${RUST_VERSION}-slim-bookworm AS builder

WORKDIR /src
# Copy only what the Rust build needs first, to maximize layer caching.
COPY evm/ evm/
RUN cargo build --release --locked --manifest-path evm/Cargo.toml

# --------------------------------------------------------------------------- #
# Stage 2: runtime
# --------------------------------------------------------------------------- #
FROM debian:bookworm-slim AS runtime

ARG EVM_VERSION
LABEL org.opencontainers.image.title="EVidenceModeler (Rust)" \
      org.opencontainers.image.description="Consensus gene structure annotation from multiple evidence sources" \
      org.opencontainers.image.version="${EVM_VERSION}" \
      org.opencontainers.image.source="https://github.com/EVidenceModeler/EVidenceModeler" \
      org.opencontainers.image.licenses="BSD-3-Clause"

# Runtime deps for the Perl/Python helper converters only.
#   perl            -> ships DB_File as a core module
#   liburi-perl     -> URI::Escape
#   libbio-perl-perl-> Bio::SearchIO (BPbtab / summarize_btab_tophits)
#   python3         -> miniprot_GFF_2_EVM_*.py
RUN apt-get update && apt-get install -y --no-install-recommends \
        perl \
        liburi-perl \
        libbio-perl-perl \
        python3 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ENV EVM_HOME=/opt/evidencemodeler
ENV PATH="/opt/evidencemodeler/EvmUtils:${PATH}"

# Rust binaries -> /usr/local/bin
COPY --from=builder \
    /src/evm/target/release/EVidenceModeler \
    /src/evm/target/release/evidence_modeler \
    /src/evm/target/release/partition_evm_inputs \
    /src/evm/target/release/recombine_evm_outputs \
    /src/evm/target/release/convert_EVM_outputs_to_GFF3 \
    /src/evm/target/release/gff3_file_to_proteins \
    /src/evm/target/release/gff3_gene_prediction_file_validator \
    /src/evm/target/release/augustus_to_evm_gff3 \
    /usr/local/bin/

# Perl/Python helper converters + Perl libraries
COPY EvmUtils/ /opt/evidencemodeler/EvmUtils/
COPY PerlLib/  /opt/evidencemodeler/PerlLib/

# Non-root user for running the tools
RUN useradd -m -s /bin/bash evm
USER evm
WORKDIR /home/evm

CMD ["EVidenceModeler", "--help"]
