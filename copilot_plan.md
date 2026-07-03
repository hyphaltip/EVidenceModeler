# EVidenceModeler Rust Rewrite — Original Design Notes

> Reformatted from the initial Visual Studio / Copilot draft. This is the
> *design sketch* that seeded the rewrite. For the current, ground-truth status
> and the forward plan, see **[PLAN.md](PLAN.md)**.

## Overview

EVM is a ~3,500-line Perl bioinformatics pipeline that weights and combines gene
predictions with protein/transcript alignment evidence via a dynamic-programming
trellis, producing consensus gene models in GFF3. The Rust rewrite should
preserve exact algorithmic behavior while delivering faster parallel execution
and eliminating the Perl + ParaFly dependency.

## Crate Structure

A single Cargo workspace:

```
evm/
├── Cargo.toml      (workspace)
├── evm-core/       (library: all algorithms)
├── evm-cli/        (binary: replaces the EVidenceModeler script)
└── evm-utils/      (binaries: utility scripts)
```

## Phase 1 — Core Data Types (`evm-core/src/types/`)

- **genome.rs** — `GenomeSequence` (wraps `Vec<u8>` uppercase DNA, 1-based access,
  `reverse_complement()`, `rev_coord(len, pos)`, substring extraction);
  `FeatureVec` (parallel `Vec<u8>` holding START=1, DONOR=2, ACCEPTOR=3, STOP=4
  per 1-indexed position; populated by scanning ATG, GT/GC, AG, TAA/TGA/TAG);
  `MaskVec` (`Vec<bool>` for repeat-masked positions).
- **evidence.rs** — `EvClass` enum {Protein, Transcript, AbinitioPrediction,
  OtherPrediction}; `EvWeight` map (parsed from weights file); `EvidenceChain`
  {accession, target, ev_type, ev_class, lend, rend, links, gaps}.
- **exon.rs** — `ExonType` {Initial, Internal, Terminal, Single, Bound};
  `Orientation` {Fwd, Rev}; `ExonPhase` u8 in [1..=6]; `Exon` struct using an
  arena/index pool (`link: Option<usize>`) to avoid self-referential lifetimes.
- **prediction.rs** — `EvmPrediction`: vector of exon indices forming one gene
  model, with `is_eliminated`, mode (Standard/Intron), span.

## Phase 2 — File Parsers (`evm-core/src/io/`)

- **fasta.rs** — streaming buffered FASTA reader.
- **gff3.rs** — streaming GFF3 parser (9-col, `#` comments, ID/Parent/Target/Query).
- **weights.rs** — 3-column weights file.
- **partitions.rs** — partitions listing (accession, base_dir, Y|N, partition_dir).
- **evm_output.rs** — intermediate EVM text format (parse/write).

## Phase 3 — Core Algorithm (`evm-core/src/algo/`)

- **splice_sites.rs** — fast substring scan (memchr); `populate_genome_features`.
- **coding_scores.rs** — per-base coding score vector (PROTEIN + ABINITIO only).
- **introns.rs** — intron scoring/evidence; forward/reverse predicted intron vecs.
- **intergenic.rs** — evidence-weighted intergenic scoring; start/stop peak augmentation.
- **peaks.rs** — sliding-window peak detection.
- **phases.rs** — `determine_good_phases`, `is_stop_codon`.
- **trellis.rs** — DP trellis: `build_trellis`, `are_compatible_exons`,
  `score_boundary_condition`, `traverse_path`.
- **filter.rs** — low-support filtering; 5'-partial→complete conversion.
- **load_predictions.rs**, **load_evidence.rs**, **process.rs** — per-strand
  pipeline coordination.
- **consensus.rs** — recursive consensus prediction capped by re-search minimums.

## Phase 4 — Partitioning (`evm-core/src/partition/`)

`partition_files_based_on_contig`, `get_range_list`, `write_genome_partition`,
`partition_gff3_range`.

## Phase 5 — Recombination (`evm-core/src/recombine/`)

`parse_and_add_predictions`, `join_intronic_preds`, `combine_predictions`
(DP over partitions maximizing non-overlapping complete-gene coding length).

## Phase 6 — GFF3 Conversion (`evm-core/src/gff3_convert/`)

`evm_to_gff3` (gene/mRNA/CDS, `evm.TU.*`/`evm.model.*` IDs), `gff3_to_proteins`
(replaces Gene_obj.pm + Nuc_translator.pm), `gff3_to_bed`.

## Phase 7 — Translation (`evm-core/src/translate/`)

Standard codon table as a compile-time perfect-hash map; alternative codes
(e.g. Tetrahymena TGA) configurable via `stop_codons`; `translate`,
`reverse_complement`.

## Phase 8 — CLI Binaries

- **evm-cli/src/main.rs** — drop-in orchestrator: clap flags matching Perl,
  checkpoints dir, partition → per-partition rayon parallelism (replaces ParaFly)
  → recombine → convert → extract proteins/CDS/BED; `.ok` checkpoint files.
- **evm-utils/src/bin/** — `evidence_modeler` (single-partition shim),
  `partition_evm_inputs`, `recombine_evm_outputs`, `gff3_file_to_proteins`.

## Phase 9 — Testing Strategy

Unit tests per module; integration test diffing full-pipeline output against the
Perl version; golden-file fixtures; optional property tests.

## Key Implementation Decisions

| Design question | Choice |
|---|---|
| Exon object graph | Arena/index pool (`Vec<Exon>`, `link: Option<usize>`) |
| Parallelism | `rayon` (replaces ParaFly) |
| GFF3 parsing | Hand-written streaming parser |
| FASTA reading | Hand-written `BufReader` |
| CLI args | clap v4 derive, Perl-compatible flag names |
| Genetic code | `phf` compile-time perfect hash |
| Output format | Byte-for-byte compatible with Perl outputs |

## Dependencies

clap, rayon, memchr, log + env_logger, anyhow, phf, serde/serde_json (optional).

## Migration / Compatibility Notes

- `evm-cli` is a drop-in replacement: identical flags and output files
  (`*.EVM.gff3`, `*.EVM.pep`, `*.EVM.cds`, `*.EVM.bed`).
- Intermediate EVM text format kept identical for downstream tools.
- Partition directory structure (`.partitions/`, `.partitions.listing`) preserved.
- ParaFly (`plugins/`) removed once Rust handles parallelism natively.

## Approximate Complexity

| Module | Complexity | Notes |
|---|---|---|
| trellis.rs + consensus.rs | ★★★★★ | Core DP; ~1,200 lines of Perl logic |
| load_predictions.rs + load_evidence.rs | ★★★★ | GFF3 parsing + exon classification |
| filter.rs | ★★★ | Low-support filtering, partial→complete |
| intergenic.rs | ★★★ | Evidence-weighted intergenic scoring |
| gff3_to_proteins.rs | ★★★ | Replaces Gene_obj.pm + Nuc_translator.pm |
| partition.rs + recombine.rs | ★★★ | File I/O heavy |
| Parsers / CLI | ★★ | Mechanical translation |
