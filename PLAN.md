# PLAN.md — Completing the EVidenceModeler Rust Rewrite

Status date: 2026-06-27 · Branch: `rust-rewrite-completion` (forked from
`copilot/rewrite-in-rust`)

This plan continues the Perl → Rust port started in `copilot_plan.md`. It is
written against the **actual measured state** of the existing Rust code, not the
original design sketch.

---

## 0. Progress log

- **2026-06-27 — Phase A complete.** Repo hygiene done (`.gitignore` added,
  1,736 `evm/target/` artifacts untracked). Fixed both startup/runtime panics:
  the clap `version` collision in `evm-cli` (now uses clap's built-in
  `version = VERSION`, `--version`/`--help` work) and the u32 coordinate
  underflows in `algo/` (saturating subtraction at the `genome_features.get(...)`
  sites in `load_evidence.rs` and `load_predictions.rs`, and `end3 - 3`). All
  compiler warnings cleared (0 warnings). `cargo test` = 22/22. The
  `evidence_modeler` single-partition shim now runs the `testing/` data to
  completion (exit 0). **Next: Phase B** — its `evm.out` (60 lines) does not yet
  match the Perl reference (170 lines); both the per-line format and the gene
  models differ. That diff is the Phase B work-list.

## 1. Where we actually are

The Rust workspace under `evm/` is **fully scaffolded** and matches the original
design. Concretely, verified on this branch:

- **~4,864 lines of Rust** across the full module tree (`types/`, `io/`, `algo/`,
  `partition/`, `recombine/`, `gff3_convert/`, `translate/`, plus `evm-cli` and
  `evm-utils`). Only **one** `TODO` remains in the source
  (`filter.rs:99` — partial→complete conversion).
- **`cargo build` succeeds** (warnings only — unused vars/`mut`, stray parens).
- **`cargo test` passes: 22/22 unit tests.**
- A Perl reference partition output (`evm.out`, 170 lines) can be generated from
  the `testing/` dataset and is our parity target.

So the skeleton is done. The remaining work is **making it run end-to-end,
proving correctness against Perl, wiring up the orchestrator, and optimizing.**

### Blocking defects found while validating (must fix first)

1. **`evm-cli` (the `EVidenceModeler` binary) panics on startup.** clap aborts
   with *"Argument names must be unique, but 'version' is in use by more than one
   argument"*. Cause: `#[command(... version ...)]` auto-generates a `--version`
   flag **and** there is a manual `version: bool` field (`evm-cli/src/main.rs:24`
   and `:126-128`). The main orchestrator currently cannot run at all.
2. **`evidence_modeler` (single-partition shim) panics on the real test data.**
   `attempt to subtract with overflow` at `load_evidence.rs:142`
   (`genome_features.get((end5 - 2) as usize)` underflows when `end5 < 2`). The
   core algorithm aborts before producing any output on `testing/`.
3. **Repo hygiene: 1,736 build artifacts under `evm/target/` are committed** and
   there is **no `.gitignore`**. This bloats the branch (multi-MB `.rlib`/`.o`
   blobs) and must be removed before any further commits.

### Known unverified risk areas

- **trellis.rs + consensus.rs (456 Rust lines) port ~1,200 lines of Perl DP.**
  This is the highest-risk translation and is completely untested against Perl
  output. Compiling ≠ correct.
- No **golden-file / integration test** exists comparing Rust output to Perl.
- The **orchestrator** (`evm-cli`) parallelism, checkpointing (`.ok` files),
  recombine, and GFF3/protein/CDS/BED emission paths have never been exercised
  end-to-end (blocked by defect #1).

---

## 2. Goal & definition of done

A single `EVidenceModeler` Rust binary that is a **drop-in replacement** for the
Perl pipeline:

- Identical CLI flags.
- **Byte-for-byte identical** `*.EVM.gff3`, `*.EVM.pep`, `*.EVM.cds`, `*.EVM.bed`
  on the `testing/` dataset (and ideally on at least one larger real genome).
- Native parallelism via `rayon` — no Perl, no ParaFly.
- Preserves the intermediate `evm.out` text format and `.partitions*` layout for
  resume-ability and downstream-tool compatibility.

"Done" = the golden-file integration test (Section 5) passes, and the binary runs
a real multi-contig genome to completion faster than the Perl+ParaFly baseline.

---

## 3. Workstream phases

Phases A–C are sequential (each unblocks the next). D–F can overlap once C is green.

### Phase A — Make it build clean & run (hygiene + the 2 panics)

- A1. Add `evm/.gitignore` (`/target`, `Cargo.lock` policy: keep for binaries),
  `git rm -r --cached evm/target`, commit the cleanup.
- A2. Fix the clap `version` collision in `evm-cli/src/main.rs` (drop the manual
  `version` field, or `disable_version_flag(true)`; keep behavior identical to
  Perl's `--version`). Add a smoke test that `--help` and `--version` don't panic.
- A3. Fix the coordinate underflow in `load_evidence.rs:142` and audit **all**
  `end5 - N` / `pos - N` / `... - 1` expressions in `algo/` for u32 underflow
  (use `checked_sub`/saturating arithmetic, or 1-based bounds guards mirroring
  the Perl). This is a class of bug, not a single line — sweep `load_evidence`,
  `load_predictions`, `phases`, `splice_sites`, `intergenic`.
- A4. Get `evidence_modeler` to run the `testing/` Contig1 partition to completion
  without panicking and emit an `evm.out`.
- A5. Clean up compiler warnings (`cargo build` warning-free; add
  `#![warn(...)]` discipline later).

**Exit:** `evidence_modeler` produces an `evm.out` on the test data; `evm-cli
--help/--version` work.

### Phase B — Single-partition correctness vs Perl (the hard part)

This is where algorithmic fidelity is won or lost.

- B1. Generate the Perl golden `evm.out` for Contig1 (already reproducible via
  `EvmUtils/evidence_modeler.pl` — see Section 5) and check it in as a fixture.
- B2. Diff Rust `evm.out` against Perl `evm.out`. Expect divergence. Drive the
  diff to zero, working module by module. Likely suspects in priority order:
  `trellis` (compatibility/boundary scoring, exon ordering), `consensus`
  (recursion limits, re-search regions), `intergenic` (peak augmentation),
  `phases` (frame math), `filter` (low-support ratio, the `filter.rs:99` TODO).
- B3. Add **focused unit tests** capturing the Perl semantics for each fixed
  discrepancy (so regressions are caught): splice-site scanning, phase
  determination, intron scoring, peak detection, exon compatibility, one full
  trellis traversal on a tiny hand-built case.
- B4. Verify both strands (`--forward-strand-only` / `--reverse-strand-only`) and
  the alternative stop-codon path.

**Exit:** Rust `evm.out` is byte-identical to Perl `evm.out` on Contig1 (both
strands).

### Phase C — Orchestrator & format parity (`evm-cli` end-to-end)

- C1. Partitioning: `partition_evm_inputs` must produce the same
  `.partitions.listing` and per-contig directory layout as
  `partition_EVM_inputs.pl` (test on `testing/`, and on a synthetic multi-contig
  + multi-partition genome to exercise `segmentSize`/`overlapSize` windowing).
- C2. Recombination: `recombine_evm_outputs` must match
  `recombine_EVM_partial_outputs.pl` (intronic-pred joining + DP combine across
  overlapping partitions).
- C3. GFF3/protein/CDS/BED conversion: match
  `convert_EVM_outputs_to_GFF3.pl` + `gff3_file_to_proteins.pl` +
  `gene_gff3_to_bed.pl`, including `evm.TU.*` / `evm.model.*` ID assignment and
  CDS phases.
- C4. Wire `evm-cli`: clap flags matching Perl exactly, checkpoint `.ok` files
  for resume, `rayon` per-partition parallelism replacing ParaFly.
- C5. Full-pipeline golden test: `EVidenceModeler ...` on `testing/` produces
  `smalltest.EVM.{gff3,pep,cds,bed}` byte-identical to Perl.

**Exit:** `make -C testing test` (re-pointed at the Rust binary) reproduces Perl
outputs exactly.

### Phase D — Utility-script parity (`evm-utils`)

Port/validate the utilities that downstream users invoke directly. Prioritize the
ones the pipeline itself depends on (already covered in C) then the standalone
converters under `EvmUtils/` and `EvmUtils/misc/` on a best-effort, test-backed
basis. Many `misc/` converters are thin format adapters and can remain Perl/Python
shims initially — **explicitly document** which are ported vs. retained.

### Phase E — Optimization

Only after correctness is locked (golden tests green):

- E1. Establish a benchmark harness (`criterion` for hot functions; wall-clock +
  peak-RSS for the full pipeline on a real genome). Record a Perl+ParaFly
  baseline.
- E2. Profile (`perf` / `cargo flamegraph`). Expected hot spots: trellis DP,
  splice-site scanning, per-base coding/intergenic score vectors.
- E3. Targeted wins: `memchr`/two-way for motif scans, avoid per-exon `String`
  keys in `exons_via_coords` (use packed integer keys), reuse score-vector
  allocations, bound `max_prev_exons_compare` sensibly, `rayon` across contigs
  and partitions. Re-run golden tests after every optimization.
- E4. `release` profile tuning (`lto = "thin"`, `codegen-units = 1`,
  `panic = "abort"` for the binary).

**Target:** correct output, and meaningfully faster than Perl+ParaFly at equal
core count, with lower memory.

### Phase F — Packaging, CI, migration

- F1. CI (GitHub Actions): `cargo fmt --check`, `cargo clippy -D warnings`,
  `cargo test`, and the golden integration test on every PR.
- F2. Update `Docker/`, `Makefile`, and `testing/runMe*.sh` to build/use the Rust
  binary; drop the ParaFly submodule once C4 lands.
- F3. Update `README.md` / `Changelog.txt`; document the migration and any
  intentional behavioral differences (there should be none in output).
- F4. Version the binary to track the Perl release line (currently 2.1.0).

---

## 4. Repository / branch hygiene (do immediately, Phase A0)

1. Add `.gitignore` (`evm/target/`, editor cruft).
2. `git rm -r --cached evm/target` — stop tracking 1,736 build artifacts.
3. Decide `Cargo.lock` policy: **commit it** (this ships binaries).
4. Keep `copilot_plan.md` (design history) and this `PLAN.md` (forward plan).

---

## 5. Testing strategy (the backbone of this rewrite)

**Golden-file parity is the primary correctness oracle.** Because the goal is a
byte-for-byte drop-in replacement, every phase is gated on matching Perl output.

- **Reference generation.** ParaFly is unbuilt, so the *full* Perl `EVidenceModeler`
  wrapper cannot run as-is, but the pieces can:
  - Partition: `EvmUtils/partition_EVM_inputs.pl` (runs).
  - Per-partition: run the command emitted in `*.partitions.evm_cmds`, i.e.
    `EvmUtils/evidence_modeler.pl -G genome.fasta -g gene_predictions.gff3 -w
    weights.txt -e transcript_alignments.gff3 -p protein_alignments.gff3
    --min_intron_length 20 --terminal_intergenic_re_search 10000` → `evm.out`
    (✅ verified working; 170-line output).
  - Recombine/convert: the respective `EvmUtils/*.pl` scripts.
  - Alternatively build the ParaFly submodule once to get an end-to-end Perl run.
- **Fixtures.** Check Perl outputs into `evm/tests/fixtures/` (`evm.out`,
  `*.EVM.gff3/pep/cds/bed`, `.partitions.listing`).
- **Integration test** (`evm/tests/golden.rs`): run each Rust stage on `testing/`
  and assert byte-equality with fixtures. Normalize only provably irrelevant
  fields (absolute paths, timestamps) and document each normalization.
- **Unit tests:** keep the existing 22; add one per discrepancy fixed in Phase B
  and per parser edge case. Aim for coverage of every `algo/` function with a
  Perl-derived expected value.
- **Property tests (optional, `proptest`):** coordinate round-trips
  (`rev_coord`), translation/reverse-complement invariants, partition window
  coverage (no gaps/overlaps beyond `overlapSize`).
- **Larger real genome:** before declaring done, run one real annotated genome
  through both pipelines and diff — `testing/` is a single contig and won't
  exercise multi-partition recombination.

---

## 6. Risks & mitigations

| Risk | Mitigation |
|---|---|
| Trellis/consensus DP subtly diverges from Perl | Golden diff at the `evm.out` level *before* touching the orchestrator; tiny hand-built trellis unit tests |
| u32 underflow / off-by-one from 1-based Perl coords | Phase-A sweep of all subtractions; bounds-guard helpers mirroring Perl |
| Floating-point score formatting differs (e.g. `%.2f`) | Match Perl's `printf` formatting exactly in `evm_output`/GFF3 writers |
| HashMap iteration-order nondeterminism affecting output | Sort before emitting; never rely on map order for output |
| `testing/` too small to catch multi-partition bugs | Add synthetic + one real multi-contig genome to the test matrix |
| Committed `target/` corrupts/bloats history | Remove in Phase A0 before any new commits |

---

## 7. Suggested milestones / sequencing

1. **M0 (hours):** Phase A0 hygiene + A1–A2 (build clean, CLI doesn't panic).
2. **M1 (days):** Phase A3–A4 — `evidence_modeler` runs `testing/` to completion.
3. **M2 (the big one):** Phase B — `evm.out` byte-identical to Perl, both strands.
4. **M3:** Phase C — full `evm-cli` pipeline reproduces all Perl outputs.
5. **M4:** Phase E — benchmarked, optimized, still byte-identical.
6. **M5:** Phase F — CI + Docker + docs; remove ParaFly.

---

## 8. Immediate next actions (proposed)

1. Phase A0 repo hygiene (`.gitignore` + untrack `target/`).
2. Fix the two panics (clap `version`, `load_evidence` underflow).
3. Run `evidence_modeler` on `testing/` and capture the first Rust-vs-Perl
   `evm.out` diff — that diff defines the real scope of Phase B.

> Awaiting your go-ahead on scope/sequencing before starting implementation.
