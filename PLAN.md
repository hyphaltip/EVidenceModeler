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

- **2026-06-27 — Phase B underway.** Golden reference captured at
  `evm/tests/fixtures/Contig1.perl.evm.out` (11 genes). Two algorithmic fixes
  landed, each cross-checked against the Perl source:
  1. **Reverse-strand transpose** (`transpose_exons_back_to_forward_strand`
     parity): the merge flipped coordinates but not `orientation`/reading frame,
     so the trellis never started a reverse gene. Now sets `orientation = Rev`
     and remaps frames 1→4/2→5/3→6. Brought back all reverse-strand genes.
  2. **Intergenic scoring** (`populate_intergenic_regions` parity): was inverting
     the coding-score vector; now sums each ABINITIO type's weight across the
     gaps between that type's neighbouring genes. Removed two spurious short
     genes caused by wrong noncoding scores.

  **Current parity on Contig1: 7 of 11 Perl genes match exactly** (was 0 reverse
  + spurious before). Remaining diffs — the Phase B3 work-list:
  - Missing: single-exon reverse gene `3632-4546` (`single-`, backed by 3
    complete ab-initio CDS predictions) — produces nothing in that region.
  - Terminal-coordinate diffs: `44377-50459` vs Perl `-50843`; `57662-59941`
    vs `57371-`; `61745-63283` vs `-63134`. All terminal-exon / start-stop-peak
    boundary selection.
  - **Output text format** still differs from Perl (`#EVM prediction mode:...`
    vs `# EVM prediction: Mode:... S-ratio ... orient ... score ...`; exon-row
    columns and evidence formatting differ) — Phase B2, not yet started.
  - Known still-divergent internals to verify next: reverse-strand coding /
    intron vectors are merged by **raw index without transposing**
    (`evidence_modeler.rs` / `main.rs` merge loops); Perl builds the fwd/rev
    intron vectors once at the end from the accumulated `%PREDICTED_INTRONS`.

- **2026-06-27 — Phase B investigation deepened (roadmap for the remaining
  diffs).** Mapped the exact Perl logic and Rust gaps for the remaining 4
  model diffs + format. Key findings:
  - **Prediction scoring is unimplemented.** `trellis::traverse_path` builds
    `EvmPrediction`s but never sets `total_score` (stays 0.0). Perl's
    `prediction_score` (= `EVM_prediction::get_score`, set to the trellis path
    total) drives both the header `score(...)` and the filter ratio. This must
    be implemented first — everything downstream depends on it.
  - **The low-support filter is a heuristic, not the Perl formula.** Port
    `filter_predictions_low_support` (Perl lines 3436-3550) faithfully:
    `noncoding = calc_intergenic_score(lend,rend)`;
    `noncoding_intron_addition = Σ FWD_PRED_INTRON_VEC[i]+REV_PRED_INTRON_VEC[i]`
    over span; per-intron `offset += calc_intergenic_score(intron) + Σ
    (strand_intron_vec[i] − predicted_intron_contrib/intron_len)`;
    `raw_noncoding = noncoding + intron_addition`;
    `noncoding_equivalent = max(raw_noncoding − offset, 0.0001·score)`;
    `score_ratio = score/noncoding_equivalent`; eliminate if
    `ratio < 0.75` (`MIN_CODING_NONCODING_SCORE_RATIO`) or
    `coding_length < (STANDARD?150:300)`. Store raw_noncoding/offset/
    noncoding_equivalent/score_ratio on the prediction for the header.
    Needs merged `predicted_introns` (the shim currently does NOT merge the
    reverse strand's `predicted_introns`).
  - **Missing single-exon gene `3632-4546` is a TRELLIS-CREATION issue, not a
    filter one** — it does not appear even with `--report-elm`, so it is never
    produced as a candidate path. Investigate single-exon (`single-`)
    candidate creation in `load_predictions` (reverse pass) and/or its
    selection/seeding in `build_trellis` (isolated single exon between two
    multi-exon genes via intergenic transitions). The 3 inputs are complete
    single-CDS reverse predictions at 3632-4546.
  - **Three terminal-boundary diffs** (`44377-50459` vs `-50843`,
    `57662-` vs `57371-`, `-63283` vs `-63134`): terminal-exon selection /
    start-stop-peak augmentation in trellis/consensus.
  - **Output format spec (Phase B2)** is fully nailed down from the Perl
    `toString` methods (Exon @4847, EVM_prediction @5112):
    - per-recursion: `!! Predictions spanning range L - R [R<n>]` (only when a
      pred survives filter or `--report-elm`).
    - header: `# EVM prediction: Mode:<m> S-ratio: <ratio> <lend>-<rend>
      orient(<o>) score(%.2f) noncoding_equivalent(%.2f) raw_noncoding(%.2f)
      offset(%.2f) ` (+` *** ELIMINATED *** ` if eliminated).
    - exon row: `end5\tend3\t<type><orient>\t<startFrame>\t<endFrame>\t` then
      `{acc;class},`-joined evidence (trailing comma chopped). Frames are 1-6,
      NOT gff 0-2. Type carries the orient suffix (`initial+`,`single-`,…).
    - intron row: `e5\te3\tINTRON\t\t\t{acc;type},…`; e5/e3 swapped for `-`.
    - components (exons+introns) sorted ascending by first coord; blank line
      after each prediction.
    Format depends on the score fields above, so do scoring+filter first.

## 0a. Live parity status (testing/Contig1)

**10 of 11 Perl genes match exactly.** Low-support filter ported faithfully and
verified against the golden header values: `raw_noncoding` matches ~exactly
(e.g. 3619.01 vs 3619.02; 6165.00, 4113.00 exact), S-ratios within ~0.5%.

**Score residuals → one remaining missing step:** `prediction_score` runs
slightly high and `offset` slightly low across genes. Traced to the still-missing
`decrement_coding_using_protein_alignment_introns` (Perl line 4206; constants
`INTRON_MEDIAN_FACTOR=2`, `MIN_ALIGNMENT_GAP_SIZE_INFER_INTRON`): it subtracts the
protein weight from coding scores over protein-alignment gaps (inferred introns)
shorter than `2×median_gap`. Needed for byte-exact scores; the `gaps` are already
on `EvidenceChain`. NOTE: this is NOT the 842 fix (it lowers `base(842)`, shifting
toward 1018 — opposite of what's needed).

**The 842 gene (1018 vs 842) — fully root-caused, fix not yet found.** Both
candidates are correctly created: genemark `initial+ 842-1127` (genemark-only
evidence) and `internal 1018-1127` (from gap2.2, a legit 12-block EST chain with
a real acceptor at 1016). It's a pure base-score tie-break the trellis resolves
to 1018 in Rust but 842 in Perl. base() should compute identically in both, so
the definitive next diagnostic is to dump Perl's actual exon base_scores for
842-1127 vs 1018-1127 (run `evidence_modeler.pl` with `$DEBUG`/`$SEE`, or build
ParaFly) and compare to Rust — that pinpoints which scoring term diverges.

 Fixes landed since the audit: byte-exact
`analyze_peaks`; faithful `augment_intergenic_from_start/stop_peaks` (peaks now
carry strand); internal-exon recovery in `recover_partial_prediction`.

**The one remaining span diff** is the first gene: Rust `1018-3150` vs Perl
`842-3150`. Root-caused: the genemark `initial+ 842-1127` exon IS created
correctly (verified: start+donor present, good_phases [1,2]); but a competing
`1018-1127` *internal* exon built from high-weight transcript evidence
(`alignAssembly`, weight 10) has a larger base score, so the trellis starts the
gene there (a 5'-partial). Perl prefers 842. The divergence is in how
transcript-alignment evidence creates/weights competing internal exons
(`load_evidence` / `instantiate_evidence_based_exons` and which exon accrues the
transcript per-base contribution in `score_exons`) — needs a `load_evidence`
deep-dive vs the Perl. NOT a creation/peak/trellis-structure bug.

## 0b. Approximation audit (directive: EXACT port first, no heuristics)

Governing rule: byte-for-byte parity with the Perl is the bar; no approximation
may stand in for a real algorithm. Status of every known approximation/gap:

**Confirmed faithful (verified vs Perl):** `score_exons`; `populate_intergenic_scores`
(prediction-span method); reverse transpose; prediction `finalize` scoring;
trellis boundary placement + `score_boundary_condition`; `are_compatible_exons`.

**Confirmed approximate / divergent — to fix in dependency order:**
1. `augment_intergenic_from_peaks` — crude ±500 smear, currently DISABLED; port
   Perl `augment_intergenic_from_start/stop_peaks` (~360 lines; needs
   `find_closest_exon_within_range`, `START_STOP_RANGE`=500). Fixes the 2 remaining
   5' diffs (1018 vs 842, 12430 vs 14459). **ACTIVE NEXT.**
2. `filter_predictions_low_support` — heuristic, not Perl `noncoding_equivalent`/
   `score_ratio` (lines 3436-3550). Unblocked (total_score exists).
3. `convert_5prime_partials` — stubbed TODO, does nothing.
4. `format_prediction` — output text not matching Perl `toString` (B2).
5. `decrement_coding_using_protein_alignment_introns` — MISSING entirely
   (Perl main flow line 726).

**Not yet verified (audit before claiming parity):** `analyze_peaks` (feeds #1);
`recombine` DP (`combine_predictions`/`join_intronic_preds`);
`recover_partial_prediction`; cross-strand intron-vector merge; `gff3_to_proteins`
translation + ID assignment; `partition::get_range_list`; PASA
`supplement_terminal_exons` (only with --terminal_exons).

**Mechanisms to guarantee exactness (to add):**
- Byte-level golden test: assert Rust `evm.out` == Perl `evm.out` (fixture exists).
- Per-function cross-check: mirror Perl `$DEBUG` intermediate dumps (intergenic,
  intron vecs, peaks, exons, trellis) and diff each vector at the source.

**Phase G (post-parity, NOT part of the port):** once bit-for-bit faithful, use
the fast Rust EVM as a trustworthy baseline to curate gold-standard annotation
and re-tune/retrain scoring — with the faithful port as the regression anchor.

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
