# PLAN.md — Completing the EVidenceModeler Rust Rewrite

Status date: 2026-06-27 · Branch: `rust-rewrite-completion` (forked from
`copilot/rewrite-in-rust`)

This plan continues the Perl → Rust port started in `copilot_plan.md`. It is
written against the **actual measured state** of the existing Rust code, not the
original design sketch.

---

## 00b. SESSION HANDOFF (2026-06-28) — Phase C END-TO-END PARITY — START HERE

**Branch `rust-rewrite-completion`. Build clean (0 warnings), `cargo test` 30/30
(26 unit + 4 golden integration).**

**Phase C is DONE on BOTH single-partition AND multi-partition `testing/Contig1`.**
The full `EVidenceModeler` Rust orchestrator runs end-to-end and its outputs match
the Perl golden in both modes:
- `smalltest.EVM.gff3` — **byte-identical**.
- `smalltest.EVM.bed` — **byte-identical**.
- `smalltest.EVM.pep` / `.cds` — **identical per-record** (Perl emits FASTA records
  in hash order = non-deterministic; oracle sorts records by header. Headers AND
  sequences byte-exact).
- `smalltest.partitions.listing` + all partitioned input files — byte-identical.

**MULTI-PARTITION verified (commit 7125c88):** ran `--segmentSize 30000
--overlapSize 10000` on Contig1 (63304 bp → 3 partitions: 1-30000, 20001-50000,
40001-63304). Per-partition evm.out (83/92/65 lines) → recombined 169 lines / 11
genes, identical to Perl golden; all four EVM outputs byte-identical (pep/cds modulo
FASTA record order). Recombine/partition bugs fixed: (1) partition dir naming
`_{lend}-{rend}` hyphen (Perl `${accession}_${lend}-${rend}`, recombine keys on
`/(\d+)-(\d+)$/`); (2) `process_prediction_text` was DROPPING the exon-type column
(`cols[2]`, e.g. `initial+`/`INTRON`) when rejoining rows → shifted every row left,
emptied GFF3 — Perl overwrites only x[0]/x[1] and rejoins ALL columns; (3) header
coordspan offset at whitespace token index 6 (was matching `S-ratio:`); (4) keep
INTRON rows; (5) blank-line separator between preds; (6) `extract_partition_lend`
mirrors `/(\d+)-(\d+)$/`. Added 2 recombine unit tests.

### What landed this session
1. **CLI flag parity (C4):** `evm-cli` clap now uses Perl's exact long names
   (`--sample_id`, `--gene_predictions`, `--segmentSize`, `--CPU`, …) via explicit
   `long = "..."`. Was kebab-case (`--sample-id`) → not a drop-in. Wired
   `--trellis_search_limit` → `max_prev_exons_compare`; added parity flags
   (`stitch_ends`, `exec_dir`, `limit_range_*`, etc.) accepted (some unused).
2. **GFF3 converter rewrite (`evm_to_gff3.rs`)** — faithful port of `EVM_to_GFF3.pl`
   + `Gene_obj::to_GFF3_format`: emits `exon` rows (were missing); correct CDS phase
   = compose `%phase_conversion` (1→0,2→1,3→2,4→0,5→1,6→2) THEN GFF3 swap 1↔2 →
   net `{1,4}→0,{2,5}→2,{3,6}→1`; `Name=EVM%20prediction%20X` (uri-escaped); exon
   order 5'→3' (asc '+' / desc '-'); model_id increments on blank lines.
3. **BED converter rewrite (`gff3_to_bed.rs`)** — `Gene_obj::to_BED_format`: name =
   `ID=<model>;<TU>;<com_name>` spaces→`_`; thick = CDS span; itemRgb `0`; blocks
   from exon features asc; final `sort -k1,1 -k2,2g -k3,3g`.
4. **PEP/CDS (`gff3_to_proteins.rs`):** header now
   `>{model} {TU}  {com_name} {contig}:lend-rend(orient)` (double space from empty
   locus_string); **protein translation trims the transcription-first CDS exon's
   GFF3 phase** (5'-partial models like gene 11, first-exon phase 1, were
   mistranslated `RLH*`). CDS/cDNA keep full sequence. Added `uri_unescape` in
   `gff3_convert/mod.rs` (decodes gff3 `Name`).
5. **De-duplicated the per-partition driver:** new
   `evm_core::pipeline::run_single_partition` is the single shared impl; both the
   `evidence_modeler` shim and the orchestrator's `run_evm_on_partition` call it
   (removed the ~120-line copy the old PLAN flagged as a drift risk).
6. **Golden integration test** `evm-core/tests/golden_convert.rs` (4 tests) +
   fixtures `evm/tests/fixtures/Contig1.perl.EVM.{gff3,bed,pep,cds}` +
   `Contig1.genome.fasta`.

### Reproduce the Perl golden (ParaFly unbuilt → run stages by hand)
In a scratch dir with the 5 `testing/` inputs:
```
perl <repo>/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta \
  --gene_predictions gene_predictions.gff3 --protein_alignments protein_alignments.gff3 \
  --transcript_alignments transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 \
  --partition_dir smalltest.partitions --partition_listing smalltest.partitions.listing
# (run evidence_modeler.pl in smalltest.partitions/Contig1 → evm.out, see §00)
perl <repo>/EvmUtils/recombine_EVM_partial_outputs.pl --partitions smalltest.partitions.listing --output_file_name evm.out
perl <repo>/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions smalltest.partitions.listing --output evm.out --genome genome.fasta
find ./smalltest.partitions -regex ".*evm.out.gff3" -exec cat {} \; > smalltest.EVM.gff3
perl <repo>/EvmUtils/gff3_file_to_proteins.pl smalltest.EVM.gff3 genome.fasta prot > smalltest.EVM.pep
perl <repo>/EvmUtils/gff3_file_to_proteins.pl smalltest.EVM.gff3 genome.fasta CDS  > smalltest.EVM.cds
bash -c "set -eou pipefail && perl <repo>/EvmUtils/gene_gff3_to_bed.pl smalltest.EVM.gff3 | sort -k1,1 -k2,2g -k3,3g > smalltest.EVM.bed"
```
Rust end-to-end: `EVidenceModeler --sample_id smalltest --genome genome.fasta
--weights ./weights.txt --gene_predictions … --segmentSize 100000 --overlapSize 10000`.

### What is NOT yet covered (next session)
- **Multi-partition: DONE** (commit 7125c88) for a single contig that splits into
  3 partitions — recombine DP (`join_intronic_preds` + `combine_predictions`) and
  the partition-`lend` offset are now exercised and byte-exact. STILL untested:
  **multi-CONTIG** (>1 record in the genome FASTA — `concatenate_gff3_outputs`
  ordering across contigs, and a case where `join_intronic_preds` actually nests
  an intron-encapsulated pred — this dataset's combine path had no nesting). To
  reproduce multi-partition golden: same stages as below but `--segmentSize 30000
  --overlapSize 10000`; run `evidence_modeler.pl` in EACH `Contig1_*-*` partition
  dir, then recombine + convert.
- **Eliminated models** (`EVM_elm` source, `*** ELIMINATED ***`) and **5'/3'
  partial** GFF3 tags (`5_prime_partial=true`) — no fixture exercises them
  (this dataset has 0). The GFF3 converter does NOT yet emit the partial tags
  (would need start/stop-codon detection on the gene model).
- **Alternate genetic code / `--forwardStrandOnly` / `--reverseStrandOnly`** end-to-end.
- The orchestrator `concatenate_gff3_outputs` uses a Rust loop; Perl uses
  `find … -regex .*evm.out.gff3 -exec cat`. Equivalent for one contig; verify
  ordering for many contigs.

---

## 00. SESSION HANDOFF (2026-06-27, latest)

**Branch `rust-rewrite-completion`. Build clean (0 warnings), `cargo test` 24/24.**

### Current parity on testing/Contig1 (single partition) — FULL PARITY
Rust `evm.out` vs Perl golden (`testing/smalltest_perl.partitions/Contig1/evm.out`,
also fixture `evm/tests/fixtures/Contig1.perl.evm.out`):
- **170 lines == 170 lines; evidence-order-normalized diff is EMPTY (byte-identical).**
- **11/11 gene spans, 79/79 exons** identical (coord+type+frame).
- **All exon/intron rows + evidence SETS identical** (0 non-header diffs with
  evidence-order normalization).
- **All 208 exon base_scores byte-exact** vs Perl `exon_list.out`.
- **ig_base byte-exact** vs Perl `intergenic.bps` (0/63304 positions differ).
- Header `score(...)`, `raw_noncoding(...)`, `offset(...)`, `noncoding_equivalent(...)`,
  `S-ratio` all byte-exact on all 11 genes.
- **Task #16 (offset residual) SOLVED.** Root cause: Rust `populate_intron_vectors`
  distributed the predicted-intron score over the RAW key span (D..A) instead of
  Perl's `intron_key_to_intron_span` span (D..A-1 for '+', A..D+1 for '-'; Perl
  `populate_forward_reverse_pred_intron_vectors` line 2974). Because the build-time
  `intron_score` covers D..A+1 but the filter offset loop iterates D..A-1
  subtracting `existing_per_base`, the wrong distribution made offset ≈ ½ of Perl's
  (residual `2·weight/intron` collapsed to `~weight/intron`). Fix: added
  `intron_key_to_intron_span()` in `introns.rs` and used it when populating the
  per-base vectors. `raw_noncoding` is unchanged (total over the full prediction
  span is conserved). Regression test in `introns.rs::tests`.

### What was completed THIS session (commits, newest first)
- `e5cade3` intergenic grouping by **full attribute column** (Perl
  `get_gene_predictions` groups CDS by `$x[8]`, not Parent; glimmerHMM/fgenesh
  CDS with a `;5_prime_partial=true` suffix split into 2 genes → intron scored
  as intergenic). Fixed ig_base tail (was 150 positions off at 62677–63076).
- `38e9792` **Phase B2 output format** — faithful `EVM_prediction::toString`
  (`!!` line, `# EVM prediction: …` header, sorted interleaved exon/INTRON rows).
- `563f115` **decrement_coding_using_protein_alignment_introns** (Perl 4206) +
  fixed `parse_evidence_chains` gap-building to sort per-link coords (median gap).
- `ed52b06`/`3a2bcd0` **convert_5prime_partials** (the "842" fix).

### KEY FINDING: Perl evm.out is NON-DETERMINISTIC
Re-running `evidence_modeler.pl` produces DIFFERENT within-row evidence ordering
each run (Perl hash randomization) — header NUMBERS are stable, only evidence
ORDER changes. So **literal byte-for-byte parity with Perl is impossible.** Rust
now **sorts evidence tokens canonically** (in `consensus.rs::format_prediction`)
so its output is reproducible. The golden-comparison oracle must normalize
evidence order on both sides (see the `norm()` awk one-liner used in this session).

### Offset residual (task #16) — SOLVED (see parity section above)
The fix landed in `algo/introns.rs::populate_intron_vectors` (use
`intron_key_to_intron_span`, not the raw key span). The diagnostic notes below are
retained for history.

`filter.rs` `offset` = Σ over the prediction's own introns of
`calc_intergenic_score(intron_span)` + Σ_i(`pred_intron_vec[i]` − `existing_per_base`).
intergenic part is now exact; `existing_per_base` = sum of ab-initio weights for
the intron (verified == Perl `get_predicted_intron_score_contribution/intron_len`).
So the residual is in the `pred_intron_vec` (FORWARD/REVERSE_PRED_INTRON_VEC, built
from `PREDICTED_INTRONS` = ab-initio introns only) summed over each intron span.
**Next step:** dump Rust `all_fwd_intron_vec`/`all_rev_intron_vec` and compare to
a Perl dump of `@FORWARD_PRED_INTRON_VEC`/`@REVERSE_PRED_INTRON_VEC` over a clean
gene (e.g. gene7 30091-32906 '+', offset 26.72 vs 54.00, raw_noncoding matches).
Perl does NOT dump the PRED vectors by default (only `introns_decomposed_to_vec.*`
from INTRONS_TO_SCORE, which is ALL evidence — different). Add a temp Perl dump
of the PRED vectors, or instrument both. Suspect: `predicted_introns` score
accumulation span vs `intron_key_to_intron_span` normalization span (off-by-bases),
or overlapping-intron contributions within a span.

### HOW TO REPRODUCE / TEST (essential commands)
Perl golden + debug dumps (run inside `testing/smalltest_perl.partitions/Contig1/`):
```
perl <repo>/EvmUtils/evidence_modeler.pl -G genome.fasta -g gene_predictions.gff3 \
  -w <repo>/testing/weights.txt -e transcript_alignments.gff3 -p protein_alignments.gff3 \
  --min_intron_length 20 --terminal_intergenic_re_search 10000 --exec_dir . --debug
```
Writes: `exon_list.out` (per-exon base_score), `intergenic.bps` (base intergenic),
`final_path` (the chosen trellis path!), `coding_vector.{+,-}.dat`,
`augment_intergenic_from_{start,stop}_peaks.dat`, `start_peaks`/`stop_peaks`.
Rust single-partition (note clap uses `--min-intron-length`, hyphens):
```
<repo>/evm/target/debug/evidence_modeler -G genome.fasta -g gene_predictions.gff3 \
  -w <repo>/testing/weights.txt -e transcript_alignments.gff3 -p protein_alignments.gff3 \
  --min-intron-length 20 --terminal-intergenic-re-search 10000 -o /tmp/rust.evm.out
```
Evidence-order-normalized diff (the real oracle, since Perl is non-deterministic):
```
norm() { awk -F'\t' '{ if ($NF ~ /^\{/) { n=split($NF,a,"},{"); for(i=1;i<=n;i++)gsub(/^\{|\}$/,"",a[i]); asort(a); s=""; for(i=1;i<=n;i++)s=s"{"a[i]"}"(i<n?",":""); $NF=s } print }' OFS='\t' "$1"; }
diff <(norm /tmp/rust.evm.out) <(norm .../evm.out)   # only header offset lines should differ
```

### AFTER offset (task #16): remaining Phase B + Phase C
- Add a golden integration test (`evm/tests/golden.rs`) using the normalized
  comparison. Decide tolerance for the offset field (or fix it first).
- Phase B4: verify `--forwardStrandOnly`/`--reverseStrandOnly` + alt stop codons.
- Phase C (next big milestone): `evm-cli` orchestrator end-to-end — partition
  parity, recombine parity, GFF3/pep/cds/bed conversion. NOTE the evm-cli
  `main.rs run_evm_on_partition` DUPLICATES the evm-utils shim merge logic — keep
  both in sync (the consensus/format fixes live in shared `evm-core` so they're
  already shared; only the merge/CLI wiring is duplicated).
- Need NEW fixtures for: multi-contig/multi-partition (recombine DP), a gene
  Perl actually ELIMINATES (filter), alternate genetic code.

### Debug instrumentation note
All temporary dumps (EVM_DUMP_EXONS / EVM_DUMP_FILTER / EVM_DUMP_IGBASE) have been
REMOVED from the committed tree. Re-add as needed for task #16.

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

**SUPERSEDED by §00 handoff — see top of file for current state.** Summary:
170/170 lines; 11/11 spans; 79/79 exons; all base_scores + ig_base byte-exact;
output format ported; evidence sets identical (Perl evidence ORDER is
non-deterministic, Rust sorts canonically). Only the filter `offset` header field
remains slightly low (task #16). Earlier notes retained below for history.

**842 SOLVED — it was a missing post-trellis step, NOT a base-score tie-break.**
The earlier diagnosis below was wrong. Perl's `--debug` `final_path` dump proves
the trellis itself picks the 5'-partial `internal 1018-1127` (identical to Rust).
The `842-1127 initial` in the golden output is produced by
`convert_5prime_partials_to_complete_genes_where_possible` (Perl line 3643),
which runs AFTER `filter_predictions_low_support` and BEFORE tail recursion: for
each multi-exon 5'-partial whose gene-start exon is `internal`, it finds the
best-scoring overlapping `initial` exon (same orient, identical `end3` +
`end_frame`) in the GLOBAL pool and replaces the gene-start exon, then re-`_init`s.
Ported in `consensus.rs` (`is_5prime_partial` + `convert_5prime_partials_to_complete_genes`,
called after the filter); `prediction.rs::finalize` now also recomputes `lend/rend`
(Perl `_init`). Confirmed via Perl debug dumps (run `evidence_modeler.pl --debug`
in `testing/smalltest_perl.partitions/Contig1`): `exon_list.out` (per-exon
base_score), `intergenic.bps` (Rust ig matched Perl exactly: 842..1017 = 352 both),
`final_path` (the trellis path).

**Score residuals → one remaining missing step:** `prediction_score` runs
slightly high and `offset` slightly low across genes. Traced to the still-missing
`decrement_coding_using_protein_alignment_introns` (Perl line 4206; constants
`INTRON_MEDIAN_FACTOR=2`, `MIN_ALIGNMENT_GAP_SIZE_INFER_INTRON`): it subtracts the
protein weight from coding scores over protein-alignment gaps (inferred introns)
shorter than `2×median_gap`. Needed for byte-exact scores; the `gaps` are already
on `EvidenceChain`. NOTE: this is NOT the 842 fix (it lowers `base(842)`, shifting
toward 1018 — opposite of what's needed).

**[RESOLVED — see "842 SOLVED" above]** The earlier root-cause hypothesis (a
base-score tie-break in `load_evidence`/`score_exons`) was WRONG. The trellis
correctly picks 1018 in both Perl and Rust; the 842 comes from the post-trellis
`convert_5prime_partials` step, now ported. The intergenic and base scores were
verified to match Perl (ig 842..1017 = 352 in both). Fixes that landed along the
way: byte-exact `analyze_peaks`; faithful `augment_intergenic_from_start/stop_peaks`
(peaks carry strand); internal-exon recovery in `recover_partial_prediction`;
`convert_5prime_partials_to_complete_genes`.

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
