//! A partition run must give the same evm.out every time. Hash iteration order
//! differs per process and per map instance, and it used to decide the order of
//! the exon pool and the float summation order of the intergenic scores.
//!
//! The unit tests in `algo::load_evidence` and `algo::intergenic` fail on the
//! old code. This end-to-end check guards the whole path, including the
//! predictions loader, which has no small entry point of its own. It passes on
//! the old code for these fixtures, so it is a regression guard, not a proof.

use std::path::PathBuf;

use evm_core::pipeline::{run_single_partition, SinglePartitionParams};

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../tests/fixtures")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn params() -> SinglePartitionParams {
    SinglePartitionParams {
        stop_codons: "TAA,TGA,TAG".to_string(),
        min_intron_length: 20,
        forward_only: false,
        reverse_only: false,
        report_elm: true,
        terminal_intergenic_re_search: 10000,
        intergenic_adjust: 1.0,
        max_prev_exons_compare: 500,
        repeats: None,
    }
}

fn run(prefix: &str) -> String {
    run_single_partition(
        &fixture(&format!("{}.genome.fasta", prefix)),
        &fixture(&format!("{}.weights.txt", prefix)),
        &fixture(&format!("{}.gene_predictions.gff3", prefix)),
        None,
        None,
        &params(),
    )
    .expect("run_single_partition")
    .concat()
}

#[test]
fn repeated_runs_give_identical_output() {
    for prefix in ["elm", "stopcodon"] {
        let first = run(prefix);
        assert!(!first.is_empty(), "{} produced no output", prefix);
        for i in 1..25 {
            assert_eq!(run(prefix), first, "{} output changed on run {}", prefix, i);
        }
    }
}
