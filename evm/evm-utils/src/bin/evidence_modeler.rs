//! evidence_modeler — single-partition EVM invocation.
//!
//! Replaces EvmUtils/evidence_modeler.pl for backward-compatible usage. The
//! actual algorithm lives in `evm_core::pipeline::run_single_partition`, shared
//! with the `EVidenceModeler` orchestrator.

use std::fs;
use std::io::Write;
use anyhow::Result;
use clap::Parser;

use evm_core::pipeline::{run_single_partition, SinglePartitionParams};

#[derive(Parser, Debug)]
#[command(name = "evidence_modeler", about = "Run EVM on a single partition")]
struct Cli {
    /// Genome FASTA (or partition FASTA)
    #[arg(long, short = 'G')]
    genome: String,

    /// Weights file
    #[arg(long, short = 'w')]
    weights: String,

    /// Gene predictions GFF3
    #[arg(long, short = 'g')]
    gene_predictions: String,

    /// Protein alignments GFF3
    #[arg(long, short = 'p')]
    protein_alignments: Option<String>,

    /// Transcript alignments GFF3
    #[arg(long, short = 'e')]
    transcript_alignments: Option<String>,

    /// Output file (default: stdout)
    #[arg(long, short = 'o')]
    output: Option<String>,

    /// Stop codons, comma-separated (default: TAA,TGA,TAG)
    #[arg(long, default_value = "TAA,TGA,TAG")]
    stop_codons: String,

    /// Minimum intron length (default: 20)
    #[arg(long, default_value_t = 20)]
    min_intron_length: u32,

    /// Forward strand only
    #[arg(long = "forwardStrandOnly")]
    forward_strand_only: bool,

    /// Reverse strand only
    #[arg(long = "reverseStrandOnly")]
    reverse_strand_only: bool,

    /// Report eliminated models
    #[arg(long = "report_ELM")]
    report_elm: bool,

    /// Intergenic score adjustment factor
    #[arg(long = "INTERGENIC_SCORE_ADJUST_FACTOR", default_value_t = 1.0)]
    intergenic_adjust: f64,

    /// Minimum intergenic size for terminal region re-search
    #[arg(long, default_value_t = 10000)]
    terminal_intergenic_re_search: u32,

    #[arg(long, default_value_t = 500)]
    max_prev_exons_compare: usize,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();
    let cli = Cli::parse();

    let params = SinglePartitionParams {
        stop_codons: cli.stop_codons,
        min_intron_length: cli.min_intron_length,
        forward_only: cli.forward_strand_only,
        reverse_only: cli.reverse_strand_only,
        report_elm: cli.report_elm,
        terminal_intergenic_re_search: cli.terminal_intergenic_re_search,
        intergenic_adjust: cli.intergenic_adjust,
        max_prev_exons_compare: cli.max_prev_exons_compare,
    };

    let output = run_single_partition(
        &cli.genome,
        &cli.weights,
        &cli.gene_predictions,
        cli.protein_alignments.as_deref(),
        cli.transcript_alignments.as_deref(),
        &params,
    )?;

    match &cli.output {
        Some(path) => {
            let mut f = fs::File::create(path)?;
            for block in &output { write!(f, "{}", block)?; }
        }
        None => {
            for block in &output { print!("{}", block); }
        }
    }

    Ok(())
}
