//! evidence_modeler — single-partition EVM invocation.
//!
//! Replaces EvmUtils/evidence_modeler.pl for backward-compatible usage.

use std::fs;
use std::io::Write;
use anyhow::Result;
use clap::Parser;

use evm_core::io::fasta::read_fasta_file;
use evm_core::io::gff3::read_gff3_file;
use evm_core::io::weights::read_weights_file;
use evm_core::types::genome::{GenomeSequence, MaskVec};
use evm_core::algo::splice_sites::parse_stop_codons;
use evm_core::algo::intergenic::{populate_intergenic_scores, augment_intergenic_from_start_stop_peaks};
use evm_core::algo::process::{process_features, ProcessConfig};
use evm_core::algo::consensus::{generate_consensus_gene_predictions, ConsensusParams};
use evm_core::types::exon::build_acceptable_linkages;
use evm_core::types::prediction::PredMode;

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
    #[arg(long, name = "forwardStrandOnly")]
    forward_strand_only: bool,

    /// Reverse strand only
    #[arg(long, name = "reverseStrandOnly")]
    reverse_strand_only: bool,

    /// Report eliminated models
    #[arg(long, name = "report_ELM")]
    report_elm: bool,

    /// Intergenic score adjustment factor
    #[arg(long, name = "INTERGENIC_SCORE_ADJUST_FACTOR", default_value_t = 1.0)]
    intergenic_adjust: f64,

    /// Minimum intergenic size for terminal region re-search
    #[arg(long, default_value_t = 10000)]
    terminal_intergenic_re_search: u32,

    #[arg(long, default_value_t = 500)]
    max_prev_exons_compare: usize,
}

/// Map a forward-strand reading frame (1,2,3) to its reverse-strand
/// equivalent (4,5,6) when transposing reverse exons back to forward coords.
fn fwd_frame_to_rev(frame: evm_core::types::exon::ExonPhase) -> evm_core::types::exon::ExonPhase {
    match frame {
        1 => 4,
        2 => 5,
        3 => 6,
        other => other,
    }
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();
    let cli = Cli::parse();

    let fasta_records = read_fasta_file(&cli.genome)?;
    if fasta_records.is_empty() {
        anyhow::bail!("No sequences found in genome file: {}", cli.genome);
    }
    let genome_seq = GenomeSequence::new(&fasta_records[0].sequence);
    let seq_len = genome_seq.len();

    let ev_weights = read_weights_file(&cli.weights)?;
    let stop_codons = parse_stop_codons(&cli.stop_codons)?;

    let gene_pred_records = read_gff3_file(&cli.gene_predictions)?;
    let protein_records = cli.protein_alignments.as_deref()
        .map(|p| read_gff3_file(p).unwrap_or_default());
    let transcript_records = cli.transcript_alignments.as_deref()
        .map(|p| read_gff3_file(p).unwrap_or_default());

    let mask = MaskVec::new(seq_len);

    let sum_pred_weights: f64 = ev_weights.values()
        .filter(|e| e.ev_class.is_prediction())
        .map(|e| e.weight)
        .sum();

    let cfg = ProcessConfig {
        genome_seq: &genome_seq,
        stop_codons: &stop_codons,
        ev_weights: &ev_weights,
        gene_pred_records: &gene_pred_records,
        protein_records: protein_records.as_deref(),
        transcript_records: transcript_records.as_deref(),
        min_intron_length: cli.min_intron_length,
        mask: &mask,
        sum_genepred_weights: sum_pred_weights,
        chain_termini_window: 250,
    };

    let fwd_state = if !cli.reverse_strand_only {
        Some(process_features('+', &cfg)?)
    } else { None };

    let rev_genome = genome_seq.to_reverse_complement();
    let rev_cfg = ProcessConfig {
        genome_seq: &rev_genome,
        ..cfg
    };
    let rev_state_raw = if !cli.forward_strand_only {
        Some(process_features('-', &rev_cfg)?)
    } else { None };

    // Merge exon pools (reversing coordinates for rev strand)
    let mut all_exons: Vec<evm_core::types::exon::Exon> = Vec::new();
    let mut all_introns_to_score = std::collections::HashMap::new();
    let mut all_introns_to_evidence: std::collections::HashMap<String, Vec<(String, String)>> = std::collections::HashMap::new();
    let mut all_predicted_introns = std::collections::HashMap::new();
    let mut all_coding_scores = evm_core::algo::coding_scores::new_coding_scores(seq_len);
    let mut all_start_peaks = Vec::new();
    let mut all_end_peaks = Vec::new();
    let mut all_fwd_intron_vec = vec![0.0f64; seq_len + 2];
    let mut all_rev_intron_vec = vec![0.0f64; seq_len + 2];

    for state in [fwd_state].into_iter().flatten() {
        all_exons.extend(state.exons);
        for (k, v) in state.introns_to_score { *all_introns_to_score.entry(k).or_insert(0.0) += v; }
        for (k, v) in state.introns_to_evidence { all_introns_to_evidence.entry(k).or_default().extend(v); }
        for (k, v) in state.predicted_introns { *all_predicted_introns.entry(k).or_insert(0.0) += v; }
        for (i, &v) in state.coding_scores.iter().enumerate() { all_coding_scores[i] += v; }
        all_start_peaks.extend(state.start_peaks);
        all_end_peaks.extend(state.end_peaks);
        for (i, &v) in state.fwd_intron_vec.iter().enumerate() { if i < all_fwd_intron_vec.len() { all_fwd_intron_vec[i] += v; } }
        for (i, &v) in state.rev_intron_vec.iter().enumerate() { if i < all_rev_intron_vec.len() { all_rev_intron_vec[i] += v; } }
    }

    if let Some(mut state) = rev_state_raw {
        // Transpose reverse-strand exons back to the forward coordinate system.
        // Mirrors Perl transpose_exons_back_to_forward_strand: revcomp the
        // coordinates, remap reading frames 1→4/2→5/3→6, and flip orientation.
        use evm_core::types::exon::Orientation;
        for exon in state.exons.iter_mut() {
            let new_e5 = seq_len as u32 - exon.end5 + 1;
            let new_e3 = seq_len as u32 - exon.end3 + 1;
            exon.end5 = new_e5;
            exon.end3 = new_e3;
            exon.start_frame = fwd_frame_to_rev(exon.start_frame);
            exon.end_frame = fwd_frame_to_rev(exon.end_frame);
            exon.orientation = Orientation::Rev;
        }
        all_exons.extend(state.exons);
        for (k, v) in state.introns_to_score { *all_introns_to_score.entry(k).or_insert(0.0) += v; }
        for (k, v) in state.introns_to_evidence { all_introns_to_evidence.entry(k).or_default().extend(v); }
        for (k, v) in state.predicted_introns { *all_predicted_introns.entry(k).or_insert(0.0) += v; }
        for (i, &v) in state.coding_scores.iter().enumerate() { all_coding_scores[i] += v; }
        all_start_peaks.extend(state.start_peaks);
        all_end_peaks.extend(state.end_peaks);
        for (i, &v) in state.fwd_intron_vec.iter().enumerate() { if i < all_fwd_intron_vec.len() { all_fwd_intron_vec[i] += v; } }
        for (i, &v) in state.rev_intron_vec.iter().enumerate() { if i < all_rev_intron_vec.len() { all_rev_intron_vec[i] += v; } }
    }

    // Base intergenic (non-augmented) for the low-support filter; augmented copy
    // (start/stop peak augmentation) for the trellis.
    let ig_base = populate_intergenic_scores(seq_len, &gene_pred_records, &ev_weights, &mask, cli.intergenic_adjust);
    let mut ig_scores = ig_base.clone();
    augment_intergenic_from_start_stop_peaks(
        &mut ig_scores, &all_start_peaks, &all_end_peaks, &all_exons, &mask,
        seq_len as u32, sum_pred_weights, 500,
    );

    let (acceptable, phased, intergenic_conns, frame_pairs) = build_acceptable_linkages();

    let mut output: Vec<String> = Vec::new();
    let mut recursion_count = 0usize;

    let mut params = ConsensusParams {
        exons: &mut all_exons,
        introns_to_score: &all_introns_to_score,
        introns_to_evidence: &all_introns_to_evidence,
        ev_weights: &ev_weights,
        mask: &mask,
        intergenic_scores: &ig_scores,
        base_intergenic_scores: &ig_base,
        acceptable_linkages: &acceptable,
        phased_connections: &phased,
        intergenic_connections: &intergenic_conns,
        frame_pairs: &frame_pairs,
        stop_codons: &stop_codons,
        coding_scores: &all_coding_scores,
        fwd_intron_vec: &all_fwd_intron_vec,
        rev_intron_vec: &all_rev_intron_vec,
        max_prev_exons_compare: cli.max_prev_exons_compare,
        min_intergenic_size_on_re_search: cli.terminal_intergenic_re_search,
        min_gene_length_size_on_re_search: 0,
        report_elm: cli.report_elm,
        recursion_limit: 200,
    };

    generate_consensus_gene_predictions(
        1, seq_len as u32, PredMode::Standard,
        &mut params, &mut recursion_count, &mut output,
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
