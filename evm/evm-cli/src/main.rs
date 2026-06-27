//! EVidenceModeler — Rust re-implementation of the main EVM orchestrator.
//!
//! Drop-in replacement for the `EVidenceModeler` Perl script.

use std::fs;
use std::io::Write;
use std::path::Path;
use anyhow::{Context, Result};
use clap::Parser;
use log::{info, warn};
use rayon::prelude::*;

use evm_core::io::partitions::read_partitions_file;
use evm_core::io::weights::read_weights_file;
use evm_core::partition::partition::{run_partition, InputFile};
use evm_core::recombine::recombine::recombine_outputs;
use evm_core::gff3_convert::evm_to_gff3::convert_all_to_gff3;
use evm_core::gff3_convert::gff3_to_proteins::{extract_sequences, SeqType};
use evm_core::gff3_convert::gff3_to_bed::gff3_to_bed;

const VERSION: &str = "EVidenceModeler-v2.1.0-rust";

#[derive(Parser, Debug)]
#[command(name = "EVidenceModeler", version = VERSION, about = "Evidence Modeler — Rust implementation")]
struct Cli {
    /// Sample ID (used for naming outputs)
    #[arg(long)]
    sample_id: String,

    /// Genome sequence in FASTA format
    #[arg(long, short = 'G')]
    genome: String,

    /// Weights file for evidence types
    #[arg(long, short = 'w')]
    weights: String,

    /// Gene predictions GFF3 file
    #[arg(long, short = 'g')]
    gene_predictions: String,

    /// Segment size for genome partitioning
    #[arg(long)]
    segment_size: Option<u32>,

    #[arg(long, name = "segmentSize")]
    segment_size2: Option<u32>,

    /// Overlap size between partitions
    #[arg(long)]
    overlap_size: Option<u32>,

    #[arg(long, name = "overlapSize")]
    overlap_size2: Option<u32>,

    /// Protein alignments GFF3 file
    #[arg(long, short = 'p')]
    protein_alignments: Option<String>,

    /// Transcript alignments GFF3 file
    #[arg(long, short = 'e')]
    transcript_alignments: Option<String>,

    /// Repeats GFF3 file
    #[arg(long, short = 'r')]
    repeats: Option<String>,

    /// Terminal exons file from PASA
    #[arg(long, short = 't')]
    terminal_exons: Option<String>,

    /// Stop codons (default: TAA,TGA,TAG)
    #[arg(long, default_value = "TAA,TGA,TAG")]
    stop_codons: String,

    /// Minimum intron length (default: 20)
    #[arg(long, default_value_t = 20)]
    min_intron_length: u32,

    /// Number of CPUs for parallel execution
    #[arg(long, name = "CPU", default_value_t = 4)]
    cpu: usize,

    /// Run only on the forward strand
    #[arg(long)]
    forward_strand_only: bool,

    #[arg(long, name = "forwardStrandOnly")]
    forward_strand_only2: bool,

    /// Run only on the reverse strand
    #[arg(long)]
    reverse_strand_only: bool,

    #[arg(long, name = "reverseStrandOnly")]
    reverse_strand_only2: bool,

    /// Report eliminated EVM predictions
    #[arg(long, name = "report_ELM")]
    report_elm: bool,

    /// Verbose output
    #[arg(short = 'S')]
    verbose: bool,

    /// Debug mode
    #[arg(long)]
    debug: bool,

    /// Search for nested genes in long introns (0 = off)
    #[arg(long, default_value_t = 0)]
    search_long_introns: u32,

    /// Re-examine intergenic regions of minimum length (0 = off)
    #[arg(long, default_value_t = 0)]
    re_search_intergenic: u32,

    /// Re-examine terminal intergenic regions of minimum length
    #[arg(long, default_value_t = 10000)]
    terminal_intergenic_re_search: u32,

    /// Intergenic score adjustment factor (default: 1.0)
    #[arg(long, name = "INTERGENIC_SCORE_ADJUST_FACTOR", default_value_t = 1.0)]
    intergenic_adjust: f64,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    let segment_size = cli.segment_size.or(cli.segment_size2)
        .context("--segmentSize is required")?;
    let overlap_size = cli.overlap_size.or(cli.overlap_size2)
        .context("--overlapSize is required")?;

    let forward_only = cli.forward_strand_only || cli.forward_strand_only2;
    let reverse_only = cli.reverse_strand_only || cli.reverse_strand_only2;
    if forward_only && reverse_only {
        anyhow::bail!("--forwardStrandOnly and --reverseStrandOnly are mutually exclusive");
    }

    // Configure rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.cpu)
        .build_global()
        .ok();

    let sample_id = &cli.sample_id;
    let checkpts_dir = format!("__{}-EVM_chckpts", sample_id);
    fs::create_dir_all(&checkpts_dir).ok();

    let genome_path = &cli.genome;
    let genome_basename = Path::new(genome_path)
        .file_name().and_then(|n| n.to_str()).unwrap_or("genome.fasta").to_string();

    let partition_dir = format!("{}.partitions", sample_id);
    let partition_listing = format!("{}.partitions.listing", sample_id);

    // ─── Step 1: Partition inputs ────────────────────────────────────────────
    let partition_ckpt = format!("{}/partition_inputs.ok", checkpts_dir);
    if !Path::new(&partition_ckpt).exists() {
        info!("Partitioning inputs...");
        let mut input_files = vec![
            InputFile::new("gene_predictions", &cli.gene_predictions),
        ];
        if let Some(p) = &cli.protein_alignments {
            input_files.push(InputFile::new("protein_alignments", p));
        }
        if let Some(t) = &cli.transcript_alignments {
            input_files.push(InputFile::new("transcript_alignments", t));
        }
        if let Some(r) = &cli.repeats {
            input_files.push(InputFile::new("repeats", r));
        }

        run_partition(
            &partition_dir, genome_path, &input_files, &genome_basename,
            segment_size, overlap_size, &partition_listing,
        )?;
        fs::write(&partition_ckpt, "")?;
        info!("Partitioning complete.");
    } else {
        info!("Partitioning already done — skipping.");
    }

    let entries = read_partitions_file(&partition_listing)?;

    // ─── Step 2: Run EVM on each partition in parallel ───────────────────────
    let evm_ckpt = format!("{}/run_evm_cmds.ok", checkpts_dir);
    if !Path::new(&evm_ckpt).exists() {
        info!("Running EVM on {} partitions...", entries.len());
        let _weights = read_weights_file(&cli.weights)?;
        let _stop_codons_parsed = evm_core::algo::splice_sites::parse_stop_codons(&cli.stop_codons)?;

        // Build list of per-partition work items
        let work_items: Vec<_> = entries.iter()
            .map(|e| {
                let data_dir = if e.is_partitioned {
                    e.partition_dir.clone().unwrap_or(e.base_dir.clone())
                } else {
                    e.base_dir.clone()
                };
                (e.accession.clone(), data_dir)
            })
            .collect();

        let results: Vec<Result<()>> = work_items.par_iter().map(|(acc, data_dir)| {
            run_evm_on_partition(
                acc, data_dir, &genome_basename, &cli.gene_predictions,
                cli.protein_alignments.as_deref(),
                cli.transcript_alignments.as_deref(),
                &cli.weights, &cli.stop_codons,
                cli.min_intron_length, forward_only, reverse_only,
                cli.report_elm, cli.search_long_introns,
                cli.re_search_intergenic, cli.terminal_intergenic_re_search,
                cli.intergenic_adjust,
            )
        }).collect();

        let mut had_error = false;
        for r in results {
            if let Err(e) = r {
                warn!("EVM partition error: {:#}", e);
                had_error = true;
            }
        }
        if had_error {
            anyhow::bail!("One or more EVM partitions failed");
        }
        fs::write(&evm_ckpt, "")?;
        info!("All EVM partitions complete.");
    }

    // ─── Step 3: Recombine partial outputs ───────────────────────────────────
    let recombine_ckpt = format!("{}/recombined_EVM_partial_outputs.ok", checkpts_dir);
    if !Path::new(&recombine_ckpt).exists() {
        info!("Recombining partial outputs...");
        recombine_outputs(&entries, "evm.out")?;
        fs::write(&recombine_ckpt, "")?;
    }

    // ─── Step 4: Convert to GFF3 ─────────────────────────────────────────────
    let gff3_ckpt = format!("{}/evm_out_to_gff3_format.ok", checkpts_dir);
    if !Path::new(&gff3_ckpt).exists() {
        info!("Converting EVM output to GFF3...");
        convert_all_to_gff3(&entries, "evm.out")?;
        fs::write(&gff3_ckpt, "")?;
    }

    // ─── Step 5: Concatenate GFF3 outputs ────────────────────────────────────
    let concat_ckpt = format!("{}/concatenate_final_output_gff3.ok", checkpts_dir);
    let final_gff3 = format!("{}.EVM.gff3", sample_id);
    if !Path::new(&concat_ckpt).exists() {
        info!("Concatenating GFF3 outputs...");
        concatenate_gff3_outputs(&entries, &partition_dir, &final_gff3)?;
        fs::write(&concat_ckpt, "")?;
    }

    // ─── Step 6: Extract protein sequences ───────────────────────────────────
    let pep_ckpt = format!("{}/make_evm_pep.ok", checkpts_dir);
    let pep_file = format!("{}.EVM.pep", sample_id);
    if !Path::new(&pep_ckpt).exists() && Path::new(&final_gff3).exists() {
        info!("Extracting protein sequences...");
        let default_stops: Vec<[u8; 3]> = vec![*b"TAA", *b"TGA", *b"TAG"];
        let seqs = extract_sequences(&final_gff3, genome_path, &SeqType::Prot, &default_stops)?;
        let mut out = fs::File::create(&pep_file)?;
        for (hdr, seq) in &seqs {
            writeln!(out, ">{}", hdr)?;
            for chunk in seq.as_bytes().chunks(60) {
                writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
            }
        }
        fs::write(&pep_ckpt, "")?;
    }

    // ─── Step 7: Extract CDS sequences ───────────────────────────────────────
    let cds_ckpt = format!("{}/make_evm_cds.ok", checkpts_dir);
    let cds_file = format!("{}.EVM.cds", sample_id);
    if !Path::new(&cds_ckpt).exists() && Path::new(&final_gff3).exists() {
        info!("Extracting CDS sequences...");
        let default_stops: Vec<[u8; 3]> = vec![*b"TAA", *b"TGA", *b"TAG"];
        let seqs = extract_sequences(&final_gff3, genome_path, &SeqType::Cds, &default_stops)?;
        let mut out = fs::File::create(&cds_file)?;
        for (hdr, seq) in &seqs {
            writeln!(out, ">{}", hdr)?;
            for chunk in seq.as_bytes().chunks(60) {
                writeln!(out, "{}", std::str::from_utf8(chunk).unwrap())?;
            }
        }
        fs::write(&cds_ckpt, "")?;
    }

    // ─── Step 8: Make BED file ────────────────────────────────────────────────
    let bed_ckpt = format!("{}/make_bed.ok", checkpts_dir);
    let bed_file = format!("{}.EVM.bed", sample_id);
    if !Path::new(&bed_ckpt).exists() && Path::new(&final_gff3).exists() {
        info!("Generating BED file...");
        let bed_lines = gff3_to_bed(&final_gff3)?;
        let mut out = fs::File::create(&bed_file)?;
        for line in &bed_lines {
            writeln!(out, "{}", line)?;
        }
        fs::write(&bed_ckpt, "")?;
    }

    info!("Done. See {}.EVM.* outputs", sample_id);
    Ok(())
}

/// Run the EVM algorithm on a single partition directory.
#[allow(clippy::too_many_arguments)]
fn run_evm_on_partition(
    _accession: &str,
    data_dir: &str,
    genome_basename: &str,
    gene_pred_global: &str,
    protein_global: Option<&str>,
    transcript_global: Option<&str>,
    weights_path: &str,
    stop_codons_str: &str,
    min_intron_length: u32,
    forward_only: bool,
    reverse_only: bool,
    report_elm: bool,
    _search_long_introns: u32,
    _re_search_intergenic: u32,
    terminal_intergenic_re_search: u32,
    intergenic_adjust: f64,
) -> Result<()> {
    let output_path = format!("{}/evm.out", data_dir);
    let _log_path = format!("{}/evm.out.log", data_dir);

    // Skip if already done
    let ckpt = format!("{}/evm.done.ok", data_dir);
    if Path::new(&ckpt).exists() {
        return Ok(());
    }

    use evm_core::io::fasta::read_fasta_file;
    use evm_core::io::gff3::read_gff3_file;
    use evm_core::io::weights::read_weights_file;
    use evm_core::types::genome::{GenomeSequence, MaskVec};
    use evm_core::algo::splice_sites::parse_stop_codons;
    use evm_core::algo::intergenic::{populate_intergenic_scores, augment_intergenic_from_start_stop_peaks};
    use evm_core::algo::process::{process_features, ProcessConfig};
    use evm_core::algo::consensus::{generate_consensus_gene_predictions, ConsensusParams};
    
    use evm_core::types::prediction::PredMode;

    let genome_path = format!("{}/{}", data_dir, genome_basename);
    let gene_pred_path = format!("{}/{}", data_dir,
        Path::new(gene_pred_global).file_name().and_then(|n| n.to_str()).unwrap_or(""));
    let protein_path = protein_global.map(|p| format!("{}/{}", data_dir,
        Path::new(p).file_name().and_then(|n| n.to_str()).unwrap_or("")));
    let transcript_path = transcript_global.map(|p| format!("{}/{}", data_dir,
        Path::new(p).file_name().and_then(|n| n.to_str()).unwrap_or("")));

    let fasta_records = read_fasta_file(&genome_path)?;
    if fasta_records.is_empty() { return Ok(()); }
    let genome_seq = GenomeSequence::new(&fasta_records[0].sequence);
    let seq_len = genome_seq.len();

    let ev_weights = read_weights_file(weights_path)?;
    let stop_codons = parse_stop_codons(stop_codons_str)?;

    let gene_pred_records = read_gff3_file(&gene_pred_path).unwrap_or_default();
    let protein_records = protein_path.as_deref().map(|p| read_gff3_file(p).unwrap_or_default());
    let transcript_records = transcript_path.as_deref().map(|p| read_gff3_file(p).unwrap_or_default());

    let mask = MaskVec::new(seq_len);

    // Sum prediction weights
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
        min_intron_length,
        mask: &mask,
        sum_genepred_weights: sum_pred_weights,
        chain_termini_window: 250,
    };

    let fwd_state = if !reverse_only {
        Some(process_features('+', &cfg)?)
    } else { None };

    let rev_genome = genome_seq.to_reverse_complement();
    let rev_cfg = ProcessConfig {
        genome_seq: &rev_genome,
        ..cfg
    };
    let rev_state = if !forward_only {
        let mut s = process_features('-', &rev_cfg)?;
        // Transpose exons back to forward strand: revcomp coords, remap
        // reading frames 1->4/2->5/3->6, flip orientation (mirrors Perl
        // transpose_exons_back_to_forward_strand).
        use evm_core::types::exon::Orientation;
        for exon in s.exons.iter_mut() {
            let new_e5 = seq_len as u32 - exon.end5 + 1;
            let new_e3 = seq_len as u32 - exon.end3 + 1;
            exon.end5 = new_e5;
            exon.end3 = new_e3;
            exon.start_frame = match exon.start_frame { 1 => 4, 2 => 5, 3 => 6, o => o };
            exon.end_frame = match exon.end_frame { 1 => 4, 2 => 5, 3 => 6, o => o };
            exon.orientation = Orientation::Rev;
        }
        Some(s)
    } else { None };

    // Merge states
    let mut all_exons: Vec<evm_core::types::exon::Exon> = Vec::new();
    let mut all_introns_to_score = std::collections::HashMap::new();
    let mut all_introns_to_evidence: std::collections::HashMap<String, Vec<(String, String)>> = std::collections::HashMap::new();
    let mut all_predicted_introns = std::collections::HashMap::new();
    let mut all_coding_scores = evm_core::algo::coding_scores::new_coding_scores(seq_len);
    let mut all_start_peaks = Vec::new();
    let mut all_end_peaks = Vec::new();
    let mut all_fwd_intron_vec = vec![0.0f64; seq_len + 2];
    let mut all_rev_intron_vec = vec![0.0f64; seq_len + 2];

    for state in [fwd_state, rev_state].into_iter().flatten() {
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

    // Base intergenic (filter) + augmented copy (trellis).
    let ig_base = populate_intergenic_scores(seq_len, &gene_pred_records, &ev_weights, &mask, intergenic_adjust);
    let mut ig_scores = ig_base.clone();
    augment_intergenic_from_start_stop_peaks(
        &mut ig_scores, &all_start_peaks, &all_end_peaks, &all_exons, &mask,
        seq_len as u32, sum_pred_weights, 500,
    );

    let (acceptable, phased, intergenic_conns, frame_pairs) =
        evm_core::types::exon::build_acceptable_linkages();

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
        max_prev_exons_compare: 500,
        min_intergenic_size_on_re_search: terminal_intergenic_re_search,
        min_gene_length_size_on_re_search: 0,
        report_elm,
        recursion_limit: 200,
    };

    generate_consensus_gene_predictions(
        1, seq_len as u32, PredMode::Standard,
        &mut params, &mut recursion_count, &mut output,
    )?;

    // Write output
    let mut out_file = fs::File::create(&output_path)?;
    for block in &output {
        write!(out_file, "{}", block)?;
    }

    fs::write(&ckpt, "")?;
    Ok(())
}

/// Concatenate all per-contig GFF3 outputs into the final output file.
fn concatenate_gff3_outputs(
    entries: &[evm_core::io::partitions::PartitionEntry],
    _partition_dir: &str,
    output_path: &str,
) -> Result<()> {
    use std::collections::HashSet;
    let mut seen_dirs: HashSet<String> = HashSet::new();
    let mut out = fs::File::create(output_path)?;

    for entry in entries {
        if seen_dirs.insert(entry.base_dir.clone()) {
            let gff3 = format!("{}/evm.out.gff3", entry.base_dir);
            if Path::new(&gff3).exists() {
                let content = fs::read_to_string(&gff3)?;
                write!(out, "{}", content)?;
            }
        }
    }
    Ok(())
}
