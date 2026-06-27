//! Convert EVM output format to GFF3.
//!
//! Mirrors `EVM_to_GFF3.pl` and `convert_EVM_outputs_to_GFF3.pl`.

use std::fs;
use std::io::{BufRead, Write};
use anyhow::{Context, Result};
use crate::io::partitions::PartitionEntry;

/// Phase conversion map: EVM phase [1..=6] → GFF3 phase [0, 1, 2].
fn evm_phase_to_gff3(phase: u8) -> u8 {
    match phase { 1 | 4 => 0, 2 | 5 => 1, 3 | 6 => 2, _ => 0 }
}

/// Convert a single EVM output file to GFF3 format, writing to `output_path`.
pub fn evm_output_to_gff3(
    evm_file: &str,
    contig_id: &str,
    output_path: &str,
) -> Result<()> {
    let f = fs::File::open(evm_file)
        .with_context(|| format!("Cannot open EVM output {}", evm_file))?;
    let reader = std::io::BufReader::new(f);
    let mut out = fs::File::create(output_path)
        .with_context(|| format!("Cannot create GFF3 output {}", output_path))?;

    // State
    let mut model_id = 1u32;
    let mut current_coords: Vec<(u32, u32, u8)> = Vec::new(); // (end5, end3, phase)
    let mut is_eliminated = false;

    let flush = |out: &mut dyn Write,
                 coords: &mut Vec<(u32, u32, u8)>,
                 model_id: &mut u32,
                 contig: &str,
                 elim: bool| -> Result<()> {
        if coords.is_empty() { return Ok(()); }
        let ev_type = if elim { "EVM_elm" } else { "EVM" };
        let (gene_lend, gene_rend, strand) = compute_span(coords, contig);
        let tu_id = format!("evm.TU.{}.{}", contig, model_id);
        let model_feat = format!("evm.model.{}.{}", contig, model_id);

        // gene feature
        writeln!(out, "{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tID={};Name={}",
            contig, ev_type, gene_lend, gene_rend, strand, tu_id, tu_id)?;

        // mRNA feature
        writeln!(out, "{}\t{}\tmRNA\t{}\t{}\t.\t{}\t.\tID={};Parent={};Name={}",
            contig, ev_type, gene_lend, gene_rend, strand, model_feat, tu_id, model_feat)?;

        // CDS features
        let mut sorted = coords.clone();
        sorted.sort_by_key(|&(e5, e3, _)| e5.min(e3));
        for (end5, end3, phase) in &sorted {
            let (cs, ce) = (end5.min(end3), end5.max(end3));
            let gff_phase = evm_phase_to_gff3(*phase);
            writeln!(out, "{}\t{}\tCDS\t{}\t{}\t.\t{}\t{}\tID=cds.{};Parent={}",
                contig, ev_type, cs, ce, strand, gff_phase, model_feat, model_feat)?;
        }
        writeln!(out)?;
        *model_id += 1;
        coords.clear();
        Ok(())
    };

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("!!") { continue; }
        if line.starts_with('#') {
            flush(&mut out, &mut current_coords, &mut model_id, contig_id, is_eliminated)?;
            is_eliminated = line.contains("ELIMINATED");
            continue;
        }
        if line.trim().is_empty() {
            flush(&mut out, &mut current_coords, &mut model_id, contig_id, is_eliminated)?;
            continue;
        }

        let cols: Vec<&str> = line.splitn(6, '\t').collect();
        if cols.len() >= 4 {
            if let (Ok(e5), Ok(e3), Ok(phase)) = (
                cols[0].parse::<u32>(), cols[1].parse::<u32>(), cols[3].parse::<u8>()
            ) {
                let etype = cols[2];
                if etype != "INTRON" {
                    current_coords.push((e5, e3, phase));
                }
            }
        }
    }
    flush(&mut out, &mut current_coords, &mut model_id, contig_id, is_eliminated)?;
    Ok(())
}

fn compute_span(coords: &[(u32, u32, u8)], _contig: &str) -> (u32, u32, char) {
    let lend = coords.iter().map(|&(e5, e3, _)| e5.min(e3)).min().unwrap_or(0);
    let rend = coords.iter().map(|&(e5, e3, _)| e5.max(e3)).max().unwrap_or(0);
    // Determine strand from first exon (end5 < end3 → '+')
    let strand = if coords[0].0 <= coords[0].1 { '+' } else { '-' };
    (lend, rend, strand)
}

/// Convert EVM outputs for all entries to GFF3.
pub fn convert_all_to_gff3(
    entries: &[PartitionEntry],
    output_file_name: &str,
) -> Result<()> {
    use std::collections::HashMap;
    let mut base_dirs: HashMap<String, String> = HashMap::new();
    for entry in entries {
        base_dirs.insert(entry.accession.clone(), entry.base_dir.clone());
    }
    for (accession, base_dir) in &base_dirs {
        let evm_file = format!("{}/{}", base_dir, output_file_name);
        let gff3_file = format!("{}/{}.gff3", base_dir, output_file_name);
        evm_output_to_gff3(&evm_file, accession, &gff3_file)?;
    }
    Ok(())
}
