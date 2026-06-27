//! Recombine partial EVM outputs from overlapping partitions.
//!
//! Mirrors the Perl `recombine_EVM_partial_outputs.pl` script.

use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, Write};
use anyhow::{Context, Result};
use crate::io::partitions::PartitionEntry;
use crate::types::prediction::{PartitionPred, PredClass};

/// Parse an EVM output file from a single partition, offsetting coordinates
/// by `partition_lend - 1` to map them back to the full-contig space.
pub fn parse_and_add_predictions(
    output_file: &str,
    partition_lend: u32,
    predictions: &mut Vec<PartitionPred>,
) -> Result<()> {
    log::debug!("Parsing {}", output_file);
    let f = match fs::File::open(output_file) {
        Ok(f) => f,
        Err(e) => {
            log::warn!("Cannot open {}: {}", output_file, e);
            return Ok(());
        }
    };
    let reader = std::io::BufReader::new(f);

    let mut current_text = String::new();
    let mut preds: Vec<PartitionPred> = Vec::new();

    let process = |text: &str, offset: u32, out: &mut Vec<PartitionPred>| {
        if text.is_empty() { return; }
        if let Some(pred) = process_prediction_text(text, offset) {
            out.push(pred);
        }
    };

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("!!") { continue; }
        if line.starts_with('#') && !line.contains("EVM") { continue; }

        if line.starts_with(|c: char| c.is_ascii_digit() || c == '#') {
            current_text.push_str(&line);
            current_text.push('\n');
        } else {
            if !current_text.is_empty() {
                process(&current_text, partition_lend, &mut preds);
                current_text.clear();
            }
        }
    }
    if !current_text.is_empty() {
        process(&current_text, partition_lend, &mut preds);
    }

    // Join nested (intronic) predictions
    let final_preds = join_intronic_preds(preds);
    predictions.extend(final_preds);
    Ok(())
}

/// Parse a single prediction text block and adjust coordinates.
fn process_prediction_text(text: &str, partition_lend: u32) -> Option<PartitionPred> {
    let offset = partition_lend.saturating_sub(1);
    let mut lines_iter = text.lines();

    // First line is the header
    let header_line = lines_iter.next()?;
    let header_parts: Vec<&str> = header_line.split_whitespace().collect();

    // Find coord span in header: look for "span:lend-rend" or positional field 6 (0-indexed)
    let coordspan = header_parts.iter()
        .find(|s| s.contains('-'))
        .and_then(|s| {
            // might be "span:100-500" or just "100-500"
            let s = s.trim_start_matches("span:");
            let parts: Vec<&str> = s.split('-').collect();
            if parts.len() == 2 {
                let l: u32 = parts[0].parse().ok()?;
                let r: u32 = parts[1].parse().ok()?;
                Some((l + offset, r + offset))
            } else { None }
        });

    let mut exon_types: Vec<String> = Vec::new();
    let mut all_coords: Vec<u32> = Vec::new();
    let mut new_text = {
        let parts: Vec<String> = header_parts.iter().map(|&s| {
            if let Some((l, r)) = coordspan {
                if s.contains('-') && s.trim_start_matches("span:").contains('-') {
                    let stripped = s.trim_start_matches("span:");
                    if stripped.split('-').count() == 2 {
                        let prefix = if s.starts_with("span:") { "span:" } else { "" };
                        return format!("{}{}-{}", prefix, l, r);
                    }
                }
            }
            s.to_string()
        }).collect();
        format!("{}\n", parts.join(" "))
    };

    for data_line in lines_iter {
        let cols: Vec<&str> = data_line.splitn(6, '\t').collect();
        if cols.len() >= 3 {
            if let (Ok(e5), Ok(e3)) = (cols[0].parse::<u32>(), cols[1].parse::<u32>()) {
                let etype = cols[2].to_string();
                if etype == "INTRON" { continue; }
                let new_e5 = e5 + offset;
                let new_e3 = e3 + offset;
                all_coords.push(new_e5);
                all_coords.push(new_e3);
                exon_types.push(etype);
                let rest: String = cols[3..].join("\t");
                new_text.push_str(&format!("{}\t{}\t{}\n", new_e5, new_e3, rest));
            }
        }
    }

    if all_coords.is_empty() { return None; }

    let gene_lend = *all_coords.iter().min()?;
    let gene_rend = *all_coords.iter().max()?;
    let length = gene_rend - gene_lend + 1;

    let type_str = exon_types.join(",");
    let class = if (type_str.contains("initial") && type_str.contains("terminal"))
        || type_str.contains("single")
    {
        PredClass::Complete
    } else {
        PredClass::Partial
    };

    Some(PartitionPred {
        lend: gene_lend,
        rend: gene_rend,
        class,
        text: new_text,
        length,
        path_score: length,
        prev_link: None,
        intronic_preds: Vec::new(),
        encaps: false,
    })
}

/// Identify predictions encapsulated within the introns of others and nest them.
pub fn join_intronic_preds(mut preds: Vec<PartitionPred>) -> Vec<PartitionPred> {
    preds.sort_by_key(|p| p.lend);

    let n = preds.len();
    let mut encaps = vec![false; n];

    for i in 0..n {
        for j in (i + 1)..n {
            if preds[j].lend > preds[i].lend && preds[j].rend < preds[i].rend {
                encaps[j] = true;
                let nested = preds[j].clone();
                let len_j = nested.length;
                let score_j = nested.path_score;
                preds[i].intronic_preds.push(nested);
                preds[i].length += len_j;
                preds[i].path_score += score_j;
            }
        }
    }

    preds.into_iter()
        .zip(encaps)
        .filter_map(|(p, enc)| if enc { None } else { Some(p) })
        .collect()
}

/// Dynamic-programming combination of predictions from multiple overlapping partitions.
/// Selects the maximal set of non-overlapping complete genes.
pub fn combine_predictions(mut preds: Vec<PartitionPred>) -> Vec<PartitionPred> {
    if preds.is_empty() { return preds; }
    preds.sort_by_key(|p| p.lend);
    let n = preds.len();

    for i in 1..n {
        let lend_i = preds[i].lend;
        let length_i = preds[i].length;
        let mut best_score = preds[i].path_score;
        let mut best_link: Option<usize> = None;

        for j in (0..i).rev() {
            let rend_j = preds[j].rend;
            if rend_j < lend_i {
                let candidate = preds[j].path_score + length_i;
                if candidate > best_score {
                    best_score = candidate;
                    best_link = Some(j);
                }
            }
        }
        preds[i].path_score = best_score;
        preds[i].prev_link = best_link;
    }

    // Find highest-scoring end
    let best_end = preds.iter()
        .enumerate()
        .max_by_key(|(_, p)| p.path_score)
        .map(|(i, _)| i);

    let mut result: Vec<PartitionPred> = Vec::new();
    let mut cur = best_end;
    while let Some(idx) = cur {
        result.push(preds[idx].clone());
        cur = preds[idx].prev_link;
    }
    result.reverse();
    result
}

/// Run recombination for all contigs listed in `entries`.
pub fn recombine_outputs(
    entries: &[PartitionEntry],
    output_file_name: &str,
) -> Result<()> {
    // Group entries by base_dir
    let mut base_to_partitions: HashMap<String, Vec<(String, u32)>> = HashMap::new();
    for entry in entries {
        if entry.is_partitioned {
            if let Some(pdir) = &entry.partition_dir {
                // Extract lend from partition dir name "..._LEND-REND"
                let lend = extract_partition_lend(pdir).unwrap_or(1);
                base_to_partitions
                    .entry(entry.base_dir.clone())
                    .or_default()
                    .push((pdir.clone(), lend));
            }
        }
    }

    for (base_dir, partition_dirs) in &base_to_partitions {
        let mut all_preds: Vec<PartitionPred> = Vec::new();
        for (pdir, lend) in partition_dirs {
            let output_path = format!("{}/{}", pdir, output_file_name);
            parse_and_add_predictions(&output_path, *lend, &mut all_preds)?;
        }

        let final_preds = combine_predictions(all_preds);
        let out_path = format!("{}/{}", base_dir, output_file_name);
        let mut out = fs::File::create(&out_path)
            .with_context(|| format!("Cannot create {}", out_path))?;

        log::debug!("Writing combined output to {}", out_path);
        for pred in &final_preds {
            write!(out, "{}", pred.text)?;
            for nested in &pred.intronic_preds {
                writeln!(out, "!! Intron-containing prediction")?;
                write!(out, "{}", nested.text)?;
            }
        }
    }
    Ok(())
}

fn extract_partition_lend(pdir: &str) -> Option<u32> {
    // Partition dir ends with _LEND-REND
    let name = std::path::Path::new(pdir).file_name()?.to_str()?;
    let last = name.split('_').last()?;
    let lend_str = last.split('-').next()?;
    lend_str.parse().ok()
}
