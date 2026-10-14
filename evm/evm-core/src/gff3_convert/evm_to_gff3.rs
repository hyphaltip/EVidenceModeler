//! Convert EVM output format to GFF3.
//!
//! Faithful port of `EVM_to_GFF3.pl` + `Gene_obj::to_GFF3_format`
//! (`convert_EVM_outputs_to_GFF3.pl` just drives `EVM_to_GFF3.pl` per
//! partition; its `join_intronic_preds`/`combine_predictions` subs are dead
//! code — never invoked — so they are intentionally not ported).

use crate::io::partitions::PartitionEntry;
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::{BufRead, BufWriter, Write};

/// EVM exon start-frame [1..=6] → GFF3 CDS phase [0,1,2].
///
/// Composes the two Perl conversions: `EVM_to_GFF3.pl`'s `%phase_conversion`
/// (1→0,2→1,3→2,4→0,5→1,6→2) followed by `Gene_obj::to_GFF3_format`'s GFF3
/// phase swap (1↔2, 0 unchanged). Net: {1,4}→0, {2,5}→2, {3,6}→1.
fn evm_startframe_to_gff3_phase(frame: u8) -> u8 {
    match frame {
        1 | 4 => 0,
        2 | 5 => 2,
        3 | 6 => 1,
        _ => 0,
    }
}

/// URI-escape matching Perl `URI::Escape::uri_escape` defaults: everything
/// except the RFC 3986 unreserved set `A-Za-z0-9-_.~` is percent-encoded.
fn uri_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for &b in s.as_bytes() {
        let unreserved = b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~');
        if unreserved {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{:02X}", b));
        }
    }
    out
}

/// One parsed EVM model: its exon rows `(end5, end3, start_frame)` in file
/// order, and the evidence-type tag for the GFF3 `source` column.
struct Model {
    coords: Vec<(u32, u32, u8)>,
    ev_type: &'static str,
}

/// Convert a single EVM output file to GFF3 format, writing to `output_path`.
pub fn evm_output_to_gff3(evm_file: &str, contig_id: &str, output_path: &str) -> Result<()> {
    let f =
        fs::File::open(evm_file).with_context(|| format!("Cannot open EVM output {}", evm_file))?;
    let bytes = evm_reader_to_gff3(std::io::BufReader::new(f), contig_id)?;
    // One write per file instead of one per GFF3 line.
    fs::write(output_path, &bytes)
        .with_context(|| format!("Cannot create GFF3 output {}", output_path))
}

/// Parse EVM output text and return the GFF3 bytes for `contig_id`.
fn evm_reader_to_gff3<R: BufRead>(reader: R, contig_id: &str) -> Result<Vec<u8>> {
    let mut out: Vec<u8> = Vec::new();

    // Mirror EVM_to_GFF3.pl's scan: model_id starts at 1 and increments on each
    // blank line; '#' header lines set the current model's ev_type; numeric
    // exon rows (6 tab cols, INTRON excluded) accumulate onto the current model.
    let mut model_id = 1u32;
    let mut models: BTreeMap<u32, Model> = BTreeMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('!') {
            continue;
        } // '^!' comment line

        if line.starts_with('#') {
            // Perl: a '#' line not matching /EVM/ is skipped; otherwise sets ev_type.
            if line.contains("EVM") {
                let ev_type = if line.contains("ELIMINATED") {
                    "EVM_elm"
                } else {
                    "EVM"
                };
                models
                    .entry(model_id)
                    .or_insert_with(|| Model {
                        coords: Vec::new(),
                        ev_type,
                    })
                    .ev_type = ev_type;
            }
            continue;
        }

        // Non-word line (blank) → end of current model.
        if !line.chars().any(|c| c.is_alphanumeric()) {
            model_id += 1;
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        // Perl requires exactly 6 columns with numeric coords.
        if cols.len() == 6 {
            if let (Ok(e5), Ok(e3), Ok(frame)) = (
                cols[0].parse::<u32>(),
                cols[1].parse::<u32>(),
                cols[3].parse::<u8>(),
            ) {
                if cols[2] != "INTRON" {
                    models
                        .entry(model_id)
                        .or_insert_with(|| Model {
                            coords: Vec::new(),
                            ev_type: "EVM",
                        })
                        .coords
                        .push((e5, e3, frame));
                }
            }
        }
    }

    for (mid, model) in &models {
        if model.coords.is_empty() {
            continue;
        }
        emit_model(&mut out, contig_id, *mid, model)?;
    }
    Ok(out)
}

/// Emit one gene model (gene + mRNA + interleaved exon/CDS rows) in GFF3.
fn emit_model(out: &mut dyn Write, contig: &str, model_id: u32, model: &Model) -> Result<()> {
    let ev_type = model.ev_type;
    let coords = &model.coords;

    let gene_lend = coords.iter().map(|&(e5, e3, _)| e5.min(e3)).min().unwrap();
    let gene_rend = coords.iter().map(|&(e5, e3, _)| e5.max(e3)).max().unwrap();
    // Strand from the first exon row (Perl Gene_obj orientation).
    let strand = if coords[0].0 <= coords[0].1 { '+' } else { '-' };

    let tu_id = format!("evm.TU.{}.{}", contig, model_id);
    let model_feat = format!("evm.model.{}.{}", contig, model_id);
    let com_name = uri_escape(&format!("EVM prediction {}.{}", contig, model_id));

    // gene + mRNA
    writeln!(
        out,
        "{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tID={};Name={}",
        contig, ev_type, gene_lend, gene_rend, strand, tu_id, com_name
    )?;
    writeln!(
        out,
        "{}\t{}\tmRNA\t{}\t{}\t.\t{}\t.\tID={};Parent={};Name={}",
        contig, ev_type, gene_lend, gene_rend, strand, model_feat, tu_id, com_name
    )?;

    // Exons in transcription order: ascending genomic position for '+',
    // descending for '-' (matches Gene_obj::get_exons 5'→3' ordering).
    let mut ordered: Vec<(u32, u32, u8)> = coords.clone();
    ordered.sort_by_key(|&(e5, e3, _)| e5.min(e3));
    if strand == '-' {
        ordered.reverse();
    }

    for (i, &(e5, e3, frame)) in ordered.iter().enumerate() {
        let (cs, ce) = (e5.min(e3), e5.max(e3));
        let exon_num = i + 1;
        writeln!(
            out,
            "{}\t{}\texon\t{}\t{}\t.\t{}\t.\tID={}.exon{};Parent={}",
            contig, ev_type, cs, ce, strand, model_feat, exon_num, model_feat
        )?;
        let phase = evm_startframe_to_gff3_phase(frame);
        writeln!(
            out,
            "{}\t{}\tCDS\t{}\t{}\t.\t{}\t{}\tID=cds.{};Parent={}",
            contig, ev_type, cs, ce, strand, phase, model_feat, model_feat
        )?;
    }
    writeln!(out)?; // blank spacer between genes
    Ok(())
}

/// Entries converted per parallel batch. Batches bound memory when a
/// concatenated output is requested and give regular progress lines.
const CHUNK_SIZE: usize = 4096;

/// Options for [`convert_all_to_gff3_opts`].
#[derive(Debug, Clone, Default)]
pub struct ConvertOptions {
    /// Write no `.gff3` for a zero-byte `evm.out` and remove a stale one.
    pub skip_empty: bool,
    /// Also write every converted GFF3, in partitions-listing order, to this file.
    pub concat_out: Option<String>,
}

/// Convert EVM outputs for all entries to GFF3.
pub fn convert_all_to_gff3(entries: &[PartitionEntry], output_file_name: &str) -> Result<()> {
    convert_all_to_gff3_opts(entries, output_file_name, &ConvertOptions::default())
}

/// Like [`convert_all_to_gff3`]. With `skip_empty`, an `evm.out` of zero bytes
/// (no models, or a failed EVM run that left only the stdout redirect) writes
/// no `.gff3`, and any stale `.gff3` from an earlier run is removed. A missing
/// `evm.out` stays an error either way.
pub fn convert_all_to_gff3_with(
    entries: &[PartitionEntry],
    output_file_name: &str,
    skip_empty: bool,
) -> Result<()> {
    let opts = ConvertOptions {
        skip_empty,
        concat_out: None,
    };
    convert_all_to_gff3_opts(entries, output_file_name, &opts)
}

enum Outcome {
    Skipped,
    /// GFF3 bytes; empty unless the caller asked for them.
    Converted(Vec<u8>),
}

fn convert_one(
    accession: &str,
    base_dir: &str,
    output_file_name: &str,
    skip_empty: bool,
    keep_bytes: bool,
) -> Result<Outcome> {
    let evm_file = format!("{}/{}", base_dir, output_file_name);
    let gff3_file = format!("{}/{}.gff3", base_dir, output_file_name);
    // One open serves the existence check, the size check and the read.
    let f = fs::File::open(&evm_file)
        .with_context(|| format!("Cannot open EVM output {}", evm_file))?;
    if skip_empty {
        let len = f
            .metadata()
            .with_context(|| format!("Cannot stat EVM output {}", evm_file))?
            .len();
        if len == 0 {
            match fs::remove_file(&gff3_file) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                Err(e) => {
                    return Err(e).with_context(|| format!("Cannot remove stale {}", gff3_file))
                }
            }
            return Ok(Outcome::Skipped);
        }
    }
    let bytes = evm_reader_to_gff3(std::io::BufReader::new(f), accession)?;
    fs::write(&gff3_file, &bytes)
        .with_context(|| format!("Cannot create GFF3 output {}", gff3_file))?;
    Ok(Outcome::Converted(if keep_bytes {
        bytes
    } else {
        Vec::new()
    }))
}

/// Convert all entries in parallel (rayon global pool). Output is identical to
/// a serial run: each contig writes only its own files, and the optional
/// concatenated file is written in partitions-listing order.
pub fn convert_all_to_gff3_opts(
    entries: &[PartitionEntry],
    output_file_name: &str,
    opts: &ConvertOptions,
) -> Result<()> {
    // One conversion per accession, in first-seen listing order; a later
    // duplicate accession replaces the base_dir, as the old HashMap did.
    let mut index: HashMap<&str, usize> = HashMap::new();
    let mut work: Vec<(&str, &str)> = Vec::new();
    for entry in entries {
        match index.get(entry.accession.as_str()) {
            Some(&i) => work[i].1 = entry.base_dir.as_str(),
            None => {
                index.insert(entry.accession.as_str(), work.len());
                work.push((entry.accession.as_str(), entry.base_dir.as_str()));
            }
        }
    }

    let mut concat = match &opts.concat_out {
        Some(path) => Some(BufWriter::new(
            fs::File::create(path)
                .with_context(|| format!("Cannot create concatenated GFF3 {}", path))?,
        )),
        None => None,
    };
    let keep_bytes = concat.is_some();

    let total = work.len();
    let start = std::time::Instant::now();
    let mut skipped = 0usize;
    let mut done = 0usize;
    log::info!(
        "Converting {} EVM outputs to GFF3 on {} threads",
        total,
        rayon::current_num_threads()
    );
    for chunk in work.chunks(CHUNK_SIZE) {
        let outcomes: Vec<Outcome> = chunk
            .par_iter()
            .map(|&(acc, dir)| convert_one(acc, dir, output_file_name, opts.skip_empty, keep_bytes))
            .collect::<Result<Vec<_>>>()?;
        for outcome in outcomes {
            match outcome {
                Outcome::Skipped => skipped += 1,
                Outcome::Converted(bytes) => {
                    if let Some(w) = concat.as_mut() {
                        w.write_all(&bytes)?;
                    }
                }
            }
        }
        done += chunk.len();
        let rate = done as f64 / start.elapsed().as_secs_f64().max(1e-9);
        log::info!(
            "GFF3 conversion: {}/{} done, {:.1} files/s, ~{:.0} s remaining",
            done,
            total,
            rate,
            (total - done) as f64 / rate
        );
    }
    if let Some(mut w) = concat {
        w.flush().context("Cannot write concatenated GFF3")?;
    }
    if opts.skip_empty {
        log::info!("Skipped {} of {} empty EVM outputs", skipped, total);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(acc: &str, dir: &std::path::Path) -> PartitionEntry {
        PartitionEntry {
            accession: acc.to_string(),
            base_dir: dir.to_string_lossy().into_owned(),
            is_partitioned: false,
            partition_dir: None,
        }
    }

    fn scratch(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("evm_gff3_{}_{}", name, std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    const ONE_MODEL: &str = "# EVM prediction: Mode:STANDARD S-ratio: 1.00 1-100 orient(+) score(1.0) noncoding_equivalent(0.0) raw_coding_score(1.0) raw_noncoding_score(0.0) offset(0)\n1\t99\tINITIAL\t1\t+\t{}\n\n";

    #[test]
    fn skip_empty_skips_zero_byte_output_and_removes_stale_gff3() {
        let d = scratch("skip");
        fs::write(d.join("evm.out"), "").unwrap();
        fs::write(d.join("evm.out.gff3"), "stale\n").unwrap();
        convert_all_to_gff3_with(&[entry("c1", &d)], "evm.out", true).unwrap();
        assert!(!d.join("evm.out.gff3").exists());
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn default_still_writes_gff3_for_empty_output() {
        let d = scratch("default");
        fs::write(d.join("evm.out"), "").unwrap();
        convert_all_to_gff3_with(&[entry("c1", &d)], "evm.out", false).unwrap();
        assert!(d.join("evm.out.gff3").exists());
        let _ = fs::remove_dir_all(&d);
    }

    #[test]
    fn skip_empty_still_converts_nonempty_and_errors_on_missing() {
        let d = scratch("nonempty");
        fs::write(d.join("evm.out"), ONE_MODEL).unwrap();
        convert_all_to_gff3_with(&[entry("c1", &d)], "evm.out", true).unwrap();
        let gff = fs::read_to_string(d.join("evm.out.gff3")).unwrap();
        assert!(gff.contains("\tgene\t"));

        let m = scratch("missing");
        assert!(convert_all_to_gff3_with(&[entry("c2", &m)], "evm.out", true).is_err());
        let _ = fs::remove_dir_all(&d);
        let _ = fs::remove_dir_all(&m);
    }
}
