//! Convert GFF3 gene models to BED format.
//!
//! Replaces gene_gff3_to_bed.pl.

use std::collections::HashMap;
use anyhow::Result;
use crate::io::gff3::read_gff3_file;

/// Convert a GFF3 file containing gene models to BED format.
///
/// Each gene is represented as a single BED12 line.
pub fn gff3_to_bed(gff3_path: &str) -> Result<Vec<String>> {
    let records = read_gff3_file(gff3_path)?;

    // Collect per-gene info
    let mut gene_records: HashMap<String, Vec<(u32, u32, char, String)>> = HashMap::new();
    let mut gene_contigs: HashMap<String, String> = HashMap::new();
    let mut gene_strands: HashMap<String, char> = HashMap::new();

    for rec in &records {
        if rec.feature == "gene" {
            let id = rec.attr("ID").unwrap_or("").to_string();
            gene_contigs.insert(id.clone(), rec.seqid.clone());
            gene_strands.insert(id, rec.strand);
        }
        if rec.feature == "CDS" {
            let parent = rec.attr("Parent").unwrap_or("").to_string();
            gene_records.entry(parent.clone()).or_default().push((
                rec.start, rec.end, rec.strand, rec.seqid.clone()
            ));
        }
    }

    let mut bed_lines = Vec::new();

    for (gene_id, exons) in &gene_records {
        if exons.is_empty() { continue; }
        let contig = &exons[0].3;
        let strand = exons[0].2;
        let chrom_start: u32 = exons.iter().map(|&(s, _, _, _)| s).min().unwrap_or(0) - 1; // BED is 0-based
        let chrom_end: u32 = exons.iter().map(|&(_, e, _, _)| e).max().unwrap_or(0);

        let mut sorted_exons: Vec<(u32, u32)> = exons.iter()
            .map(|&(s, e, _, _)| (s - 1, e)) // 0-based
            .collect();
        sorted_exons.sort();
        sorted_exons.dedup();

        let block_count = sorted_exons.len();
        let block_sizes: Vec<String> = sorted_exons.iter()
            .map(|&(s, e)| (e - s).to_string())
            .collect();
        let block_starts: Vec<String> = sorted_exons.iter()
            .map(|&(s, _)| (s - chrom_start).to_string())
            .collect();

        bed_lines.push(format!(
            "{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            contig,
            chrom_start,
            chrom_end,
            gene_id,
            strand,
            chrom_start,
            chrom_end,
            "0,0,0",
            block_count,
            block_sizes.join(","),
            block_starts.join(","),
        ));
    }

    // Sort BED lines by chromosome then start position
    bed_lines.sort_by(|a, b| {
        let a_cols: Vec<&str> = a.splitn(3, '\t').collect();
        let b_cols: Vec<&str> = b.splitn(3, '\t').collect();
        a_cols[0].cmp(b_cols[0])
            .then(a_cols[1].parse::<u32>().unwrap_or(0)
                .cmp(&b_cols[1].parse::<u32>().unwrap_or(0)))
    });

    Ok(bed_lines)
}
