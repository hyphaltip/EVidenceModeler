//! Intergenic region scoring.

use std::collections::HashMap;
use crate::types::genome::MaskVec;
use crate::types::evidence::EvWeightMap;
use crate::io::gff3::Gff3Record;

pub type IntergenicScores = Vec<f64>;

/// Compute per-base intergenic scores from gene-prediction spans.
///
/// Mirrors the Perl `populate_intergenic_regions`: for each ab-initio
/// prediction type, the per-base score in the gap between two neighbouring
/// genes of that type is incremented by that type's weight. Only
/// `ABINITIO_PREDICTION` types contribute; all genes on both strands are
/// considered. Masked positions are skipped. The genome ends [0,0] and
/// [seq_len,seq_len] are added as sentinels so the flanking regions count.
///
/// The `INTERGENIC_SCORE_ADJUST_FACTOR` is folded into the per-base value
/// here; with the default factor of 1.0 this is identical to applying it in
/// `calc_intergenic_score` as the Perl does.
pub fn populate_intergenic_scores(
    seq_len: usize,
    gene_pred_records: &[Gff3Record],
    ev_weights: &EvWeightMap,
    mask: &MaskVec,
    adjust_factor: f64,
) -> IntergenicScores {
    let mut ig = vec![0.0f64; seq_len + 2];

    // Group CDS spans per (ev_type, model id), forward-coordinate span.
    // model_spans[ev_type][model_id] = (min_coord, max_coord)
    let mut per_type: HashMap<String, HashMap<String, (u32, u32)>> = HashMap::new();

    for rec in gene_pred_records {
        if rec.feature != "CDS" { continue; }
        let entry = match ev_weights.get(&rec.source) {
            Some(e) => e,
            None => continue,
        };
        if !entry.ev_class.is_abinitio() { continue; }

        let parent = rec.raw_attributes.split("Parent=").nth(1)
            .and_then(|s| s.split(|c| c == ';' || c == ' ').next())
            .unwrap_or("")
            .to_string();
        if parent.is_empty() { continue; }

        let lend = rec.start.min(rec.end);
        let rend = rec.start.max(rec.end);
        let models = per_type.entry(rec.source.clone()).or_default();
        let span = models.entry(parent).or_insert((lend, rend));
        span.0 = span.0.min(lend);
        span.1 = span.1.max(rend);
    }

    for (ev_type, models) in &per_type {
        let weight = ev_weights[ev_type].weight * adjust_factor;

        // Collect gene spans plus genome-boundary sentinels, sorted by lend.
        let mut spans: Vec<(u32, u32)> = models.values().copied().collect();
        spans.push((0, 0));
        spans.push((seq_len as u32, seq_len as u32));
        spans.sort_by_key(|&(l, _)| l);

        for w in spans.windows(2) {
            let lend_intergenic = w[0].1;
            let rend_intergenic = w[1].0;
            if lend_intergenic > rend_intergenic { continue; }
            let mut j = lend_intergenic + 1;
            while j + 1 <= rend_intergenic {
                if !mask.get(j as usize) {
                    ig[j as usize] += weight;
                }
                j += 1;
            }
        }
    }

    ig
}

/// Compute the sum of intergenic scores over [lend, rend] (inclusive, 1-based).
pub fn calc_intergenic_score(scores: &IntergenicScores, lend: u32, rend: u32) -> f64 {
    if lend > rend { return 0.0; }
    let mut s = 0.0;
    let max = (scores.len() - 1) as u32;
    let l = lend.max(1);
    let r = rend.min(max);
    for i in l..=r {
        s += scores[i as usize];
    }
    s
}

/// Augment intergenic scores using start/stop peak positions.
/// Peaks provide evidence of gene boundaries; positions near a peak
/// receive an extra boost proportional to the peak score.
pub fn augment_intergenic_from_peaks(
    ig: &mut IntergenicScores,
    peaks: &[(u32, f64)], // (position, score)
    window: u32,
) {
    for &(pos, score) in peaks {
        let lend = pos.saturating_sub(window);
        let rend = (pos + window).min((ig.len() - 1) as u32);
        for i in lend..=rend {
            ig[i as usize] += score;
        }
    }
}

/// Get contiguous intergenic regions between adjacent predictions.
/// Returns Vec<(lend, rend)> pairs.
pub fn get_intergenic_regions(
    predictions_sorted_by_lend: &[(u32, u32)], // (lend, rend)
) -> Vec<(u32, u32)> {
    let mut regions = Vec::new();
    for pair in predictions_sorted_by_lend.windows(2) {
        let rend_prev = pair[0].1;
        let lend_next = pair[1].0;
        if lend_next > rend_prev + 1 {
            regions.push((rend_prev + 1, lend_next - 1));
        }
    }
    regions
}
