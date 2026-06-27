//! Intergenic region scoring.

use crate::types::genome::MaskVec;
use crate::types::exon::Exon;
use crate::algo::coding_scores::CodingScores;

pub type IntergenicScores = Vec<f64>;

/// Compute per-base intergenic scores by inverting the coding score vector.
///
/// Regions with zero coding score receive a positive intergenic score
/// equal to `1.0 * INTERGENIC_SCORE_ADJUST_FACTOR`; regions with coding
/// evidence receive 0.  Masked positions receive 0.
pub fn populate_intergenic_scores(
    seq_len: usize,
    coding_scores: &CodingScores,
    mask: &MaskVec,
    adjust_factor: f64,
) -> IntergenicScores {
    let mut ig = vec![0.0f64; seq_len + 2];
    for i in 1..=seq_len {
        if mask.get(i) { continue; }
        if coding_scores[i] <= 0.0 {
            ig[i] = 1.0 * adjust_factor;
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
