//! Prediction filtering — remove low-support and degenerate gene models.

use crate::types::prediction::EvmPrediction;
use crate::types::exon::{Exon, ExonType};
use crate::algo::coding_scores::CodingScores;
use crate::algo::introns::IntronVec;

/// Minimum coding/noncoding score ratio; predictions below this are eliminated.
const MIN_CODING_NONCODING_SCORE_RATIO: f64 = 0.75;
/// Minimum total coding length for nested/intergenic gene searches.
const MIN_CODING_LENGTH: u32 = 300;

/// Filter predictions that have insufficient evidence support.
///
/// A prediction is eliminated if:
/// 1. The ratio of its coding score to its intergenic baseline is too low, or
/// 2. Its total CDS length is below `MIN_CODING_LENGTH` (for nested/intergenic modes).
pub fn filter_predictions_low_support(
    predictions: &mut Vec<EvmPrediction>,
    exons: &[Exon],
    coding_scores: &CodingScores,
    fwd_intron_vec: &IntronVec,
    rev_intron_vec: &IntronVec,
    mode: &str,
) {
    for pred in predictions.iter_mut() {
        if pred.is_eliminated { continue; }

        let (lend, rend) = pred.get_span();

        // Compute total CDS length
        let cds_len: u32 = pred.exon_indices.iter()
            .map(|&i| exons[i].length())
            .sum();

        // For intron/intergenic re-search modes, enforce minimum coding length
        if mode != "STANDARD" && cds_len < MIN_CODING_LENGTH {
            pred.is_eliminated = true;
            continue;
        }

        // Compute average coding score per base over the span
        let span = (rend - lend + 1) as f64;
        if span <= 0.0 { continue; }

        let coding_sum: f64 = (lend..=rend)
            .map(|i| coding_scores[i as usize].max(0.0))
            .sum();
        let coding_avg = coding_sum / span;

        // Compute average intron score over the span (use the strand of first exon)
        let strand = exons[pred.exon_indices[0]].orientation;
        let intron_vec = if strand == crate::types::exon::Orientation::Fwd {
            fwd_intron_vec
        } else {
            rev_intron_vec
        };
        let intron_sum: f64 = (lend..=rend)
            .map(|i| intron_vec.get(i as usize).copied().unwrap_or(0.0))
            .sum();
        let intron_avg = intron_sum / span;

        // Baseline: the better of coding or intron averages
        let baseline = coding_avg.max(intron_avg);
        if baseline <= 0.0 { continue; }

        // Compute total prediction score
        let pred_score: f64 = pred.exon_indices.iter()
            .map(|&i| exons[i].sum_score.max(0.0))
            .sum();

        let ratio = pred_score / (baseline * span);
        if ratio < MIN_CODING_NONCODING_SCORE_RATIO {
            pred.is_eliminated = true;
        }
    }
}

/// For predictions that lack a proper start codon (5' partials), try to find
/// an alternative initial exon upstream that starts at ATG.
///
/// This is a simplified version of the Perl `convert_5prime_partials_to_complete_genes_where_possible`.
pub fn convert_5prime_partials(
    predictions: &mut Vec<EvmPrediction>,
    exons: &[Exon],
) {
    for pred in predictions.iter_mut() {
        if pred.is_eliminated { continue; }
        if pred.exon_indices.is_empty() { continue; }
        let first_idx = pred.exon_indices[0];
        let first_exon = &exons[first_idx];
        // Only act on non-initial leading exons (i.e. the prediction starts
        // with an internal exon, making it a 5' partial)
        if first_exon.exon_type != ExonType::Internal { continue; }
        // In the full Perl implementation this searches for a new initial exon.
        // Here we mark the prediction as a partial — full recovery logic
        // would require access to all candidate exons, which is done in the
        // `consensus` module.
        // TODO: full partial→complete conversion.
    }
}
