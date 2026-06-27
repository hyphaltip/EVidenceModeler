//! Recursive consensus gene prediction — trellis + recursion on tail/intergenic regions.

use std::collections::{HashSet, HashMap};
use anyhow::Result;
use crate::types::exon::{Exon, ExonPhase};
use crate::types::prediction::{EvmPrediction, PredMode};
use crate::algo::trellis::{build_trellis, traverse_path};
use crate::algo::filter::filter_predictions_low_support;
use crate::algo::introns::{IntronScoreMap, IntronVec};
use crate::algo::intergenic::{IntergenicScores, get_intergenic_regions};
use crate::algo::coding_scores::CodingScores;

/// Parameters controlling the recursive search.
pub struct ConsensusParams<'a> {
    pub exons: &'a mut Vec<Exon>,
    pub introns_to_score: &'a IntronScoreMap,
    pub intergenic_scores: &'a IntergenicScores,
    pub acceptable_linkages: &'a HashSet<(String, String)>,
    pub phased_connections: &'a HashSet<(String, String)>,
    pub intergenic_connections: &'a HashSet<(String, String)>,
    pub frame_pairs: &'a HashSet<(ExonPhase, ExonPhase)>,
    pub stop_codons: &'a [[u8; 3]],
    pub coding_scores: &'a CodingScores,
    pub fwd_intron_vec: &'a IntronVec,
    pub rev_intron_vec: &'a IntronVec,
    pub max_prev_exons_compare: usize,
    pub min_intergenic_size_on_re_search: u32,
    pub min_gene_length_size_on_re_search: u32,
    pub report_elm: bool,
    pub recursion_limit: usize,
}

/// Generate consensus gene predictions for [range_lend, range_rend].
///
/// Outputs formatted prediction text (matching the Perl output format) to `output`.
pub fn generate_consensus_gene_predictions(
    range_lend: u32,
    range_rend: u32,
    mode: PredMode,
    params: &mut ConsensusParams,
    recursion_count: &mut usize,
    output: &mut Vec<String>,
) -> Result<()> {
    *recursion_count += 1;

    let region_len = range_rend.saturating_sub(range_lend) + 1;
    if region_len < params.min_gene_length_size_on_re_search && mode != PredMode::Standard {
        *recursion_count -= 1;
        return Ok(());
    }

    if *recursion_count > params.recursion_limit {
        log::warn!("Recursion limit ({}) exceeded", params.recursion_limit);
        *recursion_count -= 1;
        return Ok(());
    }

    // Collect exons within range
    let exon_indices_in_range: Vec<usize> = {
        let exons_ref = &*params.exons;
        (0..exons_ref.len())
            .filter(|&i| {
                let (l, r) = exons_ref[i].coords_sorted();
                l >= range_lend && r <= range_rend
            })
            .collect()
    };

    if exon_indices_in_range.is_empty() {
        log::debug!("No exons in range {}-{}", range_lend, range_rend);
        *recursion_count -= 1;
        return Ok(());
    }

    // Build a local copy of exons for the trellis (subset in range)
    let mut local_exons: Vec<Exon> = exon_indices_in_range.iter()
        .map(|&i| params.exons[i].clone())
        .collect();
    local_exons.sort_by_key(|e| e.end5);

    let top_idx = build_trellis(
        &mut local_exons,
        range_lend,
        range_rend,
        params.acceptable_linkages,
        params.phased_connections,
        params.intergenic_connections,
        params.frame_pairs,
        params.introns_to_score,
        params.intergenic_scores,
        params.stop_codons,
        params.max_prev_exons_compare,
    );

    let top_idx = match top_idx {
        Some(i) => i,
        None => { *recursion_count -= 1; return Ok(()); }
    };

    let mut predictions = traverse_path(&local_exons, top_idx);
    if predictions.is_empty() {
        *recursion_count -= 1;
        return Ok(());
    }

    // Filter low-support predictions
    filter_predictions_low_support(
        &mut predictions,
        &local_exons,
        params.coding_scores,
        params.fwd_intron_vec,
        params.rev_intron_vec,
        mode.as_str(),
    );

    let preds_remain: Vec<&EvmPrediction> = predictions.iter()
        .filter(|p| !p.is_eliminated)
        .collect();

    if preds_remain.is_empty() && !params.report_elm {
        *recursion_count -= 1;
        return Ok(());
    }

    // Determine span of predictions
    let (pred_span_lend, pred_span_rend) = if !preds_remain.is_empty() {
        let l = preds_remain.iter().map(|p| p.lend).min().unwrap_or(range_lend);
        let r = preds_remain.iter().map(|p| p.rend).max().unwrap_or(range_rend);
        (l, r)
    } else {
        predictions.iter().fold((range_rend, range_lend), |(l, r), p| (l.min(p.lend), r.max(p.rend)))
    };

    // Emit predictions
    for pred in &predictions {
        if pred.is_eliminated && !params.report_elm { continue; }
        let text = format_prediction(pred, &local_exons, &mode, *recursion_count);
        output.push(text);
    }

    // Recursion: tail regions
    if mode != PredMode::Intron {
        let left_len = pred_span_lend.saturating_sub(range_lend);
        if left_len >= params.min_intergenic_size_on_re_search {
            generate_consensus_gene_predictions(
                range_lend, pred_span_lend - 1, mode.clone(), params, recursion_count, output,
            )?;
        }
        let right_len = range_rend.saturating_sub(pred_span_rend);
        if right_len >= params.min_intergenic_size_on_re_search {
            generate_consensus_gene_predictions(
                pred_span_rend + 1, range_rend, mode.clone(), params, recursion_count, output,
            )?;
        }

        // Recursion: intergenic regions between predictions
        if params.min_gene_length_size_on_re_search > 0 {
            let spans: Vec<(u32, u32)> = predictions.iter()
                .filter(|p| !p.is_eliminated)
                .map(|p| (p.lend, p.rend))
                .collect();
            for (ig_l, ig_r) in get_intergenic_regions(&spans) {
                let ig_len = ig_r.saturating_sub(ig_l) + 1;
                if ig_len >= params.min_gene_length_size_on_re_search {
                    generate_consensus_gene_predictions(
                        ig_l, ig_r, mode.clone(), params, recursion_count, output,
                    )?;
                }
            }
        }
    }

    *recursion_count -= 1;
    Ok(())
}

/// Format a prediction as EVM output text.
fn format_prediction(
    pred: &EvmPrediction,
    exons: &[Exon],
    mode: &PredMode,
    recursion_count: usize,
) -> String {
    use crate::types::exon::{ExonType, exon_phase_to_gff_phase};
    let prefix = if pred.is_eliminated { "#ELIMINATED EVM prediction" } else { "#EVM prediction" };
    let mut s = format!(
        "{} mode:{} span:{}-{} [R{}]\n",
        prefix, mode.as_str(), pred.lend, pred.rend, recursion_count
    );

    let mut ev_info = String::new();
    for &idx in &pred.exon_indices {
        let exon = &exons[idx];
        let (end5, end3) = (exon.end5, exon.end3);
        let phase = exon_phase_to_gff_phase(exon.start_frame);
        let etype = exon.exon_type.as_str();
        // Collect evidence strings
        let ev_str: Vec<String> = exon.evidence.iter()
            .map(|(acc, et)| format!("{}/{}", et, acc))
            .collect();
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            end5, end3, etype, phase,
            ev_str.join(";")
        ));
    }

    s
}
