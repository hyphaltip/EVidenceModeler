//! Dynamic-programming trellis: finding the highest-scoring path through exons.

use std::collections::HashSet;
use crate::types::exon::{Exon, ExonType, ExonPhase};
use crate::types::prediction::EvmPrediction;
use crate::algo::introns::IntronScoreMap;
use crate::algo::intergenic::{IntergenicScores, calc_intergenic_score};
use crate::algo::phases::is_stop_codon;

/// Result of compatibility check between two exons.
pub enum CompatResult {
    /// Incompatible — cannot link A → B.
    Incompatible,
    /// Compatible; the score bonus (intron score or intergenic score) to add.
    Compatible(f64),
}

/// Test whether exon A (left / upstream) can be followed by exon B (right / downstream)
/// in the trellis, and compute the join score.
///
/// Returns `Compatible(score)` on success, `Incompatible` otherwise.
pub fn are_compatible_exons(
    exon_a: &Exon,
    exon_b: &Exon,
    acceptable_linkages: &HashSet<(String, String)>,
    phased_connections: &HashSet<(String, String)>,
    intergenic_connections: &HashSet<(String, String)>,
    frame_pairs: &HashSet<(ExonPhase, ExonPhase)>,
    introns_to_score: &IntronScoreMap,
    intergenic_scores: &IntergenicScores,
    stop_codons: &[[u8; 3]],
) -> CompatResult {
    let key_a = exon_a.type_orient_key();
    let key_b = exon_b.type_orient_key();

    // Check linkage is allowed
    if !acceptable_linkages.contains(&(key_a.clone(), key_b.clone())) {
        return CompatResult::Incompatible;
    }

    let (a_lend, a_rend) = exon_a.coords_sorted();
    let (b_lend, b_rend) = exon_b.coords_sorted();

    // No overlap allowed
    if a_lend <= b_rend && a_rend >= b_lend {
        return CompatResult::Incompatible;
    }

    if phased_connections.contains(&(key_a.clone(), key_b.clone())) {
        // Check intron validity
        let (intron_end5, intron_end3) = if key_a.ends_with('-') {
            // Reverse strand: intron is upstream of A in genomic terms
            (b_lend - 1, a_rend + 2)
        } else {
            (a_rend + 1, b_lend - 2)
        };

        let intron_key = format!("{}_{}", intron_end5, intron_end3);
        let intron_score = match introns_to_score.get(&intron_key) {
            Some(&s) => s,
            None => return CompatResult::Incompatible,
        };

        // Check frame compatibility
        let (before, after) = if key_a.ends_with('-') { (exon_b, exon_a) } else { (exon_a, exon_b) };
        if !frame_pairs.contains(&(before.end_frame, after.start_frame)) {
            return CompatResult::Incompatible;
        }

        // Check no stop codon created across the junction
        let end_frame = before.end_frame % 3;
        let seq_junction: Vec<u8> = before.right_seq_boundary.iter()
            .chain(after.left_seq_boundary.iter())
            .copied()
            .collect();
        let potential_stop = match end_frame {
            1 => seq_junction.get(1..4),
            2 => seq_junction.get(0..3),
            _ => None, // frame 0 / 3: no partial codon at junction
        };
        if let Some(codon) = potential_stop {
            if is_stop_codon(codon, stop_codons) {
                return CompatResult::Incompatible;
            }
        }

        CompatResult::Compatible(intron_score)
    } else if intergenic_connections.contains(&(key_a.clone(), key_b.clone())) {
        let score = calc_intergenic_score(intergenic_scores, a_rend + 1, b_lend - 1);
        CompatResult::Compatible(score)
    } else {
        CompatResult::Incompatible
    }
}

/// Score the connection between a boundary node and an exon (or vice-versa).
/// Boundary conditions always return 0 (no extra score).
pub fn score_boundary_condition(_boundary: &Exon, _other: &Exon) -> f64 {
    0.0
}

/// Build the trellis over `exons` restricted to [range_lend, range_rend].
///
/// Exons must be sorted by end5 ascending before calling this function.
/// Returns the index of the highest-scoring exon (or None if empty).
pub fn build_trellis(
    exons: &mut Vec<Exon>,
    range_lend: u32,
    range_rend: u32,
    acceptable_linkages: &HashSet<(String, String)>,
    phased_connections: &HashSet<(String, String)>,
    intergenic_connections: &HashSet<(String, String)>,
    frame_pairs: &HashSet<(ExonPhase, ExonPhase)>,
    introns_to_score: &IntronScoreMap,
    intergenic_scores: &IntergenicScores,
    stop_codons: &[[u8; 3]],
    max_prev_exons_compare: usize,
) -> Option<usize> {
    if exons.is_empty() { return None; }

    // Reset link and sum_score
    for exon in exons.iter_mut() {
        exon.sum_score = exon.base_score;
        exon.link = None;
    }

    // Add boundary sentinel nodes
    let _left_bound_idx = exons.len();
    let mut left_bound = Exon::new(range_lend, range_lend);
    left_bound.exon_type = ExonType::Bound;
    left_bound.start_frame = 1;
    left_bound.end_frame = 1;
    exons.push(left_bound);

    let _right_bound_idx = exons.len();
    let mut right_bound = Exon::new(range_rend, range_rend);
    right_bound.exon_type = ExonType::Bound;
    right_bound.start_frame = 1;
    right_bound.end_frame = 1;
    exons.push(right_bound);

    let num_exons = exons.len();
    let mut highest_score = 0.0f64;
    let mut highest_idx = 0usize;

    for i in 1..num_exons {
        let base_score_i = exons[i].base_score;
        let mut best_sum = exons[i].sum_score;

        let mut compare_count = 0usize;
        let mut found_compatible = false;

        let i_type_bound = exons[i].exon_type == ExonType::Bound;

        let j_start = if i == 0 { 0 } else { i - 1 };
        let mut j = j_start as isize;

        while j >= 0
            && (compare_count < max_prev_exons_compare || !found_compatible)
        {
            let ji = j as usize;
            compare_count += 1;
            j -= 1;

            let j_type_bound = exons[ji].exon_type == ExonType::Bound;

            let join_score = if i_type_bound || j_type_bound {
                // boundary condition: always 0
                Some(score_boundary_condition(&exons[ji], &exons[i]))
            } else {
                match are_compatible_exons(
                    &exons[ji],
                    &exons[i],
                    acceptable_linkages,
                    phased_connections,
                    intergenic_connections,
                    frame_pairs,
                    introns_to_score,
                    intergenic_scores,
                    stop_codons,
                ) {
                    CompatResult::Compatible(s) => Some(s),
                    CompatResult::Incompatible => None,
                }
            };

            if let Some(js) = join_score {
                found_compatible = true;
                let candidate = base_score_i + exons[ji].sum_score + js;
                if candidate > best_sum {
                    best_sum = candidate;
                    exons[i].link = Some(ji);
                    exons[i].sum_score = best_sum;
                }
            }
        }

        if best_sum >= highest_score {
            highest_score = best_sum;
            highest_idx = i;
        }
    }

    Some(highest_idx)
}

/// Traverse the trellis from the highest-scoring exon and assemble
/// gene predictions.  Returns a Vec of `EvmPrediction` objects.
pub fn traverse_path(exons: &[Exon], top_idx: usize) -> Vec<EvmPrediction> {
    let mut chain: Vec<usize> = Vec::new();
    let mut cur = Some(top_idx);
    while let Some(idx) = cur {
        chain.push(idx);
        cur = exons[idx].link;
    }
    chain.reverse(); // left to right

    // Split chain into individual gene predictions at terminal/single boundaries
    let mut predictions: Vec<EvmPrediction> = Vec::new();
    let mut current: Vec<usize> = Vec::new();

    for &idx in &chain {
        let exon = &exons[idx];
        let type_orient = exon.type_orient_key();
        if type_orient.contains("bound") {
            if !current.is_empty() {
                let pred = EvmPrediction::new(current.clone(), exons);
                predictions.push(pred);
                current.clear();
            }
            continue;
        }
        current.push(idx);
        // Terminate prediction at terminal+ or single, or initial- (reverse strand terminal)
        if type_orient == "terminal+" || type_orient.contains("single") || type_orient == "initial-" {
            let pred = EvmPrediction::new(current.clone(), exons);
            predictions.push(pred);
            current.clear();
        }
    }
    if !current.is_empty() {
        let pred = EvmPrediction::new(current, exons);
        predictions.push(pred);
    }

    predictions
}
