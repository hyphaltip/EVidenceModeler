//! Intron scoring and intron-vector population.

use std::collections::HashMap;
use crate::types::genome::{FeatureVec, MaskVec, FEAT_DONOR, FEAT_ACCEPTOR};
use crate::types::evidence::EvClass;

/// Key for an intron: "end5_end3" stored in forward genomic coordinates.
pub type IntronKey = String;

/// Score accumulated for each intron.
pub type IntronScoreMap = HashMap<IntronKey, f64>;

/// Evidence list for each intron: Vec<(accession, ev_type)>.
pub type IntronEvidenceMap = HashMap<IntronKey, Vec<(String, String)>>;

/// Introns contributed only by ab-initio predictors.
pub type PredictedIntronMap = HashMap<IntronKey, f64>;

/// Per-base intron score vectors.
pub type IntronVec = Vec<f64>;

/// Add introns implied by a set of sorted exon coordinate pairs for one alignment chain.
///
/// `coords_list` should already be sorted ascending by end5.
/// Intron coordinates are determined from consecutive exon boundaries.
pub fn add_introns(
    accession: &str,
    coords_list: &[(u32, u32)], // (end5, end3) in forward-strand reference coords
    genomic_strand: char,
    weight: f64,
    intron_type: &str,   // ev_type string
    intron_ev_class: &EvClass,
    min_intron_length: u32,
    genome_features: &FeatureVec,
    mask: &MaskVec,
    introns_to_score: &mut IntronScoreMap,
    introns_to_evidence: &mut IntronEvidenceMap,
    predicted_introns: &mut PredictedIntronMap,
    genomic_seq_len: usize,
) {
    // Sort by first coordinate
    let mut sorted: Vec<(u32, u32)> = coords_list.to_vec();
    sorted.sort_by_key(|&(a, _)| a);

    for pair in sorted.windows(2) {
        let (_, first_end3) = pair[0];
        let (next_end5, _) = pair[1];

        if next_end5 < first_end3 {
            log::warn!(
                "ERROR adding intron for {}: next_end5 {} < first_end3 {}",
                accession, next_end5, first_end3
            );
            continue;
        }

        let potential_donor = first_end3 + 1;
        let potential_acceptor = next_end5 - 2;

        // Intron length check
        if potential_acceptor < potential_donor {
            continue;
        }
        let intron_length = potential_acceptor - potential_donor + 1;
        if intron_length < min_intron_length {
            log::warn!("Intron length ({}) < min ({})", intron_length, min_intron_length);
            continue;
        }

        // Require canonical donor and acceptor splice sites
        if genome_features.get(potential_donor as usize) != FEAT_DONOR
            || genome_features.get(potential_acceptor as usize) != FEAT_ACCEPTOR
        {
            continue;
        }

        // Score the intron region
        let mut intron_score = 0.0f64;
        for i in potential_donor..=potential_acceptor + 1 {
            if !mask.get(i as usize) {
                intron_score += weight;
            }
        }

        // Store in genomic reference coordinates
        let (intron_end5, intron_end3) = if genomic_strand == '-' {
            // Reverse-complement the coordinates back to forward strand
            let rc_end5 = genomic_seq_len as u32 - potential_donor + 1;
            let rc_end3 = genomic_seq_len as u32 - potential_acceptor + 1;
            (rc_end5, rc_end3)
        } else {
            (potential_donor, potential_acceptor)
        };

        let key = format!("{}_{}", intron_end5, intron_end3);
        *introns_to_score.entry(key.clone()).or_insert(0.0) += intron_score;
        introns_to_evidence
            .entry(key.clone())
            .or_default()
            .push((accession.to_string(), intron_type.to_string()));

        if intron_ev_class.is_abinitio() {
            *predicted_introns.entry(key).or_insert(0.0) += intron_score;
        }
    }
}

/// Parse intron key back to (end5, end3) coordinates.
pub fn intron_key_to_span(key: &str) -> Option<(u32, u32)> {
    let parts: Vec<&str> = key.split('_').collect();
    if parts.len() != 2 { return None; }
    let end5: u32 = parts[0].parse().ok()?;
    let end3: u32 = parts[1].parse().ok()?;
    Some((end5, end3))
}

/// Determine strand of an intron from its key (end5 < end3 → '+').
pub fn intron_key_strand(key: &str) -> char {
    if let Some((e5, e3)) = intron_key_to_span(key) {
        if e5 < e3 { '+' } else { '-' }
    } else {
        '+'
    }
}

/// Build per-base intron score vectors from the predicted intron map.
pub fn populate_intron_vectors(
    predicted_introns: &PredictedIntronMap,
    mask: &MaskVec,
    seq_len: usize,
) -> (IntronVec, IntronVec) {
    let mut fwd_vec = vec![0.0f64; seq_len + 2];
    let mut rev_vec = vec![0.0f64; seq_len + 2];

    for (key, &score) in predicted_introns {
        let (end5, end3) = match intron_key_to_span(key) {
            Some(v) => v,
            None => continue,
        };
        let strand = if end5 < end3 { '+' } else { '-' };
        let (lend, rend) = if end5 < end3 { (end5, end3) } else { (end3, end5) };

        // Compute adjusted length excluding masked positions
        let adj_len: f64 = (lend..=rend)
            .filter(|&i| !mask.get(i as usize))
            .count() as f64;
        if adj_len <= 0.0 { continue; }

        let score_per_bp = score / adj_len;
        let vec = if strand == '+' { &mut fwd_vec } else { &mut rev_vec };
        for i in lend..=rend {
            if !mask.get(i as usize) {
                vec[i as usize] += score_per_bp;
            }
        }
    }

    (fwd_vec, rev_vec)
}
