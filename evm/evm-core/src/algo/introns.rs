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

/// Perl `intron_key_to_intron_span`: parse the key and shift the acceptor base
/// to the exon-adjacent position — for '+' (end5 < end3) returns (end5, end3-1),
/// for '-' returns (end5, end3+1). This is the span Perl uses both when
/// populating the per-base predicted-intron vectors and in the filter offset
/// loop, and differs from the raw key span by one base at the acceptor end.
pub fn intron_key_to_intron_span(key: &str) -> Option<(u32, u32)> {
    let (end5, end3) = intron_key_to_span(key)?;
    if end5 < end3 {
        Some((end5, end3.saturating_sub(1))) // '+'
    } else {
        Some((end5, end3 + 1)) // '-'
    }
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
        // Perl distributes the score over the exon-adjacent intron span
        // (`intron_key_to_intron_span`), NOT the raw key span — this is one base
        // shorter at the acceptor end and is what the filter offset loop reads.
        let (end5, end3) = match intron_key_to_intron_span(key) {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intron_key_to_intron_span_shifts_acceptor() {
        // Perl `intron_key_to_intron_span`: '+' (end5 < end3) drops one base at
        // the acceptor (end3-1); '-' (end5 > end3) adds one (end3+1).
        assert_eq!(intron_key_to_intron_span("100_200"), Some((100, 199)));
        assert_eq!(intron_key_to_intron_span("200_100"), Some((200, 101)));
    }

    #[test]
    fn predicted_intron_vector_distributes_over_exon_adjacent_span() {
        // Regression for task #16 (filter offset residual): the per-base
        // predicted-intron vector must spread the score over the
        // `intron_key_to_intron_span` span (D..A-1 for '+'), NOT the raw key
        // span (D..A). Build a single forward intron "10_20" with score 50.
        // The exon-adjacent span is 10..=19 (10 bases) → 5.0 per base, and base
        // 20 (the raw acceptor) must receive nothing.
        let mut predicted: PredictedIntronMap = HashMap::new();
        predicted.insert("10_20".to_string(), 50.0);
        let mask = MaskVec::new(64);
        let (fwd, rev) = populate_intron_vectors(&predicted, &mask, 50);

        for i in 10..=19 {
            assert!((fwd[i] - 5.0).abs() < 1e-9, "base {i} expected 5.0, got {}", fwd[i]);
        }
        assert_eq!(fwd[20], 0.0, "raw acceptor base 20 must stay zero");
        assert_eq!(fwd[9], 0.0);
        // Total conserved == the intron score.
        let total: f64 = fwd.iter().sum();
        assert!((total - 50.0).abs() < 1e-9, "total {total} != 50.0");
        assert!(rev.iter().all(|&v| v == 0.0), "reverse vector must stay zero");
    }
}
