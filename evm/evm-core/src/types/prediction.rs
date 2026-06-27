//! EVM prediction objects.

use crate::types::exon::Exon;

/// Run mode for a prediction search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PredMode {
    Standard,
    Intron,
}

impl PredMode {
    pub fn as_str(&self) -> &'static str {
        match self { PredMode::Standard => "STANDARD", PredMode::Intron => "INTRON" }
    }
}

/// A complete or partial consensus gene prediction produced by EVM.
#[derive(Debug, Clone)]
pub struct EvmPrediction {
    /// Indices into the exon pool for the exons making up this prediction,
    /// ordered 5′ → 3′ along the genomic sequence.
    pub exon_indices: Vec<usize>,
    /// True if the prediction was eliminated by the low-support filter.
    pub is_eliminated: bool,
    /// Mode under which the prediction was generated.
    pub mode: PredMode,
    /// Leftmost coordinate of the prediction span (1-based).
    pub lend: u32,
    /// Rightmost coordinate of the prediction span (1-based).
    pub rend: u32,
    /// Total evidence-weighted score of this prediction.
    pub total_score: f64,
}

impl EvmPrediction {
    pub fn new(exon_indices: Vec<usize>, exons: &[Exon]) -> Self {
        let (lend, rend) = span_of_indices(&exon_indices, exons);
        EvmPrediction {
            exon_indices,
            is_eliminated: false,
            mode: PredMode::Standard,
            lend,
            rend,
            total_score: 0.0,
        }
    }

    pub fn is_eliminated(&self) -> bool {
        self.is_eliminated
    }

    pub fn get_span(&self) -> (u32, u32) {
        (self.lend, self.rend)
    }
}

fn span_of_indices(indices: &[usize], exons: &[Exon]) -> (u32, u32) {
    let mut lend = u32::MAX;
    let mut rend = 0u32;
    for &idx in indices {
        let (el, er) = exons[idx].coords_sorted();
        if el < lend { lend = el; }
        if er > rend { rend = er; }
    }
    (lend, rend)
}

/// Intermediate prediction structure used during partition recombination.
#[derive(Debug, Clone)]
pub struct PartitionPred {
    pub lend: u32,
    pub rend: u32,
    pub class: PredClass,
    /// Serialised EVM output text for this prediction.
    pub text: String,
    /// Span length used for DP scoring.
    pub length: u32,
    /// Cumulative path score.
    pub path_score: u32,
    /// Index of predecessor prediction in the sorted array (None = no link).
    pub prev_link: Option<usize>,
    /// Predictions nested within introns of this one.
    pub intronic_preds: Vec<PartitionPred>,
    /// Set to true if this pred is encapsulated inside another pred's intron.
    pub encaps: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PredClass {
    Complete,
    Partial,
}
