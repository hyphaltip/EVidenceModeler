//! Sliding-window peak detection for gene boundary analysis.

/// Analyse a positional score vector using a sliding window and return
/// peaks where the windowed score exceeds `threshold`.
///
/// Mimics the Perl `analyze_peaks()` subroutine:
/// - Slide a window of `window_size` positions.
/// - Track the best score/position within each non-overlapping window span.
/// - Emit a peak when `best_score > threshold` and the window has moved past.
pub fn analyze_peaks(
    vector: &[f64],     // 1-indexed; index 0 is unused
    seq_len: usize,
    window_size: usize,
    threshold: f64,
) -> Vec<(u32, f64)> {
    let mut found_peaks: Vec<(u32, f64)> = Vec::new();

    let mut best_score_so_far = 0.0f64;
    let mut best_pos_so_far: usize = 0;
    let mut current_peak_score = 0.0f64;

    for i in 1..=seq_len {
        let leading_edge = i;
        let trailing_edge = if i >= window_size { i - window_size } else { 0 };

        // Check if we've moved beyond the window
        if leading_edge.saturating_sub(best_pos_so_far) > window_size {
            if best_score_so_far > threshold {
                found_peaks.push((best_pos_so_far as u32, best_score_so_far));
            }
            // Reset
            best_score_so_far = 0.0;
            best_pos_so_far = i;
        }

        let val = vector.get(leading_edge).copied().unwrap_or(0.0);
        current_peak_score += val;

        if trailing_edge > 0 {
            current_peak_score -= vector.get(trailing_edge).copied().unwrap_or(0.0);
        }

        if current_peak_score > best_score_so_far {
            best_score_so_far = current_peak_score;
            best_pos_so_far = leading_edge;
        }
    }

    // Emit final window
    if best_score_so_far > threshold {
        found_peaks.push((best_pos_so_far as u32, best_score_so_far));
    }

    found_peaks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_peak() {
        // Build a vector with a clear spike at position 5
        let mut v = vec![0.0f64; 11]; // positions 0..=10
        v[5] = 10.0;
        let peaks = analyze_peaks(&v, 10, 3, 2.0);
        assert!(!peaks.is_empty());
        assert_eq!(peaks[0].0, 5);
    }
}
