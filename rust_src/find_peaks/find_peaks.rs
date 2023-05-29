use super::filters;
use super::utils;

pub fn find_peaks(data: &[f64], height: Option<f64>) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    // Find peaks indices
    let (mut peaks, mut left_edges, mut right_edges) = utils::find_peak_indices(data);

    // Filter down the picked peaks based on criteria
    if let Some(h) = height {
        let retain_idx = filters::filter_on_height(data, &peaks, h);

        peaks = utils::take(&peaks, &retain_idx);
        left_edges = utils::take(&left_edges, &retain_idx);
        right_edges = utils::take(&right_edges, &retain_idx);
    }

    (peaks, left_edges, right_edges)
}
