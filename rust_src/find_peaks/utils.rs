use super::super::utils;
use super::filters;

pub fn find_peak_indices(data: &[f64]) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    // Find transition points
    let transition_points: Vec<usize> = data
        .iter()
        .zip(data.iter().skip(1))
        .zip(data.iter().skip(2))
        .enumerate()
        .filter(|(p, ((&a, &b), &c))| {
            if *p == 0 || *p == data.len() - 1 {
                // first and last values are defined as transition points
                return true;
            }

            let sign_1 = (b - a).signum();
            let sign_2 = (c - b).signum();

            if sign_1 != sign_2 {
                return true;
            } else if (b - a) == 0.0 || (c - b) == 0.0 {
                return true;
            }
            false
        })
        .map(|(c, _)| c + 1)
        .collect();

    // Find peak positions - transition points where the gradient decreases
    let peak_pos: Vec<usize> = transition_points
        .iter()
        .filter(|&&x| data[x + 1] - data[x] < 0.0)
        .copied()
        .collect();

    // Select transition points to the left and right of each peak as bounds
    let left_bounds: Vec<usize> = transition_points
        .iter()
        .zip(transition_points.iter().skip(1))
        .filter(|(_, &x)| data[x + 1] - data[x] < 0.0)
        .map(|(&x, _)| x)
        .collect();

    let right_bounds: Vec<usize> = transition_points
        .iter()
        .zip(transition_points.iter().skip(1))
        .filter(|(&x, _)| data[x + 1] - data[x] < 0.0)
        .map(|(_, &x)| x)
        .collect();

    (peak_pos, left_bounds, right_bounds)
}

pub fn find_peaks(data: &[f64], height: Option<f64>) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    // Find peaks indices
    let (mut peaks, mut left_edges, mut right_edges) = find_peak_indices(data);

    // Filter down the picked peaks based on height
    if let Some(h) = height {
        let retain_idx = filters::filter_on_height(data, &peaks, h);

        peaks = utils::take(&peaks, &retain_idx);
        left_edges = utils::take(&left_edges, &retain_idx);
        right_edges = utils::take(&right_edges, &retain_idx);
    }

    (peaks, left_edges, right_edges)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_find_peak_indices() {
        let data: Vec<f64> = vec![
            0.0, 0.0, 0.1, 0.2, 0.3, 0.2, 0.3, 0.2, 0.1, 0.0, 0.0, 0.01, 0.02, 0.03, 0.01, 0.02,
            0.0, 0.1, 0.15, 0.175, 0.2, 0.3, 0.2, 0.175, 0.15, 0.1, 0.0,
        ];

        let (indices, left_bounds, right_bounds) = find_peak_indices(&data);

        assert_eq!(indices, vec![4, 6, 13, 15, 21]);
        assert_eq!(left_bounds, vec![1, 5, 10, 14, 16]);
        assert_eq!(right_bounds, vec![5, 9, 14, 16]);
    }

    #[test]
    fn test_find_peaks() {
        let data: Vec<f64> = vec![
            0.0, 0.0, 0.1, 0.2, 0.3, 0.2, 0.3, 0.2, 0.1, 0.0, 0.0, 0.01, 0.02, 0.03, 0.01, 0.02,
            0.0, 0.1, 0.15, 0.175, 0.2, 0.3, 0.2, 0.175, 0.15, 0.1, 0.0,
        ];

        let (indices, left_bounds, right_bounds) = find_peaks(&data, Some(0.2));

        assert_eq!(indices, vec![4, 6, 21]);
        assert_eq!(left_bounds, vec![1, 5, 16]);
        assert_eq!(right_bounds, vec![5, 9, 16]);
    }
}
