pub fn take<T: Copy>(data: &[T], indices: &[usize]) -> Vec<T> {
    let result = data
        .iter()
        .enumerate()
        .filter(|(c, _)| indices.contains(c))
        .map(|(_, &x)| x)
        .collect();

    result
}

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
