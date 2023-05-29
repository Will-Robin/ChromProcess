pub fn filter_on_height(data: &[f64], peak_pos: &[usize], height: f64) -> Vec<usize> {
    peak_pos
        .iter()
        .enumerate()
        .filter(|(_, &x)| data[x] > height)
        .map(|(c, _)| c)
        .collect()
}
