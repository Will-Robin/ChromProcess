pub fn filter_on_height(data: &[f64], peak_pos: &[usize], height: f64) -> Vec<usize> {
    peak_pos
        .iter()
        .enumerate()
        .filter(|(_, &x)| data[x] > height)
        .map(|(c, _)| c)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_on_height() {
        let test_data: Vec<f64> = vec![1.0, 2.0, 1.0, 3.0, 0.0];
        let peak_pos: Vec<usize> = vec![1,3];
        let result = filter_on_height(&test_data, &peak_pos, 2.0);
        assert_eq!(result, vec![1]);
    }
}
