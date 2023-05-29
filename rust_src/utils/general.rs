pub fn take<T: Copy>(data: &[T], indices: &[usize]) -> Vec<T> {
    let result = data
        .iter()
        .enumerate()
        .filter(|(c, _)| indices.contains(c))
        .map(|(_, &x)| x)
        .collect();

    result
}

