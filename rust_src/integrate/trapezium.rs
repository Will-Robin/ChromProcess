pub fn trapezium_on_points(data_y: &[f64], data_x: &[f64]) -> f64 {
    data_y
        .iter()
        .zip(data_y.iter().skip(1))
        .zip(data_x.windows(2))
        .map(|((&a, &b), c)| ((a + b) / 2.0) * (c[1] - c[0]))
        .sum::<f64>()
}

pub fn trapezium_rule_integrate(data_y: &[f64], data_x: Option<&[f64]>) -> f64 {
    match data_x {
        Some(d) => trapezium_on_points(data_y, d),
        None => trapezium_on_points(data_y, &vec![1.0; data_y.len()]),
    }
}
