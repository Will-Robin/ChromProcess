use pyo3::prelude::*;
mod find_peaks;
mod integrate;

#[pyfunction]
#[pyo3(name="find_peaks")]
fn find_peaks_wrapper(
    data: Vec<f64>,
    height: f64,
) -> PyResult<(Vec<usize>, Vec<usize>, Vec<usize>)> {
    let result = find_peaks::find_peaks(&data, Some(height));
    Ok(result)
}

#[pyfunction]
fn trapz(data_y: Vec<f64>, data_x: Vec<f64>) -> PyResult<f64> {
    let result = integrate::trapezium_rule_integrate(&data_y, Some(&data_x));
    Ok(result)
}

#[pymodule]
fn chromate(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_peaks_wrapper, m)?)?;
    m.add_function(wrap_pyfunction!(trapz, m)?)?;
    Ok(())
}
