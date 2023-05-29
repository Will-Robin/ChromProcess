use pyo3::prelude::*;
mod find_peaks;

#[pyfunction]
fn find_peaks_wrapper(
    data: Vec<f64>,
    height: f64,
) -> PyResult<(Vec<usize>, Vec<usize>, Vec<usize>)> {
    let result = find_peaks::find_peaks(&data, Some(height));
    Ok(result)
}

#[pymodule]
fn chromate(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_peaks_wrapper, m)?)?;
    Ok(())
}
