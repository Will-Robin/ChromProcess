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

#[pyfunction]
fn smooth(data: Vec<f64>, window_size: usize) -> PyResult<Vec<f64>> {
    let result = utils::adjacent_average(&data, window_size);
    Ok(result)
}

#[pymodule]
fn pyk_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_peaks_wrapper, m)?)?;
    m.add_function(wrap_pyfunction!(smooth, m)?)?;
    Ok(())
}
