use pyo3::prelude::*;
mod find_peaks;
mod integrate;
mod utils;

#[pyfunction]
#[pyo3(name = "find_peaks")]
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

#[pyfunction]
fn cluster(data: Vec<f64>, bound: f64) -> PyResult<Vec<Vec<f64>>> {
    let result = utils::agglomerative_1d_cluster(&data, bound);
    Ok(result)
}

#[pyfunction]
fn cluster_indices(data: Vec<f64>, bound: f64) -> PyResult<Vec<Vec<usize>>> {
    let result = utils::agglomerative_1d_cluster_indices(&data, bound);
    Ok(result)
}

#[pymodule]
fn chromate(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_peaks_wrapper, m)?)?;
    m.add_function(wrap_pyfunction!(trapz, m)?)?;
    m.add_function(wrap_pyfunction!(cluster, m)?)?;
    m.add_function(wrap_pyfunction!(cluster_indices, m)?)?;
    Ok(())
}
