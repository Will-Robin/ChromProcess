use pyo3::prelude::*;
mod find_peaks;
mod integrate;
mod utils;

/// This function finds peak indices and the indices of their starts and ends, given a minimum
/// height.
#[pyfunction]
#[pyo3(name = "find_peaks", signature = (data, height))]
fn find_peaks_wrapper(
    data: Vec<f64>,
    height: f64,
) -> PyResult<(Vec<usize>, Vec<usize>, Vec<usize>)> {
    let result = find_peaks::find_peaks(&data, Some(height));
    Ok(result)
}

/// This function integrates a signal using the trapezium rule.
#[pyfunction]
#[pyo3(signature = (data_y, data_x))]
fn trapz(data_y: Vec<f64>, data_x: Vec<f64>) -> PyResult<f64> {
    let result = integrate::trapezium_rule_integrate(&data_y, Some(&data_x));
    Ok(result)
}

/// This function clusters points in an ordered 1D array.
#[pyfunction]
#[pyo3(signature = (data, bound))]
fn cluster(data: Vec<f64>, bound: f64) -> PyResult<Vec<Vec<f64>>> {
    let result = utils::agglomerative_1d_cluster(&data, bound);
    Ok(result)
}

/// This function clusters points in an ordered 1D array, returning the indices of the points in each
/// cluster.
#[pyfunction]
#[pyo3(signature = (data, bound))]
fn cluster_indices(data: Vec<f64>, bound: f64) -> PyResult<Vec<Vec<usize>>> {
    let result = utils::agglomerative_1d_cluster_indices(&data, bound);
    Ok(result)
}

/// chromate is a submodule written in Rust. Detailed source code documentation will be supplied
/// seprately.
#[pymodule]
#[pyo3(name="chromate")]
fn chromate(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_peaks_wrapper, m)?)?;
    m.add_function(wrap_pyfunction!(trapz, m)?)?;
    m.add_function(wrap_pyfunction!(cluster, m)?)?;
    m.add_function(wrap_pyfunction!(cluster_indices, m)?)?;
    Ok(())
}
