use super::take;

pub fn agglomerative_1d_cluster(data: &[f64], bound: f64) -> Vec<Vec<f64>> {
    let mut output: Vec<Vec<f64>> = Vec::new();

    let mut current_cluster: Vec<f64> = Vec::new();

    for &x in data.iter() {
        let average = current_cluster.iter().sum::<f64>() / current_cluster.len() as f64;

        if !current_cluster.is_empty() && (x - average).abs() > bound {
            output.push(current_cluster.clone());
            current_cluster.clear();
        }

        current_cluster.push(x);
    }

    output
}

pub fn agglomerative_1d_cluster_indices(data: &[f64], bound: f64) -> Vec<Vec<usize>> {
    let mut output: Vec<Vec<usize>> = Vec::new();

    let mut current_cluster: Vec<usize> = Vec::new();

    for (c, &x) in data.iter().enumerate() {
        let average =
            take(data, &current_cluster).iter().sum::<f64>() / current_cluster.len() as f64;

        if !current_cluster.is_empty() && (x - average).abs() > bound {
            output.push(current_cluster.clone());
            current_cluster.clear();
        }

        current_cluster.push(c);
    }

    output
}
