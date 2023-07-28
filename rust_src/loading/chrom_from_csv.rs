use super::super::data_structs::Chromatogram;
use std::fs;

pub fn chrom_from_csv(filename: String) -> Chromatogram {
    let mut time: Vec<f64> = Vec::new();
    let mut signal: Vec<f64> = Vec::new();

    for (c, line) in fs::read_to_string(filename).unwrap().lines().enumerate() {
        if c == 0 {
            let header = line.split(",");
        } else {
            time.push(line.to_string())
        }
    }

    Chromatogram::new(time, signal);
}
