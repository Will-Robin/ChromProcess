pub struct Chromatogram {
    time: Vec<f64>,
    signal: Vec<f64>,
}

impl Chromatogram {
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Chromatogram {
        Chromatogram { time: x, signal: y }
    }
}
