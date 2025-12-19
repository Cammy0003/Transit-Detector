pub mod fits_access;
mod data_cleaner;

struct CleanLightCurve {
    time: Vec<f64>,
    flux: Vec<f64>,
    qual: Vec<i32>
}

impl CleanLightCurve {
    pub fn new(time: Vec<f64>, flux: Vec<f64>, qual: Vec<i32>) -> CleanLightCurve {
        CleanLightCurve {
            time,
            flux,
            qual
        }
    }



}