pub mod fits_access;
pub mod data_cleaner;

use fitsio::FitsFile;

struct CleanLightCurves {
    time: Vec<f64>,
    flux: Vec<f64>,
    sigma: f64
}

impl CleanLightCurves {
    pub fn new(time: Vec<f64>, flux: Vec<f64>, sigma: f64) -> CleanLightCurves {
        CleanLightCurves {
            time,
            flux,
            sigma
        }
    }
}
