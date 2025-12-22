pub mod fits_access;
pub mod data_cleaner;

use fitsio::FitsFile;

pub struct CleanLightCurve {
    time: Vec<f64>,
    flux: Vec<f64>,
    sigma: f64
}

impl CleanLightCurve {
    pub fn new(time: Vec<f64>, flux: Vec<f64>, sigma: f64) -> CleanLightCurve {
        CleanLightCurve { time, flux, sigma }
    }

    pub fn time(&self) -> &Vec<f64> {
        &self.time
    }

    pub fn flux(&self) -> &Vec<f64> {
        &self.flux
    }

    pub fn sigma(&self) -> f64 {
        self.sigma
    }
}

