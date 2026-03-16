use std::fmt::Display;
use std::marker::PhantomData;

#[derive(Clone)]
pub struct ExoplanetMarker;
#[derive(Clone)]
pub struct BinarySystemMarker;
#[derive(Clone)]
pub struct CandidateMarker;

#[derive(Copy, Clone, Default)]
pub struct OrbitalSystem<T> {
    period: f64,
    snr: f64,
    duration: f64,
    phase: f64,
    _marker: PhantomData<T>,
}

impl<T> OrbitalSystem<T> {
    pub fn new(period: f64, snr: f64, duration: f64, phase: f64) -> Self {
        Self {
            period,
            snr,
            duration,
            phase,
            _marker: PhantomData,
        }
    }

    pub fn snr(&self) -> f64 {
        self.snr
    }
    pub fn period(&self) -> f64 {
        self.period
    }
    pub fn duration(&self) -> f64 {
        self.duration
    }
    pub fn phase(&self) -> f64 {
        self.phase
    }
}

impl<T> Display for OrbitalSystem<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Period = {} days, SNR = {}, Duration = {} hrs, Phase = {}",
            self.period, self.snr, self.duration, self.phase
        )
    }
}

// Aliases for convenience
pub type Exoplanet = OrbitalSystem<ExoplanetMarker>;
pub type BinarySystem = OrbitalSystem<BinarySystemMarker>;
pub type Candidate = OrbitalSystem<CandidateMarker>;
