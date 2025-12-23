use std::fmt::Display;
use crate::data_access::{CleanLightCurve};

#[derive(Copy, Clone, Default)]
pub struct Candidate {
    period: f64,
    snr: f64,
    duration: f64,
    phase: f64,
}

impl Candidate {
    fn new(period: f64, snr: f64, duration: f64, phase: f64) -> Candidate {
        Candidate { period, snr, duration, phase }
    }

    fn snr(&self) -> f64 {
        self.snr
    }

    fn duration(&self) -> f64 {
        self.duration
    }

    fn phase(&self) -> f64 {
        self.phase
    }

    fn period(&self) -> f64 {
        self.period
    }
}

impl Display for Candidate {

    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Period = {} days, SNR = {}, Duration = {} hrs, Phase = {}",
            self.period, self.snr, self.duration, self.phase
        )
    }
}

fn signal_to_noise_ratio(deepest_dip: f64, num_points: usize, sigma: f64) -> Option<f64> {
    let snr = (- deepest_dip / sigma ) * (num_points as f64).sqrt();
    Some(snr)
}

fn binning(flux: &Vec<f64>, phases: &[f64], n_bins: usize) -> (Vec<f64>, Vec<u32>) {
    let mut count: Vec<u32> = vec![0u32; n_bins];
    let mut sum_flux: Vec<f64> = vec![0.0; n_bins];

    for (flux, phase) in flux.iter().zip(phases.iter()) {
        let bin_index = (phase * n_bins as f64).floor() as usize;
        count[bin_index] += 1;
        sum_flux[bin_index] += flux;
    }

    let mut binned_flux: Vec<f64> = Vec::with_capacity(n_bins);
    for i in 0..n_bins {
        if count[i] == 0 {
            binned_flux.push(0.0); // or some neutral value
        } else {
            binned_flux.push(sum_flux[i] / count[i] as f64);
        }
    }
    (binned_flux, count)
}

fn search_bins
(trial_period: f64, trial_duration: f64, binned_flux: &Vec<f64>, counts: &Vec<u32>,  n_bins: usize,)
    -> Option<(f64, usize, f64)> {

    let window_width: usize = (((trial_duration / 24.0) / trial_period) * n_bins as f64).round() as usize;

    if window_width == 0 || window_width > binned_flux.len() { return None; }

    let mut sum_flux: f64 = binned_flux[..window_width].iter().sum();
    let mut sum_count: u32 = counts[..window_width].iter().sum();

    if sum_count == 0 { return None;}

    let mut centre_phase_numerator: f64 = 0.0;
    let mut deepest_dip: f64 = sum_flux / window_width as f64;
    let mut best_count: u32 = sum_count;
    for right in window_width..binned_flux.len() {
        sum_flux += binned_flux[right] - binned_flux[right - window_width];
        /*
            sum_count += counts[right] - counts[right - window_width];
            will cause a panic, because u32 - u32 might be such that it is negative.
            since sum_count >= counts[right - window_with] always, it will never be negative.
        */
        sum_count += counts[right];
        sum_count -= counts[right - window_width];
        let avg = sum_flux / window_width as f64;
        if avg < deepest_dip {
            deepest_dip = avg;
            best_count = sum_count;
            centre_phase_numerator = (2.0 * right as f64 - window_width as f64) / 2.0;
        }
    }
    let centre_phase = centre_phase_numerator / n_bins as f64; // used for SNR later
    Some((deepest_dip, window_width, centre_phase))
}

fn get_phases(times: &Vec<f64>, period: f64) -> Vec<f64> {

    let mut phases = Vec::new();

    for t in times {
        let phase = (t % period) / period;
        phases.push(phase);
    }

    phases
}

pub fn trial_periods(t_size: usize) -> Vec<f64> {
    (2..=t_size).map(|i| i as f64 * 0.5).collect()
}

pub fn find_candidates
(lc: &CleanLightCurve, trial_periods: &Vec<f64>, n_bins: usize, n_dur: usize) -> Vec<Candidate> {
    let mut candidates: Vec<Candidate> = Vec::new();

    for &period in trial_periods.iter() {
        let phases: Vec<f64> = get_phases(lc.time(), period);
        let mut binned_flux: Vec<f64> = Vec::with_capacity(n_bins);
        let mut counts: Vec<u32> = Vec::with_capacity(n_bins);
        (binned_flux, counts) = binning(lc.flux(), &phases, n_bins);
        let trial_durations: Vec<f64> = (0..n_dur as i32)
            .map(|i| i as f64 * 0.5)
            .collect();
        let mut deepest_dip: f64 = f64::INFINITY;
        let mut window_length: usize = 0;
        let mut duration: f64 = 0.0;
        let mut centre_phase: f64 = 0.0;

        for &dur in trial_durations.iter() {
            if let Some((box_avg, window_len, cent_phase))
                = search_bins(period, dur, &binned_flux, &counts, n_bins)
            {
                if box_avg < deepest_dip {
                    deepest_dip = box_avg;
                    window_length = window_len;
                    duration = dur;
                    centre_phase = cent_phase;
                }
            }
        }

        if !deepest_dip.is_finite() {
            deepest_dip = 0.0;
        }

        let snr = signal_to_noise_ratio(deepest_dip, window_length, lc.sigma())
            .expect("signal_to_noise_ratio failed");

        candidates.push(Candidate::new(period, snr, duration, centre_phase));
    }

    candidates
}

pub fn transit_estimate(candidates: &Vec<Candidate>) -> Candidate {
    let mut max_snr = f64::NEG_INFINITY;
    let mut max_index: usize = 0;
    for (index, cand) in candidates.iter().enumerate() {
        let check_snr = cand.snr();
        if check_snr > max_snr {
            max_index = index;
            max_snr = check_snr;
        }
    }

    candidates[max_index].clone()
}

/*
Temporary find_candidates! I plan to design a better one soon.

pub fn find_candidates<const N_BINS: usize>
(f_all: &Vec<f64>, t_all: &Vec<f64>, trial_periods: &Vec<f64>, sigma_all: f64) -> Vec<Candidate> {
    let mut candidates: Vec<Candidate> = Vec::new();
    const N_BINS: usize = 200;
    for period in trial_periods.iter() {
        let phase: Vec<f64> = get_phases(&t_all, period);
        let mut binned_flux: Vec<f64> = Vec::with_capacity(N_BINS);
        let counts: [u32; N_BINS];
        (binned_flux, counts) = binning::<N_BINS>(&f_all, &phase);
        assert_eq!(binned_flux.len(), N_BINS);

        let trial_durations: Vec<f64> = (0..20)
            .map(|i| i as f64 * 0.5)
            .collect();
        let mut deepest_dip: f64 = f64::INFINITY;
        let mut window_length: usize = 0;
        let mut duration: f64 = 0.0;
        let mut centre_phase: f64 = 0.0;

        for d in trial_durations.iter() {
            let boxed = search_bins::<N_BINS>(period, d, &binned_flux, &counts);
            if boxed.is_some() {
                let boxed = boxed.unwrap();
                let box_average = boxed.0;
                if box_average < deepest_dip {
                    deepest_dip = box_average;
                    window_length = boxed.1;
                    duration = *d;
                    centre_phase = boxed.2;
                }
            }
        }
        if !deepest_dip.is_finite() {
            deepest_dip = 0.0;
        }

        // println!("Found deepest dip. Now onto scoring with SNR...");
        let window_length: f64 = window_length as f64;
        let snr = signal_to_noise_ratio(deepest_dip, window_length, &sigma_all)
            .expect("SNR failure");
        let potential: Candidate = Candidate {
            period: *period,
            snr,
            duration,
            phase: centre_phase,
        };
        candidates.push(potential);
    }
    candidates
}

 */

