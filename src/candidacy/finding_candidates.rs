use std::fmt::Display;

#[derive(Copy, Clone, Default)]
struct Candidate {
    period: f64,
    snr: f64,
    duration: f64,
    phase: f64,
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

fn signal_to_noise_ratio(
    deepest_dip: f64, num_points: f64, sigma: &f64) -> Option<f64> {
    let snr = (- deepest_dip / sigma ) * num_points.sqrt();
    Some(snr)
}

fn search_bins<const N_BINS: usize>(
    trial_period: &f64, trial_duration: &f64,
    binned_flux: &Vec<f64>, counts: &[u32; N_BINS])
    -> Option<(f64, usize, f64, u32)> {

    let window_width: usize = (((trial_duration / 24.0) / trial_period) * N_BINS as f64).round() as usize;

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
    let centre_phase = centre_phase_numerator / N_BINS as f64; // used for SNR later
    Some((deepest_dip, window_width, centre_phase, best_count))
}

fn signal_to_noise_ratio(
    deepest_dip: f64, num_points: f64, sigma: &f64) -> Option<f64> {
    let snr = (- deepest_dip / sigma ) * num_points.sqrt();
    Some(snr)
}

fn find_max_candidate(c_snr: Vec<f64>, c_period: Vec<f64>) -> (f64, f64) {
    let mut max_snr = f64::NEG_INFINITY;
    let mut max_index: usize = 0;
    for (index, snr) in c_snr.iter().enumerate() {
        if *snr > max_snr {
            max_snr = *snr;
            max_index = index;
        }
    }
    (max_snr, c_period[max_index])
}

fn get_phases(times: &[f64], period: &f64) -> Vec<f64> {

    let mut phases = Vec::new();

    for t in times {
        let phase = (t % period) / period;
        phases.push(phase);
    }

    phases
}

fn signal_to_noise_ratio(deepest_dip: f64, num_points: f64, sigma: &f64) -> Option<f64> {
    let snr = (- deepest_dip / sigma ) * num_points.sqrt();
    Some(snr)
}

pub fn find_candidates<const N_BINS: usize> (f_all: &Vec<f64>, t_all: &Vec<f64>) -> Vec<Candidate> {
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