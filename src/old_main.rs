/*
    Cameron Emmanuel
    Nov - 05 - 2025
 */

use fitsio::FitsFile;
use std::path::PathBuf;
use std::io::Write;
use std::process::{Command, Stdio};
use serde::Serialize;

fn plot_to_python(
    x_label: &str,
    y_label: &str,
    dat_x: &[f64],
    dat_y: &[f64],
) -> std::io::Result<()> {
    #[derive(Serialize)]
    struct Data {
        x_label: String,
        y_label: String,
        x: Vec<f64>,
        y: Vec<f64>,
    }

    let payload = Data {
        x_label: x_label.to_owned(),   // or x_label.to_string()
        y_label: y_label.to_owned(),
        x: dat_x.to_vec(),
        y: dat_y.to_vec(),
    };

    let json = serde_json::to_string(&payload).expect("failed to serialize");

    let mut child = Command::new("python3")
        .arg("plotting_data.py")
        .stdin(Stdio::piped())
        .spawn()
        .expect("failed to start python");

    {
        let stdin = child.stdin.as_mut().expect("failed to open stdin");
        stdin.write_all(json.as_bytes()).expect("failed to write json");
    }

    Ok(())
}


fn load_tess_data(f: &mut FitsFile) -> fitsio::errors::Result<(Vec<f64>, Vec<f64>, Vec<i32>)> {
    let hdu = f.hdu(1)?;
    let time: Vec<f64> = hdu.read_col(f, "TIME")?;
    let flux: Vec<f64> = hdu.read_col(f, "PDCSAP_FLUX")?;
    let qual: Vec<i32> = hdu.read_col(f, "QUALITY")?;
    Ok((time, flux, qual))
}

fn filter_and_sort_tess_data(
    flux: Vec<f64>, time: Vec<f64>, qual: Vec<i32>) -> fitsio::errors::Result<(Vec<f64>, Vec<f64>)> {
    let good_indices: Vec<usize> = (0..time.len())
        .filter(|&i| qual[i] == 0 && time[i].is_finite() && flux[i].is_finite())
        .collect();

    let mut paired: Vec<(f64, f64)> = good_indices
        .iter()
        .map(|&i| (time[i], flux[i]))
        .collect();

    paired.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let (t_sorted, f_sorted): (Vec<f64>, Vec<f64>) = paired.into_iter().unzip();

    Ok((t_sorted, f_sorted))
}

pub fn median(data: &[f64]) -> Option<f64> {
    if data.is_empty() {
        return None;
    }
    // Filter out NaNs so they don't break ordering
    let mut v: Vec<f64> = data.iter().copied().filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = v.len() / 2;

    if v.len() % 2 == 0 {
        // even length → average of the two middle values
        Some((v[mid - 1] + v[mid]) * 0.5)
    } else {
        // odd length → middle value
        Some(v[mid])
    }
}

fn median_cadence(times: &[f64]) -> Option<f64> {
    if times.len() < 2 { return None; }
    let mut dts: Vec<f64> = times3
        .windows(2)
        .map(|w| w[1] - w[0])
        .filter(|dt| dt.is_finite() && *dt > 0.0)
        .collect();
    if dts.is_empty() { return None; }
    median(&mut dts)
}

fn segment_on_gaps(t: &[f64], dt_med: f64, gap_factor: f64) -> Vec<(usize, usize)> {
    let mut segment_bounds: Vec<(usize, usize)> = Vec::new();
    let threshold = dt_med * gap_factor;

    if t.is_empty() { return segment_bounds; }

    let mut start = 0;
    for i in 1..t.len() {
        if (t[i] - t[i-1]) > threshold {
            segment_bounds.push((start, i));
            start = i;
        }
    }

    segment_bounds.push((start, t.len()));
    segment_bounds
}

struct NormSegment {
    f_norm: Vec<f64>,
    med: f64,
    sigma: f64,
}
impl NormSegment {
    fn segment_noise_clean(&mut self, k: f64){
        let threshold = k * self.sigma;

        for f in self.f_norm.iter_mut() {
            if !f.is_finite() || f.abs() > threshold {
                *f = f64::NAN;
            }
        }
    }
}

fn median_absolute_deviation(arr: &[f64]) -> Option<f64> {
    if arr.is_empty() { return None; }

    let med = median(arr)?;
    let mut deviations: Vec<f64> = arr
        .iter()
        .copied()
        .filter(|f| f.is_finite())
        .map(|f| (f - med).abs())
        .collect();

    if deviations.is_empty() { return None; }
    median(&mut deviations)
}


fn normalize_segments(f: &[f64], segments: &[(usize, usize)]) -> Vec<NormSegment> {
    let mut out = Vec::with_capacity(segments.len());

    for &(s, e) in segments {
        let f_seg = &f[s..e]; // e is exclusive
        let med = median(&f_seg).expect("empty segment");
        assert!(med.is_finite() && med != 0.0, "segment median zero or infinite; normalization blow up");
        let inv = 1.0 / med;
        let f_norm: Vec<f64> = f_seg.iter().map(|&f| inv * f - 1.0).collect();
        let mad = median_absolute_deviation(&f_norm).expect("couldn't find MAD");
        let sigma = 1.4826 * mad;

        out.push(NormSegment {
            f_norm,
            med,
            sigma,
        });
    }
    out
}

fn get_phases(times: &[f64], period: &f64) -> Vec<f64> {

    let mut phases = Vec::new();

    for t in times {
        let phase = (t % period) / period;
        phases.push(phase);
    }

    phases
}

fn binning<const N_BINS: usize>(flux: &[f64], phases: &[f64]) -> (Vec<f64>, [u32; N_BINS]) {
    let mut count: [u32; N_BINS] = [0; N_BINS];
    let mut sum_flux: [f64; N_BINS] = [0.0; N_BINS];

    for (f, p) in flux.iter().zip(phases.iter()) {
        let bin_index = (p * N_BINS as f64).floor() as usize;
        count[bin_index] += 1;
        sum_flux[bin_index] += f;
    }

    let mut binned_flux = Vec::with_capacity(N_BINS);
    for i in 0..N_BINS {
        if count[i] == 0 {
            binned_flux.push(0.0); // or some neutral value
        } else {
            binned_flux.push(sum_flux[i] / count[i] as f64);
        }
    }
    (binned_flux, count)
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

fn main() -> fitsio::errors::Result<()> {
    let fits_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join("test")
        .join("tess2025180145000-s0094-0000000006286534-0291-s_lc.fits");

    let mut fptr = FitsFile::open(fits_path)?;
    // fptr.pretty_print()?;

    let (t, f, q) = load_tess_data(&mut fptr)?;
    assert_eq!(t.len(), f.len(), "time/flux length mismatch");
    assert_eq!(t.len(), q.len(), "time/quality length mismatch");

    let (t, f) = filter_and_sort_tess_data(f, t, q)?;
    assert_eq!(t.len(), f.len(), "time/flux length mismatch");
    assert!(t.iter().all(|&val| !val.is_nan() && val.is_finite()), "NaN or infinite in t");
    assert!(f.iter().all(|&val| !val.is_nan()), "NaN detected in f");
    assert!(t.windows(2).all(|w| w[1] > w[0]), "t is not chronological");

    let dt_med = median_cadence(&t).expect("median cadence not found");
    let gap_factor = 5.0;
    let seg_bounds = segment_on_gaps(&t, dt_med, gap_factor); // dt_med loses ownership

    let mut norm_segments: Vec<NormSegment> = normalize_segments(&f, &seg_bounds);

    // cleaning all segments
    const K: f64 = 6.0;
    for seg in norm_segments.iter_mut() {
        seg.segment_noise_clean(K);
    }

    // stitch the segments back together
    let f_all: Vec<f64> = norm_segments
        .iter()
        .flat_map(|seg| seg.f_norm.iter().copied())
        .collect();

    let t_all = t; // for clarity

    let mad_all = median_absolute_deviation(&f_all)
        .expect("median cadence of f_all not found");
    let sigma_all = 1.4826 * mad_all;

    let trial_periods: Vec<f64> = (1..=250)
        .map(|i| i as f64 * 0.5)
        .collect();
    // println!("{:#?}", trial_periods);

    // goal is to get a (period, SNR)
    #[derive(Copy, Clone, Default)]
    struct Candidate {
        period: f64,
        snr: f64,
        duration: f64,
        phase: f64,
    }

    impl std::fmt::Display for Candidate {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(
                f,
                "Period = {} days, SNR = {}, Duration = {} hrs, Phase = {}",
                self.period, self.snr, self.duration, self.phase
            )
        }
    }

    println!("Testing Periods");
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

    println!("Plotting Potential candidates: SNR vs. Period");
    let candidate_periods: Vec<f64> = candidates.iter().map(|c| c.period).collect();
    let candidate_snr: Vec<f64> = candidates.iter().map(|c| c.snr).collect();
    plot_to_python("Period".into(), "SNR".into(), &candidate_periods, &candidate_snr);

    let (max_snr, max_period) = find_max_candidate(candidate_snr, candidate_periods);

    let mut transit: Candidate = Candidate::default();
    for candidate in candidates.iter() {
        if candidate.period == max_period {
            transit = *candidate;
        }
    }

    println!("Found the likely Transit:\n {}", transit);

    Ok(())
}


