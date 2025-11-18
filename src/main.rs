/*
    Cameron Emmanuel
    Nov - 05 - 2025
 */

use fitsio::FitsFile;
use fitsio::tables::Column;
use std::path::PathBuf;
use plotters::prelude::*;
use plotters::series::LineSeries;
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

fn filter_and_sort_tess_data(flux: Vec<f64>, time: Vec<f64>, qual: Vec<i32>) -> fitsio::errors::Result<(Vec<f64>, Vec<f64>)> {
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

fn median(v: &mut Vec<f64>) -> Option<f64> {
    if v.is_empty() { return None; }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = v.len() / 2;
    if v.len() % 2 == 0 {
        Some((v[mid - 1] + v[mid]) / 2.0)
    } else {
        Some(v[mid])
    }
}
fn median_cadence(times: &[f64]) -> Option<f64> {
    if times.len() < 2 { return None; }
    let mut dts: Vec<f64> = times
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
    return segment_bounds;
}

struct NormSegment<'a> {
    t: &'a [f64],
    f_raw: &'a [f64],
    f_norm: Vec<f64>,
    med: f64,
}

fn median_absolute_deviation(arr: &[f64], med: &f64) -> Option<f64> {
    if arr.is_empty() { return None; }

    let mut deviations: Vec<f64> = arr
        .iter()
        .copied()
        .filter(|f| f.is_finite())
        .map(|f| (f - med).abs())
        .collect();

    if deviations.is_empty() { return None; }

    median(&mut deviations)
}

impl<'a> NormSegment<'a> {
    fn segment_noise_clean(&mut self, k: f64){
        let Some(mad) = median_absolute_deviation(&self.f_norm, &self.med) else { return };
        let sigma = 1.4826 * mad;
        let threshold = k * sigma;

        for f in self.f_norm.iter_mut() {
            if !f.is_finite() || f.abs() > threshold {
                *f = f64::NAN;
            }
        }
    }
}

fn normalize_segments<'a>(
    t: &'a [f64],
    f: &'a [f64],
    segments: &[(usize, usize)]
) -> Vec<NormSegment<'a>> {
    let mut out = Vec::with_capacity(segments.len());

    for &(s, e) in segments {
        let f_seg = &f[s..e]; // check later if you need to have f[s..e+1]
        let mut temp: Vec<f64> = f_seg.iter().copied().collect();
        let med = median(&mut temp).expect("empty segment");
        assert!(med.is_finite() && med != 0.0, "segment median zero or infinite; normalization blow up");
        let inv = 1.0 / med;
        let f_norm: Vec<f64> = f_seg.iter().map(|&f| inv * f - 1.0).collect();

        out.push(NormSegment {
            t: &t[s..e],
            f_raw: f_seg,
            f_norm,
            med
        });
    }
    return out;
}

fn get_phases(times: &[f64], period: f64) -> Vec<f64> {

    let mut phases = Vec::new();

    for t in times {
        let phase = (t % period) / period;
        phases.push(phase);
    }

    return phases;
}

fn binning<const N_BINS: usize>(flux: &Vec<f64>, phases: &Vec<f64>) -> Vec<f64> {

    let mut binned_flux: Vec<f64> = Vec::new();

    let mut count: [i32; N_BINS] = [0; N_BINS];
    let mut sum_flux: [f64; N_BINS] = [0.0; N_BINS];

    for (f, p) in flux.iter().zip(phases.iter()) {
        let bin_index = (p * N_BINS as f64).floor() as usize;
        count[bin_index] += 1;
        sum_flux[bin_index] += f;
    }

    for (sum, count) in sum_flux.into_iter().zip(count.into_iter()) {
        // sum_flux and count don't need to retain ownership
        binned_flux.push(sum / count as f64)
    }

    binned_flux
}

fn search_bins<const N_BINS: usize>(
    trial_period: &f64,
    trial_duration: &f64,
    binned_flux: &Vec<f64>,
) -> Option<f64> {

    let window_width: usize = (((trial_duration / 24.0) / trial_period) * N_BINS as f64).round() as usize;

    /*
    println!("trial_duration = {}", trial_duration);
    println!("window_width = {}", window_width);
    println!("binned_flux.len() = {}", binned_flux.len());
     */


    if window_width == 0 || window_width > binned_flux.len() { return None; }

    let mut sum_flux: f64 = binned_flux[..window_width].iter().sum();
    let mut result: f64 = sum_flux / window_width as f64;

    for right in window_width..binned_flux.len() {
        sum_flux += binned_flux[right] - binned_flux[right - window_width];
        let avg = sum_flux / window_width as f64;

        if avg < result {
            result = avg;
        }
    }

    Some(result)
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

    let mut norm_segments: Vec<NormSegment> = normalize_segments(&t, &f, &seg_bounds);

    // cleaning all segments
    for seg in norm_segments.iter_mut() {
        seg.segment_noise_clean(6.0);
    }

    let trial_periods: Vec<f64> = (1..=4000)
        .map(|i| i as f64 * 0.5)
        .collect();
    // println!("{:#?}", trial_periods);


    // start with just one segment for now...
    let t_n = norm_segments[0].t.to_vec();
    let f_n = norm_segments[0].f_norm.clone(); // should impl copy later


    // period_phase = (period, corresponding phases)
    let mut period_phase: Vec<(f64, Vec<f64>)> = Vec::new();
    for &period in trial_periods.iter() {
        let phases = get_phases(&t_n, period);
        period_phase.push((period, phases));
    }



    /*
    println!("Plotting flux vs time");
    plot_to_python("time".into(), "flux".into(), &t_n, &f_n)?;
    println!("Plotting flux vs phase");
    plot_to_python("phase".into(), "flux".into(), &period_phase[0].1, &f_n)?;
    println!("Done both plots");
     */

    const N_BINS: usize = 200;
    let binned_flux: Vec<f64> = binning::<N_BINS>(&f_n, &period_phase[0].1);
    assert_eq!(binned_flux.len(), N_BINS);


    println!("Plotting binned_flux");
    let binned_range: Vec<f64> = (0..N_BINS).map(|i| i as f64).collect();
    plot_to_python("Bins".into(), "Binned Flux".into(), &binned_range, &binned_flux);


    let trial_durations: Vec<f64> = (1..20)
        .map(|i| i as f64 * 0.5)
        .collect();

    let mut deepest_dip: f64 = f64::INFINITY;

    for d in trial_durations.iter() {
        let box_average = search_bins::<N_BINS>(&period_phase[0].0, d, &binned_flux);
        let mut temp: f64 = 0.0;
        if box_average.is_some() {
            temp = box_average.unwrap();
            // println!("box_average = {}", temp);
            if temp < deepest_dip {
                deepest_dip = temp;
            }
        }
    }

    println!("Deepest dip: {}", deepest_dip);

    Ok(())
}


