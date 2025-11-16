/*
    Cameron Emmanuel
    Nov - 05 - 2025
 */

use fitsio::FitsFile;
use fitsio::tables::Column;
use std::path::PathBuf;
use plotters::prelude::*;
use plotters::series::LineSeries;

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
        let Some(mad) = mad_segments_norm(&self.f_norm) else { return };
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

fn binning_flux(
    trial_period: &f64,
    flux_segment: &[f64],
    time_segment: &[f64],
) -> Vec<f64> {
    let mut phase: Vec<f64> = Vec::new();

    for time in time_segment.iter() {
        phase.push( (time % trial_period) / trial_period );
    }

    let n_bins = 200;
    let start = 0.0;
    let end = 1.0;

    let divided_phase: Vec<f64> = (0..=n_bins)
        .map(|i| start + (end - start) * i as f64 / n_bins as f64)
        .collect();

    let mut count: [i32; 200] = [0; 200];
    let mut sum_flux: [f64; 200] = [0.0; 200];

    for phase in divided_phase.iter() {
        let bin_index: usize = (phase * n_bins.as_f64()).floor() as usize;
        sum_flux[bin_index] += flux_segment[bin_index];
        count[bin_index] += 1;
    }

    let binned_flux: Vec<f64> = (0..=n_bins).map(|i| sum_flux[i] / count[i] as f64).collect();

    return binned_flux;
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

    let trial_periods: Vec<f64> = (1..=40)
        .map(|i| i as f64 * 0.5)
        .collect();





    Ok(())
}
