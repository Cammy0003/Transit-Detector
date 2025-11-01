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

    let t_good: Vec<f64> = good_indices.iter().map(|&i| time[i]).collect();
    let f_good: Vec<f64> = good_indices.iter().map(|&i| flux[i]).collect();

    let mut paired: Vec<(f64, f64)> = t_good.iter().cloned().zip(f_good.iter().cloned()).collect();
    paired.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let (t_sorted, f_sorted): (Vec<f64>, Vec<f64>) = paired.into_iter().unzip();

    Ok((t_sorted, f_sorted))
}

fn median(v: &mut Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
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
    segment_bounds
}

fn main() -> fitsio::errors::Result<()> {
    let fits_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join("test")
        .join("tess2025180145000-s0094-0000000006286534-0291-s_lc.fits");

    let mut fptr = FitsFile::open(fits_path)?;
    // fptr.pretty_print()?;

    let (t, f, q) = load_tess_data(&mut fptr)?;
    let (t, f) = filter_and_sort_tess_data(f, t, q)?;
    assert!(t.iter().all(|&val| !val.is_nan() && val.is_finite()), "NaN or infinite in t");
    assert!(f.iter().all(|&val| !val.is_nan()), "NaN detected in t");
    assert!(t.windows(2).all(|w| w[1] > w[0]), "t is not chronological");

    let dt_med = median_cadence(&t).expect("median cadence not found");
    let gap_factor = 5.0;
    let seg_bounds = segment_on_gaps(&t, dt_med, gap_factor); // dt_med loses ownership

    Ok(())
}
