use crate::statistical_methods::statistics::{median, median_absolute_deviation, median_cadence};

pub struct NormSegment {
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

pub fn normalize_segments(f: &[f64], segments: &[(usize, usize)]) -> Vec<NormSegment> {
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

pub fn clean_data(flux: Vec<f64>, times: &Vec<f64>, k: f64) -> Vec<f64> {
    let dt_med = median_cadence(&times).expect("median cadence not found");
    let gap_factor = 5.0;
    let seg_bounds = segment_on_gaps(&times, dt_med, gap_factor); // dt_med loses ownership

    let mut norm_segments: Vec<NormSegment> = normalize_segments(&flux, &seg_bounds);
    for seg in norm_segments.iter_mut() {
        seg.segment_noise_clean(k);
    }
    let f_all: Vec<f64> = norm_segments
        .iter()
        .flat_map(|seg| seg.f_norm.iter().copied())
        .collect();

    f_all
}
