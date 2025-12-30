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

pub fn median_cadence(times: &[f64]) -> Option<f64> {
    if times.len() < 2 {
        return None;
    }
    let mut dts: Vec<f64> = times
        .windows(2)
        .map(|w| w[1] - w[0])
        .filter(|dt| dt.is_finite() && *dt > 0.0)
        .collect();
    if dts.is_empty() {
        return None;
    }
    median(&mut dts)
}

pub fn median_absolute_deviation(arr: &[f64]) -> Option<f64> {
    if arr.is_empty() {
        return None;
    }

    let med = median(arr)?;
    let mut deviations: Vec<f64> = arr
        .iter()
        .copied()
        .filter(|f| f.is_finite())
        .map(|f| (f - med).abs())
        .collect();

    if deviations.is_empty() {
        return None;
    }
    median(&mut deviations)
}
