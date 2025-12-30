pub mod candidacy;
pub mod data_access;
pub mod statistical_methods;

use data_access::{CleanLightCurve, data_cleaner, fits_access};

use candidacy::finding_candidates::{Candidate, find_candidates, transit_estimate, trial_periods};

fn main() {
    let fits_names = fits_access::fits_attainment().expect("Error getting fits");
    if fits_names.is_empty() {
        println!("No Fits Files");
        return;
    }
    let file_name = &fits_names[0];
    let (t, f) = fits_access::get_data(&file_name).expect("Error getting data");
    const K: f64 = 6.0;
    let (f_clean, sigma) = data_cleaner::clean_data(f, &t, K).expect("Error cleaning data");

    let light_curve: CleanLightCurve = CleanLightCurve::new(t, f_clean, sigma);

    let num_periods = 250;
    let num_bins = 200;
    let num_dur = 20;

    let t_periods = trial_periods(num_periods);

    let candidates: Vec<Candidate> = find_candidates(&light_curve, &t_periods, num_bins, num_dur);

    let transit: Candidate = transit_estimate(&candidates);

    println!("\ntransit estimate:\n {}", transit);
}
