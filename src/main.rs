pub mod candidacy;
pub mod data_access;
pub mod statistical_methods;

use data_access::{CleanLightCurve, data_cleaner, fits_access};

use candidacy::finding_candidates::{Candidate, find_candidates, transit_estimate, trial_periods};

fn main() {
    // println!("in main.rs");
    // fits_access::test::testing_fits();
    let file_name = "tess2025180145000-s0094-0000000006286534-0291-s_lc.fits";
    let (t, f) = fits_access::get_data(file_name).expect("Error getting data");
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
