pub mod data_access;
pub mod statistical_methods;

use data_access::{fits_access, data_cleaner};

fn main() {
    println!("in main.rs");
    // fits_access::test::testing_fits();
    let file_name = "tess2025180145000-s0094-0000000006286534-0291-s_lc.fits";
    let (t, f) = fits_access::get_data(file_name).expect("Error getting data");
    const K: f64 = 6.0;
    let f_clean = data_cleaner::clean_data(f, &t, K);
}