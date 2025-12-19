pub mod data_access;

use data_access::fits_access;

fn main() {
    println!("in main.rs");
    fits_access::test::testing();
}