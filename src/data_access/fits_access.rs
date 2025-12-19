use fitsio::{FitsFile, errors};

pub fn load_tess_data(fptr: &mut FitsFile) -> errors::Result<(Vec<f64>, Vec<f64>, Vec<i32>)> {
    Ok((Vec::new(), Vec::new(), Vec::new()))
}

pub fn hello() -> bool {
    return true;
}

#[cfg(test)]
pub mod test {

    use crate::fits_access;
    #[test]
    pub fn testing() {
        if fits_access::hello() {
            println!("Hello, world!");
        } else {
            println!("Not Hello, world!");
        }
    }
}