use fitsio::{FitsFile, errors};
use std::{fs, path::Path, path::PathBuf};

#[derive(Debug)]
pub enum FitsError {
    NotFitsFile,
    Fitsio(errors::Error),
    Io(std::io::Error),
}

impl From<errors::Error> for FitsError {
    fn from(err: errors::Error) -> Self {
        FitsError::Fitsio(err)
    }
}

impl From<std::io::Error> for FitsError {
    fn from(e: std::io::Error) -> Self {
        FitsError::Io(e)
    }
}

fn get_fits_ptr(file_name: &str) -> Result<FitsFile, FitsError> {
    let p = Path::new(file_name);
    if p.extension().and_then(|ext| ext.to_str()) != Some("fits") {
        return Err(FitsError::NotFitsFile);
    }

    let fits_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join("fits_files")
        .join(file_name);
    let fptr = FitsFile::open(fits_path).map_err(FitsError::Fitsio)?;
    Ok(fptr)
}

fn load_tess_data(fptr: &mut FitsFile) -> Result<(Vec<f64>, Vec<f64>, Vec<i32>), FitsError> {
    let hdu = fptr.hdu(1)?;
    let time: Vec<f64> = hdu.read_col(fptr, "TIME")?;
    let flux: Vec<f64> = hdu.read_col(fptr, "PDCSAP_FLUX")?;
    let qual: Vec<i32> = hdu.read_col(fptr, "QUALITY")?;
    Ok((time, flux, qual))
}

fn filter_and_sort_tess_data(
    flux: Vec<f64>,
    time: Vec<f64>,
    qual: Vec<i32>,
) -> Result<(Vec<f64>, Vec<f64>), FitsError> {
    let good_indices: Vec<usize> = (0..time.len())
        .filter(|&i| qual[i] == 0 && time[i].is_finite() && flux[i].is_finite())
        .collect();

    let mut paired: Vec<(f64, f64)> = good_indices.iter().map(|&i| (time[i], flux[i])).collect();

    paired.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let (t_sorted, f_sorted): (Vec<f64>, Vec<f64>) = paired.into_iter().unzip();

    Ok((t_sorted, f_sorted))
}

pub fn get_data(file_name: &str) -> Result<(Vec<f64>, Vec<f64>), FitsError> {
    let mut fptr = get_fits_ptr(file_name).expect("Got an Error!");
    let (t, f, q) = load_tess_data(&mut fptr).expect("Error loading data");

    let (t, f) = filter_and_sort_tess_data(f, t, q).expect("Error filtering data");
    Ok((t, f))
}

pub fn fits_attainment() -> Result<Vec<String>, FitsError> {
    println!("Welcome to Transit Detector By Cammy CORP!");
    println!("Please place all fits files within the data/fits_files directory");

    let mut fits_names: Vec<String> = Vec::new();

    let origin = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join("fits_files");

    for entry in fs::read_dir(origin)? {
        let entry = entry?;
        let path = entry.path();

        if path.extension().and_then(|e| e.to_str()) != Some("fits") {
            return Err(FitsError::NotFitsFile);
        }
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            fits_names.push(name.to_string());
        }
    }

    Ok(fits_names)
}

#[cfg(test)]
pub mod test {

    use crate::fits_access;

    #[test]
    pub fn testing_fits() {
        let file_name = "tess2025180145000-s0094-0000000006286534-0291-s_lc.fits";
        let mut fptr = fits_access::get_fits_ptr(file_name).expect("Got an Error!");
        let (t, f, q) = fits_access::load_tess_data(&mut fptr).expect("Error loading data");
        assert_eq!(t.len(), f.len(), "time/flux length mismatch");
        assert_eq!(t.len(), q.len(), "time/quality length mismatch");

        let (t, f) = fits_access::filter_and_sort_tess_data(f, t, q).expect("Error filtering data");
        assert_eq!(t.len(), f.len(), "time/flux length mismatch");
        assert!(
            t.iter().all(|&val| !val.is_nan() && val.is_finite()),
            "NaN or infinite in t"
        );
        assert!(f.iter().all(|&val| !val.is_nan()), "NaN detected in f");
        assert!(t.windows(2).all(|w| w[1] > w[0]), "t is not chronological");
    }
}
