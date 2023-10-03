use num::complex::{Complex64};
use rayon::prelude::*;
#[derive(Debug)]
pub struct BootstrapData {
    sample_0: Complex64,
    bootstrap_samples: Vec<Complex64>
}

impl BootstrapData {

    pub fn boot_average(&self) -> Complex64 {
        self.bootstrap_samples.par_iter().sum::<Complex64>() / (self.bootstrap_samples.len() as f64)
    }

    pub fn get_sample_0(&self) -> Complex64 {
        self.sample_0
    }

    pub fn boot_error(&self) -> Complex64 {
        let boot_average = self.boot_average();
        let sample_size: f64 = self.bootstrap_samples.len() as f64;
        let res_sum: (f64,f64) = self.bootstrap_samples
            .par_iter()
            .map(|e|(
                (e.re - boot_average.re).powi(2)/sample_size,
                (e.im - boot_average.im).powi(2)/sample_size
            ))
            .reduce(||{(0.0,0.0)}, |(a,b), (c,d)|{
                (a+c, b+d)
            });
        Complex64::new(res_sum.0.sqrt(), res_sum.1.sqrt())
    }

}

impl From<Vec<Complex64>> for BootstrapData{
    fn from(value: Vec<Complex64>) -> Self {
        let sample_0 = value.into_iter().next().expect("Empty vector while creating bootstrap data");
        let bootstrap_samples = value.collect();
        BootstrapData { sample_0, bootstrap_samples }
    }
}